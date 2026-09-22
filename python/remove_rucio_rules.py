#!/usr/bin/env python3
"""
Delete Rucio rules for the current user.

Normally this removes every rule owned by the account. If --container is
given, that container's own rule is left active (so the container keeps
being replicated to wherever it lives, e.g. the LPC T3), and instead the
contents (datasets/files) attached to the container are detached, which
empties the container without touching its rule.

Zach's container for LPC: --container user.zflowers:/Analyses/Cascades_dataset/USER#datasets
"""
import argparse
import subprocess
import sys
import re
import os

RULE_ID_RE = re.compile(r"^[a-f0-9]{32}$")
DID_RE = re.compile(r"^\S+:\S+$")

def get_account():
    try:
        return subprocess.check_output(["whoami"], text=True).strip()
    except Exception:
        print("ERROR: Could not determine user account")
        sys.exit(1)

def setup_rucio_env():
    rucio_env = os.environ.copy()
    rucio_setup = "source /cvmfs/cms.cern.ch/rucio/setup-py3.sh && env"
    out = subprocess.check_output(["bash", "-c", rucio_setup], text=True)
    for line in out.splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        rucio_env[k] = v
    for k in list(rucio_env.keys()):
        if k.startswith("BASH_FUNC_"):
            del rucio_env[k]
    return rucio_env

def _extract_rule_ids(output):
    rule_ids = []
    for line in output.splitlines():
        if not line.strip():
            continue
        first = line.split()[0]
        # Skip header and separator lines
        if not RULE_ID_RE.match(first):
            continue
        rule_ids.append(first)
    return rule_ids

def list_rules(account, env):
    """All rule IDs owned by the account."""
    cmd = ["rucio", "rule", "list", "--account", account]
    try:
        output = subprocess.check_output(cmd, text=True, env=env)
    except subprocess.CalledProcessError as e:
        print("ERROR: Failed to list rucio rules")
        print(e)
        sys.exit(1)
    return _extract_rule_ids(output)

def get_container_rule_ids(container_did, env):
    """Rule ID(s) attached directly to the container DID (to protect)."""
    cmd = ["rucio", "rule", "list", "--did", container_did]
    try:
        output = subprocess.check_output(cmd, text=True, env=env)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to list rules for container {container_did}")
        print(e)
        sys.exit(1)
    return _extract_rule_ids(output)

def list_container_contents(container_did, env):
    """DIDs currently attached to the container."""
    cmd = ["rucio", "did", "content", "list", container_did]
    try:
        output = subprocess.check_output(cmd, text=True, env=env)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to list contents of {container_did}")
        print(e)
        sys.exit(1)

    contents = []
    for line in output.splitlines():
        line = line.strip()
        # table rows look like: | scope:name | TYPE |
        if not line.startswith("|"):
            continue
        cols = [c.strip() for c in line.strip("|").split("|")]
        if not cols:
            continue
        candidate = cols[0]
        if candidate == "SCOPE:NAME":
            continue
        if not DID_RE.match(candidate):
            continue
        contents.append(candidate)
    return contents

def detach_from_container(container_did, content_did, dry_run, env):
    if dry_run:
        print(f"[DRY-RUN] Would detach {content_did} from {container_did}")
        return 0
    cmd = ["rucio", "did", "content", "remove", "--to-did", container_did, content_did]
    return subprocess.call(cmd, env=env)

def remove_rule(rule_id, dry_run, env):
    if dry_run:
        print(f"[DRY-RUN] Would remove rule {rule_id}")
        return 0
    cmd = ["rucio", "rule", "remove", rule_id]
    return subprocess.call(cmd, env=env)

def main():
    parser = argparse.ArgumentParser(
        description="Delete Rucio rules for the current user, optionally "
                     "keeping one container's rule alive and emptying its "
                     "contents instead of removing that rule."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be removed/detached without doing it",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only delete the first N non-protected rules (useful for testing)",
    )
    parser.add_argument(
        "--account",
        type=str,
        default=None,
        help="Rucio account name (will default to $USER)",
    )
    parser.add_argument(
        "--container",
        type=str,
        default=None,
        help="Container DID (scope:name) whose rule should be kept active, "
             "e.g. user.zflowers:/Analyses/Cascades_dataset/USER. Its "
             "contents are detached instead of the rule being removed.",
    )
    parser.add_argument(
        "--keep-container-contents",
        action="store_true",
        help="With --container, only protect its rule; don't detach its "
             "contents either",
    )

    args = parser.parse_args()

    account = args.account or get_account()
    rucio_env = setup_rucio_env()

    protected_rule_ids = set()
    if args.container:
        protected_rule_ids = set(get_container_rule_ids(args.container, rucio_env))
        if not protected_rule_ids:
            print(f"WARNING: No rule found for container {args.container}; "
                  "nothing to protect")
        else:
            print(f"Protecting rule(s) on container {args.container}: "
                  f"{', '.join(sorted(protected_rule_ids))}")

    rules = list_rules(account, rucio_env)

    if not rules:
        print("No Rucio rules found.")
    else:
        to_remove = [r for r in rules if r not in protected_rule_ids]
        skipped = len(rules) - len(to_remove)

        print(f"Found {len(rules)} rules")
        if skipped:
            print(f"Keeping {skipped} rule(s) tied to the container active")

        if args.limit is not None:
            to_remove = to_remove[: args.limit]
            print(f"Limiting to first {len(to_remove)} rules")

        for rid in to_remove:
            ret = remove_rule(rid, args.dry_run, rucio_env)
            if ret != 0 and not args.dry_run:
                print(f"WARNING: Failed to remove rule {rid}")

    if args.container and not args.keep_container_contents:
        contents = list_container_contents(args.container, rucio_env)
        if not contents:
            print(f"No contents found in container {args.container}.")
        else:
            print(f"Found {len(contents)} item(s) attached to {args.container}")
            for cdid in contents:
                ret = detach_from_container(args.container, cdid, args.dry_run, rucio_env)
                if ret != 0 and not args.dry_run:
                    print(f"WARNING: Failed to detach {cdid} from {args.container}")

    return 0

if __name__ == "__main__":
    sys.exit(main())