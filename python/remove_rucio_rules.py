#!/usr/bin/env python3
import argparse
import subprocess
import sys
import re
import os

RULE_ID_RE = re.compile(r"^[a-f0-9]{32}$")

def get_account():
    try:
        return subprocess.check_output(["whoami"], text=True).strip()
    except Exception:
        print("ERROR: Could not determine user account")
        sys.exit(1)

def list_rules(account, env):
    cmd = ["rucio", "rule", "list", "--account", account]
    try:
        output = subprocess.check_output(cmd, text=True, env=env)
    except subprocess.CalledProcessError as e:
        print("ERROR: Failed to list rucio rules")
        print(e)
        sys.exit(1)

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

def remove_rule(rule_id, dry_run, env):
    if dry_run:
        print(f"[DRY-RUN] Would remove rule {rule_id}")
        return 0
    cmd = ["rucio", "rule", "remove", rule_id]
    return subprocess.call(cmd, env=env)

def main():
    parser = argparse.ArgumentParser(
        description="Delete all Rucio rules for the current user"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print rules that would be deleted without removing them"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only delete the first N rules (useful for testing)"
    )
    parser.add_argument(
        "--account",
        type=str,
        default=None,
        help="Rucio account name (will default to $USER)"
    )

    args = parser.parse_args()

    account = args.account
    if not account:
        account = get_account()

    # Setup rucio environment
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

    rules = list_rules(account, rucio_env)

    if not rules:
        print("No Rucio rules found.")
        return 0

    print(f"Found {len(rules)} rules")

    if args.limit is not None:
        rules = rules[:args.limit]
        print(f"Limiting to first {len(rules)} rules")

    for rid in rules:
        ret = remove_rule(rid, args.dry_run, rucio_env)
        if ret != 0 and not args.dry_run:
            print(f"WARNING: Failed to remove rule {rid}")

if __name__ == "__main__":
    main()

