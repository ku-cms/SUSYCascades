import ROOT, os, subprocess
from collections import OrderedDict

def collect_unique_hashes(input_files):
    hashes = set()
    for f in input_files:
        tf = ROOT.TFile.Open(f)
        meta = tf.GetDirectory("meta")
        if meta:
            tag = meta.Get("GitCommitHash")
            if tag:
                hashes.add(tag.GetTitle())
        tf.Close()
    return sorted(hashes)

def write_git_hashes_to_output(output_path, hashes):
    tf = ROOT.TFile.Open(output_path, "UPDATE")
    if not tf.GetDirectory("meta"):
        tf.mkdir("meta")
    tf.cd("meta")
    for h in hashes:
        ROOT.TNamed(f"GitHash_{h[:8]}", h).Write()
    tf.Close()

def has_uncommitted_changes(repo_dir=".", sub_dir=None):
    """
    Checks if a Git repo (optionally limited to a subdirectory) has uncommitted changes.

    Args:
        repo_dir (str): Path to the Git repository.
        sub_dir (str, optional): Subdirectory inside the repo to check (e.g., 'src/' or 'include/').

    Returns:
        bool: True if there are uncommitted changes, False otherwise.
    """
    if not os.path.isdir(repo_dir):
        raise FileNotFoundError(f"Directory does not exist: {repo_dir}")

    # Verify this is a git repository
    try:
        subprocess.run(
            ["git", "-C", repo_dir, "rev-parse", "--is-inside-work-tree"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError:
        raise RuntimeError(f"{repo_dir} is not a Git repository.")

    # Build git status command
    cmd = ["git", "-C", repo_dir, "status", "--porcelain"]
    if sub_dir:
        cmd.append("--")
        cmd.append(sub_dir)

    # Run git status
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    changed_files = [line.split(maxsplit=1)[1] for line in result.stdout.strip().splitlines() if line]

    if changed_files:
        print("Changed files in:",sub_dir)
        for f in changed_files:
            print(f"  {f}")
        return True
    else:
        return False

