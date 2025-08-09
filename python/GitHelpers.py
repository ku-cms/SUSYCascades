import ROOT, os, subprocess
from collections import OrderedDict

def collect_unique_git_hashes(input_files):
    unique_hashes = OrderedDict()
    for fpath in input_files:
        tf = ROOT.TFile.Open(fpath, "READ")
        if not tf or tf.IsZombie():
            continue
        tree = tf.Get("meta/gitHash")
        if not tree:
            tf.Close()
            continue
        for entry in tree:
            h = entry.gitHash
            if h not in unique_hashes:
                unique_hashes[h] = True
        tf.Close()
    return list(unique_hashes.keys())

def overwrite_meta_git_hashes(output_path, unique_hashes):
    fout = ROOT.TFile.Open(output_path, "UPDATE")
    if not fout:
        raise RuntimeError(f"Cannot open {output_path} for update")

    # Remove old meta/gitHashes tree if exists
    meta_dir = fout.Get("meta")
    if not meta_dir:
        meta_dir = fout.mkdir("meta")
    else:
        meta_dir.Delete("gitHashes;*")
    meta_dir.cd()

    # Create new gitHashes tree
    tree = ROOT.TTree("gitHashes", "Git commit hashes")
    from ctypes import create_string_buffer, c_char_p, c_char

    # Use ROOT std::string branch to store hashes
    str_var = ROOT.std.string()
    tree.Branch("gitHash", str_var)

    for h in unique_hashes:
        str_var.replace(0, ROOT.std.string.npos, h)
        tree.Fill()

    tree.Write("", ROOT.TObject.kOverwrite)
    fout.Close()

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
        print("Changed files:")
        for f in changed_files:
            print(f"  {f}")
        return True
    else:
        print("Working directory is clean.")
        return False

