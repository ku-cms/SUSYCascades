import subprocess

class GitHashTool: 

    def __init__(self):
        self.gitHash = self.getHash()

    def getHash(self):
        try:
            result = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True)
            return result.strip()  # removes trailing newline
        except subprocess.CalledProcessError:
            raise RuntimeError("Failed to run git command")
        return ""


