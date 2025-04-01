import time
import subprocess
class CondorJobCountMonitor: 
    def __init__(self, threshold=10000, verbose=False):
        """Initialize the job monitor with a threshold value."""
        self.threshold = threshold
        if self.threshold < 1:
            self.threshold = 1
        self.verbose = verbose
    def get_total_jobs(self):
        """Returns the total number of jobs in the HTCondor scheduler for the current user."""
        try:
            output = subprocess.check_output(f"condor_q $USER -total", shell=True, text=True)
            total = 0
            for line in output.split("\n"):
                if "Total for query" in line:
                    total += int(line.split(' ')[3])  # Extract total job count
            return total
        except Exception as e:
            print(f"Error retrieving job count: {e}")
        return -1  # Return -1 if an error occurs
    
    def wait_until_jobs_below(self):
        """Waits until the number of running jobs is below the given threshold."""
        check_count = 0
        while True:
            total_jobs = self.get_total_jobs()
            if total_jobs == -1:
                print("Error retrieving job count, retrying...")
            elif total_jobs < self.threshold:
                if self.verbose: print(f"Job count ({total_jobs}) is below the threshold ({self.threshold}). Proceeding...")
                break
            else:
                if check_count % 10 == 0 and self.verbose:
                    print(f"Current jobs: {total_jobs}. Waiting for jobs to drop below {self.threshold}...")
            check_count += 1
            time.sleep(5) # Wait before checking again

