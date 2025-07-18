import subprocess, os
from collections import Counter

def count_jobs_by_status(status):
    result = subprocess.run(
        ["condor_q", "-all", "-format", "%s \n", "Owner", f"-{status}"],
        stdout=subprocess.PIPE,
        text=True
    )
    users = Counter(line.strip() for line in result.stdout.splitlines())
    return users

def print_job_summary(title, users):
    sorted_users = sorted(users.items(), key=lambda x: x[1], reverse=True)
    formatted = ', '.join(f"{user}: {count}" for user, count in sorted_users)
    print(f"{title}: {formatted}")

print("Checking users running jobs")
all_jobs = count_jobs_by_status("idle") + count_jobs_by_status("run") + count_jobs_by_status("held")
print_job_summary("Total Number Of Jobs By User", all_jobs)

idle_jobs = count_jobs_by_status("idle")
print_job_summary("Total Number Of Idle Jobs By User", idle_jobs)

held_jobs = count_jobs_by_status("held")
print_job_summary("Total Number Of Held Jobs By User", held_jobs)

running_jobs = count_jobs_by_status("run")
print_job_summary("Total Number Of Running Jobs By User", running_jobs)

user = os.environ['USER']
user_running_jobs = running_jobs.get(user, 0)
total_running_jobs = sum(running_jobs.values())

if total_running_jobs > 0:
    percent = round((user_running_jobs / total_running_jobs) * 100., 2)
    print(f"{user} owns {percent}% of running jobs")
