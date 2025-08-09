import os, argparse, glob, time
from GitHelpers import collect_unique_git_hashes, overwrite_meta_git_hashes

def run_command(command,dryrun):
    if not dryrun:
        print(f"Executing: {command}")
        os.system(command)
    else:
        print(f"Command to be executed: {command}")

def process_event_count_dirs(input_dir, output_dir, dryrun):
    if dryrun:
        print("Doing dryrun")
    # Ensure output directory exists
    os.makedirs(os.path.join(output_dir, "root", "EventCount"), exist_ok=True)
    
    excluded_dirs = [] # ["102X"] # preUL
    # List all subdirectories in input_dir that contain '_EventCount'
    # Excluding directories that we want to skip
    subdirs = [d for d in os.listdir(input_dir)
              if "_EventCount" in d and os.path.isdir(os.path.join(input_dir, d))
              and not any(excluded in d for excluded in excluded_dirs)]
    print("Subdirs to process",subdirs)
    
    for subdir in subdirs:
        # Remove '_EventCount' from the subdirectory name
        subdir_cleaned = subdir.replace("_EventCount", "")
        
        # Construct and execute hadd command
        # First hadd intermediate datasets (reduce arg list for final hadd)
        command = f"python3 scripts/DO_hadd.py -idir {input_dir}{subdir} -odir {input_dir}../HADD/{subdir} | tee HADD_logs/HADD_{subdir}.debug"
        run_command(command,dryrun)
        
        # Construct input and output paths
        input_path = os.path.join(input_dir, "../HADD/", subdir, "*.root")
        output_path = os.path.join(output_dir, f"EventCount_NANO_{subdir_cleaned}.root")
        
        # Construct and execute hadd command
        # Second hadd for entire filetag
        command = f"hadd -v 1 -f {output_path} {input_path}"
        run_command(command,dryrun)
        if not dryrun:
            unique_hashes = collect_unique_git_hashes(glob.glob(input_path))
            overwrite_meta_git_hashes(output_path, unique_hashes)

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Process EventCount directories and merge ROOT files.")
    parser.add_argument("-i", "--input", default=f"/ospool/cms-user/{os.getenv('USER')}/EventCount/root/", help="Input directory")
    parser.add_argument("-o", "--output", default=f"{os.getenv('CMSSW_BASE')}/src/SUSYCascades/root/EventCount/", help="Output directory")
    parser.add_argument("--dryrun", action='store_true', help="Do not run commands, just print to terminal")
    args = parser.parse_args()
    
    process_event_count_dirs(args.input, args.output, args.dryrun)

    print("Finished Merging All Files") 
    print("------------------------------")
    # end time
    end_time = time.time()

    # total time in seconds
    total_time_seconds = end_time - start_time
    # total time in minutes
    total_time_minutes = total_time_seconds / 60
    # total time in hours
    total_time_hours = total_time_minutes / 60
    
    print("Total for all files time: {0:.2f} seconds = {1:.2f} minutes = {2:.2f} hours".format(total_time_seconds, total_time_minutes, total_time_hours))

