import os, argparse, subprocess, time

def process_DASCheck_event_count(input_dir, working_dir, output_dir, dryrun, file, directory):
    processes = []
    if dryrun:
        print("Doing dryrun")
    if file:
        files = [d for d in os.listdir(input_dir) if "EventCount_NANO_" in d]
        for file in files:
            tag = file.replace("EventCount_NANO_", "").replace(".root", "")
            command = f"python3 python/new_countEvents.py --eventCount --filetag {tag}"
            os.makedirs("DASCheck_logs/",exist_ok=True)
            output_file = f"DASCheck_logs/DASCheck_{tag}.txt"
            if not dryrun:
                print(f"Executing: {command}")
                with open(output_file, "w") as out:
                    p = subprocess.Popen(command, shell=True, stdout=out, stderr=subprocess.STDOUT)
                    processes.append(p)
            else:
                print(f"Command that would be executed: {command}")
    
    elif directory:
        dirs = os.listdir(working_dir)
        for EC_dir in dirs:
            if EC_dir.endswith('_EventCount'):
                command = f"python3 python/CheckFiles.py -d {EC_dir} -o {output_dir}{EC_dir} --eventCount"
                os.makedirs("DASCheck_logs/",exist_ok=True)
                output_file = f"DASCheck_logs/DASCheck_{EC_dir}.txt"
                if not dryrun:
                    print(f"Executing: {command}")
                    with open(output_file, "w") as out:
                        p = subprocess.Popen(command, shell=True, stdout=out, stderr=subprocess.STDOUT)
                        processes.append(p)
                else:
                    print(f"Command that would be executed: {command}")
    
    # Wait for all processes to finish
    for p in processes:
        p.wait()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DAS Check EventCount root files.")
    parser.add_argument("-i", "--input", default=f"{os.getenv('CMSSW_BASE')}/src/SUSYCascades/root/EventCount/",
                        help="Input directory for EventCount files.")
    parser.add_argument("-w", "--working", default=f"{os.getenv('PWD')}",
                        help="Working directory containing EventCount subdirectories.")
    parser.add_argument("-o", "--output", default=f"/ospool/cms-user/{os.getenv('USER')}/EventCount/root/",
                        help="Output directory for condor jobs.")
    parser.add_argument("--dryrun", action='store_true', help="Do not run countEvents check")
    parser.add_argument("-d", action='store_true', help="Check directory of output root files from condor jobs")
    parser.add_argument("-f", action='store_true', help="Check final hadded root files")
    
    args = parser.parse_args()
    
    time.sleep(1)
    if not args.d and not args.f:
        print("NEED TO SPECIFY WHETHER RUNNING OVER FILES OR DIRECTORY OF FILES")
    else:
        process_DASCheck_event_count(args.input, args.working, args.output, args.dryrun, args.f, args.d)


