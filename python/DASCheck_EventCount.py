import os
import argparse

def process_DASCheck_event_count(input_dir,dryrun,resubmit):
    if dryrun:
        print("Doing dryrun")
    # List all filess in input_dir that contain 'EventCount_NANO_'
    files = [d for d in os.listdir(input_dir) if "EventCount_NANO_" in d]
    
    for file in files:
        # Remove 'EventCount_NANO_' and '.root' from the subdirectory name
        tag = file.replace("EventCount_NANO_", "")
        tag = tag.replace(".root", "")
        
        # Construct and execute hadd command
        command = f"python3 python/new_countEvents.py --eventCount --filetag {tag} > DASCheck_{tag}.txt"
        if not dryrun:
            print(f"Executing: {command}")
            os.system(command)
            if(resubmit):
                with open(f"DASCheck_{tag}.txt") as DAS_Result:
                  for line in DAS_Result:
                      if "failed" not in line: continue
                      dataset = line.split(' ')[1]
                      basetag = tag.replace('_SMS','')
                      command = f"condor_submit {tag}_EventCount/src/{dataset}_{basetag}.submit"
                      print(f"Resubmitting {dataset} for {tag}")
                      os.system(command) 
        else:
            print(f"Command that would be executed: {command}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DAS Check EventCount root files.")
    parser.add_argument("-i", "--input", default=f"{os.getenv('CMSSW_BASE')}/src/SUSYCascades/root/EventCount/", help="Input directory")
    parser.add_argument("--dryrun", action='store_true', help="Do not run countEvents check")
    parser.add_argument("-r","--resubmit", action='store_true', help="Do not run countEvents check")
    args = parser.parse_args()
    process_DASCheck_event_count(args.input,args.dryrun,args.resubmit)

