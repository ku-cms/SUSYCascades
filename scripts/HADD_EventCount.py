import os
import argparse

def process_event_count_dirs(input_dir, output_dir):
    # Ensure output directory exists
    os.makedirs(os.path.join(output_dir, "root", "EventCount"), exist_ok=True)
    
    # List all subdirectories in input_dir that contain '_EventCount'
    subdirs = [d for d in os.listdir(input_dir) if "_EventCount" in d and os.path.isdir(os.path.join(input_dir, d))]
    
    for subdir in subdirs:
        # Remove '_EventCount' from the subdirectory name
        subdir_cleaned = subdir.replace("_EventCount", "")
        
        # Construct input and output paths
        input_path = os.path.join(input_dir, subdir, "*", "*.root")
        output_path = os.path.join(output_dir, f"EventCount_NANO_{subdir_cleaned}.root")
        
        # Construct and execute hadd command
        command = f"hadd -j 4 -f {output_path} {input_path}"
        print(f"Executing: {command}")
        os.system(command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process EventCount directories and merge ROOT files.")
    parser.add_argument("-i", "--input", default=f"/ospool/cms-user/{os.getenv('USER')}/EventCount/root/", help="Input directory")
    parser.add_argument("-o", "--output", default=f"{os.getenv('CMSSW_BASE')}/src/SUSYCascades/root/EventCount/", help="Output directory")
    args = parser.parse_args()
    
    process_event_count_dirs(args.input, args.output)

