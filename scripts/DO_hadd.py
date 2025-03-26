import os, sys, time
from glob import glob as glob
from subprocess import Popen as pop
import subprocess, psutil

# Example submission:
#    nohup python scripts/DO_hadd.py -idir ../../../NTUPLES/Processing/Summer16_102X/ -odir ../../../NTUPLES/HADD/Summer16_102X/ > HADD_logs/HADD_Summer16_102X.debug 2>&1 &
# After hadd finishes and ready to copy to LPC:
#    nohup xrdcp --parallel 4 -f ../../../NTUPLES/HADD/Summer16_102X/* root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v1/Summer16_102X/ > xrdcp_Summer16_102X.debug 2>&1 &

def find_hadd_procs():
    current_user = psutil.Process().username()
    hadd_procs = []
    for proc in psutil.process_iter(attrs=['pid', 'name', 'username', 'cmdline']):
        try:
            if proc.info['username'] == current_user and 'hadd' in proc.info['name']:
                hadd_procs.append(proc.info)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return hadd_procs

def is_target_running(target):
    hadd_procs = find_hadd_procs()
    for proc in hadd_procs:
        cmdline = proc.get('cmdline',None)
        if any(target in arg for arg in cmdline):
            return True
    return False

def main():
    # start time
    start_time = time.time()
    print("DO_hadd.py: Start... go go go!")
    print("------------------------------")

    argv_pos = 1

    OUT_DIR = "dum"
    IN_DIR  = "dum"
    redo    = False
    SKIP_BAD_FILES = False
    SKIP_COPY = False
    SKIP_SMALL = False
    DRY_RUN = False
    
    if '-odir' in sys.argv:
        p = sys.argv.index('-odir')
        OUT_DIR = sys.argv[p+1]
        argv_pos += 2
    if '-idir' in sys.argv:
        p = sys.argv.index('-idir')
        IN_DIR = sys.argv[p+1]
        argv_pos += 2
    if '--redo' in sys.argv:
        redo = True
        argv_pos += 1
    if '--skip-bad-files' in sys.argv:
        SKIP_BAD_FILES = True
        argv_pos += 1
    if '--skip-copy' in sys.argv:
        SKIP_COPY = True
        argv_pos += 1
    if '--skip-small' in sys.argv:
        SKIP_SMALL = True
        argv_pos += 1
    if '--dry-run' in sys.argv or '--dryrun' in sys.argv:
        DRY_RUN = True
        argv_pos += 1

    if not len(sys.argv) > 1 or '-h' in sys.argv or '--help' in sys.argv or OUT_DIR == "dum" or IN_DIR == "dum":
        print(f"Usage: {sys.argv[0]} [-idir /path/input_dir] [-odir /path/output_dir] --redo --skip-bad-files")
        sys.exit(1)

    print(f"Input Directory: {IN_DIR}")
    print(f"Output Directory: {OUT_DIR}")

    # create and organize output folders
    if DRY_RUN:
        print("mkdir -p "+OUT_DIR)
        print("mkdir -p HADD_logs/")
    else:
        os.system("mkdir -p "+OUT_DIR)
        os.system("mkdir -p HADD_logs/")

    skip_list = [
        #"SMS-T2tt_mStop-400to1200_TuneCP2_13TeV-madgraphMLM-pythia8",
    ]
    redo_list = [
        #"TTTT_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X",
    ]

    if not DRY_RUN:
        if os.path.exists("scripts/startup_C.so") is False:
            os.system("cd scripts && root.exe -b -l -q startup.C+ && cd ..")
        os.environ["LD_PRELOAD"] = os.environ["PWD"]+"/scripts/startup_C.so"

    hadd_big_processes = {}
    for target in os.listdir(IN_DIR):
        skip = False
        for dataset in skip_list:
            if dataset in target:
                skip = True
        if skip:
            continue

        if redo and target not in redo_list:
            continue
        if redo and target not in os.listdir(IN_DIR):
            continue

        print(f"Target: {target}")
        hadd_sml_processes = []
        if os.path.exists("HADD_logs/"+target) is True:
            os.system("rm -r HADD_logs/"+target)
        os.system("mkdir -p HADD_logs/"+target)
        os.system(f"mkdir -p {OUT_DIR}/{target}")
        for i in range(0,10):
            os.system("mkdir -p "+OUT_DIR+"/"+target+"/"+target+"_"+str(i))
            if not DRY_RUN and not SKIP_COPY:
                for f in glob(os.path.join(IN_DIR+"/"+target+"/"+target+"_*"+str(i)+".root")):
                    os.system("cp "+f+" "+OUT_DIR+"/"+target+"/"+target+"_"+str(i)+"/") 
            command = ["hadd", "-f", f"{OUT_DIR}/{target}/{target}_{str(i)}.root"] + glob(f"{OUT_DIR}/{target}/{target}_{str(i)}/*.root")
            if SKIP_BAD_FILES:
                command.insert(2, "-k")
            if DRY_RUN:
                print(command)
            elif not SKIP_SMALL:
                hadd_sml_processes.append(pop(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE))

        for i, hadd_sml in enumerate(hadd_sml_processes):
            if not is_target_running(target+str(i)):
                out, err = hadd_sml.communicate()  # Waits for completion & captures output
                err = err.decode("utf-8")  # Convert bytes to string
                if err.strip():  # Check if error message is non-empty
                    filtered_lines = [line for line in err.splitlines() if "TFile::Cp" not in line]
                    if filtered_lines:
                        log_path = f"HADD_logs/{target}/{target}_{i}.err"
                        with open(log_path, "a") as err_log:
                            for line in filtered_lines:
                                err_log.write(line)

        # Limit big hadd running processes to 10
        while len(hadd_big_processes) >= 10:
            del_targets = []
            for current_target, hadd_big in list(hadd_big_processes.items()):
                if not is_target_running(current_target) and hadd_big.poll() is not None:
                    out, err = hadd_big.communicate()
                    err = err.decode("utf-8")  # Convert bytes to string
                    if err.strip():  # Log errors only if non-empty
                        filtered_lines = [line for line in err.splitlines() if "TFile::Cp" not in line]
                        if filtered_lines:
                            log_path = f"HADD_logs/{current_target}.err"
                            print(f"Outputting error to: {log_path}")
                            with open(log_path, "a") as err_log:
                                for line in filtered_lines:
                                    err_log.write(line)
                    del_targets.append(current_target)
                    break  # Exit loop after finding a finished process
            for del_target in del_targets:
                if del_target in hadd_big_processes:
                    del hadd_big_processes[del_target]  # Remove completed process

        # Launch new process
        command = ["hadd", "-f", f"{OUT_DIR}/{target}.root"] + glob(f"{OUT_DIR}/{target}/*.root")
        if DRY_RUN:
            print(command)
        else:
            tmp_pop = pop(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            hadd_big_processes[str(target)] = tmp_pop
    
    # Cleanup remaining hadd processes
    while len(hadd_big_processes) > 0:
        for target, hadd_big in list(hadd_big_processes.items()):
            if not is_target_running(target):
                out, err = hadd_big.communicate()
                err = err.decode("utf-8")
                if err.strip():
                    filtered_lines = [line for line in err.splitlines() if "TFile::Cp" not in line]
                    if filtered_lines:
                        log_path = f"HADD_logs/{target}.err"
                        print(f"Outputting error to: {log_path}")
                        with open(log_path, "a") as err_log:
                            for line in filtered_lines:
                                err_log.write(line)
                del hadd_big_processes[target]  # Remove finished process


    print("Finished Merging Files")    
    print("------------------------------")
    # end time
    end_time = time.time()
    print("DO_hadd.py: End... all done!")

    # total time in seconds
    total_time_seconds = end_time - start_time
    # total time in minutes
    total_time_minutes = total_time_seconds / 60
    # total time in hours
    total_time_hours = total_time_minutes / 60
    
    print("Total time: {0:.2f} seconds = {1:.2f} minutes = {2:.2f} hours".format(total_time_seconds, total_time_minutes, total_time_hours))


if __name__ == "__main__":
    main()
