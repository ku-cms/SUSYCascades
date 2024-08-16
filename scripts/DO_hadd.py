import os, sys, time
from glob import glob as glob
from subprocess import Popen as pop
import subprocess

# Example submission:
#    nohup python scripts/DO_hadd.py -idir ../../../NTUPLES/Processing/Summer16_102X/ -odir ../../../NTUPLES/HADD/Summer16_102X/ > HADD_logs/HADD_Summer16_102X.debug 2>&1 &
# After hadd finishes and ready to copy to LPC:
#    nohup xrdcp --parallel 4 -f ../../../NTUPLES/HADD/Summer16_102X/* root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v1/Summer16_102X/ > xrdcp_Summer16_102X.debug 2>&1 &

def main():
    # start time
    start_time = time.time()
    print("DO_hadd.py: Start... go go go!")
    print("------------------------------")

    argv_pos = 1

    OUT_DIR = "dum"
    IN_DIR = "dum"
    redo = False
    
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

    if not len(sys.argv) > 1 or '-h' in sys.argv or '--help' in sys.argv or OUT_DIR == "dum" or IN_DIR == "dum":
        print(f"Usage: {sys.argv[0]} [-idir /path/input_dir] [-odir /path/output_dir]")
        sys.exit(1)

    print(f"Input Directory: {IN_DIR}")
    print(f"Output Directory: {OUT_DIR}")
        
    # create and organize output folders
    os.system("mkdir -p "+OUT_DIR)
    os.system("mkdir -p HADD_logs/")

    skip_list = [
        #"SMS-T2tt_mStop-400to1200_TuneCP2_13TeV-madgraphMLM-pythia8",
    ]
    redo_list = [
        #"TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X",
        #"TTTT_TuneCP5_13TeV-amcatnlo-pythia8_Fall17_102X",
        #"TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8_Fall17_102X",
    ]

    #if os.path.exists("scripts/startup_C.so") is False:
    #    os.system("cd scripts && root.exe -b -l -q startup.C+ && cd ..")

    #os.environ["LD_PRELOAD"] = os.environ["PWD"]+"/scripts/startup_C.so"
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

        print(target)
        #haddcmd = "hadd -f "+OUT_DIR+"/"+target+".root "
        #for i in range(0,10):
            #os.system("mkdir -p "+IN_DIR+"/"+target+"/"+target+"_"+str(i))
            #os.system("mv "+IN_DIR+"/"+target+"/"+target+"_*"+str(i)+".root "+IN_DIR+"/"+target+"/"+target+"_"+str(i))
            #os.system("LD_PRELOAD=scripts/startup_C.so hadd -f -j 4 "+IN_DIR+"/"+target+"/"+target+"_"+str(i)+".root "+IN_DIR+"/"+target+"/"+target+"_"+str(i)+"/*.root")

        hadd_sml_processes = []
        if os.path.exists("HADD_logs/"+target) is True:
            os.system("rm -r HADD_logs/"+target)
        os.system("mkdir -p HADD_logs/"+target)
        for i in range(0,10):
            os.system("mkdir -p "+IN_DIR+"/"+target+"/"+target+"_"+str(i))
            for f in glob(os.path.join(IN_DIR+"/"+target+"/"+target+"_*"+str(i)+".root")):
                os.system("mv "+f+" "+IN_DIR+"/"+target+"/"+target+"_"+str(i)+"/") 
            hadd_sml_processes.append(pop("hadd -f -k "+IN_DIR+"/"+target+"/"+target+"_"+str(i)+".root "+IN_DIR+"/"+target+"/"+target+"_"+str(i)+"/*.root",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True))

        for hadd_sml in hadd_sml_processes:
            if hadd_sml.poll() is True:
                hadd_sml.wait()
            out,err = hadd_sml.communicate()
            if err != "":
                err_log = open("HADD_logs/"+"/"+target+"/"+target+"_"+str(i)+".err","a")
                err_log.write(str(err))
                err_log.close()
        os.system(f"mkdir -p {OUT_DIR}/{target}")
        os.system(f"mv {IN_DIR}/{target}/{target}_*.root {OUT_DIR}/{target}/")

    #    if len(hadd_big_processes) >= 10:
    #        for target, hadd_big in hadd_big_processes.items():
    #            if hadd_big.poll() is not None:
    #                hadd_big.wait()
    #                out,err = hadd_big.communicate()
    #                if err != "":
    #                    print("Outputting error to: HADD_logs/"+"/"+target+".err")
    #                    err_log = open("HADD_logs/"+"/"+target+".err","a")
    #                    err_log.write(err)
    #                    err_log.close()
    #                del hadd_big_processes[target]
    #            elif len(hadd_big_processes) < 10:
    #                tmp_pop = pop("hadd -f "+OUT_DIR+"/"+target+".root "+IN_DIR+"/"+target+"/*.root",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    #                hadd_big_processes[str(target)] = tmp_pop
    #        
    #    else:
    #        tmp_pop = pop("hadd -f "+OUT_DIR+"/"+target+".root "+IN_DIR+"/"+target+"/*.root",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    #        hadd_big_processes[str(target)] = tmp_pop
        

    #for target, hadd_big in hadd_big_processes.items():
    #    if hadd_big.poll() is not None:
    #        print("Waiting on big hadd job")
    #        hadd_big.wait()
    #        out,err = hadd_big.communicate()
    #        if err != "":
    #            print("Outputting error to: HADD_logs/"+"/"+target+".err")
    #            err_log = open("HADD_logs/"+"/"+target+".err","a")
    #            err_log.write(err)
    #            err_log.close()
    #        if hadd_big.poll() is None:
    #            #hadd_big_processes.pop(target,None)
    #            del hadd_big_processes[target]
    #if len(hadd_big_processes) == 0:
    #    print("Finished Merging Files")
    #else:
    #    print("Note: "+str(len(hadd_big_processes))+" hadd jobs may still be running!")
    
    print("------------------------------")
    # end time
    end_time = time.time()
    print("DO_hadd.py: End... all done!")

    # total time in seconds
    total_time_seconds = end_time - start_time
    # total time in hours
    total_time_hours = total_time_seconds / 3600
    
    print("Total time: {0:.2f} seconds = {1:.2f} hours".format(total_time_seconds, total_time_hours))


if __name__ == "__main__":
    main()
