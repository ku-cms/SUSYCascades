// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <dirent.h>
#include <vector>

// ROOT includes
#include "TSystem.h"

#include "FitInputBuilder.hh"
#include "SampleTool.hh"
#include "CategoryTool.hh"

using namespace std;

void WriteScript(const string& src_name,
		 const string& log_name,
		 const string& command);

void WriteScriptConnect(const string& src_name,
		 const string& log_name,
		 const string& command,
		 const string& OutputFold);

void MakeConnectTarball(const string& OutputFold);

int main(int argc, char* argv[]) {
  int  maxN = 1;
  bool dryRun = false;
  bool connect = false;
  bool doSigFile = false;
  string SigFile = "";
  bool sigProc = false;  
  //string NtuplePath = "root://xrootd.unl.edu//store/user/zflowers/crogan/";
  string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";
  string OutFile    = "BuildFitInput_output.root";

  bool bprint = false;
  int  year   = 2017;
  bool addBkg  = false;
  bool addSig  = false;
  bool addData = false;
  bool extrahist = false;
  vector<string> proc_to_add;

  CategoryTool CT;
  CategoryList Categories;

  bool setLumi = false;
  double lumi;

  string BuildFitInputCmd;
  
  for(int i = 0; i < argc; i++){
    if(strncmp(argv[i],"--help", 6) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-h", 2) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-path", 5) == 0){
      i++;
      NtuplePath = string(argv[i]);
      BuildFitInputCmd += "-path " + NtuplePath + " ";
    }
    if(strncmp(argv[i],"-o", 2) == 0){
      i++;
      OutFile = string(argv[i]);
    }
    if(strncmp(argv[i],"--output", 8) == 0){
      i++;
      OutFile = string(argv[i]);
    }
    if(strncmp(argv[i],"-year", 5) == 0){
      i++;
      year = std::atoi(argv[i]);
      BuildFitInputCmd += "-year " + string(argv[i]) + " ";
    }
    if(strncmp(argv[i],"+proc", 5) == 0){
      i++;
      proc_to_add.push_back(argv[i]);
    }
    if(strncmp(argv[i],"++bkg", 5) == 0){
      addBkg = true;
    }
    if(strncmp(argv[i],"++sig", 5) == 0){
      addSig = true;
    }
    if(strncmp(argv[i],"++data", 6) == 0){
      addData = true;
    }
    if(strncmp(argv[i],"++all", 5) == 0){
      addBkg  = true;
      addSig  = true;
      addData = true;
    }
    if(strncmp(argv[i],"+hist", 5) == 0){
      extrahist  = true;
      BuildFitInputCmd += "+hist ";
    }
    if(strncmp(argv[i],"+cat0L", 6) == 0){
      Categories += CT.GetCategories_0L();
      BuildFitInputCmd += "+cat0L ";
    }
    if(strncmp(argv[i],"+cat1L", 6) == 0){
      Categories += CT.GetCategories_1L();
      BuildFitInputCmd += "+cat1L ";
    }
    if(strncmp(argv[i],"+cat2L", 6) == 0){
      Categories += CT.GetCategories_2L();
      BuildFitInputCmd += "+cat2L ";
    }
    if(strncmp(argv[i],"+cat3L", 6) == 0){
      Categories += CT.GetCategories_3L();
      BuildFitInputCmd += "+cat3L ";
    }
    if(strncmp(argv[i],"++sys", 5) == 0){
      BuildFitInputCmd += "++sys ";
    }
    if(strncmp(argv[i],"++treesys",9) ==0){
      BuildFitInputCmd += "++treesys ";
    }
    if(strncmp(argv[i],"-maskSR", 7) == 0){
      BuildFitInputCmd += "-maskSR ";
    }
    if(strncmp(argv[i],"-lumi", 5) == 0){
      i++;
      setLumi = true;
      lumi = std::stof(argv[i]);
      BuildFitInputCmd += "-lumi " + string(argv[i]) + " ";
    }
    if(strncmp(argv[i],"-maxN", 5) == 0){
      i++;
      maxN = std::max(1, std::stoi(argv[i]));
    }
    if(strncmp(argv[i],"-sigfile", 8) == 0){
      i++;
      doSigFile = true;
      SigFile = argv[i];
    }
    if(strncmp(argv[i],"--dry-run", 9) == 0){
      dryRun = true;
    } 
    if(strncmp(argv[i],"--connect", 9) == 0){
      connect = true;
    }
    if(strncmp(argv[i],"-treeName", 9) == 0){
        i++;
	BuildFitInputCmd += "-treeName " + string(argv[i]) + " ";
    }
    if(strncmp(argv[i],"-treeSysName", 12) ==0){
        i++;
	BuildFitInputCmd += "-treeSysName " + string(argv[i]) + " ";
    }
    if(strncmp(argv[i],"+SMS", 4)==0){
        sigProc = true;
        BuildFitInputCmd += "+SMS ";
    }
   
 
  }
      
  if((proc_to_add.size() == 0) &&
     (addBkg  == false) &&
     (addSig  == false) &&
     (addData == false))
    bprint = true;

  if(bprint){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  Condor submission options:" << endl;
    cout << "   -maxN [number]      maximum number of processes per job" << endl;
    cout << "   -sigfile            signal filename must match this string to be included" << endl;
    cout << "   --connect           for running inside a batch job on CMS connect" << endl;
    cout << "   --dry-run           create output folders and scripts but don't submit" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)          print options" << endl;
    cout << "   -path [dest]        path to input ntuples" << endl;
    cout << "   --ouput(-o) [file]  output root file" << endl;
    cout << "   -year [year]        year to process" << endl;
    cout << "   +proc [label]       add processes matching label (can have >1 of these)" << endl;
    cout << "   ++bkg               add all background samples for year" << endl;
    cout << "   ++sig               add all signal samples" << endl;
    cout << "   ++data              add all background samples" << endl;
    cout << "   ++all               add all samples" << endl;
    cout << "   +cat0L              add 0L categories" << endl;
    cout << "   +cat1L              add 1L categories" << endl;
    cout << "   +cat2L              add 2L categories" << endl;
    cout << "   +cat3L              add 3L categories" << endl;
    cout << "   ++sys               turn on available systematics" << endl;
    cout << "   +hist               book 2D histograms also" << endl;
    cout << "   -lumi [lumi]        set luminosity to lumi" << endl;
    cout << "Example: ./BuildFitInputCondor.x -maxN 1 ++bkg +proc T2tt +cat1L -lumi 137 --connect -o name_of_BuildFitInput_output_folder/ -path root://xrootd.unl.edu//store/user/zflowers/crogan/ " << endl;

    return 0;
  }

  cout << "Initializing sample maps from path " << NtuplePath << " for year " << year << endl;
  SampleTool ST(NtuplePath, year);

  ProcessList samples;
  if(addBkg){
    cout << "Adding all background processes" << endl;
    samples += ST.Get(kBkg);
  }
  if(addSig){
    cout << "Adding all signal processes" << endl;
    samples += ST.Get(kSig);
  }
  if(addData){
    cout << "Adding all data for year " << year << endl;
    samples += ST.Get(kData);
  }
  for(int p = 0; p < int(proc_to_add.size()); p++){
    cout << "Adding processes that match \"" << proc_to_add[p] << "\"" << endl;
    //samples += ST.Get(proc_to_add[p]);
    if(sigProc){
     samples += ST.GetStrictSignalMatch(proc_to_add[p]);
    }else{
    //otherwise do the normal thing 
     samples += ST.Get(proc_to_add[p]);
    }
  }

  // prepare output folder
  string OutputFold = OutFile;
  if(OutputFold.find(".root") != string::npos){
    size_t found = OutputFold.find(".root");
    OutputFold = OutputFold.substr(0,found);
  }

  cout << "* Writing output to " << OutputFold << endl;
  gSystem->Exec(("rm -fr "+OutputFold).c_str());
  gSystem->Exec(("mkdir -p "+OutputFold).c_str());

  string SrcFold  = OutputFold+"/src/";
  string LogFold  = OutputFold+"/log/";
  string RootFold = OutputFold+"/root/";
  gSystem->Exec(("mkdir -p "+SrcFold).c_str());
  gSystem->Exec(("mkdir -p "+LogFold).c_str());
  gSystem->Exec(("mkdir -p "+RootFold).c_str());
  
  ofstream condorsubmit;
  condorsubmit.open((OutputFold+"/condor_submit.sh").c_str());
  
  int Nsample = samples.GetN();
  vector<string> procs;
  string iBFICmd;
  int Njob = 0;
  for(int p = 0; p < Nsample; p++){
    Process proc = samples[p];
//cout << "adding process " << proc.Name() << " to arguments" << endl;
    if(doSigFile && proc.Type() == kSig){
      bool keep = false;
      int Nfile = ST.NTrees(proc);
      for(int f = 0; f < Nfile; f++){
	string file = ST.FileName(proc, f);
	if(file.find(SigFile) != string::npos)
	  keep = true;
      }
      if(!keep)
	continue;
    }
    
    procs.push_back(proc.Name());
    iBFICmd += "+proc "+proc.Name()+" ";
    cout << "Adding proc " << proc.Name() << endl;
    if(procs.size() >= maxN || p == Nsample-1){
     if(connect)
     {
      int Nfile = ST.NTrees(proc);
      for(int f = 0; f < Nfile; f++){
        string new_BFICmd = BuildFitInputCmd;
        //if(proc.Type() == kBkg)
        //if(proc.Type() != kSig)
          new_BFICmd += ("-ifile "+std::to_string(f))+" ";
        if(connect)
          WriteScriptConnect(SrcFold+Form("submit_%d",Njob)+".sh",
          	  LogFold+Form("job_%d",Njob)+".log",
          	  new_BFICmd+iBFICmd+" -o "+Form("BFI_%d.root ", Njob),
          	  RootFold);
        condorsubmit << "condor_submit " << SrcFold << "submit_" << Njob << +".sh" << endl;

        Njob++;
       }
      }
      if(!connect)
      {
       WriteScript(SrcFold+Form("submit_%d",Njob)+".sh",
          	  LogFold+Form("job_%d",Njob)+".log",
          	  BuildFitInputCmd+iBFICmd+" -o "+RootFold+Form("BFI_%d.root ", Njob));
       Njob++;
      }
     iBFICmd = "";
     procs.clear();
    }
  }
  if(connect){
   MakeConnectTarball(RootFold);
  }
  condorsubmit.close();

  if(dryRun){
    cout << "Submit condor jobs by typing:" << endl;
    cout << "   source " << OutputFold << "/condor_submit.sh" << endl;
  } else {
    cout << "Submitting condor jobs" << endl;
    gSystem->Exec(("source "+OutputFold+"/condor_submit.sh").c_str());
  }
  
}
void MakeConnectTarball(const string& OutputFold){
  gSystem->Exec("mkdir -p config_BuildFitInput");
  gSystem->Exec("mkdir -p config_BuildFitInput/BtagSF");
  gSystem->Exec("cp BuildFitInput.x config_BuildFitInput/");
  gSystem->Exec("cp scripts/cmssw_setup_connect.sh config_BuildFitInput/");
  gSystem->Exec("cp scripts/setup_RestFrames_connect.sh config_BuildFitInput/");
  gSystem->Exec("cp root/BtagSF/*.root config_BuildFitInput/BtagSF/.");
  gSystem->Exec("cp csv/BtagSF/* config_BuildFitInput/BtagSF/.");
  gSystem->Exec("tar -czf config_BuildFitInput.tgz config_BuildFitInput/");
  gSystem->Exec(("mv config_BuildFitInput.tgz "+OutputFold+"/../").c_str());
  gSystem->Exec("rm -r config_BuildFitInput/");
  return;
}
void WriteScriptConnect(const string& src_name,
		 const string& log_name,
		 const string& command,
		 const string& OutputFold){
  ofstream file;
  file.open(src_name);
  string pwd = gSystem->pwd();
  string root_output = command.substr(command.find("-o")+2);
  if(root_output.find("/") != string::npos){
   int pos = root_output.find_last_of("/")+1;
   root_output = root_output.substr(pos, root_output.length()-pos);
  }
//cout << "command: " << command << endl;
/* outsource testing to speed up submission
  gSystem->Exec("mkdir -p config_BuildFitInput");
  gSystem->Exec("mkdir -p config_BuildFitInput/BtagSF");
  gSystem->Exec("cp BuildFitInput.x config_BuildFitInput/");
  gSystem->Exec("cp scripts/cmssw_setup_connect.sh config_BuildFitInput/");
  gSystem->Exec("cp scripts/setup_RestFrames_connect.sh config_BuildFitInput/");
  gSystem->Exec("cp root/BtagSF/*.root config_BuildFitInput/BtagSF/.");
  gSystem->Exec("cp csv/BtagSF/* config_BuildFitInput/BtagSF/.");
  gSystem->Exec("tar -czf config_BuildFitInput.tgz config_BuildFitInput/");
  gSystem->Exec(("mv config_BuildFitInput.tgz "+OutputFold+"/../").c_str());
  gSystem->Exec("rm -r config_BuildFitInput/");
*/ 
  file << "universe = vanilla" << endl;
  file << "executable = execute_script_BuildFitInput.sh" << endl;
  file << "getenv = True" << endl;
  file << "use_x509userproxy = true" << endl;
  file << "Arguments = " << command << endl;
  file << "output = " << log_name << ".out" << endl;
  file << "error = "  << log_name << ".err" << endl;
  file << "log = "    << log_name << ".log" << endl;
  file << "Requirements = (Machine != \"red-node000.unl.edu\") && (Machine != \"red-c2325.unl.edu\")" << endl;
  file << "request_memory = 2 GB" << endl;
  file << "transfer_input_files = /uscms/home/janguian/nobackup/CMSSW_10_6_5/src/KUEWKinoAnalysis_NewNtuples/csv/METTrigger/Parameters.csv,/uscms/home/z374f439/nobackup/whatever_you_want/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2,"+pwd+"/"+OutputFold+"/../config_BuildFitInput.tgz" << endl;
 // file << "transfer_input_files = /home/jsingera/jsingera/CMSSW_10_6_5/src/KUEWKinoAnalysis_treeSysDev/csv/METTrigger/Parameters.csv,/ospool/cms-user/zflowers/public/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2,"+pwd+"/"+OutputFold+"/../config_BuildFitInput.tgz" << endl;

  file << "should_transfer_files = YES" << endl;
  file << "when_to_transfer_output = ON_EXIT" << endl;
  file << "transfer_output_files = "+root_output << endl;
  file << "transfer_output_remaps = \""+root_output+"="+pwd+"/"+OutputFold+"/"+command.substr(command.find("-o ")+3)+"\"" << endl;
  file << "+ProjectName=\"cms.org.cern\""<< endl;
  file << "+REQUIRED_OS=\"rhel7\"" << endl;
  file << "+RequiresCVMFS=True" << endl;
  file << "queue " << endl;
  file.close();  
}

void WriteScript(const string& src_name,
		 const string& log_name,
		 const string& command){
  ofstream file;
  file.open(src_name);
  
  file << "universe = vanilla" << endl;
  file << "executable = BuildFitInput.x" << endl;
  file << "getenv = True" << endl;
  file << "use_x509userproxy = true" << endl;
  file << "Arguments = " << command << endl;
  file << "output = " << log_name << ".out" << endl;
  file << "error = "  << log_name << ".err" << endl;
  file << "log = "    << log_name << ".log" << endl;
  file << "Requirements = (Machine != \"red-node000.unl.edu\") && (Machine != \"red-c2325.unl.edu\")" << endl;
  file << "request_memory = 2 GB" << endl;
  file << "queue " << endl;
  file.close();  
}
