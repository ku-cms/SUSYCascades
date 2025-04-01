#include <iostream>
#include <string>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TColorWheel.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>

#include "NeventTool.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;
string path_to_lists = "/home/zflowers/CMSSW_13_3_1/src/SUSYCascades/samples/NANO/";

void WriteScript(const string& src_dir, const string& src_name, const string& filetag, const string& dataset, const string& proc_Name, const long double& weight_val, const int& split, const string& era, const string& extraargs);
double getTot(string input_dataset, string input_filetag, string EventCountFile);
int main(int argc, char* argv[]) {
  bool bprint = false;
  bool debug = false;
  int split = 10;
  string otag = "";
  string extraargs = "";
  bool dryrun = false;
  for(int i = 0; i < argc; i++){
    if(strncmp(argv[i],"--help", 6) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-h", 2) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"--debug", 7) == 0){
      debug = true;
    }
    if(strncmp(argv[i],"--otag", 6) == 0){
      i++;
      otag = "_"+string(argv[i]);
    }
    if(strncmp(argv[i],"--split",7) == 0){
      i++;
      sscanf(argv[i],"%d", &split);
    }
    if(strncmp(argv[i],"--dry-run",9) == 0){
      dryrun = true;
    }
    if(strncmp(argv[i],"--extraargs",11) == 0){
      for(int j = i+1; j < argc; j++){
        extraargs += string(argv[j])+" ";
      }
    }
  }

  if(bprint){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)  diplay this help menu" << endl;
    cout << "   --debug     turn on debug mode" << endl;
    cout << "   --otag      add extra tag to output dir" << endl;
    cout << "   --split     amount to split jobs by (higher = more jobs)" << endl;
    cout << "   --dry-run   do NOT submit to condor" << endl;        
    cout << "   --extraargs extra args to pass to Condor_Plot_1D_NANO.x" << endl; 
    cout << "Example: ./Submit_Plot_1D_NANO.x --otag BDT_MET --split 5 --extraargs BDT_MET" << endl;
    return 0;
  }

  RestFrames::SetStyle();

   //std::map<string,long double> dataset_list_2017 = {
   // {"/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8.txt", (1.0/getTot("WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8","Fall17_102X","root/EventCount/EventCount_NANO_Fall17_102X.root")) * 12.87 * 1.21},
   //};

   std::map<string,long double> WZ_dataset_list_2023BPix = {
    {"/WZ_TuneCP5_13p6TeV_pythia8.txt",1.},
   };
   std::map<string,long double> TTJets_dataset_list_2023BPix = {
    {"/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8.txt",1.},
    {"/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8.txt",1.},
   };
   std::map<string,long double> WJets_dataset_list_2017 = {
    {"/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
   };
   std::map<string,long double> ZDY_dataset_list_2017 = {
    {"/DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
   };
   std::map<string,long double> TTJets_dataset_list_2017 = {
    {"/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
   };
   std::map<string,long double> TTZToQQ_dataset_list_2017 = {
    {"/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8.txt",1.},
   };
   std::map<string,long double> WZ_dataset_list_2017 = {
    {"/WZ_TuneCP5_13TeV-pythia8.txt",1.},
   };
   std::map<string,long double> WW_dataset_list_2017 = {
    {"/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8.txt",1.},
    //{"/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8.txt",1.},
   };
   std::map<string,long double> WWZ_dataset_list_2017 = {
    {"/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",1.},
   };
  std::map<string,long double> TChiWZ_3_50_dataset_list_2017 = {
    {"/TChiWZ_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,long double> TChiWZ_60_90_dataset_list_2017 = {
    {"/SMS-TChiWZ_dM-60to90_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,long double> Cascades_20_dataset_list_2023BPix = {
    {"/SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,long double> Cascades_10_dataset_list_2023BPix = {
    {"/SlepSnuCascade_MN1-260_MN2-280_MC1-270_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,long double> Cascades_5_dataset_list_2023BPix = {
    {"/SlepSnuCascade_MN1-270_MN2-280_MC1-275_TuneCP5_13p6TeV_madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,long double> Cascades_10_Justin_dataset_list_2023BPix = {
    {"/SlepSnuCascade_MN1-260_MN2-280_MC1-270_TuneCP5_13p6TeV_madgraphMLM-pythia8_Justin.txt",1.},
  };

  std::map<string,std::map<string,long double>> datasets_2017 = {
    {"ttbar",TTJets_dataset_list_2017},
    //{"TTZToQQ",TTZToQQ_dataset_list_2017},
    //{"DB_WZ",WZ_dataset_list_2017},
    {"DB_WW",WW_dataset_list_2017},
    {"ZDY",ZDY_dataset_list_2017},
    {"WJets",WJets_dataset_list_2017},
    //{"TB_WWZ",WWZ_dataset_list_2017},
  };

  std::map<string,std::map<string,long double>> datasets_2017_SMS = {
    //{"TChiWZ_20",TChiWZ_3_50_dataset_list_2017},
    //{"TChiWZ_50",TChiWZ_3_50_dataset_list_2017},
    //{"TChiWZ_90",TChiWZ_60_90_dataset_list_2017},
  };
  std::map<string,std::map<string,long double>> datasets_2023_BPix_SMS = {
    {"Cascades_20",Cascades_20_dataset_list_2023BPix},
    {"Cascades_10",Cascades_10_dataset_list_2023BPix},
    //{"Cascades_10_Justin",Cascades_10_Justin_dataset_list_2023BPix},
    {"Cascades_5",Cascades_5_dataset_list_2023BPix},
  };
  std::map<string,std::map<string,long double>> datasets_2023_BPix = {
    {"ttbar",TTJets_dataset_list_2023BPix},
    //{"DB_WZ",WZ_dataset_list_2023BPix},
  };

  vector<std::map<string,std::map<string,long double>>> datasets = {
    datasets_2017,
    datasets_2017_SMS,
    datasets_2023_BPix_SMS,
    //datasets_2023_BPix,
  };

  gSystem->Exec("make cmssw -j4");
  string work_dir = "SubmitCondor_Plot_1D_NANO";
  work_dir += otag + "/"; // need '/' at end
  gSystem->Exec(("mkdir -p "+work_dir+"src/").c_str());
  gSystem->Exec(("mkdir -p "+work_dir+"root/").c_str());
  gSystem->Exec(("mkdir -p "+work_dir+"log/").c_str());
  gSystem->Exec(("mkdir -p "+work_dir+"out/").c_str());
  gSystem->Exec(("mkdir -p "+work_dir+"err/").c_str());
  gSystem->Exec(("cp Condor_Plot_1D_NANO.x "+work_dir).c_str());
  gSystem->Exec(("cp src/Condor_Plot_1D_NANO.C "+work_dir).c_str());

  string filetag = "";
  string era = "Run2";
  for(int d = 0; d < int(datasets.size()); d++){
   if(d == 0) 
   {
    cout << "   Submitting 2017..." << endl;
    filetag = "Fall17_102X";
   } 
   if(d == 1)
   {
    cout << "   Submitting 2017 SMS..." << endl;
    filetag = "Fall17_102X_SMS";
   } 
   if(d == 2)
   {
    cout << "   Submitting 2023BPix SMS..." << endl;
    filetag = "Summer23BPix_130X_SMS";
    era = "Run3";
   }
   if(d == 3)
   {
    cout << "   Submitting 2023BPix..." << endl;
    filetag = "Summer23BPix_130X";
    era = "Run3";
   }
   //setup and submit condor job
   for(std::map<string,std::map<string,long double>>::iterator iter = datasets[d].begin(); iter != datasets[d].end(); ++iter){
     string proc_Name = iter->first;
     std::map<string,long double> dataset_list = iter->second;
     for(std::map<string,long double>::iterator iter2 = dataset_list.begin(); iter2 != dataset_list.end(); ++iter2){
        string dataset = iter2->first;
        long double weight = iter2->second;
        string dataset_name = get_str_between_two_str(dataset,"/",".");
        string src_name = work_dir+"src/"+filetag+"_"+dataset_name+"_"+proc_Name+".sh";
        WriteScript(work_dir, src_name, filetag, dataset, proc_Name, weight, split, era, extraargs); 
        if(!dryrun)
          gSystem->Exec(("condor_submit "+src_name).c_str());
      } // for(std::map<string,long double>::iterator iter = datasets[d].begin(); iter != datasets[d].end(); ++iter)
    } // for(std::map<string,std::map<string,long double>>::iterator iter = datasets[d].begin(); iter != datasets[d].end(); ++iter)
  } // for(int d = 0; d < int(datasets.size()); d++)
}

void WriteScript(const string& work_dir, const string& src_name, const string& filetag, const string& dataset, const string& proc_Name, const long double& weight_val, const int& split, const string& era, const string& extraargs)
{
  ofstream file;
  file.open(src_name);
  string pwd = gSystem->pwd();
  string dataset_name = get_str_between_two_str(dataset,"/",".");
  file << "universe = vanilla" << endl;
  file << "executable = execute_script_Plot_1D_NANO.sh" << endl;
  file << "getenv = True" << endl;
  file << "use_x509userproxy = true" << endl;
  file << "Requirements = HAS_SINGULARITY == True" << endl;
  file << "MY.SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9\"" << endl;
  string args = "--proc "+proc_Name+" --weight "+std::to_string(weight_val*1.e7)+" --ifile $(Item) --ofile "+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)_$(Step).root --split=$$([$(Step)+1]),"+std::to_string(split)+" --era "+era+" "+extraargs;
  file << "Arguments = " << args << endl;
  file << "output = " << work_dir+"out/"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)" << ".out" << endl;
  file << "error = "  << work_dir+"err/"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)" << ".err" << endl;
  file << "log = "    << work_dir+"log/"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)" << ".log" << endl;
  file << "transfer_input_files = "+work_dir+"/Condor_Plot_1D_NANO.x,/ospool/cms-user/zflowers/public/sandbox-CMSSW_13_3_1-el9.tar.bz2,scripts/cmssw_setup_connect_el9.sh,scripts/setup_RestFrames_connect.sh,"+path_to_lists+filetag+dataset << endl;
  file << "should_transfer_files = YES" << endl;
  file << "when_to_transfer_output = ON_EXIT" << endl;
  file << "transfer_output_files = output_Plot_1D_NANO_"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)"+".root" << endl;
  file << "transfer_output_remaps = \"output_Plot_1D_NANO_"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)"+".root="+pwd+"/"+work_dir+"/root/output_Plot_1D_NANO_"+filetag+"_"+dataset_name+"_"+proc_Name+"_$(ItemIndex)"+"_$(Step)"+".root\"" << endl;
  file << "+ProjectName=\"cms.org.ku\"" << endl;
  file << "+REQUIRED_OS=\"rhel9\"" << endl;
  file << "job_lease_duration = 3600" << endl;
  file << "RequestCpus=ifthenelse((isUndefined(CpusUsage) || CpusUsage < 2),1,MIN(RequestCpus+2, 32))" << endl;
  file << "periodic_hold = (CpusUsage >= RequestCpus) && (JobStatus == 3)" << endl;
  file << "periodic_hold_subcode = 42" << endl;
  file << "periodic_release = (HoldReasonCode == 12 && HoldReasonSubCode == 256 || HoldReasonCode == 13 && HoldReasonSubCode == 2 || HoldReasonCode == 12 && HoldReasonSubCode == 2 || HoldReasonCode == 26 && HoldReasonSubCode == 120 || HoldReasonCode == 3 && HoldReasonSubCode == 0 || HoldReasonSubCode == 42)" << endl;
  file << "+RequiresCVMFS = True" << endl;
  file << "+RequiresSharedFS = True" << endl;
  //file << "priority = 10" << endl;
  file << "queue "+std::to_string(split)+" from "+path_to_lists+filetag+dataset << endl;
  file.close();  
}

double getTot(string input_dataset, string input_filetag, string EventCountFile)
{
  TFile* fout = new TFile(EventCountFile.c_str(),"READ");
  TBranch* b_dataset = nullptr;
  TBranch* b_filetag = nullptr;
  TBranch* b_Nweight = nullptr;
  TBranch* b_Nevent = nullptr;
  TBranch* b_MP = nullptr;
  TBranch* b_MC = nullptr;
  string* dataset = nullptr;
  string* filetag = nullptr;
  double Nweight = 0.0;
  double Nevent = 0.0;
  int MP = 0;
  int MC = 0;
  TTree* tree = nullptr;
  tree = (TTree*) fout->Get("EventCount");
  
  tree->SetMakeClass(1);
  tree->SetBranchAddress("Nevent", &Nevent,&b_Nevent);
  tree->SetBranchAddress("Nweight", &Nweight,&b_Nweight);
  tree->SetBranchAddress("filetag", &filetag,&b_filetag);
  tree->SetBranchAddress("dataset", &dataset,&b_dataset);
  tree->SetBranchAddress("MP", &MP,&b_MP);
  tree->SetBranchAddress("MC", &MC,&b_MC);

  double tot_Nevent = 0;
  double tot_Nweight = 0;

  for(int i = 0; i < tree->GetEntries(); i++)
  {
   tree->GetEntry(i);
   if(dataset->find(input_dataset) != std::string::npos && filetag->find(input_filetag) != std::string::npos)
   {
    tot_Nevent += Nevent;
    tot_Nweight += Nweight;
   }
  }
  fout->Close();
  delete fout;
  return tot_Nweight;
}
