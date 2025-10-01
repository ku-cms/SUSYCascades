// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <dirent.h>
#include <vector>
#include <variant>
#include <memory>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>

#include "ReducedNtuple.hh"
#include "SUSYNANOBase.hh"
#include "NANOULBase.hh"
#include "NANORun3.hh"

using namespace std;
using std::vector;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[400];
  char inputListName[400];
  char inputFolderName[400];
  char outputFileName[400];
  char TreeName[400];
  char DataSet[400];
  char FileTag[400];
  char EventCount[400];
  char FilterEff[400];
  char JSONFile[400];
  char PUFOLD[400];
  char BTAGFOLD[400];
  char LEPFOLD[400];
  char JMEFOLD[400];
  char METTRIGFILE[400];
  char PREFIREFILE[400];
  char XSJSONFILE[400];

  bool DO_FILE = false;
  bool DO_LIST = false;
  bool DO_FOLDER = false;
  bool DO_TREE = false;
  bool DO_SMS = false;
  bool DO_CASCADES = false;
  bool DO_PRIVATEMC = false;
  bool DO_JSON = false;
  bool IS_DATA = false;
  bool IS_FASTSIM = false;

  bool DO_SYS = false;
  bool DO_SYS_JES = false;
  bool DO_SYS_JER = false;
  bool DO_SYS_MET = false;
  bool DO_SYS_MMS = false;
  bool DO_SYS_EES = false;

  bool DO_slim = false;
  
  int ICHUNK = 1;
  int NCHUNK = 1;

  if ( argc < 2 ){
    cout << "Error at Input: please specify an input file name, a list of input ROOT files and/or a folder path"; 
    cout << " and an output filename:" << endl; 
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ifile=input.root -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ilist=input.list -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ifold=folder_path -ofile=output.root -dataset=dataset_name -filetag=sample_tag -tree=treename -eventcount=event_count --sms" << endl;
    
    return 1;
  }
  for (int i=0;i<argc;i++){
    if (strncmp(argv[i],"-ifile",6)==0){
      sscanf(argv[i],"-ifile=%s",  inputFileName);
      DO_FILE = true;
    }
    if (strncmp(argv[i],"-ilist",6)==0){
      sscanf(argv[i],"-ilist=%s",  inputListName);
      DO_LIST = true;
    }
    if (strncmp(argv[i],"-ifold",6)==0){
      sscanf(argv[i],"-ifold=%s",  inputFolderName);
      DO_FOLDER = true;
    }
    if (strncmp(argv[i],"-tree",5)==0){
      sscanf(argv[i],"-tree=%s",  TreeName);
      DO_TREE = true;
    }
    if (strncmp(argv[i],"-json",5)==0){
      sscanf(argv[i],"-json=%s",  JSONFile);
      DO_JSON = true;
    }
    if (strncmp(argv[i],"-ofile",6)==0) sscanf(argv[i],"-ofile=%s", outputFileName);
    
    if (strncmp(argv[i],"-dataset",8)==0)   sscanf(argv[i],"-dataset=%s", DataSet);
    if (strncmp(argv[i],"-filetag",8)==0)   sscanf(argv[i],"-filetag=%s", FileTag);
    if (strncmp(argv[i],"-eventcount",11)==0)   sscanf(argv[i],"-eventcount=%s", EventCount);
    if (strncmp(argv[i],"-filtereff",10)==0)   sscanf(argv[i],"-filtereff=%s", FilterEff);
    if (strncmp(argv[i],"-pu",3)==0)   sscanf(argv[i],"-pu=%s", PUFOLD);
    if (strncmp(argv[i],"-btag",5)==0)   sscanf(argv[i],"-btag=%s", BTAGFOLD);
    if (strncmp(argv[i],"-lep",4)==0)   sscanf(argv[i],"-lep=%s", LEPFOLD);
    if (strncmp(argv[i],"-jme",4)==0)   sscanf(argv[i],"-jme=%s", JMEFOLD);
    if (strncmp(argv[i],"-metfile",8)==0)   sscanf(argv[i],"-metfile=%s", METTRIGFILE);
    if (strncmp(argv[i],"-prefirefile",12)==0)   sscanf(argv[i],"-prefirefile=%s", PREFIREFILE);
    if (strncmp(argv[i],"-xsjsonfile",11)==0)   sscanf(argv[i],"-xsjsonfile=%s", XSJSONFILE);
    
    if (strncmp(argv[i],"--sms",5)==0)  DO_SMS = true;
    if (strncmp(argv[i],"--data",6)==0)  IS_DATA = true;
    if (strncmp(argv[i],"--fastsim",9)==0)  IS_FASTSIM = true;
    if (strncmp(argv[i],"--cascades",10)==0)  DO_CASCADES = true;
    if (strncmp(argv[i],"--privatemc",11)==0)  DO_PRIVATEMC = true;
    if (strncmp(argv[i],"--private",9)==0)  DO_PRIVATEMC = true;

    if (strncmp(argv[i],"--slim",6)==0) DO_slim = true;
    
    if (strncmp(argv[i],"--sys",5)==0)  DO_SYS = true;
    if (strncmp(argv[i],"--sysJES",8)==0)  DO_SYS_JES = true;
    if (strncmp(argv[i],"--sysJER",8)==0)  DO_SYS_JER = true;
    if (strncmp(argv[i],"--sysMET",8)==0)  DO_SYS_MET = true;
    if (strncmp(argv[i],"--sysMMS",8)==0)  DO_SYS_MMS = true;
    if (strncmp(argv[i],"--sysEES",8)==0)  DO_SYS_EES = true;

    if (strncmp(argv[i],"-split",6)==0)  sscanf(argv[i],"-split=%d,%d", &ICHUNK, &NCHUNK);
  }

  gROOT->ProcessLine("#include <vector>");

  vector<string> filenames;

  char Buffer[500];
  char MyRootFile[2000];  

  if(DO_FOLDER){
    DIR *dpdf;
    struct dirent *epdf;
    dpdf = opendir(inputFolderName);
    if(dpdf == NULL){
      cout << "ERROR: " << inputFolderName << " is not a directory" << endl;
      return 1;
    }
    string dname(inputFolderName);
    while ((epdf = readdir(dpdf))){
      if(string(epdf->d_name).find(".root") == string::npos)
	continue;
      string name = dname+"/"+string(epdf->d_name);
      filenames.push_back(name);
    }
  }

  if(DO_LIST){
    ifstream *inputFile = new ifstream(inputListName);
    while( !(inputFile->eof()) ){
      inputFile->getline(Buffer,500);
      if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))){
	sscanf(Buffer,"%s",MyRootFile);
	filenames.push_back(MyRootFile);
      }
    }
    inputFile->close();
    delete inputFile;
  }

  if(DO_FILE){
    filenames.push_back(inputFileName);
  }

  TChain* chain;
  if(DO_TREE)
    chain = (TChain*) new TChain(TreeName);
  else
    chain = (TChain*) new TChain("Events");
  
  // add DAS count
  int NDAS = 0;
  std::string DAS_datasetname = "PrivateMC";
  std::string DAS_filename = "PrivateMC";
  NeventTool eventTool;
  int Nfile = filenames.size();
  for(int i = 0; i < Nfile; i++){
    chain->Add(filenames[i].c_str());
    if(!DO_PRIVATEMC) NDAS += eventTool.EventsInDAS(filenames[i]);
    cout << "   Added file " << filenames[i] << endl;
  }
  if(!DO_PRIVATEMC) {
    if(NDAS == 0) return 1; // will try to resubmit job
    DAS_datasetname = eventTool.Get_DASdatasetname(filenames[0]);
    DAS_filename = eventTool.Get_DASfilename(filenames[0]);
  }

  std::cout << "Testing if input files are accessible" << std::endl;
  
  // configure redirectors to try: first entry = "" means "use filenames as-is (no replacement)"
  // subsequent entries are replacements for the known srcPrefix
  const std::string srcPrefix = "root://cmsxrootd.fnal.gov/";
  std::vector<std::string> redirectors = {
      "", // attempt 0: original redirector
      "davs://xrootd-local.unl.edu:1094/", // T2_US_Nebraska
      "root://cmsdcadisk.fnal.gov:1094//dcache/uscmsdisk/", // T1_US_FNAL_DISK
      "davs://k8s-redir-stageout.ultralight.org:1094/", // T2_US_Caltech
  };
  
  bool success = false;
  Long64_t NTOT = -1;
  int timeoutSeconds = 60;
  
  // save previous ROOT error level and then set to kFatal during GetEntries
  int prevErrLevel = gErrorIgnoreLevel;
  
  for (size_t attempt = 0; attempt < redirectors.size(); ++attempt) {
      if (attempt == 0) {
          // reuse the chain that was built during the DAS loop above
      } else {
          std::cerr << "[attempt " << attempt << "] GetEntries timed out previously; retrying with redirector: "
                    << redirectors[attempt] << std::endl;
          // delete old chain and rebuild with the new redirector
          if (chain) {
              delete chain;
              chain = nullptr;
          }
          chain = buildChain(filenames, srcPrefix, redirectors[attempt], TreeName, DO_TREE);
      }
  
      // try GetEntries with timeout
      NTOT = tryGetEntriesSubproc(chain, timeoutSeconds);
  
      if (NTOT >= 0) {
          success = true;
          if (attempt != 0)
              std::cout << "Got entries with prefix: " << redirectors[attempt] << std::endl;
          else
              std::cout << "Got entries with prefix: " << srcPrefix << std::endl;
          break;
      } else {
          std::cerr << "GetEntries() attempt " << attempt << " timed out after " << timeoutSeconds << " seconds." << std::endl;
          // loop continues to next redirector
      }
  }
  
  if (!success) {
      std::cerr << "All redirector attempts failed. Exiting with failure (return 1)." << std::endl;
      // cleanup chain pointer
      if (chain) { delete chain; chain = nullptr; }
      return 1;
  }

  using NtupleVariant = std::variant<std::unique_ptr<ReducedNtuple<SUSYNANOBase>>, std::unique_ptr<ReducedNtuple<NANOULBase>>, std::unique_ptr<ReducedNtuple<NANORun3>>>;
  NtupleVariant ntuple;
  if(string(FileTag).find("130X") != std::string::npos){ cout << "Using Run3 base" << endl; ntuple = std::make_unique<ReducedNtuple<NANORun3>>(chain); }
  else if(string(FileTag).find("UL") != std::string::npos){ cout << "Using UL base" << endl; ntuple = std::make_unique<ReducedNtuple<NANOULBase>>(chain); }
  else { cout << "Using Run2 base" << endl; ntuple = std::make_unique<ReducedNtuple<SUSYNANOBase>>(chain); }

  Long64_t N1, N0;
  std::visit([&](auto& nt) { nt->GetChunks(NDAS, N1, N0, ICHUNK, NCHUNK); }, ntuple);
  NDAS = N0 - N1;

  std::visit([&](auto& nt) { nt->AddLabels(string(DataSet),string(FileTag)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddEventCountFile(string(EventCount)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddFilterEffFile(string(FilterEff)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddPUFolder(string(PUFOLD)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddBtagFolder(string(BTAGFOLD)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddLepFolder(string(LEPFOLD)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddJMEFolder(string(JMEFOLD)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddMETTriggerFile(string(METTRIGFILE)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddPrefireFile(string(PREFIREFILE)); }, ntuple);
  std::visit([&](auto& nt) { nt->AddXSecJSON(string(XSJSONFILE)); }, ntuple);
  #ifdef _CMSSW_
  if(!DO_SMS && !IS_DATA && !DO_CASCADES)
    //ntuple->AddLHAPDF();
    std::visit([](auto& nt) { nt->AddLHAPDF(); }, ntuple);
  #endif

  if(DO_SYS)
    std::visit([](auto& nt) { nt->AddSystematics(); }, ntuple);
  if(DO_SYS_JES)
    std::visit([](auto& nt) { nt->AddJESSystematics(); }, ntuple);
  if(DO_SYS_JER)
    std::visit([](auto& nt) { nt->AddJERSystematics(); }, ntuple);
  if(DO_SYS_MET)
    std::visit([](auto& nt) { nt->AddMETSystematics(); }, ntuple);
  if(DO_SYS_EES)
    std::visit([](auto& nt) { nt->AddEESSystematics(); }, ntuple);
  if(DO_SYS_MMS)
    std::visit([](auto& nt) { nt->AddMMSSystematics(); }, ntuple);
  
  if(DO_JSON)
    std::visit([&](auto& nt) { nt->AddJSONFile(string(JSONFile)); }, ntuple);

  if(DO_SMS)
    std::visit([](auto& nt) { nt->DoSMS(); }, ntuple);

  if(IS_DATA)
    std::visit([](auto& nt) { nt->DoData(); }, ntuple);

  if(IS_FASTSIM)
    std::visit([](auto& nt) { nt->DoFastSim(); }, ntuple);

  if(DO_CASCADES)
    std::visit([](auto& nt) { nt->DoCascades(); }, ntuple);

  if(DO_PRIVATEMC)
    std::visit([](auto& nt) { nt->DoPrivateMC(); }, ntuple);

  cout << "writing output with ichunk=" << ICHUNK << " nchunk=" << NCHUNK << endl;
  bool passedDASCheck = std::visit([&](auto& nt) -> bool { return nt->WriteNtuple(string(outputFileName), ICHUNK, NCHUNK, DO_slim, NDAS, string(DAS_datasetname), string(DAS_filename)); }, ntuple);
  if(!passedDASCheck){
    std::cout << "JOB FAILED DAS CHECK!" << std::endl;
    return 1;
  }
  else return 0;

}
