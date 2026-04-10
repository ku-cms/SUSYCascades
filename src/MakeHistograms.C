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

  bool DO_FILE = false;
  bool DO_LIST = false;
  bool DO_FOLDER = false;
  bool DO_TREE = false;
  bool DO_SMS = false;
  bool DO_CASCADES = false;
  bool DO_PRIVATEMC = false;

  int ICHUNK = 1;
  int NCHUNK = 1;

  if ( argc < 2 ){
    cout << "Error at Input: please specify an input file name, a list of input ROOT files and/or a folder path"; 
    cout << " and an output filename:" << endl; 
    cout << "  Example:      ./MakeHistograms.x -ifile=input.root -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    
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
    if (strncmp(argv[i],"-ofile",6)==0) sscanf(argv[i],"-ofile=%s", outputFileName);
    
    if (strncmp(argv[i],"-dataset",8)==0)   sscanf(argv[i],"-dataset=%s", DataSet);
    if (strncmp(argv[i],"-filetag",8)==0)   sscanf(argv[i],"-filetag=%s", FileTag);
    
    if (strncmp(argv[i],"--sms",5)==0)  DO_SMS = true;
    if (strncmp(argv[i],"--cascades",10)==0)  DO_CASCADES = true;
    if (strncmp(argv[i],"--privatemc",11)==0)  DO_PRIVATEMC = true;
    if (strncmp(argv[i],"--private",9)==0)  DO_PRIVATEMC = true;

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
      "root://cmsxrootd.fnal.gov/",
      "davs://xrootd-local.unl.edu:1094/", // T2_US_Nebraska
      "root://cmsdcadisk.fnal.gov:1094//dcache/uscmsdisk/", // T1_US_FNAL_DISK
      "davs://k8s-redir-stageout.ultralight.org:1094/", // T2_US_Caltech
      "davs://cmsxrootd.hep.wisc.edu:1094", // T2_US_Wisconsin
      "davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms", // T2_UK_London_IC
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

  if(DO_SMS) std::visit([](auto& nt) { nt->DoSMS(); }, ntuple);
  if(DO_CASCADES) std::visit([](auto& nt) { nt->DoCascades(); }, ntuple);
  if(DO_PRIVATEMC) std::visit([](auto& nt) { nt->DoPrivateMC(); }, ntuple);

  std::visit([&](auto& nt) { nt->AddLabels(string(DataSet),string(FileTag),string(DAS_filename)); }, ntuple);
  std::visit([&](auto& nt) { nt->SetupBtagWP(); }, ntuple);

  cout << "writing output with ichunk=" << ICHUNK << " nchunk=" << NCHUNK << endl;
  bool passedDASCheck = std::visit([&](auto& nt) -> bool { return nt->WriteNtuple(string(outputFileName), ICHUNK, NCHUNK, false, NDAS, string(DAS_datasetname), string(DAS_filename), true, false); }, ntuple);
  if(!passedDASCheck){
    std::cout << "JOB FAILED DAS CHECK!" << std::endl;
    return 1;
  }
  else return 0;
}
