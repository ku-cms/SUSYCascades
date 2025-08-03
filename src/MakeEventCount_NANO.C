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
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>

#include "NeventTool.hh"

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
  char gitHash[400];

  bool DO_FILE = false;
  bool DO_LIST = false;
  bool DO_FOLDER = false;
  bool DO_TREE = false;
  bool DO_SMS = false;
  bool DO_CASCADES = false;
  bool DO_PRIVATEMC = false; // Private production (DAS info does not exist)

  if ( argc < 2 ){
    cout << "Error at Input: please specify an input file name, a list of input ROOT files and/or a folder path"; 
    cout << " and an output filename:" << endl; 
    cout << "  Example:      ./MakeEventCount_NANO.x -ifile=input.root -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeEventCount_NANO.x -ilist=input.list -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeEventCount_NANO.x -ifold=folder_path -ofile=output.root -dataset=dataset_name -filetag=sample_tag -tree=treename --sms" << endl;
    
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
    if (strncmp(argv[i],"-githash",8)==0)  sscanf(argv[i],"-githash=%s", gitHash);
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

  Float_t genWeight;
  UInt_t  luminosityBlock;
  
  TBranch *b_genWeight;
  TBranch *b_luminosityBlock;

  Int_t   nGenPart_Run3;
  UInt_t  nGenPart_Run2;
  Float_t GenPart_mass[400];  
  Int_t   GenPart_pdgId[400];
  TBranch *b_nGenPart;
  TBranch *b_GenPart_mass;
  TBranch *b_GenPart_pdgId;

  double Nevent = 0.;
  double Nweight = 0.;

  Double_t genEventSumw;
  Long64_t genEventCount;
  TBranch  *b_genEventSumw;
  TBranch  *b_genEventCount;

  string dataset = string(DataSet);
  string filetag = string(FileTag);

  TChain* chain;
  if(DO_TREE)
    chain = (TChain*) new TChain(TreeName);
  
  int Nfile = filenames.size();
  for(int i = 0; i < Nfile; i++){
    chain->Add(filenames[i].c_str());
    cout << "   Adding file " << filenames[i] << endl;
  }
    
  cout << "Setting Branch Addresses" << endl;

  bool status1 = false;
  bool status2 = false;
  if(strcmp(TreeName,"Runs") == 0){
    if(chain->GetListOfBranches()->FindObject("genEventSumw"))
      status1 = true;
    if(chain->GetListOfBranches()->FindObject("genEventCount"))
      status2 = true;
  }
  int NEVENT = 0;

  if(status1 && status2){
    chain->SetBranchAddress("genEventSumw", &genEventSumw, &b_genEventSumw);
    chain->SetBranchAddress("genEventCount", &genEventCount, &b_genEventCount);
  }

  // add DAS count
  int NDAS = 0;
  int Nevent_tot = 0;
  NeventTool eventTool;
  std::string DAS_datasetname = "PrivateMC";
  std::string DAS_filename = "PrivateMC";
  if(!DO_PRIVATEMC){
    cout << "Adding DAS info..." << endl;
    for(int i = 0; i < Nfile; i++)
      NDAS += eventTool.EventsInDAS(filenames[i]);
    DAS_datasetname = eventTool.Get_DASdatasetname(filenames[0]);
    DAS_filename = eventTool.Get_DASfilename(filenames[0]);
    if(NDAS == 0) return 1; // will try to resubmit job
    cout << "Added DAS info!" << endl;
  }
  bool passed_DAS = true;
  if(status1 && status2){
    chain->GetEntry(0);
    Nevent = genEventCount;
    Nweight = genEventSumw;
    if(!DO_PRIVATEMC && NDAS != Nevent){
      passed_DAS = false;
      Nevent = 0;
      Nweight = 0;
      NEVENT = 0;
    }
  }

  int MP = 0;
  int MC = 0;
  int Code = 0;
  int PDGID;
  std::vector<std::pair<int,std::pair<int,int>>> masses;
  std::map<std::pair<int,std::pair<int,int>>,double > mapNevent;
  std::map<std::pair<int,std::pair<int,int>>,double > mapNweight;
  int maxNGEN = 0;

  if(!status1 || !status2 || !passed_DAS){
    delete chain;
    chain = (TChain*) new TChain("Events");
    for(int i = 0; i < Nfile; i++)
      chain->Add(filenames[i].c_str());
    NEVENT = chain->GetEntries();
    if(NEVENT == 0) return 1;
    chain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
    chain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
    
    if(string(FileTag).find("130X") != std::string::npos) // Run3
      chain->SetBranchAddress("nGenPart", &nGenPart_Run3, &b_nGenPart);
    else
      chain->SetBranchAddress("nGenPart", &nGenPart_Run2, &b_nGenPart);

    chain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
    chain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
    chain->SetBranchStatus("*",0);
    chain->SetBranchStatus("genWeight", 1);
    chain->SetBranchStatus("luminosityBlock", 1), 
    chain->SetBranchStatus("nGenPart", 1);
    chain->SetBranchStatus("GenPart_mass", 1);
    chain->SetBranchStatus("GenPart_pdgId", 1);

    for(int e = 0; e < NEVENT; e++){
      int mymod = NEVENT/10;
      if(mymod < 1)
        mymod = 1;

      chain->GetEntry(e);
      if(e%mymod == 0)
        cout << " event = " << e << " : " << NEVENT << endl;

      Nevent += 1.;
      Nweight += genWeight;
     
      MP = 0;
      MC = 0;
      int Ngen = 0;
      if(string(FileTag).find("130X") != std::string::npos) // Run3
        Ngen = nGenPart_Run3;
      else
        Ngen = nGenPart_Run2;
      if(Ngen > 400){
        cout << "weight? " << genWeight << endl;
      }
      //cout << "NGEN " << Ngen << endl;
      if(Ngen > maxNGEN)
        maxNGEN = Ngen;
      // for cascades
      int code = 0;
      bool has_Slep = false;
      bool has_Snu = false;
      bool is_left = true;
      // minus referring to charge of e or mu (not value of PDGID)
      bool is_plus = true;
      for(int i = 0; i < Ngen; i++){
        PDGID = fabs(GenPart_pdgId[i]);
        if(PDGID > 1000000 && PDGID < 3000000){
          int mass = int(GenPart_mass[i]+0.5);
          if(PDGID == 1000022)
            MC = mass;
          else
            if(mass > MP)
              MP = mass;
        }
        // Getting 'code' for cascades
        if(DO_CASCADES && DO_SMS){
          if (abs(PDGID) % 10000 == 11 || abs(PDGID) % 10000 == 13) {
            has_Slep = true;
            if (PDGID > 0) is_plus = false;
          } 
          else if (abs(PDGID) % 10000 == 12 || abs(PDGID) % 10000 == 14) {
              has_Snu = true;
          }
          if (PDGID > 2000000) is_left = false;
        }
      } // for(int i = 0; i < Ngen; i++)
      if(DO_CASCADES && DO_SMS){
        // build code from booleans
        if(has_Slep && !has_Snu) code += 1; // SlepSlep
        if(has_Slep && has_Snu) code += 2; // SlepSnu
        if(!has_Slep && has_Snu) code += 3; // SnuSnu
        if(is_left) code += 10;
        else code += 20;
        if(is_plus) code += 100;
        else code += 200;
      }
      for(int i = 0; i < Ngen; i++){
        PDGID = abs(GenPart_pdgId[i]);
        if(PDGID > 1000000 && PDGID < 3000000){
          int mass = int(GenPart_mass[i]+0.5);
          if(PDGID == 1000022)
            MC = mass;
          else
            if(mass > MP)
              MP = mass;
        }
      }
      Code = code;
      std::pair<int,std::pair<int,int>> masspair;
      masspair = std::make_pair(Code,std::make_pair(MP,MC));
      if(mapNevent.count(masspair) == 0){
        masses.push_back(masspair);
        mapNevent[masspair]  = 0.;
        mapNweight[masspair] = 0.;
      }

      mapNevent[masspair]  += 1.;
      mapNweight[masspair] += genWeight;
    }
    cout << "MAX NGEN " << maxNGEN << endl;
  }

  TFile* fout = new TFile(string(outputFileName).c_str(),"RECREATE");
  TTree* tout = (TTree*) new TTree("EventCount", "EventCount");
  fout->SetCompressionLevel(0);
  
  tout->Branch("NDAS", &NDAS);
  tout->Branch("Nevent", &Nevent);
  tout->Branch("Nweight", &Nweight);
  tout->Branch("filetag", &filetag);
  tout->Branch("dataset", &dataset);
  tout->Branch("DAS_datasetname", &DAS_datasetname);
  tout->Branch("DAS_filename", &DAS_filename);
  tout->Branch("MP", &MP);
  tout->Branch("MC", &MC);
  tout->Branch("Code", &Code);
  if(DO_SMS){
    int Nmass = masses.size();
    for(int i = 0; i < Nmass; i++){
        Nevent_tot += mapNevent[masses[i]];
    }
  }
  else Nevent_tot = Nevent;
  passed_DAS = true;
  if(!DO_PRIVATEMC){
    if(NDAS != Nevent_tot){ 
      std::cout << "JOB FAILED DAS CHECK!" << std::endl;
      std::cout << "  NDAS: " << NDAS << std::endl;
      std::cout << "  Nevent_tot: " << Nevent_tot << std::endl;
      passed_DAS = false;
    } else {
      std::cout << "JOB PASSED DAS CHECK!" << std::endl;
    }
  }

  if(DO_SMS){
    int Nmass = masses.size();
    for(int i = 0; i < Nmass; i++){
      Nweight = mapNweight[masses[i]];
      Nevent = mapNevent[masses[i]];
      NDAS = Nevent; // since we already passed DAS check above, set NDAS to Nevent for filling tree
      Code = masses[i].first;
      MP = masses[i].second.first;
      MC = masses[i].second.second;
      tout->Fill();
    }
  } else {
    tout->Fill();
  }

  fout->cd();
  tout->Write("", TObject::kOverwrite);
  TDirectory* metaDir = fout->GetDirectory("meta");
  if(!metaDir)
    metaDir = fout->mkdir("meta");
  metaDir->cd();
  TNamed param("GitCommitHash", gitHash);
  param.Write();
  fout->cd();
  fout->Close();
  if(passed_DAS) return 0;
  else return 1;
}
