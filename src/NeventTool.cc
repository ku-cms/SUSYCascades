#include "../include/NeventTool.hh"

#include <string>

using std::string;
using std::cout;
using std::endl;

std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim)
{
 unsigned first_delim_pos = s.find(start_delim);
 unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
 unsigned last_delim_pos = s.find_first_of(stop_delim, end_pos_of_first_delim);
 return s.substr(end_pos_of_first_delim,last_delim_pos - end_pos_of_first_delim);
}

NeventTool::NeventTool(){

  m_dataset = new std::string();
  m_filetag = new std::string();
  
}

NeventTool::~NeventTool(){

}

void NeventTool::BuildMap(const std::string& rootfile){
  m_RootFile = rootfile;

  std::cout << "Opening EventCount file " << m_RootFile << std::endl;
  m_File = new TFile(m_RootFile.c_str(), "READ");
  m_Tree = nullptr;
  if(!m_File->IsOpen()){
    std::cout << "Failed to open file" << std::endl;
    return;
  }
  m_Tree = (TTree*) m_File->Get("EventCount");
  if(m_Tree == nullptr){
    std::cout << "Unable to find EventCount tree" << std::endl;
    return;
  }

  m_Tree->SetMakeClass(1);
  m_Tree->SetBranchAddress("Nevent", &m_Nevent, &b_m_Nevent);
  m_Tree->SetBranchAddress("NDAS", &m_NDAS, &b_m_NDAS);
  m_Tree->SetBranchAddress("Nweight", &m_Nweight, &b_m_Nweight);
  m_Tree->SetBranchAddress("dataset", &m_dataset, &b_m_dataset);
  m_Tree->SetBranchAddress("filetag", &m_filetag, &b_m_filetag);
  m_Tree->SetBranchAddress("MP", &m_MP, &b_m_MP);
  m_Tree->SetBranchAddress("MC", &m_MC, &b_m_MC);
  m_Tree->SetBranchAddress("Code", &m_Code, &b_m_Code);
}

void NeventTool::BuildFilterEffMap(const std::string& rootfile){
  m_RootFile_FE = rootfile;

  std::cout << "Opening FilterEff file " << m_RootFile_FE << std::endl;
  m_File_FE = new TFile(m_RootFile_FE.c_str(), "READ");
  m_Tree_FE = nullptr;
  if(!m_File_FE->IsOpen()){
    std::cout << "Failed to open file" << std::endl;
    return;
  }
  m_Tree_FE = (TTree*) m_File_FE->Get("FilterEff");
  if(m_Tree_FE == nullptr){
    std::cout << "Unable to find FilterEff tree" << std::endl;
    return;
  }

  m_Tree_FE->SetMakeClass(1);
  m_Tree_FE->SetBranchAddress("efficiency", &m_efficiency, &b_m_efficiency);
  m_Tree_FE->SetBranchAddress("lumiblock", &m_lumiblock, &b_m_lumiblock);
  m_Tree_FE->SetBranchAddress("dataset", &m_dataset, &b_m_dataset);
  m_Tree_FE->SetBranchAddress("filetag", &m_filetag, &b_m_filetag);
}

void NeventTool::Initialize_BKG(const std::string& dataset, const std::string& filetag) const {

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Tree == nullptr){
    m_Label2Nevent_BKG[label] = 0;
    m_Label2Nweight_BKG[label] = 0;
    return;
  }

  double Nevent = 0;
  double Nweight = 0;
  
  int N = m_Tree->GetEntries();
  for(int i = 0; i < N; i++){
    m_Tree->GetEntry(i);

    if((dataset == (*m_dataset)) &&
       (filetag == (*m_filetag))){
      Nevent += m_Nevent;
      Nweight += m_Nweight;
    }
  }

  m_Label2Nevent_BKG[label] = Nevent;
  m_Label2Nweight_BKG[label] = Nweight;
  
}

void NeventTool::Initialize_SMS(const std::string& dataset, const std::string& filetag) const {

  std::pair<std::string,std::string> label(dataset,filetag);

  m_Label2Nevent_SMS[label] = std::map<std::pair<int,std::pair<int,int>>,double>();
  m_Label2Nweight_SMS[label] = std::map<std::pair<int,std::pair<int,int>>,double>();
  
  if(m_Tree == nullptr){
    return;
  }

  double Nevent = 0;
  double Nweight = 0;

  //cout << "Initializing " << dataset << " " << filetag << endl;
  
  int N = m_Tree->GetEntries();
  for(int i = 0; i < N; i++){
    m_Tree->GetEntry(i);
    
    //cout << (*m_dataset) << " " << (*m_filetag) << endl;

    if((dataset == (*m_dataset)) &&
       (filetag == (*m_filetag))){

      std::pair<int,std::pair<int,int>> masses;
      masses = std::make_pair(m_Code,std::make_pair(m_MP,m_MC));
      if(m_Label2Nevent_SMS[label].count(masses) == 0){
	m_Label2Nevent_SMS[label][masses] = 0.;
	m_Label2Nweight_SMS[label][masses] = 0.;
      } else {
//cout << "HEREEREERE" << endl;
      }
      m_Label2Nevent_SMS[label][masses] += m_Nevent;
      m_Label2Nweight_SMS[label][masses] += m_Nweight;
    }
  }
}

void NeventTool::Initialize_FE(const std::string& dataset, const std::string& filetag) const {

  //std::pair<std::string,std::string> label(dataset,filetag);
  std::string label = dataset;

  m_Label2FilterEff[label] = std::map<int,double>();
  
  if(m_Tree_FE == nullptr){
    return;
  }
  
  int N = m_Tree_FE->GetEntries();
  for(int i = 0; i < N; i++){
    m_Tree_FE->GetEntry(i);
    
    if((dataset == (*m_dataset)) &&
       (filetag == (*m_filetag))){

      m_Label2FilterEff[label][m_lumiblock] = m_efficiency;
      
    }
  }
}

double NeventTool::GetNevent_BKG(const std::string& dataset, const std::string& filetag) const {
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nevent_BKG.count(label) == 0)
    Initialize_BKG(dataset, filetag);

  return m_Label2Nevent_BKG[label];
}

double NeventTool::GetNevent_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC, int Code) const {
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nevent_SMS.count(label) == 0)
    Initialize_SMS(dataset, filetag);

  std::pair<int,std::pair<int,int>> masses;
  masses = std::make_pair(Code,std::make_pair(MP,MC));

  if(m_Label2Nevent_SMS[label].count(masses) == 0)
    return 0.;
  
  return m_Label2Nevent_SMS[label][masses];
}

double NeventTool::GetNweight_BKG(const std::string& dataset, const std::string& filetag) const {
  std::cout << "Getting Nevent count for " << dataset << " " << filetag << std::endl;
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nweight_BKG.count(label) == 0)
    Initialize_BKG(dataset, filetag);

  double Nweight = m_Label2Nweight_BKG[label];
  std::cout << "Nweight from event tool is " << Nweight << std::endl;

  return Nweight;
}

double NeventTool::GetFilterEff(const std::string& dataset, const std::string& filetag, int lumiblock) const {
 
  if(m_Label2FilterEff.count(dataset) == 0)
    Initialize_FE(dataset, filetag);

  if(m_Label2FilterEff[dataset].count(lumiblock) != 0)
    return m_Label2FilterEff[dataset][lumiblock];
  
  if(m_Label2FilterEff[dataset].count(-1) == 0)
    m_Label2FilterEff[dataset][-1] = 1.;

  return m_Label2FilterEff[dataset][-1];
}

double NeventTool::GetNweight_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC, int Code) const {
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nweight_SMS.count(label) == 0)
    Initialize_SMS(dataset, filetag);
  
  std::pair<int,std::pair<int,int>> masses;
  masses = std::make_pair(Code,std::make_pair(MP,MC));

  if(m_Label2Nweight_SMS[label].count(masses) == 0)
    return 0.;
  return m_Label2Nweight_SMS[label][masses];
}

std::map<std::pair<std::string,std::string>,double> NeventTool::InitMap_Nevent_BKG(){
  std::map<std::pair<std::string,std::string>,double> Label2Nevent;

  return Label2Nevent;
}

std::map<std::pair<std::string,std::string>,double> NeventTool::InitMap_Nweight_BKG(){
  std::map<std::pair<std::string,std::string>,double> Label2Nweight;
  
  return Label2Nweight;
}

std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> NeventTool::InitMap_Nevent_SMS(){
  std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> Label2Nevent;
  
  return Label2Nevent;
}

std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> NeventTool::InitMap_Nweight_SMS(){
  std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> Label2Nweight;

  return Label2Nweight;
}

bool check_dataset_file(std::string dataset_name)
{
 std::ifstream testfile(dataset_name);
 if(testfile.peek() == std::ifstream::traits_type::eof())
 {
  testfile.close();
  return false;
 }
 testfile.close();
 return true;
}

std::string NeventTool::Get_DASdatasetname(const std::string& u_file)
{
 std::string filename = u_file;
 size_t pos = filename.find("/store/");
 if(pos != std::string::npos) filename = filename.substr(pos); // remove everything before root redirector for DAS query
 TString das_output = gSystem->GetFromPipe(("dasgoclient -query=\"dataset file="+filename+"\"").c_str());
 return std::string(das_output.Data());
}

int NeventTool::EventsInDAS(const std::string& u_file)
{
 int Events = 0;
 std::string filename = u_file;
 size_t pos = filename.find("/store/");
 if(pos != std::string::npos) filename = filename.substr(pos); // remove everything before root redirector for DAS query
 TString das_output = gSystem->GetFromPipe(("dasgoclient -query=\"file="+filename+"\" -json").c_str());
 if(das_output.First("nevents"))
 {
  string events_string = get_str_between_two_str(das_output.Data(),"nevents\":",",");
  Events = std::stoi(events_string);
 }
 return Events;
}

std::map<std::string,std::map<int,double> > NeventTool::InitMap_FilterEff(){
  std::map<std::string,std::map<int,double> > Label2FilterEff;

  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP2_13TeV-madgraphMLM-pythia8"] = std::map<int,double>();
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8"] = std::map<int,double>();
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP2_13TeV-madgraphMLM-pythia8"] = std::map<int,double>();
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8"] = std::map<int,double>();
  
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP2_13TeV-madgraphMLM-pythia8"][-1] = 0.465;
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-420_TuneCP5_13TeV-madgraphMLM-pythia8"][-1] = 0.465;
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP2_13TeV-madgraphMLM-pythia8"][-1] = 0.496;
  Label2FilterEff["SMS-T2-4bd_genMET-80_mStop-500_mLSP-490_TuneCP5_13TeV-madgraphMLM-pythia8"][-1] = 0.496;
 
  return Label2FilterEff;
}

std::map<std::pair<std::string,std::string>,double> NeventTool::m_Label2Nevent_BKG  = NeventTool::InitMap_Nevent_BKG();
std::map<std::pair<std::string,std::string>,double> NeventTool::m_Label2Nweight_BKG = NeventTool::InitMap_Nweight_BKG();
std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> NeventTool::m_Label2Nevent_SMS  = NeventTool::InitMap_Nevent_SMS();
std::map<std::pair<std::string,std::string>,std::map<std::pair<int,std::pair<int,int>>,double>> NeventTool::m_Label2Nweight_SMS = NeventTool::InitMap_Nweight_SMS();

std::map<std::string,std::map<int,double> > NeventTool::m_Label2FilterEff = NeventTool::InitMap_FilterEff();
