#include "../include/NeventTool.hh"

#include <string>

using std::string;
using std::cout;
using std::endl;

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

  m_Label2Nevent_SMS[label] = std::map<std::pair<int,int>,double>();
  m_Label2Nweight_SMS[label] = std::map<std::pair<int,int>,double>();
  
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

    //if( (*m_dataset).find("mStop") != string::npos){
    //	cout << (*m_dataset) << " " << (*m_filetag) << endl;
    //	cout << m_Nevent << " " << m_Nweight << endl;
    //	cout << m_MP << " " << m_MC << endl << endl;
    //}
    

    if((dataset == (*m_dataset)) &&
       (filetag == (*m_filetag))){

      std::pair<int,int> masses(m_MP,m_MC);
      
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

double NeventTool::GetNevent_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC) const {
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nevent_SMS.count(label) == 0)
    Initialize_SMS(dataset, filetag);

  std::pair<int,int> masses(MP,MC);

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

  int Nweight = m_Label2Nweight_BKG[label];
  std::cout << "Nweight is " << Nweight << std::endl;

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

double NeventTool::GetNweight_SMS(const std::string& dataset, const std::string& filetag, int MP, int MC) const {
  if(!m_Tree)
    return 0.;

  std::pair<std::string,std::string> label(dataset,filetag);
  
  if(m_Label2Nweight_SMS.count(label) == 0)
    Initialize_SMS(dataset, filetag);
  
  std::pair<int,int> masses(MP,MC);

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

std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > NeventTool::InitMap_Nevent_SMS(){
  std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > Label2Nevent;
  
  return Label2Nevent;
}

std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > NeventTool::InitMap_Nweight_SMS(){
  std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > Label2Nweight;

  return Label2Nweight;
}

std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim)
{
 unsigned first_delim_pos = s.find(start_delim);
 unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
 unsigned last_delim_pos = s.find_first_of(stop_delim, end_pos_of_first_delim);
 return s.substr(end_pos_of_first_delim,last_delim_pos - end_pos_of_first_delim);
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

bool NeventTool::DatasetIsFastSim(const std::string& infile_name)
{
 std::ifstream infile(infile_name);
 std::vector<string> datasets;
 std::string dataset = "";
 while (getline(infile,dataset))
  datasets.push_back(dataset);
 for(int i = 0; i < int(datasets.size()); i++)
 {
  if(datasets[i].find("FS") == std::string::npos)
  {
   infile.close();
   return false;
  }
 }
 infile.close();
 return true;
}

int NeventTool::EventsInDAS(const std::string& u_dataset, const std::string& u_filetag)
{
 std::string dataset = u_dataset;
 std::string filetag = u_filetag;
 if(filetag.find("_SMS") != std::string::npos)
   filetag.erase(filetag.length()-4);
 filetag.erase(filetag.length()-5);
 double Events = 0.;
 gSystem->Exec(("dasgoclient -query=\"dataset=/"+dataset+"/*"+filetag+"NanoAODv12*"+"/NANO*\" >> datasets_"+filetag+"_"+dataset+".txt").c_str());
 if(!check_dataset_file("datasets_"+filetag+"_"+dataset+".txt"))
   gSystem->Exec(("dasgoclient -query=\"dataset=/"+dataset+"/*"+filetag+"NanoAODv9*"+"/NANO*\" >> datasets_"+filetag+"_"+dataset+".txt").c_str());
 if(!check_dataset_file("datasets_"+filetag+"_"+dataset+".txt"))
   gSystem->Exec(("dasgoclient -query=\"dataset=/"+dataset+"/*"+filetag+"NanoAODv7*"+"/NANO*\" >> datasets_"+filetag+"_"+dataset+".txt").c_str());
 if(!check_dataset_file("datasets_"+filetag+"_"+dataset+".txt"))
   gSystem->Exec(("dasgoclient -query=\"dataset=/"+dataset+"/*"+filetag+"NanoAODv4*"+"/NANO*\" >> datasets_"+filetag+"_"+dataset+".txt").c_str());
 if(!check_dataset_file("datasets_"+filetag+"_"+dataset+".txt"))
   gSystem->Exec(("dasgoclient -query=\"dataset=/"+dataset+"/*"+filetag+"NanoAOD*"+"/NANO*\" >> datasets_"+filetag+"_"+dataset+".txt").c_str());

 std::string infile_name = "datasets_"+filetag+"_"+dataset+".txt";
 bool IsFS = DatasetIsFastSim(infile_name);
 std::ifstream infile(infile_name);
 std::string dataset_fullname = "";
 while(getline(infile,dataset_fullname))
 {
  if(dataset_fullname.find("JMENano") != std::string::npos) continue;
  if(!IsFS && dataset_fullname.find("FS") != std::string::npos) continue;
  gSystem->Exec(("dasgoclient -query=\"file dataset="+dataset_fullname+"\" -json >> "+filetag+"_"+dataset+".json").c_str());
 }
 infile.close();
 infile.open(filetag+"_"+dataset+".json");
 string line = "";
 while(getline(infile,line))
 {
  if(line.find("nevents") != std::string::npos)
  {
   string events_string = get_str_between_two_str(line,"nevents\":",",");
   Events += std::stod(events_string);
  }
 }
 //gSystem->Exec("rm datasets_*.txt");
 //gSystem->Exec("rm *.json");
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
std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > NeventTool::m_Label2Nevent_SMS  = NeventTool::InitMap_Nevent_SMS();
std::map<std::pair<std::string,std::string>,std::map<std::pair<int,int>,double> > NeventTool::m_Label2Nweight_SMS = NeventTool::InitMap_Nweight_SMS();

std::map<std::string,std::map<int,double> > NeventTool::m_Label2FilterEff = NeventTool::InitMap_FilterEff();
