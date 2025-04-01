#include <iostream>
#include <map>

#include "Process.hh"
#include "Category.hh"

///////////////////////////////////////////
////////// Process class
///////////////////////////////////////////

Process::Process(const string& title, ProcessType type){
  m_Title = title;
  m_Type = type;
}

Process::Process(const Process& proc){
  m_Title = proc.Name();
  m_Type  = proc.Type();
}

Process::~Process(){
  Clear();
}

void Process::Clear(){
  auto p = m_ProcBins.begin();
  while(p != m_ProcBins.end()){
    auto c = p->second.begin();
    while(c != p->second.end()){
      delete c->second;
      c++;
    }
    p++;
  }
}

const Process& Process::operator = (const Process& proc){
  Clear();

  m_Title = proc.Name();
  m_Type  = proc.Type();

  return *this;
}


const string& Process::Name() const {
  return m_Title;
}

ProcessType Process::Type() const {
  return m_Type;
}

bool Process::AddEvent(double weight, double Mperp, double RISR,
		       const Category& cat, const Systematic& sys, bool extrahist){

  string clabel = cat.Label()+"_"+cat.GetLabel();
  string plabel = sys.TreeName(Name());
  
  if(m_ProcBins.count(plabel) == 0)
    m_ProcBins[plabel] = map<string,FitBin*>();
    
  if(m_ProcBins[plabel].count(clabel) == 0)
    m_ProcBins[plabel][clabel] = cat.GetNewFitBin(plabel, extrahist);
      
  return m_ProcBins[plabel][clabel]->Fill(weight, Mperp, RISR);
}

////// CHECK systs includes up and down? where is map filled
void Process::AddShapeSysts(const Systematics& systs) {
  for(int s = 0; s < systs.GetN(); s++){
    string plabel = systs[s].TreeName(Name());
    
    if(m_ProcBins.count(plabel) == 0)
      m_ProcBins[plabel] = map<string,FitBin*>();
  }
}

Process Process::FakeProcess(const string& label) const {
  return Process(Name()+"_"+label, m_Type);
}

bool Process::operator <  (const Process& proc) const {
  return Name() < proc.Name();
}

bool Process::operator >  (const Process& proc) const {
  return Name() > proc.Name();
}

bool Process::operator == (const Process& proc) const {
  return Name() == proc.Name();
}

SM Process::GetSM() const {
  if(Type() != kSig)
    return SM(Name());
  
  size_t l = Name().rfind("_");
  if(l == std::string::npos)
    return SM(Name());

  SM sm(Name().substr(0, l));
  sm += Name().substr(l+1,Name().length()-l);

  return sm;
}

///////////////////////////////////////////
////////// ProcessList class
///////////////////////////////////////////

ProcessList::ProcessList(){
  m_N = 0;
}

ProcessList::ProcessList(const ProcessList& list){
  m_N = 0;
  *this += list;
}

ProcessList::~ProcessList() {}

ProcessList& ProcessList::operator = (const ProcessList& list){
  ProcessList pl(list);
  m_ProcMap.clear();
  m_Proc.clear();
  m_N = 0;
  return *this += pl;
}

ProcessList& ProcessList::operator += (const ProcessList& list){
  int N = list.GetN();
  for(int i = 0; i < N; i++)
    *this += list[i];

  return *this;
}

ProcessList& ProcessList::operator += (const Process& proc){
  if(m_ProcMap.count(proc) != 0)
    return *this;
  
  m_N++;
  m_Proc.push_back(proc);
  m_ProcMap[proc] = true;

  return *this;
}

ProcessList  ProcessList::operator + (const ProcessList& list) const {
  ProcessList ret = *this;
  ret += list;

  return ret;
}

ProcessList ProcessList::Filter(ProcessType type) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++)
    if(m_Proc[i].Type() == type)
      list += m_Proc[i];

  return list;
}

ProcessList ProcessList::Filter(const string& label) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++)
    if(m_Proc[i].Name().find(label) != std::string::npos)
      list += m_Proc[i];

  return list;
}

ProcessList ProcessList::Remove(ProcessType type) const {
  ProcessList list;
  
  for(int i = 0; i < m_N; i++)
    if(m_Proc[i].Type() != type)
      list += m_Proc[i];

  return list;
}

ProcessList ProcessList::Remove(const string& label) const {
   ProcessList list;

   for(int i = 0; i < m_N; i++)
     if(m_Proc[i].Name().find(label) == std::string::npos)
       list += m_Proc[i];
   
  return list;
}

ProcessList ProcessList::FilterOR(VS& labels) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++)
    for(auto l : labels)
      if(m_Proc[i].Name().find(l) != std::string::npos){
	list += m_Proc[i];
	break;
      }

  return list;
}

ProcessList ProcessList::FilterAND(VS& labels) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++){
    bool match = true;
    for(auto l : labels)
      if(m_Proc[i].Name().find(l) == std::string::npos){
	match = false;
	break;
      }
    if(match)
      list += m_Proc[i];
  }

  return list;
}

ProcessList ProcessList::RemoveOR(VS& labels) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++){
    bool match = false;
    for(auto l : labels)
      if(m_Proc[i].Name().find(l) != std::string::npos){
	match = true;
	break;
      }
    if(!match)
      list += m_Proc[i];
  }

  return list;  
}

ProcessList ProcessList::RemoveAND(VS& labels) const {
  ProcessList list;

  for(int i = 0; i < m_N; i++)
    for(auto l : labels)
      if(m_Proc[i].Name().find(l) == std::string::npos){
	list += m_Proc[i];
	break;
      }

  return list;
}
  
int ProcessList::GetN() const {
  return m_N;
}

Process ProcessList::operator [] (int i) const {
  if(i < 0 || i >= m_N)
    return Process("empty",kBkg);

  return m_Proc[i];
}

int ProcessList::Find(const string& label) const {
  for(int i = 0; i < m_N; i++)
    if(m_Proc[i].Name() == label)
      return i;

  return -1;
}

VS ProcessList::GetProcesses() const {
  VS list;

  for(int i = 0; i < m_N; i++)
    list += m_Proc[i].Name();

  return list;
}

VSM ProcessList::GetSignalMasses() const {
  VSM signals;

  string name;
  for (int p = 0; p < m_N; p++) {
    if (m_Proc[p].Type() != kSig)
      continue;

    name = m_Proc[p].Name();
    size_t l = name.rfind("_");
    if (l == std::string::npos)
      continue;

    // Extract the part before the last "_"
    SM sm(name.substr(0, l));

    // Extract the numeric parts after the last "_", and combine them with a 0 separator
    string combinedNumbers = "";
    size_t start = l + 1;
    while (start < name.length()) {
      size_t end = name.find("_", start);
      if (end == std::string::npos) end = name.length();
      
      // Extract the numeric part and append it to the combined string with "0" separator
      string number = name.substr(start, end - start);
      if (!combinedNumbers.empty()) {
        combinedNumbers += "0"; // Add separator before appending the next number
      }
      combinedNumbers += number;

      start = end + 1;
    }

    sm += combinedNumbers; // Append the combined number to the signal name
    signals += sm; 
  }

  return signals;
}


void ProcessList::Print() const {
   for(int p = 0; p < m_N; p++)
     std::cout << m_Proc[p].Name() << std::endl;
}

///////////////////////////////////////////
////////// ProcessBranch class
///////////////////////////////////////////

ProcessBranch::ProcessBranch() {
  m_Tree = nullptr;
}

ProcessBranch::~ProcessBranch() {}

void ProcessBranch::InitFill(TTree* tree){
  if(!tree)
    return;

  tree->Branch("proc",    &m_Proc);
  tree->Branch("subproc", &m_SubProc);
  tree->Branch("type",    &m_ProcType);

  m_Tree = tree;
}

void ProcessBranch::FillProcess(const Process& proc, const Systematic& sys){
  if(!m_Tree)
    return;

  m_Proc     = proc.Name();
  m_ProcType = proc.Type();

  if(sys.IsDefault()){
    m_SubProc  = proc.Name();
    m_Tree->Fill();
  } else {
    m_SubProc  = proc.Name()+"_"+sys.Label()+"Down";
    m_Tree->Fill();
    m_SubProc  = proc.Name()+"_"+sys.Label()+"Up";
    m_Tree->Fill();
  }
}

void ProcessBranch::FillProcess(Process& proc, TFile& file){
  if(!file.IsOpen())
    return;

  m_Proc     = proc.Name();
  m_ProcType = proc.Type();
  
  // make histogram output directories for this process
  //file.cd();
  //file.mkdir(proc.Name().c_str());
  
  // loop through all subprocesses (systematics)
  auto p = proc.m_ProcBins.begin();
  while(p != proc.m_ProcBins.end()){
    m_SubProc = p->first;
    
    // loop through all categories for a given subprocess
    auto c = p->second.begin();
    while(c != p->second.end()){
      // write FitBin to output file for each subprocess/category
      file.cd();
      if(!file.GetDirectory(c->first.c_str()))
	file.mkdir(c->first.c_str());
      file.cd();
	 
      //c->second->WriteHistogram(m_SubProc+"_"+c->first, m_Proc, file);
      c->second->WriteHistogram(m_SubProc, c->first, file);
      // clean up this bin (assuming we won't need it after writing...)
      c->second->Clear();
      
      c++;
    }

    if(m_Tree)
      m_Tree->Fill();
    
    p++;
  }
}

void ProcessBranch::InitGet(TTree* tree){
  m_ProcPtr = 0;
  m_SubProcPtr = 0;
  tree->SetBranchAddress("proc", &m_ProcPtr, &m_b_Proc);
  tree->SetBranchAddress("subproc", &m_SubProcPtr, &m_b_SubProc);
  tree->SetBranchAddress("type", &m_ProcType, &m_b_ProcType);
}
  

Process ProcessBranch::GetProcess(){
  if(!m_ProcPtr || !m_SubProcPtr)
    return Process(m_SubProc, ProcessType(m_ProcType));
  else {
    return Process(*m_SubProcPtr, ProcessType(m_ProcType));
  }
}


