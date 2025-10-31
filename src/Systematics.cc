#include "Systematics.hh"
#include "CategoryTree.hh"
///////////////////////////////////////////
////////// Systematic class
///////////////////////////////////////////

Systematic::Systematic(const string& label){
  m_Label = label;
  m_IsUp = true;
}
    
Systematic::~Systematic(){ }

std::string Systematic::Label() const {
  return m_Label;
}

const Systematic& Systematic::Up() const {
  m_IsUp = true;

  return *this;
}

const Systematic& Systematic::Down() const {
  m_IsUp = false;

  return *this;
}

bool Systematic::IsUp() const {
  return m_IsUp;
}

bool Systematic::IsDown() const {
  return !m_IsUp;
}

bool Systematic::IsDefault() const {
  return IsSame("Default");
}
 
bool Systematic::IsSame(const string& label) const {
  return m_Label == label;
}

bool Systematic::IsSame(const Systematic& sys) const {
  return IsSame(sys.Label());
}

bool Systematic::operator == (const string& label) const {
  return IsSame(label);
}

bool Systematic::operator == (const Systematic& sys) const {
  return IsSame(sys);
}

bool Systematic::operator != (const string& label) const {
  return !IsSame(label);
}

bool Systematic::operator != (const Systematic& sys) const { return
    !IsSame(sys);
}

bool Systematic::operator == (const Systematics& sys) const {
  return sys == *this;
}

bool Systematic::operator != (const Systematics& sys) const {
  return sys != *this;
}

bool Systematic::operator < (const Systematic& sys) const {
  return (Label() < sys.Label());
}

bool Systematic::operator > (const Systematic& sys) const {
  return (Label() > sys.Label());
}

std::string Systematic::TreeName(const string& name) const {
  if(m_Label.compare("METUncer_GenMET") == 0)
    return name+"_"+m_Label;
  if(this->IsDefault())
    return name;
  else
    return name+"_"+m_Label+(IsUp() ? "Up" : "Down");
}

Systematic& Systematic::Default(){
  return Systematic::m_Default;
}

Systematic Systematic::m_Default;

///////////////////////////////////////////
////////// Systematics class
///////////////////////////////////////////

Systematics::Systematics(bool include_default){
  m_N = 0;
  
  if(include_default)
    Add(Systematic::Default());
}

Systematics::Systematics(const Systematics& sys){
  m_N = 0;
  int N = sys.GetN();
  for(int i = 0; i < N; i++)
    *this += sys[i];
}

    
Systematics::~Systematics(){
  Clear();
}

void Systematics::Clear(){
  for(int i = 0; i < m_N; i++)
    delete m_Sys[i];
  m_Sys.clear();
  m_N = 0;
}

Systematics& Systematics::operator = (const Systematics& sys){
  Clear();
  
  int N = sys.GetN();
  for(int i = 0; i < N; i++)
    *this += sys[i];

  return *this;
}

Systematic& Systematics::operator [] (int i) const {
  if(i < 0 || i >= m_N)
    return Systematic::Default();

  return *m_Sys[i];
}

int Systematics::GetN() const {
  return m_N;
}

Systematics& Systematics::Add(const string& label){
  if(Contains(label))
    return *this;
      
  m_Sys.push_back(new Systematic(label));
  m_N++;
  
  return *this;
}

Systematics& Systematics::Add(const Systematic& sys){
  if(Contains(sys))
    return *this;
  
  m_Sys.push_back(new Systematic(sys.Label()));
  m_N++;
  
  return *this;
}

Systematics& Systematics::Add(const Systematics& sys){
  int N = sys.GetN();

  for(int i = 0; i < N; i++)
    Add(sys[i]);
  
  return *this;
}

Systematics& Systematics::operator += (const string& label){
  return Add(label);
}

Systematics& Systematics::operator += (const Systematic& sys){
  return Add(sys);
}

Systematics& Systematics::operator += (const Systematics& sys){
  return Add(sys);
}

bool Systematics::Contains(const string& label) const {
  for(int i = 0; i < m_N; i++)
    if(*m_Sys[i] == label)
      return true;

  return false;
}

bool Systematics::Contains(const Systematic& sys) const {
  for(int i = 0; i < m_N; i++)
    if(*m_Sys[i] == sys)
      return true;

  return false;
}

bool Systematics::operator == (const string& label) const {
  return Contains(label);
}

bool Systematics::operator == (const Systematic& sys) const {
  return Contains(sys);
}

bool Systematics::operator != (const string& label) const {
  return !Contains(label);
}

bool Systematics::operator != (const Systematic& sys) const {
  return !Contains(sys);
}

Systematics Systematics::Filter(const string& label) const {
  Systematics list;
  for(auto s : m_Sys)
    if(s->Label().find(label) != std::string::npos)
      list += *s;

  return list;
}

Systematics Systematics::Remove(const string& label) const {
  Systematics list;
  for(auto s : m_Sys)
    if(s->Label().find(label) == std::string::npos)
      list += *s;

  return list;
}

Systematics Systematics::FilterOR(const VS& labels) const {
  Systematics list;
  for(auto s : m_Sys)
    for(auto l : labels)
      if(s->Label().find(l) != std::string::npos){
	list += *s;
	break;
      }
  
  return list;
}

Systematics Systematics::FilterAND(const VS& labels) const {
  Systematics list;
  for(auto s : m_Sys){
    bool match = true;
    for(auto l : labels)
      if(s->Label().find(l) == std::string::npos){
	match = false;
	break;
      }
    if(match)
      list += *s;
  }
  
  return list;
}

Systematics Systematics::RemoveOR(const VS& labels) const {
  Systematics list;
  for(auto s : m_Sys){
    bool match = false;
    for(auto l : labels)
      if(s->Label().find(l) != std::string::npos){
	match = true;
	break;
      }
    if(!match)
      list += *s;
  }
  
  return list;
}


Systematics Systematics::RemoveAND(const VS& labels) const {
  Systematics list;
  for(auto s : m_Sys)
    for(auto l : labels)
      if(s->Label().find(l) == std::string::npos){
	list += *s;
	break;
      }
  
  return list;
}

///////////////////////////////////////////
////////// SystematicsTool class
///////////////////////////////////////////

SystematicsTool::SystematicsTool(){
  Init();
}

SystematicsTool::~SystematicsTool() {}


//get list to convert shape systematics to normalizations
Systematics SystematicsTool::GetConvertedSystematics() const {
   Systematics list;
   
  list += "BTAGHF_SF";
  list += "BTAGLF_SF";
 
return list;
}

// Default weight/SF systematics
Systematics SystematicsTool::GetWeightSystematics() const {
  Systematics list;
  
  //list += "PU_SF"; // turn off for now
  list += "BTAGHF_SF";
  list += "BTAGLF_SF";
  //list += "MET_TRIG_SF"; // default
  ////list += "MET_TRIG_el";
  ////list += "MET_TRIG_mu";

  return list;
}

//pass this list of real procs ie ttbar, Wjets
//Systematics SystematicsTool::GetFakeShapeSystematics(CategoryTree CT, VS procs) const {
//        Systematics list;
//	vector<const CategoryTree*> catTrees;
//	CT.GetListDepth(catTrees,1);
//
//	VS procGroups;
//	if(std::count(procs.begin(),procs.end(),"ttbar") || std::count(procs.begin(),procs.end(),"ST")) procGroups += "ttbarST";
//	if(std::count(procs.begin(),procs.end(),"Wjets") || std::count(procs.begin(),procs.end(),"TB") || std::count(procs.begin(),procs.end(),"DB")) procGroups += "WjetsDBTB";
//	if(std::count(procs.begin(),procs.end(),"ZDY")) procGroups += "ZDY";
//	if(std::count(procs.begin(),procs.end(),"QCD")) procGroups += "QCD";
//
//
//	for(int c = 0; c < int(catTrees.size()); c++){
//		for(int p = 0; p < procGroups.size(); p++){
//			list += Systematic(procGroups[p]+"_"+catTrees[c]->GetSpecLabel()+"f0_RISR").Up();
//			list += Systematic(procGroups[p]+"_"+catTrees[c]->GetSpecLabel()+"f1_RISR").Down();
//			list += Systematic(procGroups[p]+"_"+catTrees[c]->GetSpecLabel()+"f0_Mperp").Up();
//			list += Systematic(procGroups[p]+"_"+catTrees[c]->GetSpecLabel()+"f1_Mperp").Down();
//		}
//	}
//  
//
//  return list;
//}

// Default alternative tree systematics
Systematics SystematicsTool::GetTreeSystematics() const {
  
  Systematics list;
  
  //list += "JESUncer_Full";
  //list += "JESUncer_Reduced";
  list += "JESUncer_Total";
  list += "JERUncer_Total";
  list += "METUncer_UnClust";
  list += "METUncer_GenMET";
  // list += "MMSUncer_Total";
  // list += "EESUncer_Total";

  return list;
}

const Systematics& SystematicsTool::JESSystematics() const {
  return m_JESSys;
}

const Systematics& SystematicsTool::JERSystematics() const {
  return m_JERSys;
}

const Systematics& SystematicsTool::MMSSystematics() const {
  return m_MMSSys;
}

const Systematics& SystematicsTool::EESSystematics() const {
  return m_EESSys;
}

const Systematics& SystematicsTool::METSystematics() const {
  return m_METSys;
}
  
void SystematicsTool::Init(){

  m_JESSys += "JESUncer_Total";
  m_JESSys += "JESUncer_AbsoluteStat";
  m_JESSys += "JESUncer_AbsoluteScale";
  m_JESSys += "JESUncer_AbsoluteFlavMap";
  m_JESSys += "JESUncer_AbsoluteMPFBias";
  m_JESSys += "JESUncer_Fragmentation";
  m_JESSys += "JESUncer_SinglePionECAL";
  m_JESSys += "JESUncer_SinglePionHCAL";
  m_JESSys += "JESUncer_FlavorQCD";
  m_JESSys += "JESUncer_TimePtEta";
  m_JESSys += "JESUncer_RelativeJEREC1";
  m_JESSys += "JESUncer_RelativeJEREC2";
  m_JESSys += "JESUncer_RelativeJERHF";
  m_JESSys += "JESUncer_RelativePtBB";
  m_JESSys += "JESUncer_RelativePtEC1";
  m_JESSys += "JESUncer_RelativePtEC2";
  m_JESSys += "JESUncer_RelativePtHF";
  m_JESSys += "JESUncer_RelativeBal";
  m_JESSys += "JESUncer_RelativeSample";
  m_JESSys += "JESUncer_RelativeFSR";
  m_JESSys += "JESUncer_RelativeStatFSR";
  m_JESSys += "JESUncer_RelativeStatEC";
  m_JESSys += "JESUncer_RelativeStatHF";
  m_JESSys += "JESUncer_PileUpDataMC";
  m_JESSys += "JESUncer_PileUpPtRef";
  m_JESSys += "JESUncer_PileUpPtBB";
  m_JESSys += "JESUncer_PileUpPtEC1";
  m_JESSys += "JESUncer_PileUpPtEC2";
  m_JESSys += "JESUncer_PileUpPtHF";
  m_JESSys += "JESUncer_PileUpMuZero";
  m_JESSys += "JESUncer_PileUpEnvelope";
  m_JESSys += "JESUncer_SubTotalPileUp";
  m_JESSys += "JESUncer_SubTotalRelative";
  m_JESSys += "JESUncer_SubTotalPt";
  m_JESSys += "JESUncer_SubTotalScale";
  m_JESSys += "JESUncer_SubTotalAbsolute";
  m_JESSys += "JESUncer_SubTotalMC";
  m_JESSys += "JESUncer_TotalNoFlavor";
  m_JESSys += "JESUncer_TotalNoTime";
  m_JESSys += "JESUncer_TotalNoFlavorNoTime";
  m_JESSys += "JESUncer_FlavorZJet";
  m_JESSys += "JESUncer_FlavorPhotonJet";
  m_JESSys += "JESUncer_FlavorPureGluon";
  m_JESSys += "JESUncer_FlavorPureQuark";
  m_JESSys += "JESUncer_FlavorPureCharm";
  m_JESSys += "JESUncer_FlavorPureBottom";
  m_JESSys += "JESUncer_TimeRunBCD";
  m_JESSys += "JESUncer_TimeRunEF";
  m_JESSys += "JESUncer_TimeRunGH";
  m_JESSys += "JESUncer_CorrelationGroupMPFInSitu";
  m_JESSys += "JESUncer_CorrelationGroupIntercalibration";
  m_JESSys += "JESUncer_CorrelationGroupbJES";
  m_JESSys += "JESUncer_CorrelationGroupFlavor";
  m_JESSys += "JESUncer_CorrelationGroupUncorrelated";

  m_MMSSys += "MMSUncer_Total";
  
  m_EESSys += "EESUncer_Total";

  m_JERSys += "JERUncer_Total";
  
  m_METSys += "METUncer_UnClust";
  m_METSys += "METUncer_GenMET";
  m_METSys.Add(m_JESSys);

}
