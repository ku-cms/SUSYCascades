#include <iostream>
#include <map>

#include "../include/Leptonic.hh"

LepIDs::LepIDs() {}

LepIDs::~LepIDs() {}

int LepIDs::GetN() const {
  return m_IDs.size();
}

LepID LepIDs::Get(int i) const {
  if(i < 0 || i >= GetN())
    return m_IDs[0];

  return m_IDs[i];
}

LepID LepIDs::operator [] (int i) const {
  if(i < 0 || i >= GetN())
    return m_IDs[0];

  return m_IDs[i];
}

void LepIDs::Add(LepID id){
  m_IDs.push_back(id);
}

void LepIDs::operator += (LepID id){
  m_IDs.push_back(id);
}

LepIDs& LepIDs::a(LepID id){
  m_IDs.push_back(id);

  return *this;
}

bool LepIDs::IsEqual(const LepIDs& ids) const {
  if(GetN() != ids.GetN())
    return false;

  int N = GetN();
  for(int i = 0; i < N; i++)
    if(m_IDs[i] != ids[i])
      return false;

  return true;
}

bool LepIDs::operator == (const LepIDs& ids) const {
  return IsEqual(ids);
}

void LepIDs::swap(int index_a, int index_b){
  int N = GetN();
  if(index_a < 0   || index_b < 0 ||
     index_a > N-1 || index_b > N-1)
    return;

  LepID temp = m_IDs[index_a];
  m_IDs[index_a] = m_IDs[index_b];
  m_IDs[index_b] = temp;
}

///////////////////////////////////////////
////////// LepIDsList class
///////////////////////////////////////////

LepIDsList::LepIDsList(const string& label){
  m_Label = label;
}

LepIDsList::~LepIDsList(){ }

string LepIDsList::Label() const {
  return m_Label;
}

void LepIDsList::Add(const LepIDs& ids){
  if(!Contains(ids))
     m_IDs.push_back(ids);
}

void LepIDsList::operator += (const LepIDs& ids){
  Add(ids);
}

void LepIDsList::Add(const LepIDsList& ids){
  int N = ids.GetN();
  for(int i = 0; i < N; i++)
    Add(ids[i]);
}

void LepIDsList::operator += (const LepIDsList& ids){
  Add(ids);
}

int LepIDsList::GetN() const {
  return m_IDs.size();
}

const LepIDs& LepIDsList::Get(int i) const {
  if(i < 0 || i >= GetN())
    return m_IDs[0];

  return m_IDs[i];
}

const LepIDs& LepIDsList::operator [] (int i) const {
  if(i < 0 || i >= GetN())
    return m_IDs[0];

  return m_IDs[i];
}

bool LepIDsList::IsEqual(const LepIDsList& ids) const {
  if(GetN() != ids.GetN())
    return false;

  int N = GetN();
  for(int i = 0; i < N; i++)
    if(!Contains(ids[i]))
       return false;

  return true;
}

bool LepIDsList::operator == (const LepIDsList& ids) const {
  return IsEqual(ids);
}

bool LepIDsList::Contains(const LepIDs& ids) const {
  int N = GetN();

  for(int i = 0; i < N; i++)
    if(m_IDs[i] == ids)
      return true;

  return false;
}

bool LepIDsList::operator >= (const LepIDs& ids) const {
  return Contains(ids);
}

bool LepIDsList::Contains(const LepIDsList& ids) const {
  int N = ids.GetN();
  if(N > GetN())
    return false;

  for(int i = 0; i < N; i++)
    if(!Contains(ids[i]))
       return false;

  return true;
}
bool LepIDsList::operator >= (const LepIDsList& ids) const {
  return Contains(ids);
}

bool LepIDsList::operator <= (const LepIDsList& ids) const {
  int N = GetN();
  if(N > ids.GetN())
    return false;

  for(int i = 0; i < N; i++)
    if(!ids.Contains(Get(i)))
      return false;

  return true;
}

LepIDsList LepIDsList::Combinatorics(int Nlep) const {
  LepIDsList comb(Label());

  int N = GetN();
  for(int i = 0; i < N; i++)
    comb += Combinatorics(Nlep, i);
  
  return comb;
}

LepIDsList LepIDsList::Combinatorics(int Nlep, int index) const {
  LepIDsList comb(Label());

  LepIDs id = m_IDs[index];
  int N = id.GetN(); 

  while(N < Nlep){
    id += kGold;
    N++;
  }

  heapPermutation(id, N, N, comb);

  return comb;
}

void LepIDsList::heapPermutation(LepIDs& id, int size, int N, LepIDsList& comb) const {
  if(size == 1){
    comb += id;
    // int N = id.GetN();
    // cout << "Adding" << endl;
    // for(int i = 0; i < N; i++)
    //   cout << id[i] << " ";
    // cout << endl;
    // return;
  }

  for(int i = 0; i < size; i++){
    heapPermutation(id, size - 1, N, comb);

    if(size % 2 == 1)
      id.swap(0, size-1);
    else
      id.swap(i, size-1);
  }
}

///////////////////////////////////////////
////////// Lep class
///////////////////////////////////////////

Lep::Lep(LepFlavor flavor, LepCharge charge,
	 LepID id, LepSource source){
  m_Flavor = flavor;
  m_Charge = charge;
  m_ID = id;
  m_Source = source;
}

Lep::Lep(const Lep& lep){
  m_Flavor = lep.Flavor();
  m_Charge = lep.Charge();
  m_ID = lep.ID();
  m_Source = lep.Source();
}

Lep::~Lep(){ }

// for sorting leptons in category
bool Lep::operator < (const Lep& lep) const {
  if(m_Flavor < lep.Flavor())
    return true;
  if(m_Flavor > lep.Flavor())
    return false;

  if(m_ID < lep.ID())
    return true;
  if(m_ID > lep.ID())
    return false;

  if(m_Charge < lep.Charge())
    return true;
  if(m_Charge > lep.Charge())
    return false;
  
  if(m_Source < lep.Source())
    return true;

  return false;
  
}

// for matching lepton categories
bool Lep::operator == (const Lep& lep) const {
  return (m_Flavor == lep.Flavor() &&
	  m_ID     == lep.ID() &&
	  m_Charge == lep.Charge());
}

LepFlavor Lep::Flavor() const {
  return m_Flavor;
}

LepCharge Lep::Charge() const {
  return m_Charge;
}

LepID Lep::ID() const {
  return m_ID;
}

LepSource Lep::Source() const {
  return m_Source;
}

Lep& Lep::SetFlavor(LepFlavor flavor){
  m_Flavor = flavor;

  return *this;
}

Lep& Lep::SetCharge(LepCharge charge){
  m_Charge = charge;
  
  return *this;
}

Lep& Lep::SetID(LepID id){
  m_ID = id;
  
  return *this;
}

Lep& Lep::SetSource(LepSource source){
  m_Source = source;
  
  return *this;
}

std::string Lep::IDLabel() const {
  string ret = string(m_Flavor == LepFlavor::kElec ? "el" : "mu");
  if(m_ID == kGold)
    ret += "G";
  if(m_ID == kSilver)
    ret += "S";
  if(m_ID == kBronze)
    ret += "B";
  
  return ret;
}

LepSource GetLepSource(int PDGID, int genPDGID, int momPDGID){
  if(abs(genPDGID) != abs(PDGID)){ // wrong type
    return kFake;
    
  } else if(genPDGID != PDGID){ // charge flip
    return kSignal;
    
  } else if((abs(momPDGID) == abs(PDGID)) ||
	    (abs(momPDGID) == 23) ||
	    (abs(momPDGID) == 24) ||
	    ((abs(momPDGID) > 1000000) &&
	     (abs(momPDGID) < 3000000))){ // "signal" processes
    return kSignal;
    
  } else if(abs(momPDGID) == 15){ // Taus
    return kTau;
    
  } else if((abs(momPDGID) == 4) ||
	    ((abs(momPDGID) > 400) && (abs(momPDGID) < 500)) ||
	    ((abs(momPDGID) > 4000) && (abs(momPDGID) < 5000))){
    return kHFC;
    
  } else if((abs(momPDGID) == 5) ||
	    ((abs(momPDGID) > 500) && (abs(momPDGID) < 600)) ||
	    ((abs(momPDGID) > 5000) && (abs(momPDGID) < 6000))){
    return kHFB;
    
  } else {
    return kLF;
  }
}

///////////////////////////////////////////
////////// LepList class
///////////////////////////////////////////

LepList::LepList(){
  m_N = 0;
}

LepList::LepList(const LepList& leps){
  m_N = 0;
  int N = leps.GetN();
  for(int i = 0; i < N; i++)
    *this += leps[i];
}

LepList::~LepList(){
  Clear();
}

void LepList::Clear(){
  for(int i = 0; i < m_N; i++)
    delete m_Lep[i];
    
  m_Lep.clear();
  m_N = 0;
}

LepList& LepList::operator = (const LepList& lep){
  Clear();
  int N = lep.GetN();
  for(int i = 0; i < N; i++)
    *this += lep[i];
  
  return *this;
}

int LepList::GetN() const {
  return m_N;
}

const Lep& LepList::operator [] (int i) const {
  // no protection for array out-of-bounds
  return *m_Lep[i];
}

LepList& LepList::Add(const Lep& lep){
  return (*this) += lep;
}

LepList& LepList::operator += (const Lep& lep){
  m_Lep.push_back(new Lep(lep));
  m_N++;
  
  if(m_N == 1)
    return *this;
  
  // keep list ordered
  int index = m_N-2;
  while(index >= 0){
    if(*m_Lep[index+1] < *m_Lep[index]){
      Lep* temp = m_Lep[index+1];
      m_Lep[index+1] = m_Lep[index];
      m_Lep[index] = temp;
    } else {
      break;
    }
    index--;
  }

  return *this;
}

LepList& LepList::operator += (const LepList& leps){
  int N = leps.GetN();
  for(int i = 0; i < N; i++)
    *this += leps[i];
  
  return *this;
}

// for matching hemispheres
bool LepList::operator == (const LepList& leps) const {
  if(m_N != leps.GetN())
    return false;

  for(int i = 0; i < m_N; i++)
    if(!((*this)[i] == leps[i]))
      return false;

  return true;
}

// for sorting hemispheres, < means "more important"
bool LepList::operator < (const LepList& leps) const {
  if(m_N < leps.GetN())
    return false;
  if(m_N > leps.GetN())
    return true;

  for(int i = 0; i < m_N; i++){
    if((*this)[i] < leps[i])
      return false;
    if(leps[i] < (*this)[i])
      return true;
  }

  return false;
}

std::string LepList::GetIDLabel() const {
  if(m_N == 0)
    return "0";
  string label = "";
  for(int i = 0; i < m_N; i++)
    label += (*this)[i].IDLabel();

  return label;
}

VS LepList::GetFakeLabels(int Nfake) const {
  std::string label0 = "Fakes_";
  std::string label;
  vector<std::string> labels;
  
  for(int i = 0; i < m_N; i++){
    LepSource source = (*this)[i].Source();
    if(source > kTau){
      label = ((*this)[i].Flavor() == LepFlavor::kElec ? "elf" : "muf");
      if(source == kHFB || source == kHFC)
	label += "0";
      else
	label += "1";
      if( (*this)[i].ID() == kBronze )
	label += "B";
      else if( (*this)[i].ID() == kSilver)
	label += "S";
      else
	label += "G";
      
      labels.push_back(label);
    }
  }

  if(Nfake > 1)
    labels = Combinatorics(labels, Nfake);

  Nfake = labels.size();
  for(int i = 0; i < Nfake; i++)
    labels[i] = label0 + labels[i];
  
  return labels;
}

VS LepList::Combinatorics(const VS& labels, int N) const {
  if(N <= 1)
    return labels;

  VS labels_ret;
  VS labels_m1;
  int Nl = labels.size();
  
  if(Nl <= N){
    string lab = labels[0];
    for(int i = 1; i < Nl; i++)
      //lab += "_" + labels[i];
      lab += labels[i];
    labels_ret += lab;
    return labels_ret;
  }
    
  for(int i = 0; i < Nl; i++){
    labels_m1.clear();
    for(int j = i+1; j < Nl; j++)
      labels_m1 += labels[j];
    labels_m1 = Combinatorics(labels_m1, N-1);
    for(auto l : labels_m1)
      //labels_ret += labels[i] + "_" + l;
      labels_ret += labels[i] + l;
    }

  return labels_ret;
}

LepList LepList::GetFakes() const {
  LepList list;
  for(int i = 0; i < m_N; i++){
    if((*this)[i].Source() > kTau){
      Lep lep = (*this)[i];
      //lep.SetID(kGold);
      lep.SetCharge(kPos);
      list += lep;
    }
  }

  return list;
}

///////////////////////////////////////////
////////// Leptonic class
///////////////////////////////////////////

Leptonic::Leptonic(const LepList& Hem_a, const LepList& Hem_b, const std::string& label)
  : Criteria(label) {
  if(Hem_b < Hem_a){
    m_Lep_a = Hem_b;
    m_Lep_b = Hem_a;
  } else {
    m_Lep_a = Hem_a;
    m_Lep_b = Hem_b;
  }

   m_Nleptonic = 1;
   m_Leptonic.push_back(this);
}

Leptonic::Leptonic(const Leptonic& lep)
  : Criteria(lep) {
  m_Lep_a = lep.GetHemA();
  m_Lep_b = lep.GetHemB();
  int N = lep.GetN();
  m_Nleptonic = 1;
  m_Leptonic.push_back(this);
  
  for(int i = 1; i < N; i++)
    *this += lep[i];
}
    
Leptonic::~Leptonic(){
  for(int i = 1; i < m_Nleptonic; i++)
    delete m_Leptonic[i];
}

// see if there is a match between leptonic systems
bool Leptonic::operator == (const Criteria& val) const {
  const Leptonic* lep = dynamic_cast<const Leptonic*>(&val);
  if(!lep)
    return false;
  
  int N = lep->GetN();
  for(int i = 0; i < m_Nleptonic; i++)
    for(int j = 0; j < N; j++)
      if(((*this)[i].GetHemA() == (*lep)[j].GetHemA()) &&
	 ((*this)[i].GetHemB() == (*lep)[j].GetHemB()))
	return true;
  
  return false;
}

int Leptonic::GetN() const{
  return m_Nleptonic;
}

const Leptonic& Leptonic::operator [] (int i) const {
  // no check for array-out-of-bounds
  return *m_Leptonic[i];
}

Leptonic& Leptonic::operator += (const Leptonic& leptonic){
  if((m_Lep_a == leptonic.GetHemA()) &&
     (m_Lep_b == leptonic.GetHemB())){
    return *this;
  }
  m_Nleptonic++;
  m_Leptonic.push_back(new Leptonic(leptonic));
  return *this;
}

const LepList& Leptonic::GetHemA() const {
  return m_Lep_a;
}

const LepList& Leptonic::GetHemB() const {
  return m_Lep_b;
}

