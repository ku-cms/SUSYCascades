#include "V_Cand.hh"

V_Cand::V_Cand() {}

V_Cand::V_Cand(ParticleList PL, ConstRestFrameList RL){
  init(PL,RL);
}

V_Cand::V_Cand(ConstRestFrameList RL, ParticleList PL){
  init(PL,RL);
}

V_Cand::~V_Cand() {}
  
void V_Cand::init(ParticleList PL, ConstRestFrameList RL){
  if(PL.size() != RL.GetN()){
    std::cout << "Can't Make Candidate with different sized lists! \n ParticleList size: " << PL.size() << " \n RestFrameList size: " << RL.GetN() << std::endl;
    return;
  }
  m_pair = std::make_pair(PL,RL);
  m_TLV.SetPtEtaPhiM(0.,0.,0.,0.);
  for(int i = 0; i < int(PL.size()); i++){
    TLorentzVector dummy_TLV;
    dummy_TLV.SetPtEtaPhiM(PL[i].Pt(),PL[i].Eta(),PL[i].Phi(),PL[i].M());
    m_TLV += dummy_TLV;
  }
}

const ParticleList V_Cand::PL(){
  return m_pair.first;
}

const ConstRestFrameList V_Cand::RL(){
  return m_pair.second;
}

Particle V_Cand::Cand_Part(int index){
  return PL()[index];
}

V_CandMatch V_Cand::Match(){
  return m_Match;
}

void V_Cand::SetMatch(V_CandMatch match){
  m_Match = match;
}

V_CandSide V_Cand::Side(){
  return m_Side;
}

void V_Cand::SetSide(V_CandSide side){
  m_Side = side;
}

V_CandType V_Cand::Type(){
  return m_Type;
}

void V_Cand::SetType(V_CandType type){
  m_Type = type;
}

Particle V_Cand::operator[](int index){
  return Cand_Part(index);
}

int V_Cand::size(){
  return PL().size();
}

int V_Cand::prongs(){
  return size();
}

double V_Cand::Pt(){
  return m_TLV.Pt();
}

double V_Cand::Eta(){
  return m_TLV.Eta();
}

double V_Cand::Phi(){
  return m_TLV.Phi();
}

double V_Cand::M(){
  return m_TLV.M();
}

double V_Cand::Mass(){
  return M();
}

double V_Cand::P(){
  return m_TLV.P();
}

double V_Cand::E(){
  return m_TLV.E();
}

double V_Cand::ProngDeltaPhi(){ 
  if(prongs() < 1) return 0.;
  else if(prongs() == 1) return Phi();
  double maxDeltaPhi = PL()[0].DeltaPhi(PL()[1]);
  if(prongs() == 2) return maxDeltaPhi;
  for(int i = 1; i < prongs(); i++){
    double DeltaPhi = PL()[i].DeltaPhi(PL()[i-1]);
    if(DeltaPhi > maxDeltaPhi)
      maxDeltaPhi = DeltaPhi;
  }
  return maxDeltaPhi;
}

double V_Cand::ProngDeltaEta(){ 
  if(prongs() < 1) return 0.;
  else if(prongs() == 1) return Eta();
  double maxDeltaEta = PL()[0].Eta()-PL()[1].Eta();
  if(prongs() == 2) return maxDeltaEta;
  for(int i = 1; i < prongs(); i++){
    double DeltaEta = PL()[i].Eta()-PL()[i-1].Eta();
    if(DeltaEta > maxDeltaEta)
      maxDeltaEta = DeltaEta;
  }
  return maxDeltaEta;
}

double V_Cand::ProngAbsDeltaEta(){ 
  if(prongs() < 1) return 0.;
  else if(prongs() == 1) return Eta();
  double maxDeltaEta = fabs(PL()[0].Eta())-fabs(PL()[1].Eta());
  if(prongs() == 2) return maxDeltaEta;
  for(int i = 1; i < prongs(); i++){
    double DeltaEta = fabs(PL()[i].Eta())-fabs(PL()[i-1].Eta());
    if(DeltaEta > maxDeltaEta)
      maxDeltaEta = DeltaEta;
  }
  return maxDeltaEta;
}

double V_Cand::ProngDeltaR(){ 
  if(prongs() < 2) return 0.;
  double maxDeltaR = PL()[0].DeltaR(PL()[1]);
  for(int i = 1; i < prongs(); i++){
    double DeltaR = PL()[i].DeltaR(PL()[i-1]);
    if(DeltaR > maxDeltaR)
      maxDeltaR = DeltaR;
  }
  return maxDeltaR;
}

double V_Cand::ProngMassRatio(){ 
  if(prongs() != 2) return -1.;
  if(PL()[1].M() < 1.e-5) return -1.;
  double prongMassRatio = PL()[0].M()/PL()[1].M();
  return prongMassRatio;
}

double V_Cand::PMR(){
  return ProngMassRatio();
}

const RestFrame& V_Cand::CandFrame(){
  const RestFrame& LAB = RL().Get(0).GetLabFrame();
  ConstRestFrameList children = LAB.GetListVisibleFrames();
  int minDepth = LAB.GetFrameDepth(RL().Get(0));
  int minIndex = 0;
  for(int i = 1; i < RL().GetN(); i++){
    int newDepth = LAB.GetFrameDepth(RL().Get(i));
    if(newDepth < minDepth){
      minDepth = LAB.GetFrameDepth(RL().Get(i));
      minIndex = i;
    }
  }
  return RL().Get(minIndex).GetProductionFrame();
}

double V_Cand::CosDecayAngle(const RestFrame& Frame){
  return CandFrame().GetCosDecayAngle(Frame);
}

bool V_Cand::IsProng(const RestFrame& frame){
  for(int i = 0; i < prongs(); i++)
    if(frame == RL().Get(i)) return true;
  return false;
}

