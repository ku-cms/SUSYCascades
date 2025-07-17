#include "L_Cand.hh"

L_Cand::L_Cand() {}

L_Cand::L_Cand(ParticleList PL){
  init(PL);
}

L_Cand::L_Cand(ParticleList PL, ConstRestFrameList RL){
  init(PL,RL);
}

L_Cand::L_Cand(ConstRestFrameList RL, ParticleList PL){
  init(PL,RL);
}

L_Cand::~L_Cand() {}
  
void L_Cand::init(ParticleList PL){
  m_hemi = false;
  m_PL = PL;
  m_TLV.SetPtEtaPhiM(0.,0.,0.,0.);
  for(int i = 0; i < int(PL.size()); i++){
    TLorentzVector dummy_TLV;
    dummy_TLV.SetPtEtaPhiM(PL[i].Pt(),PL[i].Eta(),PL[i].Phi(),PL[i].M());
    m_TLV += dummy_TLV;
  }
}

void L_Cand::init(ParticleList PL, ConstRestFrameList RL){
  if(PL.size() != RL.GetN()){
    std::cout << "Can't Make Candidate with different sized lists! \n ParticleList size: " << PL.size() << " \n RestFrameList size: " << RL.GetN() << std::endl;
    return;
  }
  init(PL);
  m_pair = std::make_pair(PL,RL);
}

const ParticleList L_Cand::PL(){
  return m_PL;
}

const ConstRestFrameList L_Cand::RL(){
  return m_pair.second;
}

TLorentzVector L_Cand::TLV(int index){
  if(index < 0)
    return m_TLV;
  else return PL()[index];
}

L_CandMatch L_Cand::Match(){
  return m_Match;
}

void L_Cand::SetMatch(L_CandMatch match){
  m_Match = match;
}

LepFlavor L_Cand::Flavor(){
  return m_Flav;
}

void L_Cand::SetFlavor(LepFlavor flav){
  m_Flav = flav;
}

void L_Cand::SetSameHemi(bool hemi){
  m_hemi = hemi;
}

bool L_Cand::IsSameHemi(){
  return m_hemi;
}

Particle L_Cand::Cand_Part(int index){
  return PL()[index];
}

Particle L_Cand::operator[](int index){
  return Cand_Part(index);
}

Particle L_Cand::Cand_PartPlus(){
  if(Cand_Part(0).Charge() > 0) return Cand_Part(0);
  else return Cand_Part(1);
}

Particle L_Cand::Cand_PartMinus(){
  if(Cand_Part(0).Charge() < 0) return Cand_Part(0);
  else return Cand_Part(1);
}

TLorentzVector L_Cand::TLV_Part(int index){
  return TLorentzVector(PL()[index].X(), PL()[index].Y(), PL()[index].Z(), PL()[index].T());
}

TLorentzVector L_Cand::Cand_TLVPlus(){
  if(Cand_Part(0).Charge() > 0) return Cand_Part(0);
  else return Cand_Part(1);
}

TLorentzVector L_Cand::Cand_TLVMinus(){
  if(Cand_Part(0).Charge() < 0) return Cand_Part(0);
  else return Cand_Part(1);
}

L_Cand::L_Cand(const L_Cand& other) {
    // Copy the basic members directly
    m_Match = other.m_Match;
    m_Side = other.m_Side;
    m_TLV = other.m_TLV;
    m_Flav = other.m_Flav;
    m_hemi = other.m_hemi;

    // Copy the ParticleList
    m_PL = other.m_PL;

    // Now call init() to recalculate m_TLV based on the copied ParticleList
    init(m_PL);

    m_pair = other.m_pair;  // If m_pair contains RestFrameList, copy it
}

L_Cand& L_Cand::operator=(const L_Cand& other){
  if (this != &other) {  // Avoid self-assignment
    // Clean up any resources if needed, or reset some members
    m_Match = other.m_Match;
    m_Side = other.m_Side;
    m_TLV = other.m_TLV;
    m_Flav = other.m_Flav;
    m_hemi = other.m_hemi;

    // Call init() to ensure correct initialization of m_PL and m_TLV
    m_PL = other.m_PL; // Copy the ParticleList (shallow copy by default)
    init(m_PL);  // Re-initialize m_TLV based on the new ParticleList

    // If RestFrameList is also part of the copy and needs to be handled:
    if (other.m_pair.second.GetN() > 0) {
      // If you need to copy RestFrameList too, call the second init overload:
      m_pair = other.m_pair;
      init(m_PL, m_pair.second);  // Initialize with both ParticleList and RestFrameList
    }
  }
  return *this; 
}

double L_Cand::Pt(){
  return m_TLV.Pt();
}

double L_Cand::Eta(){
  return m_TLV.Eta();
}

double L_Cand::Phi(){
  return m_TLV.Phi();
}

double L_Cand::M(){
  return m_TLV.M();
}

double L_Cand::Mass(){
  return M();
}

double L_Cand::P(){
  return m_TLV.P();
}

double L_Cand::E(){
  return m_TLV.E();
}

double L_Cand::ProngDeltaPhi(){ 
  return PL()[0].DeltaPhi(PL()[1]);
}

double L_Cand::DeltaPhi(const TLorentzVector& v){
  return m_TLV.DeltaPhi(v);
}

double L_Cand::DeltaPhi(const TVector3& v){
  return m_TLV.Vect().DeltaPhi(v);
}

double L_Cand::ProngDeltaEta(){ 
  return PL()[0].Eta()-PL()[1].Eta();
}

double L_Cand::ProngAbsDeltaEta(){ 
  return fabs(PL()[0].Eta())-fabs(PL()[1].Eta());
}

double L_Cand::ProngDeltaR(){ 
  return PL()[0].DeltaR(PL()[1]);
}

double L_Cand::ProngMassRatio(){ 
  return PL()[0].M()/PL()[1].M();
}

double L_Cand::PMR(){
  return ProngMassRatio();
}

double L_Cand::Beta(){ 
  return m_TLV.P()/m_TLV.E();
}

const RestFrame& L_Cand::CandFrame(){
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

double L_Cand::CosDecayAngleRF(const RestFrame& Frame){
  return CandFrame().GetCosDecayAngle(Frame);
}

double L_Cand::CosDecayAngle(){
  TVector3 boost = m_TLV.BoostVector();
  Particle part_Cand_Child = Cand_PartPlus();
  TLorentzVector TLV_Cand_Child;
  TLV_Cand_Child.SetPtEtaPhiM(part_Cand_Child.Pt(), part_Cand_Child.Eta(), part_Cand_Child.Phi(), part_Cand_Child.M());
  TLV_Cand_Child.Boost(-boost);
  return TLV_Cand_Child.Vect().Unit().Dot(boost.Unit());
}

double L_Cand::MCon(TVector3 boost){
  TLorentzVector TLV_Cand_Child_a = TLV_Part(0);
  TLorentzVector TLV_Cand_Child_b = TLV_Part(1);
  TLV_Cand_Child_a.Boost(boost);
  TLV_Cand_Child_b.Boost(boost);
  return sqrt(TLV_Cand_Child_a.E()*TLV_Cand_Child_b.E() + TLV_Cand_Child_a.Vect().Dot(TLV_Cand_Child_b.Vect()) + TLV_Cand_Child_a.Mag2() + TLV_Cand_Child_b.Mag2());
}

double L_Cand::MCon(TLorentzVector frame){
  TVector3 boost = -frame.BoostVector();
  return MCon(boost);
}

