#include "H_Cand.hh"

H_Cand::H_Cand() {}

H_Cand::H_Cand(ParticleList PL, ConstRestFrameList RL){
  init(PL,RL);
}

H_Cand::H_Cand(ConstRestFrameList RL, ParticleList PL){
  init(PL,RL);
}

H_Cand::~H_Cand() {}
  
void H_Cand::init(ParticleList PL, ConstRestFrameList RL){
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

const ParticleList H_Cand::PL(){
  return m_pair.first;
}

const ConstRestFrameList H_Cand::RL(){
  return m_pair.second;
}

Particle H_Cand::Cand_Part(int index){
  return PL()[index];
}

H_CandMatch H_Cand::Match(){
  return m_Match;
}

void H_Cand::SetMatch(H_CandMatch match){
  m_Match = match;
}

H_CandSide H_Cand::Side(){
  return m_Side;
}

void H_Cand::SetSide(H_CandSide side){
  m_Side = side;
}

H_CandType H_Cand::Type(){
  return m_Type;
}

void H_Cand::SetType(H_CandType type){
  m_Type = type;
}

Particle H_Cand::operator[](int index){
  return Cand_Part(index);
}

int H_Cand::size(){
  return PL().size();
}

int H_Cand::prongs(){
  return size();
}

double H_Cand::Pt(){
  return m_TLV.Pt();
}

double H_Cand::Eta(){
  return m_TLV.Eta();
}

double H_Cand::Phi(){
  return m_TLV.Phi();
}

double H_Cand::M(){
  return m_TLV.M();
}

double H_Cand::Mass(){
  return M();
}

double H_Cand::P(){
  return m_TLV.P();
}

double H_Cand::E(){
  return m_TLV.E();
}

double H_Cand::ProngDeltaPhi(){ 
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

double H_Cand::ProngDeltaEta(){ 
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

double H_Cand::ProngAbsDeltaEta(){ 
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

double H_Cand::ProngDeltaR(){ 
  if(prongs() < 2) return 0.;
  double maxDeltaR = PL()[0].DeltaR(PL()[1]);
  for(int i = 1; i < prongs(); i++){
    double DeltaR = PL()[i].DeltaR(PL()[i-1]);
    if(DeltaR > maxDeltaR)
      maxDeltaR = DeltaR;
  }
  return maxDeltaR;
}

double H_Cand::ProngMassRatio(){ 
  if(prongs() != 2) return -1.;
  if(PL()[1].M() < 1.e-5) return -1.;
  double prongMassRatio = PL()[0].M()/PL()[1].M();
  return prongMassRatio;
}

double H_Cand::PMR(){
  return ProngMassRatio();
}

const RestFrame& H_Cand::CandFrame(){
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

double H_Cand::CosDecayAngle(const RestFrame& Frame){
  return CandFrame().GetCosDecayAngle(Frame);
}

bool H_Cand::IsProng(const RestFrame& frame){
  for(int i = 0; i < prongs(); i++)
    if(frame == RL().Get(i)) return true;
  return false;
}

void H_Cand::cand_matching(std::vector<H_Cand>& cand_list){ 
  int N_cands = cand_list.size();
    for(int i = 0; i < N_cands; i++){
      bool unmatched = true; // both jets are radiative
      bool matched = false; // both jets come from same boson
      for(int j = 0; j < int(cand_list[i].size()); j++){
        for(int k = j+1; k < int(cand_list[i].size()); k++){
          if((abs(cand_list[i][j].MomPDGID()) == 23 || abs(cand_list[i][j].MomPDGID()) == 24) && cand_list[i][j].GenMomIndex() == cand_list[i][k].GenMomIndex()){
            cand_list[i].SetMatch(kMatched);
            unmatched = false;
            matched = true;
            break;
          }
        }
        if(!matched){
          if(cand_list[i][j].MomPDGID() == 24){
            cand_list[i].SetMatch(kW);
            unmatched = false;
          }
          else if(cand_list[i][j].MomPDGID() == 23){
            cand_list[i].SetMatch(kZ);
            unmatched = false;
          }
          else if(cand_list[i][j].MomPDGID() == 6){
            cand_list[i].SetMatch(kB);
            unmatched = false;
          }
        } // if(!matched)
      } // for(int j = 0; j < int(cand_list[i].size()); j++)
      if(unmatched){
        cand_list[i].SetMatch(kUnmatched);
      }
    } // for(int i = 0; i < N_V_had; i++)
} // cand_matching(const std::vector<H_Cand>& cand_list)
