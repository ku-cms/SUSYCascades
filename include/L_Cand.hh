#ifndef L_Cand_HH
#define L_Cand_HH

#include "ParticleList.hh"
#include "Leptonic.hh"
#include "RestFrames/RestFrames.hh"
using namespace RestFrames;

enum L_CandMatch { kMatched, kB, kUnmatched, kW, kZ };
enum L_CandSide { kAside, kBside };
 
class L_Cand {
  private:
    std::pair<ParticleList,RestFrameList> m_pair;
    void init(ParticleList PL, ConstRestFrameList RL);
    L_CandMatch m_Match;
    L_CandSide m_Side;
    TLorentzVector m_TLV;
    void init(ParticleList PL);
    LepFlavor m_Flav;
    ParticleList m_PL;
    const ConstRestFrameList m_RL;
  
  public:
    L_Cand();
    L_Cand(ParticleList PL);
    L_Cand(ParticleList PL, ConstRestFrameList RL);
    L_Cand(ConstRestFrameList RL, ParticleList PL);
    virtual ~L_Cand();
  
    const ParticleList PL();
    const ConstRestFrameList RL();
  
    L_CandMatch Match();
    void SetMatch(L_CandMatch match);
    LepFlavor Flavor();
    void SetFlavor(LepFlavor flav);
    LepID CandQual() const;
    void SetCandQual(LepID qual);
    const RestFrame& CandFrame(); // candidate's frame
    TLorentzVector TLV(int index = -1);
  
    Particle Cand_Part(int index);
    
    Particle operator [] (int index);
  
    double Pt();
    double Eta();
    double Phi();
    double M();
    double Mass();
    double P();
    double E();
    double ProngDeltaPhi();
    double DeltaPhi(const TLorentzVector& v);
    double DeltaPhi(const TVector3& v);
    double ProngDeltaEta();
    double ProngAbsDeltaEta();
    double ProngDeltaR();
    double ProngMassRatio(); // PMR
    double PMR();
    double Beta();
    double CosDecayAngle(const RestFrame& Frame = RestFrame::Empty());
  
};

inline void cand_matching(std::vector<L_Cand>& cand_list){ 
  int N_cands = cand_list.size();
  for(int i = 0; i < N_cands; i++){
    bool unmatched = true; // both jets are radiative
    bool matched = false; // both jets come from same boson
    //if((cand_list[i][0].MomPDGID() == 23) && cand_list[i][0].GenMomIndex() == cand_list[i][1].GenMomIndex()){
    if((cand_list[i][0].MomPDGID() == 23 && cand_list[i][1].MomPDGID() == 23)){
      cand_list[i].SetMatch(kMatched);
      unmatched = false;
      matched = true;
      break;
    }
    if(!matched){
      if(cand_list[i][0].MomPDGID() == 23 || cand_list[i][1].MomPDGID() == 23){
        cand_list[i].SetMatch(kZ);
        unmatched = false;
      }
      else if(cand_list[i][0].MomPDGID() == 24 || cand_list[i][1].MomPDGID() == 24){
        cand_list[i].SetMatch(kW);
        unmatched = false;
      }
      else if(cand_list[i][0].MomPDGID() == 6 || cand_list[i][1].MomPDGID() == 6){
        cand_list[i].SetMatch(kB);
        unmatched = false;
      }
    } // if(!matched)
    if(unmatched){
      cand_list[i].SetMatch(kUnmatched);
    }
  } // for(int i = 0; i < N_V_had; i++)
} // cand_matching(const std::vector<L_Cand>& cand_list)

#endif
