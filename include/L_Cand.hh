#ifndef L_Cand_HH
#define L_Cand_HH

#include "ParticleList.hh"
#include "Leptonic.hh"
#include "RestFrames/RestFrames.hh"
using namespace RestFrames;

enum L_CandMatch { kMatched, kB, kUnmatched, kW, kZ };
enum L_CandSide { kAside, kBside };
enum L_CandQF {kOSSF, kSSSF, kOSOF, kSSOF};
 
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
    //const ConstRestFrameList m_RL;
    bool m_hemi;
    int m_N = 0;
  
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
    void SetSameHemi(bool hemi);
    bool IsSameHemi();
  
    L_Cand(const L_Cand& other);
    L_Cand& operator=(const L_Cand& other);
    Particle Cand_Part(int index);
    Particle operator [] (int index);
    Particle Cand_PartPlus();
    Particle Cand_PartMinus();
    TLorentzVector TLV_Part(int index);
    TLorentzVector Cand_TLVPlus();
    TLorentzVector Cand_TLVMinus();
    int GetN();
  
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
    double CosDecayAngleRF(const RestFrame& Frame = RestFrame::Empty());
    double CosDecayAngle();
    double MCon(TVector3 boost = TVector3(0.,0.,0.)); // contravariant mass evaluated in frame set by boost
    double MCon(TLorentzVector frame = TLorentzVector(0.,0.,0.,0.)); // contravariant mass evaluated in frame
  
};

inline void cand_matching(L_Cand& cand){ 
    bool unmatched = true; // both leps are radiative
    bool matched = false; // both leps come from same boson
    if(cand[0].GenMomIndex() == cand[1].GenMomIndex() && cand[0].GenMomIndex() >= 0){
      cand.SetMatch(kMatched);
      unmatched = false;
      matched = true;
    }
    if(!matched){
      if(cand[0].MomPDGID() == 23 || cand[1].MomPDGID() == 23){
        cand.SetMatch(kZ);
        unmatched = false;
      }
      else if(cand[0].MomPDGID() == 24 || cand[1].MomPDGID() == 24){
        cand.SetMatch(kW);
        unmatched = false;
      }
      else if(cand[0].MomPDGID() == 6 || cand[1].MomPDGID() == 6){
        cand.SetMatch(kB);
        unmatched = false;
      }
    } // if(!matched)
    if(unmatched){
      cand.SetMatch(kUnmatched);
    }
} // cand_matching(L_Cand& cand)

inline void cand_matching(std::vector<L_Cand>& cand_list){ 
  int N_cands = cand_list.size();
  for(int i = 0; i < N_cands; i++){
    cand_matching(cand_list[i]);
  } // for(int i = 0; i < N_cands; i++)
} // cand_matching(const std::vector<L_Cand>& cand_list)

#endif
