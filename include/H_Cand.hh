#ifndef H_Cand_HH
#define H_Cand_HH

#include "ParticleList.hh"
#include "RestFrames/RestFrames.hh"
using namespace RestFrames;

enum H_CandMatch { kMatched, kW, kZ, kB, kUnmatched };
enum H_CandSide { kAside, kBside };
enum H_CandType { kSib, kAunt, kLep }; // kSib = 2prong, kAunt = Aunt and Niece pair, kLep = Aunt and Niece pair and sib is lep
  
class H_Cand {
private:
  std::pair<ParticleList,RestFrameList> m_pair;
  void init(ParticleList PL, ConstRestFrameList RL);
  H_CandMatch m_Match;
  H_CandSide m_Side;
  H_CandType m_Type;
  TLorentzVector m_TLV;

public:
  H_Cand();
  H_Cand(ParticleList PL, ConstRestFrameList RL);
  H_Cand(ConstRestFrameList RL, ParticleList PL);
  virtual ~H_Cand();

  const ParticleList PL();
  const ConstRestFrameList RL();

  H_CandMatch Match();
  void SetMatch(H_CandMatch match);
  H_CandSide Side();
  void SetSide(H_CandSide side);
  H_CandType Type();
  void SetType(H_CandType type);

  Particle Cand_Part(int index);
  
  Particle operator [] (int index);

  int size();
  int prongs();

  double Pt();
  double Eta();
  double Phi();
  double M();
  double Mass();
  double P();
  double E();
  double CosDecayAngle(const RestFrame& Frame = RestFrame::Empty());
  double ProngDeltaPhi();
  double ProngDeltaEta();
  double ProngAbsDeltaEta();
  double ProngDeltaR();
  double ProngMassRatio(); // PMR
  double PMR();
  void cand_matching(std::vector<H_Cand>& cand_list); 

  const RestFrame& CandFrame(); // candidate's frame
  bool IsProng(const RestFrame& frame);

};

  

#endif
