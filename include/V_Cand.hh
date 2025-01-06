#ifndef V_Cand_HH
#define V_Cand_HH

#include "ParticleList.hh"
#include "RestFrames/RestFrames.hh"
using namespace RestFrames;

enum V_CandMatch { kMatched, kW, kZ, kB, kUnmatched };
enum V_CandSide { kAside, kBside };
enum V_CandType { kSib, kAunt, kLep }; // kSib = 2prong, kAunt = Aunt and Niece pair, kLep = Aunt and Niece pair and sib is lep
  
class V_Cand {
private:
  std::pair<ParticleList,RestFrameList> m_pair;
  void init(ParticleList PL, ConstRestFrameList RL);
  V_CandMatch m_Match;
  V_CandSide m_Side;
  V_CandType m_Type;
  TLorentzVector m_TLV;

public:
  V_Cand();
  V_Cand(ParticleList PL, ConstRestFrameList RL);
  V_Cand(ConstRestFrameList RL, ParticleList PL);
  virtual ~V_Cand();

  const ParticleList PL();
  const ConstRestFrameList RL();

  V_CandMatch Match();
  void SetMatch(V_CandMatch match);
  V_CandSide Side();
  void SetSide(V_CandSide side);
  V_CandType Type();
  void SetType(V_CandType type);

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

  const RestFrame& CandFrame(); // candidate's frame
  bool IsProng(const RestFrame& frame);

};

  

#endif
