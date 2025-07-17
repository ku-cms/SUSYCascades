#ifndef Particle_HH
#define Particle_HH

#include "TLorentzVector.h"

class ParticleList;

/// Particle ID level
enum ParticleIDType { kNothing, kVeryLoose, kLoose, kMedium, kTight, kVeryTight, kVeryVeryTight };
enum LepID { kGold, kSilver, kBronze };

class Particle : public TLorentzVector {
public:
  Particle();
    
  virtual ~Particle();

  int Charge() const;
  void SetCharge(int charge);

  int PDGID() const;
  void SetPDGID(int pdgid);

  int MomPDGID() const;
  void SetMomPDGID(int pdgid);

  int SourceID() const;
  void SetSourceID(int sourceid);

  int Index() const;
  void SetIndex(int index);

  int GenIndex() const;
  void SetGenIndex(int genindex);

  int GenMomIndex() const;
  void SetGenMomIndex(int genmomindex);

  ParticleIDType ParticleID() const;
  void SetParticleID(ParticleIDType id);

  double RelIso() const;
  double MiniIso() const;
  void SetRelIso(double iso);
  void SetMiniIso(double iso);

  double Dxy() const;
  double DxyErr() const;
  double Dz() const;
  double DzErr() const;
  double IP3D() const;
  double SIP3D() const;
  void SetDxy(double val);
  void SetDxyErr(double val);
  void SetDz(double val);
  void SetDzErr(double val);
  void SetIP3D(double val);
  void SetSIP3D(double val);

  double D3d() const;
  double D3dSig() const;
  double CosTheta() const;
  double Ndof() const;
  double ProbB() const;
  double ProbC() const;
  void SetD3d(double val);
  void SetD3dSig(double val);
  void SetCosTheta(double val);
  void SetNdof(double val);  
  void SetProbB(double val);
  void SetProbC(double val);

  double Btag() const;
  void SetBtag(double btag);

  ParticleIDType BtagID() const;
  void SetBtagID(ParticleIDType id);

  LepID LepQual() const;
  void SetLepQual(LepID qual);

  bool IsLowPt() const;
  void SetIsLowPt(bool is_low_pt);

  int TightCharge() const;
  void SetTightCharge(int tightCharge);

  Particle Merge(const Particle&) const;
  
  operator ParticleList() const;
  ParticleList operator + (const Particle& part) const; 
  ParticleList operator + (const ParticleList& parts) const; 
    
private:
  int m_Charge;
  int m_PDGID;
  int m_MomPDGID;
  int m_SourceID;
  int m_Index;
  int m_GenIndex; // index of particle in gen record, for reco objects, corresponds to matched gen particle index
  int m_GenMomIndex; // index of mom of particle in gen record, for reco objects, corresponds to mom of matched gen particle index
  ParticleIDType m_ParticleID;
  
  double m_Btag;
  ParticleIDType m_BtagID;

  // lepton stuff
  LepID m_LepQual;
  bool m_IsLowPt;
  int m_TightCharge;
  
  double m_RelIso;
  double m_MiniIso;

  double m_Dxy;
  double m_DxyErr;
  double m_Dz;
  double m_DzErr;
  double m_IP3D;
  double m_SIP3D;

  double m_D3d;
  double m_D3dSig;
  double m_CosTheta;
  double m_Ndof;
  double m_ProbB; 
  double m_ProbC;
  
};

bool sortbypt(const Particle& p1, const Particle& p2);

template <typename V>
inline bool inVec(const std::vector<V>& vect, const V& value){
  return std::find(vect.begin(), vect.end(), value) != vect.end();
}

#endif
