#ifndef Leptonic_HH
#define Leptonic_HH

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <iostream>
#include <vector>

#include "Criteria.hh"
#include "Particle.hh"

using std::vector;
using std::string;

enum LepFlavor { kElec, kMuon };
enum LepCharge { kPos, kNeg };
enum LepSource { kSignal, kTau, kHFB, kHFC, kLF, kFake };

class LepIDs {
public:
  LepIDs();
  virtual ~LepIDs();

  int GetN() const;
  LepID Get(int i) const;
  LepID operator [] (int i) const;
  
  void Add(LepID id);
  void operator += (LepID id);
  LepIDs& a(LepID id);

  bool IsEqual(const LepIDs& ids) const;
  bool operator == (const LepIDs& ids) const;

  void swap(int index_a, int index_b);
  
private:
  std::vector<LepID> m_IDs;
  
};

class LepIDsList {
public:
  LepIDsList(const string& label);
  virtual ~LepIDsList();

  string Label() const;
  
  void Add(const LepIDs& ids);
  void operator += (const LepIDs& ids);

  void Add(const LepIDsList& ids);
  void operator += (const LepIDsList& ids);

  int GetN() const;
  const LepIDs& Get(int i) const;
  const LepIDs& operator [] (int i) const;

  bool IsEqual(const LepIDsList& ids) const;
  bool operator == (const LepIDsList& ids) const;
  
  bool Contains(const LepIDs& ids) const;
  bool operator >= (const LepIDs& ids) const;

  bool Contains(const LepIDsList& ids) const;
  bool operator >= (const LepIDsList& ids) const;
  bool operator <= (const LepIDsList& ids) const;

  LepIDsList Combinatorics(int Nlep) const;
  
private:
  std::vector<LepIDs> m_IDs;
  string m_Label;

  LepIDsList Combinatorics(int Nlep, int index) const;
  void heapPermutation(LepIDs& id, int size, int N, LepIDsList& comb) const; 
};

///////////////////////////////////////////
////////// Lep class
///////////////////////////////////////////

class Lep {
public:
  Lep(LepFlavor flavor, LepCharge charge,
      LepID id = kGold, LepSource source = kSignal);
  Lep(const Lep&);
  virtual ~Lep();

  // for sorting leptons in category
  bool operator < (const Lep& lep) const;

  // for matching lepton categories
  bool operator == (const Lep& lep) const;

  LepFlavor Flavor() const;
  LepCharge Charge() const;
  LepID     ID() const;
  LepSource Source() const;

  Lep& SetFlavor(LepFlavor flavor);
  Lep& SetCharge(LepCharge charge);
  Lep& SetID(LepID id);
  Lep& SetSource(LepSource source);

  std::string IDLabel() const;
    
private:
  LepFlavor m_Flavor;
  LepCharge m_Charge;
  LepID     m_ID;
  LepSource m_Source;
};

LepSource GetLepSource(int PDGID, int genPDGID = 0, int momPDGID = 0);

///////////////////////////////////////////
////////// LepList class
///////////////////////////////////////////

class LepList {
public:
  LepList();
  LepList(const LepList& leps);
  virtual ~LepList();

  void Clear();
  
  LepList& operator = (const LepList& lep);
  
  int GetN() const;
  const Lep& operator [] (int i) const;

  LepList& Add(const Lep& lep);
  LepList& operator += (const Lep& lep);
  LepList& operator += (const LepList& leps);

  // for matching hemispheres
  bool operator == (const LepList& leps) const;

  // for sorting hemispheres
  bool operator < (const LepList& leps) const;

  std::string GetIDLabel() const;
  VS GetFakeLabels(int Nfake = 1) const;

  LepList GetFakes() const;
  
private:
  int m_N;
  vector<Lep*> m_Lep;

  VS Combinatorics(const VS& labels, int N) const;
};

///////////////////////////////////////////
////////// Leptonic class
///////////////////////////////////////////

class Leptonic : public Criteria {
public:
  Leptonic(const LepList& Hem_a, const LepList& Hem_b, const std::string& label = "none");
  Leptonic(const Leptonic&);
  
  virtual ~Leptonic();
  
  bool operator == (const Criteria& leptonic) const;

  int GetN() const;
  const Leptonic& operator [] (int i) const;
  Leptonic& operator += (const Leptonic& leptonic);

  const LepList& GetHemA() const;
  const LepList& GetHemB() const;
    
private:
  int m_Nleptonic;
  vector<Leptonic*> m_Leptonic;

  LepList m_Lep_a;
  LepList m_Lep_b;

};

#endif

