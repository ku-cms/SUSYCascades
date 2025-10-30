#ifndef AnalysisBase_h
#define AnalysisBase_h

#include <iostream>
#include <filesystem>
#include <regex>
#include <stdexcept>
#include <optional>
#include <algorithm>
#include <cctype>
#include <variant>
#include <cmath>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <string>
#include <functional>

#include "NeventTool.hh"
#include "XsecTool.hh"
#include "JSONTool.hh"
#include "PUTool.hh"
#include "LHETool.hh"
#include "BtagSFTool.hh"
#include "LepSFTool.hh"
#include "JMETool.hh"
#include "METTriggerTool.hh"
#include "PrefireTool.hh"

#include "Particle.hh"
#include "Systematics.hh"
#include "CascadesTreeEncoder.hh"
#include "L_Cand.hh"

#include "correction.h"
#include "mt2_bisect.hh"

using namespace std;
namespace fs = std::filesystem;

class ParticleList;

template <class Base>
class AnalysisBase : public Base {

public:
  AnalysisBase(TTree* tree = 0);
  virtual ~AnalysisBase();

  void AddLabels(const string& dataset, const string& filetag);
  void AddEventCountFile(const string& rootfile);
  void AddFilterEffFile(const string& rootfile);
  void AddJSONFile(const string& jsonfile);
  void AddPUFolder(const string& pufold);
  void AddBtagFolder(const string& btagfold, const string& proc_rootfile="", int year=1);
  void AddLepFolder(const string& lepfold);
  void AddJMEFolder(const string& jmefold);
  void AddMETTriggerFile(const string& csvfile);
  void AddPrefireFile(const string& prefirefile);
  void AddXSecJSON(const string& XSjsonfile);
  #ifdef _CMSSW_
  void AddLHAPDF();
  #endif
  void DoSMS();
  void DoData();
  void DoFastSim();
  void DoPrivateMC();
  void DoCascades();
  void AddSystematics();
  void AddJESSystematics();
  void AddJERSystematics();
  void AddMETSystematics();
  void AddEESSystematics();
  void AddMMSSystematics();
  
  void InitializeHistograms(vector<TH1D*>& histos);
  void BookHistograms(vector<TH1D*>& histos);

  string GetEntry(int entry);

  // event functions
  virtual int GetRunNum();
  virtual int GetLumiNum();
  virtual long GetEventNum();
  virtual int GetNpartons();

  virtual bool PassEventFilter();
  virtual bool FastSimEventVeto(const ParticleList& GenJets);
  virtual double EGvalue(int jetIndex, int updown);
  virtual double GetPrefireWeight(int updown);
  virtual bool PassesJVM(const ParticleList& jets) const;
  
  // analysis functions
  virtual int GetNPV();
  virtual int GetNPUtrue();

  virtual bool GetMETtrigger();
  virtual bool GetMETORtrigger();
  virtual bool GetMETDoubleMutrigger();
  virtual bool GetSingleElectrontrigger();
  virtual bool GetSingleMuontrigger();
  virtual bool GetDoubleElectrontrigger();
  virtual bool GetDoubleMuontrigger();
  virtual bool GetTripleElectrontrigger();
  virtual bool GetTripleMuonLowPTtrigger();
  virtual bool GetTripleMuonHighPTtrigger();
  virtual bool GetDiMuEleLowPTtrigger();
  virtual bool GetDiMuEleHighPTtrigger();
  virtual bool GetDiEleMutrigger();
  virtual bool GetEMutrigger(); 
  virtual bool GetEMuMutrigger(); // High MuPt Low ElePt
  virtual bool GetEMuEtrigger(); // Low MuPt High ElePt
  
  virtual void ApplyMETPhiCorrections(TLorentzVector &metOut, int runNumber, int npvGood, const std::string &variation);
  virtual TVector3 GetPV(bool& good);
  virtual TVector3 GetMET();
  virtual ParticleList GetJets(int id = -1);
  virtual ParticleList GetJetsMET(TVector3& MET, int id = -1);
  virtual ParticleList GetGenJets();
  virtual ParticleList GetElectrons();
  virtual ParticleList GetLowPtElectrons();
  virtual ParticleList GetMuons();

  virtual TVector3 GetAltMET();
  
  virtual double Get_LHE_HT();
  virtual double Get_LHE_HTIncoming();

  virtual bool IsHEM(Particle part);

  virtual TVector3 GetGenMET();
  virtual ParticleList GetGenElectrons();
  virtual ParticleList GetGenMuons();
  virtual ParticleList GetGenNeutrinos();
  virtual ParticleList GetGenBosons();
  virtual ParticleList GetGenSparticles();
 
  double DeltaPhiMin(const vector<TLorentzVector>& JETs, const TVector3& MET, int N = -1);
  double DeltaPhiMin(const vector<pair<TLorentzVector, bool> >& JETs, const TVector3& MET, int N = -1);
  
  void MomTensorCalc(vector<TLorentzVector>& input, vector<double>& eigenvalues, double pow = 1., bool threeD = true);

  virtual std::pair<int,int> GetSUSYMasses();
  virtual int GetGenMass(const int& u_PDGID);
  virtual int GetGenSUSYNBosons(const int& u_PDGID);

  virtual uint16_t GetGenCascadesTree();
  virtual int GetGenCascadesProduction();
  virtual int GetGenCascadesProduction(int& firstSpart, int& secondSpart);
  virtual std::pair<int,int> GetGenCascadesDecayMode(const int& GenIndex);
  virtual int GetGenCascadesIndex(const int& u_momIndex, const int& u_PDGID);
  virtual vector<int> GetLSPParents();

  bool IsSMS(){ return m_IsSMS; }
  bool IsData(){ return m_IsData; }
  bool IsPrivateMC(){ return m_IsPrivateMC; }
  bool IsCascades(){ return m_IsCascades; }
  bool IsFastSim(){ return m_IsFastSim; }
  
  string GetDataSet(){ return m_DataSet; }
  string GetFileTag(){ return m_FileTag; }
  int    GetYear(){ return m_year; }
  
protected:
  bool m_IsSMS;
  bool m_IsData;
  bool m_IsFastSim;
  bool m_IsPrivateMC;
  bool m_IsCascades;
  
  virtual double GetEventWeight();
  virtual double GetGenEventWeight();
  virtual double GetPUWeight(int updown = 0);
  virtual double GetPDFWeight(int updown = 0);
  virtual double GetMuFWeight(int updown = 0);
  virtual double GetMuRWeight(int updown = 0);
  virtual double GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown = 0, ParticleIDType tag = kMedium);

  virtual double GetElIDSFWeight(const ParticleList& els, int updown = 0);
  virtual double GetElISOSFWeight(const ParticleList& els, int updown = 0);
  virtual double GetElSIPSFWeight(const ParticleList& els, int updown = 0);
  virtual double GetElVLIDSFWeight(const ParticleList& els, int updown = 0);
  virtual double GetMuIDSFWeight(const ParticleList& mus, int updown = 0);
  virtual double GetMuISOSFWeight(const ParticleList& mus, int updown = 0);
  virtual double GetMuSIPSFWeight(const ParticleList& mus, int updown = 0);
  virtual double GetMuVLIDSFWeight(const ParticleList& mus, int updown = 0);
 
  virtual double GetMETTriggerSFWeight(double MET, double HT, int Nele, int Nmu, int updown = 0);
  virtual int GetMETTriggerSFCurve(double HT, int Nele, int Nmu);
  virtual double GetXsec();
  virtual double GetFilterEff();
  virtual double GetNevent();
  virtual double GetNweight();
  virtual bool   IsGoodEvent();

  void SetSystematic(const Systematic& sys);

  string m_DataSet;
  string m_FileTag;
  int m_year;
  bool m_IsAPV = false;
  bool m_IsUL = false;
  bool m_IsEE = false;
  bool m_IsBPix = false;
  bool m_IsRun3 = false;

  Systematics m_Systematics;
  std::string m_JMEYearOverride{};
  std::string m_JMEEraOverride{};
  
private:

  NeventTool      m_NeventTool;
  XsecTool        m_XsecTool;
  JSONTool        m_JSONTool;
  PUTool          m_PUTool;
  LHETool         m_LHETool;
  BtagSFTool      m_BtagSFTool;
  LepSFTool       m_LepSFTool;
  JMETool         m_JMETool;
  SystematicsTool m_SysTool;
  METTriggerTool  m_METTriggerTool;
  PrefireTool     m_PrefireTool;

  int m_SampleIndex;
  virtual int GetSampleIndex();
  int m_Nsample;
  std::map<long long,int>   m_HashToIndex;
  std::map<int,std::string> m_IndexToSample;
  std::map<int,double>      m_IndexToXsec;
  std::map<int,double>      m_IndexToNevent;
  std::map<int,double>      m_IndexToNweight;

  const Systematic* m_CurSys;
  const Systematic& CurrentSystematic() const;
  std::unique_ptr<correction::CorrectionSet> m_cset_Btag;
  std::unique_ptr<correction::CorrectionSet> m_cset_METPhi;
  double m_BtagLooseWP;
  double m_BtagMediumWP;
  double m_BtagTightWP;
  double m_BtagVeryTightWP;
  double m_BtagVeryVeryTightWP;
  virtual void clip_string(string& str, const string& clip){
    size_t pos = str.find(clip);
    if (pos != std::string::npos)
      str.erase(pos, clip.length());
  }
  virtual void wrapPhi(double& phi) {
    while (phi > M_PI)  phi -= 2.0*M_PI;
    while (phi <= -M_PI) phi += 2.0*M_PI;
  }
  bool minus_iso_hoe(int WPBitMap, int value, std::function<bool(int, int)> comp);
  int extract_nano_version(const std::string& dirname);
  std::string normalize_tag(const std::string& filetag);
  std::string normalize_filetag(const std::string& tag);
  std::string find_clib_file(const std::string& fold, const std::string& filename);

  std::string getJMEYearKey() const;
  std::string m_jmeYearKey;
  void setJMEYearOverride(const std::string& s) { m_JMEYearOverride = s; }
  void setJMEEraOverride(const std::string& s)  { m_JMEEraOverride = s;  }
  bool ci_find_substr(const std::string &hay, const std::string &needle);
  std::optional<std::string> findJESKeyForLabel(const std::string &label);
   
};

#endif
