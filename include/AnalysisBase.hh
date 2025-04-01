#ifndef AnalysisBase_h
#define AnalysisBase_h

#include <iostream>

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

#include "correction.h"

using namespace std;

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

  virtual bool PassEventFilter();
  virtual bool FastSimEventVeto(const ParticleList& GenJets);
  virtual double EGvalue(int jetIndex, int updown);
  virtual double GetPrefireWeight(int updown);
  
  // analysis functions
  virtual int GetNPV();
  virtual int GetNPUtrue();

  virtual bool GetMETtrigger();
  virtual bool GetMETORtrigger();

  virtual bool GetSingleElectrontrigger();
  virtual bool GetSingleMuontrigger();
  virtual bool GetDoubleElectrontrigger();
  virtual bool GetDoubleMuontrigger();
  virtual bool GetEMutrigger(); 
  
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

  bool IsSMS(){ return m_DoSMS; }
  bool IsData(){ return m_IsData; }
  bool IsFastSim(){ return m_IsFastSim; }
  
  string GetDataSet(){ return m_DataSet; }
  string GetFileTag(){ return m_FileTag; }
  int    GetYear(){ return m_year; }
  
protected:
  bool m_DoSMS;
  bool m_IsData;
  bool m_IsFastSim;
  
  virtual double GetEventWeight();
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
  std::map<long long,int>         m_HashToIndex;
  std::map<int,std::string> m_IndexToSample;
  std::map<int,double>      m_IndexToXsec;
  std::map<int,double>      m_IndexToNevent;
  std::map<int,double>      m_IndexToNweight;

  const Systematic* m_CurSys;
  const Systematic& CurrentSystematic() const;
  std::unique_ptr<correction::CorrectionSet> m_cset_Btag;
  virtual void clip_string(string& str, const string& clip){
    size_t pos = str.find(clip);
    if (pos != std::string::npos)
      str.erase(pos, clip.length());
  }
  bool minus_iso_hoe(int WPBitMap, int value, std::function<bool(int, int)> comp);
  
};


#endif
