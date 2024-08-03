#include <TH1D.h>
#include <iostream>

#include "AnalysisBase.hh"
#include "ParticleList.hh"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "NANORun3.hh"
#include "Leptonic.hh"

using namespace std;

template <class Base>
AnalysisBase<Base>::AnalysisBase(TTree* tree)
  : Base(tree), m_Systematics(true)
{
  m_Nsample = 0;
  m_SampleIndex = 0;
  m_DoSMS = false;
  m_IsData = false;
  m_IsFastSim = false;

  m_CurSys = &Systematic::Default();
}

template <class Base>
AnalysisBase<Base>::~AnalysisBase() {}

template <class Base>
void AnalysisBase<Base>::SetSystematic(const Systematic& sys){
  m_CurSys = &sys;
}

template <class Base>
const Systematic& AnalysisBase<Base>::CurrentSystematic() const {
  return *m_CurSys;
}

template <class Base>
void AnalysisBase<Base>::AddSystematics(){
  m_Systematics.Add(m_SysTool.GetTreeSystematics());
}

template <class Base>
void AnalysisBase<Base>::AddJESSystematics(){
  m_Systematics.Add(m_SysTool.JESSystematics());
}

template <class Base>
void AnalysisBase<Base>::AddJERSystematics(){
  m_Systematics.Add(m_SysTool.JERSystematics());
}

template <class Base>
void AnalysisBase<Base>::AddMETSystematics(){
  m_Systematics.Add(m_SysTool.METSystematics());
}

template <class Base>
void AnalysisBase<Base>::AddEESSystematics(){
  m_Systematics.Add(m_SysTool.EESSystematics());
}

template <class Base>
void AnalysisBase<Base>::AddMMSSystematics(){
  m_Systematics.Add(m_SysTool.MMSSystematics());
}

template <class Base>
void AnalysisBase<Base>::DoSMS(){
  m_DoSMS = true;
  m_IsData = false;
}

template <class Base>
void AnalysisBase<Base>::DoData(){
  m_IsData = true;
  m_IsFastSim = false;
  m_DoSMS = false;
}

template <class Base>
void AnalysisBase<Base>::DoFastSim(){
  m_IsFastSim = true;
  m_IsData = false;
}

template <class Base>
string AnalysisBase<Base>::GetEntry(int entry){
  if (!Base::fChain) return 0;

  Base::ClearEvent();
  Base::fChain->GetEntry(entry);
  m_SampleIndex = GetSampleIndex();
  
  return m_IndexToSample[m_SampleIndex];
}

template <class Base>
int AnalysisBase<Base>::GetSampleIndex(){
  if(m_Nsample == 0){
    m_IndexToSample[0] = "KUAnalysis";
    m_IndexToXsec[0] = 1.;
    m_IndexToNevent[0] = 1.;
    m_IndexToNweight[0] = 1.;
    m_Nsample++;
  }

  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetXsec(){
  if(m_Nsample)
    return m_IndexToXsec[m_SampleIndex];
  else
    return 0.;
}
  
template <class Base>
void AnalysisBase<Base>::AddLabels(const string& dataset, const string& filetag){
  m_DataSet = dataset;
  m_FileTag = filetag;
}

template <class Base>
void AnalysisBase<Base>::AddEventCountFile(const string& rootfile){
  m_NeventTool.BuildMap(rootfile);
}

template <class Base>
void AnalysisBase<Base>::AddFilterEffFile(const string& rootfile){
  m_NeventTool.BuildFilterEffMap(rootfile);
}

template <class Base>
void AnalysisBase<Base>::AddJSONFile(const string& jsonfile){
  m_JSONTool.BuildMap(jsonfile);
}

template <class Base>
void AnalysisBase<Base>::AddPUFolder(const string& pufold){
  m_PUTool.BuildMap(pufold);
}

template <class Base>
void AnalysisBase<Base>::AddBtagFolder(const string& btagfold){
  m_BtagSFTool.BuildMap(btagfold);
}

template <class Base>
void AnalysisBase<Base>::AddLepFolder(const string& lepfold){
  m_LepSFTool.BuildMap(lepfold);
}

template <class Base>
void AnalysisBase<Base>::AddJMEFolder(const string& jmefold){
  m_JMETool.BuildMap(jmefold);
  m_JMETool.BuildJERMap(jmefold);
}

template <class Base>
void AnalysisBase<Base>::AddMETTriggerFile(const string& csvfile){
  m_METTriggerTool.BuildMap(csvfile);
}

template <class Base>
void AnalysisBase<Base>::AddPrefireFile(const string& prefirefile){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  bool UseEMpT = false; // seems to always be false from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  m_PrefireTool = PrefireTool(year,UseEMpT,prefirefile);
}

template <class Base>
void AnalysisBase<Base>::InitializeHistograms(vector<TH1D*>& histos){}

template <class Base>
void AnalysisBase<Base>::BookHistograms(vector<TH1D*>& histos){}

template <class Base>
double AnalysisBase<Base>::DeltaPhiMin(const vector<TLorentzVector>& JETs, const TVector3& MET, int N){
  double dphimin = acos(-1);
  int Njet = JETs.size();
  for(int i = 0; i < Njet; i++){
    if(N > 0 && i >= N) break;
    if(fabs(JETs[i].Vect().DeltaPhi(MET)) < dphimin) dphimin = fabs(JETs[i].Vect().DeltaPhi(MET));
  }
  return dphimin;
}

template <class Base>
double AnalysisBase<Base>::DeltaPhiMin(const vector<pair<TLorentzVector, bool> >& JETs, const TVector3& MET, int N){
  double dphimin = acos(-1);
  int Njet = JETs.size();
  for(int i = 0; i < Njet; i++){
    if(N > 0 && i >= N) break;
    if(fabs(JETs[i].first.Vect().DeltaPhi(MET)) < dphimin) dphimin = fabs(JETs[i].first.Vect().DeltaPhi(MET));
  }
  return dphimin;
}

template <class Base>
int AnalysisBase<Base>::GetRunNum(){
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetLumiNum(){
  return 0;
}

template <class Base>
long AnalysisBase<Base>::GetEventNum(){
  return 0;
}

template <class Base>
bool AnalysisBase<Base>::PassEventFilter(){
  return true;
}

template <class Base>
bool AnalysisBase<Base>::FastSimEventVeto(const ParticleList& GenJets){
  return true;
}

template <class Base>
double AnalysisBase<Base>::GetPrefireWeight(int updown){
  return 1.;
}

template <class Base>
double AnalysisBase<Base>::EGvalue(int jetIndex, int updown){
  return 1.;
}

template <class Base>
bool AnalysisBase<Base>::IsHEM(Particle part){
  if(part.Eta() > -3.2 && part.Eta() < -1.2 && part.Phi() > -1.77 && part.Phi() < -0.67)
    return true;

  return false;
}

template <class Base>
double AnalysisBase<Base>::GetEventWeight(){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetPUWeight(int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuFWeight(int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuRWeight(int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetPDFWeight(int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown, ParticleIDType tag){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetElIDSFWeight(const ParticleList& els, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetElISOSFWeight(const ParticleList& els, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetElSIPSFWeight(const ParticleList& els, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetElVLIDSFWeight(const ParticleList& els, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuIDSFWeight(const ParticleList& mus, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuISOSFWeight(const ParticleList& mus, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuSIPSFWeight(const ParticleList& mus, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMuVLIDSFWeight(const ParticleList& mus, int updown){
  return 0;
}

template <class Base>
double AnalysisBase<Base>::GetMETTriggerSFWeight(double MET, double HT, int Nele, int Nmu, int updown){
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetMETTriggerSFCurve(double HT, int Nele, int Nmu){
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetNPV(){
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetNPUtrue(){
  return 0;
}

template <class Base>
bool AnalysisBase<Base>::IsGoodEvent(){
  return false;
}

template <class Base>
TVector3 AnalysisBase<Base>::GetMET(){
  return TVector3(0.,0.,0.);
}

template <class Base>
TVector3 AnalysisBase<Base>::GetAltMET(){
  return TVector3(0.,0.,0.);
}

template <class Base>
double AnalysisBase<Base>::Get_LHE_HT(){
  return 0.;
}

template <class Base>
double AnalysisBase<Base>::Get_LHE_HTIncoming(){
  return 0.;
}

template <class Base>
TVector3 AnalysisBase<Base>::GetPV(bool& good){
  good = false;
  return TVector3();
}

template <class Base>
bool AnalysisBase<Base>::GetMETtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetMETHTtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetMETORtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetSingleElectrontrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetSingleMuontrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetDoubleElectrontrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetDoubleMuontrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetEMutrigger(){
  return false;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetJetsMET(TVector3& MET, int id){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetJets(int id){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenJets(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetElectrons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetMuons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenElectrons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenMuons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenNeutrinos(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenBosons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenSparticles(){
  return ParticleList();
}

template <class Base>
std::pair<int,int> AnalysisBase<Base>::GetSUSYMasses(){

  return std::pair<int,int>(0,0);
}

template <class Base>
void AnalysisBase<Base>::MomTensorCalc(vector<TLorentzVector>& input, vector<double>& eigenvalues, double power, bool threeD){

  eigenvalues.clear();
  
  int N = input.size();

  if(threeD){
    if(N <= 0){
      for(int i = 0; i < 3; i++) eigenvalues.push_back(0.);
      return;
    }
    if(N == 1){
      eigenvalues.push_back(1.);
      for(int i = 0; i < 2; i++) eigenvalues.push_back(0.);
      return;
    }
    
    TMatrixDSym momTensor(3);
    momTensor.Zero();

    double norm = 0.;
    double P = 0.;
    double pnorm = 0.;
    for(int i = 0; i < N; i++){
      P = input[i].P();
      if( P > 0. ){
	norm += pow(P, power);
	pnorm = pow(P, power - 2.);
	momTensor(0,0) += pnorm*input[i].Px()*input[i].Px();
	momTensor(0,1) += pnorm*input[i].Px()*input[i].Py();
	momTensor(0,2) += pnorm*input[i].Px()*input[i].Pz();
	momTensor(1,0) += pnorm*input[i].Py()*input[i].Px();
	momTensor(1,1) += pnorm*input[i].Py()*input[i].Py();
	momTensor(1,2) += pnorm*input[i].Py()*input[i].Pz();
	momTensor(2,0) += pnorm*input[i].Pz()*input[i].Px();
	momTensor(2,1) += pnorm*input[i].Pz()*input[i].Py();
	momTensor(2,2) += pnorm*input[i].Pz()*input[i].Pz();
      }
    }
    if(norm > 0.){
      momTensor = (1./norm)*momTensor;
      TVectorD evalues(3);
      momTensor.EigenVectors(evalues);
      for(int i = 0; i < 3; i++) eigenvalues.push_back(evalues(i));
      return;
    } else {
      for(int i = 0; i < 3; i++) eigenvalues.push_back(0.);
      return;
    }

  } else { // transverse
    if(N <= 0){
      for(int i = 0; i < 2; i++) eigenvalues.push_back(0.);
      return;
    }
    if(N == 1){
      eigenvalues.push_back(1.);
      eigenvalues.push_back(0.);
      return;
    }

    TMatrixDSym momTensor(2);
    momTensor.Zero();

    double norm = 0.;
    double P = 0.;
    double pnorm = 0.;
    for(int i = 0; i < N; i++){
      P = input[i].Pt();
      if( P > 0. ){
	norm += pow(P, power);
	pnorm = pow(P, power - 2.);
	momTensor(0,0) += pnorm*input[i].Px()*input[i].Px();
	momTensor(0,1) += pnorm*input[i].Px()*input[i].Py();
	momTensor(1,0) += pnorm*input[i].Py()*input[i].Px();
	momTensor(1,1) += pnorm*input[i].Py()*input[i].Py();
      }
    }
    if(norm > 0.){
      momTensor = (1./norm)*momTensor;
      TVectorD evalues(2);
      momTensor.EigenVectors(evalues);
      for(int i = 0; i < 2; i++) eigenvalues.push_back(evalues(i));
      return;
    } else{
      for(int i = 0; i < 2; i++) eigenvalues.push_back(0.);
      return;
    }

  }
} 

/////////////////////////////////////////////////
// NANORun3 specific methods
/////////////////////////////////////////////////

template <>
int AnalysisBase<NANORun3>::GetRunNum(){
  return run;
}

template <>
int AnalysisBase<NANORun3>::GetLumiNum(){
  return luminosityBlock;
}

template <>
long AnalysisBase<NANORun3>::GetEventNum(){
  return event;
}

template <>
bool AnalysisBase<NANORun3>::PassEventFilter(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
//
//  if(year == 2016){
//    return Flag_goodVertices &&
//      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
//      Flag_HBHENoiseFilter &&
//      Flag_HBHENoiseIsoFilter &&
//      Flag_EcalDeadCellTriggerPrimitiveFilter &&
//      Flag_BadPFMuonFilter;
//  }
//  if(year == 2017){
//    return Flag_goodVertices &&
//      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
//      Flag_HBHENoiseFilter &&
//      Flag_HBHENoiseIsoFilter &&
//      Flag_EcalDeadCellTriggerPrimitiveFilter &&
//      Flag_BadPFMuonFilter;
//  }
//  if(year == 2018){
//    return Flag_goodVertices &&
//      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
//      Flag_HBHENoiseFilter &&
//      Flag_HBHENoiseIsoFilter &&
//      Flag_EcalDeadCellTriggerPrimitiveFilter &&
//      Flag_BadPFMuonFilter;
//  }
  
  return true;
}

template<>
bool AnalysisBase<NANORun3>::FastSimEventVeto(const ParticleList& GenJets){

  ParticleList jets;
  for(int i = 0; i < nJet; i++){
    if(Jet_pt[i] < 20. || fabs(Jet_eta[i]) > 2.5)
      continue;  
    if(Jet_chEmEF[i] > 0.1)
      continue;

    Particle jet;
    float mass = Jet_mass[i];
    if(std::isnan(mass))
      mass = 0;
    if(std::isinf(mass))
      mass = 0;
    if(mass < 0.)
      mass = 0.;
    jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i],
		     Jet_phi[i], mass);
    jets.push_back(jet);
  }
  jets.RemoveOverlap(GenJets, 0.3);
  if(jets.size() > 0)
    return false;

  return true;
}

template<>
double AnalysisBase<NANORun3>::EGvalue(int jetIndex, int updown){
  double PhotonMinPt = 20.;
  double PhotonMaxPt = 500.;
  double PhotonMinEta = 2.;
  double PhotonMaxEta = 3.;
  double phopf = 1.;

  vector<int> PhotonInJet;

  for(int p = 0; p < nPhoton; p++)
    {
      if(Photon_jetIdx[p] == jetIndex){
	if(Photon_pt[p] >= PhotonMinPt && fabs(Photon_eta[p]) <= PhotonMaxEta && fabs(Photon_eta[p]) >= PhotonMinEta){
	  double phopf_temp = 1. - m_PrefireTool.GetPrefireProbability(false, Photon_eta[p], Photon_pt[p], PhotonMaxPt, updown);
	  double elepf_temp = 1.;
	  if(Photon_electronIdx[p] > -1){
	    if(Electron_pt[Photon_electronIdx[p]] >= PhotonMinPt && fabs(Electron_eta[Photon_electronIdx[p]]) <= PhotonMaxEta && fabs(Electron_eta[Photon_electronIdx[p]]) >= PhotonMinEta){
	      elepf_temp = 1. - m_PrefireTool.GetPrefireProbability(false, Electron_eta[Photon_electronIdx[p]], Electron_pt[Photon_electronIdx[p]], PhotonMaxPt, updown);
	    }
	  }
	  phopf *= min(phopf_temp,elepf_temp);
	  PhotonInJet.push_back(p);
	}   
      }
    }
  for(int e = 0; e < nElectron; e++)
    {
      if(Electron_jetIdx[e] == jetIndex && std::count(PhotonInJet.begin(), PhotonInJet.end(), Electron_photonIdx[e]) == 0){
	if(Electron_pt[e] >= PhotonMinPt && fabs(Electron_eta[e]) <= PhotonMaxEta && fabs(Electron_eta[e]) >= PhotonMinEta){
	  phopf *= 1. - m_PrefireTool.GetPrefireProbability(false, Electron_eta[e], Electron_pt[e], PhotonMaxPt, updown);
	}
      }
    }
  return phopf;
}

template<>
double AnalysisBase<NANORun3>::GetPrefireWeight(int updown){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    return 1.; // no prefire weight for 2018

  double JetMinPt = 20.;
  double JetMaxPt = 500.;
  double JetMinEta = 2.;
  double JetMaxEta = 3.;


  double prefw = 1.;

  // loop over jets
  for(int j = 0; j < nJet; j++){
    double jetpf = 1.0;
    double jetpt = Jet_pt[j];
    if(m_PrefireTool.Get_UseEMpT())
      jetpt *= Jet_chEmEF[j] + Jet_neEmEF[j];
    if(jetpt >= JetMinPt && fabs(Jet_eta[j]) <= JetMaxEta && fabs(Jet_eta[j]) >= JetMinEta)
      jetpf *= 1. - m_PrefireTool.GetPrefireProbability(true, Jet_eta[j], jetpt, JetMaxPt, updown);
    double phopf = EGvalue(j,updown);
    prefw *= std::min(jetpf,phopf);
  }
  prefw *= EGvalue(-1,updown); //loop over all photons/electrons not associated to jets
 
  return prefw;
}

template <>
int AnalysisBase<NANORun3>::GetNPV(){
  return nOtherPV+1;
}

template <>
int AnalysisBase<NANORun3>::GetNPUtrue(){
  if(!IsData())
    return Pileup_nPU;
  
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetMETtrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  if(m_FileTag.find("22") != std::string::npos)
    year = 2022;
  if(m_FileTag.find("23") != std::string::npos)
    year = 2023;
  if(m_FileTag.find("24") != std::string::npos)
    year = 2024;

  if(year == 2016)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  if(year == 2017 || year == 2018 || year == 2022 || year == 2023 || year == 2024)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

  return 0;
}
  
template <>
bool AnalysisBase<NANORun3>::GetMETHTtrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  if(m_FileTag.find("22") != std::string::npos)
    year = 2022;
  if(m_FileTag.find("23") != std::string::npos)
    year = 2023;
  if(m_FileTag.find("24") != std::string::npos)
    year = 2024;

  if(year == 2016)
    return (HLT_PFHT300_PFMET100 ||
	    HLT_PFHT300_PFMET110);
  if(year == 2017 ||
     year == 2018)
    return (HLT_PFHT500_PFMET100_PFMHT100_IDTight ||
	    HLT_PFHT500_PFMET110_PFMHT110_IDTight ||
	    HLT_PFHT700_PFMET85_PFMHT85_IDTight ||
	    HLT_PFHT700_PFMET95_PFMHT95_IDTight ||
	    HLT_PFHT800_PFMET75_PFMHT75_IDTight ||
	    HLT_PFHT800_PFMET85_PFMHT85_IDTight);
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetMETORtrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_PFMETNoMu90_PFMHTNoMu90_IDTight ||
	    HLT_PFMETNoMu100_PFMHTNoMu100_IDTight ||
	    HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET90_PFMHT90_IDTight ||
	    HLT_PFMET100_PFMHT100_IDTight ||
	    HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned ||
	    HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight);
  if(year == 2017 ||
     year == 2018)
    return (HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMET130_PFMHT130_IDTight ||
	    HLT_PFMET140_PFMHT140_IDTight ||
	    HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMETNoMu130_PFMHTNoMu130_IDTight ||
	    HLT_PFMETNoMu140_PFMHTNoMu140_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

*/
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetSingleElectrontrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_Ele27_WPTight_Gsf);
  if(year == 2017 ||
     year == 2018)
    return (HLT_Ele32_WPTight_Gsf);
*/
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetSingleMuontrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_IsoMu24 ||
            HLT_IsoTkMu24);
  if(year == 2017 ||
     year == 2018)
    return (HLT_IsoMu24 ||
            HLT_IsoTkMu24);
*/
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetDoubleElectrontrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
  if(year == 2017 ||
     year == 2018)
    return (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
*/
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetDoubleMuontrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
  if(year == 2017 ||
     year == 2018)
    return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
*/
  return 0;
}

template <>
bool AnalysisBase<NANORun3>::GetEMutrigger(){
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

/*
  if(year == 2016)
    return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
  if(year == 2017 ||
     year == 2018)
    return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
*/
  return 0;
}

template <>
std::pair<int,int> AnalysisBase<NANORun3>::GetSUSYMasses(){
  if(!IsData()){
    int MP = 0;
    int MC = 0;
    int Ngen = nGenPart;
    for(int i = 0; i < Ngen; i++){
      int PDGID = abs(GenPart_pdgId[i]);
      if(PDGID > 1000000 && PDGID < 3000000){
	int mass = int(GenPart_mass[i]+0.5);
	if(PDGID == 1000022)
	  MC = mass;
	else
	  if(mass > MP)
	    MP = mass;
      }
    }
  
    return std::pair<int,int>(MP,MC);
  } else {
    return std::pair<int,int>(0,0);
  }
}

template <>
int AnalysisBase<NANORun3>::GetSampleIndex(){
  if(!m_DoSMS){
    if(m_Nsample == 0){
      m_IndexToSample[0]  = "KUAnalysis";
      m_IndexToXsec[0]    = m_XsecTool.GetXsec_BKG(m_DataSet);
      m_IndexToNevent[0]  = m_NeventTool.GetNevent_BKG(m_DataSet, m_FileTag);
      m_IndexToNweight[0] = m_NeventTool.GetNweight_BKG(m_DataSet, m_FileTag);
      m_Nsample++;
    }
    return 0;
  }
  
  int MP = 0;
  int MC = 0;
  int Ngen = nGenPart;
  int PDGID;
  for(int i = 0; i < nGenPart; i++){
    PDGID = fabs(GenPart_pdgId[i]);
    if(PDGID > 1000000 && PDGID < 3000000){
      int mass = int(GenPart_mass[i]+0.5);
      if(PDGID == 1000022)
	MC = mass;
      else
	if(mass > MP)
	  MP = mass;
    }
  }
  
  int hash = 100000*MP + MC;
  if(m_HashToIndex.count(hash) == 0){
    m_HashToIndex[hash] = m_Nsample;
    m_IndexToSample[m_Nsample]  = std::string(Form("SMS_%d_%d", MP, MC));
    m_IndexToXsec[m_Nsample]    = m_XsecTool.GetXsec_SMS(m_DataSet, MP);
    m_IndexToNevent[m_Nsample]  = m_NeventTool.GetNevent_SMS(m_DataSet, m_FileTag, MP, MC);
    m_IndexToNweight[m_Nsample] = m_NeventTool.GetNweight_SMS(m_DataSet, m_FileTag, MP, MC);
  
    m_Nsample++;
  }

  return m_HashToIndex[hash];
}


template <>
double AnalysisBase<NANORun3>::GetEventWeight(){
  if(IsData())
    return 1.;
  
  if(m_IndexToNweight[m_SampleIndex] > 0.){
    if(!m_DoSMS)
      return genWeight*m_IndexToXsec[m_SampleIndex]/m_IndexToNweight[m_SampleIndex];

    //cout << "genWeight " << genWeight << endl;
    //cout << "Xsec " << m_IndexToXsec[m_SampleIndex] << endl;
    //cout << "Nweight " << m_IndexToNweight[m_SampleIndex] << endl;
    //cout << "Filter eff " << m_NeventTool.GetFilterEff(m_DataSet,m_FileTag,luminosityBlock) << endl;
    return genWeight*m_IndexToXsec[m_SampleIndex]/m_IndexToNweight[m_SampleIndex]*m_NeventTool.GetFilterEff(m_DataSet,m_FileTag,luminosityBlock);
  } else
    return 0.;
}

template <>
double AnalysisBase<NANORun3>::GetPUWeight(int updown){
  if(IsData())
    return 1.;

  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  return m_PUTool.GetWeight(Pileup_nPU, year, updown);
}

// [0] is muR=0.5 muF=0.5 ; [1] is muR=0.5 muF=1.0 ; [2] is muR=0.5 muF=2.0 ;
// [3] is muR=0.1 muF=0.5 ; [4] is muR=1.0 muF=1.0 ; [5] is muR=1.0 muF=2.0 ;
// [6] is muR=2.0 muF=0.5 ; [7] is muR=2.0 muF=1.0 ; [8] is muR=2.0 muF=2.0 ;

template <>
double AnalysisBase<NANORun3>::GetMuFWeight(int updown){
  if(IsData())
    return 1.;
  if(nLHEScaleWeight == 0) return 1.;
  if(updown > 0)
    return LHEScaleWeight[5];
  else if(updown < 0) 
    return LHEScaleWeight[3];
  else
    return LHEScaleWeight[4]; //nominal
}

template <>
double AnalysisBase<NANORun3>::GetMuRWeight(int updown){
  if(IsData())
    return 1.;
  if(nLHEScaleWeight == 0) return 1.;
  if(updown > 0)
    return LHEScaleWeight[7];
  else if(updown < 0) 
    return LHEScaleWeight[1];
  else
    return LHEScaleWeight[4]; //nominal
  return 1.;
}

#ifdef _CMSSW_

template<>
void AnalysisBase<NANORun3>::AddLHAPDF(){
  if(IsData() || IsSMS())
    return;
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
 
  m_LHETool.AddLHAPDF(year);
}

#endif

template <>
double AnalysisBase<NANORun3>::GetPDFWeight(int updown){
  if(IsData() || IsSMS())
    return 1.;
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  return m_LHETool.GetWeight(nLHEPdfWeight,LHEPdfWeight,Generator_id1,Generator_id2,Generator_x1,Generator_x2,Generator_scalePDF,year,updown);
}

template <>
double AnalysisBase<NANORun3>::GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown, ParticleIDType tag){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Njet = jets.size();
  int iflavor;
  double EFF, SF;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Njet; i++){
    if(abs(jets[i].PDGID()) == 5)
      iflavor = 0;
    else if(abs(jets[i].PDGID()) == 4)
      iflavor = 1;
    else
      iflavor = 2;

    if(HForLF && iflavor == 2)
      continue;
    if(!HForLF && iflavor != 2)
      continue;
    
    EFF = m_BtagSFTool.EFF(jets[i].Pt(), year, iflavor, FastSim);
    SF  = m_BtagSFTool.SF(jets[i].Pt(), year, iflavor, updown);
    if(FastSim)
      SF *= m_BtagSFTool.SF(jets[i].Pt(), year, iflavor, updown, FastSim);

    if(jets[i].BtagID() >= tag){
      probMC   *= EFF;
      probDATA *= SF*EFF;
    } else {
      probMC   *= (1.-EFF);
      probDATA *= (1.-SF*EFF);
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetElIDSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(els[i].Pt(), els[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataIdError(els[i].Pt(), els[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(els[i].ParticleID() >= kMedium){
      probMC   *= EFFMC;
      probDATA *= SF*EFFMC;
    } else {
      probMC   *= EFFMC;
      probDATA *= EFFMC/SF;
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetElISOSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(els[i].Pt(), els[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataIsoError(els[i].Pt(), els[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(els[i].Pt(), els[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(els[i].RelIso() < 4. && els[i].MiniIso() < 4.){
      probMC   *= EFFMC;
      probDATA *= SF*EFFMC;
    } else {
      probMC   *= EFFMC;
      probDATA *= EFFMC/SF;
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetElSIPSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(els[i].Pt(), els[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataSipError(els[i].Pt(), els[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(els[i].Pt(), els[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(els[i].LepQual() != kBronze){
      if(els[i].SIP3D() < 2.){
	probMC   *= EFFMC;
	probDATA *= SF*EFFMC;
      } else {
	probMC   *= EFFMC;
	probDATA *= EFFMC/SF;
      }
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetElVLIDSFWeight(const ParticleList& els, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(els[i].Pt(), els[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(els[i].Pt(), els[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // apply VL SF to all leptons
    probDATA *= SF;
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetMuIDSFWeight(const ParticleList& mus, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(mus[i].Pt(), mus[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataIdError(mus[i].Pt(), mus[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(mus[i].ParticleID() >= kMedium){
      probMC   *= EFFMC;
      probDATA *= SF*EFFMC;
    } else {
      probMC   *= EFFMC;
      probDATA *= EFFMC/SF;
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetMuISOSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(mus[i].Pt(), mus[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataIsoError(mus[i].Pt(), mus[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(mus[i].RelIso() < 4. && mus[i].MiniIso() < 4.){
      probMC   *= EFFMC;
      probDATA *= SF*EFFMC;
    } else {
      probMC   *= EFFMC;
      probDATA *= EFFMC/SF;
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetMuSIPSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(mus[i].Pt(), mus[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataSipError(mus[i].Pt(), mus[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // Evaluate cut
    if(mus[i].LepQual() != kBronze){
      if(mus[i].SIP3D() < 2.){
	probMC   *= EFFMC;
	probDATA *= SF*EFFMC;
      } else {
	probMC   *= EFFMC;
	probDATA *= EFFMC/SF;
      }
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetMuVLIDSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, year);
      if(FSSF > 0){
	SF *= FSSF;
	EFFMC /= FSSF;
      }
    }

    // apply VL SF to all leptons
    probDATA *= SF;
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;

  return probDATA/probMC;
}

template <>
double AnalysisBase<NANORun3>::GetMETTriggerSFWeight(double MET, double HT, int Nele, int Nmu, int updown){
  if(IsData())
    return 1.;

  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  if(IsFastSim()){
    return m_METTriggerTool.Get_EFF(MET, HT, year,
				    (Nele > 0), (Nmu > 0),
				    false, updown)*
      m_METTriggerTool.Get_SF(MET, HT, year,
			      (Nele > 0), (Nmu > 0),
			      false, updown);
  } else {
    return m_METTriggerTool.Get_SF(MET, HT, year,
				   (Nele > 0), (Nmu > 0),
				   false, updown);
  }
}

template <>
int AnalysisBase<NANORun3>::GetMETTriggerSFCurve(double HT, int Nele, int Nmu){

  int year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  return m_METTriggerTool.Get_Curve_Index(HT, year, (Nele > 0), (Nmu > 0), IsData());
}

template <>
bool AnalysisBase<NANORun3>::IsGoodEvent(){
  return m_JSONTool.IsGood(run, luminosityBlock);
}

template <>
TVector3 AnalysisBase<NANORun3>::GetGenMET(){
  if(IsData())
    return TVector3();
  
  TVector3 vmet;
  vmet.SetPtEtaPhi(GenMET_pt,0.0,GenMET_phi);
  return vmet;
}

template <>
TVector3 AnalysisBase<NANORun3>::GetPV(bool& good){
  good = false;
  TVector3 PV(0.,0.,0.);
  
  if(PV_chi2 < 0.)
    return PV;
  
  if(PV_ndof < 5)
    return PV;
  
  if(fabs(PV_z) > 24.)
    return PV;
  
  if(PV_x*PV_x + PV_y*PV_y > 4.)
    return PV;
  
  good = true;
  PV.SetXYZ(PV_x,PV_y,PV_z);
  
  return PV;
}

template <>
void AnalysisBase<NANORun3>::InitializeHistograms(vector<TH1D*>& histos){
  // nPU
  TH1D* h_nPU = new TH1D("hist_NPU", "hist_NPU", 75, 0., 75.);
  histos.push_back(h_nPU);

  // Btag efficiencies
  vector<double> bin_edges;
  bin_edges.push_back(20.);
  bin_edges.push_back(30.);
  bin_edges.push_back(40.);
  bin_edges.push_back(50.);
  bin_edges.push_back(60.);
  bin_edges.push_back(70.);
  bin_edges.push_back(85.);
  bin_edges.push_back(100.);
  bin_edges.push_back(120.);
  bin_edges.push_back(140.);
  bin_edges.push_back(170.);
  bin_edges.push_back(200.);
  bin_edges.push_back(250.);
  bin_edges.push_back(300.);
  bin_edges.push_back(400.);
  bin_edges.push_back(600.);
  bin_edges.push_back(800.);
  bin_edges.push_back(1000.);

  TH1D* h_btag[3][2]; // [flavor][den/num]
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      h_btag[i][j] = (TH1D*) new TH1D(Form("hist_btag_flavor%d_%s", i, (j == 0 ? "den" : "num")),
				      Form("hist_btag_flavor%d_%s", i, (j == 0 ? "den" : "num")),
				      17, &bin_edges[0]);
      histos.push_back(h_btag[i][j]);
    }
  }
}

template <>
void AnalysisBase<NANORun3>::BookHistograms(vector<TH1D*>& histos){
  int year = 2017;
  if(m_FileTag.find("16") != std::string::npos)
    year = 2016;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int ihist = 0;

  // nPU
  histos[ihist]->Fill(GetNPUtrue());

  ihist++;

  // Btag efficiencies
  int Njet = nJet;
  for(int i = 0; i < Njet; i++){
    if(Jet_pt[i] < 20. || fabs(Jet_eta[i] > 2.4))
      continue;

    bool btag = false;
    if(year == 2016)
      if(Jet_btagDeepFlavB[i] > 0.3093)
	btag = true;
    if(year == 2017)
      if(Jet_btagDeepFlavB[i] > 0.3033)
	btag = true;
    if(year == 2018)
      if(Jet_btagDeepFlavB[i] > 0.2770)
	btag = true;

    int flavor;
    if(abs(Jet_partonFlavour[i]) == 5)
      flavor = 0;
    else if(abs(Jet_partonFlavour[i]) == 4)
      flavor = 1;
    else
      flavor = 2;

    histos[ihist+2*flavor]->Fill(Jet_pt[i]);
    if(btag) histos[ihist+2*flavor+1]->Fill(Jet_pt[i]);
  }
}

template <>
ParticleList AnalysisBase<NANORun3>::GetGenJets(){
  ParticleList list;
  if(IsData())
    return list;
  int year = 2017;
  if(m_FileTag.find("16") != std::string::npos)
    year = 2016;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  int NGenjet = nGenJet;

  for(int i = 0; i < NGenjet; i++){
    if(GenJet_pt[i] < 15. || fabs(GenJet_eta[i]) > 5.)
      continue;
    Particle jet;
    float mass = GenJet_mass[i];
    if(std::isnan(mass))
      mass = 0;
    if(std::isinf(mass))
      mass = 0;
    if(mass < 0.)
      mass = 0.;
    jet.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i],
		     GenJet_phi[i], mass);
    //jet.SetPDGID( GenPart_pdgId[i] );

    list.push_back(jet);
  }
  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetJetsMET(TVector3& MET, int id){
  int year = 2017;
  if(m_FileTag.find("16") != std::string::npos)
    year = 2016;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  
  ParticleList list;
  bool passID = true;
  int Njet = nJet;

  double delta  = (CurrentSystematic().IsUp() ? 1. : -1.);
  TVector3 deltaMET(0.,0.,0.);
  bool DO_JES = false;
  if(m_SysTool.JESSystematics() == CurrentSystematic())
    DO_JES = true;

  bool DO_JER = false;
  if(m_SysTool.JERSystematics() == CurrentSystematic())
    DO_JER = true;
  
  for(int i = 0; i < Njet; i++){
    bool failID = false;
    //if(Jet_pt[i] < 15. || fabs(Jet_eta[i]) > 5.)
    //  continue;
    //if(Jet_jetId[i] < id)
    //  continue;
    
    Particle jet;
    float mass = Jet_mass[i];
    if(std::isnan(mass))
      mass = 0;
    if(std::isinf(mass))
      mass = 0;
    if(mass < 0.)
      mass = 0.;
    jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i],
		     Jet_phi[i], mass);
    
    if(DO_JES){
      double uncer = m_JMETool.GetJESFactor(year, CurrentSystematic().Label(),
					    Jet_pt[i], Jet_eta[i]);
      
      deltaMET -= delta*uncer*jet.Vect();
      jet.SetPtEtaPhiM((1.+delta*uncer)*jet.Pt(),
		       jet.Eta(), jet.Phi(),
		       (1.+delta*uncer)*jet.M());
    }
    
    // recalibrate jets
    // double raw = 1. - Jet_rawFactor[i];
    // double L1 = m_JMETool.GetJESFactor(year, "L1FastJet",
    // 				     raw*Jet_pt[i], Jet_eta[i],
    // 				     Jet_area[i], fixedGridRhoFastjetAll);
    // double L2 = m_JMETool.GetJESFactor(year, "L2Relative",
    // 				     raw*L1*Jet_pt[i], Jet_eta[i],
    // 				     Jet_area[i], fixedGridRhoFastjetAll);

    // cout << raw << " " << L1 << " " << L2 << " " << raw*L1*L2 << endl;

    if(!IsData()){

      // JER recipe based on https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
     
      double smearFactor = 1.;
      //double JER = m_JMETool.GetJERFactor(year, Jet_pt[i], Jet_eta[i], fixedGridRhoFastjetAll); // using this for rho based on: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/0127d46a973e894d97e9a16bd3939f421b2b689e/python/postprocessing/modules/jme/jetmetUncertainties.py#L49
      double JER = m_JMETool.GetJERFactor(year, Jet_pt[i], Jet_eta[i], Rho_fixedGridRhoFastjetAll); // need to check this rho for UL and Run3
      double SF = m_JMETool.GetJERSFFactor(year,Jet_eta[i],0);

      if(DO_JER)
        SF = m_JMETool.GetJERSFFactor(year,Jet_eta[i],delta);
       
      //cout << SF << " " << JER << " " << Jet_pt[i] << " " << Jet_eta[i] << " " << Njet << " " << nGenJet <<  endl;
      
      // check for gen jet matching:
      bool gen_match = false;
      Particle genJet;
      genJet.SetPtEtaPhiM(0.,0.,0.,0.);
      
      for(int g = 0; g < nGenJet; g++){
	genJet.SetPtEtaPhiM(GenJet_pt[g],GenJet_eta[g],GenJet_phi[g],GenJet_mass[g]);
        if(fabs(Jet_pt[i] - GenJet_pt[g]) < 3.*JER*Jet_pt[i] && jet.DeltaR(genJet) < 0.2){
          gen_match = true;
          break;
        }
      }

      // 3 different cases to consider
      // Case 1: we have a "good" gen level jet matched to reco jet
      if(gen_match){
        double dPt = jet.Pt() - genJet.Pt();
        smearFactor = 1. + (SF - 1.)*dPt/jet.Pt();
      }

      // Case 2: Smear jet pT using a random Gaussian variation
      else if(!gen_match && SF > 1.){
        TRandom3 rand3;
        rand3.SetSeed(event);
        double rand_val = rand3.Gaus(0.,JER);
        smearFactor = 1.+rand_val*sqrt(SF*SF-1.);
      }

      // Case 3: Resolution in data is better than res in sim so do nothing
      else
        smearFactor = 1.;
      
      if(smearFactor*jet.E() < 1.e-2)
        smearFactor = 1.e-2/jet.E();
      
      Particle oldJet = jet;
      jet.SetPtEtaPhiM(jet.Pt()*smearFactor,jet.Eta(),jet.Phi(),jet.M()*smearFactor);
      //deltaMET += (oldJet-jet).Vect(); //should be default?
      deltaMET -= (oldJet-jet).Vect();

    } //end JER

    if(Jet_pt[i] < 15. || fabs(Jet_eta[i]) > 5.)
      continue;
    if(Jet_jetId[i] < id)
      continue;
    
    if(Jet_jetId[i] >= 3)
      jet.SetParticleID(kTight);
    else if(Jet_jetId[i] >= 2) 
      jet.SetParticleID(kMedium);
    else if(Jet_jetId[i] >= 1)
      jet.SetParticleID(kLoose);
    
    // DeepCSV tagger
    // jet.SetBtag(Jet_btagDeepB[i]);

    // DeepFlavour tagger
    jet.SetBtag(Jet_btagDeepFlavB[i]);

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
    if(year == 2016){
      // Deep CSV
      // if(jet.Btag() > 0.8953)
      // 	jet.SetBtagID(kTight);
      // else if(jet.Btag() > 0.6321) 
      // 	jet.SetBtagID(kMedium);
      // else if(jet.Btag() > 0.2217)
      // 	jet.SetBtagID(kLoose);

      // Deep Flavor
      if(jet.Btag() > 0.7221)
	jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.3093) 
	jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0614)
	jet.SetBtagID(kLoose);
    }

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    if(year == 2017){
      // DeepCSV
      // if(jet.Btag() > 0.8001)
      // 	jet.SetBtagID(kTight);
      // else if(jet.Btag() > 0.4941) 
      // 	jet.SetBtagID(kMedium);
      // else if(jet.Btag() > 0.1522)
      // 	jet.SetBtagID(kLoose);

      // Deep Flavor
      if(jet.Btag() > 0.7489)
	jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.3033) 
	jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0521)
	jet.SetBtagID(kLoose);
    }

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    if(year == 2018){
      // DeepCSV
      // if(jet.Btag() > 0.7527)
      // 	jet.SetBtagID(kTight);
      // else if(jet.Btag() > 0.4184) 
      // 	jet.SetBtagID(kMedium);
      // else if(jet.Btag() > 0.1241)
      // 	jet.SetBtagID(kLoose);

      // DeepFlavor
      if(jet.Btag() > 0.7264)
	jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.2770) 
	jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0494)
	jet.SetBtagID(kLoose);
    }

    jet.SetPDGID( Jet_partonFlavour[i] );
      
    list.push_back(jet);
  }

  // If one jet fails jet ID, 
  if(!passID)
    return ParticleList();
  
  //if(year == 2017)
  //  MET.SetPtEtaPhi(METFixEE2017_pt,0.0,METFixEE2017_phi);
  //else
    MET.SetPtEtaPhi(MET_pt,0.0,MET_phi);
  
  deltaMET.SetZ(0.);
  MET += deltaMET;
  
  if(CurrentSystematic() == Systematic("METUncer_UnClust")){
    //if(year == 2017)
    //  deltaMET.SetXYZ(delta*METFixEE2017_MetUnclustEnUpDeltaX,
    //    	      delta*METFixEE2017_MetUnclustEnUpDeltaY, 0.);
    //else
      deltaMET.SetXYZ(delta*MET_MetUnclustEnUpDeltaX,
		      delta*MET_MetUnclustEnUpDeltaY, 0.);
    MET += deltaMET;
  }

  if(CurrentSystematic() == Systematic("METUncer_GenMET"))
    MET.SetPtEtaPhi(GenMET_pt,0.,GenMET_phi);
  
  return list;
}

template <>
TVector3 AnalysisBase<NANORun3>::GetMET(){
  TVector3 MET;
  GetJetsMET(MET);

  return MET;
}

template <>
TVector3 AnalysisBase<NANORun3>::GetAltMET(){
  TVector3 MET;
  MET.SetPtEtaPhi(MET_pt,0.0,MET_phi);

  return MET;
}

template <>
double AnalysisBase<NANORun3>::Get_LHE_HT(){
  return LHE_HT;
}

template <>
double AnalysisBase<NANORun3>::Get_LHE_HTIncoming(){
  return LHE_HTIncoming;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetJets(int id){
  TVector3 dum;
  return GetJetsMET(dum, id);
}

template <>
ParticleList AnalysisBase<NANORun3>::GetElectrons(){
  int year = 2017;
  if(m_FileTag.find("16") != std::string::npos)
    year = 2016;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;

  ParticleList list;

  int N = nElectron;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Electron_pt[i] < 5. || fabs(Electron_eta[i]) > 2.5)
      continue;
    // if(fabs(Electron_dxy[i]) >= 0.05 || fabs(Electron_dz[i]) >= 0.1 ||
    //    Electron_ip3d[i] >= 0.0175 || Electron_sip3d[i] >= 2.5)
    if(fabs(Electron_dxy[i]) >= 0.05 || fabs(Electron_dz[i]) >= 0.1 ||
       Electron_sip3d[i] >= 8)
      continue;
    if(Electron_pfRelIso03_all[i]*Electron_pt[i] >= 20. + 300./Electron_pt[i])
      continue;

    Particle lep;
    lep.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i],
		     Electron_phi[i], std::max(Electron_mass[i],float(1.e-6)));
    lep.SetPDGID( (Electron_charge[i] < 0. ? 11 : -11) );
    lep.SetCharge( (Electron_charge[i] < 0. ? -1 : 1) );

    lep.SetDxy(Electron_dxy[i]);
    lep.SetDxyErr(Electron_dxyErr[i]);
    lep.SetDz(Electron_dz[i]);
    lep.SetDzErr(Electron_dzErr[i]);
    lep.SetIP3D(Electron_ip3d[i]);
    lep.SetSIP3D(Electron_sip3d[i]);
    lep.SetIsLowPt(false);

    lep.SetRelIso(Electron_pfRelIso03_all[i]);
    lep.SetMiniIso(Electron_miniPFRelIso_all[i]);

    // https://twiki.cern.ch/twiki/pub/CMS/SUSLeptonSF/Run2_SUSYwp_EleCB_MVA_8Jan19.pdf
    
    // FO baseline criteria
    if(Electron_lostHits[i] == 0 && Electron_convVeto[i]){

//      double mva = Electron_mvaFall17V1noIso[i];
//      if(year == 2016 || year == 2018)
//	mva = Electron_mvaFall17V2noIso[i];
	double mva = Electron_mvaNoIso[i];


      // convert to raw MVA output
      if(mva == -1.)
	mva = -999.;
      else if(mva == 1.)
	mva = 999.;
      else
	mva = -0.5*log((1.-mva)/(1.+mva));
      
      // FO VLoose
      if(year == 2016){ // Summer16_94X legacy
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > -0.259)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.388 + 0.109*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.388)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > -0.256)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.696 + 0.106*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.696)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > -1.630)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -1.219 + 0.148*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -1.219)
	      lep.SetParticleID(kVeryLoose);
	  }
	}
      }

      if(year == 2017){ // Fall17_94X
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > -0.135)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.93 + (0.043/15.)*(lep.Pt()-10.)))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.887)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > -0.417)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.93 + (0.04/15.)*(lep.Pt()-10.)))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.89)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > -0.470)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.942 + (0.032/15.)*(lep.Pt()-10.)))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.91)
	      lep.SetParticleID(kVeryLoose);
	  }
	}
      }

      if(year == 2018){ // Autumn18_102X
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > 0.053)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.106 + 0.062*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.106)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > -0.434)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.769 + 0.038*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -0.769)
	      lep.SetParticleID(kVeryLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > -0.956)
	      lep.SetParticleID(kVeryLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -1.461 + 0.042*(lep.Pt() - 25.))
	      lep.SetParticleID(kVeryLoose);
	  } else {
	    if(mva > -1.461)
	      lep.SetParticleID(kVeryLoose);
	  }
	}
      }

      // VLoose electron
      if(year == 2016){ // Summer16_94X legacy
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > 1.309)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > 0.887 + 0.088*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > 0.887)
	      lep.SetParticleID(kLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > 0.373)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > 0.112 + 0.099*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > 0.112)
	      lep.SetParticleID(kLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > 0.071)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.017 + 0.137*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > -0.017)
	      lep.SetParticleID(kLoose);
	  }
	}
      }

      if(year == 2017){ // Fall17_94X
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > 0.488)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.788 + (0.148/15.)*(lep.Pt()-10.)) )
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > -0.64)
	      lep.SetParticleID(kLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > -0.045)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.85 + (0.075/15.)*(lep.Pt()-10.)))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > -0.775)
	      lep.SetParticleID(kLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > 0.176)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > (-0.81 + (0.077/15.)*(lep.Pt()-10.)))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > -0.733)
	      lep.SetParticleID(kLoose);
	  }
	}
      }

      if(year == 2018){ // Autumn18_102X
	if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	  if(lep.Pt() < 10.){
	    if(mva > 1.320)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > 1.204 + 0.066*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > 1.204)
	      lep.SetParticleID(kLoose);
	  }
	} else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	  if(lep.Pt() < 10.){
	    if(mva > 0.192)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > 0.084 + 0.033*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > 0.084)
	      lep.SetParticleID(kLoose);
	  }
	} else { // eta < 2.5
	  if(lep.Pt() < 10.){
	    if(mva > 0.362)
	      lep.SetParticleID(kLoose);
	  } else if(lep.Pt() < 25.) {
	    if(mva > -0.123 + 0.053*(lep.Pt() - 25.))
	      lep.SetParticleID(kLoose);
	  } else {
	    if(mva > -0.123)
	      lep.SetParticleID(kLoose);
	  }
	}
      }
	    
      // signal lepton IDs (only Tight for now) baseline criteria
      // if(lep.IP3D() < 0.01 && lep.SIP3D() < 2.){
      if(true){
	// Tight electron
	if(year == 2016){ // Summer16_94X legacy
	  if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 1.309)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 40.) {
	      if(mva > 3.447 + 0.063*(lep.Pt()- 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 4.392)
		lep.SetParticleID(kTight);
	    }
	  } else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.373)
		lep.SetParticleID(kTight); // just changed me
	    } else if(lep.Pt() < 40.) {
	      if(mva > 2.522 + 0.058*(lep.Pt() - 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 3.392)
		lep.SetParticleID(kTight);
	    }
	  } else { // eta < 2.5
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.071)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 40.) {
	      if(mva > 1.555 + 0.075*(lep.Pt() - 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 2.680)
		lep.SetParticleID(kTight);
	    }
	  }
	}

	if(year == 2017){ // Fall17_94X
	  if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.488)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > 0.2+0.032*(lep.Pt() - 10.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 0.68)
		lep.SetParticleID(kTight);
	    }
	  } else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > -0.045)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > 0.1+0.025*(lep.Pt() - 10.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 0.475)
		lep.SetParticleID(kTight);
	    }
	  } else { // eta < 2.5
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.176)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > -0.1+0.028*(lep.Pt() - 10.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 0.32)
		lep.SetParticleID(kTight);
	    }
	  }
	}

	if(year == 2018){ // Autumn18_102X
	  if(fabs(lep.Eta()) < 0.8){ // eta < 0.8
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 1.320)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > 4.277 + 0.112*(lep.Pt() - 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 4.277)
		lep.SetParticleID(kTight);
	    }
	  } else if(fabs(lep.Eta()) < 1.479){ // eta < 1.479
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.192)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > 3.152 + 0.060*(lep.Pt() - 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 3.152)
		lep.SetParticleID(kTight);
	    }
	  } else { // eta < 2.5
	    if(lep.Pt() < 10.){ // using VLoose ID for low pT
	      if(mva > 0.362)
		lep.SetParticleID(kTight);
	    } else if(lep.Pt() < 25.) {
	      if(mva > 2.359 + 0.087*(lep.Pt() - 25.))
		lep.SetParticleID(kTight);
	    } else {
	      if(mva > 2.359)
		lep.SetParticleID(kTight);
	    }
	  }
	}
	
      }
    } // end lepton id

    // Lepton quality
    // matching old BFI qual logic exactly
    /*
    if(lep.ParticleID() >= kMedium &&
       (lep.RelIso()*lep.Pt() < 4. && lep.MiniIso()*lep.Pt() < 4.)){
      if(lep.SIP3D() < 2.)
	lep.SetLepQual(kGold);
      else
	lep.SetLepQual(kSilver);
    } else {
      lep.SetLepQual(kBronze);
    }
    */
     //its okay to cut on kMedium for elecs, they are either tight(4) or loose(2&1) - this makes the logic identical to old BFI
     if(lep.ParticleID() < kMedium || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
	    lep.SetLepQual(kBronze);
	  else if(lep.SIP3D() > 2.)
	    lep.SetLepQual(kSilver);
	  else
	    lep.SetLepQual(kGold);

    list.push_back(lep);
  }

  // Adding the lowpt electrons here

  int N1 = nLowPtElectron;
  for(int i = 0; i < N1; i++){
    // baseline lepton definition
    if(LowPtElectron_pt[i] < 1. || LowPtElectron_pt[i] >= 5. || fabs(LowPtElectron_eta[i]) > 2.5)
      continue;
    if(LowPtElectron_convVeto[i] == 0)
      continue;
    if(fabs(LowPtElectron_dxy[i]) >= 0.05 || fabs(LowPtElectron_dz[i]) >= 0.1)
      continue;
    if(LowPtElectron_ID[i] < 1.4) // need to tune after feedback from Brady
      continue;
    
    //if(LowPtElectron_dxyErr[i] < 1.e-8 || LowPtElectron_dzErr[i] < 1.e-8)
    //  continue;
    
    // Calculate IP_3D and SIP_3D = IP_3D / IP_3D_err for the LowPtElectron collection.
    // - IP_3D:     3D impact parameter wrt first PV, in cm
    // - SIP_3D:    3D impact parameter significance wrt first PV, in cm
    float dxy       = LowPtElectron_dxy[i];
    float dz        = LowPtElectron_dz[i];
    float dxy_err   = LowPtElectron_dxyErr[i];
    float dz_err    = LowPtElectron_dzErr[i];
    float IP_3D     = sqrt(dxy*dxy + dz*dz);
    float SIP_3D    = IP_3D*IP_3D / sqrt((dxy*dxy)*(dxy_err*dxy_err) + (dz*dz)*(dz_err*dz_err));

    // TODO: Review and tune SIP_3D cut.
    //if (SIP_3D >= 8)
    //  continue;

    // FIXME: fix PFIso requirement for the LowPtElectron collection.
    //if(Electron_pfRelIso03_all[i]*Electron_pt[i] >= 20. + 300./Electron_pt[i])
    //  continue;

    Particle lep;
    lep.SetPtEtaPhiM(LowPtElectron_pt[i], LowPtElectron_eta[i],
		     LowPtElectron_phi[i], std::max(LowPtElectron_mass[i],float(1.e-6)));
    lep.SetPDGID( (LowPtElectron_charge[i] < 0. ? 11 : -11) );
    lep.SetCharge( (LowPtElectron_charge[i] < 0. ? -1 : 1) );

    lep.SetDxy(LowPtElectron_dxy[i]);
    lep.SetDxyErr(LowPtElectron_dxyErr[i]);
    lep.SetDz(LowPtElectron_dz[i]);
    lep.SetDzErr(LowPtElectron_dzErr[i]);
    lep.SetIP3D(IP_3D);
    lep.SetSIP3D(SIP_3D);
    lep.SetIsLowPt(true);

    lep.SetRelIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetMiniIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetParticleID(kMedium);
    
    //if (LowPtElectron_pt[i] < 3.0 || (LowPtElectron_pt[i] >= 3.0 && LowPtElectron_embeddedID[i] >=6))
    //  lep.SetParticleID(kTight);
    if (LowPtElectron_pt[i] < 3.0 || (LowPtElectron_pt[i] >= 3.0 && LowPtElectron_ID[i] >=1.4)) // need to tune after feedback from Brady
      lep.SetParticleID(kTight);
 
    list.push_back(lep);
  }

  return list;

}

template <>
ParticleList AnalysisBase<NANORun3>::GetMuons(){
  int year = 2017;
  if(m_FileTag.find("16") != std::string::npos)
    year = 2016;
  if(m_FileTag.find("18") != std::string::npos)
    year = 2018;
  
  ParticleList list;

  int N = nMuon;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Muon_pt[i] < 3. || fabs(Muon_eta[i]) > 2.4)
      continue;
    // if(fabs(Muon_dxy[i]) >= 0.05 || fabs(Muon_dz[i]) >= 0.1 ||
    //    Muon_ip3d[i] >= 0.0175 || Muon_sip3d[i] >= 2.5)
    if(fabs(Muon_dxy[i]) >= 0.05 || fabs(Muon_dz[i]) >= 0.1 || Muon_sip3d[i] >= 8.)
      continue;
    if(Muon_pfRelIso03_all[i]*Muon_pt[i] >= 20. + 300./Muon_pt[i])
      continue;
    
    Particle lep;
    lep.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i],
		     Muon_phi[i], std::max(float(0.),Muon_mass[i]));
    lep.SetPDGID( (Muon_charge[i] < 0. ? 13 : -13) );
    lep.SetCharge( (Muon_charge[i] < 0. ? -1 : 1) );	
    lep.SetDxy(Muon_dxy[i]);
    lep.SetDxyErr(Muon_dxyErr[i]);
    lep.SetDz(Muon_dz[i]);
    lep.SetDzErr(Muon_dzErr[i]);
    lep.SetIP3D(Muon_ip3d[i]);
    lep.SetSIP3D(Muon_sip3d[i]);
    lep.SetIsLowPt(false);

    lep.SetRelIso(Muon_pfRelIso03_all[i]);
    lep.SetMiniIso(Muon_miniPFRelIso_all[i]);

    // FO baseline criteria
    if(true){
      lep.SetParticleID(kLoose);

      // signal lep criteria
      //if(lep.IP3D() < 0.01 && lep.SIP3D() < 2.){
      if(true){
	if(Muon_tightId[i])
	  lep.SetParticleID(kTight);
	else if(lep.Pt() < 0.){
	  if(Muon_softId[i])
	    lep.SetParticleID(kMedium);
	} else {
	  if(Muon_mediumId[i])
	    lep.SetParticleID(kMedium);
	}
      }
    }

    // Lepton quality
    /*
    if(lep.ParticleID() >= kMedium &&
       (lep.RelIso()*lep.Pt() < 4. && lep.MiniIso()*lep.Pt() < 4.)){
      if(lep.SIP3D() < 2.)
	lep.SetLepQual(kGold);
      else
	lep.SetLepQual(kSilver);
    } else {
      lep.SetLepQual(kBronze);
    }
    */
    if(lep.ParticleID() < kMedium || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
            lep.SetLepQual(kBronze);
          else if(lep.SIP3D() > 2.)
            lep.SetLepQual(kSilver);
          else
            lep.SetLepQual(kGold);

    
    list.push_back(lep);
  }

  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetGenElectrons(){
  ParticleList list;

  if(IsData())
    return list;
  
  int N = nGenPart;
  int PDGID;
  for(int i = 0; i < N; i++){
    PDGID = GenPart_pdgId[i];
    if(abs(PDGID) == 11 && GenPart_pt[i] > 2. && GenPart_status[i] == 1){
      Particle lep;
      
      lep.SetPDGID(PDGID);
      int mom = GenPart_genPartIdxMother[i];
      if(mom >= 0 && mom < N){

	int momID = GenPart_pdgId[mom];
	int momStatus = GenPart_status[mom];

	while(abs(momID) == 11){

	  if(momStatus == 23){
	    lep.SetMomPDGID(PDGID);
	    lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
	    break;
	  }

	  mom = GenPart_genPartIdxMother[mom];
	  if(mom < 0 || mom >= N)
	    continue;
	  momID = GenPart_pdgId[mom];
	  momStatus = GenPart_status[mom];

	}
	
	lep.SetMomPDGID(momID);
	lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
      }
     
      lep.SetCharge( (PDGID > 0 ? -1 : 1) );
      lep.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i],
		       GenPart_phi[i], max(float(0.),GenPart_mass[i]));
      
      list.push_back(lep);
    }
  }

  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetGenMuons(){
  ParticleList list;

  if(IsData())
    return list;
  
  int N = nGenPart;
  int PDGID;
  for(int i = 0; i < N; i++){
    PDGID = GenPart_pdgId[i];
    if(abs(PDGID) == 13 && GenPart_pt[i] > 2. && GenPart_status[i] == 1){
      Particle lep;
      
      lep.SetPDGID(PDGID);
      int mom = GenPart_genPartIdxMother[i];
      if(mom >= 0 && mom < N){

	int momID = GenPart_pdgId[mom];
	int momStatus = GenPart_status[mom];

	while(abs(momID) == 13){

	  if(momStatus == 23){
	    lep.SetMomPDGID(PDGID);
	    lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
	    break;
	  }

	  mom = GenPart_genPartIdxMother[mom];
	  if(mom < 0 || mom >= N)
	    continue;
	  momID = GenPart_pdgId[mom];
	  momStatus = GenPart_status[mom];

	}
	
	lep.SetMomPDGID(momID);
	lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
      }
     
      lep.SetCharge( (PDGID > 0 ? -1 : 1) );
      lep.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i],
		       GenPart_phi[i], max(float(0.),GenPart_mass[i]));
      
      list.push_back(lep);
    }
  }

  return list;
}



template <>
ParticleList AnalysisBase<NANORun3>::GetGenNeutrinos(){
  ParticleList list;

  if(IsData())
    return list;
  
  int N = nGenPart;
  int PDGID;
  for(int i = 0; i < N; i++){
    PDGID = GenPart_pdgId[i];
    if(abs(PDGID) == 12 || abs(PDGID) == 14 || abs(PDGID) == 16){
      Particle lep;
      
      lep.SetPDGID(PDGID);
      int mom = GenPart_genPartIdxMother[i];
      if(mom >= 0 && mom < N)
	lep.SetMomPDGID(GenPart_pdgId[mom]);
      lep.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i],
		       GenPart_phi[i], max(float(0.),GenPart_mass[i]));

      list.push_back(lep);
    }
  }

  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetGenBosons(){
  ParticleList list;

  if(IsData())
    return list;
  
  int N = nGenPart;
  int PDGID;
  for(int i = 0; i < N; i++){
    PDGID = GenPart_pdgId[i];
    if(abs(PDGID) == 23 || abs(PDGID) == 24 || abs(PDGID) == 25){
      Particle p;
      
      p.SetPDGID(PDGID);
      int mom = GenPart_genPartIdxMother[i];
      if(mom >= 0 && mom < N)
	p.SetMomPDGID(GenPart_pdgId[mom]);
      p.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i],
		     GenPart_phi[i], max(float(0.),GenPart_mass[i]));

      list.push_back(p);
    }
  }

  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetGenSparticles(){
  ParticleList list;

  if(IsData())
    return list;
  
  int N = nGenPart;
  int PDGID;
  for(int i = 0; i < N; i++){
    PDGID = GenPart_pdgId[i];
    if(abs(PDGID) >= 1000000 && abs(PDGID) < 3000000){
      Particle p;
      
      p.SetPDGID(PDGID);
      int mom = GenPart_genPartIdxMother[i];
      if(mom >= 0 && mom < N)
	p.SetMomPDGID(GenPart_pdgId[mom]);
      p.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i],
		     GenPart_phi[i], max(float(0.),GenPart_mass[i]));

      list.push_back(p);
    }
  }

  return list;
}

template class AnalysisBase<NANORun3>;

