#include <TH1D.h>
#include <iostream>

#include "AnalysisBase.hh"
#include "ParticleList.hh"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "SUSYNANOBase.hh"
#include "NANOULBase.hh"
#include "NANORun3.hh"
#include "Leptonic.hh"

using namespace std;

template <class Base>
AnalysisBase<Base>::AnalysisBase(TTree* tree)
  : Base(tree), m_Systematics(true)
{
  m_Nsample = 0;
  m_SampleIndex = 0;
  m_IsSMS = false;
  m_IsData = false;
  m_IsFastSim = false;

  m_CurSys = &Systematic::Default();
}

template <class Base>
AnalysisBase<Base>::~AnalysisBase() {}

/////////////////////////////////////////////////
// Start Base generic methods
/////////////////////////////////////////////////

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
  m_IsSMS = true;
  m_IsData = false;
}

template <class Base>
void AnalysisBase<Base>::DoData(){
  m_IsData = true;
  m_IsFastSim = false;
  m_IsSMS = false;
}

template <class Base>
void AnalysisBase<Base>::DoFastSim(){
  m_IsFastSim = true;
  m_IsData = false;
}

template <class Base>
void AnalysisBase<Base>::DoPrivateMC(){
  m_IsPrivateMC = true;
  m_IsData = false;
}

template <class Base>
void AnalysisBase<Base>::DoCascades(){
  m_IsCascades = true;
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
double AnalysisBase<Base>::GetXsec(){
  if(m_Nsample)
    return m_IndexToXsec[m_SampleIndex];
  else
    return 0.;
}

template <class Base>
double AnalysisBase<Base>::GetNevent(){
  if(m_Nsample)
    return m_IndexToNevent[m_SampleIndex];
  else
    return 0.;
}

template <class Base>
double AnalysisBase<Base>::GetNweight(){
  if(m_Nsample)
    return m_IndexToNweight[m_SampleIndex];
  else
    return 0.;
}

template <class Base>
double AnalysisBase<Base>::GetFilterEff(){
  if(m_IsSMS && std::is_member_object_pointer<decltype(&Base::luminosityBlock)>::value)
    return m_NeventTool.GetFilterEff(m_DataSet,m_FileTag,this->luminosityBlock);
  else
    return 1.;
}

template <class Base>
void AnalysisBase<Base>::AddLabels(const string& dataset, const string& filetag){
  m_DataSet = dataset;
  m_FileTag = filetag;
  m_year = 2016;
       if(m_FileTag.find("17") != std::string::npos) m_year = 2017;
  else if(m_FileTag.find("18") != std::string::npos) m_year = 2018;
  else if(m_FileTag.find("22") != std::string::npos) m_year = 2022;
  else if(m_FileTag.find("23") != std::string::npos) m_year = 2023;
  else if(m_FileTag.find("24") != std::string::npos) m_year = 2024;
  else if(m_FileTag.find("25") != std::string::npos) m_year = 2025;
  else if(m_FileTag.find("26") != std::string::npos) m_year = 2026;
  if(m_FileTag.find("APV") != std::string::npos) m_IsAPV = true;
  if(m_FileTag.find("UL") != std::string::npos) m_IsUL = true;
  if(m_FileTag.find("EE") != std::string::npos) m_IsEE = true;
  if(m_FileTag.find("BPix") != std::string::npos) m_IsBPix = true;
  if(m_FileTag.find("130X") != std::string::npos) m_IsRun3 = true;
  if(m_FileTag.find("Cascades") != std::string::npos) DoCascades();
  if(m_FileTag.find("SMS") != std::string::npos) DoSMS();
  if(m_FileTag.find("Data") != std::string::npos) DoData();
  m_XsecTool.SetFileTag(filetag);
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
void AnalysisBase<Base>::AddBtagFolder(const string& btagfold, const string& proc_rootfile, int year){
  if(!m_IsUL and m_year < 2019)
    m_BtagSFTool.BuildMap(btagfold, proc_rootfile, year);
  else{
    std::string Btag_file = find_clib_file(btagfold, "btagging.json.gz");
    m_cset_Btag = correction::CorrectionSet::from_file(Btag_file);
    m_BtagLooseWP = m_cset_Btag->at("deepJet_wp_values")->evaluate({"L"});
    m_BtagMediumWP = m_cset_Btag->at("deepJet_wp_values")->evaluate({"M"});
    m_BtagTightWP = m_cset_Btag->at("deepJet_wp_values")->evaluate({"T"});
    if(m_year > 2018){
      m_BtagVeryTightWP = m_cset_Btag->at("deepJet_wp_values")->evaluate({"XT"});
      m_BtagVeryVeryTightWP = m_cset_Btag->at("deepJet_wp_values")->evaluate({"XXT"});
    }
  }
}

template <class Base>
void AnalysisBase<Base>::AddLepFolder(const string& lepfold){
  m_LepSFTool.BuildMap(lepfold);
}

template <class Base>
bool AnalysisBase<Base>::ci_find_substr(const std::string &hay, const std::string &needle){
  auto it = std::search(
    hay.begin(), hay.end(),
    needle.begin(), needle.end(),
    [](char a, char b){ return std::tolower(a) == std::tolower(b); }
  );
  return it != hay.end();
}

template <class Base>
std::optional<std::string> AnalysisBase<Base>::findJESKeyForLabel(const std::string &label){
  using nlohmann::json;
  const auto &cfg = JecConfigReader::getCfgAK4();
  if (!cfg.contains(m_jmeYearKey)) return std::nullopt;
  const auto &yearObj = cfg.at(m_jmeYearKey);
  if (!yearObj.contains("ForUncertaintyJES")) return std::nullopt;
  const auto &fJES = yearObj.at("ForUncertaintyJES");

  const std::vector<std::string> sections = {
    "ForUncertaintyJESReduced",
    "ForUncertaintyJESFull",
    "ForUncertaintyJESTotal"
  };

  // normalise label: if legacy prefix exists, strip it for matching convenience
  std::string core = label;
  const std::string legacyPrefix = "JESUncer_";
  if (core.rfind(legacyPrefix, 0) == 0) core = core.substr(legacyPrefix.size());

  for (const auto &sec : sections) {
    if (!fJES.contains(sec)) continue;
    const auto &arr = fJES.at(sec);
    for (const auto &entry : arr) {
      if (!entry.is_array() || entry.size() < 2) continue;
      std::string corrKey = entry.at(0).get<std::string>();   // CLIB key
      std::string humanLabel = entry.at(1).get<std::string>(); // descriptive label

      // Try strict-ish matches first: contain core as token or exact human label match
      if (ci_find_substr(humanLabel, core) || ci_find_substr(corrKey, core)) {
        return corrKey;
      }
      // fallback: match full label (in case label is something like "CMS_scale_j_Total")
      if (ci_find_substr(humanLabel, label) || ci_find_substr(corrKey, label)) {
        return corrKey;
      }
    }
  }

  return std::nullopt;
}

template <class Base>
bool AnalysisBase<Base>::PassesJVM(const ParticleList& jets) const {
  for (const auto& jet : jets) {
    if (m_JMETool.JetInJVM(jet.Pt(), jet.Eta(), jet.Phi(), jet.jetID(), jet.ChEmEF(), jet.NeEmEF()))
      return false; // veto event if any jet inside veto region
  }
  return true;
}

template <class Base>
std::string AnalysisBase<Base>::getJMEYearKey() const {
  // Respect explicit override first
  if (!m_JMEYearOverride.empty()) return m_JMEYearOverride;

  // Map integer year + flags -> helper year-key
  if (m_year == 2016) {
      // APV == "preVFP"
      return m_IsAPV ? "2016Pre" : "2016Post";
  }
  if (m_year == 2017) return "2017";
  if (m_year == 2018) return "2018";

  // For 2022: try to guess Pre vs Post from file tag
  if (m_year == 2022) {
    // Heuristic: if filename/tag mentions "RunC" or "RunD" prefer Pre, otherwise Post.
    std::string tag = m_FileTag;
    if (tag.find("RunC") != std::string::npos || tag.find("RunD") != std::string::npos
     || tag.find("C") != std::string::npos || tag.find("D") != std::string::npos) {
        return "2022Pre";
    }
    return "2022Post";
  }

  // 2023 has run-dependent payloads (Pre/Post). Heuristic below:
  if (m_year == 2023) {
    // If the input tag explicitly references "BPix" the original code used a BPix-specific tag;
    // use that as a weak indicator for the "Post" payload in some workflows.
    if (m_IsBPix) return "2023Post";
    // Otherwise default to 2023Pre (safe fallback; override if needed)
    return "2023Pre";
  }

  if (m_year == 2024) return "2024";
  if (m_year == 2025) return "2025";
  if (m_year == 2026) return "2026";

  // Unknown -> empty string
  return "";
}

template <class Base>
void AnalysisBase<Base>::AddJMEFolder(const std::string& jmefold) {
  if (!m_IsUL && m_year < 2019) {
    m_JMETool.BuildMap(jmefold);
    m_JMETool.BuildJERMap(jmefold);
  } else {
    std::string METPhi_file;
    if(m_year < 2019) METPhi_file = find_clib_file(jmefold, "met.json.gz");
    else if(m_year == 2022 && !m_IsEE) METPhi_file = find_clib_file(jmefold, "met_xyCorrections_2022_2022.json.gz");
    else if(m_year == 2022 && m_IsEE) METPhi_file = find_clib_file(jmefold, "met_xyCorrections_2022_2022EE.json.gz");
    else if(m_year == 2023 && !m_IsBPix) METPhi_file = find_clib_file(jmefold, "met_xyCorrections_2023_2023.json.gz");
    else if(m_year == 2023 && m_IsBPix) METPhi_file = find_clib_file(jmefold, "met_xyCorrections_2023_2023BPix.json.gz");
    m_cset_METPhi = correction::CorrectionSet::from_file(METPhi_file);
    m_jmeYearKey = getJMEYearKey();
    // --- Set up Jet Veto Map ---
    try {
      auto jvm = JvmConfigReader::getJvmForYear(m_jmeYearKey);
      if (jvm.use) {
        m_JMETool.SetupJVM(jvm);
        std::cout << "[JVM] Loaded veto map for year " << m_year
                  << " (key=" << jvm.key << ")" << std::endl;
      } else {
        std::cout << "[JVM] No veto map configured for year " << m_year << std::endl;
      }
    } catch (const std::exception& e) {
      std::cerr << "[JVM] Error loading veto map for year " << m_year << ": " << e.what() << std::endl;
    }
  }
}

template <class Base>
void AnalysisBase<Base>::AddMETTriggerFile(const string& csvfile){
  m_METTriggerTool.BuildMap(csvfile);
}

template <class Base>
void AnalysisBase<Base>::AddXSecJSON(const string& XSjsonfile){
  m_XsecTool.UpdateXsecFromJSON(XSjsonfile);
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
bool AnalysisBase<Base>::PassEventFilter(){
  return true;
}

template <class Base>
void AnalysisBase<Base>::AddPrefireFile(const string& prefirefile){
  if(m_year > 2018) return; // prefire only need for Run2
  bool UseEMpT = false; // seems to always be false from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  m_PrefireTool = PrefireTool(m_year,UseEMpT,prefirefile);
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
TVector3 AnalysisBase<Base>::GetMET(){
  return TVector3(0.,0.,0.);
}

template <class Base>
TVector3 AnalysisBase<Base>::GetAltMET(){
  return TVector3(0.,0.,0.);
}

template <class Base>
bool AnalysisBase<Base>::GetMETtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetMETORtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetMETDoubleMutrigger(){
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
bool AnalysisBase<Base>::GetTripleElectrontrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetTripleMuonLowPTtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetTripleMuonHighPTtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetDiMuEleLowPTtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetDiMuEleHighPTtrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetDiEleMutrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetEMutrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetEMuMutrigger(){
  return false;
}

template <class Base>
bool AnalysisBase<Base>::GetEMuEtrigger(){
  return false;
}

template <class Base>
void AnalysisBase<Base>::ApplyMETPhiCorrections(TLorentzVector &metOut, int runNumber, int npvGood, const std::string &variation) {}

template <class Base>
ParticleList AnalysisBase<Base>::GetJetsMET(TVector3& MET, int id){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetJets(int id){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetElectrons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetLowPtElectrons(){
  return ParticleList();
}

template <class Base>
ParticleList AnalysisBase<Base>::GetMuons(){
  return ParticleList();
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

template <class Base>
bool AnalysisBase<Base>::minus_iso_hoe(int WPBitMap, int threshold, std::function<bool(int, int)> comp){
  if(!m_IsRun3) return true; // not needed for run2
  // Define the bit shifts corresponding to each cut
  std::pair<std::string, int> cuts[] = {
    {"MinPtCut", 0},
    {"GsfEleSCEtaMultiRangeCut", 3},
    {"GsfEleDEtaInSeedCut", 6},
    {"GsfEleDPhiInCut", 9},
    {"GsfEleFull5x5SigmaIEtaIEtaCut", 12},
    {"GsfEleEInverseMinusPInverseCut", 18},
    {"GsfEleConversionVetoCut", 24},
    {"GsfEleMissingHitsCut", 27}
  };

  // Check if all required bits pass threshold check
  for (const auto& cut : cuts) {
    int value = (WPBitMap >> cut.second) & 0b111;  // Extract 3-bit value
    if (!comp(value, threshold))
      return false;
  }
  return true; 
}

// Map UL campaign tokens to the Run2-style directory substrings
template <class Base>
std::string AnalysisBase<Base>::normalize_tag(const std::string& filetag) {
    if (filetag.find("Summer20UL16APV") != std::string::npos)
        return "Run2-2016preVFP-UL";
    if (filetag.find("Summer20UL16") != std::string::npos)
        return "Run2-2016postVFP-UL";
    if (filetag.find("Summer20UL17") != std::string::npos)
        return "Run2-2017-UL";
    if (filetag.find("Summer20UL18") != std::string::npos)
        return "Run2-2018-UL";

    // Keep as-is for Run 3 (Summer22, Summer23, etc.)
    return filetag;
}

// Strip known build suffixes (_106X, _130X, _102X) then apply UL->Run2 mapping
template <class Base>
std::string AnalysisBase<Base>::normalize_filetag(const std::string& filetag_in) {
    std::string tag = filetag_in;
    static const std::vector<std::string> suffixes = {"_106X", "_130X", "_102X"};
    for (const auto& suf : suffixes) {
        if (tag.size() > suf.size() &&
            tag.compare(tag.size() - suf.size(), suf.size(), suf) == 0) {
            tag = tag.substr(0, tag.size() - suf.size());
            break;
        }
    }
    return normalize_tag(tag);
}

// member function to extract NanoAOD version from a directory name
template <class Base>
int AnalysisBase<Base>::extract_nano_version(const std::string& dirname) {
    static const std::regex re("NanoAODv([0-9]+)");
    std::smatch m;
    if (std::regex_search(dirname, m, re)) {
        return std::stoi(m[1]);
    }
    return -1; // if not found
}

template <class Base>
std::string AnalysisBase<Base>::find_clib_file(const std::string& fold, const std::string& filename) {
    // normalize_filetag makes a modifiable copy and applies UL->Run2 mapping
    string filetag = m_FileTag;
    clip_string(filetag, "_Data");
    clip_string(filetag, "_SMS");
    clip_string(filetag, "_Cascades");
    filetag = normalize_filetag(filetag);
    std::vector<std::pair<int, fs::path>> matches;

    // Decide how to search to avoid accidental matches (e.g. Summer22 vs Summer22EE).
    // For Run2-mapped tags (start with "Run2-"), search for the tag as-is.
    // For campaign tokens like "Summer22", prefer "-Summer22-" to avoid "Summer22EE".
    std::string needle;
    if (filetag.rfind("Run2-", 0) == 0 || filetag.rfind("Run3-", 0) == 0) {
        needle = filetag; // match substring anywhere in dirname
    } else {
        needle = "-" + filetag + "-";
    }

    for (const auto& entry : fs::directory_iterator(fold)) {
        if (!entry.is_directory()) continue;
        const std::string dirname = entry.path().filename().string();
        if (dirname.find(needle) != std::string::npos) {
            int version = extract_nano_version(dirname);
            if (version > 0) {
                fs::path candidate = entry.path() / "latest" / filename;
                if (fs::exists(candidate)) {
                    matches.emplace_back(version, candidate);
                }
            }
        }
    }

    if (matches.empty()) {
        throw std::runtime_error("No matching clib file found for tag=" + filetag + " " + filename);
    }

    // Pick entry with highest NanoAOD version
    auto best = std::max_element(matches.begin(), matches.end(),
                                 [](const auto& a, const auto& b){ return a.first < b.first; });

    // Check if multiple entries have the same highest version (avoid silent ambiguity)
    int max_version = best->first;
    int count = std::count_if(matches.begin(), matches.end(),
                              [&](const auto& p){ return p.first == max_version; });

    if (count > 1) {
        std::string err = "Multiple matches for tag=" + filetag +
                          " with same NanoAODv" + std::to_string(max_version) + ":\n";
        for (const auto& m : matches) {
            if (m.first == max_version) err += "  " + m.second.string() + "\n";
        }
        throw std::runtime_error(err);
    }

    return best->second.string();
}

/////////////////////////////////////////////////
// End Base generic methods
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// Start NANOAOD generic methods
/////////////////////////////////////////////////

template <class Base>
int AnalysisBase<Base>::GetRunNum(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::run)>::value)
    return this->run;
  else
    return 0;
}

template <class Base>
int AnalysisBase<Base>::GetLumiNum(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::luminosityBlock)>::value)
    return this->luminosityBlock;
  else
    return 0;
}

template <class Base>
long AnalysisBase<Base>::GetEventNum(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::event)>::value)
    return this->event;
  else
    return 0;
}

template <class Base>
double AnalysisBase<Base>::Get_LHE_HT(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::LHE_HT)>::value) return this->LHE_HT;
  else return 0.;
}

template <class Base>
double AnalysisBase<Base>::Get_LHE_HTIncoming(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::LHE_HTIncoming)>::value) return this->LHE_HTIncoming;
  else return 0.;
}

template <class Base>
int AnalysisBase<Base>::GetNPV(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::nOtherPV)>::value) return this->nOtherPV+1;
  else return 0;
}

template <class Base>
int AnalysisBase<Base>::GetNPUtrue(){
  if(!IsData())
    if constexpr (std::is_member_object_pointer<decltype(&Base::Pileup_nPU)>::value) return this->Pileup_nPU;
  return 0;
}

template <class Base>
bool AnalysisBase<Base>::IsGoodEvent(){
  if constexpr (std::is_member_object_pointer<decltype(&Base::run)>::value &&
                std::is_member_object_pointer<decltype(&Base::luminosityBlock)>::value)
    return m_JSONTool.IsGood(this->run, this->luminosityBlock);
  else return false;
}

template <class Base>
int AnalysisBase<Base>::GetNpartons(){
  if(IsData()) return 0;
  if constexpr (std::is_member_object_pointer<decltype(&Base::LHE_Njets)>::value)
    return this->LHE_Njets;
  return 0;
}

template <class Base>
bool AnalysisBase<Base>::FastSimEventVeto(const ParticleList& GenJets){
  if constexpr (std::is_member_object_pointer<decltype(&Base::nJet)>::value &&
                std::is_member_object_pointer<decltype(&Base::Jet_pt)>::value &&
                std::is_member_object_pointer<decltype(&Base::Jet_eta)>::value &&
                std::is_member_object_pointer<decltype(&Base::Jet_phi)>::value &&
                std::is_member_object_pointer<decltype(&Base::Jet_mass)>::value &&
                std::is_member_object_pointer<decltype(&Base::Jet_chEmEF)>::value){
    ParticleList jets;
    for(int i = 0; i < this->nJet; i++){
      if(this->Jet_pt[i] < 20. || fabs(this->Jet_eta[i]) > 2.5)
        continue;  
      if(this->Jet_chEmEF[i] > 0.1)
        continue;
  
      Particle jet;
      float mass = this->Jet_mass[i];
      if(std::isnan(mass))
        mass = 0;
      if(std::isinf(mass))
        mass = 0;
      if(mass < 0.)
        mass = 0.;
      jet.SetPtEtaPhiM(this->Jet_pt[i], this->Jet_eta[i],
  		     this->Jet_phi[i], mass);
      jets.push_back(jet);
    }
    jets.RemoveOverlap(GenJets, 0.3);
    if(jets.size() > 0)
      return false;
  }
  return true;
}

template <class Base>
TVector3 AnalysisBase<Base>::GetPV(bool& good){
  good = false;
  TVector3 PV(0.,0.,0.);
  if constexpr (
    std::is_member_object_pointer<decltype(&Base::PV_chi2)>::value &&
    std::is_member_object_pointer<decltype(&Base::PV_ndof)>::value && 
    std::is_member_object_pointer<decltype(&Base::PV_x)>::value && 
    std::is_member_object_pointer<decltype(&Base::PV_y)>::value && 
    std::is_member_object_pointer<decltype(&Base::PV_z)>::value){ 

    int PV_chi2 = this->PV_chi2; int PV_ndof = this->PV_ndof;
    int PV_x = this->PV_x; int PV_y = this->PV_y; int PV_z = this->PV_z;
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
  } 
  return PV;
}

template <class Base>
TVector3 AnalysisBase<Base>::GetGenMET(){
  if(IsData())
    return TVector3();
  TVector3 vmet;
  if constexpr (
    std::is_member_object_pointer<decltype(&Base::GenMET_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenMET_phi)>::value)
    vmet.SetPtEtaPhi(this->GenMET_pt,0.0,this->GenMET_phi);
  return vmet;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenJets(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenJet)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenJet_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenJet_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenJet_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenJet_mass)>::value ){
    int NGenjet = this->nGenJet;

    for(int i = 0; i < NGenjet; i++){
      if(this->GenJet_pt[i] < 15. || fabs(this->GenJet_eta[i]) > 5.)
        continue;
      Particle jet;
      float mass = this->GenJet_mass[i];
      if(std::isnan(mass))
        mass = 0;
      if(std::isinf(mass))
        mass = 0;
      if(mass < 0.)
        mass = 0.;
      jet.SetPtEtaPhiM(this->GenJet_pt[i], this->GenJet_eta[i],
          	       this->GenJet_phi[i], mass);

      list.push_back(jet);
    }
  }
  return list;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenElectrons(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_status)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) == 11 && this->GenPart_pt[i] > 0.5 && this->GenPart_status[i] == 1){
        Particle lep;
        lep.SetPDGID(PDGID);
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          int momID = this->GenPart_pdgId[mom];
          int momStatus = this->GenPart_status[mom];
          while(abs(momID) == 11 && (mom >= 0 && mom < N)){
            if(momStatus == 23){
              lep.SetMomPDGID(PDGID);
              lep.SetGenMomIndex(mom);
              lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
              break;
            }
            mom = this->GenPart_genPartIdxMother[mom];
            if(mom < 0 || mom >= N)
              continue;
            momID = this->GenPart_pdgId[mom];
            momStatus = this->GenPart_status[mom];
          }
          lep.SetMomPDGID(momID);
          lep.SetGenMomIndex(mom);
          lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
        }
        lep.SetCharge( (PDGID > 0 ? -1 : 1) );
        lep.SetPtEtaPhiM(this->GenPart_pt[i], this->GenPart_eta[i],
          	         this->GenPart_phi[i], max(float(0.),this->GenPart_mass[i]));
        list.push_back(lep);
      }
    }
  }
  return list;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenMuons(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_status)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) == 13 && this->GenPart_pt[i] > 2. && this->GenPart_status[i] == 1){
        Particle lep;
        lep.SetPDGID(PDGID);
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          int momID = this->GenPart_pdgId[mom];
          int momStatus = this->GenPart_status[mom];
          while(abs(momID) == 13 && (mom >= 0 && mom < N)){
            if(momStatus == 23){
              lep.SetMomPDGID(PDGID);
              lep.SetGenMomIndex(mom);
              lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
              break;
            }
            mom = this->GenPart_genPartIdxMother[mom];
            if(mom < 0 || mom >= N)
              continue;
            momID = this->GenPart_pdgId[mom];
            momStatus = this->GenPart_status[mom];
          }
          lep.SetMomPDGID(momID);
          lep.SetGenMomIndex(mom);
          lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
        }
        lep.SetCharge( (PDGID > 0 ? -1 : 1) );
        lep.SetPtEtaPhiM(this->GenPart_pt[i], this->GenPart_eta[i],
          	       this->GenPart_phi[i], max(float(0.),this->GenPart_mass[i]));
        list.push_back(lep);
      }
    }
  }
  return list;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenNeutrinos(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) == 12 || abs(PDGID) == 14 || abs(PDGID) == 16){
        Particle lep;
        lep.SetPDGID(PDGID);
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          lep.SetMomPDGID(this->GenPart_pdgId[mom]);
          lep.SetGenMomIndex(mom);
        }
        lep.SetPtEtaPhiM(this->GenPart_pt[i], this->GenPart_eta[i],
          	       this->GenPart_phi[i], max(float(0.),this->GenPart_mass[i]));
        list.push_back(lep);
      }
    }
  }
  return list;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenBosons(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) == 23 || abs(PDGID) == 24 || abs(PDGID) == 25){
        Particle p;
        p.SetPDGID(PDGID);
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          p.SetMomPDGID(this->GenPart_pdgId[mom]);
          p.SetGenMomIndex(mom);
        }
        p.SetPtEtaPhiM(this->GenPart_pt[i], this->GenPart_eta[i],
          	     this->GenPart_phi[i], max(float(0.),this->GenPart_mass[i]));
        list.push_back(p);
      }
    }
  }
  return list;
}

template <class Base>
ParticleList AnalysisBase<Base>::GetGenSparticles(){
  ParticleList list;
  if(IsData())
    return list;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pt)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_eta)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_phi)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) >= 1000000 && abs(PDGID) < 3000000){
        Particle p;
        p.SetPDGID(PDGID);
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          p.SetMomPDGID(this->GenPart_pdgId[mom]);
          p.SetGenMomIndex(mom);
        }
        p.SetPtEtaPhiM(this->GenPart_pt[i], this->GenPart_eta[i],
          	     this->GenPart_phi[i], max(float(0.),this->GenPart_mass[i]));
        list.push_back(p);
      }
    }
  }
  return list;
}

template <class Base>
std::pair<int,int> AnalysisBase<Base>::GetSUSYMasses(){
  if(!IsData()){
      if constexpr (
        std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
        std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
        std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value ){
      int MP = 0;
      int MC = 0;
      int Ngen = this->nGenPart;
      for(int i = 0; i < Ngen; i++){
        int PDGID = abs(this->GenPart_pdgId[i]);
        if(PDGID > 1000000 && PDGID < 3000000){
          int mass = int(this->GenPart_mass[i]+0.5);
          if(PDGID == 1000022)
            MC = mass;
          else
            if(mass > MP)
              MP = mass;
        }
      }
    return std::pair<int,int>(MP,MC);
    } // if !constexpr
  } // if(!IsData)
  return std::pair<int,int>(0,0);
}

template <class Base>
int AnalysisBase<Base>::GetGenMass(const int& u_PDGID){
  if(!IsData())
      if constexpr (
        std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
        std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
        std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value )
      for(int i = 0; i < int(this->nGenPart); i++)
        if(abs(this->GenPart_pdgId[i]) == u_PDGID)
          return int(this->GenPart_mass[i]+0.5);
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetGenSUSYNBosons(const int& u_PDGID){
  int NBosons = 0;
  if(IsData())
    return 0;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID;
    for(int i = 0; i < N; i++){
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) == u_PDGID){
        int mom = this->GenPart_genPartIdxMother[i];
        if(mom >= 0 && mom < N){
          int momPDGID = this->GenPart_pdgId[mom];
          if(momPDGID > 1000000 && momPDGID < 3000000)
            NBosons++;
        }
      }
    }
  }
  return NBosons;
}

template <class Base>
int AnalysisBase<Base>::GetGenCascadesProduction(int& firstSpart, int& secondSpart){
  if(!IsCascades())
    return 0;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_status)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value){

    int N = this->nGenPart;
    int PDGID;
    bool SlepEP = false; // Positive Selectron
    bool SlepEM = false; // Negative Selectron
    bool SlepMP = false; // Positive Smuon
    bool SlepMM = false; // Negative Smuon
    bool SneuE  = false; // Selectron Sneutrino
    bool SneuM  = false; // Smuon Sneutrino
    bool found1st = false;
    for(int i = 0; i < N; i++){
      bool found = false;
      PDGID = this->GenPart_pdgId[i];
      if(abs(PDGID) > 1000010 && abs(PDGID) < 1000015 && this->GenPart_status[i] == 22){ // is Slep or Sneu
        if(PDGID == -1000011)          { SlepEP = true; found = true; }
        else if(PDGID == 1000011)      { SlepEM = true; found = true; }
        else if(PDGID == -1000013)     { SlepMP = true; found = true; }
        else if(PDGID == 1000013)      { SlepMM = true; found = true; }
        else if(abs(PDGID) == 1000012) { SneuE = true;  found = true; } 
        else if(abs(PDGID) == 1000014) { SneuM = true;  found = true; }
        if(found && !found1st) { firstSpart = i; found1st = true; }
        else if(found &&  found1st) { secondSpart = i; break; }
      }
    }
    if(SlepEP && SlepEM) return 1; // SlepSlep Electron
    if(SlepEP && SneuE)  return 2; // SlepSneu + Electron
    if(SlepEM && SneuE)  return 3; // SlepSneu - Electron
    if(SneuE && !(SlepEP || SlepEP)) return 4; // SneuSneu Electron
    if(SlepMP && SlepMM) return 5; // SlepSlep Muon
    if(SlepMP && SneuM)  return 6; // SlepSneu + Muon
    if(SlepMM && SneuM)  return 7; // SlepSneu - Muon
    if(SneuM && !(SlepMP || SlepMP)) return 8; // SneuSneu Muon
  }
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetGenCascadesProduction(){
  int dum1 = 0;
  int dum2 = 0;
  return GetGenCascadesProduction(dum1, dum2);
}

template <class Base>
std::pair<int,int> AnalysisBase<Base>::GetGenCascadesDecayMode(const int& GenIndex){ // get decay of given particle by index
  if(!IsCascades())
    return std::pair<int,int>(0,0);

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int mom, momPDGID;
    int D1PDGID = 0, D2PDGID = 0; // Daughter PDG IDs
    for(int i = 0; i < N; i++){
      mom = this->GenPart_genPartIdxMother[i];
      if(this->GenPart_pdgId[i] == this->GenPart_pdgId[GenIndex]) continue;
      if(mom >= 0 && mom < N){
        momPDGID = this->GenPart_pdgId[mom];
        while(momPDGID == this->GenPart_pdgId[GenIndex] && (mom >= 0 && mom < N)){
          if(mom == GenIndex){
            if(D1PDGID == 0) D1PDGID = this->GenPart_pdgId[i];
            else             D2PDGID = this->GenPart_pdgId[i];
          }
          mom = this->GenPart_genPartIdxMother[mom];
        }
      }
    }
    if(abs(D1PDGID) > abs(D2PDGID)) // convention that the 'D1' is always the lower abs(PDGID)
       std::swap(D1PDGID,D2PDGID);  // swap value of vars to fit convention
    return std::pair<int,int>(D1PDGID, D2PDGID);
  }
  return std::pair<int,int>(0,0);
}

template <class Base>
// get index of particle with PDGID that comes from given mother
// Example: Have a particle with PDGID of LSP and the index of the mom. Want to know the index of the LSP
int AnalysisBase<Base>::GetGenCascadesIndex(const int& u_momIndex, const int& u_PDGID){
  if(!IsCascades()) return -1;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){
    
    int N = this->nGenPart;
    int mom;
    for(int i = 0; i < N; i++){
      if(this->GenPart_pdgId[i] != u_PDGID) continue;
      mom = this->GenPart_genPartIdxMother[i];
      while(mom > 0 && mom < N){
        if(mom == u_momIndex) return i;
        mom = this->GenPart_genPartIdxMother[mom];
      }
    }

  }
  return -1;
}

template <class Base>
uint16_t AnalysisBase<Base>::GetGenCascadesTree(){
  if(!IsCascades()) return 0;

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value){

    int SlepSneu1st = -1;
    int SlepSneu2nd = -1;
    int prod_mode = GetGenCascadesProduction(SlepSneu1st, SlepSneu2nd);
    if(prod_mode == 0) return 0;

    if(prod_mode == 1 || prod_mode == 5) // SlepSlep prods
      if(this->GenPart_pdgId[SlepSneu1st] < 0) // Slep^{-} is always '1st' or 'left' side of tree
         std::swap(SlepSneu1st,SlepSneu2nd);
    if(prod_mode == 2 || prod_mode == 3 || prod_mode == 6 || prod_mode == 7) // SlepSneu prods
      if(abs(this->GenPart_pdgId[SlepSneu1st]) % 2 == 0) // Slep is always '1st' or 'left' side of tree
         std::swap(SlepSneu1st,SlepSneu2nd);

    int SlepSneu1stDecayMode = 0;
    std::pair<int,int> SlepSneu1stDecays = GetGenCascadesDecayMode(SlepSneu1st);
    if(SlepSneu1stDecays.second == 1000022) SlepSneu1stDecayMode = 1; // lep + N1
    else if(SlepSneu1stDecays.second == 1000023) SlepSneu1stDecayMode = 2; // lep + N2
    else if(abs(SlepSneu1stDecays.second) == 1000024) SlepSneu1stDecayMode = 3; // lep + C1
    int SlepSneu2ndDecayMode = 0;
    std::pair<int,int> SlepSneu2ndDecays = GetGenCascadesDecayMode(SlepSneu2nd);
    if(SlepSneu2ndDecays.second == 1000022) SlepSneu2ndDecayMode = 1; // lep + N1
    else if(SlepSneu2ndDecays.second == 1000023) SlepSneu2ndDecayMode = 2; // lep + N2
    else if(abs(SlepSneu2ndDecays.second) == 1000024) SlepSneu2ndDecayMode = 3; // lep + C1
    
    int N2_1st_Index = -1;
    int N2_1stDecayMode = 0;
    if(SlepSneu1stDecayMode == 2){ // need to figure out N2 decay
      N2_1st_Index = GetGenCascadesIndex(SlepSneu1st, 1000023); // get index of N2 coming from 1st Slep/Sneu
      std::pair<int,int> N2_1stDecays = GetGenCascadesDecayMode(N2_1st_Index); // get decays of N2 from 1st Slep/Sneu
      if(N2_1stDecays.first == 22) N2_1stDecayMode = 1; // N2 -> photon + LSP
      else if(N2_1stDecays.first == 23) N2_1stDecayMode = 2; // N2 -> Z + LSP
      else if(N2_1stDecays.second == 1000024) N2_1stDecayMode = 3; // N2 -> W- + C1+
      else if(N2_1stDecays.second == -1000024) N2_1stDecayMode = 4; // N2 -> W+ + C1-
      else if(abs(N2_1stDecays.first) == 24 && abs(N2_1stDecays.second) == 24) N2_1stDecayMode = 5; // N2 -> W+ + W-
    }
    int N2_2nd_Index = -1;
    int N2_2ndDecayMode = 0;
    if(SlepSneu2ndDecayMode == 2){ // need to figure out N2 decay
      N2_2nd_Index = GetGenCascadesIndex(SlepSneu2nd, 1000023); // get index of N2 coming from 2nd Slep/Sneu
      std::pair<int,int> N2_2ndDecays = GetGenCascadesDecayMode(N2_2nd_Index); // get decays of N2 from 2nd Slep/Sneu
      if(N2_2ndDecays.first == 22) N2_2ndDecayMode = 1; // N2 -> photon + LSP
      else if(N2_2ndDecays.first == 23) N2_2ndDecayMode = 2; // N2 -> Z + LSP
      else if(N2_2ndDecays.second == 1000024) N2_2ndDecayMode = 3; // N2 -> W- + C1+
      else if(N2_2ndDecays.second == -1000024) N2_2ndDecayMode = 4; // N2 -> W+ + C1-
      else if(abs(N2_2ndDecays.first) == 24 && abs(N2_2ndDecays.second) == 24) N2_2ndDecayMode = 5; // N2 -> W+ + W-
    }

    uint16_t packed = CascadesTreeEncoder::Encode(prod_mode, SlepSneu1stDecayMode, SlepSneu2ndDecayMode, N2_1stDecayMode, N2_2ndDecayMode);
    return packed;
  }
  return 0;
}

template <class Base>
int AnalysisBase<Base>::GetSampleIndex(){
  if(!m_IsSMS && !m_IsCascades){
    if(m_Nsample == 0){
      m_IndexToSample[0]  = "KUAnalysis";
      m_IndexToXsec[0]    = m_XsecTool.GetXsec_BKG(m_DataSet);
      m_IndexToNevent[0]  = m_NeventTool.GetNevent_BKG(m_DataSet, m_FileTag);
      m_IndexToNweight[0] = m_NeventTool.GetNweight_BKG(m_DataSet, m_FileTag);
      m_Nsample++;
    }
    return 0;
  }

  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_mass)>::value ){
  
    int MP = 0;
    int MC = 0;
    int Ngen = this->nGenPart;
    int PDGID;
    // for cascades as SMS
    int code = 0;
    bool has_Slep = false;
    bool has_Snu = false;
    bool is_left = true;
    bool is_plus = true; // plus referring to charge of e or mu (not value of PDGID)
    for(int i = 0; i < Ngen; i++){
      PDGID = fabs(this->GenPart_pdgId[i]);
      if(PDGID > 1000000 && PDGID < 3000000){
        int mass = int(this->GenPart_mass[i]+0.5);
        if(PDGID == 1000022)
          MC = mass;
        else
          if(mass > MP)
            MP = mass;
      }
      // Getting 'code' for cascades
      if(m_IsCascades && m_IsSMS){
        if (abs(PDGID) % 10000 == 11 || abs(PDGID) % 10000 == 13) {
          has_Slep = true;
          if (PDGID > 0) is_plus = false;
        } 
        else if (abs(PDGID) % 10000 == 12 || abs(PDGID) % 10000 == 14) {
            has_Snu = true;
        }
        if (PDGID > 2000000) is_left = false;
      }
    } // for(int i = 0; i < Ngen; i++)
    if(m_IsCascades && m_IsSMS){
      // build code from booleans
      if(has_Slep && !has_Snu) code += 1; // SlepSlep
      if(has_Slep && has_Snu) code += 2; // SlepSnu
      if(!has_Slep && has_Snu) code += 3; // SnuSnu
      if(is_left) code += 10;
      else code += 20;
      if(is_plus) code += 100;
      else code += 200;
    }
    
    //int hash = 100000*MP + MC;
    long long hash = ((long long)MP << 28) | ((long long)MC << 14) | code;
    if(m_IsCascades){
      if(m_IsSMS){ // should not be needed
        if(m_HashToIndex.count(hash) == 0){
          m_HashToIndex[hash] = m_Nsample;
          m_IndexToSample[m_Nsample]  = std::string(Form("SMS_%d_%d_%d", MP, MC, code));
          m_IndexToXsec[m_Nsample]    = m_XsecTool.GetXsec_SMS_code(m_DataSet, MP, code, m_IsRun3);
          m_IndexToNevent[m_Nsample]  = m_NeventTool.GetNevent_SMS_code(m_DataSet, m_FileTag, MP, MC, code);
          m_IndexToNweight[m_Nsample] = m_NeventTool.GetNweight_SMS_code(m_DataSet, m_FileTag, MP, MC, code);
          m_Nsample++;
        }
        return m_HashToIndex[hash];
      } else {
        if(m_Nsample == 0){
          m_IndexToSample[0]  = "KUAnalysis";
          m_IndexToXsec[0]    = m_XsecTool.GetXsec_Cascades(m_DataSet, MP, m_IsRun3);
          m_IndexToNevent[0]  = m_NeventTool.GetNevent_Cascades(m_DataSet, m_FileTag);
          m_IndexToNweight[0] = m_NeventTool.GetNweight_Cascades(m_DataSet, m_FileTag);
          m_Nsample++;
        }
        return 0.;
      }
    } else {
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
  }
  return 0.;
}

template <class Base>
vector<int> AnalysisBase<Base>::GetLSPParents(){
  std::vector<int> LSPParents;
  if constexpr (
    std::is_member_object_pointer<decltype(&Base::nGenPart)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_pdgId)>::value &&
    std::is_member_object_pointer<decltype(&Base::GenPart_genPartIdxMother)>::value){

    int N = this->nGenPart;
    int PDGID, mom;
    for(int i = 0; i < N; i++){
      if(this->GenPart_pdgId[i] != 1000022) continue;
      mom = this->GenPart_genPartIdxMother[i];
      LSPParents.push_back(this->GenPart_pdgId[mom]);
    }

  }
  return LSPParents;
}

template <class Base>
double AnalysisBase<Base>::GetEventWeight(){
  if(IsData())
    return 1.;
  
  if(m_IndexToNweight[m_SampleIndex] > 0.){
    if constexpr (std::is_member_object_pointer<decltype(&Base::genWeight)>::value && 
                  std::is_member_object_pointer<decltype(&Base::luminosityBlock)>::value){
      double weight = m_IndexToXsec[m_SampleIndex];
      if(m_IsPrivateMC)
        weight *= 1./m_IndexToNevent[m_SampleIndex];
      else 
        weight *= this->genWeight/m_IndexToNweight[m_SampleIndex];
      if(m_IsSMS)
        weight *= m_NeventTool.GetFilterEff(m_DataSet,m_FileTag,this->luminosityBlock);
      if(weight == 0. || weight == 1.){
        cout << "weight " << weight << endl;
        cout << "genWeight " << this->genWeight << endl;
        cout << "Xsec " << m_IndexToXsec[m_SampleIndex] << endl;
        cout << "Nevent " << m_IndexToNevent[m_SampleIndex] << endl;
        cout << "Nweight " << m_IndexToNweight[m_SampleIndex] << endl;
        if(m_IsSMS)
          cout << "Filter eff " << m_NeventTool.GetFilterEff(m_DataSet,m_FileTag,this->luminosityBlock) << endl;
      }
      return weight;

  } } else return 0.;
}

template <class Base>
double AnalysisBase<Base>::GetGenEventWeight(){
  if(IsData())
    return 1.;
  if constexpr (std::is_member_object_pointer<decltype(&Base::genWeight)>::value)
    return this->genWeight;
  else
    return 0;
}

// [0] is muR=0.5 muF=0.5 ; [1] is muR=0.5 muF=1.0 ; [2] is muR=0.5 muF=2.0 ;
// [3] is muR=0.1 muF=0.5 ; [4] is muR=1.0 muF=1.0 ; [5] is muR=1.0 muF=2.0 ;
// [6] is muR=2.0 muF=0.5 ; [7] is muR=2.0 muF=1.0 ; [8] is muR=2.0 muF=2.0 ;

template <class Base>
double AnalysisBase<Base>::GetMuFWeight(int updown){
  if(IsData())
    return 1.;
  if constexpr (std::is_member_object_pointer<decltype(&Base::LHEScaleWeight)>::value &&
                std::is_member_object_pointer<decltype(&Base::nLHEScaleWeight)>::value){
    if(this->nLHEScaleWeight == 0) return 1.;
    if(updown > 0)
      return this->LHEScaleWeight[5];
    else if(updown < 0) 
      return this->LHEScaleWeight[3];
    else
      return this->LHEScaleWeight[4]; //nominal
  }
}

template <class Base>
double AnalysisBase<Base>::GetMuRWeight(int updown){
  if(IsData())
    return 1.;
  if constexpr (std::is_member_object_pointer<decltype(&Base::LHEScaleWeight)>::value &&
                std::is_member_object_pointer<decltype(&Base::nLHEScaleWeight)>::value){
    if(this->nLHEScaleWeight == 0) return 1.;
    if(updown > 0)
      return this->LHEScaleWeight[7];
    else if(updown < 0) 
      return this->LHEScaleWeight[1];
    else
      return this->LHEScaleWeight[4]; //nominal
  }
}

template <class Base>
double AnalysisBase<Base>::GetPUWeight(int updown){
  if(IsData())
    return 1.;
  if constexpr (std::is_member_object_pointer<decltype(&Base::Pileup_nPU)>::value)
    return m_PUTool.GetWeight(this->Pileup_nPU, m_year, updown);
  return 1.;
}

#ifdef _CMSSW_
template <class Base>
void AnalysisBase<Base>::AddLHAPDF(){
  if(IsData() || IsSMS() || IsCascades())
    return;
  m_LHETool.AddLHAPDF(m_year);
}
#endif

template <class Base>
double AnalysisBase<Base>::GetPDFWeight(int updown){
  if(IsData() || IsSMS() || IsCascades())
    return 1.;
  if constexpr (std::is_member_object_pointer<decltype(&Base::nLHEPdfWeight)>::value &&
                std::is_member_object_pointer<decltype(&Base::LHEPdfWeight)>::value &&
                std::is_member_object_pointer<decltype(&Base::Generator_id1)>::value &&
                std::is_member_object_pointer<decltype(&Base::Generator_id2)>::value &&
                std::is_member_object_pointer<decltype(&Base::Generator_x1)>::value &&
                std::is_member_object_pointer<decltype(&Base::Generator_x2)>::value &&
                std::is_member_object_pointer<decltype(&Base::Generator_scalePDF)>::value){
    return m_LHETool.GetWeight(this->nLHEPdfWeight,this->LHEPdfWeight,this->Generator_id1,this->Generator_id2,this->Generator_x1,this->Generator_x2,this->Generator_scalePDF,this->m_year,updown);
  }
  return 1.;
}

/////////////////////////////////////////////////
// End NANOAOD generic methods
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// Start SUSYNANOBase specific methods
/////////////////////////////////////////////////

template<>
double AnalysisBase<SUSYNANOBase>::EGvalue(int jetIndex, int updown){
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
double AnalysisBase<SUSYNANOBase>::GetPrefireWeight(int updown){
  if(m_year == 2018)
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
bool AnalysisBase<SUSYNANOBase>::PassEventFilter(){
  if(m_year == 2016){
    return Flag_goodVertices &&
      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
      Flag_HBHENoiseFilter &&
      Flag_HBHENoiseIsoFilter &&
      Flag_EcalDeadCellTriggerPrimitiveFilter &&
      Flag_BadPFMuonFilter;
  }
  if(m_year == 2017){
    return Flag_goodVertices &&
      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
      Flag_HBHENoiseFilter &&
      Flag_HBHENoiseIsoFilter &&
      Flag_EcalDeadCellTriggerPrimitiveFilter &&
      Flag_BadPFMuonFilter;
  }
  if(m_year == 2018){
    return Flag_goodVertices &&
      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
      Flag_HBHENoiseFilter &&
      Flag_HBHENoiseIsoFilter &&
      Flag_EcalDeadCellTriggerPrimitiveFilter &&
      Flag_BadPFMuonFilter;
  }
  return true;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetMETtrigger(){
  if(m_year == 2016)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  if(m_year == 2017 || m_year == 2018)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  return 0;
}
  
template <>
bool AnalysisBase<SUSYNANOBase>::GetMETORtrigger(){
  if(m_year == 2016)
    return (HLT_PFMETNoMu90_PFMHTNoMu90_IDTight ||
	    HLT_PFMETNoMu100_PFMHTNoMu100_IDTight ||
	    //HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET90_PFMHT90_IDTight ||
	    HLT_PFMET100_PFMHT100_IDTight ||
	    HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned ||
	    //HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight);
  if(m_year == 2017 || m_year == 2018)
    return (//HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMET130_PFMHT130_IDTight ||
	    HLT_PFMET140_PFMHT140_IDTight ||
	    //HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMETNoMu130_PFMHTNoMu130_IDTight ||
	    HLT_PFMETNoMu140_PFMHTNoMu140_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetMETDoubleMutrigger(){
  if(m_year == 2016)
     return HLT_DoubleMu3_PFMET50;
  else return HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetSingleElectrontrigger(){
  if(m_year == 2016)
    return (HLT_Ele27_WPTight_Gsf);
  if(m_year == 2017 ||
     m_year == 2018)
    return (HLT_Ele32_WPTight_Gsf);
  return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetSingleMuontrigger(){
  if(m_year == 2016 ||
     m_year == 2017 ||
     m_year == 2018  )
    return (HLT_IsoMu24 ||
            HLT_IsoTkMu24);
  return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetDoubleElectrontrigger(){
  if(m_year == 2016 ||
     m_year == 2017 ||
     m_year == 2018  )
    return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetDoubleMuontrigger(){
  if(m_year == 2016 ||
     m_year == 2017 ||
     m_year == 2018  )
    return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetTripleElectrontrigger(){
  return HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetTripleMuonLowPTtrigger(){
  if(m_year == 2016)
    return HLT_TripleMu_5_3_3_DZ_Mass3p8;
  else if(m_year == 2017)
    return (HLT_TripleMu_5_3_3_Mass3p8to60_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8to60_DCA
           );
  else if(m_year == 2018)
    return (HLT_TripleMu_5_3_3_Mass3p8to60_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8to60_DCA ||
            HLT_TripleMu_5_3_3_Mass3p8_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8_DCA
           );
  else return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetTripleMuonHighPTtrigger(){
  if(m_year == 2016)
    return HLT_TripleMu_12_10_5;
  else if(m_year == 2017)
    return (HLT_TripleMu_12_10_5 ||
            HLT_TripleMu_10_5_5_DZ
           );
  else if(m_year == 2018)
    return (HLT_TripleMu_12_10_5 ||
            HLT_TripleMu_10_5_5_DZ
           );
  else return 0;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetDiMuEleLowPTtrigger(){
  if(m_year == 2018) return HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
  else return false;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetDiMuEleHighPTtrigger(){
  if(m_year == 2016) return HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  else return HLT_DiMu9_Ele9_CaloIdL_TrackIdL || HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetDiEleMutrigger(){
  if(m_year == 2016) return HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  else return HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ || HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetEMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetEMuMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<SUSYNANOBase>::GetEMuEtrigger(){
  return HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
double AnalysisBase<SUSYNANOBase>::GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown, ParticleIDType tag){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

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
    
    EFF = m_BtagSFTool.EFF(jets[i].Pt(), m_year, iflavor, FastSim);
    SF  = m_BtagSFTool.SF(jets[i].Pt(), m_year, iflavor, updown);
    if(FastSim)
      SF *= m_BtagSFTool.SF(jets[i].Pt(), m_year, iflavor, updown, FastSim);

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
double AnalysisBase<SUSYNANOBase>::GetElIDSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetElISOSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetElSIPSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetElVLIDSFWeight(const ParticleList& els, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int m_year = 2016;
  if(m_FileTag.find("17") != std::string::npos)
    m_year = 2017;
  if(m_FileTag.find("18") != std::string::npos)
    m_year = 2018;

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetMuIDSFWeight(const ParticleList& mus, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetMuISOSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetMuSIPSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetMuVLIDSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<SUSYNANOBase>::GetMETTriggerSFWeight(double MET, double HT, int Nele, int Nmu, int updown){
  if(IsData())
    return 1.;

  if(IsFastSim()){
    return m_METTriggerTool.Get_EFF(MET, HT, m_year,
				    (Nele > 0), (Nmu > 0),
				    false, updown)*
      m_METTriggerTool.Get_SF(MET, HT, m_year,
			      (Nele > 0), (Nmu > 0),
			      false, updown);
  } else {
    return m_METTriggerTool.Get_SF(MET, HT, m_year,
				   (Nele > 0), (Nmu > 0),
				   false, updown);
  }
}

template <>
int AnalysisBase<SUSYNANOBase>::GetMETTriggerSFCurve(double HT, int Nele, int Nmu){
  return m_METTriggerTool.Get_Curve_Index(HT, m_year, (Nele > 0), (Nmu > 0), IsData());
}

template <>
void AnalysisBase<SUSYNANOBase>::InitializeHistograms(vector<TH1D*>& histos){
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
void AnalysisBase<SUSYNANOBase>::BookHistograms(vector<TH1D*>& histos){
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
    if(m_year == 2016)
      if(Jet_btagDeepFlavB[i] > 0.3093)
	btag = true;
    if(m_year == 2017)
      if(Jet_btagDeepFlavB[i] > 0.3033)
	btag = true;
    if(m_year == 2018)
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
void AnalysisBase<SUSYNANOBase>::ApplyMETPhiCorrections(TLorentzVector &metOut, int runNumber, int npvGood, const std::string &variation) {}

template <>
ParticleList AnalysisBase<SUSYNANOBase>::GetJetsMET(TVector3& MET, int id){
  ParticleList list;
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
      double uncer = m_JMETool.GetJESFactor(m_year, CurrentSystematic().Label(),
					    Jet_pt[i], Jet_eta[i]);
      
      deltaMET -= delta*uncer*jet.Vect();
      jet.SetPtEtaPhiM((1.+delta*uncer)*jet.Pt(),
		       jet.Eta(), jet.Phi(),
		       (1.+delta*uncer)*jet.M());
    }
    
    if(!IsData()){

      // JER recipe based on https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
     
      double smearFactor = 1.;
      double JER = m_JMETool.GetJERFactor(m_year, Jet_pt[i], Jet_eta[i], fixedGridRhoFastjetAll); // using this for rho based on: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/0127d46a973e894d97e9a16bd3939f421b2b689e/python/postprocessing/modules/jme/jetmetUncertainties.py#L49
      double SF = m_JMETool.GetJERSFFactor(m_year,Jet_eta[i],0);

      if(DO_JER)
        SF = m_JMETool.GetJERSFFactor(m_year,Jet_eta[i],delta);
       
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
    if(m_year == 2016){
      // Deep Flavor
      if(jet.Btag() > 0.7221)
	jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.3093) 
	jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0614)
	jet.SetBtagID(kLoose);
    }

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    if(m_year == 2017){
      // Deep Flavor
      if(jet.Btag() > 0.7489)
	jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.3033) 
	jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0521)
	jet.SetBtagID(kLoose);
    }

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    if(m_year == 2018){
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
  
  if(m_year == 2017)
    MET.SetPtEtaPhi(METFixEE2017_pt,0.0,METFixEE2017_phi);
  else
    MET.SetPtEtaPhi(MET_pt,0.0,MET_phi);
  
  deltaMET.SetZ(0.);
  MET += deltaMET;
  
  if(CurrentSystematic() == Systematic("METUncer_UnClust")){
    if(m_year == 2017)
      deltaMET.SetXYZ(delta*METFixEE2017_MetUnclustEnUpDeltaX,
        	      delta*METFixEE2017_MetUnclustEnUpDeltaY, 0.);
    else
      deltaMET.SetXYZ(delta*MET_MetUnclustEnUpDeltaX,
		      delta*MET_MetUnclustEnUpDeltaY, 0.);
    MET += deltaMET;
  }

  if(CurrentSystematic() == Systematic("METUncer_GenMET"))
    MET.SetPtEtaPhi(GenMET_pt,0.,GenMET_phi);
  
  return list;
}

template <>
TVector3 AnalysisBase<SUSYNANOBase>::GetMET(){
  TVector3 MET;
  GetJetsMET(MET);

  return MET;
}

template <>
TVector3 AnalysisBase<SUSYNANOBase>::GetAltMET(){
  TVector3 MET;
  MET.SetPtEtaPhi(MET_pt,0.0,MET_phi);

  return MET;
}

template <>
ParticleList AnalysisBase<SUSYNANOBase>::GetJets(int id){
  TVector3 dum;
  return GetJetsMET(dum, id);
}

template <>
ParticleList AnalysisBase<SUSYNANOBase>::GetElectrons(){
  ParticleList list;

  int N = nElectron;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Electron_pt[i] < 5. || fabs(Electron_eta[i]) > 2.5)
      continue;
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
    lep.SetTightCharge(Electron_tightCharge[i]);

    lep.SetRelIso(Electron_pfRelIso03_all[i]);
    lep.SetMiniIso(Electron_miniPFRelIso_all[i]);

    // https://twiki.cern.ch/twiki/pub/CMS/SUSLeptonSF/Run2_SUSYwp_EleCB_MVA_8Jan19.pdf
    
    // FO baseline criteria
    if(Electron_lostHits[i] == 0 && Electron_convVeto[i]){

      double mva = Electron_mvaFall17V1noIso[i];
      if(m_year == 2016 || m_year == 2018)
	mva = Electron_mvaFall17V2noIso[i];


      // convert to raw MVA output
      if(mva == -1.)
	mva = -999.;
      else if(mva == 1.)
	mva = 999.;
      else
	mva = -0.5*log((1.-mva)/(1.+mva));
      
      // FO VLoose
      if(m_year == 2016){ // Summer16_94X legacy
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

      if(m_year == 2017){ // Fall17_94X
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

      if(m_year == 2018){ // Autumn18_102X
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
      if(m_year == 2016){ // Summer16_94X legacy
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

      if(m_year == 2017){ // Fall17_94X
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

      if(m_year == 2018){ // Autumn18_102X
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
	if(m_year == 2016){ // Summer16_94X legacy
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

	if(m_year == 2017){ // Fall17_94X
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

	if(m_year == 2018){ // Autumn18_102X
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
ParticleList AnalysisBase<SUSYNANOBase>::GetMuons(){
  ParticleList list;

  int N = nMuon;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Muon_pt[i] < 3. || fabs(Muon_eta[i]) > 2.4)
      continue;
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
    lep.SetTightCharge(Muon_tightCharge[i]);

    lep.SetRelIso(Muon_pfRelIso03_all[i]);
    lep.SetMiniIso(Muon_miniPFRelIso_all[i]);

    // FO baseline criteria
    lep.SetParticleID(kLoose);

    if(Muon_tightId[i])
      lep.SetParticleID(kTight);
    else if(lep.Pt() < 0.){
      if(Muon_softId[i])
        lep.SetParticleID(kMedium);
    } else {
      if(Muon_mediumId[i])
        lep.SetParticleID(kMedium);
    }
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

/////////////////////////////////////////////////
// End SUSYNANOBase specific methods
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// Start NANOULBase specific methods
/////////////////////////////////////////////////

template <>
bool AnalysisBase<NANOULBase>::PassEventFilter(){
  if(m_year == 2016){
    return Flag_goodVertices &&
      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
      Flag_HBHENoiseFilter &&
      Flag_HBHENoiseIsoFilter &&
      Flag_EcalDeadCellTriggerPrimitiveFilter &&
      Flag_BadPFMuonFilter &&
      Flag_BadPFMuonDzFilter &&
      Flag_eeBadScFilter;
  }
  else if(m_year == 2017 || m_year == 2018){
    return Flag_goodVertices &&
      (IsFastSim() ? true : Flag_globalSuperTightHalo2016Filter) &&
      Flag_HBHENoiseFilter &&
      Flag_HBHENoiseIsoFilter &&
      Flag_EcalDeadCellTriggerPrimitiveFilter &&
      Flag_BadPFMuonFilter &&
      Flag_BadPFMuonDzFilter &&
      Flag_eeBadScFilter && 
      Flag_ecalBadCalibFilter;
  }
  return true;
}

template<>
double AnalysisBase<NANOULBase>::EGvalue(int jetIndex, int updown){
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

template <>
bool AnalysisBase<NANOULBase>::GetMETtrigger(){
  if(m_year == 2016)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
  if(m_year == 2017 || m_year == 2018)
    return (HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);

  return 0;
}

template <>
bool AnalysisBase<NANOULBase>::GetMETORtrigger(){
  if(m_year == 2016)
    return (HLT_PFMETNoMu90_PFMHTNoMu90_IDTight ||
	    HLT_PFMETNoMu100_PFMHTNoMu100_IDTight ||
	    //HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMET90_PFMHT90_IDTight ||
	    HLT_PFMET100_PFMHT100_IDTight ||
	    //HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight);
  if(m_year == 2017 ||
     m_year == 2018)
    return (//HLT_PFMET110_PFMHT110_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight ||
	    HLT_PFMET130_PFMHT130_IDTight ||
	    HLT_PFMET140_PFMHT140_IDTight ||
	    //HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
	    HLT_PFMETNoMu130_PFMHTNoMu130_IDTight ||
	    HLT_PFMETNoMu140_PFMHTNoMu140_IDTight ||
	    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
	    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
  return 0;
}

template <>
bool AnalysisBase<NANOULBase>::GetMETDoubleMutrigger(){
  if(m_year == 2016)
     return HLT_DoubleMu3_PFMET50;
  else return HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
}

template <>
bool AnalysisBase<NANOULBase>::GetSingleElectrontrigger(){
  if(m_year == 2016)
    return (HLT_Ele27_WPTight_Gsf);
  if(m_year == 2017)
    return (HLT_Ele35_WPTight_Gsf);
  if(m_year == 2018)
    return (HLT_Ele32_WPTight_Gsf);
  return 0;
}

template <>
bool AnalysisBase<NANOULBase>::GetSingleMuontrigger(){
  return HLT_IsoMu24;
}

template <>
bool AnalysisBase<NANOULBase>::GetDoubleElectrontrigger(){
  return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
}

template <>
bool AnalysisBase<NANOULBase>::GetDoubleMuontrigger(){
  return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
}

template <>
bool AnalysisBase<NANOULBase>::GetTripleElectrontrigger(){
  return HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<NANOULBase>::GetTripleMuonLowPTtrigger(){
  if(m_year == 2016)
    return HLT_TripleMu_5_3_3_DZ_Mass3p8;
  else if(m_year == 2017)
    return (HLT_TripleMu_5_3_3_Mass3p8to60_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8to60_DCA
           );
  else if(m_year == 2018)
    return (HLT_TripleMu_5_3_3_Mass3p8to60_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8to60_DCA ||
            HLT_TripleMu_5_3_3_Mass3p8_DZ ||
            HLT_TripleMu_5_3_3_Mass3p8_DCA
           );
  else return 0;
}

template <>
bool AnalysisBase<NANOULBase>::GetTripleMuonHighPTtrigger(){
  if(m_year == 2016)
    return HLT_TripleMu_12_10_5;
  else if(m_year == 2017)
    return (HLT_TripleMu_12_10_5 ||
            HLT_TripleMu_10_5_5_DZ
           );
  else if(m_year == 2018)
    return (HLT_TripleMu_12_10_5 ||
            HLT_TripleMu_10_5_5_DZ
           );
  else return 0;
}

template <>
bool AnalysisBase<NANOULBase>::GetDiMuEleLowPTtrigger(){
  if(m_year == 2018) return HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
  else return false;
}

template <>
bool AnalysisBase<NANOULBase>::GetDiMuEleHighPTtrigger(){
  if(m_year == 2016) return HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  else return HLT_DiMu9_Ele9_CaloIdL_TrackIdL || HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
}

template <>
bool AnalysisBase<NANOULBase>::GetDiEleMutrigger(){
  if(m_year == 2016) return HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  else return HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ || HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<NANOULBase>::GetEMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<NANOULBase>::GetEMuMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<NANOULBase>::GetEMuEtrigger(){
  return HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>  
double AnalysisBase<NANOULBase>::GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown, ParticleIDType tag) {
  if (IsData()) 
      return 1.;

  bool FastSim = IsFastSim();
  int Njet = jets.size();
  int iflavor = 0;
  double probMC = 1.;
  double probDATA = 1.;
  std::string syst = "central";
  if(updown > 0) syst = "up";
  else if(updown < 0) syst = "down";
  
  for (int i = 0; i < Njet; i++) {
      if(abs(jets[i].PDGID()) == 5)
        iflavor = 5;
      else if(abs(jets[i].PDGID()) == 4)
        iflavor = 4;
      if(HForLF && iflavor == 0)
        continue;
      if(!HForLF && iflavor != 0)
        continue;
      std::vector<std::variant<int, double, std::string>> evalArgs;
      evalArgs.push_back(syst);
      evalArgs.push_back("M"); // Working Point ('M' for medium)
      evalArgs.push_back(iflavor);
      evalArgs.push_back(abs(jets[i].Eta()));
      evalArgs.push_back(jets[i].Pt());

      double SF = 1.;
      double EFF = 1.; // need to measure the efficiencies
      if(iflavor == 0)
        SF = m_cset_Btag->at("deepJet_incl")->evaluate(evalArgs);
      else
        SF = m_cset_Btag->at("deepJet_comb")->evaluate(evalArgs);

      if (jets[i].BtagID() >= tag) {
          probMC *= EFF;
          probDATA *= SF * EFF;
      } else {
          probMC *= (1. - EFF);
          probDATA *= (1. - SF * EFF);
      }
  }

  if (probMC <= 0. || probDATA <= 0.)
      return 1.;

  return probDATA / probMC;
}

template <>
double AnalysisBase<NANOULBase>::GetElIDSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetElISOSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetElSIPSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetElVLIDSFWeight(const ParticleList& els, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetMuIDSFWeight(const ParticleList& mus, int updown){
   if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetMuISOSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetMuSIPSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetMuVLIDSFWeight(const ParticleList& mus, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
double AnalysisBase<NANOULBase>::GetMETTriggerSFWeight(double MET, double HT, int Nele, int Nmu, int updown){
  if(IsData())
    return 1.;

  if(IsFastSim()){
    return m_METTriggerTool.Get_EFF(MET, HT, m_year,
				    (Nele > 0), (Nmu > 0),
				    false, updown)*
      m_METTriggerTool.Get_SF(MET, HT, m_year,
			      (Nele > 0), (Nmu > 0),
			      false, updown);
  } else {
    return m_METTriggerTool.Get_SF(MET, HT, m_year,
				   (Nele > 0), (Nmu > 0),
				   false, updown);
  }
}

template <>
int AnalysisBase<NANOULBase>::GetMETTriggerSFCurve(double HT, int Nele, int Nmu){
  return m_METTriggerTool.Get_Curve_Index(HT, m_year, (Nele > 0), (Nmu > 0), IsData());
}

template <>
void AnalysisBase<NANOULBase>::InitializeHistograms(vector<TH1D*>& histos){
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
void AnalysisBase<NANOULBase>::BookHistograms(vector<TH1D*>& histos){
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
    if(m_year == 2016)
      if(Jet_btagDeepFlavB[i] > 0.3093)
	btag = true;
    if(m_year == 2017)
      if(Jet_btagDeepFlavB[i] > 0.3033)
	btag = true;
    if(m_year == 2018)
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
void AnalysisBase<NANOULBase>::ApplyMETPhiCorrections(TLorentzVector &metOut, int runNumber, int npvGood, const std::string &variation){
  if (!m_cset_METPhi) {
    std::cerr << "[ApplyMETPhiCorrections] ERROR: m_cset_METPhi not loaded.\n";
    return;
  }

  // read the inputs from the outgoing TLorentzVector (should be Type-1 propagated)
  double in_pt  = metOut.Pt();
  double in_phi = metOut.Phi();

  // helper to check key presence using at() (CorrectionSet throws out_of_range if missing)
  auto has_key = [&](const std::string& name) -> bool {
    try {
      m_cset_METPhi->at(name);
      return true;
    } catch (const std::out_of_range&) {
      return false;
    }
  };

  try {
    std::string metLabel = "pfmet";
    std::string dt = m_IsData ? "data" : "mc"; // assume you have m_IsData member
    std::string phiKey = "phi_metphicorr_" + metLabel + "_" + dt;
    std::string ptKey  = "pt_metphicorr_"  + metLabel + "_" + dt;

    // Build args: met_pt (real), met_phi (real), npvs (real/int), run (real/int)
    std::vector<std::variant<int,double,std::string>> args;
    args.push_back(in_pt);
    args.push_back(in_phi);
    args.push_back(static_cast<double>(npvGood));
    args.push_back(static_cast<double>(runNumber));

    double corrected_phi = in_phi;
    double corrected_pt  = in_pt;

    try {
      // access + evaluate phi
      if (!has_key(phiKey)) { 
        std::cerr << "[ApplyMETPhiCorrections] WARNING: Run-2 phi key '" << phiKey
                  << "' not found in MET cset. Skipping phi correction.\n";
      } else {
        corrected_phi = m_cset_METPhi->at(phiKey)->evaluate(args);
      }    
    } catch (const std::exception &e) {
      std::cerr << "[ApplyMETPhiCorrections] ERROR evaluating phi key '" << phiKey
                << "': " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ApplyMETPhiCorrections] Unknown error evaluating phi key '"
                << phiKey << "'.\n";
    }    

    try {
      if (!has_key(ptKey)) { 
        // some Run-2 configs do not provide pt-correction; that's fine
      } else {
        corrected_pt = m_cset_METPhi->at(ptKey)->evaluate(args);
      }    
    } catch (const std::exception &e) {
      std::cerr << "[ApplyMETPhiCorrections] ERROR evaluating pt key '" << ptKey
                << "': " << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ApplyMETPhiCorrections] Unknown error evaluating pt key '"
                << ptKey << "'.\n";
    }    

    wrapPhi(corrected_phi);
    metOut.SetPtEtaPhiM(corrected_pt, 0.0, corrected_phi, 0.); 
    return;

  } catch (const std::exception &e) {
    std::cerr << "[ApplyMETPhiCorrections] Caught exception: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ApplyMETPhiCorrections] Caught unknown exception\n";
  }
}

template <>
ParticleList AnalysisBase<NANOULBase>::GetJetsMET(TVector3& MET, int id){
  ParticleList list;
  int Njet = nJet;

  // Determine whether current Systematic is a JES/JER
  bool DO_JES = (m_SysTool.JESSystematics() == CurrentSystematic());
  bool DO_JER = (m_SysTool.JERSystematics() == CurrentSystematic());

  // Prepare vectors for MET propagation (RAW AK4)
  const int nAk4Jets = std::max(int(nJet), 0);
  std::vector<JecApplication::JetForMet> jetsForMet;
  jetsForMet.reserve(static_cast<size_t>(nAk4Jets));
  std::vector<JecApplication::JerInputs> jersForMet;
  jersForMet.reserve(static_cast<size_t>(nAk4Jets));
  long long EventNum = static_cast<long long>(GetEventNum());
  int RunNum = static_cast<int>(GetRunNum());
  int LumiBlock = static_cast<int>(GetLumiNum());

  for (int j = 0; j < nAk4Jets; ++j) {
    JecApplication::JetForMet v;
    v.phi  = Jet_phi[j];
    v.eta  = Jet_eta[j];
    v.area = Jet_area[j];
    v.rawPt = (1.0 - Jet_rawFactor[j]) * Jet_pt[j]; // reconstruct RAW from NanoAOD
    v.muonSubtrFactor = Jet_muonSubtrFactor[j];
    v.chEmEf = Jet_chEmEF[j];
    v.neEmEf = Jet_neEmEF[j];
    jetsForMet.push_back(v);

    // JER inputs for MET propagation
    JecApplication::JerInputs jerIn;
    jerIn.event = EventNum;
    jerIn.run = RunNum;
    jerIn.lumi = LumiBlock;
    const int genIdx = Jet_genJetIdx[j];
    if (genIdx > -1 && static_cast<UInt_t>(genIdx) < static_cast<UInt_t>(nGenJet)) {
      jerIn.hasGen = true;
      jerIn.genPt  = GenJet_pt[genIdx];
      jerIn.genEta = GenJet_eta[genIdx];
      jerIn.genPhi = GenJet_phi[genIdx];
      jerIn.maxDr  = 0.2; // AK4 match radius
    } else {
      jerIn.hasGen = false;
    }
    jersForMet.push_back(jerIn);
  }

  // Build SystematicOptions for MET helper (JES & JER)
  JecApplication::SystematicOptions systOpts{};
  if (!m_IsData) {
    // JES: resolve CLIB key (single resolution for MET propagation)
    if (DO_JES) {
      const std::string systLabel = CurrentSystematic().Label();
      static std::map<std::string, std::optional<std::string>> jesKeyCache;
      const std::string cacheKey = m_jmeYearKey + "|" + systLabel;
      std::optional<std::string> optClibKey;
      auto it = jesKeyCache.find(cacheKey);
      if (it != jesKeyCache.end()) {
        optClibKey = it->second;
      } else {
        optClibKey = findJESKeyForLabel(systLabel);
        jesKeyCache.emplace(cacheKey, optClibKey);
        if (!optClibKey)
          std::cerr << "[JME] Warning: could not resolve JES systematic '" << systLabel
                    << "' for year key " << m_jmeYearKey << ". MET will use nominal JES.\n";
      }

      if (optClibKey) {
        systOpts.jesSystName = *optClibKey;
        systOpts.jesSystVar  = (CurrentSystematic().IsUp() ? "Up" : "Down");
      } else {
        systOpts.jesSystName.clear();
        systOpts.jesSystVar.clear();
      }
    }

    if (DO_JER) systOpts.jerVar = (CurrentSystematic().IsUp() ? "up" : "down");
    else        systOpts.jerVar = "nom";
  }

  // --- Process jets for analysis (apply nominal JES, JES syst, JER) ---
  for(int i = 0; i < Njet; ++i){
    Particle jet;
    float mass = Jet_mass[i];
    if(std::isnan(mass) || std::isinf(mass) || mass < 0.) mass = 0.;

    // --- Nominal JES ---
    double ptRaw = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
    const double jesNominal = m_JMETool.GetJESFactorCLIB(
      m_jmeYearKey, m_IsData, RunNum, ptRaw,
      Jet_eta[i], Jet_phi[i], Jet_area[i], fixedGridRhoFastjetAll
    );
    jet.SetPtEtaPhiM(ptRaw * jesNominal, Jet_eta[i], Jet_phi[i], mass);

    // --- JES syst (per-jet multiplicative) ---
    if (DO_JES && !m_IsData && !systOpts.jesSystName.empty()) {
      const std::string var = CurrentSystematic().IsUp() ? "Up" : "Down";
      double jesSystFactor = m_JMETool.GetJESFactorSystCLIB(m_jmeYearKey, jet.Pt(), Jet_eta[i], Jet_phi[i], Jet_area[i], fixedGridRhoFastjetAll, systOpts.jesSystName, var);
      jet.SetPtEtaPhiM(jet.Pt() * jesSystFactor, jet.Eta(), jet.Phi(), mass);
    }

    // --- Jet selection ---
    if(Jet_pt[i] < 15. || fabs(Jet_eta[i]) > 5.) continue;
    //if(Jet_jetId[i] < id) continue;
    // Can't use jetId in NanoAODv12: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV#nanoAOD_Flags
    if (!(Jet_jetId[i] & (1 << 1)) || Jet_muEF[i] >= 0.8 || Jet_chEmEF[i] >= 0.8)
      continue;

    if(Jet_jetId[i] >= 3) jet.SetParticleID(kTight);
    else if(Jet_jetId[i] >= 2) jet.SetParticleID(kMedium);
    else if(Jet_jetId[i] >= 1) jet.SetParticleID(kLoose);

    // --- B-tagging ---
    jet.SetBtag(Jet_btagDeepFlavB[i]);
    if(jet.Btag() > m_BtagTightWP) jet.SetBtagID(kTight);
    else if(jet.Btag() > m_BtagMediumWP) jet.SetBtagID(kMedium);
    else if(jet.Btag() > m_BtagLooseWP) jet.SetBtagID(kLoose);

    jet.SetPDGID(Jet_partonFlavour[i]);

    // --- JER: apply per-jet (nominal or up/down) via JMETool wrapper that uses JecApplication ---
    if(!m_IsData) {
      JecApplication::JesInputs jAfter{jet.Pt(), Jet_eta[i], Jet_phi[i], Jet_area[i], fixedGridRhoFastjetAll, 0.0};

      JecApplication::JerInputs jerIn = jersForMet[i];
      std::string jerVar = (DO_JER ? (CurrentSystematic().IsUp() ? "up" : "down") : "nom");
      const double sJer = m_JMETool.GetJERScaleFactor(m_jmeYearKey, jAfter, jerIn, jerVar, m_IsData);

      jet.SetPtEtaPhiM(jet.Pt() * sJer, jet.Eta(), jet.Phi(), jet.M() * sJer);
      jet.SetjetID(Jet_jetId[i]);
      jet.SetChEmEF(Jet_chEmEF[i]);
      jet.SetNeEmEF(Jet_neEmEF[i]);
    }

    list.push_back(jet);
  } // end jet loop

  // Call helper to get corrected MET 4-vector (JES+JER applied, MET propagation)
  TLorentzVector p4CorrectedMET = m_JMETool.GetCorrectedMET(m_jmeYearKey, m_IsData, RunNum, RawMET_pt, RawMET_phi, jetsForMet, jersForMet, fixedGridRhoFastjetAll, systOpts, false);

  int NPV = static_cast<int>(GetNPV());
  ApplyMETPhiCorrections(p4CorrectedMET, RunNum, NPV, "nom");

  MET.SetPtEtaPhi(p4CorrectedMET.Pt(), 0.0, p4CorrectedMET.Phi());

  // Unclustered MET systematic
  if(CurrentSystematic() == Systematic("METUncer_UnClust")){
    const double delta = (CurrentSystematic().IsUp() ? 1. : -1.);
    TVector3 dUncl(delta*MET_MetUnclustEnUpDeltaX, delta*MET_MetUnclustEnUpDeltaY, 0.);
    MET += dUncl;
  }

  if(CurrentSystematic() == Systematic("METUncer_GenMET"))
    MET.SetPtEtaPhi(GenMET_pt,0.,GenMET_phi);

  return list;
}

template <>
TVector3 AnalysisBase<NANOULBase>::GetMET(){
  TVector3 MET;
  GetJetsMET(MET);

  return MET;
}

template <>
TVector3 AnalysisBase<NANOULBase>::GetAltMET(){
  TVector3 MET;
  MET.SetPtEtaPhi(MET_pt,0.0,MET_phi);

  return MET;
}

template <>
ParticleList AnalysisBase<NANOULBase>::GetJets(int id){
  TVector3 dum;
  return GetJetsMET(dum, id);
}


template <>
ParticleList AnalysisBase<NANOULBase>::GetElectrons(){
  ParticleList list;

  int N = nElectron;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Electron_pt[i] < 7. || fabs(Electron_eta[i]) > 2.5)
      continue;
    if(Electron_pt[i] < 10. && fabs(Electron_eta[i]) > 1.4442) // remove low pt endcap ele
      continue;
    if(fabs(Electron_eta[i]) >= 1.4442 && fabs(Electron_eta[i]) <= 1.566)
      continue;
    if(fabs(Electron_dxy[i]) >= 0.05 || fabs(Electron_dz[i]) >= 0.1 ||
       Electron_sip3d[i] >= 6.)
      continue;
    if(Electron_pfRelIso03_all[i]*Electron_pt[i] >= 20. + 300./Electron_pt[i])
      continue;
    if(!minus_iso_hoe(Electron_vidNestedWPBitmap[i], 2, std::greater_equal<int>()))
      continue;
    if(Electron_lostHits[i] != 0)
      continue;
    if(Electron_convVeto[i] == 0)
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
    lep.SetTightCharge(Electron_tightCharge[i]);
    lep.SetParticleID(kVeryLoose); // need to set to something

    lep.SetRelIso(Electron_pfRelIso03_all[i]);
    lep.SetMiniIso(Electron_miniPFRelIso_all[i]);
    if( ( Electron_pt[i] < 20. && (
          lep.MiniIso()*lep.Pt() >= 4. || 
          lep.RelIso()*lep.Pt() >= 4. ||
          !minus_iso_hoe(Electron_vidNestedWPBitmap[i], 4, std::greater_equal<int>())
        )) ||
        ( Electron_pt[i] >= 20. && (
          !Electron_mvaFall17V2Iso_WP90[i]
        ))
      ) 
        lep.SetLepQual(kBronze);
    else if(( Electron_pt[i] < 20. && (
          lep.MiniIso()*lep.Pt() < 4. && 
          lep.RelIso()*lep.Pt() < 4. &&
          minus_iso_hoe(Electron_vidNestedWPBitmap[i], 4, std::greater_equal<int>())
        )) ||
        ( Electron_pt[i] >= 20. && (
          Electron_mvaFall17V2Iso_WP90[i]
        ))
      )
    { 
      if(lep.SIP3D() >= 3.)
        lep.SetLepQual(kSilver);
      else
        lep.SetLepQual(kGold);
    }
    list.push_back(lep);

  } // for(int i = 0; i < N; i++)
  return list;
}

template <>
ParticleList AnalysisBase<NANOULBase>::GetLowPtElectrons(){
  ParticleList list;
  int N1 = nLowPtElectron;
  for(int i = 0; i < N1; i++){
    // baseline lepton definition
    if(LowPtElectron_pt[i] < 2. || LowPtElectron_pt[i] >= 7. || fabs(LowPtElectron_eta[i]) > 2.5)
      continue;
    if(fabs(LowPtElectron_eta[i]) >= 1.4442 && fabs(LowPtElectron_eta[i]) <= 1.566)
      continue;
    if(LowPtElectron_pt[i] < 10. && fabs(LowPtElectron_eta[i]) > 1.4442) // remove low pt endcap ele
      continue;
    if(LowPtElectron_convVeto[i] == 0)
      continue;
    if(fabs(LowPtElectron_dxy[i]) >= 0.05 || fabs(LowPtElectron_dz[i]) >= 0.1)
      continue;
    if(LowPtElectron_ID[i] < 2.)
      continue;
    if(LowPtElectron_lostHits[i] != 0)
      continue;
    
    if(LowPtElectron_dxyErr[i] < 1.e-8 || LowPtElectron_dzErr[i] < 1.e-8)
      continue;
    
    // Calculate IP_3D and SIP_3D = IP_3D / IP_3D_err for the LowPtElectron collection.
    // - IP_3D:     3D impact parameter wrt first PV, in cm
    // - SIP_3D:    3D impact parameter significance wrt first PV, in cm
    float dxy       = LowPtElectron_dxy[i];
    float dz        = LowPtElectron_dz[i];
    float dxy_err   = LowPtElectron_dxyErr[i];
    float dz_err    = LowPtElectron_dzErr[i];
    // the SIP_3D defined below is coming from Suyash's talk, its something we defined.
    // The sigmas are needed for that calculation.
    float sigma_xy  = dxy/dxy_err;
    float sigma_z   = dz/dz_err;
    float IP_3D     = sqrt(dxy*dxy + dz*dz);
    float SIP_3D    = sqrt(sigma_xy*sigma_xy + sigma_z*sigma_z);

    if (SIP_3D >= 6.)
      continue;

    if(LowPtElectron_miniPFRelIso_all[i]*LowPtElectron_pt[i] >= 20. + 300./LowPtElectron_pt[i])
     continue;

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
    lep.SetTightCharge(0);

    lep.SetRelIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetMiniIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetParticleID(kVeryLoose); // need to set to something for later on
    
    float pt = lep.Pt();
    float absEta = std::abs(lep.Eta());
    float id = LowPtElectron_ID[i];

    bool failsIso = lep.MiniIso() * pt >= 4. || lep.RelIso() * pt >= 4.;
    if (failsIso) {
      lep.SetLepQual(kBronze);
    } else {
      bool passesID = false;

      if (pt >= 2. && pt < 4.) {
        if (absEta >= 0.8 && absEta < 1.442)
          passesID = id >= 3.;
        else if (absEta < 0.8)
          passesID = id >= 2.6;
      } else if (pt >= 4. && pt < 7.) {
        if (absEta >= 0.8 && absEta < 1.442)
          passesID = id >= 3.2;
        else if (absEta < 0.8)
          passesID = id >= 2.8;
      }

      if (!passesID) {
        lep.SetLepQual(kBronze);
      } else if (lep.SIP3D() >= 3.) {
        lep.SetLepQual(kSilver);
      } else {
        lep.SetLepQual(kGold);
      }
    } // else
 
    list.push_back(lep);
  }
  return list;

}

template <>
ParticleList AnalysisBase<NANOULBase>::GetMuons(){
  ParticleList list;

  int N = nMuon;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Muon_pt[i] < 3. || fabs(Muon_eta[i]) > 2.4)
      continue;
    if(fabs(Muon_dxy[i]) >= 0.05 || fabs(Muon_dz[i]) >= 0.1 || Muon_sip3d[i] >= 6.)
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
    lep.SetTightCharge(Muon_tightCharge[i]);

    lep.SetRelIso(Muon_pfRelIso03_all[i]);
    lep.SetMiniIso(Muon_miniPFRelIso_all[i]);

    // FO baseline criteria
    lep.SetParticleID(kLoose);

    // signal lep criteria
    if(Muon_tightId[i])
      lep.SetParticleID(kTight);
    else if(Muon_mediumId[i])
      lep.SetParticleID(kMedium);
    if(lep.ParticleID() < kTight || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
      lep.SetLepQual(kBronze);
    else if(lep.SIP3D() >= 3.)
      lep.SetLepQual(kSilver);
    else
      lep.SetLepQual(kGold);
    list.push_back(lep);
  }
  return list;
}

/////////////////////////////////////////////////
// End NANOULBase specific methods
/////////////////////////////////////////////////

/////////////////////////////////////////////////
// Start NANORun3 specific methods
/////////////////////////////////////////////////

template <>
bool AnalysisBase<NANORun3>::PassEventFilter(){
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC
   return Flag_goodVertices && 
   Flag_globalSuperTightHalo2016Filter && 
   Flag_EcalDeadCellTriggerPrimitiveFilter && 
   Flag_BadPFMuonFilter && 
   Flag_BadPFMuonDzFilter && 
   Flag_hfNoisyHitsFilter && 
   Flag_eeBadScFilter && 
   Flag_ecalBadCalibFilter;
}

template <>
bool AnalysisBase<NANORun3>::GetMETtrigger(){
// https://cmshltinfo.app.cern.ch/summary?search=HLT_PFMET&year=2023&paths=true&prescaled=true&stream-types=Physics
  return (HLT_PFMET120_PFMHT120_IDTight ||
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
}

template <>
bool AnalysisBase<NANORun3>::GetMETORtrigger(){
  return (//HLT_PFMET110_PFMHT110_IDTight ||
    HLT_PFMET120_PFMHT120_IDTight ||
    HLT_PFMET130_PFMHT130_IDTight ||
    HLT_PFMET140_PFMHT140_IDTight ||
    //HLT_PFMETNoMu110_PFMHTNoMu110_IDTight ||
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight ||
    HLT_PFMETNoMu130_PFMHTNoMu130_IDTight ||
    HLT_PFMETNoMu140_PFMHTNoMu140_IDTight ||
    HLT_PFMET120_PFMHT120_IDTight_PFHT60 ||
    HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
}

template <>
bool AnalysisBase<NANORun3>::GetMETDoubleMutrigger(){
  return HLT_DoubleMu3_DZ_PFMET50_PFMHT60;
}

template <>
bool AnalysisBase<NANORun3>::GetSingleElectrontrigger(){
  return HLT_Ele30_WPTight_Gsf;
}

template <>
bool AnalysisBase<NANORun3>::GetSingleMuontrigger(){
  return HLT_IsoMu24;
}

template <>
bool AnalysisBase<NANORun3>::GetDoubleElectrontrigger(){
  return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
}

template <>
bool AnalysisBase<NANORun3>::GetDoubleMuontrigger(){
  return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
}

template <>
bool AnalysisBase<NANORun3>::GetTripleElectrontrigger(){
  return HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<NANORun3>::GetTripleMuonLowPTtrigger(){
  return (HLT_TripleMu_5_3_3_Mass3p8_DZ ||
          HLT_TripleMu_5_3_3_Mass3p8_DCA
         );
}

template <>
bool AnalysisBase<NANORun3>::GetTripleMuonHighPTtrigger(){
  return (HLT_TripleMu_12_10_5 ||
          HLT_TripleMu_10_5_5_DZ
         );
}

template <>
bool AnalysisBase<NANORun3>::GetDiMuEleLowPTtrigger(){
  return HLT_DiMu4_Ele9_CaloIdL_TrackIdL_DZ_Mass3p8;
}

template <>
bool AnalysisBase<NANORun3>::GetDiMuEleHighPTtrigger(){
  return HLT_DiMu9_Ele9_CaloIdL_TrackIdL || HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ;
}

template <>
bool AnalysisBase<NANORun3>::GetDiEleMutrigger(){
  return HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ || HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
}

template <>
bool AnalysisBase<NANORun3>::GetEMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<NANORun3>::GetEMuMutrigger(){
  return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>
bool AnalysisBase<NANORun3>::GetEMuEtrigger(){
  return HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
}

template <>  
double AnalysisBase<NANORun3>::GetBtagSFWeight(const ParticleList& jets, bool HForLF, int updown, ParticleIDType tag) {
  if (IsData()) 
      return 1.;

  bool FastSim = IsFastSim();
  int Njet = jets.size();
  int iflavor = 0;
  double probMC = 1.;
  double probDATA = 1.;
  std::string syst = "central";
  if(updown > 0) syst = "up";
  else if(updown < 0) syst = "down";
  
  for (int i = 0; i < Njet; i++) {
      if(abs(jets[i].PDGID()) == 5)
        iflavor = 5;
      else if(abs(jets[i].PDGID()) == 4)
        iflavor = 4;
      if(HForLF && iflavor == 0)
        continue;
      if(!HForLF && iflavor != 0)
        continue;
      std::vector<std::variant<int, double, std::string>> evalArgs;
      evalArgs.push_back(syst);
      evalArgs.push_back("M"); // Working Point ('M' for medium)
      evalArgs.push_back(iflavor);
      evalArgs.push_back(abs(jets[i].Eta()));
      evalArgs.push_back(jets[i].Pt());

      double SF = 1.;
      double EFF = 1.; // need to measure the efficiencies
      if(iflavor == 0)
        SF = m_cset_Btag->at("deepJet_light")->evaluate(evalArgs);
      else
        SF = m_cset_Btag->at("deepJet_comb")->evaluate(evalArgs);

      if (jets[i].BtagID() >= tag) {
          probMC *= EFF;
          probDATA *= SF * EFF;
      } else {
          probMC *= (1. - EFF);
          probDATA *= (1. - SF * EFF);
      }
  }

  if (probMC <= 0. || probDATA <= 0.)
      return 1.;

  return probDATA / probMC;
}

template <>
double AnalysisBase<NANORun3>::GetElIDSFWeight(const ParticleList& els, int updown){
  if(IsData())
    return 1.;

  bool FastSim = IsFastSim();

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
  
  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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

  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
  
  int Nlep = els.size();
  int pdg = 11;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(els[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(els[i].Pt(), els[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(els[i].Pt(), els[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(els[i].Pt(), els[i].Eta(), pdg, m_year);
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
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataIsoEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataIsoError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getIsoFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataSipEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0 && EFFData > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataSipError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getSipFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
  
  int Nlep = mus.size();
  int pdg = 13;
  double EFFMC, EFFData, SF;
  double EFFMCErr, EFFDataErr, SFErr;

  double probMC   = 1.;
  double probDATA = 1.;
  
  for(int i = 0; i < Nlep; i++){
    if(mus[i].SourceID() > 0)
      continue;
    
    EFFMC = m_LepSFTool.getMCVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
    EFFData = m_LepSFTool.getDataVLIdEfficiency(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

    if(EFFMC > 0)
      SF = EFFData/EFFMC;
    else
      SF = 1;

    SFErr = 0;
    if(updown != 0){
      EFFMCErr = m_LepSFTool.getMCVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
      EFFDataErr = m_LepSFTool.getDataVLIdError(mus[i].Pt(), mus[i].Eta(), pdg, m_year);

      if(EFFMC > 0 && EFFData > 0)
	SFErr = updown * SF * sqrt( EFFMCErr*EFFMCErr/EFFMC/EFFMC + EFFDataErr*EFFDataErr/EFFData/EFFData );
    }

    SF += SFErr;
    
    if(FastSim){
      double FSSF = m_LepSFTool.getVLIdFastSimSF(mus[i].Pt(), mus[i].Eta(), pdg, m_year);
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
    if(m_year == 2016)
      if(Jet_btagDeepFlavB[i] > 0.3093)
	btag = true;
    if(m_year == 2017)
      if(Jet_btagDeepFlavB[i] > 0.3033)
	btag = true;
    if(m_year == 2018)
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
void AnalysisBase<NANORun3>::ApplyMETPhiCorrections(TLorentzVector &metOut, int runNumber, int npvGood, const std::string &variation){
  if (!m_cset_METPhi) {
    std::cerr << "[ApplyMETPhiCorrections] ERROR: m_cset_METPhi not loaded.\n";
    return;
  }

  // read the inputs from the outgoing TLorentzVector (should be Type-1 propagated)
  double in_pt  = metOut.Pt();
  double in_phi = metOut.Phi();

  // helper to check key presence using at() (CorrectionSet throws out_of_range if missing)
  auto has_key = [&](const std::string& name) -> bool {
    try {
      m_cset_METPhi->at(name);
      return true;
    } catch (const std::out_of_range&) {
      return false;
    }
  };

  try {
    // Build met_type and epoch strings using existing flags (mirror your find_clib_file logic)
    std::string met_type = "PuppiMET";
    std::string dtmc = m_IsData ? "DATA" : "MC";
    std::string epoch;
    if (m_year == 2022) epoch = m_IsEE ? "2022EE" : "2022";
    else if (m_year == 2023) epoch = m_IsBPix ? "2023BPix" : "2023";
    else epoch = std::to_string(m_year); // fallback; most likely only 2022/2023 are used

    std::string var = variation.empty() ? "nom" : variation;

    // Build args: pt_phi (string), met_type (string), epoch (string), dtmc (string),
    //             variation (string), met_pt (real), met_phi (real), npvGood (real)
    std::vector<std::variant<int,double,std::string>> args;
    args.push_back(std::string("phi"));     // request phi
    args.push_back(met_type);
    args.push_back(epoch);
    args.push_back(dtmc);
    args.push_back(var);
    args.push_back(in_pt);
    args.push_back(in_phi);
    args.push_back(static_cast<double>(npvGood));

    std::vector<std::variant<int,double,std::string>> args_pt = args;
    args_pt[0] = std::string("pt"); // request pt

    double corrected_phi = in_phi;
    double corrected_pt  = in_pt;

    try {
      const std::string nodeName = "met_xy_corrections";
      if (!has_key(nodeName)) {
        std::cerr << "[ApplyMETPhiCorrections] WARNING: Run-3 node '" << nodeName
                  << "' not found in MET cset. Skipping Run-3 corrections.\n";
      } else {
        corrected_phi = m_cset_METPhi->at(nodeName)->evaluate(args);
      }
    } catch (const std::exception &e) {
      std::cerr << "[ApplyMETPhiCorrections] ERROR evaluating Run-3 phi node: "
                << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ApplyMETPhiCorrections] Unknown error evaluating Run-3 phi node.\n";
    }

    try {
      const std::string nodeName = "met_xy_corrections";
      if (!has_key(nodeName)) {
        // already warned above
      } else {
        corrected_pt = m_cset_METPhi->at(nodeName)->evaluate(args_pt);
      }
    } catch (const std::exception &e) {
      std::cerr << "[ApplyMETPhiCorrections] ERROR evaluating Run-3 pt node: "
                << e.what() << "\n";
    } catch (...) {
      std::cerr << "[ApplyMETPhiCorrections] Unknown error evaluating Run-3 pt node.\n";
    }

    wrapPhi(corrected_phi);
    metOut.SetPtEtaPhiM(corrected_pt, 0.0, corrected_phi, 0.);
    return;

  } catch (const std::exception &e) {
    std::cerr << "[ApplyMETPhiCorrections] Caught exception: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "[ApplyMETPhiCorrections] Caught unknown exception\n";
  }
}

template <>
ParticleList AnalysisBase<NANORun3>::GetJetsMET(TVector3& MET, int id){
  ParticleList list;
  int Njet = nJet;

  // Determine whether current Systematic is a JES/JER
  bool DO_JES = (m_SysTool.JESSystematics() == CurrentSystematic());
  bool DO_JER = (m_SysTool.JERSystematics() == CurrentSystematic());

  // Prepare vectors for MET propagation (RAW AK4)
  const int nAk4Jets = std::max(int(nJet), 0);
  std::vector<JecApplication::JetForMet> jetsForMet;
  jetsForMet.reserve(static_cast<size_t>(nAk4Jets));
  std::vector<JecApplication::JerInputs> jersForMet;
  jersForMet.reserve(static_cast<size_t>(nAk4Jets));
  long long EventNum = static_cast<long long>(GetEventNum());
  int RunNum = static_cast<int>(GetRunNum());
  int LumiBlock = static_cast<int>(GetLumiNum());

  for (int j = 0; j < nAk4Jets; ++j) {
    JecApplication::JetForMet v;
    v.phi  = Jet_phi[j];
    v.eta  = Jet_eta[j];
    v.area = Jet_area[j];
    v.rawPt = (1.0 - Jet_rawFactor[j]) * Jet_pt[j]; // reconstruct RAW from NanoAOD
    v.muonSubtrFactor = Jet_muonSubtrFactor[j];
    v.chEmEf = Jet_chEmEF[j];
    v.neEmEf = Jet_neEmEF[j];
    jetsForMet.push_back(v);

    // JER inputs for MET propagation
    JecApplication::JerInputs jerIn;
    jerIn.event = EventNum;
    jerIn.run = RunNum;
    jerIn.lumi = LumiBlock;
    const int genIdx = Jet_genJetIdx[j];
    if (genIdx > -1 && static_cast<UInt_t>(genIdx) < static_cast<UInt_t>(nGenJet)) {
      jerIn.hasGen = true;
      jerIn.genPt  = GenJet_pt[genIdx];
      jerIn.genEta = GenJet_eta[genIdx];
      jerIn.genPhi = GenJet_phi[genIdx];
      jerIn.maxDr  = 0.2; // AK4 match radius
    } else {
      jerIn.hasGen = false;
    }
    jersForMet.push_back(jerIn);
  }

  // Build SystematicOptions for MET helper (JES & JER)
  JecApplication::SystematicOptions systOpts{};
  if (!m_IsData) {
    // JES: resolve CLIB key (single resolution for MET propagation)
    if (DO_JES) {
      const std::string systLabel = CurrentSystematic().Label();
      static std::map<std::string, std::optional<std::string>> jesKeyCache;
      const std::string cacheKey = m_jmeYearKey + "|" + systLabel;
      std::optional<std::string> optClibKey;
      auto it = jesKeyCache.find(cacheKey);
      if (it != jesKeyCache.end()) {
        optClibKey = it->second;
      } else {
        optClibKey = findJESKeyForLabel(systLabel);
        jesKeyCache.emplace(cacheKey, optClibKey);
        if (!optClibKey)
          std::cerr << "[JME] Warning: could not resolve JES systematic '" << systLabel
                    << "' for year key " << m_jmeYearKey << ". MET will use nominal JES.\n";
      }

      if (optClibKey) {
        systOpts.jesSystName = *optClibKey;
        systOpts.jesSystVar  = (CurrentSystematic().IsUp() ? "Up" : "Down");
      } else {
        systOpts.jesSystName.clear();
        systOpts.jesSystVar.clear();
      }
    }

    if (DO_JER) systOpts.jerVar = (CurrentSystematic().IsUp() ? "up" : "down");
    else        systOpts.jerVar = "nom";
  }

  // --- Process jets for analysis (apply nominal JES, JES syst, JER) ---
  for(int i = 0; i < Njet; ++i){
    Particle jet;
    float mass = Jet_mass[i];
    if(std::isnan(mass) || std::isinf(mass) || mass < 0.) mass = 0.;

    // --- Nominal JES ---
    double ptRaw = Jet_pt[i] * (1.0 - Jet_rawFactor[i]);
    const double jesNominal = m_JMETool.GetJESFactorCLIB(
      m_jmeYearKey, m_IsData, RunNum, ptRaw,
      Jet_eta[i], Jet_phi[i], Jet_area[i], Rho_fixedGridRhoFastjetAll
    );
    jet.SetPtEtaPhiM(ptRaw * jesNominal, Jet_eta[i], Jet_phi[i], mass);

    // --- JES syst (per-jet multiplicative) ---
    if (DO_JES && !m_IsData && !systOpts.jesSystName.empty()) {
      const std::string var = CurrentSystematic().IsUp() ? "Up" : "Down";
      double jesSystFactor = m_JMETool.GetJESFactorSystCLIB(m_jmeYearKey, jet.Pt(), Jet_eta[i], Jet_phi[i], Jet_area[i], Rho_fixedGridRhoFastjetAll, systOpts.jesSystName, var);
      jet.SetPtEtaPhiM(jet.Pt() * jesSystFactor, jet.Eta(), jet.Phi(), mass);
    }

    // --- Jet selection ---
    if(Jet_pt[i] < 15. || fabs(Jet_eta[i]) > 5.) continue;
    if(Jet_jetId[i] < id) continue;

    if(Jet_jetId[i] >= 3) jet.SetParticleID(kTight);
    else if(Jet_jetId[i] >= 2) jet.SetParticleID(kMedium);
    else if(Jet_jetId[i] >= 1) jet.SetParticleID(kLoose);

    // --- B-tagging ---
    jet.SetBtag(Jet_btagDeepFlavB[i]);
    if(jet.Btag() > m_BtagVeryVeryTightWP) jet.SetBtagID(kVeryVeryTight);
    else if(jet.Btag() > m_BtagVeryTightWP) jet.SetBtagID(kVeryTight);
    else if(jet.Btag() > m_BtagTightWP) jet.SetBtagID(kTight);
    else if(jet.Btag() > m_BtagMediumWP) jet.SetBtagID(kMedium);
    else if(jet.Btag() > m_BtagLooseWP) jet.SetBtagID(kLoose);

    jet.SetPDGID(Jet_partonFlavour[i]);

    // --- JER: apply per-jet (nominal or up/down) via JMETool wrapper that uses JecApplication ---
    if(!m_IsData) {
      JecApplication::JesInputs jAfter{jet.Pt(), Jet_eta[i], Jet_phi[i], Jet_area[i], Rho_fixedGridRhoFastjetAll, 0.0};

      JecApplication::JerInputs jerIn = jersForMet[i];
      std::string jerVar = (DO_JER ? (CurrentSystematic().IsUp() ? "up" : "down") : "nom");
      const double sJer = m_JMETool.GetJERScaleFactor(m_jmeYearKey, jAfter, jerIn, jerVar, m_IsData);

      jet.SetPtEtaPhiM(jet.Pt() * sJer, jet.Eta(), jet.Phi(), jet.M() * sJer);
      jet.SetjetID(Jet_jetId[i]);
      jet.SetChEmEF(Jet_chEmEF[i]);
      jet.SetNeEmEF(Jet_neEmEF[i]);
    }

    list.push_back(jet);
  } // end jet loop

  // Call helper to get corrected MET 4-vector (JES+JER applied, MET propagation)
  TLorentzVector p4CorrectedMET = m_JMETool.GetCorrectedMET(m_jmeYearKey, m_IsData, RunNum, RawPuppiMET_pt, RawPuppiMET_phi, jetsForMet, jersForMet, Rho_fixedGridRhoFastjetAll, systOpts, false);

  int NPV = static_cast<int>(GetNPV());
  ApplyMETPhiCorrections(p4CorrectedMET, RunNum, NPV, "nom");

  MET.SetPtEtaPhi(p4CorrectedMET.Pt(), 0.0, p4CorrectedMET.Phi());

  // Unclustered MET systematic
  if(CurrentSystematic() == Systematic("METUncer_UnClust")){
    const double delta = (CurrentSystematic().IsUp() ? 1. : -1.);
    TVector3 dUncl(0., 0., 0.);
    if(delta >= 0.)
      dUncl.SetPtEtaPhi(PuppiMET_ptUnclusteredUp, 0., PuppiMET_phiUnclusteredUp);
    else
      dUncl.SetPtEtaPhi(PuppiMET_ptUnclusteredDown, 0., PuppiMET_phiUnclusteredDown);
    MET += dUncl;
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
ParticleList AnalysisBase<NANORun3>::GetJets(int id){
  TVector3 dum;
  return GetJetsMET(dum, id);
}

template <>
ParticleList AnalysisBase<NANORun3>::GetElectrons(){
  ParticleList list;

  int N = nElectron;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Electron_pt[i] < 7. || fabs(Electron_eta[i]) > 2.5)
      continue;
    if(Electron_pt[i] < 10. && fabs(Electron_eta[i]) > 1.4442) // remove low pt endcap ele
      continue;
    if(fabs(Electron_eta[i]) >= 1.4442 && fabs(Electron_eta[i]) <= 1.566)
      continue;
    if(fabs(Electron_dxy[i]) >= 0.05 || fabs(Electron_dz[i]) >= 0.1 ||
       Electron_sip3d[i] >= 6.)
      continue;
    if(Electron_pfRelIso03_all[i]*Electron_pt[i] >= 20. + 300./Electron_pt[i])
      continue;
    if(!minus_iso_hoe(Electron_vidNestedWPBitmap[i], 2, std::greater_equal<int>()))
      continue;
    if(Electron_lostHits[i] != 0)
      continue;
    if(Electron_convVeto[i] == 0)
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
    lep.SetTightCharge(Electron_tightCharge[i]);
    lep.SetParticleID(kVeryLoose); // need to set to something

    lep.SetRelIso(Electron_pfRelIso03_all[i]);
    lep.SetMiniIso(Electron_miniPFRelIso_all[i]);
    if( ( Electron_pt[i] < 20. && (
          lep.MiniIso()*lep.Pt() >= 4. || 
          lep.RelIso()*lep.Pt() >= 4. ||
          !minus_iso_hoe(Electron_vidNestedWPBitmap[i], 4, std::greater_equal<int>())
        )) ||
        ( Electron_pt[i] >= 20. && (
          !Electron_mvaIso_WP90[i]
        ))
      ) 
        lep.SetLepQual(kBronze);
    else if(( Electron_pt[i] < 20. && (
          lep.MiniIso()*lep.Pt() < 4. && 
          lep.RelIso()*lep.Pt() < 4. &&
          minus_iso_hoe(Electron_vidNestedWPBitmap[i], 4, std::greater_equal<int>())
        )) ||
        ( Electron_pt[i] >= 20. && (
          Electron_mvaIso_WP90[i]
        ))
      )
    { 
      if(lep.SIP3D() >= 3.)
        lep.SetLepQual(kSilver);
      else
        lep.SetLepQual(kGold);
    }
    list.push_back(lep);

  } // for(int i = 0; i < N; i++)
  return list;
}

template <>
ParticleList AnalysisBase<NANORun3>::GetLowPtElectrons(){
  ParticleList list;
  int N1 = nLowPtElectron;
  for(int i = 0; i < N1; i++){
    // baseline lepton definition
    if(LowPtElectron_pt[i] < 2. || LowPtElectron_pt[i] >= 7. || fabs(LowPtElectron_eta[i]) > 2.5)
      continue;
    if(fabs(LowPtElectron_eta[i]) >= 1.4442 && fabs(LowPtElectron_eta[i]) <= 1.566)
      continue;
    if(LowPtElectron_pt[i] < 10. && fabs(LowPtElectron_eta[i]) > 1.4442) // remove low pt endcap ele
      continue;
    if(LowPtElectron_convVeto[i] == 0)
      continue;
    if(fabs(LowPtElectron_dxy[i]) >= 0.05 || fabs(LowPtElectron_dz[i]) >= 0.1)
      continue;
    if(LowPtElectron_ID[i] < 1.5)
      continue;
    if(LowPtElectron_lostHits[i] != 0)
      continue;
    
    if(LowPtElectron_dxyErr[i] < 1.e-8 || LowPtElectron_dzErr[i] < 1.e-8)
      continue;
    
    // Calculate IP_3D and SIP_3D = IP_3D / IP_3D_err for the LowPtElectron collection.
    // - IP_3D:     3D impact parameter wrt first PV, in cm
    // - SIP_3D:    3D impact parameter significance wrt first PV, in cm
    float dxy       = LowPtElectron_dxy[i];
    float dz        = LowPtElectron_dz[i];
    float dxy_err   = LowPtElectron_dxyErr[i];
    float dz_err    = LowPtElectron_dzErr[i];
    // the SIP_3D defined below is coming from Suyash's talk, its something we defined.
    // The sigmas are needed for that calculation.
    float sigma_xy  = dxy/dxy_err;
    float sigma_z   = dz/dz_err;
    float IP_3D     = sqrt(dxy*dxy + dz*dz);
    float SIP_3D    = sqrt(sigma_xy*sigma_xy + sigma_z*sigma_z);

    if (SIP_3D >= 6.)
      continue;

    if(LowPtElectron_miniPFRelIso_all[i]*LowPtElectron_pt[i] >= 20. + 300./LowPtElectron_pt[i])
     continue;

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
    lep.SetTightCharge(0);

    lep.SetRelIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetMiniIso(LowPtElectron_miniPFRelIso_all[i]);
    lep.SetParticleID(kVeryLoose); // need to set to something for later on
    
    float pt = lep.Pt();
    float absEta = std::abs(lep.Eta());
    float id = LowPtElectron_ID[i];

    bool failsIso = lep.MiniIso() * pt >= 4. || lep.RelIso() * pt >= 4.;
    if (failsIso) {
      lep.SetLepQual(kBronze);
    } else {
      bool passesID = false;

      if (pt >= 2. && pt < 5.) {
        if (absEta >= 0.8 && absEta < 1.442)
          passesID = id >= 3;
        else if (absEta < 0.8)
          passesID = id >= 2.3;
      } else if (pt >= 5. && pt < 7.) {
        if (absEta >= 0.8 && absEta < 1.442)
          passesID = id >= 3;
        else if (absEta < 0.8)
          passesID = id >= 2.3;
      }

      if (!passesID) {
        lep.SetLepQual(kBronze);
      } else if (lep.SIP3D() >= 3.) {
        lep.SetLepQual(kSilver);
      } else {
        lep.SetLepQual(kGold);
      }
    } // else

    list.push_back(lep);
  }
  return list;

}

template <>
ParticleList AnalysisBase<NANORun3>::GetMuons(){
  ParticleList list;

  int N = nMuon;
  for(int i = 0; i < N; i++){
    // baseline lepton definition
    if(Muon_pt[i] < 3. || fabs(Muon_eta[i]) > 2.4)
      continue;
    if(fabs(Muon_dxy[i]) >= 0.05 || fabs(Muon_dz[i]) >= 0.1 || Muon_sip3d[i] >= 6.)
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
    lep.SetTightCharge(Muon_tightCharge[i]);

    lep.SetRelIso(Muon_pfRelIso03_all[i]);
    lep.SetMiniIso(Muon_miniPFRelIso_all[i]);

    // FO baseline criteria
    lep.SetParticleID(kLoose);

    // signal lep criteria
    if(Muon_tightId[i])
      lep.SetParticleID(kTight);
    else if(Muon_mediumId[i])
      lep.SetParticleID(kMedium);
    if(lep.ParticleID() < kTight || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
      lep.SetLepQual(kBronze);
    else if(lep.SIP3D() >= 3.)
      lep.SetLepQual(kSilver);
    else
      lep.SetLepQual(kGold);
    list.push_back(lep);
  }
  return list;
}

/////////////////////////////////////////////////
// End NANORun3 specific methods
/////////////////////////////////////////////////

template class AnalysisBase<SUSYNANOBase>; // preUL Run2 NANOAODv7
template class AnalysisBase<NANOULBase>; // UL Run2 NANOAODv9
template class AnalysisBase<NANORun3>; // Run3 NANOAODv12

