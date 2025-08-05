#ifndef ReducedNtuple_h
#define ReducedNtuple_h

#include "NtupleBase.hh"
#include "RestFrames/RestFrames.hh"

template class std::vector<std::vector<int> >;

using namespace RestFrames;

template <class Base>
class ReducedNtuple : public NtupleBase<Base> {

public:
  ReducedNtuple(TTree* tree = 0);
  virtual ~ReducedNtuple();

private:
  // Number of RJR analysis trees
  const static int m_aTrees=4; // [JetsLeps_ISR, Leps_ISR, Jets_ISR, Jets]

  vector<bool> m_treeSkipped;
  bool m_library_generated;
  TTree* InitOutputTree(const string& sample, bool do_slim = false);
  void FillOutputTree(TTree* tree, const Systematic& sys = Systematic::Default(), bool do_slim = false);

  void ClearVariables();

  bool m_event_skipped;
  
  // common variables for output tree
  double m_weight;
  double m_PUweight;
  double m_PUweight_up;
  double m_PUweight_down;
  double m_MuFweight;
  double m_MuFweight_up;
  double m_MuFweight_down;
  double m_MuRweight;
  double m_MuRweight_up;
  double m_MuRweight_down;
  double m_PDFweight;
  double m_PDFweight_up;
  double m_PDFweight_down;
  double m_BtagHFSFweight;
  double m_BtagHFSFweight_up;
  double m_BtagHFSFweight_down;
  double m_BtagLFSFweight;
  double m_BtagLFSFweight_up;
  double m_BtagLFSFweight_down;

  double m_elIDSFweight;
  double m_elIDSFweight_up;
  double m_elIDSFweight_down;
  double m_elISOSFweight;
  double m_elISOSFweight_up;
  double m_elISOSFweight_down;
  double m_elSIPSFweight;
  double m_elSIPSFweight_up;
  double m_elSIPSFweight_down;
  double m_elVLSFweight;
  double m_elVLSFweight_up;
  double m_elVLSFweight_down;
  double m_muIDSFweight;
  double m_muIDSFweight_up;
  double m_muIDSFweight_down;
  double m_muISOSFweight;
  double m_muISOSFweight_up;
  double m_muISOSFweight_down;
  double m_muSIPSFweight;
  double m_muSIPSFweight_up;
  double m_muSIPSFweight_down;
  double m_muVLSFweight;
  double m_muVLSFweight_up;
  double m_muVLSFweight_down;
  
  double m_MetTrigSFweight;
  double m_MetTrigSFweight_up;
  double m_MetTrigSFweight_down;
  int    m_MetTrigSFCurveIndex;

  int m_runnum;
  int m_luminum;
  Long64_t m_eventnum;

  // additional vars for signal
  int m_cascades_tree;
  int m_cascades_prod;
  int m_cascades_SlepSneu_1stDecay;
  int m_cascades_SlepSneu_2ndDecay;
  int m_cascades_N2_1stDecay;
  int m_cascades_N2_2ndDecay;
  int m_Npartons;
  int m_MSlepL;
  int m_MSneu;
  int m_MN2;
  int m_MC1;
  int m_MN1;
  int m_MP;
  int m_NSparticleW;
  int m_NSparticleZ;
  int m_NSparticlePhoton;
  vector<int> m_LSPParents;

  // Additional vars for gen weighting
  double m_XSec;
  double m_genweight;
  double m_Nweight;
  double m_Nevent;

  int m_NPV;
  int m_NPU;
  
  double m_MET;
  double m_MET_phi;

  double m_altMET;
  double m_altMET_phi;

  double m_genMET;
  double m_genMET_phi;

  double m_LHE_HT;
  double m_LHE_HTIncoming;

  double m_HT_eta24;
  double m_HT_eta24_id;
  double m_HT_eta3;
  double m_HT_eta3_id;
  double m_HT_eta5;
  double m_HT_eta5_id;
  
  bool m_EventFilter;
  bool m_FastSimEventVeto;

  bool m_EventFlag_FailJetID;
  bool m_EventFlag_JetInHEM;
  bool m_EventFlag_JetInHEM_Pt20;
  bool m_EventFlag_JetInHEM_Pt20_JetID;
  bool m_HEM_Veto;

  double m_PrefireWeight;
  double m_PrefireWeight_up;
  double m_PrefireWeight_down;
  
  bool m_METtrigger;
  bool m_METORtrigger;

  bool m_SingleElectrontrigger;
  bool m_SingleMuontrigger;
  bool m_DoubleElectrontrigger;
  bool m_DoubleMuontrigger;
  bool m_TripleElectrontrigger;
  bool m_TripleMuontrigger;
  bool m_EMutrigger;  
  
  int m_Nele;
  int m_Nlele;
  int m_Nmu;
  
  int m_Nlep;
  vector<double> m_PT_lep;
  vector<double> m_Eta_lep;
  vector<double> m_Phi_lep;
  vector<double> m_M_lep;
  vector<int>    m_Charge_lep;
  vector<int>    m_PDGID_lep;
  vector<double> m_RelIso_lep;
  vector<double> m_MiniIso_lep;
  vector<double> m_Dxy_lep;
  vector<double> m_DxyErr_lep;
  vector<double> m_Dz_lep;
  vector<double> m_DzErr_lep;
  vector<double> m_IP3D_lep;
  vector<double> m_SIP3D_lep;
  vector<int>    m_ID_lep;
  vector<int>    m_SourceID_lep;
  vector<int>    m_LepQual_lep;
  vector<int>    m_IsLowPt_lep;
  vector<int>    m_TightCharge_lep;
  vector<int>    m_Index_lep;

  int m_Njet;
  int m_Nbjet;
  vector<double> m_PT_jet;
  vector<double> m_Eta_jet;
  vector<double> m_Phi_jet;
  vector<double> m_M_jet;
  vector<double> m_Btag_jet;
  vector<int>    m_BtagID_jet;
  vector<double> m_Flavor_jet;
  vector<double> m_ProbB_SV;
  vector<double> m_ProbC_SV;

  int m_NGenjet;
  vector<double> m_PT_Genjet;
  vector<double> m_Eta_Genjet;
  vector<double> m_Phi_Genjet;
  vector<double> m_M_Genjet;
  vector<int>    m_Index_jet;

  int m_NSV;
  vector<double> m_PT_SV;
  vector<double> m_Eta_SV;
  vector<double> m_Phi_SV;
  vector<double> m_M_SV;

  int m_genNele;
  int m_genNmu;

  int m_genNlep;
  vector<double> m_genPT_lep;
  vector<double> m_genEta_lep;
  vector<double> m_genPhi_lep;
  vector<double> m_genM_lep;
  vector<int>    m_genCharge_lep;
  vector<int>    m_genPDGID_lep;
  vector<int>    m_genMomPDGID_lep;
  vector<int>    m_genSourceID_lep;
  vector<int>    m_genIndex_lep;
  vector<int>    m_genMomIndex_lep;

  int m_genNnu;
  vector<double> m_genPT_nu;
  vector<double> m_genEta_nu;
  vector<double> m_genPhi_nu;
  vector<int>    m_genPDGID_nu;
  vector<int>    m_genMomPDGID_nu;
  
  int m_genNboson;
  vector<double> m_genPT_boson;
  vector<double> m_genEta_boson;
  vector<double> m_genPhi_boson;
  vector<double> m_genM_boson;
  vector<int>    m_genPDGID_boson;
  vector<int>    m_genMomPDGID_boson;
  
  int m_genNsusy;
  vector<double> m_genPT_susy;
  vector<double> m_genEta_susy;
  vector<double> m_genPhi_susy;
  vector<double> m_genM_susy;
  vector<int>    m_genPDGID_susy;
  vector<int>    m_genMomPDGID_susy;

  //////////////////////
  // derived observables
  //////////////////////

  // Object Counting Variables

  vector<double> m_dphi_lep_S;
  vector<double> m_cos_lep_S;
  vector<double> m_dphi_SV_S;
  vector<double> m_cos_SV_S;
  vector<double> m_dphi_jet_S;
  vector<double> m_cos_jet_S;
  
  vector<double> m_dphiMET_lep_S;
  vector<double> m_dphiMET_SV_S; 
  vector<double> m_dphiMET_jet_S; 
    
  int m_Njet_ISR;
  int m_Njet_S;
  int m_Nbjet_ISR;
  int m_Nbjet_S;
  int m_Nlep_ISR;
  int m_Nlep_S;
  int m_NSV_ISR;
  int m_NSV_S;
  vector<int> m_index_jet_ISR;
  vector<int> m_index_jet_S;
  vector<int> m_index_SV_ISR;
  vector<int> m_index_SV_S;
  vector<int> m_index_lep_ISR;
  vector<int> m_index_lep_S;

  int m_Njet_a;
  int m_Njet_b;
  int m_Nbjet_a;
  int m_Nbjet_b;
  int m_Nlep_a;
  int m_Nlep_b;
  int m_NSV_a;
  int m_NSV_b;
 
  vector<int> m_index_jet_a;
  vector<int> m_index_jet_b;
  vector<int> m_index_lep_a;
  vector<int> m_index_lep_b;
  vector<int> m_index_SV_a;
  vector<int> m_index_SV_b;

  // LEP ISR tree
  int m_Nlep_a_LEP;
  int m_Nlep_b_LEP;
 
  vector<int> m_index_lep_a_LEP;
  vector<int> m_index_lep_b_LEP;

  // JET ISR tree
  int m_Njet_ISR_JET_ISR;
  int m_Njet_S_JET_ISR;
  int m_Nbjet_ISR_JET_ISR;
  int m_Nbjet_S_JET_ISR;
  vector<int> m_index_jet_ISR_JET_ISR;
  vector<int> m_index_jet_S_JET_ISR;

  int m_Njet_a_JET_ISR;
  int m_Njet_b_JET_ISR;
  int m_Nbjet_a_JET_ISR;
  int m_Nbjet_b_JET_ISR;
 
  vector<int> m_index_jet_a_JET_ISR;
  vector<int> m_index_jet_b_JET_ISR;

  // JET tree
  int m_Njet_a_JET;
  int m_Njet_b_JET;
  int m_Nbjet_a_JET;
  int m_Nbjet_b_JET;
 
  vector<int> m_index_jet_a_JET;
  vector<int> m_index_jet_b_JET;
  
  // Kinematic Variables

  double m_PTCM;
  double m_PzCM;
  double m_cosCM;
  double m_dphiCM;
  double m_dphiCMI;
  
  double m_MS;
  double m_PS;
  double m_cosS;
  double m_dphiS;
  double m_dphiSI;
  double m_PTS;
  double m_PzS;

  double m_EtaCM;
  double m_PhiCM;
  double m_MCM;
  double m_EtaS;
  double m_PhiS;
  double m_LAB_Pt;
  double m_LAB_Eta;
  double m_LAB_Phi;
  double m_LAB_M; 

  double m_EVa;
  double m_EVb;
  double m_PVa;
  double m_PVb;
  double m_EJa;
  double m_EJb;
  double m_PJa;
  double m_PJb;

  double m_MX2a;
  double m_MX2b;
  double m_ELa;
  double m_ELb;
  double m_PLa;
  double m_PLb;

  double m_MV;
  double m_PV;
  double m_MVa;
  double m_MVb;

  double m_PV_lab;
  double m_dphiMET_V;

  double m_MJ;
  double m_ML;
  double m_EJ;
  double m_EL;
  double m_PJ;
  double m_PL;
  
  double m_PX2;
  double m_PX2_BoostT;
  double m_MX2a_BoostT;
  double m_MX2b_BoostT;
  double m_Mperp;
  double m_gammaT;

  double m_PV_BoostT;
  
  double m_EVa_BoostT;
  double m_EVb_BoostT;
  double m_PVa_BoostT;
  double m_PVb_BoostT;

  double m_EJ_BoostT;
  double m_EL_BoostT;
  double m_PJ_BoostT;
  double m_PL_BoostT;
  
  double m_MJa;
  double m_MJb;
  double m_MLa;
  double m_MLb;
  double m_cosJa;
  double m_cosJb;
  double m_cosLa;
  double m_cosLb;

  double m_MT2;

  // ISR related variables
  double m_PISR;
  double m_PTISR;
  double m_MISR;
  double m_RISR;
  double m_RISRT;

  // New Observables for Run3/Cascades
  // analysis using ISR boosted tree
  // variables using 4-vects constructed by evaluating in CM frame and making perp to boost
  double m_MSperpCM0; 
  double m_MaPerpCM0;
  double m_MbPerpCM0;
  double m_MaVPerpCM0;
  double m_MbVPerpCM0;
  double m_MQperpCM0; // sqrt(Ma*Ma + Mb*Mb)/sqrt(2)
  double m_gammaPerpCM0; // 2*MQ/MS
  double m_MVisAperpCM0;
  double m_MVisBperpCM0;
  // variables using 4-vects constructed by evaluating in CM frame
  double m_MSCM0; 
  double m_MaCM0;
  double m_MbCM0;
  double m_MaVCM0;
  double m_MbVCM0;
  double m_MQCM0; // sqrt(Ma*Ma + Mb*Mb)/sqrt(2)
  double m_gammaCM0; // 2*MQ/MS
  double m_MVisACM0;
  double m_MVisBCM0;

  // New info for RJR 'LEP' tree
  double m_RISR_LEP;
  double m_PTISR_LEP;
  double m_MS_LEP;
  double m_MSV_LEP;
  double m_MQ_LEP;
  double m_gamma_LEP;
  double m_MPa_LEP;
  double m_MPb_LEP;
  double m_MVa_LEP;
  double m_MVb_LEP;
  double m_PTS_CM_LEP;
  double m_MS_S0_LEP;
  double m_MV_S0_LEP;
  double m_MQ_S0_LEP;
  double m_gamma_S0_LEP;
  double m_MPTilde_LEP;
  double m_MSTilde_LEP;
  double m_gammaTilde_LEP;
  double m_CosDecayAngle_Pa_LEP;
  double m_CosDecayAngle_Pb_LEP;
  double m_CosDecayAngle_Va_LEP;
  double m_CosDecayAngle_Vb_LEP;
  double m_CosDecayAngle_S_LEP;
  double m_RZPara_LEP;
  double m_MINV_LEP;

  double m_PTCM_LEP;
  double m_PzCM_LEP;
  double m_cosCM_LEP;
  double m_dphiCM_LEP;
  double m_dphiCMI_LEP;
  double m_dphiMET_V_LEP;

  // New info for RJR 'JET ISR' tree
  double m_gammaT_JET_ISR;
  double m_Mperp_JET_ISR;
  double m_RISR_JET_ISR;
  double m_PTISR_JET_ISR;
  double m_MS_JET_ISR;
  double m_MSV_JET_ISR;
  double m_MQ_JET_ISR;
  double m_gamma_JET_ISR;
  double m_MPa_JET_ISR;
  double m_MPb_JET_ISR;
  double m_MVa_JET_ISR;
  double m_MVb_JET_ISR;
  double m_PTS_CM_JET_ISR;
  double m_MS_S0_JET_ISR;
  double m_MV_S0_JET_ISR;
  double m_MQ_S0_JET_ISR;
  double m_gamma_S0_JET_ISR;
  double m_MPTilde_JET_ISR;
  double m_MSTilde_JET_ISR;
  double m_gammaTilde_JET_ISR;
  double m_CosDecayAngle_Pa_JET_ISR;
  double m_CosDecayAngle_Pb_JET_ISR;
  double m_CosDecayAngle_Va_JET_ISR;
  double m_CosDecayAngle_Vb_JET_ISR;
  double m_CosDecayAngle_S_JET_ISR;
  double m_RZPara_JET_ISR;
  double m_MINV_JET_ISR;

  double m_PTCM_JET_ISR;
  double m_PzCM_JET_ISR;
  double m_cosCM_JET_ISR;
  double m_dphiCM_JET_ISR;
  double m_dphiCMI_JET_ISR;
  double m_dphiMET_V_JET_ISR;

  // New info for RJR 'JET' tree
  double m_MS_JET;
  double m_MSV_JET;
  double m_MQ_JET;
  double m_gamma_JET;
  double m_MPa_JET;
  double m_MPb_JET;
  double m_MVa_JET;
  double m_MVb_JET;
  double m_MS_S0_JET;
  double m_MV_S0_JET;
  double m_MQ_S0_JET;
  double m_gamma_S0_JET;
  double m_MPTilde_JET;
  double m_MSTilde_JET;
  double m_gammaTilde_JET;
  double m_CosDecayAngle_Pa_JET;
  double m_CosDecayAngle_Pb_JET;
  double m_CosDecayAngle_Va_JET;
  double m_CosDecayAngle_Vb_JET;
  double m_CosDecayAngle_S_JET;

  double m_PTCM_JET;
  double m_PzCM_JET;
  double m_cosCM_JET;
  double m_dphiCM_JET;
  double m_dphiCMI_JET;
  double m_dphiMET_V_JET;

  // RestFrames frames and friends
  LabRecoFrame*     LAB[m_aTrees];
  DecayRecoFrame*   CM[m_aTrees];
  DecayRecoFrame*   S[m_aTrees];
  DecayRecoFrame*   X2a[m_aTrees];
  DecayRecoFrame*   X2b[m_aTrees];
  VisibleRecoFrame*   ISR[m_aTrees];
  VisibleRecoFrame*   Ja[m_aTrees];
  VisibleRecoFrame*   Jb[m_aTrees];
  VisibleRecoFrame*   La[m_aTrees];
  VisibleRecoFrame*   Lb[m_aTrees];
  InvisibleRecoFrame* X1a[m_aTrees];
  InvisibleRecoFrame* X1b[m_aTrees];

  InvisibleGroup*       INV[m_aTrees];
  SetMassInvJigsaw*     InvM[m_aTrees];
  SetRapidityInvJigsaw* InvEta[m_aTrees];
  MinMassesSqInvJigsaw* InvSplit[m_aTrees];
  
  CombinatoricGroup*   COMB_J[m_aTrees];
  MinMassesSqCombJigsaw* CombSplit_ISR[m_aTrees];
  MinMassesSqCombJigsaw* CombSplit_J[m_aTrees];

  CombinatoricGroup*   COMB_L[m_aTrees];
  MinMassesSqCombJigsaw* CombSplit_L[m_aTrees];

};

#endif
