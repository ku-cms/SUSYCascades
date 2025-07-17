//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  5 17:12:21 2025 by ROOT version 6.26/11
// from TTree KUAnalysis/KUAnalysis
// found on file: /local-scratch/zflowers/NTUPLES/HADD/Summer23BPix_130X/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Summer23BPix_130X.root
//////////////////////////////////////////////////////////

#ifndef ReducedBase_V2_h
#define ReducedBase_V2_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using std::vector;

class ReducedBase_V2 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          event_skipped;
   vector<bool>    *treeSkipped;
   Double_t        weight;
   Double_t        genweight;
   Double_t        PUweight;
   Double_t        PUweight_up;
   Double_t        PUweight_down;
   Double_t        MuFweight;
   Double_t        MuFweight_up;
   Double_t        MuFweight_down;
   Double_t        MuRweight;
   Double_t        MuRweight_up;
   Double_t        MuRweight_down;
   Double_t        PDFweight;
   Double_t        PDFweight_up;
   Double_t        PDFweight_down;
   Double_t        BtagHFSFweight;
   Double_t        BtagHFSFweight_up;
   Double_t        BtagHFSFweight_down;
   Double_t        BtagLFSFweight;
   Double_t        BtagLFSFweight_up;
   Double_t        BtagLFSFweight_down;
   Double_t        elIDSFweight;
   Double_t        elIDSFweight_up;
   Double_t        elIDSFweight_down;
   Double_t        elISOSFweight;
   Double_t        elISOSFweight_up;
   Double_t        elISOSFweight_down;
   Double_t        elSIPSFweight;
   Double_t        elSIPSFweight_up;
   Double_t        elSIPSFweight_down;
   Double_t        elVLSFweight;
   Double_t        elVLSFweight_up;
   Double_t        elVLSFweight_down;
   Double_t        muIDSFweight;
   Double_t        muIDSFweight_up;
   Double_t        muIDSFweight_down;
   Double_t        muISOSFweight;
   Double_t        muISOSFweight_up;
   Double_t        muISOSFweight_down;
   Double_t        muSIPSFweight;
   Double_t        muSIPSFweight_up;
   Double_t        muSIPSFweight_down;
   Double_t        muVLSFweight;
   Double_t        muVLSFweight_up;
   Double_t        muVLSFweight_down;
   Double_t        MetTrigSFweight;
   Double_t        MetTrigSFweight_up;
   Double_t        MetTrigSFweight_down;
   Int_t           MetTrigSFCurveIndex;
   Int_t           runnum;
   Int_t           luminum;
   Long64_t        eventnum;
   Int_t           NPV;
   Bool_t          EventFilter;
   Bool_t          FastSimEventVeto;
   Double_t        PrefireWeight;
   Double_t        PrefireWeight_up;
   Double_t        PrefireWeight_down;
   Bool_t          METtrigger;
   Bool_t          METORtrigger;
   Bool_t          DoubleElectrontrigger;
   Bool_t          DoubleMuontrigger;
   Bool_t          EventFlag_FailJetID;
   Bool_t          EventFlag_JetInHEM;
   Bool_t          EventFlag_JetInHEM_Pt20;
   Bool_t          EventFlag_JetInHEM_Pt20_JetID;
   Bool_t          HEM_Veto;
   Bool_t          SingleElectrontrigger;
   Bool_t          SingleMuontrigger;
   Bool_t          EMutrigger;
   Double_t        MET;
   Double_t        MET_phi;
   Double_t        altMET;
   Double_t        altMET_phi;
   Double_t        LHE_HT;
   Double_t        LHE_HTIncoming;
   Double_t        HT_eta24;
   Double_t        HT_eta24_id;
   Double_t        HT_eta3;
   Double_t        HT_eta3_id;
   Double_t        HT_eta5;
   Double_t        HT_eta5_id;
   Int_t           Nele;
   Int_t           Nlele;
   Int_t           Nmu;
   Int_t           Nlep;
   vector<double>  *PT_lep;
   vector<double>  *Eta_lep;
   vector<double>  *Phi_lep;
   vector<double>  *M_lep;
   vector<int>     *Charge_lep;
   vector<int>     *PDGID_lep;
   vector<int>     *ID_lep;
   vector<int>     *SourceID_lep;
   vector<int>     *LepQual_lep;
   vector<int>     *IsLowPt_lep;
   vector<int>     *TightCharge_lep;
   vector<int>     *Index_lep;
   vector<double>  *RelIso_lep;
   vector<double>  *MiniIso_lep;
   vector<double>  *Dxy_lep;
   vector<double>  *DxyErr_lep;
   vector<double>  *Dz_lep;
   vector<double>  *DzErr_lep;
   vector<double>  *IP3D_lep;
   vector<double>  *SIP3D_lep;
   Int_t           Njet;
   Int_t           Nbjet;
   vector<double>  *PT_jet;
   vector<double>  *Eta_jet;
   vector<double>  *Phi_jet;
   vector<double>  *M_jet;
   vector<double>  *Btag_jet;
   vector<int>     *BtagID_jet;
   vector<double>  *Flavor_jet;
   vector<int>     *index_jet_a;
   vector<int>     *index_jet_b;
   vector<int>     *index_jet_ISR;
   vector<int>     *index_jet_S;
   vector<double>  *PT_Genjet;
   vector<double>  *Eta_Genjet;
   vector<double>  *Phi_Genjet;
   vector<double>  *M_Genjet;
   vector<int>     *Index_jet;
   Int_t           Njet_ISR;
   Int_t           Njet_S;
   Int_t           Nbjet_ISR;
   Int_t           Nbjet_S;
   Int_t           Nlep_ISR;
   Int_t           Nlep_S;
   vector<int>     *index_lep_ISR;
   vector<int>     *index_lep_S;
   vector<double>  *dphi_lep_S;
   vector<double>  *cos_lep_S;
   vector<double>  *dphi_jet_S;
   vector<double>  *cos_jet_S;
   vector<double>  *dphiMET_lep_S;
   vector<double>  *dphiMET_jet_S;
   Int_t           Njet_a;
   Int_t           Njet_b;
   Int_t           Nbjet_a;
   Int_t           Nbjet_b;
   Int_t           Nlep_a;
   Int_t           Nlep_b;
   vector<int>     *index_lep_a;
   vector<int>     *index_lep_b;
   Double_t        PTCM;
   Double_t        PzCM;
   Double_t        cosCM;
   Double_t        dphiCM;
   Double_t        dphiCMI;
   Double_t        dphiMET_V;
   Double_t        Mperp;
   Double_t        gammaT;
   Double_t        EJ_BoostT;
   Double_t        EL_BoostT;
   Double_t        PTISR;
   Double_t        RISR;
   Double_t        EtaCM;
   Double_t        PhiCM;
   Double_t        MCM;
   Double_t        EtaS;
   Double_t        PhiS;
   Double_t        LAB_Pt;
   Double_t        LAB_Eta;
   Double_t        LAB_Phi;
   Double_t        LAB_M;
   Double_t        MSperpCM0;
   Double_t        MaPerpCM0;
   Double_t        MbPerpCM0;
   Double_t        MaVPerpCM0;
   Double_t        MbVPerpCM0;
   Double_t        MQperpCM0;
   Double_t        gammaPerpCM0;
   Double_t        MVisAperpCM0;
   Double_t        MVisBperpCM0;
   Double_t        MSCM0;
   Double_t        MaCM0;
   Double_t        MbCM0;
   Double_t        MaVCM0;
   Double_t        MbVCM0;
   Double_t        MQCM0;
   Double_t        gammaCM0;
   Double_t        MVisACM0;
   Double_t        MVisBCM0;
   Double_t        MS;
   Double_t        PS;
   Double_t        cosS;
   Double_t        dphiS;
   Double_t        dphiSI;
   Double_t        PTS;
   Double_t        PzS;
   Double_t        EVa;
   Double_t        EVb;
   Double_t        PVa;
   Double_t        PVb;
   Double_t        EJa;
   Double_t        EJb;
   Double_t        PJa;
   Double_t        PJb;
   Double_t        MX2a;
   Double_t        MX2b;
   Double_t        ELa;
   Double_t        ELb;
   Double_t        PLa;
   Double_t        PLb;
   Double_t        MV;
   Double_t        PV;
   Double_t        MVa;
   Double_t        MVb;
   Double_t        PV_lab;
   Double_t        MJa;
   Double_t        MJb;
   Double_t        MLa;
   Double_t        MLb;
   Double_t        cosJa;
   Double_t        cosJb;
   Double_t        cosLa;
   Double_t        cosLb;
   Double_t        MJ;
   Double_t        ML;
   Double_t        EJ;
   Double_t        EL;
   Double_t        PJ;
   Double_t        PL;
   Double_t        PX2;
   Double_t        PX2_BoostT;
   Double_t        MX2a_BoostT;
   Double_t        MX2b_BoostT;
   Double_t        PV_BoostT;
   Double_t        EVa_BoostT;
   Double_t        EVb_BoostT;
   Double_t        PVa_BoostT;
   Double_t        PVb_BoostT;
   Double_t        PJ_BoostT;
   Double_t        PL_BoostT;
   Double_t        PISR;
   Double_t        RISRT;
   Double_t        MISR;
   Int_t           NPU;
   Double_t        genMET;
   Double_t        genMET_phi;
   Int_t           genNele;
   Int_t           genNmu;
   Int_t           genNlep;
   vector<double>  *genPT_lep;
   vector<double>  *genEta_lep;
   vector<double>  *genPhi_lep;
   vector<double>  *genM_lep;
   vector<int>     *genCharge_lep;
   vector<int>     *genPDGID_lep;
   vector<int>     *genMomPDGID_lep;
   vector<int>     *genSourceID_lep;
   vector<int>     *genIndex_lep;
   vector<int>     *genMomIndex_lep;
   Int_t           genNnu;
   vector<double>  *genPT_nu;
   vector<double>  *genEta_nu;
   vector<double>  *genPhi_nu;
   vector<int>     *genPDGID_nu;
   vector<int>     *genMomPDGID_nu;
   Int_t           genNboson;
   vector<double>  *genPT_boson;
   vector<double>  *genEta_boson;
   vector<double>  *genPhi_boson;
   vector<double>  *genM_boson;
   vector<int>     *genPDGID_boson;
   vector<int>     *genMomPDGID_boson;
   Int_t           genNsusy;
   vector<double>  *genPT_susy;
   vector<double>  *genEta_susy;
   vector<double>  *genPhi_susy;
   vector<double>  *genM_susy;
   vector<int>     *genPDGID_susy;
   vector<int>     *genMomPDGID_susy;
   Double_t        MT2;
   Double_t        RISR_LEP;
   Double_t        PTISR_LEP;
   Double_t        MS_LEP;
   Double_t        MSV_LEP;
   Double_t        MQ_LEP;
   Double_t        gamma_LEP;
   Double_t        MPa_LEP;
   Double_t        MPb_LEP;
   Double_t        MVa_LEP;
   Double_t        MVb_LEP;
   Double_t        PTS_CM_LEP;
   Double_t        MS_S0_LEP;
   Double_t        MV_S0_LEP;
   Double_t        MQ_S0_LEP;
   Double_t        gamma_S0_LEP;
   Double_t        MPTilde_LEP;
   Double_t        MSTilde_LEP;
   Double_t        gammaTilde_LEP;
   Double_t        CosDecayAngle_Pa_LEP;
   Double_t        CosDecayAngle_Pb_LEP;
   Double_t        CosDecayAngle_S_LEP;
   Double_t        RZPara_LEP;
   Double_t        MINV_LEP;
   Int_t           cascades_tree;
   Int_t           cascades_prod;
   Int_t           cascades_SlepSneu_1stDecay;
   Int_t           cascades_SlepSneu_2ndDecay;
   Int_t           cascades_N2_1stDecay;
   Int_t           cascades_N2_2ndDecay;
   Int_t           MSlepL;
   Int_t           MSneu;
   Int_t           MN2;
   Int_t           MC1;
   Int_t           MN1;
   Int_t           MP;
   Int_t           NSparticleW;
   Int_t           NSparticleZ;
   Int_t           NSparticlePhoton;
   vector<int>     *LSPParents;
   Int_t           Npartons;
   Double_t        XSec;
   Double_t        Nweight;
   Double_t        Nevent;
   Double_t        Mperp_JET;
   Double_t        gammaT_JET;
   Double_t        RISR_JET;
   Double_t        PTISR_JET;
   Double_t        MS_JET;
   Double_t        MSV_JET;
   Double_t        MQ_JET;
   Double_t        gamma_JET;
   Double_t        MPa_JET;
   Double_t        MPb_JET;
   Double_t        MVa_JET;
   Double_t        MVb_JET;
   Double_t        PTS_CM_JET;
   Double_t        MS_S0_JET;
   Double_t        MV_S0_JET;
   Double_t        MQ_S0_JET;
   Double_t        gamma_S0_JET;
   Double_t        MPTilde_JET;
   Double_t        MSTilde_JET;
   Double_t        gammaTilde_JET;
   Double_t        CosDecayAngle_Pa_JET;
   Double_t        CosDecayAngle_Pb_JET;
   Double_t        CosDecayAngle_S_JET;
   Double_t        RZPara_JET;
   Double_t        MINV_JET;

   // List of branches
   TBranch        *b_event_skipped;   //!
   TBranch        *b_treeSkipped;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_PUweight_up;   //!
   TBranch        *b_PUweight_down;   //!
   TBranch        *b_MuFweight;   //!
   TBranch        *b_MuFweight_up;   //!
   TBranch        *b_MuFweight_down;   //!
   TBranch        *b_MuRweight;   //!
   TBranch        *b_MuRweight_up;   //!
   TBranch        *b_MuRweight_down;   //!
   TBranch        *b_PDFweight;   //!
   TBranch        *b_PDFweight_up;   //!
   TBranch        *b_PDFweight_down;   //!
   TBranch        *b_BtagHFSFweight;   //!
   TBranch        *b_BtagHFSFweight_up;   //!
   TBranch        *b_BtagHFSFweight_down;   //!
   TBranch        *b_BtagLFSFweight;   //!
   TBranch        *b_BtagLFSFweight_up;   //!
   TBranch        *b_BtagLFSFweight_down;   //!
   TBranch        *b_elIDSFweight;   //!
   TBranch        *b_elIDSFweight_up;   //!
   TBranch        *b_elIDSFweight_down;   //!
   TBranch        *b_elISOSFweight;   //!
   TBranch        *b_elISOSFweight_up;   //!
   TBranch        *b_elISOSFweight_down;   //!
   TBranch        *b_elSIPSFweight;   //!
   TBranch        *b_elSIPSFweight_up;   //!
   TBranch        *b_elSIPSFweight_down;   //!
   TBranch        *b_elVLSFweight;   //!
   TBranch        *b_elVLSFweight_up;   //!
   TBranch        *b_elVLSFweight_down;   //!
   TBranch        *b_muIDSFweight;   //!
   TBranch        *b_muIDSFweight_up;   //!
   TBranch        *b_muIDSFweight_down;   //!
   TBranch        *b_muISOSFweight;   //!
   TBranch        *b_muISOSFweight_up;   //!
   TBranch        *b_muISOSFweight_down;   //!
   TBranch        *b_muSIPSFweight;   //!
   TBranch        *b_muSIPSFweight_up;   //!
   TBranch        *b_muSIPSFweight_down;   //!
   TBranch        *b_muVLSFweight;   //!
   TBranch        *b_muVLSFweight_up;   //!
   TBranch        *b_muVLSFweight_down;   //!
   TBranch        *b_MetTrigSFweight;   //!
   TBranch        *b_MetTrigSFweight_up;   //!
   TBranch        *b_MetTrigSFweight_down;   //!
   TBranch        *b_MetTrigSFCurveIndex;   //!
   TBranch        *b_runnum;   //!
   TBranch        *b_luminum;   //!
   TBranch        *b_eventnum;   //!
   TBranch        *b_NPV;   //!
   TBranch        *b_EventFilter;   //!
   TBranch        *b_FastSimEventVeto;   //!
   TBranch        *b_PrefireWeight;   //!
   TBranch        *b_PrefireWeight_up;   //!
   TBranch        *b_PrefireWeight_down;   //!
   TBranch        *b_METtrigger;   //!
   TBranch        *b_METORtrigger;   //!
   TBranch        *b_DoubleElectrontrigger;   //!
   TBranch        *b_DoubleMuontrigger;   //!
   TBranch        *b_EventFlag_FailJetID;   //!
   TBranch        *b_EventFlag_JetInHEM;   //!
   TBranch        *b_EventFlag_JetInHEM_Pt20;   //!
   TBranch        *b_EventFlag_JetInHEM_Pt20_JetID;   //!
   TBranch        *b_HEM_Veto;   //!
   TBranch        *b_SingleElectrontrigger;   //!
   TBranch        *b_SingleMuontrigger;   //!
   TBranch        *b_EMutrigger;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_altMET;   //!
   TBranch        *b_altMET_phi;   //!
   TBranch        *b_LHE_HT;   //!
   TBranch        *b_LHE_HTIncoming;   //!
   TBranch        *b_HT_eta24;   //!
   TBranch        *b_HT_eta24_id;   //!
   TBranch        *b_HT_eta3;   //!
   TBranch        *b_HT_eta3_id;   //!
   TBranch        *b_HT_eta5;   //!
   TBranch        *b_HT_eta5_id;   //!
   TBranch        *b_Nele;   //!
   TBranch        *b_Nlele;   //!
   TBranch        *b_Nmu;   //!
   TBranch        *b_Nlep;   //!
   TBranch        *b_PT_lep;   //!
   TBranch        *b_Eta_lep;   //!
   TBranch        *b_Phi_lep;   //!
   TBranch        *b_M_lep;   //!
   TBranch        *b_Charge_lep;   //!
   TBranch        *b_PDGID_lep;   //!
   TBranch        *b_ID_lep;   //!
   TBranch        *b_SourceID_lep;   //!
   TBranch        *b_LepQual_lep;   //!
   TBranch        *b_IsLowPt_lep;   //!
   TBranch        *b_TightCharge_lep;   //!
   TBranch        *b_Index_lep;   //!
   TBranch        *b_RelIso_lep;   //!
   TBranch        *b_MiniIso_lep;   //!
   TBranch        *b_Dxy_lep;   //!
   TBranch        *b_DxyErr_lep;   //!
   TBranch        *b_Dz_lep;   //!
   TBranch        *b_DzErr_lep;   //!
   TBranch        *b_IP3D_lep;   //!
   TBranch        *b_SIP3D_lep;   //!
   TBranch        *b_Njet;   //!
   TBranch        *b_Nbjet;   //!
   TBranch        *b_PT_jet;   //!
   TBranch        *b_Eta_jet;   //!
   TBranch        *b_Phi_jet;   //!
   TBranch        *b_M_jet;   //!
   TBranch        *b_Btag_jet;   //!
   TBranch        *b_BtagID_jet;   //!
   TBranch        *b_Flavor_jet;   //!
   TBranch        *b_index_jet_a;   //!
   TBranch        *b_index_jet_b;   //!
   TBranch        *b_index_jet_ISR;   //!
   TBranch        *b_index_jet_S;   //!
   TBranch        *b_PT_Genjet;   //!
   TBranch        *b_Eta_Genjet;   //!
   TBranch        *b_Phi_Genjet;   //!
   TBranch        *b_M_Genjet;   //!
   TBranch        *b_Index_jet;   //!
   TBranch        *b_Njet_ISR;   //!
   TBranch        *b_Njet_S;   //!
   TBranch        *b_Nbjet_ISR;   //!
   TBranch        *b_Nbjet_S;   //!
   TBranch        *b_Nlep_ISR;   //!
   TBranch        *b_Nlep_S;   //!
   TBranch        *b_index_lep_ISR;   //!
   TBranch        *b_index_lep_S;   //!
   TBranch        *b_dphi_lep_S;   //!
   TBranch        *b_cos_lep_S;   //!
   TBranch        *b_dphi_jet_S;   //!
   TBranch        *b_cos_jet_S;   //!
   TBranch        *b_dphiMET_lep_S;   //!
   TBranch        *b_dphiMET_jet_S;   //!
   TBranch        *b_Njet_a;   //!
   TBranch        *b_Njet_b;   //!
   TBranch        *b_Nbjet_a;   //!
   TBranch        *b_Nbjet_b;   //!
   TBranch        *b_Nlep_a;   //!
   TBranch        *b_Nlep_b;   //!
   TBranch        *b_index_lep_a;   //!
   TBranch        *b_index_lep_b;   //!
   TBranch        *b_PTCM;   //!
   TBranch        *b_PzCM;   //!
   TBranch        *b_cosCM;   //!
   TBranch        *b_dphiCM;   //!
   TBranch        *b_dphiCMI;   //!
   TBranch        *b_dphiMET_V;   //!
   TBranch        *b_Mperp;   //!
   TBranch        *b_gammaT;   //!
   TBranch        *b_EJ_BoostT;   //!
   TBranch        *b_EL_BoostT;   //!
   TBranch        *b_PTISR;   //!
   TBranch        *b_RISR;   //!
   TBranch        *b_EtaCM;   //!
   TBranch        *b_PhiCM;   //!
   TBranch        *b_MCM;   //!
   TBranch        *b_EtaS;   //!
   TBranch        *b_PhiS;   //!
   TBranch        *b_LAB_Pt;   //!
   TBranch        *b_LAB_Eta;   //!
   TBranch        *b_LAB_Phi;   //!
   TBranch        *b_LAB_M;   //!
   TBranch        *b_MSperpCM0;   //!
   TBranch        *b_MaPerpCM0;   //!
   TBranch        *b_MbPerpCM0;   //!
   TBranch        *b_MaVPerpCM0;   //!
   TBranch        *b_MbVPerpCM0;   //!
   TBranch        *b_MQperpCM0;   //!
   TBranch        *b_gammaPerpCM0;   //!
   TBranch        *b_MVisAperpCM0;   //!
   TBranch        *b_MVisBperpCM0;   //!
   TBranch        *b_MSCM0;   //!
   TBranch        *b_MaCM0;   //!
   TBranch        *b_MbCM0;   //!
   TBranch        *b_MaVCM0;   //!
   TBranch        *b_MbVCM0;   //!
   TBranch        *b_MQCM0;   //!
   TBranch        *b_gammaCM0;   //!
   TBranch        *b_MVisACM0;   //!
   TBranch        *b_MVisBCM0;   //!
   TBranch        *b_MS;   //!
   TBranch        *b_PS;   //!
   TBranch        *b_cosS;   //!
   TBranch        *b_dphiS;   //!
   TBranch        *b_dphiSI;   //!
   TBranch        *b_PTS;   //!
   TBranch        *b_PzS;   //!
   TBranch        *b_EVa;   //!
   TBranch        *b_EVb;   //!
   TBranch        *b_PVa;   //!
   TBranch        *b_PVb;   //!
   TBranch        *b_EJa;   //!
   TBranch        *b_EJb;   //!
   TBranch        *b_PJa;   //!
   TBranch        *b_PJb;   //!
   TBranch        *b_MX2a;   //!
   TBranch        *b_MX2b;   //!
   TBranch        *b_ELa;   //!
   TBranch        *b_ELb;   //!
   TBranch        *b_PLa;   //!
   TBranch        *b_PLb;   //!
   TBranch        *b_MV;   //!
   TBranch        *b_PV;   //!
   TBranch        *b_MVa;   //!
   TBranch        *b_MVb;   //!
   TBranch        *b_PV_lab;   //!
   TBranch        *b_MJa;   //!
   TBranch        *b_MJb;   //!
   TBranch        *b_MLa;   //!
   TBranch        *b_MLb;   //!
   TBranch        *b_cosJa;   //!
   TBranch        *b_cosJb;   //!
   TBranch        *b_cosLa;   //!
   TBranch        *b_cosLb;   //!
   TBranch        *b_MJ;   //!
   TBranch        *b_ML;   //!
   TBranch        *b_EJ;   //!
   TBranch        *b_EL;   //!
   TBranch        *b_PJ;   //!
   TBranch        *b_PL;   //!
   TBranch        *b_PX2;   //!
   TBranch        *b_PX2_BoostT;   //!
   TBranch        *b_MX2a_BoostT;   //!
   TBranch        *b_MX2b_BoostT;   //!
   TBranch        *b_PV_BoostT;   //!
   TBranch        *b_EVa_BoostT;   //!
   TBranch        *b_EVb_BoostT;   //!
   TBranch        *b_PVa_BoostT;   //!
   TBranch        *b_PVb_BoostT;   //!
   TBranch        *b_PJ_BoostT;   //!
   TBranch        *b_PL_BoostT;   //!
   TBranch        *b_PISR;   //!
   TBranch        *b_RISRT;   //!
   TBranch        *b_MISR;   //!
   TBranch        *b_NPU;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMET_phi;   //!
   TBranch        *b_genNele;   //!
   TBranch        *b_genNmu;   //!
   TBranch        *b_genNlep;   //!
   TBranch        *b_genPT_lep;   //!
   TBranch        *b_genEta_lep;   //!
   TBranch        *b_genPhi_lep;   //!
   TBranch        *b_genM_lep;   //!
   TBranch        *b_genCharge_lep;   //!
   TBranch        *b_genPDGID_lep;   //!
   TBranch        *b_genMomPDGID_lep;   //!
   TBranch        *b_genSourceID_lep;   //!
   TBranch        *b_genIndex_lep;   //!
   TBranch        *b_genMomIndex_lep;   //!
   TBranch        *b_genNnu;   //!
   TBranch        *b_genPT_nu;   //!
   TBranch        *b_genEta_nu;   //!
   TBranch        *b_genPhi_nu;   //!
   TBranch        *b_genPDGID_nu;   //!
   TBranch        *b_genMomPDGID_nu;   //!
   TBranch        *b_genNboson;   //!
   TBranch        *b_genPT_boson;   //!
   TBranch        *b_genEta_boson;   //!
   TBranch        *b_genPhi_boson;   //!
   TBranch        *b_genM_boson;   //!
   TBranch        *b_genPDGID_boson;   //!
   TBranch        *b_genMomPDGID_boson;   //!
   TBranch        *b_genNsusy;   //!
   TBranch        *b_genPT_susy;   //!
   TBranch        *b_genEta_susy;   //!
   TBranch        *b_genPhi_susy;   //!
   TBranch        *b_genM_susy;   //!
   TBranch        *b_genPDGID_susy;   //!
   TBranch        *b_genMomPDGID_susy;   //!
   TBranch        *b_MT2;   //!
   TBranch        *b_RISR_LEP;   //!
   TBranch        *b_PTISR_LEP;   //!
   TBranch        *b_MS_LEP;   //!
   TBranch        *b_MSV_LEP;   //!
   TBranch        *b_MQ_LEP;   //!
   TBranch        *b_gamma_LEP;   //!
   TBranch        *b_MPa_LEP;   //!
   TBranch        *b_MPb_LEP;   //!
   TBranch        *b_MVa_LEP;   //!
   TBranch        *b_MVb_LEP;   //!
   TBranch        *b_PTS_CM_LEP;   //!
   TBranch        *b_MS_S0_LEP;   //!
   TBranch        *b_MV_S0_LEP;   //!
   TBranch        *b_MQ_S0_LEP;   //!
   TBranch        *b_gamma_S0_LEP;   //!
   TBranch        *b_MPTilde_LEP;   //!
   TBranch        *b_MSTilde_LEP;   //!
   TBranch        *b_gammaTilde_LEP;   //!
   TBranch        *b_CosDecayAngle_Pa_LEP;   //!
   TBranch        *b_CosDecayAngle_Pb_LEP;   //!
   TBranch        *b_CosDecayAngle_S_LEP;   //!
   TBranch        *b_RZPara_LEP;   //!
   TBranch        *b_MINV_LEP;   //!
   TBranch        *b_MP;   //!
   TBranch        *b_Mperp_JET;   //!
   TBranch        *b_gammaT_JET;   //!
   TBranch        *b_RISR_JET;   //!
   TBranch        *b_PTISR_JET;   //!
   TBranch        *b_MS_JET;   //!
   TBranch        *b_MSV_JET;   //!
   TBranch        *b_MQ_JET;   //!
   TBranch        *b_gamma_JET;   //!
   TBranch        *b_MPa_JET;   //!
   TBranch        *b_MPb_JET;   //!
   TBranch        *b_MVa_JET;   //!
   TBranch        *b_MVb_JET;   //!
   TBranch        *b_PTS_CM_JET;   //!
   TBranch        *b_MS_S0_JET;   //!
   TBranch        *b_MV_S0_JET;   //!
   TBranch        *b_MQ_S0_JET;   //!
   TBranch        *b_gamma_S0_JET;   //!
   TBranch        *b_MPTilde_JET;   //!
   TBranch        *b_MSTilde_JET;   //!
   TBranch        *b_gammaTilde_JET;   //!
   TBranch        *b_CosDecayAngle_Pa_JET;   //!
   TBranch        *b_CosDecayAngle_Pb_JET;   //!
   TBranch        *b_CosDecayAngle_S_JET;   //!
   TBranch        *b_RZPara_JET;   //!
   TBranch        *b_MINV_JET;   //!
   TBranch        *b_cascades_tree;   //!
   TBranch        *b_cascades_prod;   //!
   TBranch        *b_cascades_SlepSneu_1stDecay;   //!
   TBranch        *b_cascades_SlepSneu_2ndDecay;   //!
   TBranch        *b_cascades_N2_1stDecay;   //!
   TBranch        *b_cascades_N2_2ndDecay;   //!
   TBranch        *b_MSlepL;   //!
   TBranch        *b_MSneu;   //!
   TBranch        *b_MN2;   //!
   TBranch        *b_MC1;   //!
   TBranch        *b_MN1;   //!
   TBranch        *b_NSparticleW;   //!
   TBranch        *b_NSparticleZ;   //!
   TBranch        *b_NSparticlePhoton;   //!
   TBranch        *b_LSPParents;   //!
   TBranch        *b_Npartons;   //!
   TBranch        *b_XSec;   //!
   TBranch        *b_Nweight;   //!
   TBranch        *b_Nevent;   //!

   ReducedBase_V2(TTree *tree=0);
   virtual ~ReducedBase_V2();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline ReducedBase_V2::ReducedBase_V2(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/local-scratch/zflowers/NTUPLES/HADD/Summer23BPix_130X/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Summer23BPix_130X.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/local-scratch/zflowers/NTUPLES/HADD/Summer23BPix_130X/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Summer23BPix_130X.root");
      }
      f->GetObject("KUAnalysis",tree);

   }
   Init(tree);
}

inline ReducedBase_V2::~ReducedBase_V2()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t ReducedBase_V2::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t ReducedBase_V2::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

inline void ReducedBase_V2::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   treeSkipped = 0;
   PT_lep = 0;
   Eta_lep = 0;
   Phi_lep = 0;
   M_lep = 0;
   Charge_lep = 0;
   PDGID_lep = 0;
   ID_lep = 0;
   SourceID_lep = 0;
   LepQual_lep = 0;
   IsLowPt_lep = 0;
   TightCharge_lep = 0;
   Index_lep = 0;
   RelIso_lep = 0;
   MiniIso_lep = 0;
   Dxy_lep = 0;
   DxyErr_lep = 0;
   Dz_lep = 0;
   DzErr_lep = 0;
   IP3D_lep = 0;
   SIP3D_lep = 0;
   PT_jet = 0;
   Eta_jet = 0;
   Phi_jet = 0;
   M_jet = 0;
   Btag_jet = 0;
   BtagID_jet = 0;
   Flavor_jet = 0;
   index_jet_a = 0;
   index_jet_b = 0;
   index_jet_ISR = 0;
   index_jet_S = 0;
   PT_Genjet = 0;
   Eta_Genjet = 0;
   Phi_Genjet = 0;
   M_Genjet = 0;
   Index_jet = 0;
   index_lep_ISR = 0;
   index_lep_S = 0;
   dphi_lep_S = 0;
   cos_lep_S = 0;
   dphi_jet_S = 0;
   cos_jet_S = 0;
   dphiMET_lep_S = 0;
   dphiMET_jet_S = 0;
   index_lep_a = 0;
   index_lep_b = 0;
   genPT_lep = 0;
   genEta_lep = 0;
   genPhi_lep = 0;
   genM_lep = 0;
   genCharge_lep = 0;
   genPDGID_lep = 0;
   genMomPDGID_lep = 0;
   genSourceID_lep = 0;
   genIndex_lep = 0;
   genMomIndex_lep = 0;
   genPT_nu = 0;
   genEta_nu = 0;
   genPhi_nu = 0;
   genPDGID_nu = 0;
   genMomPDGID_nu = 0;
   genPT_boson = 0;
   genEta_boson = 0;
   genPhi_boson = 0;
   genM_boson = 0;
   genPDGID_boson = 0;
   genMomPDGID_boson = 0;
   genPT_susy = 0;
   genEta_susy = 0;
   genPhi_susy = 0;
   genM_susy = 0;
   genPDGID_susy = 0;
   genMomPDGID_susy = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_skipped", &event_skipped, &b_event_skipped);
   fChain->SetBranchAddress("treeSkipped", &treeSkipped, &b_treeSkipped);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweight_up", &PUweight_up, &b_PUweight_up);
   fChain->SetBranchAddress("PUweight_down", &PUweight_down, &b_PUweight_down);
   fChain->SetBranchAddress("MuFweight", &MuFweight, &b_MuFweight);
   fChain->SetBranchAddress("MuFweight_up", &MuFweight_up, &b_MuFweight_up);
   fChain->SetBranchAddress("MuFweight_down", &MuFweight_down, &b_MuFweight_down);
   fChain->SetBranchAddress("MuRweight", &MuRweight, &b_MuRweight);
   fChain->SetBranchAddress("MuRweight_up", &MuRweight_up, &b_MuRweight_up);
   fChain->SetBranchAddress("MuRweight_down", &MuRweight_down, &b_MuRweight_down);
   fChain->SetBranchAddress("PDFweight", &PDFweight, &b_PDFweight);
   fChain->SetBranchAddress("PDFweight_up", &PDFweight_up, &b_PDFweight_up);
   fChain->SetBranchAddress("PDFweight_down", &PDFweight_down, &b_PDFweight_down);
   fChain->SetBranchAddress("BtagHFSFweight", &BtagHFSFweight, &b_BtagHFSFweight);
   fChain->SetBranchAddress("BtagHFSFweight_up", &BtagHFSFweight_up, &b_BtagHFSFweight_up);
   fChain->SetBranchAddress("BtagHFSFweight_down", &BtagHFSFweight_down, &b_BtagHFSFweight_down);
   fChain->SetBranchAddress("BtagLFSFweight", &BtagLFSFweight, &b_BtagLFSFweight);
   fChain->SetBranchAddress("BtagLFSFweight_up", &BtagLFSFweight_up, &b_BtagLFSFweight_up);
   fChain->SetBranchAddress("BtagLFSFweight_down", &BtagLFSFweight_down, &b_BtagLFSFweight_down);
   fChain->SetBranchAddress("elIDSFweight", &elIDSFweight, &b_elIDSFweight);
   fChain->SetBranchAddress("elIDSFweight_up", &elIDSFweight_up, &b_elIDSFweight_up);
   fChain->SetBranchAddress("elIDSFweight_down", &elIDSFweight_down, &b_elIDSFweight_down);
   fChain->SetBranchAddress("elISOSFweight", &elISOSFweight, &b_elISOSFweight);
   fChain->SetBranchAddress("elISOSFweight_up", &elISOSFweight_up, &b_elISOSFweight_up);
   fChain->SetBranchAddress("elISOSFweight_down", &elISOSFweight_down, &b_elISOSFweight_down);
   fChain->SetBranchAddress("elSIPSFweight", &elSIPSFweight, &b_elSIPSFweight);
   fChain->SetBranchAddress("elSIPSFweight_up", &elSIPSFweight_up, &b_elSIPSFweight_up);
   fChain->SetBranchAddress("elSIPSFweight_down", &elSIPSFweight_down, &b_elSIPSFweight_down);
   fChain->SetBranchAddress("elVLSFweight", &elVLSFweight, &b_elVLSFweight);
   fChain->SetBranchAddress("elVLSFweight_up", &elVLSFweight_up, &b_elVLSFweight_up);
   fChain->SetBranchAddress("elVLSFweight_down", &elVLSFweight_down, &b_elVLSFweight_down);
   fChain->SetBranchAddress("muIDSFweight", &muIDSFweight, &b_muIDSFweight);
   fChain->SetBranchAddress("muIDSFweight_up", &muIDSFweight_up, &b_muIDSFweight_up);
   fChain->SetBranchAddress("muIDSFweight_down", &muIDSFweight_down, &b_muIDSFweight_down);
   fChain->SetBranchAddress("muISOSFweight", &muISOSFweight, &b_muISOSFweight);
   fChain->SetBranchAddress("muISOSFweight_up", &muISOSFweight_up, &b_muISOSFweight_up);
   fChain->SetBranchAddress("muISOSFweight_down", &muISOSFweight_down, &b_muISOSFweight_down);
   fChain->SetBranchAddress("muSIPSFweight", &muSIPSFweight, &b_muSIPSFweight);
   fChain->SetBranchAddress("muSIPSFweight_up", &muSIPSFweight_up, &b_muSIPSFweight_up);
   fChain->SetBranchAddress("muSIPSFweight_down", &muSIPSFweight_down, &b_muSIPSFweight_down);
   fChain->SetBranchAddress("muVLSFweight", &muVLSFweight, &b_muVLSFweight);
   fChain->SetBranchAddress("muVLSFweight_up", &muVLSFweight_up, &b_muVLSFweight_up);
   fChain->SetBranchAddress("muVLSFweight_down", &muVLSFweight_down, &b_muVLSFweight_down);
   fChain->SetBranchAddress("MetTrigSFweight", &MetTrigSFweight, &b_MetTrigSFweight);
   fChain->SetBranchAddress("MetTrigSFweight_up", &MetTrigSFweight_up, &b_MetTrigSFweight_up);
   fChain->SetBranchAddress("MetTrigSFweight_down", &MetTrigSFweight_down, &b_MetTrigSFweight_down);
   fChain->SetBranchAddress("MetTrigSFCurveIndex", &MetTrigSFCurveIndex, &b_MetTrigSFCurveIndex);
   fChain->SetBranchAddress("runnum", &runnum, &b_runnum);
   fChain->SetBranchAddress("luminum", &luminum, &b_luminum);
   fChain->SetBranchAddress("eventnum", &eventnum, &b_eventnum);
   fChain->SetBranchAddress("NPV", &NPV, &b_NPV);
   fChain->SetBranchAddress("EventFilter", &EventFilter, &b_EventFilter);
   fChain->SetBranchAddress("FastSimEventVeto", &FastSimEventVeto, &b_FastSimEventVeto);
   fChain->SetBranchAddress("PrefireWeight", &PrefireWeight, &b_PrefireWeight);
   fChain->SetBranchAddress("PrefireWeight_up", &PrefireWeight_up, &b_PrefireWeight_up);
   fChain->SetBranchAddress("PrefireWeight_down", &PrefireWeight_down, &b_PrefireWeight_down);
   fChain->SetBranchAddress("METtrigger", &METtrigger, &b_METtrigger);
   fChain->SetBranchAddress("METORtrigger", &METORtrigger, &b_METORtrigger);
   fChain->SetBranchAddress("DoubleElectrontrigger", &DoubleElectrontrigger, &b_DoubleElectrontrigger);
   fChain->SetBranchAddress("DoubleMuontrigger", &DoubleMuontrigger, &b_DoubleMuontrigger);
   fChain->SetBranchAddress("EventFlag_FailJetID", &EventFlag_FailJetID, &b_EventFlag_FailJetID);
   fChain->SetBranchAddress("EventFlag_JetInHEM", &EventFlag_JetInHEM, &b_EventFlag_JetInHEM);
   fChain->SetBranchAddress("EventFlag_JetInHEM_Pt20", &EventFlag_JetInHEM_Pt20, &b_EventFlag_JetInHEM_Pt20);
   fChain->SetBranchAddress("EventFlag_JetInHEM_Pt20_JetID", &EventFlag_JetInHEM_Pt20_JetID, &b_EventFlag_JetInHEM_Pt20_JetID);
   fChain->SetBranchAddress("HEM_Veto", &HEM_Veto, &b_HEM_Veto);
   fChain->SetBranchAddress("SingleElectrontrigger", &SingleElectrontrigger, &b_SingleElectrontrigger);
   fChain->SetBranchAddress("SingleMuontrigger", &SingleMuontrigger, &b_SingleMuontrigger);
   fChain->SetBranchAddress("EMutrigger", &EMutrigger, &b_EMutrigger);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("altMET", &altMET, &b_altMET);
   fChain->SetBranchAddress("altMET_phi", &altMET_phi, &b_altMET_phi);
   fChain->SetBranchAddress("LHE_HT", &LHE_HT, &b_LHE_HT);
   fChain->SetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming, &b_LHE_HTIncoming);
   fChain->SetBranchAddress("HT_eta24", &HT_eta24, &b_HT_eta24);
   fChain->SetBranchAddress("HT_eta24_id", &HT_eta24_id, &b_HT_eta24_id);
   fChain->SetBranchAddress("HT_eta3", &HT_eta3, &b_HT_eta3);
   fChain->SetBranchAddress("HT_eta3_id", &HT_eta3_id, &b_HT_eta3_id);
   fChain->SetBranchAddress("HT_eta5", &HT_eta5, &b_HT_eta5);
   fChain->SetBranchAddress("HT_eta5_id", &HT_eta5_id, &b_HT_eta5_id);
   fChain->SetBranchAddress("Nele", &Nele, &b_Nele);
   fChain->SetBranchAddress("Nlele", &Nlele, &b_Nlele);
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("Nlep", &Nlep, &b_Nlep);
   fChain->SetBranchAddress("PT_lep", &PT_lep, &b_PT_lep);
   fChain->SetBranchAddress("Eta_lep", &Eta_lep, &b_Eta_lep);
   fChain->SetBranchAddress("Phi_lep", &Phi_lep, &b_Phi_lep);
   fChain->SetBranchAddress("M_lep", &M_lep, &b_M_lep);
   fChain->SetBranchAddress("Charge_lep", &Charge_lep, &b_Charge_lep);
   fChain->SetBranchAddress("PDGID_lep", &PDGID_lep, &b_PDGID_lep);
   fChain->SetBranchAddress("ID_lep", &ID_lep, &b_ID_lep);
   fChain->SetBranchAddress("SourceID_lep", &SourceID_lep, &b_SourceID_lep);
   fChain->SetBranchAddress("LepQual_lep", &LepQual_lep, &b_LepQual_lep);
   fChain->SetBranchAddress("IsLowPt_lep", &IsLowPt_lep, &b_IsLowPt_lep);
   fChain->SetBranchAddress("TightCharge_lep", &TightCharge_lep, &b_TightCharge_lep);
   fChain->SetBranchAddress("Index_lep", &Index_lep, &b_Index_lep);
   fChain->SetBranchAddress("RelIso_lep", &RelIso_lep, &b_RelIso_lep);
   fChain->SetBranchAddress("MiniIso_lep", &MiniIso_lep, &b_MiniIso_lep);
   fChain->SetBranchAddress("Dxy_lep", &Dxy_lep, &b_Dxy_lep);
   fChain->SetBranchAddress("DxyErr_lep", &DxyErr_lep, &b_DxyErr_lep);
   fChain->SetBranchAddress("Dz_lep", &Dz_lep, &b_Dz_lep);
   fChain->SetBranchAddress("DzErr_lep", &DzErr_lep, &b_DzErr_lep);
   fChain->SetBranchAddress("IP3D_lep", &IP3D_lep, &b_IP3D_lep);
   fChain->SetBranchAddress("SIP3D_lep", &SIP3D_lep, &b_SIP3D_lep);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("Nbjet", &Nbjet, &b_Nbjet);
   fChain->SetBranchAddress("PT_jet", &PT_jet, &b_PT_jet);
   fChain->SetBranchAddress("Eta_jet", &Eta_jet, &b_Eta_jet);
   fChain->SetBranchAddress("Phi_jet", &Phi_jet, &b_Phi_jet);
   fChain->SetBranchAddress("M_jet", &M_jet, &b_M_jet);
   fChain->SetBranchAddress("Btag_jet", &Btag_jet, &b_Btag_jet);
   fChain->SetBranchAddress("BtagID_jet", &BtagID_jet, &b_BtagID_jet);
   fChain->SetBranchAddress("Flavor_jet", &Flavor_jet, &b_Flavor_jet);
   fChain->SetBranchAddress("index_jet_a", &index_jet_a, &b_index_jet_a);
   fChain->SetBranchAddress("index_jet_b", &index_jet_b, &b_index_jet_b);
   fChain->SetBranchAddress("index_jet_ISR", &index_jet_ISR, &b_index_jet_ISR);
   fChain->SetBranchAddress("index_jet_S", &index_jet_S, &b_index_jet_S);
   fChain->SetBranchAddress("PT_Genjet", &PT_Genjet, &b_PT_Genjet);
   fChain->SetBranchAddress("Eta_Genjet", &Eta_Genjet, &b_Eta_Genjet);
   fChain->SetBranchAddress("Phi_Genjet", &Phi_Genjet, &b_Phi_Genjet);
   fChain->SetBranchAddress("M_Genjet", &M_Genjet, &b_M_Genjet);
   fChain->SetBranchAddress("Index_jet", &Index_jet, &b_Index_jet);
   fChain->SetBranchAddress("Njet_ISR", &Njet_ISR, &b_Njet_ISR);
   fChain->SetBranchAddress("Njet_S", &Njet_S, &b_Njet_S);
   fChain->SetBranchAddress("Nbjet_ISR", &Nbjet_ISR, &b_Nbjet_ISR);
   fChain->SetBranchAddress("Nbjet_S", &Nbjet_S, &b_Nbjet_S);
   fChain->SetBranchAddress("Nlep_ISR", &Nlep_ISR, &b_Nlep_ISR);
   fChain->SetBranchAddress("Nlep_S", &Nlep_S, &b_Nlep_S);
   fChain->SetBranchAddress("index_lep_ISR", &index_lep_ISR, &b_index_lep_ISR);
   fChain->SetBranchAddress("index_lep_S", &index_lep_S, &b_index_lep_S);
   fChain->SetBranchAddress("dphi_lep_S", &dphi_lep_S, &b_dphi_lep_S);
   fChain->SetBranchAddress("cos_lep_S", &cos_lep_S, &b_cos_lep_S);
   fChain->SetBranchAddress("dphi_jet_S", &dphi_jet_S, &b_dphi_jet_S);
   fChain->SetBranchAddress("cos_jet_S", &cos_jet_S, &b_cos_jet_S);
   fChain->SetBranchAddress("dphiMET_lep_S", &dphiMET_lep_S, &b_dphiMET_lep_S);
   fChain->SetBranchAddress("dphiMET_jet_S", &dphiMET_jet_S, &b_dphiMET_jet_S);
   fChain->SetBranchAddress("Njet_a", &Njet_a, &b_Njet_a);
   fChain->SetBranchAddress("Njet_b", &Njet_b, &b_Njet_b);
   fChain->SetBranchAddress("Nbjet_a", &Nbjet_a, &b_Nbjet_a);
   fChain->SetBranchAddress("Nbjet_b", &Nbjet_b, &b_Nbjet_b);
   fChain->SetBranchAddress("Nlep_a", &Nlep_a, &b_Nlep_a);
   fChain->SetBranchAddress("Nlep_b", &Nlep_b, &b_Nlep_b);
   fChain->SetBranchAddress("index_lep_a", &index_lep_a, &b_index_lep_a);
   fChain->SetBranchAddress("index_lep_b", &index_lep_b, &b_index_lep_b);
   fChain->SetBranchAddress("PTCM", &PTCM, &b_PTCM);
   fChain->SetBranchAddress("PzCM", &PzCM, &b_PzCM);
   fChain->SetBranchAddress("cosCM", &cosCM, &b_cosCM);
   fChain->SetBranchAddress("dphiCM", &dphiCM, &b_dphiCM);
   fChain->SetBranchAddress("dphiCMI", &dphiCMI, &b_dphiCMI);
   fChain->SetBranchAddress("dphiMET_V", &dphiMET_V, &b_dphiMET_V);
   fChain->SetBranchAddress("Mperp", &Mperp, &b_Mperp);
   fChain->SetBranchAddress("gammaT", &gammaT, &b_gammaT);
   fChain->SetBranchAddress("EJ_BoostT", &EJ_BoostT, &b_EJ_BoostT);
   fChain->SetBranchAddress("EL_BoostT", &EL_BoostT, &b_EL_BoostT);
   fChain->SetBranchAddress("PTISR", &PTISR, &b_PTISR);
   fChain->SetBranchAddress("RISR", &RISR, &b_RISR);
   fChain->SetBranchAddress("EtaCM", &EtaCM, &b_EtaCM);
   fChain->SetBranchAddress("PhiCM", &PhiCM, &b_PhiCM);
   fChain->SetBranchAddress("MCM", &MCM, &b_MCM);
   fChain->SetBranchAddress("EtaS", &EtaS, &b_EtaS);
   fChain->SetBranchAddress("PhiS", &PhiS, &b_PhiS);
   fChain->SetBranchAddress("LAB_Pt", &LAB_Pt, &b_LAB_Pt);
   fChain->SetBranchAddress("LAB_Eta", &LAB_Eta, &b_LAB_Eta);
   fChain->SetBranchAddress("LAB_Phi", &LAB_Phi, &b_LAB_Phi);
   fChain->SetBranchAddress("LAB_M", &LAB_M, &b_LAB_M);
   fChain->SetBranchAddress("MSperpCM0", &MSperpCM0, &b_MSperpCM0);
   fChain->SetBranchAddress("MaPerpCM0", &MaPerpCM0, &b_MaPerpCM0);
   fChain->SetBranchAddress("MbPerpCM0", &MbPerpCM0, &b_MbPerpCM0);
   fChain->SetBranchAddress("MaVPerpCM0", &MaVPerpCM0, &b_MaVPerpCM0);
   fChain->SetBranchAddress("MbVPerpCM0", &MbVPerpCM0, &b_MbVPerpCM0);
   fChain->SetBranchAddress("MQperpCM0", &MQperpCM0, &b_MQperpCM0);
   fChain->SetBranchAddress("gammaPerpCM0", &gammaPerpCM0, &b_gammaPerpCM0);
   fChain->SetBranchAddress("MVisAperpCM0", &MVisAperpCM0, &b_MVisAperpCM0);
   fChain->SetBranchAddress("MVisBperpCM0", &MVisBperpCM0, &b_MVisBperpCM0);
   fChain->SetBranchAddress("MSCM0", &MSCM0, &b_MSCM0);
   fChain->SetBranchAddress("MaCM0", &MaCM0, &b_MaCM0);
   fChain->SetBranchAddress("MbCM0", &MbCM0, &b_MbCM0);
   fChain->SetBranchAddress("MaVCM0", &MaVCM0, &b_MaVCM0);
   fChain->SetBranchAddress("MbVCM0", &MbVCM0, &b_MbVCM0);
   fChain->SetBranchAddress("MQCM0", &MQCM0, &b_MQCM0);
   fChain->SetBranchAddress("gammaCM0", &gammaCM0, &b_gammaCM0);
   fChain->SetBranchAddress("MVisACM0", &MVisACM0, &b_MVisACM0);
   fChain->SetBranchAddress("MVisBCM0", &MVisBCM0, &b_MVisBCM0);
   fChain->SetBranchAddress("MS", &MS, &b_MS);
   fChain->SetBranchAddress("PS", &PS, &b_PS);
   fChain->SetBranchAddress("cosS", &cosS, &b_cosS);
   fChain->SetBranchAddress("dphiS", &dphiS, &b_dphiS);
   fChain->SetBranchAddress("dphiSI", &dphiSI, &b_dphiSI);
   fChain->SetBranchAddress("PTS", &PTS, &b_PTS);
   fChain->SetBranchAddress("PzS", &PzS, &b_PzS);
   fChain->SetBranchAddress("EVa", &EVa, &b_EVa);
   fChain->SetBranchAddress("EVb", &EVb, &b_EVb);
   fChain->SetBranchAddress("PVa", &PVa, &b_PVa);
   fChain->SetBranchAddress("PVb", &PVb, &b_PVb);
   fChain->SetBranchAddress("EJa", &EJa, &b_EJa);
   fChain->SetBranchAddress("EJb", &EJb, &b_EJb);
   fChain->SetBranchAddress("PJa", &PJa, &b_PJa);
   fChain->SetBranchAddress("PJb", &PJb, &b_PJb);
   fChain->SetBranchAddress("MX2a", &MX2a, &b_MX2a);
   fChain->SetBranchAddress("MX2b", &MX2b, &b_MX2b);
   fChain->SetBranchAddress("ELa", &ELa, &b_ELa);
   fChain->SetBranchAddress("ELb", &ELb, &b_ELb);
   fChain->SetBranchAddress("PLa", &PLa, &b_PLa);
   fChain->SetBranchAddress("PLb", &PLb, &b_PLb);
   fChain->SetBranchAddress("MV", &MV, &b_MV);
   fChain->SetBranchAddress("PV", &PV, &b_PV);
   fChain->SetBranchAddress("MVa", &MVa, &b_MVa);
   fChain->SetBranchAddress("MVb", &MVb, &b_MVb);
   fChain->SetBranchAddress("PV_lab", &PV_lab, &b_PV_lab);
   fChain->SetBranchAddress("MJa", &MJa, &b_MJa);
   fChain->SetBranchAddress("MJb", &MJb, &b_MJb);
   fChain->SetBranchAddress("MLa", &MLa, &b_MLa);
   fChain->SetBranchAddress("MLb", &MLb, &b_MLb);
   fChain->SetBranchAddress("cosJa", &cosJa, &b_cosJa);
   fChain->SetBranchAddress("cosJb", &cosJb, &b_cosJb);
   fChain->SetBranchAddress("cosLa", &cosLa, &b_cosLa);
   fChain->SetBranchAddress("cosLb", &cosLb, &b_cosLb);
   fChain->SetBranchAddress("MJ", &MJ, &b_MJ);
   fChain->SetBranchAddress("ML", &ML, &b_ML);
   fChain->SetBranchAddress("EJ", &EJ, &b_EJ);
   fChain->SetBranchAddress("EL", &EL, &b_EL);
   fChain->SetBranchAddress("PJ", &PJ, &b_PJ);
   fChain->SetBranchAddress("PL", &PL, &b_PL);
   fChain->SetBranchAddress("PX2", &PX2, &b_PX2);
   fChain->SetBranchAddress("PX2_BoostT", &PX2_BoostT, &b_PX2_BoostT);
   fChain->SetBranchAddress("MX2a_BoostT", &MX2a_BoostT, &b_MX2a_BoostT);
   fChain->SetBranchAddress("MX2b_BoostT", &MX2b_BoostT, &b_MX2b_BoostT);
   fChain->SetBranchAddress("PV_BoostT", &PV_BoostT, &b_PV_BoostT);
   fChain->SetBranchAddress("EVa_BoostT", &EVa_BoostT, &b_EVa_BoostT);
   fChain->SetBranchAddress("EVb_BoostT", &EVb_BoostT, &b_EVb_BoostT);
   fChain->SetBranchAddress("PVa_BoostT", &PVa_BoostT, &b_PVa_BoostT);
   fChain->SetBranchAddress("PVb_BoostT", &PVb_BoostT, &b_PVb_BoostT);
   fChain->SetBranchAddress("PJ_BoostT", &PJ_BoostT, &b_PJ_BoostT);
   fChain->SetBranchAddress("PL_BoostT", &PL_BoostT, &b_PL_BoostT);
   fChain->SetBranchAddress("PISR", &PISR, &b_PISR);
   fChain->SetBranchAddress("RISRT", &RISRT, &b_RISRT);
   fChain->SetBranchAddress("MISR", &MISR, &b_MISR);
   fChain->SetBranchAddress("NPU", &NPU, &b_NPU);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMET_phi", &genMET_phi, &b_genMET_phi);
   fChain->SetBranchAddress("genNele", &genNele, &b_genNele);
   fChain->SetBranchAddress("genNmu", &genNmu, &b_genNmu);
   fChain->SetBranchAddress("genNlep", &genNlep, &b_genNlep);
   fChain->SetBranchAddress("genPT_lep", &genPT_lep, &b_genPT_lep);
   fChain->SetBranchAddress("genEta_lep", &genEta_lep, &b_genEta_lep);
   fChain->SetBranchAddress("genPhi_lep", &genPhi_lep, &b_genPhi_lep);
   fChain->SetBranchAddress("genM_lep", &genM_lep, &b_genM_lep);
   fChain->SetBranchAddress("genCharge_lep", &genCharge_lep, &b_genCharge_lep);
   fChain->SetBranchAddress("genPDGID_lep", &genPDGID_lep, &b_genPDGID_lep);
   fChain->SetBranchAddress("genMomPDGID_lep", &genMomPDGID_lep, &b_genMomPDGID_lep);
   fChain->SetBranchAddress("genSourceID_lep", &genSourceID_lep, &b_genSourceID_lep);
   fChain->SetBranchAddress("genIndex_lep", &genIndex_lep, &b_genIndex_lep);
   fChain->SetBranchAddress("genMomIndex_lep", &genMomIndex_lep, &b_genMomIndex_lep);
   fChain->SetBranchAddress("genNnu", &genNnu, &b_genNnu);
   fChain->SetBranchAddress("genPT_nu", &genPT_nu, &b_genPT_nu);
   fChain->SetBranchAddress("genEta_nu", &genEta_nu, &b_genEta_nu);
   fChain->SetBranchAddress("genPhi_nu", &genPhi_nu, &b_genPhi_nu);
   fChain->SetBranchAddress("genPDGID_nu", &genPDGID_nu, &b_genPDGID_nu);
   fChain->SetBranchAddress("genMomPDGID_nu", &genMomPDGID_nu, &b_genMomPDGID_nu);
   fChain->SetBranchAddress("genNboson", &genNboson, &b_genNboson);
   fChain->SetBranchAddress("genPT_boson", &genPT_boson, &b_genPT_boson);
   fChain->SetBranchAddress("genEta_boson", &genEta_boson, &b_genEta_boson);
   fChain->SetBranchAddress("genPhi_boson", &genPhi_boson, &b_genPhi_boson);
   fChain->SetBranchAddress("genM_boson", &genM_boson, &b_genM_boson);
   fChain->SetBranchAddress("genPDGID_boson", &genPDGID_boson, &b_genPDGID_boson);
   fChain->SetBranchAddress("genMomPDGID_boson", &genMomPDGID_boson, &b_genMomPDGID_boson);
   fChain->SetBranchAddress("genNsusy", &genNsusy, &b_genNsusy);
   fChain->SetBranchAddress("genPT_susy", &genPT_susy, &b_genPT_susy);
   fChain->SetBranchAddress("genEta_susy", &genEta_susy, &b_genEta_susy);
   fChain->SetBranchAddress("genPhi_susy", &genPhi_susy, &b_genPhi_susy);
   fChain->SetBranchAddress("genM_susy", &genM_susy, &b_genM_susy);
   fChain->SetBranchAddress("genPDGID_susy", &genPDGID_susy, &b_genPDGID_susy);
   fChain->SetBranchAddress("genMomPDGID_susy", &genMomPDGID_susy, &b_genMomPDGID_susy);
   fChain->SetBranchAddress("cascades_tree", &cascades_tree, &b_cascades_tree);
   fChain->SetBranchAddress("cascades_prod", &cascades_prod, &b_cascades_prod);
   fChain->SetBranchAddress("cascades_SlepSneu_1stDecay", &cascades_SlepSneu_1stDecay, &b_cascades_SlepSneu_1stDecay);
   fChain->SetBranchAddress("cascades_SlepSneu_2ndDecay", &cascades_SlepSneu_2ndDecay, &b_cascades_SlepSneu_2ndDecay);
   fChain->SetBranchAddress("cascades_N2_1stDecay", &cascades_N2_1stDecay, &b_cascades_N2_1stDecay);
   fChain->SetBranchAddress("cascades_N2_2ndDecay", &cascades_N2_2ndDecay, &b_cascades_N2_2ndDecay);
   fChain->SetBranchAddress("MSlepL", &MSlepL, &b_MSlepL);
   fChain->SetBranchAddress("MSneu", &MSneu, &b_MSneu);
   fChain->SetBranchAddress("MN2", &MN2, &b_MN2);
   fChain->SetBranchAddress("MC1", &MC1, &b_MC1);
   fChain->SetBranchAddress("MN1", &MN1, &b_MN1);
   fChain->SetBranchAddress("NSparticleW", &NSparticleW, &b_NSparticleW);
   fChain->SetBranchAddress("NSparticleZ", &NSparticleZ, &b_NSparticleZ);
   fChain->SetBranchAddress("NSparticlePhoton", &NSparticlePhoton, &b_NSparticlePhoton);
   fChain->SetBranchAddress("LSPParents", &LSPParents, &b_LSPParents);
   fChain->SetBranchAddress("Npartons", &Npartons, &b_Npartons);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("XSec", &XSec, &b_XSec);
   fChain->SetBranchAddress("Nweight", &Nweight, &b_Nweight);
   fChain->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
   fChain->SetBranchAddress("MT2", &MT2, &b_MT2);
   fChain->SetBranchAddress("RISR_LEP", &RISR_LEP, &b_RISR_LEP);
   fChain->SetBranchAddress("PTISR_LEP", &PTISR_LEP, &b_PTISR_LEP);
   fChain->SetBranchAddress("MS_LEP", &MS_LEP, &b_MS_LEP);
   fChain->SetBranchAddress("MSV_LEP", &MSV_LEP, &b_MSV_LEP);
   fChain->SetBranchAddress("MQ_LEP", &MQ_LEP, &b_MQ_LEP);
   fChain->SetBranchAddress("gamma_LEP", &gamma_LEP, &b_gamma_LEP);
   fChain->SetBranchAddress("MPa_LEP", &MPa_LEP, &b_MPa_LEP);
   fChain->SetBranchAddress("MPb_LEP", &MPb_LEP, &b_MPb_LEP);
   fChain->SetBranchAddress("MVa_LEP", &MVa_LEP, &b_MVa_LEP);
   fChain->SetBranchAddress("MVb_LEP", &MVb_LEP, &b_MVb_LEP);
   fChain->SetBranchAddress("PTS_CM_LEP", &PTS_CM_LEP, &b_PTS_CM_LEP);
   fChain->SetBranchAddress("MS_S0_LEP", &MS_S0_LEP, &b_MS_S0_LEP);
   fChain->SetBranchAddress("MV_S0_LEP", &MV_S0_LEP, &b_MV_S0_LEP);
   fChain->SetBranchAddress("MQ_S0_LEP", &MQ_S0_LEP, &b_MQ_S0_LEP);
   fChain->SetBranchAddress("gamma_S0_LEP", &gamma_S0_LEP, &b_gamma_S0_LEP);
   fChain->SetBranchAddress("MPTilde_LEP", &MPTilde_LEP, &b_MPTilde_LEP);
   fChain->SetBranchAddress("MSTilde_LEP", &MSTilde_LEP, &b_MSTilde_LEP);
   fChain->SetBranchAddress("gammaTilde_LEP", &gammaTilde_LEP, &b_gammaTilde_LEP);
   fChain->SetBranchAddress("CosDecayAngle_Pa_LEP", &CosDecayAngle_Pa_LEP, &b_CosDecayAngle_Pa_LEP);
   fChain->SetBranchAddress("CosDecayAngle_Pb_LEP", &CosDecayAngle_Pb_LEP, &b_CosDecayAngle_Pb_LEP);
   fChain->SetBranchAddress("CosDecayAngle_S_LEP", &CosDecayAngle_S_LEP, &b_CosDecayAngle_S_LEP);
   fChain->SetBranchAddress("RZPara_LEP", &RZPara_LEP, &b_RZPara_LEP);
   fChain->SetBranchAddress("MINV_LEP", &MINV_LEP, &b_MINV_LEP);
   fChain->SetBranchAddress("MP", &MP, &b_MP);
   fChain->SetBranchAddress("Mperp_JET", &Mperp_JET, &b_Mperp_JET);
   fChain->SetBranchAddress("gammaT_JET", &gammaT_JET, &b_gammaT_JET);
   fChain->SetBranchAddress("RISR_JET", &RISR_JET, &b_RISR_JET);
   fChain->SetBranchAddress("PTISR_JET", &PTISR_JET, &b_PTISR_JET);
   fChain->SetBranchAddress("MS_JET", &MS_JET, &b_MS_JET);
   fChain->SetBranchAddress("MSV_JET", &MSV_JET, &b_MSV_JET);
   fChain->SetBranchAddress("MQ_JET", &MQ_JET, &b_MQ_JET);
   fChain->SetBranchAddress("gamma_JET", &gamma_JET, &b_gamma_JET);
   fChain->SetBranchAddress("MPa_JET", &MPa_JET, &b_MPa_JET);
   fChain->SetBranchAddress("MPb_JET", &MPb_JET, &b_MPb_JET);
   fChain->SetBranchAddress("MVa_JET", &MVa_JET, &b_MVa_JET);
   fChain->SetBranchAddress("MVb_JET", &MVb_JET, &b_MVb_JET);
   fChain->SetBranchAddress("PTS_CM_JET", &PTS_CM_JET, &b_PTS_CM_JET);
   fChain->SetBranchAddress("MS_S0_JET", &MS_S0_JET, &b_MS_S0_JET);
   fChain->SetBranchAddress("MV_S0_JET", &MV_S0_JET, &b_MV_S0_JET);
   fChain->SetBranchAddress("MQ_S0_JET", &MQ_S0_JET, &b_MQ_S0_JET);
   fChain->SetBranchAddress("gamma_S0_JET", &gamma_S0_JET, &b_gamma_S0_JET);
   fChain->SetBranchAddress("MPTilde_JET", &MPTilde_JET, &b_MPTilde_JET);
   fChain->SetBranchAddress("MSTilde_JET", &MSTilde_JET, &b_MSTilde_JET);
   fChain->SetBranchAddress("gammaTilde_JET", &gammaTilde_JET, &b_gammaTilde_JET);
   fChain->SetBranchAddress("CosDecayAngle_Pa_JET", &CosDecayAngle_Pa_JET, &b_CosDecayAngle_Pa_JET);
   fChain->SetBranchAddress("CosDecayAngle_Pb_JET", &CosDecayAngle_Pb_JET, &b_CosDecayAngle_Pb_JET);
   fChain->SetBranchAddress("CosDecayAngle_S_JET", &CosDecayAngle_S_JET, &b_CosDecayAngle_S_JET);
   fChain->SetBranchAddress("RZPara_JET", &RZPara_JET, &b_RZPara_JET);
   fChain->SetBranchAddress("MINV_JET", &MINV_JET, &b_MINV_JET);
   Notify();
}

inline Bool_t ReducedBase_V2::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void ReducedBase_V2::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t ReducedBase_V2::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
