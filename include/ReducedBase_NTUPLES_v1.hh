//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 11 21:47:02 2023 by ROOT version 6.14/09
// from TTree KUAnalysis/KUAnalysis
// found on file: root://cmseos.fnal.gov//store/group/lpcsusylep/NTUPLES_v1/Fall17_102X/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X.root
//////////////////////////////////////////////////////////

#ifndef ReducedBase_h
#define ReducedBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using std::vector;

class ReducedBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Bool_t          event_skipped;
   Double_t        weight;
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
   Bool_t          EventFlag_FailJetID;
   Bool_t          EventFlag_JetInHEM;
   Bool_t          EventFlag_JetInHEM_Pt20;
   Bool_t          EventFlag_JetInHEM_Pt20_JetID;
   Bool_t          HEM_Veto;
   Bool_t          METtrigger;
   Bool_t          METHTtrigger;
   Bool_t          METORtrigger;
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
   Int_t           Nmu;
   Int_t           Nlep;
   vector<double>  *PT_lep;
   vector<double>  *Eta_lep;
   vector<double>  *Phi_lep;
   vector<double>  *M_lep;
   vector<int>     *Charge_lep;
   vector<int>     *PDGID_lep;
   vector<double>  *RelIso_lep;
   vector<double>  *MiniIso_lep;
   vector<double>  *Dxy_lep;
   vector<double>  *DxyErr_lep;
   vector<double>  *Dz_lep;
   vector<double>  *DzErr_lep;
   vector<double>  *IP3D_lep;
   vector<double>  *SIP3D_lep;
   vector<int>     *ID_lep;
   vector<int>     *SourceID_lep;
   vector<int>     *LepQual_lep;
   vector<int>     *Index_lep;
   Int_t           Njet;
   Int_t           Nbjet;
   Int_t           NSV;
   vector<double>  *PT_SV;
   vector<double>  *Eta_SV;
   vector<double>  *Phi_SV;
   vector<double>  *M_SV;
   vector<double>  *ProbB_SV;
   vector<double>  *ProbC_SV;
   Int_t           Njet_ISR;
   Int_t           Njet_S;
   Int_t           Nbjet_ISR;
   Int_t           Nbjet_S;
   Int_t           Nlep_ISR;
   Int_t           Nlep_S;
   Int_t           NSV_ISR;
   Int_t           NSV_S;
   vector<int>     *index_SV_ISR;
   vector<int>     *index_SV_S;
   vector<int>     *index_lep_ISR;
   vector<int>     *index_lep_S;
   Int_t           Njet_a;
   Int_t           Njet_b;
   Int_t           Nbjet_a;
   Int_t           Nbjet_b;
   Int_t           Nlep_a;
   Int_t           Nlep_b;
   Int_t           NSV_a;
   Int_t           NSV_b;
   vector<int>     *index_lep_a;
   vector<int>     *index_lep_b;
   vector<int>     *index_SV_a;
   vector<int>     *index_SV_b;
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

   // List of branches
   TBranch        *b_event_skipped;   //!
   TBranch        *b_weight;   //!
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
   TBranch        *b_EventFlag_FailJetID;   //!
   TBranch        *b_EventFlag_JetInHEM;   //!
   TBranch        *b_EventFlag_JetInHEM_Pt20;   //!
   TBranch        *b_EventFlag_JetInHEM_Pt20_JetID;   //!
   TBranch        *b_HEM_Veto;   //!
   TBranch        *b_METtrigger;   //!
   TBranch        *b_METHTtrigger;   //!
   TBranch        *b_METORtrigger;   //!
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
   TBranch        *b_Nmu;   //!
   TBranch        *b_Nlep;   //!
   TBranch        *b_PT_lep;   //!
   TBranch        *b_Eta_lep;   //!
   TBranch        *b_Phi_lep;   //!
   TBranch        *b_M_lep;   //!
   TBranch        *b_Charge_lep;   //!
   TBranch        *b_PDGID_lep;   //!
   TBranch        *b_RelIso_lep;   //!
   TBranch        *b_MiniIso_lep;   //!
   TBranch        *b_Dxy_lep;   //!
   TBranch        *b_DxyErr_lep;   //!
   TBranch        *b_Dz_lep;   //!
   TBranch        *b_DzErr_lep;   //!
   TBranch        *b_IP3D_lep;   //!
   TBranch        *b_SIP3D_lep;   //!
   TBranch        *b_ID_lep;   //!
   TBranch        *b_SourceID_lep;   //!
   TBranch        *b_LepQual_lep;   //!
   TBranch        *b_Index_lep;   //!
   TBranch        *b_Njet;   //!
   TBranch        *b_Nbjet;   //!
   TBranch        *b_NSV;   //!
   TBranch        *b_PT_SV;   //!
   TBranch        *b_Eta_SV;   //!
   TBranch        *b_Phi_SV;   //!
   TBranch        *b_M_SV;   //!
   TBranch        *b_ProbB_SV;   //!
   TBranch        *b_ProbC_SV;   //!
   TBranch        *b_Njet_ISR;   //!
   TBranch        *b_Njet_S;   //!
   TBranch        *b_Nbjet_ISR;   //!
   TBranch        *b_Nbjet_S;   //!
   TBranch        *b_Nlep_ISR;   //!
   TBranch        *b_Nlep_S;   //!
   TBranch        *b_NSV_ISR;   //!
   TBranch        *b_NSV_S;   //!
   TBranch        *b_index_SV_ISR;   //!
   TBranch        *b_index_SV_S;   //!
   TBranch        *b_index_lep_ISR;   //!
   TBranch        *b_index_lep_S;   //!
   TBranch        *b_Njet_a;   //!
   TBranch        *b_Njet_b;   //!
   TBranch        *b_Nbjet_a;   //!
   TBranch        *b_Nbjet_b;   //!
   TBranch        *b_Nlep_a;   //!
   TBranch        *b_Nlep_b;   //!
   TBranch        *b_NSV_a;   //!
   TBranch        *b_NSV_b;   //!
   TBranch        *b_index_lep_a;   //!
   TBranch        *b_index_lep_b;   //!
   TBranch        *b_index_SV_a;   //!
   TBranch        *b_index_SV_b;   //!
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

   ReducedBase(TTree *tree=0);
   virtual ~ReducedBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

inline ReducedBase::ReducedBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://cmseos.fnal.gov//store/group/lpcsusylep/NTUPLES_v2/Fall17_102X/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://cmseos.fnal.gov//store/group/lpcsusylep/NTUPLES_v2/Fall17_102X/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8_Fall17_102X.root");
      }
      f->GetObject("KUAnalysis",tree);

   }
   Init(tree);
}

inline ReducedBase::~ReducedBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t ReducedBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t ReducedBase::LoadTree(Long64_t entry)
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

inline void ReducedBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PT_lep = 0;
   Eta_lep = 0;
   Phi_lep = 0;
   M_lep = 0;
   Charge_lep = 0;
   PDGID_lep = 0;
   RelIso_lep = 0;
   MiniIso_lep = 0;
   Dxy_lep = 0;
   DxyErr_lep = 0;
   Dz_lep = 0;
   DzErr_lep = 0;
   IP3D_lep = 0;
   SIP3D_lep = 0;
   ID_lep = 0;
   SourceID_lep = 0;
   LepQual_lep = 0;
   Index_lep = 0;
   PT_SV = 0;
   Eta_SV = 0;
   Phi_SV = 0;
   M_SV = 0;
   ProbB_SV = 0;
   ProbC_SV = 0;
   index_SV_ISR = 0;
   index_SV_S = 0;
   index_lep_ISR = 0;
   index_lep_S = 0;
   index_lep_a = 0;
   index_lep_b = 0;
   index_SV_a = 0;
   index_SV_b = 0;
   genPT_lep = 0;
   genEta_lep = 0;
   genPhi_lep = 0;
   genM_lep = 0;
   genCharge_lep = 0;
   genPDGID_lep = 0;
   genMomPDGID_lep = 0;
   genSourceID_lep = 0;
   genIndex_lep = 0;
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
   fChain->SetBranchAddress("weight", &weight, &b_weight);
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
   fChain->SetBranchAddress("EventFlag_FailJetID", &EventFlag_FailJetID, &b_EventFlag_FailJetID);
   fChain->SetBranchAddress("EventFlag_JetInHEM", &EventFlag_JetInHEM, &b_EventFlag_JetInHEM);
   fChain->SetBranchAddress("EventFlag_JetInHEM_Pt20", &EventFlag_JetInHEM_Pt20, &b_EventFlag_JetInHEM_Pt20);
   fChain->SetBranchAddress("EventFlag_JetInHEM_Pt20_JetID", &EventFlag_JetInHEM_Pt20_JetID, &b_EventFlag_JetInHEM_Pt20_JetID);
   fChain->SetBranchAddress("HEM_Veto", &HEM_Veto, &b_HEM_Veto);
   fChain->SetBranchAddress("METtrigger", &METtrigger, &b_METtrigger);
   fChain->SetBranchAddress("METHTtrigger", &METHTtrigger, &b_METHTtrigger);
   fChain->SetBranchAddress("METORtrigger", &METORtrigger, &b_METORtrigger);
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
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("Nlep", &Nlep, &b_Nlep);
   fChain->SetBranchAddress("PT_lep", &PT_lep, &b_PT_lep);
   fChain->SetBranchAddress("Eta_lep", &Eta_lep, &b_Eta_lep);
   fChain->SetBranchAddress("Phi_lep", &Phi_lep, &b_Phi_lep);
   fChain->SetBranchAddress("M_lep", &M_lep, &b_M_lep);
   fChain->SetBranchAddress("Charge_lep", &Charge_lep, &b_Charge_lep);
   fChain->SetBranchAddress("PDGID_lep", &PDGID_lep, &b_PDGID_lep);
   fChain->SetBranchAddress("RelIso_lep", &RelIso_lep, &b_RelIso_lep);
   fChain->SetBranchAddress("MiniIso_lep", &MiniIso_lep, &b_MiniIso_lep);
   fChain->SetBranchAddress("Dxy_lep", &Dxy_lep, &b_Dxy_lep);
   fChain->SetBranchAddress("DxyErr_lep", &DxyErr_lep, &b_DxyErr_lep);
   fChain->SetBranchAddress("Dz_lep", &Dz_lep, &b_Dz_lep);
   fChain->SetBranchAddress("DzErr_lep", &DzErr_lep, &b_DzErr_lep);
   fChain->SetBranchAddress("IP3D_lep", &IP3D_lep, &b_IP3D_lep);
   fChain->SetBranchAddress("SIP3D_lep", &SIP3D_lep, &b_SIP3D_lep);
   fChain->SetBranchAddress("ID_lep", &ID_lep, &b_ID_lep);
   fChain->SetBranchAddress("SourceID_lep", &SourceID_lep, &b_SourceID_lep);
   fChain->SetBranchAddress("LepQual_lep", &LepQual_lep, &b_LepQual_lep);
   fChain->SetBranchAddress("Index_lep", &Index_lep, &b_Index_lep);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("Nbjet", &Nbjet, &b_Nbjet);
   fChain->SetBranchAddress("NSV", &NSV, &b_NSV);
   fChain->SetBranchAddress("PT_SV", &PT_SV, &b_PT_SV);
   fChain->SetBranchAddress("Eta_SV", &Eta_SV, &b_Eta_SV);
   fChain->SetBranchAddress("Phi_SV", &Phi_SV, &b_Phi_SV);
   fChain->SetBranchAddress("M_SV", &M_SV, &b_M_SV);
   fChain->SetBranchAddress("ProbB_SV", &ProbB_SV, &b_ProbB_SV);
   fChain->SetBranchAddress("ProbC_SV", &ProbC_SV, &b_ProbC_SV);
   fChain->SetBranchAddress("Njet_ISR", &Njet_ISR, &b_Njet_ISR);
   fChain->SetBranchAddress("Njet_S", &Njet_S, &b_Njet_S);
   fChain->SetBranchAddress("Nbjet_ISR", &Nbjet_ISR, &b_Nbjet_ISR);
   fChain->SetBranchAddress("Nbjet_S", &Nbjet_S, &b_Nbjet_S);
   fChain->SetBranchAddress("Nlep_ISR", &Nlep_ISR, &b_Nlep_ISR);
   fChain->SetBranchAddress("Nlep_S", &Nlep_S, &b_Nlep_S);
   fChain->SetBranchAddress("NSV_ISR", &NSV_ISR, &b_NSV_ISR);
   fChain->SetBranchAddress("NSV_S", &NSV_S, &b_NSV_S);
   fChain->SetBranchAddress("index_SV_ISR", &index_SV_ISR, &b_index_SV_ISR);
   fChain->SetBranchAddress("index_SV_S", &index_SV_S, &b_index_SV_S);
   fChain->SetBranchAddress("index_lep_ISR", &index_lep_ISR, &b_index_lep_ISR);
   fChain->SetBranchAddress("index_lep_S", &index_lep_S, &b_index_lep_S);
   fChain->SetBranchAddress("Njet_a", &Njet_a, &b_Njet_a);
   fChain->SetBranchAddress("Njet_b", &Njet_b, &b_Njet_b);
   fChain->SetBranchAddress("Nbjet_a", &Nbjet_a, &b_Nbjet_a);
   fChain->SetBranchAddress("Nbjet_b", &Nbjet_b, &b_Nbjet_b);
   fChain->SetBranchAddress("Nlep_a", &Nlep_a, &b_Nlep_a);
   fChain->SetBranchAddress("Nlep_b", &Nlep_b, &b_Nlep_b);
   fChain->SetBranchAddress("NSV_a", &NSV_a, &b_NSV_a);
   fChain->SetBranchAddress("NSV_b", &NSV_b, &b_NSV_b);
   fChain->SetBranchAddress("index_lep_a", &index_lep_a, &b_index_lep_a);
   fChain->SetBranchAddress("index_lep_b", &index_lep_b, &b_index_lep_b);
   fChain->SetBranchAddress("index_SV_a", &index_SV_a, &b_index_SV_a);
   fChain->SetBranchAddress("index_SV_b", &index_SV_b, &b_index_SV_b);
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
   Notify();
}

inline Bool_t ReducedBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void ReducedBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t ReducedBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif
