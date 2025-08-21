#include "ReducedNtuple.hh"
#include "ParticleList.hh"
#include "TInterpreter.h"

#include "SUSYNANOBase.hh"
#include "NANOULBase.hh"
#include "NANORun3.hh"
#include "Leptonic.hh"

using namespace RestFrames;

template <class Base>
ReducedNtuple<Base>::ReducedNtuple(TTree* tree)
  : NtupleBase<Base>(tree)
{
  m_library_generated = false;

  for(int t = 0; t < m_aTrees; t++){
    // first tree is ISR Boosted tree
    LAB[t] = new LabRecoFrame("LAB","LAB");
    CM[t]  = new DecayRecoFrame("CM","CM");
    S[t]   = new DecayRecoFrame("S","#tilde{S}");
    X2a[t] = new DecayRecoFrame("X2a","#tilde{#chi}_{2a}");
    X2b[t] = new DecayRecoFrame("X2b","#tilde{#chi}_{2b}");
    ISR[t] = new VisibleRecoFrame("ISR","ISR");
    sJa[t] = new SelfAssemblingRecoFrame("sJa","jets_{a}");
    sJb[t] = new SelfAssemblingRecoFrame("sJb","jets_{b}");
    sLa[t] = new SelfAssemblingRecoFrame("sLa","leps_{a}");
    sLb[t] = new SelfAssemblingRecoFrame("sLb","leps_{b}");
    Ja[t] = new VisibleRecoFrame("Ja","jets_{a}");
    Jb[t] = new VisibleRecoFrame("Jb","jets_{b}");
    La[t] = new VisibleRecoFrame("La","leps_{a}");
    Lb[t] = new VisibleRecoFrame("Lb","leps_{b}");
    X1a[t] = new InvisibleRecoFrame("X1a","#tilde{#chi}_{1a}");
    X1b[t] = new InvisibleRecoFrame("X1b","#tilde{#chi}_{1b}");
    
    if(t==0){ // ISR
      LAB[t]->SetChildFrame(*CM[t]);
      CM[t]->AddChildFrame(*S[t]);
      CM[t]->AddChildFrame(*ISR[t]);
      S[t]->AddChildFrame(*X2a[t]);
      S[t]->AddChildFrame(*X2b[t]);
      X2a[t]->AddChildFrame(*sJa[t]);
      X2b[t]->AddChildFrame(*sJb[t]);
      X2a[t]->AddChildFrame(*sLa[t]);
      X2b[t]->AddChildFrame(*sLb[t]);
      sJa[t]->AddChildFrame(*Ja[t]);
      sJb[t]->AddChildFrame(*Jb[t]);
      sLa[t]->AddChildFrame(*La[t]);
      sLb[t]->AddChildFrame(*Lb[t]);
      X2a[t]->AddChildFrame(*X1a[t]);
      X2b[t]->AddChildFrame(*X1b[t]);
    }

    if(t==1){ // LEP only
      LAB[t]->SetChildFrame(*CM[t]);
      CM[t]->AddChildFrame(*S[t]);
      CM[t]->AddChildFrame(*ISR[t]);
      S[t]->AddChildFrame(*X2a[t]);
      S[t]->AddChildFrame(*X2b[t]);
      X2a[t]->AddChildFrame(*sLa[t]);
      X2b[t]->AddChildFrame(*sLb[t]);
      sLa[t]->AddChildFrame(*La[t]);
      sLb[t]->AddChildFrame(*Lb[t]);
      X2a[t]->AddChildFrame(*X1a[t]);
      X2b[t]->AddChildFrame(*X1b[t]);
    }

    if(t==2){ // JET & ISR only
      LAB[t]->SetChildFrame(*CM[t]);
      CM[t]->AddChildFrame(*S[t]);
      CM[t]->AddChildFrame(*ISR[t]);
      S[t]->AddChildFrame(*X2a[t]);
      S[t]->AddChildFrame(*X2b[t]);
      X2a[t]->AddChildFrame(*sJa[t]);
      X2b[t]->AddChildFrame(*sJb[t]);
      sJa[t]->AddChildFrame(*Ja[t]);
      sJb[t]->AddChildFrame(*Jb[t]);
      X2a[t]->AddChildFrame(*X1a[t]);
      X2b[t]->AddChildFrame(*X1b[t]);
    }

    if(t==3){ // JET only
      LAB[t]->SetChildFrame(*S[t]);
      S[t]->AddChildFrame(*X2a[t]);
      S[t]->AddChildFrame(*X2b[t]);
      X2a[t]->AddChildFrame(*sJa[t]);
      X2b[t]->AddChildFrame(*sJb[t]);
      sJa[t]->AddChildFrame(*Ja[t]);
      sJb[t]->AddChildFrame(*Jb[t]);
      X2a[t]->AddChildFrame(*X1a[t]);
      X2b[t]->AddChildFrame(*X1b[t]);
    }

    if(!LAB[t]->InitializeTree()){
      cout << "Problem initializing tree" << endl;
    }

    INV[t] = new InvisibleGroup("INV","Invisible System");
    INV[t]->AddFrame(*X1a[t]);
    INV[t]->AddFrame(*X1b[t]);
    
    InvM[t] = new SetMassInvJigsaw("InvM", "Set inv. system mass");
    INV[t]->AddJigsaw(*InvM[t]);
    
    InvEta[t] = new SetRapidityInvJigsaw("InvEta", "Set inv. system rapidity");
    INV[t]->AddJigsaw(*InvEta[t]);
    InvEta[t]->AddVisibleFrames(S[t]->GetListVisibleFrames());
    
    InvSplit[t] = new MinMassesSqInvJigsaw("InvSplit", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}", 2);
    INV[t]->AddJigsaw(*InvSplit[t]);
    InvSplit[t]->AddVisibleFrames(X2a[t]->GetListVisibleFrames(),0);
    InvSplit[t]->AddVisibleFrames(X2b[t]->GetListVisibleFrames(),1);
    InvSplit[t]->AddInvisibleFrame(*X1a[t], 0);
    InvSplit[t]->AddInvisibleFrame(*X1b[t], 1);
    
    COMB_J[t] =  new CombinatoricGroup("COMB_J", "Combinatoric System of Jets");
    CombSplit_ISR[t] = new MinMassesSqCombJigsaw("CombSplitSq_ISR", "Minimize M_{ISR}_{T}^{2} + M_{S}_{T}^{2}",2,2);
    CombSplit_J[t] = new MinMassesSqCombJigsaw("CombSplitSq_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
    
    COMB_L[t] =  new CombinatoricGroup("COMB_L", "Combinatoric System of Leps");
    CombSplit_L[t] = new MinMassesSqCombJigsaw("CombSplit_L", "Minimize M_{Va}^{2} + M_{Vb}^{2}",2,2);
    
    if(t==0 || t==2){
      COMB_J[t]->AddFrame(*ISR[t]);
      COMB_J[t]->SetNElementsForFrame(*ISR[t], 1);
      COMB_J[t]->AddFrame(*Ja[t]);
      COMB_J[t]->SetNElementsForFrame(*Ja[t], 0);
      COMB_J[t]->AddFrame(*Jb[t]);
      COMB_J[t]->SetNElementsForFrame(*Jb[t], 0);
    }
    
    if(t==0 || t==2){
      COMB_J[t]->AddJigsaw(*CombSplit_ISR[t]);
      CombSplit_ISR[t]->SetTransverse();
      CombSplit_ISR[t]->AddCombFrame(*ISR[t], 0);
      CombSplit_ISR[t]->AddCombFrame(*Ja[t], 1);
      CombSplit_ISR[t]->AddCombFrame(*Jb[t], 1);
      CombSplit_ISR[t]->AddObjectFrames(ISR[t]->GetListVisibleFrames(), 0);
      CombSplit_ISR[t]->AddObjectFrames(S[t]->GetListVisibleFrames(), 1);
      
      COMB_J[t]->AddJigsaw(*CombSplit_J[t]);
      CombSplit_J[t]->AddCombFrame(*Ja[t], 0);
      CombSplit_J[t]->AddCombFrame(*Jb[t], 1);
      CombSplit_J[t]->AddObjectFrames(X2a[t]->GetListVisibleFrames(), 0);
      CombSplit_J[t]->AddObjectFrames(X2b[t]->GetListVisibleFrames(), 1);
    }

    if(t == 3){
      COMB_J[t]->AddFrame(*Ja[t]);
      COMB_J[t]->SetNElementsForFrame(*Ja[t], 0);
      COMB_J[t]->AddFrame(*Jb[t]);
      COMB_J[t]->SetNElementsForFrame(*Jb[t], 0);

      COMB_J[t]->AddJigsaw(*CombSplit_J[t]);
      CombSplit_J[t]->AddCombFrame(*Ja[t], 0);
      CombSplit_J[t]->AddCombFrame(*Jb[t], 1);
      CombSplit_J[t]->AddObjectFrames(X2a[t]->GetListVisibleFrames(), 0);
      CombSplit_J[t]->AddObjectFrames(X2b[t]->GetListVisibleFrames(), 1);
    }
      
    if(t < 2){
      COMB_L[t]->AddFrame(*La[t]);
      COMB_L[t]->SetNElementsForFrame(*La[t], 0);
      COMB_L[t]->AddFrame(*Lb[t]);
      COMB_L[t]->SetNElementsForFrame(*Lb[t], 0);
          
      COMB_L[t]->AddJigsaw(*CombSplit_L[t]);
      CombSplit_L[t]->AddCombFrame(*La[t], 0);
      CombSplit_L[t]->AddCombFrame(*Lb[t], 1);
    }

    if(t==0){
      CombSplit_L[t]->AddObjectFrames(X2a[t]->GetListVisibleFrames(), 0);
      CombSplit_L[t]->AddObjectFrames(X2b[t]->GetListVisibleFrames(), 1);
    }
    if(t==1){
      CombSplit_L[t]->AddObjectFrame(*La[t], 0);
      CombSplit_L[t]->AddObjectFrame(*Lb[t], 1);
    }

    if(!LAB[t]->InitializeAnalysis())
      cout << "Problem initializing analysis tree" << endl;
    /*
      TreePlot tree_plot("TreePlot","TreePlot");

      tree_plot.SetTree(*LAB);
      tree_plot.Draw("ANA_tree", "Reconstruction Tree");

      tree_plot.SetTree(*COMB_J);
      tree_plot.Draw("ANA_comb", "Combinatoric Jigsaws for jets");

      tree_plot.SetTree(*COMB_L);
      tree_plot.Draw("ANA_comb_L", "Combinatoric Jigsaws for leps");

      tree_plot.SetTree(*INV);
      tree_plot.Draw("ANA_inv", "Invisible Jigsaws");
 
      tree_plot.WriteOutput("trees.root");
    */
    m_treeSkipped.push_back(false);
  }
}

template <class Base>
ReducedNtuple<Base>::~ReducedNtuple() {
 
  for(int t = 0; t < m_aTrees; t++){
    delete LAB[t];
    delete CM[t];
    delete S[t];
    delete X2a[t];
    delete X2b[t];
    delete ISR[t];
    delete sJa[t];
    delete sJb[t];
    delete sLa[t];
    delete sLb[t];
    delete Ja[t];
    delete Jb[t];
    delete La[t];
    delete Lb[t];
    delete X1a[t];
    delete X1b[t];

    delete INV[t];
    delete InvM[t];
    delete InvEta[t];
    delete InvSplit[t];
      
    delete COMB_J[t];
    delete CombSplit_ISR[t];
    delete CombSplit_J[t];
      
    delete COMB_L[t];
    delete CombSplit_L[t];
  }
  
}

template <class Base>
TTree* ReducedNtuple<Base>::InitOutputTree(const string& sample, bool do_slim){
  
  TTree* tree = (TTree*) new TTree(sample.c_str(), sample.c_str());

  tree->Branch("event_skipped", &m_event_skipped);
  tree->Branch("treeSkipped", &m_treeSkipped);
  
  tree->Branch("weight", &m_weight);
  tree->Branch("weight2", &m_weight2);

  tree->Branch("PUweight", &m_PUweight);
  tree->Branch("PUweight_up", &m_PUweight_up);
  tree->Branch("PUweight_down", &m_PUweight_down);

  tree->Branch("MuFweight", &m_MuFweight);
  tree->Branch("MuFweight_up", &m_MuFweight_up);
  tree->Branch("MuFweight_down", &m_MuFweight_down);

  tree->Branch("MuRweight", &m_MuRweight);
  tree->Branch("MuRweight_up", &m_MuRweight_up);
  tree->Branch("MuRweight_down", &m_MuRweight_down);

  tree->Branch("PDFweight", &m_PDFweight);
  tree->Branch("PDFweight_up", &m_PDFweight_up);
  tree->Branch("PDFweight_down", &m_PDFweight_down);

  tree->Branch("BtagHFSFweight", &m_BtagHFSFweight);
  tree->Branch("BtagHFSFweight_up", &m_BtagHFSFweight_up);
  tree->Branch("BtagHFSFweight_down", &m_BtagHFSFweight_down);
  tree->Branch("BtagLFSFweight", &m_BtagLFSFweight);
  tree->Branch("BtagLFSFweight_up", &m_BtagLFSFweight_up);
  tree->Branch("BtagLFSFweight_down", &m_BtagLFSFweight_down);

  tree->Branch("elIDSFweight", &m_elIDSFweight);
  tree->Branch("elIDSFweight_up", &m_elIDSFweight_up);
  tree->Branch("elIDSFweight_down", &m_elIDSFweight_down);
  tree->Branch("elISOSFweight", &m_elISOSFweight);
  tree->Branch("elISOSFweight_up", &m_elISOSFweight_up);
  tree->Branch("elISOSFweight_down", &m_elISOSFweight_down);
  tree->Branch("elSIPSFweight", &m_elSIPSFweight);
  tree->Branch("elSIPSFweight_up", &m_elSIPSFweight_up);
  tree->Branch("elSIPSFweight_down", &m_elSIPSFweight_down);
  tree->Branch("elVLSFweight", &m_elVLSFweight);
  tree->Branch("elVLSFweight_up", &m_elVLSFweight_up);
  tree->Branch("elVLSFweight_down", &m_elVLSFweight_down);
  tree->Branch("muIDSFweight", &m_muIDSFweight);
  tree->Branch("muIDSFweight_up", &m_muIDSFweight_up);
  tree->Branch("muIDSFweight_down", &m_muIDSFweight_down);
  tree->Branch("muISOSFweight", &m_muISOSFweight);
  tree->Branch("muISOSFweight_up", &m_muISOSFweight_up);
  tree->Branch("muISOSFweight_down", &m_muISOSFweight_down);
  tree->Branch("muSIPSFweight", &m_muSIPSFweight);
  tree->Branch("muSIPSFweight_up", &m_muSIPSFweight_up);
  tree->Branch("muSIPSFweight_down", &m_muSIPSFweight_down);
  tree->Branch("muVLSFweight", &m_muVLSFweight);
  tree->Branch("muVLSFweight_up", &m_muVLSFweight_up);
  tree->Branch("muVLSFweight_down", &m_muVLSFweight_down);
   
  tree->Branch("MetTrigSFweight", &m_MetTrigSFweight);
  tree->Branch("MetTrigSFweight_up", &m_MetTrigSFweight_up);
  tree->Branch("MetTrigSFweight_down", &m_MetTrigSFweight_down);
  tree->Branch("MetTrigSFCurveIndex", &m_MetTrigSFCurveIndex);

  tree->Branch("runnum", &m_runnum);
  tree->Branch("luminum", &m_luminum);
  tree->Branch("eventnum", &m_eventnum);

  tree->Branch("NPV", &m_NPV);

  tree->Branch("EventFilter", &m_EventFilter);
  tree->Branch("FastSimEventVeto", &m_FastSimEventVeto);

  tree->Branch("PrefireWeight", &m_PrefireWeight);
  tree->Branch("PrefireWeight_up", &m_PrefireWeight_up);
  tree->Branch("PrefireWeight_down", &m_PrefireWeight_down);

  if(!do_slim){
    tree->Branch("EventFlag_FailJetID", &m_EventFlag_FailJetID);
    tree->Branch("EventFlag_JetInHEM", &m_EventFlag_JetInHEM);
    tree->Branch("EventFlag_JetInHEM_Pt20", &m_EventFlag_JetInHEM_Pt20);
    tree->Branch("EventFlag_JetInHEM_Pt20_JetID", &m_EventFlag_JetInHEM_Pt20_JetID);
  }
  tree->Branch("HEM_Veto", &m_HEM_Veto);

  tree->Branch("METtrigger", &m_METtrigger);
  tree->Branch("METORtrigger", &m_METORtrigger);
  tree->Branch("SingleElectrontrigger", &m_SingleElectrontrigger);
  tree->Branch("SingleMuontrigger", &m_SingleMuontrigger);
  tree->Branch("EMutrigger", &m_EMutrigger);
  tree->Branch("EMuMutrigger", &m_EMuMutrigger);
  tree->Branch("EMuEtrigger", &m_EMuEtrigger);
  tree->Branch("DoubleElectrontrigger", &m_DoubleElectrontrigger);
  tree->Branch("DoubleMuontrigger", &m_DoubleMuontrigger);
  tree->Branch("TripleElectrontrigger", &m_TripleElectrontrigger);
  tree->Branch("TripleMuonLowPTtrigger", &m_TripleMuonLowPTtrigger);
  tree->Branch("TripleMuonHighPTtrigger", &m_TripleMuonHighPTtrigger);
  tree->Branch("DiMuEleLowPTtrigger", &m_DiMuEleLowPTtrigger);
  tree->Branch("DiMuEleHighPTtrigger", &m_DiMuEleHighPTtrigger);
  tree->Branch("DiEleMutrigger", &m_DiEleMutrigger);
  
  tree->Branch("MET", &m_MET);
  tree->Branch("MET_phi", &m_MET_phi);

  tree->Branch("altMET", &m_altMET);
  tree->Branch("altMET_phi", &m_altMET_phi);

  tree->Branch("LHE_HT", &m_LHE_HT);
  tree->Branch("LHE_HTIncoming", &m_LHE_HTIncoming);

  tree->Branch("HT_eta24", &m_HT_eta24);
  tree->Branch("HT_eta24_id", &m_HT_eta24_id);
  tree->Branch("HT_eta3", &m_HT_eta3);
  tree->Branch("HT_eta3_id", &m_HT_eta3_id);
  tree->Branch("HT_eta5", &m_HT_eta5);
  tree->Branch("HT_eta5_id", &m_HT_eta5_id);

  tree->Branch("Nele", &m_Nele);
  tree->Branch("Nlele", &m_Nlele);
  tree->Branch("Nmu", &m_Nmu);

  tree->Branch("Nlep", &m_Nlep);
  tree->Branch("PT_lep",  &m_PT_lep);
  
  tree->Branch("Eta_lep", &m_Eta_lep);
  tree->Branch("Phi_lep", &m_Phi_lep);
  tree->Branch("M_lep",   &m_M_lep);
  
  tree->Branch("Charge_lep",  &m_Charge_lep);
  tree->Branch("PDGID_lep",   &m_PDGID_lep);
  tree->Branch("ID_lep", &m_ID_lep);
  tree->Branch("SourceID_lep", &m_SourceID_lep);
  tree->Branch("LepQual_lep", &m_LepQual_lep);
  tree->Branch("IsLowPt_lep", &m_IsLowPt_lep);
  tree->Branch("TightCharge_lep", &m_TightCharge_lep);
  tree->Branch("Index_lep", &m_Index_lep);

  tree->Branch("RelIso_lep",  &m_RelIso_lep);
  tree->Branch("MiniIso_lep", &m_MiniIso_lep);
  tree->Branch("Dxy_lep", &m_Dxy_lep);
  tree->Branch("DxyErr_lep", &m_DxyErr_lep);
  tree->Branch("Dz_lep", &m_Dz_lep);
  tree->Branch("DzErr_lep", &m_DzErr_lep);
  tree->Branch("IP3D_lep", &m_IP3D_lep);
  tree->Branch("SIP3D_lep", &m_SIP3D_lep);
  
  tree->Branch("Njet", &m_Njet);
  tree->Branch("Nbjet", &m_Nbjet);

  if(!do_slim){
    tree->Branch("PT_jet",  &m_PT_jet);
    tree->Branch("Eta_jet", &m_Eta_jet);
    tree->Branch("Phi_jet", &m_Phi_jet);
    tree->Branch("M_jet",   &m_M_jet);
    tree->Branch("Btag_jet",   &m_Btag_jet);
    tree->Branch("BtagID_jet",   &m_BtagID_jet);
    tree->Branch("Flavor_jet",   &m_Flavor_jet);
    tree->Branch("index_jet_a", &m_index_jet_a);
    tree->Branch("index_jet_b", &m_index_jet_b);
    tree->Branch("index_jet_ISR", &m_index_jet_ISR);
    tree->Branch("index_jet_S", &m_index_jet_S);
  }

  if(!do_slim){
    tree->Branch("PT_Genjet",  &m_PT_Genjet);
    tree->Branch("Eta_Genjet", &m_Eta_Genjet);
    tree->Branch("Phi_Genjet", &m_Phi_Genjet);
    tree->Branch("M_Genjet",   &m_M_Genjet);
    tree->Branch("Index_jet",   &m_Index_jet);
  }

  if(!do_slim){
    tree->Branch("dphi_lep_S", &m_dphi_lep_S);
    tree->Branch("cos_lep_S", &m_cos_lep_S);
    tree->Branch("dphi_jet_S", &m_dphi_jet_S);
    tree->Branch("cos_jet_S", &m_cos_jet_S);

    tree->Branch("dphiMET_lep_S", &m_dphiMET_lep_S);
    tree->Branch("dphiMET_jet_S", &m_dphiMET_jet_S);
  }
  
  // Object counting vars
  // Base ISR+JET+LEP tree
  tree->Branch("Njet_ISR", &m_Njet_ISR);
  tree->Branch("Njet_S", &m_Njet_S);
  tree->Branch("Nbjet_ISR", &m_Nbjet_ISR);
  tree->Branch("Nbjet_S", &m_Nbjet_S);
  tree->Branch("Nlep_ISR", &m_Nlep_ISR);
  tree->Branch("Nlep_S", &m_Nlep_S);
  tree->Branch("index_lep_ISR", &m_index_lep_ISR);
  tree->Branch("index_lep_S", &m_index_lep_S);
  tree->Branch("Njet_a", &m_Njet_a);
  tree->Branch("Njet_b", &m_Njet_b);
  tree->Branch("Nbjet_a", &m_Nbjet_a);
  tree->Branch("Nbjet_b", &m_Nbjet_b);
  tree->Branch("Nlep_a", &m_Nlep_a);
  tree->Branch("Nlep_b", &m_Nlep_b);
  tree->Branch("index_lep_a", &m_index_lep_a);
  tree->Branch("index_lep_b", &m_index_lep_b);

  // ISR+LEP tree
  tree->Branch("Nlep_a_LEP", &m_Nlep_a_LEP);
  tree->Branch("Nlep_b_LEP", &m_Nlep_b_LEP);
  tree->Branch("index_lep_a_LEP", &m_index_lep_a_LEP);
  tree->Branch("index_lep_b_LEP", &m_index_lep_b_LEP);

  // ISR+JET tree
  tree->Branch("Njet_ISR_JET_ISR", &m_Njet_ISR_JET_ISR);
  tree->Branch("Njet_S_JET_ISR", &m_Njet_S_JET_ISR);
  tree->Branch("Nbjet_ISR_JET_ISR", &m_Nbjet_ISR_JET_ISR);
  tree->Branch("Nbjet_S_JET_ISR", &m_Nbjet_S_JET_ISR);
  tree->Branch("Njet_a_JET_ISR", &m_Njet_a_JET_ISR);
  tree->Branch("Njet_b_JET_ISR", &m_Njet_b_JET_ISR);
  tree->Branch("Nbjet_a_JET_ISR", &m_Nbjet_a_JET_ISR);
  tree->Branch("Nbjet_b_JET_ISR", &m_Nbjet_b_JET_ISR);
  tree->Branch("index_jet_a_JET_ISR", &m_index_jet_a_JET_ISR);
  tree->Branch("index_jet_b_JET_ISR", &m_index_jet_b_JET_ISR);
  tree->Branch("index_jet_ISR_JET_ISR", &m_index_jet_ISR_JET_ISR);
  tree->Branch("index_jet_S_JET_ISR", &m_index_jet_S_JET_ISR);

  // JET tree
  tree->Branch("Njet_a_JET", &m_Njet_a_JET);
  tree->Branch("Njet_b_JET", &m_Njet_b_JET);
  tree->Branch("Nbjet_a_JET", &m_Nbjet_a_JET);
  tree->Branch("Nbjet_b_JET", &m_Nbjet_b_JET);
  tree->Branch("index_jet_a_JET", &m_index_jet_a_JET);
  tree->Branch("index_jet_b_JET", &m_index_jet_b_JET);

  // Kinematics
  tree->Branch("PTCM", &m_PTCM);
  tree->Branch("PzCM", &m_PzCM);
  tree->Branch("cosCM", &m_cosCM);
  tree->Branch("dphiCM", &m_dphiCM);
  tree->Branch("dphiCMI", &m_dphiCMI);
  tree->Branch("dphiMET_V", &m_dphiMET_V);

  tree->Branch("CosDecayAngle_Pa", &m_CosDecayAngle_Pa);
  tree->Branch("CosDecayAngle_Pb", &m_CosDecayAngle_Pb);
  tree->Branch("CosDecayAngle_Va", &m_CosDecayAngle_Va);
  tree->Branch("CosDecayAngle_Vb", &m_CosDecayAngle_Vb);
  tree->Branch("CosDecayAngle_S", &m_CosDecayAngle_S);

  tree->Branch("PTCM_LEP", &m_PTCM_LEP);
  tree->Branch("PzCM_LEP", &m_PzCM_LEP);
  tree->Branch("cosCM_LEP", &m_cosCM_LEP);
  tree->Branch("dphiCM_LEP", &m_dphiCM_LEP);
  tree->Branch("dphiCMI_LEP", &m_dphiCMI_LEP);
  tree->Branch("dphiMET_V_LEP", &m_dphiMET_V_LEP);

  tree->Branch("PTCM_JET_ISR", &m_PTCM_JET_ISR);
  tree->Branch("PzCM_JET_ISR", &m_PzCM_JET_ISR);
  tree->Branch("cosCM_JET_ISR", &m_cosCM_JET_ISR);
  tree->Branch("dphiCM_JET_ISR", &m_dphiCM_JET_ISR);
  tree->Branch("dphiCMI_JET_ISR", &m_dphiCMI_JET_ISR);
  tree->Branch("dphiMET_V_JET_ISR", &m_dphiMET_V_JET_ISR);

  tree->Branch("PTCM_JET", &m_PTCM_JET);
  tree->Branch("PzCM_JET", &m_PzCM_JET);
  tree->Branch("cosCM_JET", &m_cosCM_JET);
  tree->Branch("dphiCM_JET", &m_dphiCM_JET);
  tree->Branch("dphiCMI_JET", &m_dphiCMI_JET);
  tree->Branch("dphiMET_V_JET", &m_dphiMET_V_JET);

  tree->Branch("Mperp", &m_Mperp);
  tree->Branch("gammaT", &m_gammaT);
  tree->Branch("EJ_BoostT", &m_EJ_BoostT);
  tree->Branch("EL_BoostT", &m_EL_BoostT);
  tree->Branch("PTISR", &m_PTISR);
  tree->Branch("RISR", &m_RISR);

  tree->Branch("EtaCM", &m_EtaCM);
  tree->Branch("PhiCM", &m_PhiCM);
  tree->Branch("MCM", &m_MCM);
  tree->Branch("EtaS", &m_EtaS);
  tree->Branch("PhiS", &m_PhiS);

  tree->Branch("LAB_Pt", &m_LAB_Pt);
  tree->Branch("LAB_Eta", &m_LAB_Eta);
  tree->Branch("LAB_Phi", &m_LAB_Phi);
  tree->Branch("LAB_M", &m_LAB_M);

  // New Observables for Run3/Cascades
  // analysis using ISR boosted tree
  tree->Branch("MSperpCM0", &m_MSperpCM0); 
  tree->Branch("MaPerpCM0", &m_MaPerpCM0);
  tree->Branch("MbPerpCM0", &m_MbPerpCM0);
  tree->Branch("MaVPerpCM0", &m_MaVPerpCM0);
  tree->Branch("MbVPerpCM0", &m_MbVPerpCM0);
  tree->Branch("MQperpCM0", &m_MQperpCM0); // sqrt(Ma*Ma + Mb*Mb)/sqrt(2)
  tree->Branch("gammaPerpCM0", &m_gammaPerpCM0); // 2*MQ/MS
  tree->Branch("MSCM0", &m_MSCM0); 
  tree->Branch("MaCM0", &m_MaCM0);
  tree->Branch("MbCM0", &m_MbCM0);
  tree->Branch("MaVCM0", &m_MaVCM0);
  tree->Branch("MbVCM0", &m_MbVCM0);
  tree->Branch("MQCM0", &m_MQCM0); // sqrt(Ma*Ma + Mb*Mb)/sqrt(2)
  tree->Branch("gammaCM0", &m_gammaCM0); // 2*MQ/MS

  tree->Branch("MT2", &m_MT2);
  tree->Branch("RISR_LEP", &m_RISR_LEP);
  tree->Branch("PTISR_LEP", &m_PTISR_LEP);
  tree->Branch("MS_LEP", &m_MS_LEP);
  tree->Branch("MSV_LEP", &m_MSV_LEP);
  tree->Branch("MQ_LEP", &m_MQ_LEP);
  tree->Branch("gamma_LEP", &m_gamma_LEP);
  tree->Branch("MPa_LEP", &m_MPa_LEP);
  tree->Branch("MPb_LEP", &m_MPb_LEP);
  tree->Branch("MVa_LEP", &m_MVa_LEP);
  tree->Branch("MVb_LEP", &m_MVb_LEP);
  tree->Branch("PTS_CM_LEP", &m_PTS_CM_LEP);
  tree->Branch("MS_S0_LEP", &m_MS_S0_LEP);
  tree->Branch("MV_S0_LEP", &m_MV_S0_LEP);
  tree->Branch("MQ_S0_LEP", &m_MQ_S0_LEP);
  tree->Branch("gamma_S0_LEP", &m_gamma_S0_LEP);
  tree->Branch("MPTilde_LEP", &m_MPTilde_LEP);
  tree->Branch("MSTilde_LEP", &m_MSTilde_LEP);
  tree->Branch("gammaTilde_LEP", &m_gammaTilde_LEP);
  tree->Branch("CosDecayAngle_Pa_LEP", &m_CosDecayAngle_Pa_LEP);
  tree->Branch("CosDecayAngle_Pb_LEP", &m_CosDecayAngle_Pb_LEP);
  tree->Branch("CosDecayAngle_Va_LEP", &m_CosDecayAngle_Va_LEP);
  tree->Branch("CosDecayAngle_Vb_LEP", &m_CosDecayAngle_Vb_LEP);
  tree->Branch("CosDecayAngle_S_LEP", &m_CosDecayAngle_S_LEP);
  tree->Branch("RZPara_LEP", &m_RZPara_LEP);
  tree->Branch("MINV_LEP", &m_MINV_LEP);

  tree->Branch("Mperp_JET_ISR", &m_Mperp_JET_ISR);
  tree->Branch("gammaT_JET_ISR", &m_gammaT_JET_ISR);
  tree->Branch("RISR_JET_ISR", &m_RISR_JET_ISR);
  tree->Branch("PTISR_JET_ISR", &m_PTISR_JET_ISR);
  tree->Branch("MS_JET_ISR", &m_MS_JET_ISR);
  tree->Branch("MSV_JET_ISR", &m_MSV_JET_ISR);
  tree->Branch("MQ_JET_ISR", &m_MQ_JET_ISR);
  tree->Branch("gamma_JET_ISR", &m_gamma_JET_ISR);
  tree->Branch("MPa_JET_ISR", &m_MPa_JET_ISR);
  tree->Branch("MPb_JET_ISR", &m_MPb_JET_ISR);
  tree->Branch("MVa_JET_ISR", &m_MVa_JET_ISR);
  tree->Branch("MVb_JET_ISR", &m_MVb_JET_ISR);
  tree->Branch("PTS_CM_JET_ISR", &m_PTS_CM_JET_ISR);
  tree->Branch("MS_S0_JET_ISR", &m_MS_S0_JET_ISR);
  tree->Branch("MV_S0_JET_ISR", &m_MV_S0_JET_ISR);
  tree->Branch("MQ_S0_JET_ISR", &m_MQ_S0_JET_ISR);
  tree->Branch("gamma_S0_JET_ISR", &m_gamma_S0_JET_ISR);
  tree->Branch("MPTilde_JET_ISR", &m_MPTilde_JET_ISR);
  tree->Branch("MSTilde_JET_ISR", &m_MSTilde_JET_ISR);
  tree->Branch("gammaTilde_JET_ISR", &m_gammaTilde_JET_ISR);
  tree->Branch("CosDecayAngle_Pa_JET_ISR", &m_CosDecayAngle_Pa_JET_ISR);
  tree->Branch("CosDecayAngle_Pb_JET_ISR", &m_CosDecayAngle_Pb_JET_ISR);
  tree->Branch("CosDecayAngle_Va_JET_ISR", &m_CosDecayAngle_Va_JET_ISR);
  tree->Branch("CosDecayAngle_Vb_JET_ISR", &m_CosDecayAngle_Vb_JET_ISR);
  tree->Branch("CosDecayAngle_S_JET_ISR", &m_CosDecayAngle_S_JET_ISR);
  tree->Branch("RZPara_JET_ISR", &m_RZPara_JET_ISR);
  tree->Branch("MINV_JET_ISR", &m_MINV_JET_ISR);

  tree->Branch("MS_JET", &m_MS_JET);
  tree->Branch("MSV_JET", &m_MSV_JET);
  tree->Branch("MQ_JET", &m_MQ_JET);
  tree->Branch("gamma_JET", &m_gamma_JET);
  tree->Branch("MPa_JET", &m_MPa_JET);
  tree->Branch("MPb_JET", &m_MPb_JET);
  tree->Branch("MVa_JET", &m_MVa_JET);
  tree->Branch("MVb_JET", &m_MVb_JET);
  tree->Branch("MS_S0_JET", &m_MS_S0_JET);
  tree->Branch("MV_S0_JET", &m_MV_S0_JET);
  tree->Branch("MQ_S0_JET", &m_MQ_S0_JET);
  tree->Branch("gamma_S0_JET", &m_gamma_S0_JET);
  tree->Branch("MPTilde_JET", &m_MPTilde_JET);
  tree->Branch("MSTilde_JET", &m_MSTilde_JET);
  tree->Branch("gammaTilde_JET", &m_gammaTilde_JET);
  tree->Branch("CosDecayAngle_Pa_JET", &m_CosDecayAngle_Pa_JET);
  tree->Branch("CosDecayAngle_Pb_JET", &m_CosDecayAngle_Pb_JET);
  tree->Branch("CosDecayAngle_Va_JET", &m_CosDecayAngle_Va_JET);
  tree->Branch("CosDecayAngle_Vb_JET", &m_CosDecayAngle_Vb_JET);
  tree->Branch("CosDecayAngle_S_JET", &m_CosDecayAngle_S_JET);

  tree->Branch("MQV", &m_MQV);
  tree->Branch("gammaV", &m_gammaV);
  tree->Branch("MQV_LEP", &m_MQV_LEP);
  tree->Branch("gammaV_LEP", &m_gammaV_LEP);
  tree->Branch("MQV_JET_ISR", &m_MQV_JET_ISR);
  tree->Branch("gammaV_JET_ISR", &m_gammaV_JET_ISR);
  tree->Branch("MQV_JET", &m_MQV_JET);
  tree->Branch("gammaV_JET", &m_gammaV_JET);

  if(!do_slim){
    tree->Branch("MS", &m_MS);
    tree->Branch("PS", &m_PS);
    tree->Branch("cosS", &m_cosS);
    tree->Branch("dphiS", &m_dphiS);
    tree->Branch("dphiSI", &m_dphiSI);
    tree->Branch("PTS", &m_PTS);
    tree->Branch("PzS", &m_PzS);
    
    tree->Branch("EVa", &m_EVa);
    tree->Branch("EVb", &m_EVb);
    tree->Branch("PVa", &m_PVa);
    tree->Branch("PVb", &m_PVb);
    tree->Branch("EJa", &m_EJa);
    tree->Branch("EJb", &m_EJb);
    tree->Branch("PJa", &m_PJa);
    tree->Branch("PJb", &m_PJb);

    tree->Branch("MX2a", &m_MX2a);
    tree->Branch("MX2b", &m_MX2b);
    tree->Branch("ELa", &m_ELa);
    tree->Branch("ELb", &m_ELb);
    tree->Branch("PLa", &m_PLa);
    tree->Branch("PLb", &m_PLb);

    tree->Branch("MV", &m_MV);
    tree->Branch("PV", &m_PV);
    tree->Branch("MVa", &m_MVa);
    tree->Branch("MVb", &m_MVb);
    tree->Branch("PV_lab", &m_PV_lab);

    tree->Branch("MJa", &m_MJa);
    tree->Branch("MJb", &m_MJb);
    tree->Branch("MLa", &m_MLa);
    tree->Branch("MLb", &m_MLb);
    tree->Branch("cosJa", &m_cosJa);
    tree->Branch("cosJb", &m_cosJb);
    tree->Branch("cosLa", &m_cosLa);
    tree->Branch("cosLb", &m_cosLb);
    
    tree->Branch("MJ", &m_MJ);
    tree->Branch("ML", &m_ML);
    tree->Branch("EJ", &m_EJ);
    tree->Branch("EL", &m_EL);
    tree->Branch("PJ", &m_PJ);
    tree->Branch("PL", &m_PL);
    
    tree->Branch("PX2", &m_PX2);
    tree->Branch("PX2_BoostT", &m_PX2_BoostT);
    tree->Branch("MX2a_BoostT", &m_MX2a_BoostT);
    tree->Branch("MX2b_BoostT", &m_MX2b_BoostT);

    tree->Branch("PV_BoostT", &m_PV_BoostT);
    
    tree->Branch("EVa_BoostT", &m_EVa_BoostT);
    tree->Branch("EVb_BoostT", &m_EVb_BoostT);
    tree->Branch("PVa_BoostT", &m_PVa_BoostT);
    tree->Branch("PVb_BoostT", &m_PVb_BoostT);
    
    tree->Branch("PJ_BoostT", &m_PJ_BoostT);
    tree->Branch("PL_BoostT", &m_PL_BoostT);
  
    tree->Branch("PISR", &m_PISR);
    tree->Branch("RISRT", &m_RISRT);
    tree->Branch("MISR", &m_MISR);
  }

  if(!AnalysisBase<Base>::IsData()){
    tree->Branch("NPU", &m_NPU);
    tree->Branch("genMET", &m_genMET);
    tree->Branch("genMET_phi", &m_genMET_phi);
    
    tree->Branch("genNele", &m_genNele);
    tree->Branch("genNmu", &m_genNmu);
    
    tree->Branch("genNlep", &m_genNlep);
    tree->Branch("genPT_lep",  &m_genPT_lep);
    tree->Branch("genEta_lep", &m_genEta_lep);
    tree->Branch("genPhi_lep", &m_genPhi_lep);
    tree->Branch("genM_lep",   &m_genM_lep);
    tree->Branch("genCharge_lep",  &m_genCharge_lep);
    tree->Branch("genPDGID_lep",   &m_genPDGID_lep);
    tree->Branch("genMomPDGID_lep",   &m_genMomPDGID_lep);
    tree->Branch("genSourceID_lep",   &m_genSourceID_lep);
    tree->Branch("genIndex_lep",   &m_genIndex_lep);
    tree->Branch("genMomIndex_lep",   &m_genMomIndex_lep);
    
    tree->Branch("genNnu", &m_genNnu);
    tree->Branch("genPT_nu",  &m_genPT_nu);
    tree->Branch("genEta_nu", &m_genEta_nu);
    tree->Branch("genPhi_nu", &m_genPhi_nu);
    tree->Branch("genPDGID_nu",   &m_genPDGID_nu);
    tree->Branch("genMomPDGID_nu",   &m_genMomPDGID_nu);
    
    tree->Branch("genNboson", &m_genNboson);
    tree->Branch("genPT_boson",  &m_genPT_boson);
    tree->Branch("genEta_boson", &m_genEta_boson);
    tree->Branch("genPhi_boson", &m_genPhi_boson);
    tree->Branch("genM_boson",   &m_genM_boson);
    tree->Branch("genPDGID_boson",   &m_genPDGID_boson);
    tree->Branch("genMomPDGID_boson",   &m_genMomPDGID_boson);
    
    tree->Branch("genNsusy", &m_genNsusy);
    tree->Branch("genPT_susy",  &m_genPT_susy);
    tree->Branch("genEta_susy", &m_genEta_susy);
    tree->Branch("genPhi_susy", &m_genPhi_susy);
    tree->Branch("genM_susy",   &m_genM_susy);
    tree->Branch("genPDGID_susy", &m_genPDGID_susy);
    tree->Branch("genMomPDGID_susy", &m_genMomPDGID_susy);
  }

  if(AnalysisBase<Base>::IsCascades()){
    tree->Branch("cascades_tree", &m_cascades_tree);
    tree->Branch("cascades_prod", &m_cascades_prod);
    tree->Branch("cascades_SlepSneu_1stDecay", &m_cascades_SlepSneu_1stDecay);
    tree->Branch("cascades_SlepSneu_2ndDecay", &m_cascades_SlepSneu_2ndDecay);
    tree->Branch("cascades_N2_1stDecay", &m_cascades_N2_1stDecay);
    tree->Branch("cascades_N2_2ndDecay", &m_cascades_N2_2ndDecay);
  }

  if(AnalysisBase<Base>::IsCascades() || AnalysisBase<Base>::IsSMS()){
    tree->Branch("MSlepL", &m_MSlepL);
    tree->Branch("MSneu", &m_MSneu);
    tree->Branch("MN2", &m_MN2);
    tree->Branch("MC1", &m_MC1);
    tree->Branch("MN1", &m_MN1);
    tree->Branch("MP", &m_MP);
    tree->Branch("NSparticleW", &m_NSparticleW);
    tree->Branch("NSparticleZ", &m_NSparticleZ);
    tree->Branch("NSparticlePhoton", &m_NSparticlePhoton);
    tree->Branch("LSPParents", &m_LSPParents);
  }

  if(!AnalysisBase<Base>::IsData()){
    tree->Branch("Npartons", &m_Npartons);
    tree->Branch("genweight", &m_genweight);
    tree->Branch("XSec", &m_XSec);
    tree->Branch("FilterEff", &m_FilterEff);
    tree->Branch("Nweight", &m_Nweight);
    tree->Branch("Nevent", &m_Nevent);
  }

  return tree;
}

template <class Base>
void ReducedNtuple<Base>::ClearVariables(){

  m_event_skipped = false;
  m_treeSkipped.clear();
  for(int t = 0; t < m_aTrees; t++)
    m_treeSkipped.push_back(false);

  m_dphi_lep_S.clear();
  m_cos_lep_S.clear();
  m_dphi_jet_S.clear();
  m_cos_jet_S.clear();
  m_dphi_SV_S.clear();
  m_cos_SV_S.clear();

  m_dphiMET_lep_S.clear();
  m_dphiMET_jet_S.clear();
  m_dphiMET_SV_S.clear();

  m_Njet_ISR = 0;
  m_Njet_S = 0;
  m_Nbjet_ISR = 0;
  m_Nbjet_S = 0;
  m_Nlep_ISR = 0;
  m_Nlep_S = 0;
  m_NSV_ISR = 0;
  m_NSV_S = 0;
  m_index_jet_ISR.clear();
  m_index_jet_S.clear();
  m_index_SV_ISR.clear();
  m_index_SV_S.clear();
  m_index_lep_ISR.clear();
  m_index_lep_S.clear();
 
  m_Njet_a = 0;
  m_Njet_b = 0;
  m_Nbjet_a = 0;
  m_Nbjet_b = 0;
  m_Nlep_a = 0;
  m_Nlep_b = 0;
  m_NSV_a = 0;
  m_NSV_b = 0;
   
  m_index_jet_a.clear();
  m_index_jet_b.clear();
  m_index_lep_a.clear();
  m_index_lep_b.clear();
  m_index_SV_a.clear();
  m_index_SV_b.clear();

  m_Nlep_a_LEP = 0;
  m_Nlep_b_LEP = 0;
  m_index_lep_a_LEP.clear();
  m_index_lep_b_LEP.clear();
  m_Njet_ISR_JET_ISR = 0;
  m_Njet_S_JET_ISR = 0;
  m_Nbjet_ISR_JET_ISR = 0;
  m_Nbjet_S_JET_ISR = 0;
  m_Njet_a_JET_ISR = 0;
  m_Njet_b_JET_ISR = 0;
  m_Nbjet_a_JET_ISR = 0;
  m_Nbjet_b_JET_ISR = 0;
  m_index_jet_a_JET_ISR.clear();
  m_index_jet_b_JET_ISR.clear();
  m_index_jet_ISR_JET_ISR.clear();
  m_index_jet_S_JET_ISR.clear();
  m_Njet_a_JET = 0;
  m_Njet_b_JET = 0;
  m_Nbjet_a_JET = 0;
  m_Nbjet_b_JET = 0;
  m_index_jet_a_JET.clear();
  m_index_jet_b_JET.clear();

}

template <class Base>
void ReducedNtuple<Base>::FillOutputTree(TTree* tree, const Systematic& sys, bool do_slim){
  
  AnalysisBase<Base>::SetSystematic(sys);

  m_EventFilter = AnalysisBase<Base>::PassEventFilter();

  if(AnalysisBase<Base>::IsData())
    if(!AnalysisBase<Base>::IsGoodEvent())
      return;

  bool good_PV;
  TVector3 PV = AnalysisBase<Base>::GetPV(good_PV);

  if(!good_PV)
    return;

  TVector3 ETMiss;
  ParticleList Jets_noID = AnalysisBase<Base>::GetJetsMET(ETMiss);
  ParticleList Jets      = AnalysisBase<Base>::GetJetsMET(ETMiss, 3); // jet ID 3
  ParticleList GenJets   = AnalysisBase<Base>::GetGenJets();

  m_FastSimEventVeto = AnalysisBase<Base>::FastSimEventVeto(GenJets);
  
  m_PrefireWeight = AnalysisBase<Base>::GetPrefireWeight(0);
  m_PrefireWeight_up = AnalysisBase<Base>::GetPrefireWeight(1);
  m_PrefireWeight_down = AnalysisBase<Base>::GetPrefireWeight(-1);
  
  if(std::isnan(ETMiss.Mag())) return;
  if(ETMiss.Mag() < 0.)
    return;
  
  ClearVariables();

  if(Jets_noID.size() != Jets.size())
    m_EventFlag_FailJetID = true;
  else
    m_EventFlag_FailJetID = false;

  int Njet      = Jets.size();
  int Njet_noID = Jets_noID.size();

  // any jets in HEM region?
  m_EventFlag_JetInHEM = false;
  for(int j = 0; j < Njet_noID; j++)
    if(AnalysisBase<Base>::IsHEM(Jets_noID[j]))
      m_EventFlag_JetInHEM = true;

  Jets      = Jets.PtEtaCut(20., 5.);
  Jets_noID = Jets_noID.PtEtaCut(20., 5.);
  GenJets   = GenJets.PtEtaCut(20., 5.);
  Njet      = Jets.size();
  Njet_noID = Jets_noID.size();

  m_EventFlag_JetInHEM_Pt20 = false;
  m_HT_eta5 = 0.;
  for(int j = 0; j < Njet_noID; j++){
    if(AnalysisBase<Base>::IsHEM(Jets_noID[j]))
      m_EventFlag_JetInHEM_Pt20 = true;
    m_HT_eta5 += Jets_noID[j].Pt();
  }

  m_EventFlag_JetInHEM_Pt20_JetID = false;
  m_HT_eta5_id = 0.;
  for(int j = 0; j < Njet; j++){
    if(AnalysisBase<Base>::IsHEM(Jets[j]))
      m_EventFlag_JetInHEM_Pt20_JetID = true;
    m_HT_eta5_id += Jets[j].Pt();
  }

  Jets      = Jets.PtEtaCut(20., 3.);
  Jets_noID = Jets_noID.PtEtaCut(20., 3.);
  Njet      = Jets.size();
  Njet_noID = Jets_noID.size();

  m_HT_eta3 = 0.;
  for(int j = 0; j < Njet_noID; j++){
    m_HT_eta3 += Jets_noID[j].Pt();
  }

  m_HT_eta3_id = 0.;
  for(int j = 0; j < Njet; j++){
    m_HT_eta3_id += Jets[j].Pt();
  }

  Jets      = Jets.PtEtaCut(20., 2.4);
  Jets_noID = Jets_noID.PtEtaCut(20., 2.4);
  Njet      = Jets.size();
  Njet_noID = Jets_noID.size();

  m_HT_eta24 = 0.;
  for(int j = 0; j < Njet_noID; j++){
    m_HT_eta24 += Jets_noID[j].Pt();
  }

  m_HT_eta24_id = 0.;
  for(int j = 0; j < Njet; j++){
    m_HT_eta24_id += Jets[j].Pt();
  }

  m_LHE_HT = AnalysisBase<Base>::Get_LHE_HT();
  m_LHE_HTIncoming = AnalysisBase<Base>::Get_LHE_HTIncoming();
  
  ParticleList Muons = AnalysisBase<Base>::GetMuons();
  Muons = Muons.ParticleIDCut(kVeryLoose);
  Muons = Muons.PtEtaCut(3.0);

  ParticleList Electrons = AnalysisBase<Base>::GetElectrons();
  Electrons = Electrons.ParticleIDCut(kVeryLoose);
  Electrons = Electrons.PtEtaCut(5.0);

  ParticleList LowPtElectrons = AnalysisBase<Base>::GetLowPtElectrons();
  LowPtElectrons = LowPtElectrons.ParticleIDCut(kVeryLoose);
  LowPtElectrons = LowPtElectrons.PtEtaCut(2.0); // upper cut in AnalysisBase functions

  ParticleList Leptons = Electrons+LowPtElectrons+Muons;
  Leptons.SortByPt();

  Jets = Jets.RemoveOverlap(Leptons, 0.2);

  ParticleList GenMuons = AnalysisBase<Base>::GetGenMuons();
  ParticleList GenElectrons = AnalysisBase<Base>::GetGenElectrons();
  ParticleList GenLeptons = GenElectrons+GenMuons;
  GenLeptons.SortByPt();

  m_genNele = GenElectrons.size();
  m_genNmu  = GenMuons.size();
  m_genNlep = GenLeptons.size();

  GenJets = GenJets.RemoveOverlap(GenLeptons, 0.2);

  // merge jets until total number is combinatorically manageable
  Jets = Jets.BinaryMerge(13-Leptons.size());
  Jets = Jets.PtEtaCut(20., 2.4);
  Jets.SortByPt();
  GenJets.SortByPt();
  
  m_Njet = Jets.size();
  m_NGenjet = GenJets.size();  

  ParticleList BJets;
  vector<int> BJets_index;
  for(int i = 0; i < m_Njet; i++){
    if(Jets[i].BtagID() >= kMedium){
      BJets += Jets[i];
      BJets_index.push_back(i);
    }
  }
  m_Nbjet = BJets.size();

  m_Nele = Electrons.size();
  m_Nlele = LowPtElectrons.size();
  m_Nmu  = Muons.size();
  m_Nlep = Leptons.size();
  
  // NTUPLE Selection
  if(m_Nlep < 2 && ETMiss.Mag() < 150)
    return;
  if(AnalysisBase<Base>::IsCascades() || AnalysisBase<Base>::IsSMS()){
    std::pair<int,int> temp_masses = AnalysisBase<Base>::GetSUSYMasses();
    m_MP = temp_masses.first;
    m_MN1 = temp_masses.second;
    if(1.*m_MN1/m_MP < 0.45) return; // only keep compressed signals
  }

  m_HEM_Veto = m_EventFlag_JetInHEM_Pt20;
  for(int l = 0; l < m_Nlep; l++)
    if(AnalysisBase<Base>::IsHEM(Leptons[l]))
      m_HEM_Veto = true;

  TVector3 altETMiss = AnalysisBase<Base>::GetAltMET();
  m_altMET     = altETMiss.Pt();
  m_altMET_phi = altETMiss.Phi();

  if(m_Nlep >= 2){
    double pa[3] = {Leptons[0].M(), Leptons[0].Px(), Leptons[0].Py()};
    double pb[3] = {Leptons[1].M(), Leptons[1].Px(), Leptons[1].Py()};
    double pmiss[3] = {0.0, ETMiss.Px(), ETMiss.Py()};

    mt2_bisect::mt2 mt2calc;
    mt2calc.set_momenta(pa, pb, pmiss);
    mt2calc.set_mn(0.); // assume massless invisible particles like neutrinos
    m_MT2 = mt2calc.get_mt2();
  }
  
  // Sparticle pair-production trees analysis
  for(int t = 0; t < m_aTrees; t++){
    
    LAB[t]->ClearEvent(); 

    INV[t]->SetLabFrameThreeVector(ETMiss);
    
    // ISR analysis    
    if(t==0){
      if(m_Njet == 0){ // can't do ISR analysis without at least one jet
        m_treeSkipped[t] = true;
        continue;
      }
      std::vector<RFKey> lepID;
      for(int i = 0; i < m_Nlep; i++){
        lepID.push_back(COMB_L[t]->AddLabFrameFourVector(Leptons[i]));
      }
      std::vector<RFKey> jetID;
      for(int i = 0; i < m_Njet; i++){
        jetID.push_back(COMB_J[t]->AddLabFrameFourVector(Jets[i]));
      }
      if(!LAB[t]->AnalyzeEvent())
        cout << "Something went wrong with tree event analysis" << endl;

      bool BSideIsA = false;
      if((Lb[t]->GetFourVector() + Jb[t]->GetFourVector()).M() > (La[t]->GetFourVector() + Ja[t]->GetFourVector()).M()) BSideIsA = true;

      // jet counting in ISR/S, hemispheres
      for(int i = 0; i < m_Njet; i++){
        if(COMB_J[t]->GetFrame(jetID[i]) == *ISR[t]){
          m_Njet_ISR++;
          if(Jets[i].BtagID() >= kMedium)
            m_Nbjet_ISR++;
          m_index_jet_ISR.push_back(i);
        }
        if(COMB_J[t]->GetFrame(jetID[i]) == *Ja[t]){
          if(!BSideIsA){
            m_Njet_S++;
            m_Njet_a++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S++;
              m_Nbjet_a++;
            }
            m_index_jet_S.push_back(i);
            m_index_jet_a.push_back(i);
          }
          else{
            m_Njet_S++;
            m_Njet_b++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S++;
              m_Nbjet_b++;
            }
            m_index_jet_S.push_back(i);
            m_index_jet_b.push_back(i);
          }
        }
        if(COMB_J[t]->GetFrame(jetID[i]) == *Jb[t]){
          if(!BSideIsA){
            m_Njet_S++;
            m_Njet_b++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S++;
              m_Nbjet_b++;
            }
            m_index_jet_S.push_back(i);
            m_index_jet_b.push_back(i);
          }
          else{
            m_Njet_S++;
            m_Njet_a++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S++;
              m_Nbjet_a++;
            }
            m_index_jet_S.push_back(i);
            m_index_jet_a.push_back(i);
          }
        }
      }
          
      // lepton counting in ISR/S, hemispheres
      for(int i = 0; i < m_Nlep; i++){
        m_Nlep_S++;
        if(COMB_L[t]->GetFrame(lepID[i]) == *La[t]){
          if(!BSideIsA){
            m_Nlep_a++;
            m_index_lep_a.push_back(i);
            m_index_lep_S.push_back(i);
          }
          else{
            m_Nlep_b++;
            m_index_lep_b.push_back(i);
            m_index_lep_S.push_back(i);
          }
        }
        if(COMB_L[t]->GetFrame(lepID[i]) == *Lb[t]){
          if(!BSideIsA){
            m_Nlep_b++;
            m_index_lep_b.push_back(i);
            m_index_lep_S.push_back(i);
          }
          else{
            m_Nlep_a++;
            m_index_lep_a.push_back(i);
            m_index_lep_S.push_back(i);
          }
        }
      }

      if(m_Nlep + m_Njet_S < 1){
        m_treeSkipped[t] = true;
        continue;
      }
      
      // Fill Observable Branches
      m_PTCM = CM[t]->GetFourVector().Pt();
      m_EtaCM = CM[t]->GetFourVector().Eta();
      m_PhiCM = CM[t]->GetFourVector().Phi();
      m_MCM = CM[t]->GetFourVector().M();
      m_PzCM = CM[t]->GetFourVector().Pz();
      m_cosCM = CM[t]->GetCosDecayAngle();
      m_dphiCM = CM[t]->GetDeltaPhiDecayAngle();
      m_dphiCMI = CM[t]->GetDeltaPhiBoostVisible();
      m_PS = S[t]->GetMomentum(*CM[t]);
      
      m_MS = S[t]->GetMass();
      m_cosS = S[t]->GetCosDecayAngle();
      m_dphiS = S[t]->GetDeltaPhiDecayAngle();
      m_dphiSI = S[t]->GetDeltaPhiBoostVisible();
      m_PTS = S[t]->GetFourVector().Pt();
      m_EtaS = S[t]->GetFourVector().Eta();
      m_PhiS = S[t]->GetFourVector().Phi();
      m_PzS = S[t]->GetFourVector().Pz();

      m_LAB_Pt = LAB[t]->GetFourVector().Pt();
      m_LAB_Eta = LAB[t]->GetFourVector().Eta();
      m_LAB_Phi = LAB[t]->GetFourVector().Phi();
      m_LAB_M = LAB[t]->GetFourVector().M();
      
      m_EVa = X2a[t]->GetListVisibleFrames().GetFourVector(*X2a[t]).E();
      m_EVb = X2b[t]->GetListVisibleFrames().GetFourVector(*X2b[t]).E();
      m_PVa = X2a[t]->GetListVisibleFrames().GetFourVector(*X2a[t]).P();
      m_PVb = X2b[t]->GetListVisibleFrames().GetFourVector(*X2b[t]).P();
      m_EJa = Ja[t]->GetFourVector(*X2a[t]).E();
      m_EJb = Jb[t]->GetFourVector(*X2b[t]).E();
      m_PJa = Ja[t]->GetFourVector(*X2a[t]).P();
      m_PJb = Jb[t]->GetFourVector(*X2b[t]).P();

      if(BSideIsA){
        m_EVa = X2b[t]->GetListVisibleFrames().GetFourVector(*X2b[t]).E();
        m_EVb = X2a[t]->GetListVisibleFrames().GetFourVector(*X2a[t]).E();
        m_PVa = X2b[t]->GetListVisibleFrames().GetFourVector(*X2b[t]).P();
        m_PVb = X2a[t]->GetListVisibleFrames().GetFourVector(*X2a[t]).P();
        m_EJa = Jb[t]->GetFourVector(*X2b[t]).E();
        m_EJb = Ja[t]->GetFourVector(*X2a[t]).E();
        m_PJa = Jb[t]->GetFourVector(*X2b[t]).P();
        m_PJb = Ja[t]->GetFourVector(*X2a[t]).P();
      }
      
      m_MX2a = X2a[t]->GetMass();
      //m_cosX2a = X2a[t]->GetCosDecayAngle();
      m_MX2b = X2b[t]->GetMass();
      //m_cosX2b = X2b[t]->GetCosDecayAngle();
      m_ELa = La[t]->GetFourVector(*X2a[t]).E();
      m_ELb = Lb[t]->GetFourVector(*X2b[t]).E();;
      m_PLa = La[t]->GetFourVector(*X2a[t]).P();
      m_PLb = Lb[t]->GetFourVector(*X2b[t]).P();

      if(BSideIsA){
        m_MX2a = X2b[t]->GetMass();
        //m_cosX2a = X2b[t]->GetCosDecayAngle();
        m_MX2b = X2a[t]->GetMass();
        //m_cosX2b = X2a[t]->GetCosDecayAngle();
        m_ELa = Lb[t]->GetFourVector(*X2b[t]).E();
        m_ELb = La[t]->GetFourVector(*X2a[t]).E();;
        m_PLa = Lb[t]->GetFourVector(*X2b[t]).P();
        m_PLb = La[t]->GetFourVector(*X2a[t]).P();
      }
      
      m_MV = S[t]->GetListVisibleFrames().GetMass();
      m_PV = S[t]->GetListVisibleFrames().GetFourVector(*S[t]).P();
      m_MVa = X2a[t]->GetListVisibleFrames().GetMass();
      m_MVb = X2b[t]->GetListVisibleFrames().GetMass();
      if(BSideIsA){
        m_MVa = X2b[t]->GetListVisibleFrames().GetMass();
        m_MVb = X2a[t]->GetListVisibleFrames().GetMass();
      }

      m_PV_lab    = S[t]->GetListVisibleFrames().GetFourVector().P();
      m_dphiMET_V = S[t]->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);
      
      m_MJa = Ja[t]->GetMass();
      m_MJb = Jb[t]->GetMass();
      m_MLa = La[t]->GetMass();
      m_MLb = Lb[t]->GetMass();
      m_cosJa = sJa[t]->GetCosDecayAngle(*LAB[t]);
      m_cosJb = sJb[t]->GetCosDecayAngle(*LAB[t]);
      m_cosLa = sLa[t]->GetCosDecayAngle(*LAB[t]);
      m_cosLb = sLb[t]->GetCosDecayAngle(*LAB[t]);

      if(BSideIsA){
        m_MJa = Jb[t]->GetMass();
        m_MJb = Ja[t]->GetMass();
        m_MLa = Lb[t]->GetMass();
        m_MLb = La[t]->GetMass();
        m_cosJa = sJb[t]->GetCosDecayAngle(*LAB[t]);
        m_cosJb = sJa[t]->GetCosDecayAngle(*LAB[t]);
        m_cosLa = sLb[t]->GetCosDecayAngle(*LAB[t]);
        m_cosLb = sLa[t]->GetCosDecayAngle(*LAB[t]);
      }
      
      TLorentzVector vP_S_CM  = S[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Ja_S  = Ja[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Jb_S  = Jb[t]->GetFourVector(*S[t]);
      TLorentzVector vP_La_S  = La[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Lb_S  = Lb[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Ia_S  = X1a[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Ib_S  = X1b[t]->GetFourVector(*S[t]);
      
      TLorentzVector vP_Ja_X2a  = Ja[t]->GetFourVector(*X2a[t]);
      TLorentzVector vP_Jb_X2b  = Jb[t]->GetFourVector(*X2b[t]);
      TLorentzVector vP_La_X2a  = La[t]->GetFourVector(*X2a[t]);
      TLorentzVector vP_Lb_X2b  = Lb[t]->GetFourVector(*X2b[t]);
      TLorentzVector vP_Ia_X2a  = X1a[t]->GetFourVector(*X2a[t]);
      TLorentzVector vP_Ib_X2b  = X1b[t]->GetFourVector(*X2b[t]);

      m_MJ = (vP_Ja_S+vP_Jb_S).M();
      m_ML = (vP_La_S+vP_Lb_S).M();
      m_EJ = (vP_Ja_S+vP_Jb_S).E();
      m_EL = (vP_La_S+vP_Lb_S).E();
      m_PJ = (vP_Ja_S+vP_Jb_S).P();
      m_PL = (vP_La_S+vP_Lb_S).P();
      
      m_PX2 = (vP_Ja_S+vP_La_S+vP_Ia_S).P();
      
      // removing momentum components parallel to CM->S boost
      TVector3 boostVis = (vP_Ja_S+vP_La_S+vP_Jb_S+vP_Lb_S).BoostVector();
      TVector3 boostInv = (vP_Ia_S+vP_Ib_S).BoostVector();
      TVector3 daBoost = vP_S_CM.Vect().Unit();
      
      if((std::isnan(boostInv.Mag()) || std::isnan(boostVis.Mag())))
        cout << "boost NAN " << boostInv.Mag() << " " << boostVis.Mag() << " Ja=" << vP_Ja_S.P() << " Jb=" << vP_Jb_S.P() << " La=" << vP_La_S.P() << " Lb=" << vP_Lb_S.P() << " Inva=" << vP_Ia_S.P() << " Invb=" << vP_Ib_S.P() << " MET=" << ETMiss.Mag() << endl;
      
      boostVis = (boostVis.Dot(daBoost))*daBoost;
      boostInv = (boostInv.Dot(daBoost))*daBoost;

      if((!std::isnan(boostVis.Mag())) &&
         (boostVis.Mag() < 1)){
        vP_Ja_S.Boost(-boostVis);
        vP_Jb_S.Boost(-boostVis);
        vP_La_S.Boost(-boostVis);
        vP_Lb_S.Boost(-boostVis);
      } else {
        vP_Ja_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ja_S.M()));
        vP_Jb_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Jb_S.M()));
        vP_La_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_La_S.M()));
        vP_Lb_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Lb_S.M()));
      }
      if((!std::isnan(boostInv.Mag())) &&
         (boostInv.Mag() < 1)){
        vP_Ia_S.Boost(-boostInv);
        vP_Ib_S.Boost(-boostInv);
      } else {
        vP_Ia_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ia_S.M()));
        vP_Ib_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ib_S.M()));
      }
        
      m_PX2_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).P();
      m_MX2a_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).M();
      m_MX2b_BoostT = (vP_Jb_S+vP_Lb_S+vP_Ib_S).M();
      if(BSideIsA){
        m_MX2a_BoostT = (vP_Jb_S+vP_Lb_S+vP_Ib_S).M();
        m_MX2b_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).M();
      }
      m_Mperp = sqrt(m_MX2a_BoostT*m_MX2a_BoostT+m_MX2b_BoostT*m_MX2b_BoostT)/sqrt(2.);
      m_gammaT = 2*m_Mperp/(sqrt(m_MX2a_BoostT*m_MX2a_BoostT+m_PX2_BoostT*m_PX2_BoostT) +
            		sqrt(m_MX2b_BoostT*m_MX2b_BoostT+m_PX2_BoostT*m_PX2_BoostT));

      m_PV_BoostT = (vP_Ja_S+vP_La_S+vP_Jb_S+vP_Lb_S).P();
      
      m_EVa_BoostT = (vP_Ja_S+vP_La_S).E();
      m_EVb_BoostT = (vP_Jb_S+vP_Lb_S).E();
      m_PVa_BoostT = (vP_Ja_S+vP_La_S).P();
      m_PVb_BoostT = (vP_Jb_S+vP_Lb_S).P();

      if(BSideIsA){
        m_EVa_BoostT = (vP_Jb_S+vP_Lb_S).E();
        m_EVb_BoostT = (vP_Ja_S+vP_La_S).E();
        m_PVa_BoostT = (vP_Jb_S+vP_Lb_S).P();
        m_PVb_BoostT = (vP_Ja_S+vP_La_S).P();
      }

      m_EJ_BoostT = (vP_Ja_S+vP_Jb_S).E();
      m_EL_BoostT = (vP_La_S+vP_Lb_S).E();
      m_PJ_BoostT = (vP_Ja_S+vP_Jb_S).P();
      m_PL_BoostT = (vP_La_S+vP_Lb_S).P();

      // ISR related variables
      if(m_Njet_ISR > 0){
        TVector3 vPISR = S[t]->GetFourVector(*CM[t]).Vect();
        TVector3 vPINV = (X1a[t]->GetFourVector(*CM[t])+X1b[t]->GetFourVector(*CM[t])).Vect();

        m_MISR = ISR[t]->GetMass();
        m_RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();

        m_PISR = vPISR.Mag();
        // transverse
        TVector3 vPTISR = S[t]->GetTransverseFourVector(*CM[t]).Vect();
        TVector3 vPTINV = (X1a[t]->GetTransverseFourVector(*CM[t])+X1b[t]->GetTransverseFourVector(*CM[t])).Vect();

        m_PTISR  = vPTISR.Mag();
        m_RISRT  = fabs(vPTINV.Dot(vPTISR.Unit())) / vPTISR.Mag();
        
        if((std::isnan(m_RISR) || std::isnan(m_Mperp))){
          cout << "VAR NAN " << vPTISR.Mag() << " " << vPTINV.Mag() << " NjetS=" << m_Njet_S << " Njeta=" << m_Njet_a << " Njetb=" << m_Njet_b << " Nlep=" << m_Nlep << " Nlepa=" << m_Nlep_a << " Nlepb=" << m_Nlep_b << " Mperp=" << m_Mperp << " RISR=" << m_RISR << endl;
        }

        m_PISR = vPISR.Mag();
        TVector3 isr_t = ISR[t]->GetTransverseFourVector(*S[t]).Vect();
        TVector3 isr   = ISR[t]->GetFourVector(*S[t]).Vect();
        
        for(int l = 0; l < m_Nlep_S; l++){
          TVector3 lep   = S[t]->GetFourVector(Leptons[m_index_lep_S[l]]).Vect();
          TVector3 lep_t = S[t]->GetTransverseFourVector(Leptons[m_index_lep_S[l]]).Vect();
          m_dphi_lep_S.push_back( lep_t.Angle(isr_t) );
          m_cos_lep_S.push_back( lep.Unit().Dot(isr.Unit()) );
          m_dphiMET_lep_S.push_back( Leptons[m_index_lep_S[l]].Vect().DeltaPhi(ETMiss) );
        }
        m_PISR = vPISR.Mag();
        for(int l = 0; l < m_Njet_S; l++){
          TVector3 jet   = S[t]->GetFourVector(Jets[m_index_jet_S[l]]).Vect();
          TVector3 jet_t = S[t]->GetTransverseFourVector(Jets[m_index_jet_S[l]]).Vect();
          m_dphi_jet_S.push_back( jet_t.Angle(isr_t) );
          m_cos_jet_S.push_back( jet.Unit().Dot(isr.Unit()) );
          m_dphiMET_jet_S.push_back( Jets[m_index_jet_S[l]].Vect().DeltaPhi(ETMiss) );
        }
      }
      // New Observables for Run3/Cascades
      // analysis using ISR boosted tree
      // variables using 4-vects constructed by evaluating in CM frame and making perp to boost
      TLorentzVector vP_Ja_CM = Ja[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Jb_CM = Jb[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_La_CM = La[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Lb_CM = Lb[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Ia_CM = X1a[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Ib_CM = X1b[t]->GetFourVector(*CM[t]);

      TVector3 vP_Ja_CM_Perp = vP_Ja_CM.Vect();
      TVector3 vP_Jb_CM_Perp = vP_Jb_CM.Vect();
      TVector3 vP_La_CM_Perp = vP_La_CM.Vect();
      TVector3 vP_Lb_CM_Perp = vP_Lb_CM.Vect();
      TVector3 vP_Ia_CM_Perp = vP_Ia_CM.Vect();
      TVector3 vP_Ib_CM_Perp = vP_Ib_CM.Vect();

      vP_Ja_CM_Perp = vP_Ja_CM_Perp - vP_Ja_CM_Perp.Dot(daBoost)*daBoost;
      vP_Jb_CM_Perp = vP_Jb_CM_Perp - vP_Jb_CM_Perp.Dot(daBoost)*daBoost;
      vP_La_CM_Perp = vP_La_CM_Perp - vP_La_CM_Perp.Dot(daBoost)*daBoost;
      vP_Lb_CM_Perp = vP_Lb_CM_Perp - vP_Lb_CM_Perp.Dot(daBoost)*daBoost;
      vP_Ia_CM_Perp = vP_Ia_CM_Perp - vP_Ia_CM_Perp.Dot(daBoost)*daBoost;
      vP_Ib_CM_Perp = vP_Ib_CM_Perp - vP_Ib_CM_Perp.Dot(daBoost)*daBoost;

      TLorentzVector vP_Ja_CM_PerpM0_TLV;
      TLorentzVector vP_Jb_CM_PerpM0_TLV;
      TLorentzVector vP_La_CM_PerpM0_TLV;
      TLorentzVector vP_Lb_CM_PerpM0_TLV;
      TLorentzVector vP_Ia_CM_PerpM0_TLV;
      TLorentzVector vP_Ib_CM_PerpM0_TLV;
      vP_Ja_CM_PerpM0_TLV.SetPtEtaPhiM(vP_Ja_CM_Perp.Pt(),vP_Ja_CM_Perp.Eta(),vP_Ja_CM_Perp.Phi(),0.);
      vP_Jb_CM_PerpM0_TLV.SetPtEtaPhiM(vP_Jb_CM_Perp.Pt(),vP_Jb_CM_Perp.Eta(),vP_Jb_CM_Perp.Phi(),0.);
      vP_La_CM_PerpM0_TLV.SetPtEtaPhiM(vP_La_CM_Perp.Pt(),vP_La_CM_Perp.Eta(),vP_La_CM_Perp.Phi(),0.);
      vP_Lb_CM_PerpM0_TLV.SetPtEtaPhiM(vP_Lb_CM_Perp.Pt(),vP_Lb_CM_Perp.Eta(),vP_Lb_CM_Perp.Phi(),0.);
      vP_Ia_CM_PerpM0_TLV.SetPtEtaPhiM(vP_Ia_CM_Perp.Pt(),vP_Ia_CM_Perp.Eta(),vP_Ia_CM_Perp.Phi(),0.);
      vP_Ib_CM_PerpM0_TLV.SetPtEtaPhiM(vP_Ib_CM_Perp.Pt(),vP_Ib_CM_Perp.Eta(),vP_Ib_CM_Perp.Phi(),0.);

      TLorentzVector vP_Ja_CM_Perp_TLV;
      TLorentzVector vP_Jb_CM_Perp_TLV;
      TLorentzVector vP_La_CM_Perp_TLV;
      TLorentzVector vP_Lb_CM_Perp_TLV;
      TLorentzVector vP_Ia_CM_Perp_TLV;
      TLorentzVector vP_Ib_CM_Perp_TLV;
      vP_Ja_CM_Perp_TLV.SetPtEtaPhiM(vP_Ja_CM_Perp.Pt(),vP_Ja_CM_Perp.Eta(),vP_Ja_CM_Perp.Phi(),vP_Ja_CM.M());
      vP_Jb_CM_Perp_TLV.SetPtEtaPhiM(vP_Jb_CM_Perp.Pt(),vP_Jb_CM_Perp.Eta(),vP_Jb_CM_Perp.Phi(),vP_Jb_CM.M());
      vP_La_CM_Perp_TLV.SetPtEtaPhiM(vP_La_CM_Perp.Pt(),vP_La_CM_Perp.Eta(),vP_La_CM_Perp.Phi(),vP_La_CM.M());
      vP_Lb_CM_Perp_TLV.SetPtEtaPhiM(vP_Lb_CM_Perp.Pt(),vP_Lb_CM_Perp.Eta(),vP_Lb_CM_Perp.Phi(),vP_Lb_CM.M());
      vP_Ia_CM_Perp_TLV.SetPtEtaPhiM(vP_Ia_CM_Perp.Pt(),vP_Ia_CM_Perp.Eta(),vP_Ia_CM_Perp.Phi(),vP_Ia_CM.M());
      vP_Ib_CM_Perp_TLV.SetPtEtaPhiM(vP_Ib_CM_Perp.Pt(),vP_Ib_CM_Perp.Eta(),vP_Ib_CM_Perp.Phi(),vP_Ib_CM.M());

      m_MSperpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_Jb_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV + vP_Ia_CM_PerpM0_TLV + vP_Ib_CM_PerpM0_TLV).Mag();
      m_MaPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV + vP_Ia_CM_PerpM0_TLV).Mag();
      m_MaVPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV).Mag();
      m_MbPerpCM0 = (vP_Jb_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV + vP_Ib_CM_PerpM0_TLV).Mag();
      m_MbVPerpCM0 = (vP_Jb_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV).Mag();
      if(BSideIsA){
        m_MaPerpCM0 = (vP_Jb_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV + vP_Ib_CM_PerpM0_TLV).Mag();
        m_MaVPerpCM0 = (vP_Jb_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV).Mag();
        m_MbPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV + vP_Ia_CM_PerpM0_TLV).Mag();
        m_MbVPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV).Mag();
      }
      m_MQperpCM0 = sqrt((m_MaPerpCM0*m_MaPerpCM0 + m_MbPerpCM0*m_MbPerpCM0)/2.); 
      m_gammaPerpCM0 = 2*m_MQperpCM0/m_MSperpCM0;
      // variables using 4-vects constructed by evaluating in CM frame
      TLorentzVector vP_Ja_CM0;
      TLorentzVector vP_Jb_CM0;
      TLorentzVector vP_La_CM0;
      TLorentzVector vP_Lb_CM0;
      TLorentzVector vP_Ia_CM0;
      TLorentzVector vP_Ib_CM0;

      vP_Ja_CM0.SetPtEtaPhiM(vP_Ja_CM.Pt(), vP_Ja_CM.Eta(), vP_Ja_CM.Phi(), 0.);
      vP_Jb_CM0.SetPtEtaPhiM(vP_Jb_CM.Pt(), vP_Jb_CM.Eta(), vP_Jb_CM.Phi(), 0.);
      vP_La_CM0.SetPtEtaPhiM(vP_La_CM.Pt(), vP_La_CM.Eta(), vP_La_CM.Phi(), 0.);
      vP_Lb_CM0.SetPtEtaPhiM(vP_Lb_CM.Pt(), vP_Lb_CM.Eta(), vP_Lb_CM.Phi(), 0.);
      vP_Ia_CM0.SetPtEtaPhiM(vP_Ia_CM.Pt(), vP_Ia_CM.Eta(), vP_Ia_CM.Phi(), 0.);
      vP_Ib_CM0.SetPtEtaPhiM(vP_Ib_CM.Pt(), vP_Ib_CM.Eta(), vP_Ib_CM.Phi(), 0.);

      m_MSCM0 = (vP_Ja_CM0 + vP_Jb_CM0 + vP_La_CM0 + vP_Lb_CM0 + vP_Ia_CM0 + vP_Ib_CM0).Mag(); // this is Mr or 'ass mass' in the LLP analysis
      m_MaCM0 = (vP_Ja_CM0 + vP_La_CM0 + vP_Ia_CM0).Mag();
      m_MbCM0 = (vP_Jb_CM0 + vP_Lb_CM0 + vP_Ib_CM0).Mag();
      m_MaVCM0 = (vP_Ja_CM0 + vP_La_CM0).Mag();
      m_MbVCM0 = (vP_Jb_CM0 + vP_Lb_CM0).Mag();
      if(BSideIsA){
        m_MaCM0 = (vP_Jb_CM0 + vP_Lb_CM0 + vP_Ib_CM0).Mag();
        m_MbCM0 = (vP_Ja_CM0 + vP_La_CM0 + vP_Ia_CM0).Mag();
        m_MaVCM0 = (vP_Jb_CM0 + vP_Lb_CM0).Mag();
        m_MbVCM0 = (vP_Ja_CM0 + vP_La_CM0).Mag();
      }
      m_MQCM0 = sqrt((m_MaCM0*m_MaCM0 + m_MbCM0*m_MbCM0)/2.); 
      m_gammaCM0 = 2.*m_MQCM0/m_MSCM0; // this is 'R' in the LLP analysis
      m_MQV = sqrt(((Ja[t]->GetFourVector()+La[t]->GetFourVector()).M2() + (Jb[t]->GetFourVector()+Lb[t]->GetFourVector()).M2())/2.);
      m_gammaV = 2.*m_MQV/m_MSCM0; // this is 'Rv' in the LLP analysis

      m_CosDecayAngle_S = S[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pa = X2a[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pb = X2b[t]->GetCosDecayAngle();
      if(BSideIsA){
        m_CosDecayAngle_Pa = X2b[t]->GetCosDecayAngle();
        m_CosDecayAngle_Pb = X2a[t]->GetCosDecayAngle();
      }

      ParticleList cand_list_parts_a;
      for(int i = 0; i < m_Nlep_a; i++){
        Particle part = Leptons[m_index_lep_a[i]];
        cand_list_parts_a.push_back(part);
      }
      for(int i = 0; i < m_Njet_a; i++){
        Particle part = Jets[m_index_jet_a[i]];
        cand_list_parts_a.push_back(part);
      }
      L_Cand cand_a(cand_list_parts_a);

      ParticleList cand_list_parts_b;
      for(int i = 0; i < m_Nlep_b; i++){
        Particle part = Leptons[m_index_lep_b[i]];
        cand_list_parts_b.push_back(part);
      }
      for(int i = 0; i < m_Njet_b; i++){
        Particle part = Jets[m_index_jet_b[i]];
        cand_list_parts_b.push_back(part);
      }
      L_Cand cand_b(cand_list_parts_b);

      if(cand_a.GetN() + cand_b.GetN() < 3){
        ParticleList dum_cand_list_parts = cand_a.PL() + cand_b.PL();
        cand_a = L_Cand(dum_cand_list_parts);
      }
      
      if(cand_a.GetN() > 1)
        m_CosDecayAngle_Va = cand_a.CosDecayAngle();
      else
        m_CosDecayAngle_Va = 0.;
      if(cand_b.GetN() > 1)
        m_CosDecayAngle_Vb = cand_b.CosDecayAngle();
      else
        m_CosDecayAngle_Vb = 0.;

    }
    if(t==1){
      if(m_Nlep < 2){
        m_treeSkipped[t] = true;
        continue;
      }
      std::vector<RFKey> lepID;
      for(int i = 0; i < m_Nlep; i++){
        lepID.push_back(COMB_L[t]->AddLabFrameFourVector(Leptons[i]));
      }
      TLorentzVector TLV_ISR;
      for(int i = 0; i < m_Njet; i++){
        TLorentzVector jet; jet.SetPtEtaPhiM( Jets[i].Pt(),
                                Jets[i].Eta(),
                                Jets[i].Phi(),
                                std::max(0.,Jets[i].M()) );
        TLV_ISR += jet;
      }
      ISR[t]->SetLabFrameFourVector(TLV_ISR);
      if(!LAB[t]->AnalyzeEvent())
        cout << "Something went wrong with tree event analysis" << endl;
          
      bool BSideIsA = false;
      if(Lb[t]->GetFourVector().M() > La[t]->GetFourVector().M()) BSideIsA = true;

      for(int i = 0; i < m_Nlep; i++){
        if(COMB_L[t]->GetFrame(lepID[i]) == *La[t]){
          if(!BSideIsA){
            m_Nlep_a_LEP++;
            m_index_lep_a_LEP.push_back(i);
          }
          else{
            m_Nlep_b_LEP++;
            m_index_lep_b_LEP.push_back(i);
          }
        }
        if(COMB_L[t]->GetFrame(lepID[i]) == *Lb[t]){
          if(!BSideIsA){
            m_Nlep_b_LEP++;
            m_index_lep_b_LEP.push_back(i);
          }
          else{
            m_Nlep_a_LEP++;
            m_index_lep_a_LEP.push_back(i);
          }
        }
      }

      m_PTCM_LEP = CM[t]->GetFourVector().Pt();
      m_PzCM_LEP = CM[t]->GetFourVector().Pz();
      m_cosCM_LEP = CM[t]->GetCosDecayAngle();
      m_dphiCM_LEP = CM[t]->GetDeltaPhiDecayAngle();
      m_dphiCMI_LEP = CM[t]->GetDeltaPhiBoostVisible();
      m_dphiMET_V_LEP = S[t]->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

      TVector3 vPISR = S[t]->GetFourVector(*CM[t]).Vect();
      m_PTISR_LEP = vPISR.Pt();
      TVector3 vPINV = (X1a[t]->GetFourVector(*CM[t])+X1b[t]->GetFourVector(*CM[t])).Vect();
      m_RISR_LEP = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();

      m_MPa_LEP = X2a[t]->GetFourVector().M();
      m_MPb_LEP = X2b[t]->GetFourVector().M();
      if(BSideIsA){ m_MPa_LEP = X2b[t]->GetFourVector().M(); m_MPb_LEP = X2a[t]->GetFourVector().M(); }
      m_MS_LEP = S[t]->GetFourVector().M();
      m_MSV_LEP = (La[t]->GetFourVector()+Lb[t]->GetFourVector()).M();
      m_MQ_LEP = sqrt(m_MPa_LEP*m_MPa_LEP+m_MPb_LEP*m_MPb_LEP)/sqrt(2.);
      m_gamma_LEP = 2.*m_MQ_LEP/m_MS_LEP;
      m_CosDecayAngle_S_LEP = S[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pa_LEP = X2a[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pb_LEP = X2b[t]->GetCosDecayAngle();
      if(BSideIsA){ m_CosDecayAngle_Pa_LEP = X2b[t]->GetCosDecayAngle(); m_CosDecayAngle_Pb_LEP = X2a[t]->GetCosDecayAngle(); } 

      ParticleList cand_list_parts_a;
      for(int i = 0; i < m_Nlep_a_LEP; i++){
        Particle part = Leptons[m_index_lep_a_LEP[i]];
        cand_list_parts_a.push_back(part);
      }
      L_Cand cand_a(cand_list_parts_a);

      ParticleList cand_list_parts_b;
      for(int i = 0; i < m_Nlep_b_LEP; i++){
        Particle part = Leptons[m_index_lep_b_LEP[i]];
        cand_list_parts_b.push_back(part);
      }
      L_Cand cand_b(cand_list_parts_b);

      if(cand_a.GetN() + cand_b.GetN() < 3){
        ParticleList dum_cand_list_parts = cand_a.PL() + cand_b.PL();
        cand_a = L_Cand(dum_cand_list_parts);
      }
      
      if(cand_a.GetN() > 1)
        m_CosDecayAngle_Va_LEP = cand_a.CosDecayAngle();
      else
        m_CosDecayAngle_Va_LEP = 0.;
      if(cand_b.GetN() > 1)
        m_CosDecayAngle_Vb_LEP = cand_b.CosDecayAngle();
      else
        m_CosDecayAngle_Vb_LEP = 0.;

      m_MVa_LEP = La[t]->GetFourVector().M();
      m_MVb_LEP = Lb[t]->GetFourVector().M();
      if(BSideIsA){ m_MVa = Lb[t]->GetFourVector().M(); m_MVb = La[t]->GetFourVector().M(); }
      m_PTS_CM_LEP = S[t]->GetFourVector(*CM[t]).Pt();
      TLorentzVector TLV_L_CMLEP = La[t]->GetFourVector(*CM[t]) + Lb[t]->GetFourVector(*CM[t]);
      m_RZPara_LEP = TLV_L_CMLEP.Vect().Dot(S[t]->GetFourVector(*CM[t]).Vect().Unit())/S[t]->GetFourVector(*CM[t]).Vect().Mag();
      m_MINV_LEP = TLV_L_CMLEP.M()*m_RISR_LEP/m_RZPara_LEP;

      TLorentzVector X1a_S = X1a[t]->GetFourVector(*S[t]);
      TLorentzVector X1b_S = X1b[t]->GetFourVector(*S[t]);
      TLorentzVector La_S = La[t]->GetFourVector(*S[t]);
      TLorentzVector Lb_S = Lb[t]->GetFourVector(*S[t]);
      TLorentzVector Ia_S0;
      TLorentzVector Ib_S0;
      TLorentzVector La_S0;
      TLorentzVector Lb_S0;
      Ia_S0.SetPtEtaPhiM(X1a_S.Pt(), X1a_S.Eta(), X1a_S.Phi(), 0.);
      Ib_S0.SetPtEtaPhiM(X1b_S.Pt(), X1b_S.Eta(), X1b_S.Phi(), 0.);
      La_S0.SetPtEtaPhiM(La_S.Pt(), La_S.Eta(), La_S.Phi(), 0.);
      Lb_S0.SetPtEtaPhiM(Lb_S.Pt(), Lb_S.Eta(), Lb_S.Phi(), 0.);
      m_MS_S0_LEP = (Ia_S0+Ib_S0+La_S0+Lb_S0).M(); // this is Mr or 'ass mass' in the LLP analysis
      m_MV_S0_LEP = (La_S0+Lb_S0).M();
      m_MQ_S0_LEP = sqrt(((Ia_S0+La_S0).M2()+(Ib_S0+Lb_S0).M2())/2.);
      m_gamma_S0_LEP = 2.*m_MQ_S0_LEP/m_MS_S0_LEP; // this is 'R' in the LLP analysis
      m_MQV_LEP = sqrt(((La[t]->GetMass()*La[t]->GetMass())+(Lb[t]->GetMass()*Lb[t]->GetMass()))/2.);
      m_gammaV_LEP = 2.*m_MQV_LEP/m_MS_S0_LEP; // this is 'Rv' in the LLP analysis

      TLorentzVector X1a_X2a = X1a[t]->GetFourVector(*X2a[t]);
      TLorentzVector X1b_X2b = X1b[t]->GetFourVector(*X2b[t]);
      TLorentzVector X2a_S = X2a[t]->GetFourVector(*S[t]);
      m_MPTilde_LEP = X1a_X2a.Vect().Mag() + X1b_X2b.Vect().Mag();
      m_MSTilde_LEP = sqrt(m_MPTilde_LEP*m_MPTilde_LEP+X2a_S.Vect().Mag2());
      m_gammaTilde_LEP = m_MPTilde_LEP/m_MSTilde_LEP;
    }
    if(t==2){ // JET + ISR tree
      if(m_Njet < 2){
        m_treeSkipped[t] = true;
        continue;
      }
      std::vector<RFKey> jetID;
      for(int i = 0; i < m_Njet; i++){
        jetID.push_back(COMB_J[t]->AddLabFrameFourVector(Jets[i]));
      }
      if(!LAB[t]->AnalyzeEvent())
        cout << "Something went wrong with tree event analysis" << endl;

      bool BSideIsA = false;
      if(Jb[t]->GetFourVector().M() > Ja[t]->GetFourVector().M()) BSideIsA = true;

      // jet counting in ISR/S, hemispheres
      for(int i = 0; i < m_Njet; i++){
        if(COMB_J[t]->GetFrame(jetID[i]) == *ISR[t]){
          m_Njet_ISR_JET_ISR++;
          if(Jets[i].BtagID() >= kMedium)
            m_Nbjet_ISR_JET_ISR++;
          m_index_jet_ISR_JET_ISR.push_back(i);
        }
        if(COMB_J[t]->GetFrame(jetID[i]) == *Ja[t]){
          if(!BSideIsA){
            m_Njet_S_JET_ISR++;
            m_Njet_a_JET_ISR++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S_JET_ISR++;
              m_Nbjet_a_JET_ISR++;
            }
            m_index_jet_S_JET_ISR.push_back(i);
            m_index_jet_a_JET_ISR.push_back(i);
          }
          else{
            m_Njet_S_JET_ISR++;
            m_Njet_b_JET_ISR++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S_JET_ISR++;
              m_Nbjet_b_JET_ISR++;
            }
            m_index_jet_S_JET_ISR.push_back(i);
            m_index_jet_b_JET_ISR.push_back(i);
          }
        }
        if(COMB_J[t]->GetFrame(jetID[i]) == *Jb[t]){
          if(!BSideIsA){
            m_Njet_S_JET_ISR++;
            m_Njet_b_JET_ISR++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S_JET_ISR++;
              m_Nbjet_b_JET_ISR++;
            }
            m_index_jet_S_JET_ISR.push_back(i);
            m_index_jet_b_JET_ISR.push_back(i);
          }
          else{
            m_Njet_S_JET_ISR++;
            m_Njet_a_JET_ISR++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_S_JET_ISR++;
              m_Nbjet_a_JET_ISR++;
            }
            m_index_jet_S_JET_ISR.push_back(i);
            m_index_jet_a_JET_ISR.push_back(i);
          }
        }
      }

      m_PTCM_JET_ISR = CM[t]->GetFourVector().Pt();
      m_PzCM_JET_ISR = CM[t]->GetFourVector().Pz();
      m_cosCM_JET_ISR = CM[t]->GetCosDecayAngle();
      m_dphiCM_JET_ISR = CM[t]->GetDeltaPhiDecayAngle();
      m_dphiCMI_JET_ISR = CM[t]->GetDeltaPhiBoostVisible();
      m_dphiMET_V_JET_ISR = S[t]->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);
          
      TVector3 vPISR = S[t]->GetFourVector(*CM[t]).Vect();
      m_PTISR_JET_ISR = vPISR.Pt();
      TVector3 vPINV = (X1a[t]->GetFourVector(*CM[t])+X1b[t]->GetFourVector(*CM[t])).Vect();
      m_RISR_JET_ISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();

      m_MPa_JET_ISR = X2a[t]->GetFourVector().M();
      m_MPb_JET_ISR = X2b[t]->GetFourVector().M();
      if(BSideIsA){ m_MPa_JET_ISR = X2b[t]->GetFourVector().M(); m_MPb_JET_ISR = X2a[t]->GetFourVector().M(); }
      m_MS_JET_ISR = S[t]->GetFourVector().M();
      m_MSV_JET_ISR = (Ja[t]->GetFourVector()+Jb[t]->GetFourVector()).M();
      m_MQ_JET_ISR = sqrt(m_MPa_JET_ISR*m_MPa_JET_ISR+m_MPb_JET_ISR*m_MPb_JET_ISR)/sqrt(2.);
      m_gamma_JET_ISR = 2.*m_MQ_JET_ISR/m_MS_JET_ISR;
      m_CosDecayAngle_S_JET_ISR = S[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pa_JET_ISR = X2a[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pb_JET_ISR = X2b[t]->GetCosDecayAngle();
      if(BSideIsA){ m_CosDecayAngle_Pa_JET_ISR = X2b[t]->GetCosDecayAngle(); m_CosDecayAngle_Pb_JET_ISR = X2a[t]->GetCosDecayAngle(); }

      ParticleList cand_list_parts_a;
      for(int i = 0; i < m_Njet_a_JET_ISR; i++){
        Particle part = Jets[m_index_jet_a_JET_ISR[i]];
        cand_list_parts_a.push_back(part);
      }
      L_Cand cand_a(cand_list_parts_a);

      ParticleList cand_list_parts_b;
      for(int i = 0; i < m_Njet_b_JET_ISR; i++){
        Particle part = Jets[m_index_jet_b_JET_ISR[i]];
        cand_list_parts_b.push_back(part);
      }
      L_Cand cand_b(cand_list_parts_b);

      if(cand_a.GetN() + cand_b.GetN() < 3){
        ParticleList dum_cand_list_parts = cand_a.PL() + cand_b.PL();
        cand_a = L_Cand(dum_cand_list_parts);
      }

      if(cand_a.GetN() > 1)
        m_CosDecayAngle_Va_JET_ISR = cand_a.CosDecayAngle();
      else
        m_CosDecayAngle_Va_JET_ISR = 0.;
      if(cand_b.GetN() > 1)
        m_CosDecayAngle_Vb_JET_ISR = cand_b.CosDecayAngle();
      else
        m_CosDecayAngle_Vb_JET_ISR = 0.;

      m_MVa_JET_ISR = Ja[t]->GetFourVector().M();
      m_MVb_JET_ISR = Jb[t]->GetFourVector().M();
      if(BSideIsA){ m_MVa = Jb[t]->GetFourVector().M(); m_MVb = Ja[t]->GetFourVector().M(); }
      m_PTS_CM_JET_ISR = S[t]->GetFourVector(*CM[t]).Pt();
      TLorentzVector TLV_L_CMJET_ISR = Ja[t]->GetFourVector(*CM[t]) + Jb[t]->GetFourVector(*CM[t]);
      m_RZPara_JET_ISR = TLV_L_CMJET_ISR.Vect().Dot(S[t]->GetFourVector(*CM[t]).Vect().Unit())/S[t]->GetFourVector(*CM[t]).Vect().Mag();
      m_MINV_JET_ISR = TLV_L_CMJET_ISR.M()*m_RISR_JET_ISR/m_RZPara_JET_ISR;

      TLorentzVector X1a_S = X1a[t]->GetFourVector(*S[t]);
      TLorentzVector X1b_S = X1b[t]->GetFourVector(*S[t]);
      TLorentzVector Ja_S = Ja[t]->GetFourVector(*S[t]);
      TLorentzVector Jb_S = Jb[t]->GetFourVector(*S[t]);
      TLorentzVector Ia_S0;
      TLorentzVector Ib_S0;
      TLorentzVector Ja_S0;
      TLorentzVector Jb_S0;
      Ia_S0.SetPtEtaPhiM(X1a_S.Pt(), X1a_S.Eta(), X1a_S.Phi(), 0.);
      Ib_S0.SetPtEtaPhiM(X1b_S.Pt(), X1b_S.Eta(), X1b_S.Phi(), 0.);
      Ja_S0.SetPtEtaPhiM(Ja_S.Pt(), Ja_S.Eta(), Ja_S.Phi(), 0.);
      Jb_S0.SetPtEtaPhiM(Jb_S.Pt(), Jb_S.Eta(), Jb_S.Phi(), 0.);
      m_MS_S0_JET_ISR = (Ia_S0+Ib_S0+Ja_S0+Jb_S0).M(); // this is Mr or 'ass mass' in the LLP analysis
      m_MV_S0_JET_ISR = (Ja_S0+Jb_S0).M();
      m_MQ_S0_JET_ISR = sqrt(((Ia_S0+Ja_S0).M2()+(Ib_S0+Jb_S0).M2())/2.);
      m_gamma_S0_JET_ISR = 2.*m_MQ_S0_JET_ISR/m_MS_S0_JET_ISR; // this is 'R' in the LLP analysis
      m_MQV_JET_ISR = sqrt(((Ja[t]->GetMass()*Ja[t]->GetMass())+(Jb[t]->GetMass()*Jb[t]->GetMass()))/2.);
      m_gammaV_JET_ISR = 2.*m_MQV_JET_ISR/m_MS_S0_JET_ISR; // this is 'Rv' in the LLP analysis

      TLorentzVector X1a_X2a = X1a[t]->GetFourVector(*X2a[t]);
      TLorentzVector X1b_X2b = X1b[t]->GetFourVector(*X2b[t]);
      TLorentzVector X2a_S = X2a[t]->GetFourVector(*S[t]);
      m_MPTilde_JET_ISR = X1a_X2a.Vect().Mag() + X1b_X2b.Vect().Mag();
      m_MSTilde_JET_ISR = sqrt(m_MPTilde_JET_ISR*m_MPTilde_JET_ISR+X2a_S.Vect().Mag2());
      m_gammaTilde_JET_ISR = m_MPTilde_JET_ISR/m_MSTilde_JET_ISR;

      // removing momentum components parallel to CM->S boost
      TLorentzVector vP_S_CM  = S[t]->GetFourVector(*CM[t]);
      TLorentzVector vP_Ja_S  = Ja[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Jb_S  = Jb[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Ia_S  = X1a[t]->GetFourVector(*S[t]);
      TLorentzVector vP_Ib_S  = X1b[t]->GetFourVector(*S[t]);
      
      TLorentzVector vP_Ja_X2a  = Ja[t]->GetFourVector(*X2a[t]);
      TLorentzVector vP_Jb_X2b  = Jb[t]->GetFourVector(*X2b[t]);
      TLorentzVector vP_Ia_X2a  = X1a[t]->GetFourVector(*X2a[t]);
      TLorentzVector vP_Ib_X2b  = X1b[t]->GetFourVector(*X2b[t]);

      TVector3 boostVis = (vP_Ja_S+vP_Jb_S).BoostVector();
      TVector3 boostInv = (vP_Ia_S+vP_Ib_S).BoostVector();
      TVector3 daBoost = vP_S_CM.Vect().Unit();
      
      if((std::isnan(boostInv.Mag()) || std::isnan(boostVis.Mag())))
        cout << "boost NAN " << boostInv.Mag() << " " << boostVis.Mag() << " Ja=" << vP_Ja_S.P() << " Jb=" << vP_Jb_S.P() << " Inva=" << vP_Ia_S.P() << " Invb=" << vP_Ib_S.P() << " MET=" << ETMiss.Mag() << endl;
      
      boostVis = (boostVis.Dot(daBoost))*daBoost;
      boostInv = (boostInv.Dot(daBoost))*daBoost;

      if((!std::isnan(boostVis.Mag())) &&
         (boostVis.Mag() < 1)){
        vP_Ja_S.Boost(-boostVis);
        vP_Jb_S.Boost(-boostVis);
      } else {
        vP_Ja_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ja_S.M()));
        vP_Jb_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Jb_S.M()));
      }
      if((!std::isnan(boostInv.Mag())) &&
         (boostInv.Mag() < 1)){
        vP_Ia_S.Boost(-boostInv);
        vP_Ib_S.Boost(-boostInv);
      } else {
        vP_Ia_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ia_S.M()));
        vP_Ib_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ib_S.M()));
      }
        
      double PX2_BoostT = (vP_Ja_S+vP_Ia_S).P();
      double MX2a_BoostT = (vP_Ja_S+vP_Ia_S).M();
      double MX2b_BoostT = (vP_Jb_S+vP_Ib_S).M();
      m_Mperp_JET_ISR = sqrt(MX2a_BoostT*MX2a_BoostT+MX2b_BoostT*MX2b_BoostT)/sqrt(2.);
      m_gammaT_JET_ISR = 2*m_Mperp_JET_ISR/(sqrt(MX2a_BoostT*MX2a_BoostT+PX2_BoostT*PX2_BoostT) +
            		sqrt(MX2b_BoostT*MX2b_BoostT+PX2_BoostT*PX2_BoostT));
    } // End JET_ISR
    if(t==3){ // JET only (no ISR) tree
      if(m_Njet < 2){
        m_treeSkipped[t] = true;
        continue;
      }
      std::vector<RFKey> jetID;
      for(int i = 0; i < m_Njet; i++){
        jetID.push_back(COMB_J[t]->AddLabFrameFourVector(Jets[i]));
      }
      if(!LAB[t]->AnalyzeEvent())
        cout << "Something went wrong with tree event analysis" << endl;

      bool BSideIsA = false;
      if(Jb[t]->GetFourVector().M() > Ja[t]->GetFourVector().M()) BSideIsA = true;

      // jet counting in ISR/S, hemispheres
      for(int i = 0; i < m_Njet; i++){
        if(COMB_J[t]->GetFrame(jetID[i]) == *Ja[t]){
          if(!BSideIsA){
            m_Njet_a_JET++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_a_JET++;
            }
            m_index_jet_a_JET.push_back(i);
          }
          else{
            m_Njet_b_JET++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_b_JET++;
            }
            m_index_jet_b_JET.push_back(i);
          }
        }
        if(COMB_J[t]->GetFrame(jetID[i]) == *Jb[t]){
          if(!BSideIsA){
            m_Njet_b_JET++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_b_JET++;
            }
            m_index_jet_b_JET.push_back(i);
          }
          else{
            m_Njet_a_JET++;
            if(Jets[i].BtagID() >= kMedium){
              m_Nbjet_a_JET++;
            }
            m_index_jet_a_JET.push_back(i);
          }
        }
      }

      // no 'CM' system since everything is in S
      m_PTCM_JET = S[t]->GetFourVector().Pt();
      m_PzCM_JET = S[t]->GetFourVector().Pz();
      m_cosCM_JET = S[t]->GetCosDecayAngle();
      m_dphiCM_JET = S[t]->GetDeltaPhiDecayAngle();
      m_dphiCMI_JET = S[t]->GetDeltaPhiBoostVisible();
      m_dphiMET_V_JET = S[t]->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

      m_MPa_JET = X2a[t]->GetFourVector().M();
      m_MPb_JET = X2b[t]->GetFourVector().M();
      if(BSideIsA){ m_MPa_JET = X2b[t]->GetFourVector().M(); m_MPb_JET = X2a[t]->GetFourVector().M(); }
      m_MS_JET = S[t]->GetFourVector().M();
      m_MSV_JET = (Ja[t]->GetFourVector()+Jb[t]->GetFourVector()).M();
      m_MQ_JET = sqrt(m_MPa_JET*m_MPa_JET+m_MPb_JET*m_MPb_JET)/sqrt(2.);
      m_gamma_JET = 2.*m_MQ_JET/m_MS_JET;
      m_CosDecayAngle_S_JET = S[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pa_JET = X2a[t]->GetCosDecayAngle();
      m_CosDecayAngle_Pb_JET = X2b[t]->GetCosDecayAngle();
      if(BSideIsA){ m_CosDecayAngle_Pa_JET = X2b[t]->GetCosDecayAngle(); m_CosDecayAngle_Pb_JET = X2a[t]->GetCosDecayAngle(); }

      ParticleList cand_list_parts_a;
      for(int i = 0; i < m_Njet_a_JET; i++){
        Particle part = Jets[m_index_jet_a_JET[i]];
        cand_list_parts_a.push_back(part);
      }
      L_Cand cand_a(cand_list_parts_a);

      ParticleList cand_list_parts_b;
      for(int i = 0; i < m_Njet_b_JET; i++){
        Particle part = Jets[m_index_jet_b_JET[i]];
        cand_list_parts_b.push_back(part);
      }
      L_Cand cand_b(cand_list_parts_b);

      if(cand_a.GetN() + cand_b.GetN() < 3){
        ParticleList dum_cand_list_parts = cand_a.PL() + cand_b.PL();
        cand_a = L_Cand(dum_cand_list_parts);
      }

      if(cand_a.GetN() > 1)
        m_CosDecayAngle_Va_JET = cand_a.CosDecayAngle();
      else
        m_CosDecayAngle_Va_JET = 0.;
      if(cand_b.GetN() > 1)
        m_CosDecayAngle_Vb_JET = cand_b.CosDecayAngle();
      else
        m_CosDecayAngle_Vb_JET = 0.;

      m_MVa_JET = Ja[t]->GetFourVector().M();
      m_MVb_JET = Jb[t]->GetFourVector().M();
      if(BSideIsA){ m_MVa = Jb[t]->GetFourVector().M(); m_MVb = Ja[t]->GetFourVector().M(); }

      TLorentzVector X1a_S = X1a[t]->GetFourVector(*S[t]);
      TLorentzVector X1b_S = X1b[t]->GetFourVector(*S[t]);
      TLorentzVector Ja_S = Ja[t]->GetFourVector(*S[t]);
      TLorentzVector Jb_S = Jb[t]->GetFourVector(*S[t]);
      TLorentzVector Ia_S0;
      TLorentzVector Ib_S0;
      TLorentzVector Ja_S0;
      TLorentzVector Jb_S0;
      Ia_S0.SetPtEtaPhiM(X1a_S.Pt(), X1a_S.Eta(), X1a_S.Phi(), 0.);
      Ib_S0.SetPtEtaPhiM(X1b_S.Pt(), X1b_S.Eta(), X1b_S.Phi(), 0.);
      Ja_S0.SetPtEtaPhiM(Ja_S.Pt(), Ja_S.Eta(), Ja_S.Phi(), 0.);
      Jb_S0.SetPtEtaPhiM(Jb_S.Pt(), Jb_S.Eta(), Jb_S.Phi(), 0.);
      m_MS_S0_JET = (Ia_S0+Ib_S0+Ja_S0+Jb_S0).M(); // this is Mr or 'ass mass' in the LLP analysis
      m_MV_S0_JET = (Ja_S0+Jb_S0).M();
      m_MQ_S0_JET = sqrt(((Ia_S0+Ja_S0).M2()+(Ib_S0+Jb_S0).M2())/2.);
      m_gamma_S0_JET = 2.*m_MQ_S0_JET/m_MS_S0_JET; // this is 'R' in the LLP analysis
      m_MQV_JET = sqrt(((Ja[t]->GetMass()*Ja[t]->GetMass())+(Jb[t]->GetMass()*Jb[t]->GetMass()))/2.);
      m_gammaV_JET = 2.*m_MQV_JET/m_MS_S0_JET; // this is 'Rv' in the LLP analysis

      TLorentzVector X1a_X2a = X1a[t]->GetFourVector(*X2a[t]);
      TLorentzVector X1b_X2b = X1b[t]->GetFourVector(*X2b[t]);
      TLorentzVector X2a_S = X2a[t]->GetFourVector(*S[t]);
      m_MPTilde_JET = X1a_X2a.Vect().Mag() + X1b_X2b.Vect().Mag();
      m_MSTilde_JET = sqrt(m_MPTilde_JET*m_MPTilde_JET+X2a_S.Vect().Mag2());
      m_gammaTilde_JET = m_MPTilde_JET/m_MSTilde_JET;

    } // End JET
  } // End of RJR trees analysis

  if(!AnalysisBase<Base>::IsData()){
    m_weight = AnalysisBase<Base>::GetEventWeight();
    m_weight2 = m_weight*m_weight;
    //if(m_weight == 0.) { cout << "Event Weight: " << tree->GetName() << " " << m_weight << endl; } // debug statement
    m_genweight = AnalysisBase<Base>::GetGenEventWeight();
    m_XSec = AnalysisBase<Base>::GetXsec();
    m_FilterEff = AnalysisBase<Base>::GetFilterEff();
    m_Nevent = AnalysisBase<Base>::GetNevent();
    m_Nweight = AnalysisBase<Base>::GetNweight();
    m_Npartons = AnalysisBase<Base>::GetNpartons();
    
    m_PUweight = AnalysisBase<Base>::GetPUWeight(0);
    m_PUweight_up = AnalysisBase<Base>::GetPUWeight(1);
    m_PUweight_down = AnalysisBase<Base>::GetPUWeight(-1);

    m_MuFweight = AnalysisBase<Base>::GetMuFWeight(0);
    m_MuFweight_up = AnalysisBase<Base>::GetMuFWeight(1);
    m_MuFweight_down = AnalysisBase<Base>::GetMuFWeight(-1);

    m_MuRweight = AnalysisBase<Base>::GetMuRWeight(0);
    m_MuRweight_up = AnalysisBase<Base>::GetMuRWeight(1);
    m_MuRweight_down = AnalysisBase<Base>::GetMuRWeight(-1);

    m_PDFweight = AnalysisBase<Base>::GetPDFWeight(0);
    m_PDFweight_up = AnalysisBase<Base>::GetPDFWeight(1);
    m_PDFweight_down = AnalysisBase<Base>::GetPDFWeight(-1);
    
    m_BtagHFSFweight = AnalysisBase<Base>::GetBtagSFWeight(Jets, true, 0, kMedium);
    m_BtagHFSFweight_up = AnalysisBase<Base>::GetBtagSFWeight(Jets, true, 1, kMedium);
    m_BtagHFSFweight_down = AnalysisBase<Base>::GetBtagSFWeight(Jets, true, -1, kMedium);
    m_BtagLFSFweight = AnalysisBase<Base>::GetBtagSFWeight(Jets, false, 0, kMedium);
    m_BtagLFSFweight_up = AnalysisBase<Base>::GetBtagSFWeight(Jets, false, 1, kMedium);
    m_BtagLFSFweight_down = AnalysisBase<Base>::GetBtagSFWeight(Jets, false, -1, kMedium);

    m_elIDSFweight = AnalysisBase<Base>::GetElIDSFWeight(Electrons, 0);
    m_elIDSFweight_up = AnalysisBase<Base>::GetElIDSFWeight(Electrons, 1);
    m_elIDSFweight_down = AnalysisBase<Base>::GetElIDSFWeight(Electrons, -1);
    m_elISOSFweight = AnalysisBase<Base>::GetElISOSFWeight(Electrons, 0);
    m_elISOSFweight_up = AnalysisBase<Base>::GetElISOSFWeight(Electrons, 1);
    m_elISOSFweight_down = AnalysisBase<Base>::GetElISOSFWeight(Electrons, -1);
    m_elSIPSFweight = AnalysisBase<Base>::GetElSIPSFWeight(Electrons, 0);
    m_elSIPSFweight_up = AnalysisBase<Base>::GetElSIPSFWeight(Electrons, 1);
    m_elSIPSFweight_down = AnalysisBase<Base>::GetElSIPSFWeight(Electrons, -1);
    m_elVLSFweight = AnalysisBase<Base>::GetElVLIDSFWeight(Electrons, 0);
    m_elVLSFweight_up = AnalysisBase<Base>::GetElVLIDSFWeight(Electrons, 1);
    m_elVLSFweight_down = AnalysisBase<Base>::GetElVLIDSFWeight(Electrons, -1);
    m_muIDSFweight = AnalysisBase<Base>::GetMuIDSFWeight(Muons, 0);
    m_muIDSFweight_up = AnalysisBase<Base>::GetMuIDSFWeight(Muons, 1);
    m_muIDSFweight_down = AnalysisBase<Base>::GetMuIDSFWeight(Muons, -1);
    m_muISOSFweight = AnalysisBase<Base>::GetMuISOSFWeight(Muons, 0);
    m_muISOSFweight_up = AnalysisBase<Base>::GetMuISOSFWeight(Muons, 1);
    m_muISOSFweight_down = AnalysisBase<Base>::GetMuISOSFWeight(Muons, -1);
    m_muSIPSFweight = AnalysisBase<Base>::GetMuSIPSFWeight(Muons, 0);
    m_muSIPSFweight_up = AnalysisBase<Base>::GetMuSIPSFWeight(Muons, 1);
    m_muSIPSFweight_down = AnalysisBase<Base>::GetMuSIPSFWeight(Muons, -1);
    m_muVLSFweight = AnalysisBase<Base>::GetMuVLIDSFWeight(Muons, 0);
    m_muVLSFweight_up = AnalysisBase<Base>::GetMuVLIDSFWeight(Muons, 1);
    m_muVLSFweight_down = AnalysisBase<Base>::GetMuVLIDSFWeight(Muons, -1);

    m_MetTrigSFweight = AnalysisBase<Base>::GetMETTriggerSFWeight(m_MET, m_HT_eta5, m_Nele, m_Nmu, 0);
    m_MetTrigSFweight_up = AnalysisBase<Base>::GetMETTriggerSFWeight(m_MET, m_HT_eta5, m_Nele, m_Nmu, 1);
    m_MetTrigSFweight_down = AnalysisBase<Base>::GetMETTriggerSFWeight(m_MET, m_HT_eta5, m_Nele, m_Nmu, -1);
    m_MetTrigSFCurveIndex = AnalysisBase<Base>::GetMETTriggerSFCurve(m_HT_eta5, m_Nele, m_Nmu);
   
    m_NPU = AnalysisBase<Base>::GetNPUtrue();

    TVector3 genETMiss = AnalysisBase<Base>::GetGenMET();
    m_genMET     = genETMiss.Pt();
    m_genMET_phi = genETMiss.Phi();
  } else {
    m_weight = 1;
    m_PUweight = 1;
    m_PUweight_up = 1;
    m_PUweight_down = 1;
    m_MuFweight = 1;
    m_MuFweight_up = 1;
    m_MuFweight_down = 1;
    m_MuRweight = 1;
    m_MuRweight_up = 1;
    m_MuRweight_down = 1;
    m_PDFweight = 1;
    m_PDFweight_up = 1;
    m_PDFweight_down = 1;
    m_BtagHFSFweight = 1;
    m_BtagHFSFweight_up = 1;
    m_BtagHFSFweight_down = 1;
    m_BtagLFSFweight = 1;
    m_BtagLFSFweight_up = 1;
    m_BtagLFSFweight_down = 1;
    m_elIDSFweight = 1;
    m_elIDSFweight_up = 1;
    m_elIDSFweight_down = 1;
    m_elISOSFweight = 1;
    m_elISOSFweight_up = 1;
    m_elISOSFweight_down = 1;
    m_elSIPSFweight = 1;
    m_elSIPSFweight_up = 1;
    m_elSIPSFweight_down = 1;
    m_elVLSFweight = 1;
    m_elVLSFweight_up = 1;
    m_elVLSFweight_down = 1;
    m_muIDSFweight = 1;
    m_muIDSFweight_up = 1;
    m_muIDSFweight_down = 1;
    m_muISOSFweight = 1;
    m_muISOSFweight_up = 1;
    m_muISOSFweight_down = 1;
    m_muSIPSFweight = 1;
    m_muSIPSFweight_up = 1;
    m_muSIPSFweight_down = 1;
    m_muVLSFweight = 1;
    m_muVLSFweight_up = 1;
    m_muVLSFweight_down = 1;
    m_MetTrigSFweight = 1.;
    m_MetTrigSFweight_up = 1.;
    m_MetTrigSFweight_down = 1.;
    m_MetTrigSFCurveIndex = 0;
    m_NPU = 0.;
    m_genweight = 1;
    m_cascades_tree = 0;
    m_cascades_prod = -1;
    m_cascades_SlepSneu_1stDecay = -1;
    m_cascades_SlepSneu_2ndDecay = -1;
    m_cascades_N2_1stDecay = -1;
    m_cascades_N2_2ndDecay = -1;
    m_MSlepL = 0;
    m_MSneu = 0;
    m_MN2 = 0;
    m_MC1 = 0;
    m_MN1 = 0;
    m_MP = 0;
    m_NSparticleW = 0;
    m_NSparticleZ = 0;
    m_NSparticlePhoton = 0;
    m_Npartons = 0;
    m_XSec = 0;
    m_FilterEff = 1;
    m_Nweight = 0;
    m_Nevent = 0;
  }
  
  m_runnum   = AnalysisBase<Base>::GetRunNum();
  m_luminum  = AnalysisBase<Base>::GetLumiNum();
  m_eventnum = AnalysisBase<Base>::GetEventNum();

  m_NPV = AnalysisBase<Base>::GetNPV();

  m_METtrigger   = AnalysisBase<Base>::GetMETtrigger();
  m_METORtrigger = AnalysisBase<Base>::GetMETORtrigger();

  m_SingleElectrontrigger = AnalysisBase<Base>::GetSingleElectrontrigger();
  m_SingleMuontrigger = AnalysisBase<Base>::GetSingleMuontrigger();
  m_DoubleElectrontrigger = AnalysisBase<Base>::GetDoubleElectrontrigger();
  m_DoubleMuontrigger = AnalysisBase<Base>::GetDoubleMuontrigger();
  m_TripleElectrontrigger = AnalysisBase<Base>::GetTripleElectrontrigger();
  m_TripleMuonLowPTtrigger = AnalysisBase<Base>::GetTripleMuonLowPTtrigger();
  m_TripleMuonHighPTtrigger = AnalysisBase<Base>::GetTripleMuonHighPTtrigger();
  m_DiMuEleLowPTtrigger = AnalysisBase<Base>::GetDiMuEleLowPTtrigger();
  m_DiMuEleHighPTtrigger = AnalysisBase<Base>::GetDiMuEleHighPTtrigger();
  m_DiEleMutrigger = AnalysisBase<Base>::GetDiEleMutrigger();
  m_EMutrigger = AnalysisBase<Base>::GetEMutrigger(); 
  m_EMuMutrigger = AnalysisBase<Base>::GetEMuMutrigger(); 
  m_EMuEtrigger = AnalysisBase<Base>::GetEMuEtrigger(); 
  
  m_MET     = ETMiss.Pt();
  m_MET_phi = ETMiss.Phi();

  // Fill Jets
  if(!do_slim){
    m_PT_jet.clear();
    m_Eta_jet.clear();
    m_Phi_jet.clear();
    m_M_jet.clear();
    m_Btag_jet.clear();
    m_BtagID_jet.clear();
    m_Flavor_jet.clear();
    for(int i = 0; i < m_Njet; i++){
      m_PT_jet.push_back(Jets[i].Pt());
      m_Eta_jet.push_back(Jets[i].Eta());
      m_Phi_jet.push_back(Jets[i].Phi());
      m_M_jet.push_back(Jets[i].M());
      m_Btag_jet.push_back(Jets[i].Btag());
      m_BtagID_jet.push_back(Jets[i].BtagID());
      m_Flavor_jet.push_back(Jets[i].PDGID());
    }

    // Fill GenJets
    vector<int> genmatch_jet;
    for(int i = 0; i < m_NGenjet; i++)
      genmatch_jet.push_back(-1);
    m_PT_Genjet.clear();
    m_Eta_Genjet.clear();
    m_Phi_Genjet.clear();
    m_M_Genjet.clear();
    m_Index_jet.clear();
    for(int i = 0; i < m_NGenjet; i++){
      m_PT_Genjet.push_back(GenJets[i].Pt());
      m_Eta_Genjet.push_back(GenJets[i].Eta());
      m_Phi_Genjet.push_back(GenJets[i].Phi());
      m_M_Genjet.push_back(GenJets[i].M());
      int index = -1;
      double minDR = 0.1;
      for(int g = 0; g < m_Njet; g++)
	if(Jets[g].DeltaR(GenJets[i]) < minDR){
	  minDR = Jets[g].DeltaR(GenJets[i]);
	  index = i;
	  genmatch_jet[i] = g;
	}
      m_Index_jet.push_back(index);
    }
  }

  // Fill reconstructed lepton branches
  m_PT_lep.clear();
  m_Eta_lep.clear();
  m_Phi_lep.clear();
  m_M_lep.clear();
  m_Charge_lep.clear();
  m_PDGID_lep.clear();
  m_RelIso_lep.clear();
  m_MiniIso_lep.clear();
  m_Dxy_lep.clear();
  m_DxyErr_lep.clear();
  m_Dz_lep.clear();
  m_DzErr_lep.clear();
  m_IP3D_lep.clear();
  m_SIP3D_lep.clear();
  m_ID_lep.clear();
  m_SourceID_lep.clear();
  m_LepQual_lep.clear();
  m_IsLowPt_lep.clear();
  m_TightCharge_lep.clear();
  m_Index_lep.clear();
  vector<int> genmatch;
  for(int i = 0; i < m_genNlep; i++)
    genmatch.push_back(-1);
  for(int r = 0; r < m_Nlep; r++){
    m_PT_lep.push_back(Leptons[r].Pt());
    m_Eta_lep.push_back(Leptons[r].Eta());
    m_Phi_lep.push_back(Leptons[r].Phi());
    m_M_lep.push_back(Leptons[r].M());
    m_Charge_lep.push_back(Leptons[r].Charge());
    m_PDGID_lep.push_back(Leptons[r].PDGID());
    m_RelIso_lep.push_back(Leptons[r].RelIso());
    m_MiniIso_lep.push_back(Leptons[r].MiniIso());
    m_Dxy_lep.push_back(Leptons[r].Dxy());
    m_DxyErr_lep.push_back(Leptons[r].DxyErr());
    m_Dz_lep.push_back(Leptons[r].Dz());
    m_DzErr_lep.push_back(Leptons[r].DzErr());
    m_IP3D_lep.push_back(Leptons[r].IP3D());
    m_SIP3D_lep.push_back(Leptons[r].SIP3D());
    m_ID_lep.push_back(Leptons[r].ParticleID());
    m_LepQual_lep.push_back(Leptons[r].LepQual());
    m_IsLowPt_lep.push_back(Leptons[r].IsLowPt());
    m_TightCharge_lep.push_back(Leptons[r].TightCharge());
    int index = -1;
    double minDR = 0.1;
    for(int g = 0; g < m_genNlep; g++)
      if(Leptons[r].DeltaR(GenLeptons[g]) < minDR){
	minDR = Leptons[r].DeltaR(GenLeptons[g]);
	index = g;
	genmatch[g] = r;
      }
 
    m_Index_lep.push_back(index);
    if(index >= 0)
      Leptons[r].SetSourceID(GenLeptons[index].SourceID());
    else
      Leptons[r].SetSourceID(kFake);
    
    // m_ID_lep.push_back(Leptons[r].SourceID());
    m_SourceID_lep.push_back(Leptons[r].SourceID());
  }

  if(!AnalysisBase<Base>::IsData()){
    // Fill gen lepton branches
    m_genPT_lep.clear();
    m_genEta_lep.clear();
    m_genPhi_lep.clear();
    m_genM_lep.clear();
    m_genCharge_lep.clear();
    m_genPDGID_lep.clear();
    m_genMomPDGID_lep.clear();
    m_genSourceID_lep.clear();
    m_genIndex_lep.clear();
    m_genMomIndex_lep.clear();
    for(int g = 0; g < m_genNlep; g++){
      m_genPT_lep.push_back(GenLeptons[g].Pt());
      m_genEta_lep.push_back(GenLeptons[g].Eta());
      m_genPhi_lep.push_back(GenLeptons[g].Phi());
      m_genM_lep.push_back(GenLeptons[g].M());
      m_genCharge_lep.push_back(GenLeptons[g].Charge());
      m_genPDGID_lep.push_back(GenLeptons[g].PDGID());
      m_genMomPDGID_lep.push_back(GenLeptons[g].MomPDGID());
      m_genSourceID_lep.push_back(GenLeptons[g].SourceID());
      m_genIndex_lep.push_back(genmatch[g]);
      m_genMomIndex_lep.push_back(GenLeptons[g].GenMomIndex());
    }
  
    // Fill gen neutrino branches
    ParticleList GenNus = AnalysisBase<Base>::GetGenNeutrinos();
    m_genNnu = GenNus.size();
    m_genPT_nu.clear();
    m_genEta_nu.clear();
    m_genPhi_nu.clear();
    m_genPDGID_nu.clear();
    m_genMomPDGID_nu.clear();
    for(int i = 0; i < m_genNnu; i++){
      m_genPT_nu.push_back(GenNus[i].Pt());
      m_genEta_nu.push_back(GenNus[i].Eta());
      m_genPhi_nu.push_back(GenNus[i].Phi());
      m_genPDGID_nu.push_back(GenNus[i].PDGID());
      m_genMomPDGID_nu.push_back(GenNus[i].MomPDGID());
    }
  
    // Fill gen boson branches
    ParticleList GenBosons = AnalysisBase<Base>::GetGenBosons();
    m_genNboson = GenBosons.size();
    m_genPT_boson.clear();
    m_genEta_boson.clear();
    m_genPhi_boson.clear();
    m_genM_boson.clear();
    m_genPDGID_boson.clear();
    m_genMomPDGID_boson.clear();
    for(int i = 0; i < m_genNboson; i++){
      m_genPT_boson.push_back(GenBosons[i].Pt());
      m_genEta_boson.push_back(GenBosons[i].Eta());
      m_genPhi_boson.push_back(GenBosons[i].Phi());
      m_genM_boson.push_back(GenBosons[i].M());
      m_genPDGID_boson.push_back(GenBosons[i].PDGID());
      m_genMomPDGID_boson.push_back(GenBosons[i].MomPDGID());
    }

    // Fill gen sparticle branches
    ParticleList GenSparticles = AnalysisBase<Base>::GetGenSparticles();
    m_genNsusy = GenSparticles.size();
    m_genPT_susy.clear();
    m_genEta_susy.clear();
    m_genPhi_susy.clear();
    m_genM_susy.clear();
    m_genPDGID_susy.clear();
    m_genMomPDGID_susy.clear();
    for(int i = 0; i < m_genNsusy; i++){
      m_genPT_susy.push_back(GenSparticles[i].Pt());
      m_genEta_susy.push_back(GenSparticles[i].Eta());
      m_genPhi_susy.push_back(GenSparticles[i].Phi());
      m_genM_susy.push_back(GenSparticles[i].M());
      m_genPDGID_susy.push_back(GenSparticles[i].PDGID());
      m_genMomPDGID_susy.push_back(GenSparticles[i].MomPDGID());
    }
  }

  if(AnalysisBase<Base>::IsCascades()){
    m_cascades_tree = AnalysisBase<Base>::GetGenCascadesTree();
    std::tie(
     m_cascades_prod,
     m_cascades_SlepSneu_1stDecay,
     m_cascades_SlepSneu_2ndDecay,
     m_cascades_N2_1stDecay,
     m_cascades_N2_2ndDecay
    ) = CascadesTreeEncoder::Decode(m_cascades_tree);
  }

  if(AnalysisBase<Base>::IsCascades() || AnalysisBase<Base>::IsSMS()){
    if(m_cascades_prod > 0 && m_cascades_prod < 5){
      m_MSlepL = AnalysisBase<Base>::GetGenMass(1000011);
      m_MSneu  = AnalysisBase<Base>::GetGenMass(1000012);
    }
    else{
      m_MSlepL = AnalysisBase<Base>::GetGenMass(1000013);
      m_MSneu  = AnalysisBase<Base>::GetGenMass(1000014);
    }
    m_MN2 = AnalysisBase<Base>::GetGenMass(1000023);
    m_MC1 = AnalysisBase<Base>::GetGenMass(1000024);
    m_NSparticleW = AnalysisBase<Base>::GetGenSUSYNBosons(24);
    m_NSparticleZ = AnalysisBase<Base>::GetGenSUSYNBosons(23);
    m_NSparticlePhoton = AnalysisBase<Base>::GetGenSUSYNBosons(22);
    m_LSPParents = AnalysisBase<Base>::GetLSPParents();
  }
  
  // Fill output tree
  if(tree)
    tree->Fill();
}

template class ReducedNtuple<SUSYNANOBase>;
template class ReducedNtuple<NANOULBase>;
template class ReducedNtuple<NANORun3>;
