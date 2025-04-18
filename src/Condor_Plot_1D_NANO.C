#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TColorWheel.h>
#include <TH1.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TError.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include "SUSYNANOBase.hh"
#include "CategoryTool.hh"
#include "ParticleList.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

string load_file_from_list(string fileListName, int line_number);
string double_to_string(double val);

int main(int argc, char* argv[]) {

  string plot_folder = "";
  bool use_gen_jets = false;
  bool boson_acceptance_cut = false;
  bool gen_lepton_cut = false;
  int min_lep_cut = -1;
  int max_lep_cut = -1;
  int min_Sjet = -1;
  int max_Sjet = -1;
  int min_La = -1;
  int max_La = -1;
  int min_Lb = -1;
  int max_Lb = -1;
  double MET_cut = 0.;
  double PTISR_cut = 0.;
  double RISR_cut = 0.;

  bool bprint = false;
  bool debug = false;
  string ofile = "";
  string ifile = "";
  string proc = "";
  long double weight = 1.;
  int ichunk = 1;
  int nchunk = 1;
  bool treeplot = false;
  bool treeplot_INV = true; // whether to flip treeplots to dark mode
  string era = "Run2";

  if(argc == 0) bprint = true;
  for(int i = 0; i < argc; i++){
    if(strncmp(argv[i],"--help", 6) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-h", 2) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"--debug", 7) == 0){
      debug = true;
      cout << "turned on debug mode" << endl;
    }
    if(strncmp(argv[i],"--ofile", 7) == 0){
      i++;
      ofile = string(argv[i]);
      if(debug) cout << "ofile: " << ofile << endl;
    }
    if(strncmp(argv[i],"--ifile", 7) == 0){
      i++;
      ifile = string(argv[i]);
      if(debug) cout << "ifile: " << ifile << endl;
    }
    if(strncmp(argv[i],"--proc", 6) == 0){
      i++;
      proc = string(argv[i]);
      if(debug) cout << "proc: " << proc << endl;
    }
    if(strncmp(argv[i],"--treeplot", 10) == 0){
      treeplot = true;
        if(!gSystem->AccessPathName("trees.root"))
          gSystem->Exec("rm trees.root");
    }
    if(strncmp(argv[i],"--weight", 8) == 0){
      i++;
      weight = std::atoll(argv[i]);
    }
    if(strncmp(argv[i],"--split",7) == 0){
      sscanf(argv[i],"--split=%d,%d", &ichunk, &nchunk);
    }
    if(strncmp(argv[i],"--era", 5) == 0){
      i++;
      era = string(argv[i]);
    }
    if(strncmp(argv[i],"--genjets", 9) == 0){
      use_gen_jets = true;
      if(debug) cout << "using gen jets" << endl;
    }
    if(strncmp(argv[i],"--genlep", 8) == 0){
      gen_lepton_cut = true;
      if(debug) cout << "requiring one gen lep" << endl;
    }
    if(strncmp(argv[i],"--MET", 5) == 0){
      i++;
      MET_cut = std::atof(argv[i]);
      if(debug) cout << "using MET cut at: " << MET_cut << endl;
    }
    if(strncmp(argv[i],"--PTISR", 7) == 0){
      i++;
      PTISR_cut = std::atof(argv[i]);
      if(debug) cout << "using PTISR cut at: " << PTISR_cut << endl;
    }
    if(strncmp(argv[i],"--RISR", 6) == 0){
      i++;
      RISR_cut = std::atof(argv[i]);
      if(debug) cout << "using RISR cut at: " << RISR_cut << endl;
    }
    if(strncmp(argv[i],"--min_lep_cut", 13) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &min_lep_cut);
      if(debug) cout << "Requiring min " << min_lep_cut << " reco leps" << endl;
    }
    if(strncmp(argv[i],"--max_lep_cut", 13) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &max_lep_cut);
      if(debug) cout << "Requiring max " << max_lep_cut << " reco leps" << endl;
    }
    if(strncmp(argv[i],"--min_Sjet", 10) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &min_Sjet);
      if(debug) cout << "Requiring min " << min_Sjet << " jets in S" << endl;
    }
    if(strncmp(argv[i],"--max_Sjet", 10) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &max_Sjet);
      if(debug) cout << "Requiring max " << max_Sjet << " jets in S" << endl;
    }
    if(strncmp(argv[i],"--min_La", 8) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &min_La);
      if(debug) cout << "Requiring min " << min_La << " leps in Sa" << endl;
    }
    if(strncmp(argv[i],"--max_La", 8) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &max_La);
      if(debug) cout << "Requiring max " << max_La << " leps in Sa" << endl;
    }
    if(strncmp(argv[i],"--min_Lb", 8) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &min_Lb);
      if(debug) cout << "Requiring min " << min_Lb << " leps in Sb" << endl;
    }
    if(strncmp(argv[i],"--max_Lb", 8) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &max_Lb);
      if(debug) cout << "Requiring max " << max_Lb << " leps in Sb" << endl;
    }
  } // end arg loop

  if(bprint){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)          print options" << endl;
    cout << "Example: ./Condor_Plot_1D_NANO.x --proc ttbar --weight 10000000.000000 --ifile root://cmsxrootd.fnal.gov//store/mc/RunIIFall17NanoAODv7/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/PU2017_12Apr2018_Nano02Apr2020_102X_mc2017_realistic_v8-v1/100000/F66A024B-54F2-B340-BECD-9EFB7DFCF723.root --ofile Fall17_102X_TTJets_TuneCP5_13TeV-madgraphMLM-pythia8_ttbar_0_0.root --split=1,10 --era Run2" << endl;
    return 0;
  }

  RestFrames::SetStyle();
  TreePlot tree_plot("TreePlot","TreePlot");

  plot_folder += "plots";
  if(use_gen_jets) plot_folder += "_genjets";
  else plot_folder += "_recojets";
  if(boson_acceptance_cut) plot_folder += "_qqbar";
  else plot_folder += "_XX";
  if(gen_lepton_cut) plot_folder += "_l";
  else plot_folder += "_XX";

  // S ISR RJR tree
  LabRecoFrame LAB("LAB","lab");
  DecayRecoFrame CM("CM","CM");
  DecayRecoFrame S("S","S");
  DecayRecoFrame X2a("X2a","X2a");
  DecayRecoFrame X2b("X2b","X2b");
  VisibleRecoFrame J_X2a("J_X2a","J_{a}");
  VisibleRecoFrame J_X2b("J_X2b","J_{b}");
  VisibleRecoFrame L_X2a("L_X2a","L_{a}");
  VisibleRecoFrame L_X2b("L_X2b","L_{b}");
  VisibleRecoFrame ISR("ISR","ISR");
  InvisibleRecoFrame Ia("Ia","Ia");
  InvisibleRecoFrame Ib("Ib","Ib");

  LAB.SetChildFrame(CM);
  CM.AddChildFrame(S);
  CM.AddChildFrame(ISR);
  S.AddChildFrame(X2a);
  S.AddChildFrame(X2b);
  X2a.AddChildFrame(J_X2a);
  X2b.AddChildFrame(J_X2b);
  X2a.AddChildFrame(L_X2a);
  X2a.AddChildFrame(Ia);
  X2b.AddChildFrame(L_X2b);
  X2b.AddChildFrame(Ib);

  LAB.InitializeTree();

  CombinatoricGroup COMB_J("COMB_J", "Combinatoric System of Jets");
  MinMassesSqCombJigsaw CombSplitSq_J_Transverse("CombSplitSq_ISR", "Minimize M_{T}_{ISR}^{2} + M_{T}_{S}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_JS("CombSplitSq_S", "Minimize M_{X2a}^{2} + M_{X2b}^{2}",2,2);
  
  CombinatoricGroup COMB_L("COMB_L", "Combinatoric System of Leps");
  MinMassesSqCombJigsaw CombSplitSq_L("CombSplitSq_L", "Minimize M_{Va}^{2} + M_{Vb}^{2}",2,2);

  COMB_L.AddFrame(L_X2a);
  COMB_L.SetNElementsForFrame(L_X2a, 0);
  COMB_L.AddFrame(L_X2b);
  COMB_L.SetNElementsForFrame(L_X2b, 0);
          
  COMB_L.AddJigsaw(CombSplitSq_L);
  CombSplitSq_L.AddCombFrame(L_X2a, 0);
  CombSplitSq_L.AddCombFrame(L_X2b, 1);
  CombSplitSq_L.AddObjectFrames(X2a.GetListVisibleFrames(), 0);
  CombSplitSq_L.AddObjectFrames(X2b.GetListVisibleFrames(), 1);
        
  COMB_J.AddFrame(ISR);
  COMB_J.SetNElementsForFrame(ISR, 1);
  COMB_J.AddFrame(J_X2a);
  COMB_J.SetNElementsForFrame(J_X2a, 0);
  COMB_J.AddFrame(J_X2b);
  COMB_J.SetNElementsForFrame(J_X2b, 0);
      
  COMB_J.AddJigsaw(CombSplitSq_J_Transverse);
  CombSplitSq_J_Transverse.SetTransverse();
  CombSplitSq_J_Transverse.AddCombFrame(ISR, 0);
  CombSplitSq_J_Transverse.AddCombFrame(J_X2a, 1);
  CombSplitSq_J_Transverse.AddCombFrame(J_X2b, 1);
  CombSplitSq_J_Transverse.AddObjectFrame(ISR, 0);
  CombSplitSq_J_Transverse.AddObjectFrame(S, 1);

  COMB_J.AddJigsaw(CombSplitSq_JS);
  CombSplitSq_JS.AddCombFrame(J_X2a, 0);
  CombSplitSq_JS.AddCombFrame(J_X2b, 1);
  CombSplitSq_JS.AddObjectFrame(X2a, 0);
  CombSplitSq_JS.AddObjectFrame(X2b, 1);

  InvisibleGroup INV("INV","Invisible System");
  INV.AddFrame(Ia);
  INV.AddFrame(Ib);
  
  SetMassInvJigsaw InvM("InvM", "Set inv. system mass");
  INV.AddJigsaw(InvM);
  
  SetRapidityInvJigsaw InvEta("InvEta", "Set inv. system rapidity");
  INV.AddJigsaw(InvEta);
  InvEta.AddVisibleFrames(S.GetListVisibleFrames());
  
  MinMassesSqInvJigsaw InvSplit("InvSplit", "INV -> I_{a} + I_{b}", 2);
  INV.AddJigsaw(InvSplit);
  InvSplit.AddVisibleFrame(J_X2a, 0);
  InvSplit.AddVisibleFrame(J_X2b, 1);
  InvSplit.AddVisibleFrame(L_X2a, 0);
  InvSplit.AddVisibleFrame(L_X2b, 1);
  InvSplit.AddInvisibleFrame(Ia, 0);
  InvSplit.AddInvisibleFrame(Ib, 1);

  LAB.InitializeAnalysis();
  if(treeplot){
    tree_plot.SetTree(LAB);
    tree_plot.Draw("ISR_Sparticle2_tree", "Analysis Tree", treeplot_INV);
    tree_plot.SetTree(INV);
    tree_plot.Draw("ISR_Sparticle2_inv", "Invisible Jigsaws", treeplot_INV);
    tree_plot.SetTree(COMB_J);
    tree_plot.Draw("ISR_Sparticle2_jets", "Jet Jigsaws", treeplot_INV);
    tree_plot.SetTree(COMB_L);
    tree_plot.Draw("ISR_Sparticle2_leps", "Lep Jigsaws", treeplot_INV);

    tree_plot.WriteOutput("trees.root");
    std::cout << "Writing trees to trees.root" << endl;
    return 0;
  }

  int quark_cut = 0;
  int tot_quark = 0;

  string g_Label = "";
  double g_NX = 128;
  vector<TH1*> hists1;
  vector<TH2*> hists2;
  vector<TEfficiency*> effs;
  string title = proc;

  double DM = 0.;
  if(proc.find("TChi") != std::string::npos)
    DM = std::stod(proc.substr(proc.find("_") + 1));

  TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET;").c_str(), g_NX, 0., 500.);
  hists1.push_back(hist_MET);
  TH1D* hist_NQuarks = new TH1D((title+"_NQuarks").c_str(), (title+"_NQuarks;NQuarks").c_str(), 6, 0, 6);
  hists1.push_back(hist_NQuarks);

  TEfficiency* eff_quarkMatchedGenJet_eta = new TEfficiency((title+"_eff_quarkMatchedGenJet_eta").c_str(),"Efficiency of Quark Getting Matched to GenJet;quarkEta;Quark Matching Efficiency", g_NX, -5., 5.);
  effs.push_back(eff_quarkMatchedGenJet_eta);

  TH1D* hist_quarkPt = new TH1D((title+"_quarkPt").c_str(), (title+"_quarkPt;quarkPt").c_str(), g_NX, 0., 200.);
  //hists1.push_back(hist_quarkPt);
  TH1D* hist_quarkEta = new TH1D((title+"_quarkEta").c_str(), (title+"_quarkEta;quarkEta").c_str(), g_NX, -5., 5.);
  //hists1.push_back(hist_quarkEta);

  TH1D* hist_DeltaR_quark_genjet = new TH1D((title+"_DeltaR_quark_genjet").c_str(), (title+"_DeltaR_quark_genjet;DeltaR(quark,genjet)").c_str(), g_NX, 0., 1.);
  //hists1.push_back(hist_DeltaR_quark_genjet);
  TH1D* hist_DeltaR_quark_jet = new TH1D((title+"_DeltaR_quark_jet").c_str(), (title+"_DeltaR_quark_jet;DeltaR(quark,jet)").c_str(), g_NX, 0., 1.);
  //hists1.push_back(hist_DeltaR_quark_jet);

  TH1D* hist_genW_pt = new TH1D((title+"_genW_pt").c_str(), (title+"_genW_pt;genW Pt").c_str(), g_NX, 0., 250.);
  //hists1.push_back(hist_genW_pt);
  TH1D* hist_genW_eta = new TH1D((title+"_genW_eta").c_str(), (title+"_genW_eta;genW Eta").c_str(), g_NX, -5., 5.);
  //hists1.push_back(hist_genW_eta);
  TH1D* hist_genW_mass = new TH1D((title+"_genW_mass").c_str(), (title+"_genW_mass;genW Mass").c_str(), g_NX, 50., 110.);
  //hists1.push_back(hist_genW_mass);

  TH1D* hist_Count_ISR_Sparticle = new TH1D((title+"_Count_ISR_Sparticle").c_str(), (title+"_Count_ISR_Sparticle;").c_str(), 20, 0, 20); // generic counting
  
  TH2D* hist_matched_Njets_S_unmatched_Njets_S = new TH2D((title+"_matched_Njets_S_unmatched_Njets_S").c_str(), (title+"_matched_Njets_S_unmatched_Njets_S;matched Njets S;unmatched Njets S").c_str(), 5, 0, 5, 5, 0, 5);
  //hists2.push_back(hist_matched_Njets_S_unmatched_Njets_S);
  TH2D* hist_matched_Njets_ISR_unmatched_Njets_ISR = new TH2D((title+"_matched_Njets_ISR_unmatched_Njets_ISR").c_str(), (title+"_matched_Njets_ISR_unmatched_Njets_ISR;matched Njets ISR;unmatched Njets ISR").c_str(), 5, 0, 5, 5, 0, 5);
  //hists2.push_back(hist_matched_Njets_ISR_unmatched_Njets_ISR);
  TH2D* hist_matched_Njets_Sa_matched_Njets_Sb = new TH2D((title+"_matched_Njets_Sa_matched_Njets_Sb").c_str(), (title+"_matched_Njets_Sa_matched_Njets_Sb;matched Njets Sa;matched Njets Sb").c_str(), 5, 0, 5, 5, 0, 5);
  //hists2.push_back(hist_matched_Njets_Sa_matched_Njets_Sb);
  TH2D* hist_unmatched_Njets_Sa_unmatched_Njets_Sb = new TH2D((title+"_unmatched_Njets_Sa_unmatched_Njets_Sb").c_str(), (title+"_unmatched_Njets_Sa_unmatched_Njets_Sb;unmatched Njets Sa;unmatched Njets Sb").c_str(), 5, 0, 5, 5, 0, 5);
  //hists2.push_back(hist_unmatched_Njets_Sa_unmatched_Njets_Sb);

  // Sparticle2 KIN Plots
  TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  hists2.push_back(hist_RISR_PTISR);
  TH2D* hist_dphiCMI_PTCM = new TH2D((title+"_dphiCMI_PTCM").c_str(), (title+"_dphiCMI_PTCM;#Delta #phi_{(CM,I)};p_{T}^{CM}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 0., 500.);
  hists2.push_back(hist_dphiCMI_PTCM);
  TH2D* hist_dphiMETV_PTISR = new TH2D((title+"_dphiMETV_PTISR").c_str(), (title+"_dphiMETV_PTISR;#Delta #phi_{(I,V)};p_{T}^{ISR}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 200., 800.);
  hists2.push_back(hist_dphiMETV_PTISR);
  TH2D* hist_gammaPerp_RISR = new TH2D((title+"_gammaPerp_RISR").c_str(), (title+"_gammaPerp_RISR;#gamma_{#perp};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_gammaPerp_RISR);
  TH2D* hist_MQperp_RISR = new TH2D((title+"_MQperp_RISR").c_str(), (title+"_MQperp_RISR;M_{#perp};RISR").c_str(), g_NX/2., 0., 100., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_MQperp_RISR);
  TH2D* hist_RISR_mllLEAD = new TH2D((title+"_RISR_mllLEAD").c_str(), (title+"_RISR_mllLEAD;R_{ISR};m_{ll}LEAD").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
  hists2.push_back(hist_RISR_mllLEAD);
  TH2D* hist_RISR_mL = new TH2D((title+"_RISR_mL").c_str(), (title+"_RISR_mL;R_{ISR};mL").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.); // mass of leptonic system
  hists2.push_back(hist_RISR_mL);

  TH2D* hist_MSperpCM0_RISR = new TH2D((title+"_MSperpCM0_RISR").c_str(), (title+"_MSperpCM0_RISR;MS_{#perp CM0};RISR").c_str(), g_NX/2., 0., 500., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_MSperpCM0_RISR);
  TH2D* hist_MQperpCM0_RISR = new TH2D((title+"_MQperpCM0_RISR").c_str(), (title+"_MQperpCM0_RISR;MQ_{#perp CM0};RISR").c_str(), g_NX/2., 0., 300., g_NX/2., 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
  hists2.push_back(hist_MQperpCM0_RISR);
  TH2D* hist_gammaPerpCM0_RISR = new TH2D((title+"_gammaPerpCM0_RISR").c_str(), (title+"_gammaPerpCM0_RISR;#gamma_{#perp CM0};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_gammaPerpCM0_RISR);
  TH2D* hist_MVisAperpCM_MVisBperpCM = new TH2D((title+"_MVisAperpCM_MVisBperpCM").c_str(), (title+"_MVisAperpCM_MVisBperpCM;M^{VisA}_{#perp CM0};M^{VisB}_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MVisAperpCM_MVisBperpCM);

  TH2D* hist_MSperpCM0_gammaPerpCM0 = new TH2D((title+"_MSperpCM0_gammaPerpCM0").c_str(), (title+"_MSperpCM0_gammaPerpCM0;MS_{#perp CM0};#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 1.);
  hists2.push_back(hist_MSperpCM0_gammaPerpCM0);
  TH2D* hist_MSperpCM0_MQperpCM0 = new TH2D((title+"_MSperpCM0_MQperpCM0").c_str(), (title+"_MSperpCM0_MQperpCM0;MS_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 300.);
  hists2.push_back(hist_MSperpCM0_MQperpCM0);
  TH2D* hist_gammaPerpCM0_MQperpCM0 = new TH2D((title+"_gammaPerpCM0_MQperpCM0").c_str(), (title+"_gammaPerpCM0_MQperpCM0;#gamma_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 300.);
  hists2.push_back(hist_gammaPerpCM0_MQperpCM0);
  
  TH2D* hist_MQperpCM0_MVisAperpCM = new TH2D((title+"_MQperpCM0_MVisAperpCM").c_str(), (title+"_MQperpCM0_MVisAperpCM;MQ_{#perp CM0};M^{VisA}_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MQperpCM0_MVisAperpCM);
  TH2D* hist_MQperpCM0_MVisBperpCM = new TH2D((title+"_MQperpCM0_MVisBperpCM").c_str(), (title+"_MQperpCM0_MVisBperpCM;MQ_{#perp CM0};M^{VisB}_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MQperpCM0_MVisBperpCM);

  TH2D* hist_MJ_MQperpCM0 = new TH2D((title+"_MJ_MQperpCM0").c_str(), (title+"_MJ_MQperpCM0;MJ;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MJ_MQperpCM0);
  TH2D* hist_MLa_MQperpCM0 = new TH2D((title+"_MLa_MQperpCM0").c_str(), (title+"_MLa_MQperpCM0;MLa;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLa_MQperpCM0);
  TH2D* hist_MLb_MQperpCM0 = new TH2D((title+"_MLb_MQperpCM0").c_str(), (title+"_MLb_MQperpCM0;MLb;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLb_MQperpCM0);
  TH2D* hist_MJ_MSperpCM0 = new TH2D((title+"_MJ_MSperpCM0").c_str(), (title+"_MJ_MSperpCM0;MJ;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MJ_MSperpCM0);
  TH2D* hist_MLa_MSperpCM0 = new TH2D((title+"_MLa_MSperpCM0").c_str(), (title+"_MLa_MSperpCM0;MLa;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MLa_MSperpCM0);
  TH2D* hist_MLb_MSperpCM0 = new TH2D((title+"_MLb_MSperpCM0").c_str(), (title+"_MLb_MSperpCM0;MLb;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MLb_MSperpCM0);
  TH2D* hist_MJ_gammaPerpCM0 = new TH2D((title+"_MJ_gammaPerpCM0").c_str(), (title+"_MJ_gammaPerpCM0;MJ;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MJ_gammaPerpCM0);
  TH2D* hist_MLa_gammaPerpCM0 = new TH2D((title+"_MLa_gammaPerpCM0").c_str(), (title+"_MLa_gammaPerpCM0;MLa;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MLa_gammaPerpCM0);
  TH2D* hist_MLb_gammaPerpCM0 = new TH2D((title+"_MLb_gammaPerpCM0").c_str(), (title+"_MLb_gammaPerpCM0;MLb;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MLb_gammaPerpCM0);

  TH2D* hist_MSCM0_RISR = new TH2D((title+"_MSCM0_RISR").c_str(), (title+"_MSCM0_RISR;MS_{CM0};RISR").c_str(), g_NX/2., 0., 500., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_MSCM0_RISR);
  TH2D* hist_MQCM0_RISR = new TH2D((title+"_MQCM0_RISR").c_str(), (title+"_MQCM0_RISR;MQ_{CM0};RISR").c_str(), g_NX/2., 0., 300., g_NX/2., 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
  hists2.push_back(hist_MQCM0_RISR);
  TH2D* hist_gammaCM0_RISR = new TH2D((title+"_gammaCM0_RISR").c_str(), (title+"_gammaCM0_RISR;#gamma_{CM0};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
  hists2.push_back(hist_gammaCM0_RISR);
  TH2D* hist_MVisACM_MVisBCM = new TH2D((title+"_MVisACM_MVisBCM").c_str(), (title+"_MVisACM_MVisBCM;M^{VisA}_{CM};M^{VisB}_{CM}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MVisACM_MVisBCM);

  TH2D* hist_MSCM0_gammaCM0 = new TH2D((title+"_MSCM0_gammaCM0").c_str(), (title+"_MSCM0_gammaCM0;MS_{CM0};#gamma_{CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 1.);
  hists2.push_back(hist_MSCM0_gammaCM0);
  TH2D* hist_MSCM0_MQCM0 = new TH2D((title+"_MSCM0_MQCM0").c_str(), (title+"_MSCM0_MQCM0;MS_{CM0};MQ_{CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 300.);
  hists2.push_back(hist_MSCM0_MQCM0);
  TH2D* hist_gammaCM0_MQCM0 = new TH2D((title+"_gammaCM0_MQCM0").c_str(), (title+"_gammaCM0_MQCM0;#gamma_{CM0};MQ_{CM0}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 300.);
  hists2.push_back(hist_gammaCM0_MQCM0);
  
  TH2D* hist_MQCM0_MVisACM = new TH2D((title+"_MQCM0_MVisACM").c_str(), (title+"_MQCM0_MVisACM;MQ_{CM0};M^{VisA}_{CM}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MQCM0_MVisACM);
  TH2D* hist_MQCM0_MVisBCM = new TH2D((title+"_MQCM0_MVisBCM").c_str(), (title+"_MQCM0_MVisBCM;MQ_{CM0};M^{VisB}_{CM}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MQCM0_MVisBCM);

  TH2D* hist_MJ_MQCM0 = new TH2D((title+"_MJ_MQCM0").c_str(), (title+"_MJ_MQCM0;MJ;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MJ_MQCM0);
  TH2D* hist_MLa_MQCM0 = new TH2D((title+"_MLa_MQCM0").c_str(), (title+"_MLa_MQCM0;MLa;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLa_MQCM0);
  TH2D* hist_MLb_MQCM0 = new TH2D((title+"_MLb_MQCM0").c_str(), (title+"_MLb_MQCM0;MLb;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLb_MQCM0);
  TH2D* hist_MJ_MSCM0 = new TH2D((title+"_MJ_MSCM0").c_str(), (title+"_MJ_MSCM0;MJ;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MJ_MSCM0);
  TH2D* hist_MLa_MSCM0 = new TH2D((title+"_MLa_MSCM0").c_str(), (title+"_MLa_MSCM0;MLa;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MLa_MSCM0);
  TH2D* hist_MLb_MSCM0 = new TH2D((title+"_MLb_MSCM0").c_str(), (title+"_MLb_MSCM0;MLb;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
  hists2.push_back(hist_MLb_MSCM0);
  TH2D* hist_MJ_gammaCM0 = new TH2D((title+"_MJ_gammaCM0").c_str(), (title+"_MJ_gammaCM0;MJ;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MJ_gammaCM0);
  TH2D* hist_MLa_gammaCM0 = new TH2D((title+"_MLa_gammaCM0").c_str(), (title+"_MLa_gammaCM0;MLa;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MLa_gammaCM0);
  TH2D* hist_MLb_gammaCM0 = new TH2D((title+"_MLb_gammaCM0").c_str(), (title+"_MLb_gammaCM0;MLb;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
  hists2.push_back(hist_MLb_gammaCM0);

  TH2D* hist_MLb_MJ = new TH2D((title+"_MLb_MJ").c_str(), (title+"_MLb_MJ;MLb;MJ").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLb_MJ);
  TH2D* hist_MLa_MJ = new TH2D((title+"_MLa_MJ").c_str(), (title+"_MLa_MJ;MLa;MJ").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLa_MJ);
  TH2D* hist_MLa_MLb = new TH2D((title+"_MLa_MLb").c_str(), (title+"_MLa_MLb;MLa;MLb").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
  hists2.push_back(hist_MLa_MLb);

  TH2D* hist_RISR_MJ = new TH2D((title+"_RISR_MJ").c_str(), (title+"_RISR_MJ;R_{ISR};MJ").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
  hists2.push_back(hist_RISR_MJ);
  TH2D* hist_RISR_MLa = new TH2D((title+"_RISR_MLa").c_str(), (title+"_RISR_MLa;R_{ISR};MLa").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
  hists2.push_back(hist_RISR_MLa);
  TH2D* hist_RISR_MLb = new TH2D((title+"_RISR_MLb").c_str(), (title+"_RISR_MLb;R_{ISR};MLb").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
  hists2.push_back(hist_RISR_MLb);

  TH2D* hist_EventCount = new TH2D("EventCount", "EventCount", 12, 0., 12., 20, 0., 20.);
   
  if(debug) cout << "Loading TChains" << endl;
  TChain* chain = new TChain("Events");
  chain->Add(ifile.c_str());

  if(debug) cout << "Making base class" << endl;
  SUSYNANOBase* base = new SUSYNANOBase(chain, era == "Run3"); // if era is equal to Run3 then pass true and use Run3 branches
  int processed_events = 0;
  int kept_file_events = 0;
  int SKIP = 1; // 1 is default (no skipping of events)
  
  if(!debug) gErrorIgnoreLevel = kFatal;
  int NTOT = 0;
  if(debug) cout << "Getting nentries" << endl;
  NTOT = base->fChain->GetEntries();
  if(debug) cout << "Got " << NTOT << " entries" << endl;
  gErrorIgnoreLevel = 0;
  if(nchunk < 1 || ichunk < 1 || ichunk > nchunk){
    ichunk = 1;
    nchunk = 1;
  }
  Long64_t N1, N0;
  if(nchunk >= NTOT){
    N1 = ichunk;
    N0 = ichunk-1;
  } else {
    N1 = NTOT/nchunk;
    if(NTOT%nchunk > 0)
      N1++;
    N0 = (ichunk-1)*N1;
    N1 = N0 + N1;
  }
  
  if(debug) cout << "RUNNING EVENT LOOP" << endl;
  // event loop
  for(Long64_t e = N0; e < N1 && e < NTOT; e+=SKIP){
    processed_events++;
    
    base->GetEntry(e);
    int mymod = (N1-N0)/10;
    if(mymod < 1)
      mymod = 1;
    if(e%mymod == 0)
      cout << " event = " << e << " : [" << N0 << " , " << N1 << "]" << endl;
    double hweight = base->genWeight*weight*double(SKIP)/1.e7; // need 1.e7 factor to normalize from input!
    int N = base->nGenPart;
    int PDGID;
    int MP = 0;
    int MC = 0;
    if(proc.find("TChi") != std::string::npos || proc.find("Cascade") != std::string::npos){
      for(int i = 0; i < N; i++){
        PDGID = abs(base->GenPart_pdgId[i]);
        if(PDGID > 1000000 && PDGID < 3000000){
          int mass = int(base->GenPart_mass[i]+0.5);
          if(PDGID == 1000022)
            MC = mass;
          else
            if(mass > MP)
              MP = mass;
        }
      }
    }
    int e_DM = MP - MC;

    if(DM > 0.)
      if(MP != 0 && MC != 0 && (e_DM > (DM + 1.) || e_DM < (DM - 1.))) continue;

    int slep = 0;
    int snu = 0;
    if(proc.find("Cascade") != std::string::npos){
      for(int i = 0; i < N; i++){
        PDGID = abs(base->GenPart_pdgId[i]);
        int status = base->GenPart_status[i];
        if(status == 22 && (PDGID == 1000011 || PDGID == 1000013))
          slep++;
        if(status == 22 && (PDGID == 1000012 || PDGID == 1000014))
          snu++;
      }
    }

    int year = 2017;
    // loop over jets 
    ParticleList jets;
    int Njet = base->nJet;
    for(int i = 0; i < Njet; i++){
      Particle jet;
      float mass = base->Jet_mass[i];
      if(std::isnan(mass))
        mass = 0;
      if(std::isinf(mass))
        mass = 0;
      if(mass < 0.)
        mass = 0.;
      jet.SetPtEtaPhiM(base->Jet_pt[i], base->Jet_eta[i],
          	     base->Jet_phi[i], mass);
      int jetId;
      if(era == "Run3")
        jetId = base->Run3Jet_jetId[i];
      else
        jetId = base->Jet_jetId[i];
      if(jetId < 3) continue;
      // DeepFlavour tagger
      jet.SetBtag(base->Jet_btagDeepFlavB[i]);
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
      // Deep Flavor
      if(jet.Btag() > 0.7489)
        jet.SetBtagID(kTight);
      else if(jet.Btag() > 0.3033) 
        jet.SetBtagID(kMedium);
      else if(jet.Btag() > 0.0521)
        jet.SetBtagID(kLoose);
      if(era == "Run3")
        jet.SetPDGID( base->Run3Jet_partonFlavour[i] );
      else
        jet.SetPDGID( base->Jet_partonFlavour[i] );
      jet.SetGenIndex(-1);
      jet.SetMomPDGID(0);
      jet.SetGenMomIndex(-1);
      jets.push_back(jet);
    }

     ParticleList electrons;

     int Nele = base->nElectron;
     for(int i = 0; i < Nele; i++){
       // baseline lepton definition
       //if(base->Electron_pt[i] < 5. || fabs(base->Electron_eta[i]) > 2.5)
       //  continue;
       if(fabs(base->Electron_dxy[i]) >= 0.05 || fabs(base->Electron_dz[i]) >= 0.1 ||
          base->Electron_sip3d[i] >= 8)
         continue;
       if(base->Electron_pfRelIso03_all[i]*base->Electron_pt[i] >= 20. + 300./base->Electron_pt[i])
         continue;

       Particle lep;
       lep.SetPtEtaPhiM(base->Electron_pt[i], base->Electron_eta[i],
           	     base->Electron_phi[i], std::max(base->Electron_mass[i],float(1.e-6)));
       lep.SetPDGID( (base->Electron_charge[i] < 0. ? 11 : -11) );
       lep.SetCharge( (base->Electron_charge[i] < 0. ? -1 : 1) );

       lep.SetDxy(base->Electron_dxy[i]);
       lep.SetDxyErr(base->Electron_dxyErr[i]);
       lep.SetDz(base->Electron_dz[i]);
       lep.SetDzErr(base->Electron_dzErr[i]);
       lep.SetIP3D(base->Electron_ip3d[i]);
       lep.SetSIP3D(base->Electron_sip3d[i]);

       lep.SetRelIso(base->Electron_pfRelIso03_all[i]);
       lep.SetMiniIso(base->Electron_miniPFRelIso_all[i]);

       // https://twiki.cern.ch/twiki/pub/CMS/SUSLeptonSF/Run2_SUSYwp_EleCB_MVA_8Jan19.pdf
       
       // FO baseline criteria
       if(base->Electron_lostHits[i] == 0 && base->Electron_convVeto[i]){
         if(era == "Run3"){
           if(base->Electron_mvaIso_WP80[i]){ // run3
             if(lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
                   lep.SetLepQual(kBronze);
                 else if(lep.SIP3D() > 2.)
                   lep.SetLepQual(kSilver);
                 else
                   lep.SetLepQual(kGold);
             electrons.push_back(lep);
           }
         } // Run3
         else{
           if(base->Electron_mvaFall17V2Iso_WP80[i]){ // run2 2017
             if(lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
                   lep.SetLepQual(kBronze);
                 else if(lep.SIP3D() > 2.)
                   lep.SetLepQual(kSilver);
                 else
                   lep.SetLepQual(kGold);
             electrons.push_back(lep);
           }
         } // Run2
       } // if(base->Electron_lostHits[i] == 0 && base->Electron_convVeto[i])
     } // for(int i = 0; i < Nele; i++)

     ParticleList muons;
     int Nmu = base->nMuon;
     for(int i = 0; i < Nmu; i++){
       // baseline lepton definition
       if(base->Muon_pt[i] < 3. || fabs(base->Muon_eta[i]) > 2.5)
         continue;
       if(fabs(base->Muon_dxy[i]) >= 0.05 || fabs(base->Muon_dz[i]) >= 0.1 || base->Muon_sip3d[i] >= 8.)
         continue;
       if(base->Muon_pfRelIso03_all[i]*base->Muon_pt[i] >= 20. + 300./base->Muon_pt[i])
         continue;
       
       Particle lep;
       lep.SetPtEtaPhiM(base->Muon_pt[i], base->Muon_eta[i],
           	     base->Muon_phi[i], std::max(float(0.),base->Muon_mass[i]));
       lep.SetPDGID( (base->Muon_charge[i] < 0. ? 13 : -13) );
       lep.SetCharge( (base->Muon_charge[i] < 0. ? -1 : 1) );	
       lep.SetDxy(base->Muon_dxy[i]);
       lep.SetDxyErr(base->Muon_dxyErr[i]);
       lep.SetDz(base->Muon_dz[i]);
       lep.SetDzErr(base->Muon_dzErr[i]);
       lep.SetIP3D(base->Muon_ip3d[i]);
       lep.SetSIP3D(base->Muon_sip3d[i]);
       lep.SetRelIso(base->Muon_pfRelIso03_all[i]);
       lep.SetMiniIso(base->Muon_miniPFRelIso_all[i]);
       // FO baseline criteria
       if(true){
         lep.SetParticleID(kLoose);

         // signal lep criteria
         //if(lep.IP3D() < 0.01 && lep.SIP3D() < 2.){
         if(true){
           if(base->Muon_tightId[i])
             lep.SetParticleID(kTight);
           else if(lep.Pt() < 0.){
             if(base->Muon_softId[i])
               lep.SetParticleID(kMedium);
           } else {
             if(base->Muon_mediumId[i])
               lep.SetParticleID(kMedium);
           }
         }
       }
       if(lep.ParticleID() < kMedium || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
               lep.SetLepQual(kBronze);
             else if(lep.SIP3D() > 2.)
               lep.SetLepQual(kSilver);
             else
               lep.SetLepQual(kGold);
       muons.push_back(lep);
     }



     ParticleList gen_jets;
     int NGenjet = base->nGenJet;

     for(int i = 0; i < NGenjet; i++){
       //if(base->GenJet_pt[i] < 15. || fabs(base->GenJet_eta[i]) > 5.)
       //  continue;
       Particle jet;
       float mass = base->GenJet_mass[i];
       if(std::isnan(mass))
         mass = 0;
       if(std::isinf(mass))
         mass = 0;
       if(mass < 0.)
         mass = 0.;
       jet.SetPtEtaPhiM(base->GenJet_pt[i], base->GenJet_eta[i],
           	     base->GenJet_phi[i], mass);
        jet.SetGenIndex(-1);
        jet.SetMomPDGID(0);
        jet.SetGenMomIndex(-1);

       gen_jets.push_back(jet);
     }

     ParticleList gen_electrons;
     int gen_Nleps_A = 0;
     int gen_Nleps_B = 0;
    
     for(int i = 0; i < N; i++){
       PDGID = base->GenPart_pdgId[i];
       if(abs(PDGID) == 11 && base->GenPart_pt[i] > 2. && base->GenPart_status[i] == 1){
         Particle lep;
         
         lep.SetPDGID(PDGID);
         int mom = 0;
         if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[i];
         else mom = base->GenPart_genPartIdxMother[i];
         if(mom >= 0 && mom < N){
    
          int momID = base->GenPart_pdgId[mom];
          int momStatus = base->GenPart_status[mom];
    
          while(abs(momID) == 11){
    
            if(momStatus == 23){
              lep.SetMomPDGID(PDGID);
              lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
              break;
            }
    
            if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[mom];
            else mom = base->GenPart_genPartIdxMother[mom];
            if(mom < 0 || mom >= N)
              continue;
            momID = base->GenPart_pdgId[mom];
            momStatus = base->GenPart_status[mom];
    
          }
          
          lep.SetMomPDGID(momID);
          lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
         }
        
         lep.SetCharge( (PDGID > 0 ? -1 : 1) );
         lep.SetPtEtaPhiM(base->GenPart_pt[i], base->GenPart_eta[i],
          	       base->GenPart_phi[i], max(float(0.),base->GenPart_mass[i]));
         
         gen_electrons.push_back(lep);
       }
     }
    
     ParticleList gen_muons;
     for(int i = 0; i < N; i++){
       PDGID = base->GenPart_pdgId[i];
       if(abs(PDGID) == 13 && base->GenPart_pt[i] > 2. && base->GenPart_status[i] == 1){
         Particle lep;
         lep.SetPDGID(PDGID);
         int mom = 0;
         if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[i];
         else mom = base->GenPart_genPartIdxMother[i];
         if(mom >= 0 && mom < N){
          int momID = base->GenPart_pdgId[mom];
          int momStatus = base->GenPart_status[mom];
          while(abs(momID) == 13){
            if(momStatus == 23){
              lep.SetMomPDGID(PDGID);
              lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
              break;
            }
            if(mom < 0 || mom >= N)
              continue;
            if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[mom];
            else mom = base->GenPart_genPartIdxMother[mom];
            momID = base->GenPart_pdgId[mom];
            momStatus = base->GenPart_status[mom];
          }
          lep.SetMomPDGID(momID);
          lep.SetSourceID(GetLepSource(PDGID, PDGID, momID));
         }
         lep.SetCharge( (PDGID > 0 ? -1 : 1) );
         lep.SetPtEtaPhiM(base->GenPart_pt[i], base->GenPart_eta[i],
          	       base->GenPart_phi[i], max(float(0.),base->GenPart_mass[i]));
         gen_muons.push_back(lep);
       }
     }
    
     ParticleList gen_bosons;
     for(int i = 0; i < N; i++){
       PDGID = base->GenPart_pdgId[i];
       if(abs(PDGID) == 23 || abs(PDGID) == 24 || abs(PDGID) == 25){
         if(base->GenPart_status[i] != 22) continue;
         Particle p;
         p.SetPDGID(PDGID);
         p.SetGenIndex(i);
         int mom = 0;
         if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[i];
         else mom = base->GenPart_genPartIdxMother[i];
         if(mom >= 0 && mom < N)
          p.SetMomPDGID(base->GenPart_pdgId[mom]);
          p.SetGenMomIndex(mom);
          p.SetPtEtaPhiM(base->GenPart_pt[i], base->GenPart_eta[i],
          	     base->GenPart_phi[i], max(float(0.),base->GenPart_mass[i]));
         gen_bosons.push_back(p);
       }
    }

    ParticleList gen_quarks;
    ParticleList gen_bquarks;
    ParticleList gen_hadbosons;
    for(int i = 0; i < N; i++){
      PDGID = base->GenPart_pdgId[i];
      int status = base->GenPart_status[i];
      int mom = 0;
      if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[i];
      else mom = base->GenPart_genPartIdxMother[i];
      int momID = base->GenPart_pdgId[mom];
      if(abs(PDGID) > 5 || status != 23) continue;
      bool found = false;
      Particle p;
      while(mom >= 0 && mom < N){
        if(found) break;
        momID = base->GenPart_pdgId[mom];
        if(abs(PDGID) == 5 && abs(momID) == 6){ // b-quark from top
          p.SetGenIndex(i);
          p.SetMomPDGID(momID);
          p.SetGenMomIndex(mom);
          p.SetPtEtaPhiM(base->GenPart_pt[i], base->GenPart_eta[i],
               	       base->GenPart_phi[i], max(float(0.),base->GenPart_mass[i]));
          gen_bquarks.push_back(p);
          found = true;
        }
      // check quarks are in acceptance
      if(boson_acceptance_cut && (base->GenPart_pt[i] < 5. || abs(base->GenPart_eta[i]) > 2.7)) break;
      if((abs(momID) == 24 || abs(momID) == 23) && base->GenPart_status[mom] == 22){
        int grandmom = 0;
        if(era == "Run3") grandmom = base->Run3GenPart_genPartIdxMother[i];
        else base->GenPart_genPartIdxMother[mom];
        while(grandmom >= 0 && grandmom < N){
          int grandmomID = base->GenPart_pdgId[grandmom];
          int grandmomStatus = base->GenPart_status[grandmom];
          if((abs(grandmomID) == 6 || abs(grandmomID) > 1000000) || (abs(grandmomID) < 6 && grandmomStatus == 21)){
            p.SetGenIndex(i);
            p.SetMomPDGID(momID);
            p.SetGenMomIndex(mom);
            p.SetPtEtaPhiM(base->GenPart_pt[i], base->GenPart_eta[i],
             	     base->GenPart_phi[i], max(float(0.),base->GenPart_mass[i]));
            gen_quarks.push_back(p);
            found = true;
            // gen W info
            // check if we already have found this mom from another quark
            bool already_found_mom = false;
            for(int mi = 0; mi < int(gen_hadbosons.size()); mi++)
              if(gen_hadbosons[mi].GenIndex() == mom) already_found_mom = true;
            if(!already_found_mom){
              hist_genW_pt->Fill(base->GenPart_pt[mom], hweight);
              hist_genW_eta->Fill(base->GenPart_eta[mom], hweight);
              hist_genW_mass->Fill(base->GenPart_mass[mom], hweight);
              Particle m;
              m.SetPDGID(base->GenPart_pdgId[mom]);
              m.SetGenIndex(mom);
              m.SetGenMomIndex(grandmom);
              m.SetMomPDGID(base->GenPart_pdgId[grandmom]);
              m.SetPtEtaPhiM(base->GenPart_pt[mom], base->GenPart_eta[mom],
                             base->GenPart_phi[mom], max(float(0.),base->GenPart_mass[mom]));
              gen_hadbosons.push_back(m);
            }
            break;
          }
          if(era == "Run3") grandmom = base->Run3GenPart_genPartIdxMother[grandmom];
          else grandmom = base->GenPart_genPartIdxMother[grandmom];
        } // while(grandmom >= 0 && grandmom < N)
      }
      if(era == "Run3") mom = base->Run3GenPart_genPartIdxMother[mom];
      else mom = base->GenPart_genPartIdxMother[mom];
      }
    }
          
    ParticleList allquarks = gen_quarks+gen_bquarks;
    int NQuarks = gen_quarks.size();
    int NBQuarks = gen_bquarks.size();
    tot_quark += NQuarks;
    quark_cut += (NQuarks - gen_quarks.PtEtaCut(20.,2.4).size());
    int Nhadbosons = gen_hadbosons.size();
    int bosons_in_acceptance = 0;
    for(int i = 0; i < NQuarks; i++){
      hist_quarkPt->Fill(gen_quarks[i].Pt(), hweight);
      hist_quarkEta->Fill(gen_quarks[i].Eta(), hweight);
      int GenMomIndex = gen_quarks[i].GenMomIndex();
      for(int j = i + 1; j < NQuarks; j++){
        if(GenMomIndex == gen_quarks[j].GenMomIndex()){
          bosons_in_acceptance++;
          break;
        }
      }
    }
    if(boson_acceptance_cut)
      Nhadbosons = bosons_in_acceptance;
    if(boson_acceptance_cut && Nhadbosons < 1) continue;

    muons = muons.ParticleIDCut(kVeryLoose);
    muons = muons.PtEtaCut(3.0,2.5);
    electrons = electrons.ParticleIDCut(kVeryLoose);
    electrons = electrons.PtEtaCut(5.0,2.5); 
    ParticleList leptons = electrons+muons;
    leptons.SortByPt();
    jets = jets.PtEtaCut(20., 2.4);
    jets = jets.RemoveOverlap(leptons, 0.2);
    jets.SortByPt();
    gen_jets = gen_jets.PtEtaCut(20.,2.4);
    gen_jets.SortByPt();
    gen_electrons = gen_electrons.PtEtaCut(1.,2.5);
    gen_muons = gen_muons.PtEtaCut(1.,2.5);
    ParticleList gen_leptons = gen_electrons+gen_muons;
    gen_leptons.SortByPt();

    int Njets = int(jets.size());
    int Nleps = int(leptons.size());
    int Nobjs = Njets + Nleps;
    if(Nobjs < 2 || Njets < 1) continue;
    if(min_lep_cut > -1 && Nleps < min_lep_cut) continue;
    if(max_lep_cut > -1 && Nleps > max_lep_cut) continue;
    int Ngen_jets = int(gen_jets.size());
    int Ngen_leps = int(gen_leptons.size());

    int Ngen_sig_leps = 0; // leps coming from something interesting
    for(int i = 0; i < Ngen_leps; i++){
      if(gen_leptons[i].SourceID() == kSignal)
        Ngen_sig_leps++;
    }
    if(gen_lepton_cut && Ngen_sig_leps < 1) continue;

    for(int i = 0; i < NQuarks; i++)
      for(int j = 0; j < Ngen_jets; j++)
          hist_DeltaR_quark_genjet->Fill(gen_quarks[i].DeltaR(gen_jets[j]), hweight);

    for(int i = 0; i < NQuarks; i++)
      for(int j = 0; j < Njets; j++)
          hist_DeltaR_quark_jet->Fill(gen_quarks[i].DeltaR(jets[j]), hweight);

    // gen jet matching
    for(int i = 0; i < int(allquarks.size()); i++){
      double minDR = 0.35;
      int index = -1;
      bool matched = false;
      for(int j = 0; j < Ngen_jets; j++){
        if(allquarks[i].DeltaR(gen_jets[j]) < minDR){
          index = j;
          minDR = allquarks[i].DeltaR(gen_jets[j]);
        }
      }
      if(index > -1){
        matched = true;
        gen_jets[index].SetMomPDGID(allquarks[i].MomPDGID());
        gen_jets[index].SetGenIndex(allquarks[i].GenIndex());
        gen_jets[index].SetGenMomIndex(allquarks[i].GenMomIndex());
      }
      eff_quarkMatchedGenJet_eta->Fill(matched, allquarks[i].Eta(), hweight);
    }

    ParticleList unmatchedGenJets;
    ParticleList matchedGenJets; // matched to jet coming from W
    ParticleList matchedBGenJets; // matched to b quark
    for(int i = 0; i < Ngen_jets; i++){
      if(abs(gen_jets[i].MomPDGID()) == 6)
        matchedBGenJets.push_back(gen_jets[i]);
      else if(abs(gen_jets[i].MomPDGID()) == 23)
        matchedGenJets.push_back(gen_jets[i]);
      else if(abs(gen_jets[i].MomPDGID()) == 24)
        matchedGenJets.push_back(gen_jets[i]);
      else
        unmatchedGenJets.push_back(gen_jets[i]);
    }

    // reco jet matching
    for(int i = 0; i < int(allquarks.size()); i++){
      double minDR = 0.35;
      int index = -1;
      bool matched = false;
      for(int j = 0; j < Njets; j++){
        if(allquarks[i].DeltaR(jets[j]) < minDR){
          index = j;
          minDR = allquarks[i].DeltaR(jets[j]);
        }
      }
      if(index > -1){
        matched = true;
        jets[index].SetMomPDGID(allquarks[i].MomPDGID());
        jets[index].SetGenIndex(allquarks[i].GenIndex());
        jets[index].SetGenMomIndex(allquarks[i].GenMomIndex());
      }
      eff_quarkMatchedGenJet_eta->Fill(matched, allquarks[i].Eta(), hweight);
    }

    ParticleList unmatchedJets;
    ParticleList matchedJets; // matched to jet coming from W
    ParticleList matchedBJets; // matched to b quark
    for(int i = 0; i < Njets; i++){
      if(abs(jets[i].MomPDGID()) == 6)
        matchedBJets.push_back(jets[i]);
      else if(abs(jets[i].MomPDGID()) == 23)
        matchedJets.push_back(jets[i]);
      else if(abs(jets[i].MomPDGID()) == 24)
        matchedJets.push_back(jets[i]);
      else
        unmatchedJets.push_back(jets[i]);
    }

    hist_NQuarks->Fill(NQuarks, hweight);

    TVector3 MET(0.,0.,0.);
    MET.SetPtEtaPhi(base->MET_pt,0.,base->MET_phi);
    MET.SetZ(0.);

    // binary decay tree with S & ISR splitting #2
    LAB.ClearEvent();
    vector<RFKey> jetID_BDT;
    vector<RFKey> lepID_BDT;
    vector<RFKey> objID_BDT;
    INV.SetLabFrameThreeVector(MET);

    // use reco jets for BDT
    if(!use_gen_jets){
      for(int i = 0; i < Njets; i++){
        RFKey key = COMB_J.AddLabFrameFourVector(jets[i]);
        jetID_BDT.push_back(key);
        objID_BDT.push_back(key);
      }
    }

    for(int i = 0; i < Nleps; i++){
      RFKey key = COMB_L.AddLabFrameFourVector(leptons[i]);
      lepID_BDT.push_back(key);
      objID_BDT.push_back(key);
    }
    if(!LAB.AnalyzeEvent()) cout << "Problem with RF ISR Sparticle Analyze Event #2 \n";
    TVector3 vPISR = S.GetFourVector(CM).Vect();
    double PTISR = S.GetTransverseFourVector(CM).Vect().Mag();
    TVector3 vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();
    double RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
    hist_RISR_PTISR->Fill(RISR, PTISR, hweight);
    if(RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut) continue;

    int Nleps_A = 0;
    int Nleps_B = 0;
    int Njets_S = 0;
    int Njets_ISR = 0;
    int Njets_S_matched = 0;
    int Njets_S_unmatched = 0;
    int Njets_Sa_matched = 0;
    int Njets_Sa_unmatched = 0;
    int Njets_Sb_matched = 0;
    int Njets_Sb_unmatched = 0;
    int Njets_ISR_matched = 0;
    int Njets_ISR_unmatched = 0;
    bool reco_sameside = false;
    for(int i = 0; i < Njets; i++){
      //if(COMB_J.GetFrame(jetID_BDT[i]) == J_X2a){
      if(COMB_J.GetFrame(jetID_BDT[i]) == J_X2a
      || COMB_J.GetFrame(jetID_BDT[i]) == J_X2b){
        Njets_S++;
        if(abs(jets[i].MomPDGID()) == 23 || abs(jets[i].MomPDGID()) == 24)
          Njets_S_matched++;
        else
          Njets_S_unmatched++;
      }
      else{
        Njets_ISR++;
        if(abs(jets[i].MomPDGID()) == 23 || abs(jets[i].MomPDGID()) == 24)
          Njets_ISR_matched++;
        else
          Njets_ISR_unmatched++;
      }
    }
    for(int i = 0; i < Nleps; i++){
      if(COMB_L.GetFrame(lepID_BDT[i]) == L_X2a)
        Nleps_A++;
      else
        Nleps_B++;
    }
    if(min_Sjet >= 0 && Njets_S < min_Sjet) continue;
    if(max_Sjet >= 0 && Njets_S > max_Sjet) continue;
    if(min_La >= 0 && Nleps_A < min_La) continue;
    if(max_La >= 0 && Nleps_A > max_La) continue;
    if(min_Lb >= 0 && Nleps_B < min_Lb) continue;
    if(max_Lb >= 0 && Nleps_B > max_Lb) continue;
    if(Nleps_A == 0 || Nleps_B == 0) reco_sameside = true;

    hist_matched_Njets_S_unmatched_Njets_S->Fill(Njets_S_matched, Njets_S_unmatched, hweight);
    hist_matched_Njets_ISR_unmatched_Njets_ISR->Fill(Njets_ISR_matched, Njets_ISR_unmatched, hweight);
    hist_matched_Njets_Sa_matched_Njets_Sb->Fill(Njets_Sa_matched, Njets_Sb_matched, hweight);
    hist_unmatched_Njets_Sa_unmatched_Njets_Sb->Fill(Njets_Sa_unmatched, Njets_Sb_unmatched, hweight);
    if(Nleps >= 2){
      TLorentzVector l1l2 = leptons[0]+leptons[1];
      double mllLEAD = l1l2.M();
      hist_RISR_mllLEAD->Fill(RISR, mllLEAD, hweight);
    }
    if(Nleps >= 2) hist_RISR_mL->Fill(RISR, (L_X2a.GetFourVector() + L_X2b.GetFourVector()).M(), hweight);
    hist_dphiCMI_PTCM->Fill(CM.GetDeltaPhiBoostVisible(), CM.GetFourVector().Pt(), hweight);
    hist_dphiMETV_PTISR->Fill(S.GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(MET), PTISR, hweight);
    double MVis = (J_X2a.GetFourVector() + J_X2b.GetFourVector() + L_X2a.GetFourVector() + L_X2b.GetFourVector()).M();
    double MLa = L_X2a.GetFourVector().Mag();
    double MLb = L_X2b.GetFourVector().Mag();
    double MJ =  (J_X2a.GetFourVector()+J_X2b.GetFourVector()).Mag();

    TLorentzVector vP_S_CM = S.GetFourVector(CM);
    TVector3 daBoost = vP_S_CM.Vect().Unit();

    TLorentzVector vP_Ja_CM = J_X2a.GetFourVector(CM);
    TLorentzVector vP_Jb_CM = J_X2b.GetFourVector(CM);
    TLorentzVector vP_La_CM = L_X2a.GetFourVector(CM);
    TLorentzVector vP_Lb_CM = L_X2b.GetFourVector(CM);
    TLorentzVector vP_Ia_CM = Ia.GetFourVector(CM);
    TLorentzVector vP_Ib_CM = Ib.GetFourVector(CM);

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

    double MSperpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_Jb_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV + vP_Ia_CM_PerpM0_TLV + vP_Ib_CM_PerpM0_TLV).Mag();
    double MaPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV + vP_Ia_CM_PerpM0_TLV).Mag();
    double MbPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV + vP_Ib_CM_PerpM0_TLV).Mag();
    double MaVPerpCM0 = (vP_Ja_CM_PerpM0_TLV + vP_La_CM_PerpM0_TLV).Mag();
    double MbVPerpCM0 = (vP_Jb_CM_PerpM0_TLV + vP_Lb_CM_PerpM0_TLV).Mag();
    double MQperpCM0 = sqrt((MaPerpCM0*MaPerpCM0 + MbPerpCM0*MbPerpCM0)/2.); 
    double gammaPerpCM0 = 2*MQperpCM0/MSperpCM0;
    double MVisAperpCM = (vP_Ja_CM_Perp_TLV + vP_La_CM_Perp_TLV).Mag();
    double MVisBperpCM = (vP_Jb_CM_Perp_TLV + vP_Lb_CM_Perp_TLV).Mag();

    hist_MSperpCM0_RISR->Fill(MSperpCM0, RISR, hweight);
    hist_MQperpCM0_RISR->Fill(MQperpCM0, RISR, hweight);
    hist_gammaPerpCM0_RISR->Fill(gammaPerpCM0, RISR, hweight);
    hist_MVisAperpCM_MVisBperpCM->Fill(MVisAperpCM, MVisBperpCM, hweight);

    hist_MSperpCM0_gammaPerpCM0->Fill(MSperpCM0, gammaPerpCM0, hweight);
    hist_MSperpCM0_MQperpCM0->Fill(MSperpCM0, MQperpCM0, hweight);
    hist_gammaPerpCM0_MQperpCM0->Fill(gammaPerpCM0, MQperpCM0, hweight);
    hist_MQperpCM0_MVisAperpCM->Fill(MQperpCM0, MVisAperpCM, hweight);
    hist_MQperpCM0_MVisBperpCM->Fill(MQperpCM0, MVisBperpCM, hweight);

    hist_MJ_MQperpCM0->Fill(MJ, MQperpCM0, hweight);
    hist_MLa_MQperpCM0->Fill(MLa, MQperpCM0, hweight);
    hist_MLb_MQperpCM0->Fill(MLb, MQperpCM0, hweight);
    hist_MJ_MSperpCM0->Fill(MJ, MSperpCM0, hweight);
    hist_MLa_MSperpCM0->Fill(MLa, MSperpCM0, hweight);
    hist_MLb_MSperpCM0->Fill(MLb, MSperpCM0, hweight);
    hist_MJ_gammaPerpCM0->Fill(MJ, gammaPerpCM0, hweight);
    hist_MLa_gammaPerpCM0->Fill(MLa, gammaPerpCM0, hweight);
    hist_MLb_gammaPerpCM0->Fill(MLb, gammaPerpCM0, hweight);

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

    double MSCM0 = (vP_Ja_CM0 + vP_Jb_CM0 + vP_La_CM0 + vP_Lb_CM0 + vP_Ia_CM0 + vP_Ib_CM0).Mag();
    double MaCM0 = (vP_Ja_CM0 + vP_La_CM0 + vP_Ia_CM0).Mag();
    double MbCM0 = (vP_Jb_CM0 + vP_Lb_CM0 + vP_Ib_CM0).Mag();
    double MaVCM0 = (vP_Ja_CM0 + vP_La_CM0).Mag();
    double MbVCM0 = (vP_Jb_CM0 + vP_Lb_CM0).Mag();
    double MQCM0 = sqrt((MaCM0*MaCM0 + MbCM0*MbCM0)/2.); 
    double gammaCM0 = 2*MQCM0/MSCM0;
    double MVisACM = (vP_Ja_CM + vP_La_CM).Mag();
    double MVisBCM = (vP_Jb_CM + vP_Lb_CM).Mag();

    hist_MSCM0_RISR->Fill(MSCM0, RISR, hweight);
    hist_MQCM0_RISR->Fill(MQCM0, RISR, hweight);
    hist_gammaCM0_RISR->Fill(gammaCM0, RISR, hweight);
    hist_MVisACM_MVisBCM->Fill(MVisACM, MVisBCM, hweight);

    hist_MSCM0_gammaCM0->Fill(MSCM0, gammaCM0, hweight);
    hist_MSCM0_MQCM0->Fill(MSCM0, MQCM0, hweight);
    hist_gammaCM0_MQCM0->Fill(gammaCM0, MQCM0, hweight);
    hist_MQCM0_MVisACM->Fill(MQCM0, MVisACM, hweight);
    hist_MQCM0_MVisBCM->Fill(MQCM0, MVisBCM, hweight);

    hist_MJ_MQCM0->Fill(MJ, MQCM0, hweight);
    hist_MLa_MQCM0->Fill(MLa, MQCM0, hweight);
    hist_MLb_MQCM0->Fill(MLb, MQCM0, hweight);
    hist_MJ_MSCM0->Fill(MJ, MSCM0, hweight);
    hist_MLa_MSCM0->Fill(MLa, MSCM0, hweight);
    hist_MLb_MSCM0->Fill(MLb, MSCM0, hweight);
    hist_MJ_gammaCM0->Fill(MJ, gammaCM0, hweight);
    hist_MLa_gammaCM0->Fill(MLa, gammaCM0, hweight);
    hist_MLb_gammaCM0->Fill(MLb, gammaCM0, hweight);

    hist_MLa_MJ->Fill(MLa, MJ, hweight);
    hist_MLb_MJ->Fill(MLb, MJ, hweight);
    hist_MLa_MLb->Fill(MLa, MLb, hweight);
    hist_RISR_MJ->Fill(RISR, MJ, hweight);
    hist_RISR_MLa->Fill(RISR, MLa, hweight);
    hist_RISR_MLb->Fill(RISR, MLb, hweight);

    TLorentzVector vP_Ja_S = J_X2a.GetFourVector(S);
    TLorentzVector vP_Jb_S = J_X2b.GetFourVector(S);
    TLorentzVector vP_La_S = L_X2a.GetFourVector(S);
    TLorentzVector vP_Lb_S = L_X2b.GetFourVector(S);
    TLorentzVector vP_Ia_S = Ia.GetFourVector(S);
    TLorentzVector vP_Ib_S = Ib.GetFourVector(S);
    TVector3 boostVis = (vP_Ja_S+vP_Jb_S+vP_La_S+vP_Lb_S).BoostVector();
    TVector3 boostInv = (vP_Ia_S+vP_Ib_S).BoostVector();
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
    double PX3_BoostT = (vP_Ja_S+vP_Jb_S+vP_La_S+vP_Ia_S).P();
    double MX3a_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).M();
    double MX3b_BoostT = (vP_Lb_S+vP_Jb_S+vP_Ib_S).M();
    double MQperp = sqrt(MX3a_BoostT*MX3a_BoostT+MX3b_BoostT*MX3b_BoostT)/sqrt(2.);
    double gammaPerp = 2*MQperp/(sqrt(MX3a_BoostT*MX3a_BoostT+PX3_BoostT*PX3_BoostT) +
		                sqrt(MX3b_BoostT*MX3b_BoostT+PX3_BoostT*PX3_BoostT));
    hist_MQperp_RISR->Fill(MQperp, RISR, hweight);
    hist_gammaPerp_RISR->Fill(gammaPerp, RISR, hweight);

    int EC_X = 0; // root hists have underflow in bin 0
    int cat_Nleps = Nleps;
    int cat_Njets_S = Njets_S;
    if(cat_Nleps > 4) cat_Nleps = 4; // ge4L is upper limit
    if(cat_Njets_S > 2) cat_Njets_S = 2; // ge2J is upper limit
    EC_X += 3*cat_Nleps-5;
    EC_X += cat_Njets_S;
    int EC_Y = 0; // root hists have underflow in bin 0
    if(proc.find("ttbar") != std::string::npos) EC_Y = 1;
    if(proc.find("DB") != std::string::npos) EC_Y = 2;
    if(proc.find("WJets") != std::string::npos) EC_Y = 3;
    if(proc.find("ZDY") != std::string::npos) EC_Y = 4;
    int Nbkg_proc = 4;
    if(proc.find("Cascades") != std::string::npos){
        EC_Y = Nbkg_proc+1 + slep;
      if(proc.find("_10") != std::string::npos)
        EC_Y += 3;
      if(proc.find("_5") != std::string::npos)
        EC_Y += 6;
    }
    hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+1);
    hist_EventCount->SetBinContent(0,EC_Y,hist_EventCount->GetBinContent(0,EC_Y)+1);

    kept_file_events++;
    hist_MET->Fill(base->MET_pt, hweight);
  }
  
  string output_root_filename = "output_Plot_1D_NANO_"+ofile;
  TFile* output_file = new TFile((output_root_filename).c_str(),"RECREATE");
  output_file->mkdir(plot_folder.c_str());
  output_file->cd(plot_folder.c_str());
  for(int hist1 = 0; hist1 < int(hists1.size()); hist1++) hists1[hist1]->Write();
  for(int hist2 = 0; hist2 < int(hists2.size()); hist2++) hists2[hist2]->Write();
  for(int eff = 0; eff < int(effs.size()); eff++) effs[eff]->Write();
  hist_EventCount->Write();
  output_file->Close(); 
  delete base;
  delete chain;
}

string load_file_from_list(string fileListName, int line_number)
{
 std::ifstream infile(fileListName);
 string line = "";
 int counter = 0;
 while(getline(infile,line))
 {
  if(counter == line_number) return line;
  counter++;
 }
 return line;
}
