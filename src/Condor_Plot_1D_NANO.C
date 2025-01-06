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
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include "SUSYNANOBase.hh"
#include "V_Cand.hh"
#include "CategoryTool.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

bool use_prong2 = false;
bool use_aunt_prong = false;
bool use_lep_prong = false;
double min_Cand_CosDecayAngle_CM = -1.; // default: -1
double max_Cand_CosDecayAngle_CM = 1.; // default: 1.

string load_file_from_list(string fileListName, int line_number);
string double_to_string(double val);
int FrameDistance(const RestFrame& frame, const RestFrame& sibling);
double getX2(const double& pMass, const double& mass, const double& width);
void getJetsMinX2(int& j1_index, int& j2_index, ParticleList& jets);
template <typename V>
bool inVec(const std::vector<V>& vect, const V& value){
  return std::find(vect.begin(), vect.end(), value) != vect.end();
}
std::vector<V_Cand> cand_list(const ParticleList& jets, const ParticleList& leps,
  const CombinatoricGroup& Comb, const SelfAssemblingRecoFrame& BDT, const VisibleRecoFrame& VIS, const RestFrame& CM,
  const std::vector<RFKey>& jetID, const std::vector<RFKey>& lepID);
void cand_matching(std::vector<V_Cand>& cand_list);
void cand_side(std::vector<V_Cand>& cand_list, const SelfAssemblingRecoFrame& BDT);

int main(int argc, char* argv[]) {

  string plot_folder = "";
  bool use_gen_jets = false;
  bool boson_acceptance_cut = false;
  bool gen_lepton_cut = false;
  int reco_lep_cut = 0;
  int min_Sjet = -1;
  int max_Sjet = -1;
  int min_La = -1;
  int max_La = -1;
  int min_Lb = -1;
  int max_Lb = -1;
  bool KIN = false; // kinematic cuts like PTISR, RISR, and MET
  double MET_cut = 0.;
  double PTISR_cut = 0.;
  double RISR_cut = 0.;

  bool Zero_RJR_BDT_JetMass = false; // set mass of jets to zero before adding into RJR tree
  bool BDT_MET = false; // whether MET is in BDT
  bool S_MET = false; // whether MET is in S system
  bool transverseOBJ = false; // whether or not only transverse component of objs should be passed to BDT
  bool transverseJIG = false; // whether or not jigsaw that splits S & ISR should be transverse
  bool SqJIG = false; // whether or not jigsaw that splits S & ISR should be squared masses (for using leptons to decided on S/ISR splitting)
  bool S_leps = false; // whether or not to force leps into Sparticle system
  bool Sparticle1 = false; // whether or not to save hists for "ISR_Sparticle" RJR tree
  bool Sparticle2 = false; // whether or not to save hists for "ISR_Sparticle2" RJR tree
  bool HEM = false; // whether or not to "force" HEM splitting in S system
  bool Cand = false; // whether or not to do things with BDT and Cands
  double min_Cand_CosDecayAngle_CM = -1.; // default: -1
  double max_Cand_CosDecayAngle_CM = 1.; // default: 1.

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
    if(strncmp(argv[i],"--KIN", 5) == 0){
      KIN = true;
      if(debug) cout << "using KIN cuts" << endl;
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
    if(strncmp(argv[i],"--ZeroJetMass", 13) == 0){
      Zero_RJR_BDT_JetMass = true;
    }
    if(strncmp(argv[i],"--BDT_MET", 9) == 0){
      BDT_MET = true;
      if(debug) cout << "turning on BDT_MET" << endl;
    }
    if(strncmp(argv[i],"--S_MET", 7) == 0){
      S_MET = true;
      if(debug) cout << "turning on S_MET" << endl;
    }
    if(strncmp(argv[i],"--transverseJIG", 15) == 0){
      transverseJIG = true;
      if(debug) cout << "turning on transverseJIG" << endl;
    }
    if(strncmp(argv[i],"--transverseOBJ", 15) == 0){
      transverseOBJ = true;
      if(debug) cout << "turning on transverseOBJ" << endl;
    }
    if(strncmp(argv[i],"--SqJIG", 7) == 0){
      SqJIG = true;
      if(debug) cout << "turning on SqJIG" << endl;
    }
    if(strncmp(argv[i],"--S_leps", 8) == 0){
      S_leps = true;
      if(debug) cout << "turning on S_leps" << endl;
    }
    if(strncmp(argv[i],"--Sparticle1", 12) == 0){
      Sparticle1 = true;
    }
    if(strncmp(argv[i],"--Sparticle2", 12) == 0){
      Sparticle2 = true;
    }
    if(strncmp(argv[i],"--Cand", 6) == 0){
      Cand = true;
    }
    if(strncmp(argv[i],"--HEM", 5) == 0){
      HEM = true;
      if(debug) cout << "Forcing hemisphere split in S system" << endl;
    }
    if(strncmp(argv[i],"--reco_lep_cut", 14) == 0){
      i++;
      sscanf(argv[i],"%d,%d", &reco_lep_cut);
      if(debug) cout << "Requiring " << reco_lep_cut << " reco leps" << endl;
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
  if(Zero_RJR_BDT_JetMass) plot_folder += "_ZeroJetMass";
  if(BDT_MET) plot_folder += "_BDTMET";
  if(S_MET) plot_folder += "_SMET";
  if(transverseOBJ) plot_folder += "_transverseOBJ";
  if(transverseJIG) plot_folder += "_transverseJIG";
  if(SqJIG) plot_folder += "_SqJIG";
  if(KIN) plot_folder += "_KIN";
  if(use_prong2) plot_folder += "_prong2";
  if(use_aunt_prong) plot_folder += "_auntprong";
  if(use_lep_prong) plot_folder += "_lepprong";
  if(boson_acceptance_cut) plot_folder += "_qqbar";
  if(HEM) plot_folder += "_HEM";
  else plot_folder += "_XX";
  if(gen_lepton_cut) plot_folder += "_l";
  else plot_folder += "_XX";
  if(min_Cand_CosDecayAngle_CM > -1. || max_Cand_CosDecayAngle_CM < 1.){
    plot_folder += "_CosThetaCMCut";
    if(min_Cand_CosDecayAngle_CM > -1.){
      plot_folder += "_G"+double_to_string(min_Cand_CosDecayAngle_CM);
    }
    if(max_Cand_CosDecayAngle_CM < 1.){
      plot_folder += "_L"+double_to_string(max_Cand_CosDecayAngle_CM);
    }
  }

  // BDT
  LabRecoFrame LAB_BDT("LAB_BDT","lab");
  SelfAssemblingRecoFrame BDT_CM("BDT_CM","CM");
  VisibleRecoFrame V1("V1","V");
  CombinatoricGroup VIS("VIS","Visible Object Jigsaws");
  VIS.AddFrame(V1);
  VIS.SetNElementsForFrame(V1,1,false);
  InvisibleRecoFrame I1("I1","#vec{E}_{T}^{miss}");
  InvisibleGroup INV_BDT("INV_BDT","Invisible State Jigsaws");
  INV_BDT.AddFrame(I1);
  LAB_BDT.SetChildFrame(BDT_CM);
  if(BDT_MET)
    BDT_CM.AddChildFrame(I1);
  BDT_CM.AddChildFrame(V1);
  LAB_BDT.InitializeTree();
  SetMassInvJigsaw MinMassJigsaw_BDT("MINMASS_JIGSAW_BDT", "Invisible system mass Jigsaw");
  SetRapidityInvJigsaw RapidityJigsaw_BDT("RAPIDITY_JIGSAW_BDT", "Invisible system rapidity Jigsaw");
  if(BDT_MET){
    INV_BDT.AddJigsaw(MinMassJigsaw_BDT);
    INV_BDT.AddJigsaw(RapidityJigsaw_BDT);
    RapidityJigsaw_BDT.AddVisibleFrames((LAB_BDT.GetListVisibleFrames()));
  }
  LAB_BDT.InitializeAnalysis();
  if(treeplot){
    tree_plot.SetTree(LAB_BDT);
    tree_plot.Draw("BDT_tree", "BDT Tree", treeplot_INV);
    if(BDT_MET){
      tree_plot.SetTree(INV_BDT);
      tree_plot.Draw("BDT_inv", "Invisible Jigsaws", treeplot_INV);
    }
  }

  // BDT ISR Sparticle Splitting
  LabRecoFrame LAB_ISR_Sparticle("LAB_ISR_Sparticle","lab");
  DecayRecoFrame CM_ISR_Sparticle("CM_ISR_Sparticle","CM");
  DecayRecoFrame S_ISR_Sparticle("S_ISR_Sparticle","S");
  DecayRecoFrame Sa_ISR_Sparticle("Sa_ISR_Sparticle","Sa");
  DecayRecoFrame Sb_ISR_Sparticle("Sb_ISR_Sparticle","Sb");
  SelfAssemblingRecoFrame BDT_S_ISR_Sparticle("BDT_S_ISR_Sparticle","S");
  SelfAssemblingRecoFrame BDT_Sa_ISR_Sparticle("BDT_S_ISR_Sparticle","Sa");
  SelfAssemblingRecoFrame BDT_Sb_ISR_Sparticle("BDT_S_ISR_Sparticle","Sb");
  SelfAssemblingRecoFrame BDT_ISR_ISR_Sparticle("BDT_ISR_ISR_Sparticle","ISR");
  VisibleRecoFrame V_S_ISR_Sparticle("V_S_ISR_Sparticle","V");
  VisibleRecoFrame V_Sa_ISR_Sparticle("V_Sa_ISR_Sparticle","V_{a}");
  VisibleRecoFrame V_Sb_ISR_Sparticle("V_Sb_ISR_Sparticle","V_{b}");
  VisibleRecoFrame V_ISR_ISR_Sparticle("V_ISR_ISR_Sparticle","V_{ISR}");
  VisibleRecoFrame J_S_ISR_Sparticle("J_S_ISR_Sparticle","J");
  VisibleRecoFrame L_S_ISR_Sparticle("L_S_ISR_Sparticle","L");
  VisibleRecoFrame J_Sa_ISR_Sparticle("J_Sa_ISR_Sparticle","J_{a}");
  VisibleRecoFrame L_Sa_ISR_Sparticle("L_Sa_ISR_Sparticle","L_{a}");
  VisibleRecoFrame J_Sb_ISR_Sparticle("J_Sb_ISR_Sparticle","J_{b}");
  VisibleRecoFrame L_Sb_ISR_Sparticle("L_Sb_ISR_Sparticle","L_{b}");
  InvisibleRecoFrame I_ISR_Sparticle("I_ISR_Sparticle","#vec{E}_{T}^{miss}");
  InvisibleRecoFrame Ia_ISR_Sparticle("Ia_ISR_Sparticle","#vec{Ea}_{T}^{miss}");
  InvisibleRecoFrame Ib_ISR_Sparticle("Ib_ISR_Sparticle","#vec{Eb}_{T}^{miss}");

  InvisibleGroup INV_ISR_Sparticle("INV_ISR_Sparticle","Invisible State Jigsaws");

  LAB_ISR_Sparticle.SetChildFrame(CM_ISR_Sparticle);
  CM_ISR_Sparticle.AddChildFrame(S_ISR_Sparticle);
  CM_ISR_Sparticle.AddChildFrame(BDT_ISR_ISR_Sparticle);
  BDT_ISR_ISR_Sparticle.AddChildFrame(V_ISR_ISR_Sparticle);
  if(!HEM){
    if(BDT_MET || S_MET){
      INV_ISR_Sparticle.AddFrame(I_ISR_Sparticle);
    }
    S_ISR_Sparticle.AddChildFrame(BDT_S_ISR_Sparticle);
    if(BDT_MET){
      BDT_S_ISR_Sparticle.AddChildFrame(I_ISR_Sparticle);
    }
    else if(S_MET){
      S_ISR_Sparticle.AddChildFrame(I_ISR_Sparticle);
    }
    if(S_leps){
      BDT_S_ISR_Sparticle.AddChildFrame(J_S_ISR_Sparticle);
      BDT_S_ISR_Sparticle.AddChildFrame(L_S_ISR_Sparticle);
    }
    else{
      BDT_S_ISR_Sparticle.AddChildFrame(V_S_ISR_Sparticle);
    }
  }
  else if(HEM){
    if(BDT_MET || S_MET){
      INV_ISR_Sparticle.AddFrame(Ia_ISR_Sparticle);
      INV_ISR_Sparticle.AddFrame(Ib_ISR_Sparticle);
    }
    S_ISR_Sparticle.AddChildFrame(Sa_ISR_Sparticle);
    S_ISR_Sparticle.AddChildFrame(Sb_ISR_Sparticle);
    Sa_ISR_Sparticle.AddChildFrame(BDT_Sa_ISR_Sparticle);
    Sb_ISR_Sparticle.AddChildFrame(BDT_Sb_ISR_Sparticle);
    if(BDT_MET){
      BDT_Sa_ISR_Sparticle.AddChildFrame(Ia_ISR_Sparticle);
      BDT_Sb_ISR_Sparticle.AddChildFrame(Ib_ISR_Sparticle);
    }
    else if(S_MET){
      Sa_ISR_Sparticle.AddChildFrame(Ia_ISR_Sparticle);
      Sb_ISR_Sparticle.AddChildFrame(Ib_ISR_Sparticle);
    }
    if(S_leps){
      BDT_Sa_ISR_Sparticle.AddChildFrame(J_Sa_ISR_Sparticle);
      BDT_Sa_ISR_Sparticle.AddChildFrame(L_Sa_ISR_Sparticle);
      BDT_Sb_ISR_Sparticle.AddChildFrame(J_Sb_ISR_Sparticle);
      BDT_Sb_ISR_Sparticle.AddChildFrame(L_Sb_ISR_Sparticle);
    }
    else{
      BDT_Sa_ISR_Sparticle.AddChildFrame(V_Sa_ISR_Sparticle);
      BDT_Sb_ISR_Sparticle.AddChildFrame(V_Sb_ISR_Sparticle);
    }
  }
  LAB_ISR_Sparticle.InitializeTree();

  CombinatoricGroup VIS_ISR_Sparticle("VIS","Visible Object Jigsaws");
  CombinatoricGroup COMB_J_ISR_Sparticle("VIS","Visible Object Jigsaws");
  CombinatoricGroup COMB_L_ISR_Sparticle("VIS","Visible Object Jigsaws");
  SetMassInvJigsaw MinMassJigsaw_BDT_ISR_Sparticle("MINMASS_JIGSAW_BDT", "Invisible system mass Jigsaw");
  SetRapidityInvJigsaw RapidityJigsaw_BDT_ISR_Sparticle("RAPIDITY_JIGSAW_BDT", "Invisible system rapidity Jigsaw");
  MinMassesCombJigsaw CombSplit_ISR_ISR_Sparticle_Transverse("CombSplit_ISR", "Minimize M_{T}^{ISR} and M_{T}^{S}");
  MinMassesCombJigsaw CombSplit_ISR_ISR_Sparticle("CombSplit_ISR", "Minimize M^{ISR} and M^{S}");
  MinMassesSqCombJigsaw CombSplitSq_ISR_ISR_Sparticle("CombSplitSq_ISR", "Minimize M_{ISR}^{2} + M_{S}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_ISR_Sparticle("CombSplitSq_S_V", "Minimize M_{Sa}^{2} + M_{Sb}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_L_ISR_Sparticle("CombSplitSq_S_L", "Minimize M_{Sa}^{2} + M_{Sb}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_J_ISR_Sparticle("CombSplitSq_S_J", "Minimize M_{Sa}^{2} + M_{Sb}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_ISR_ISR_Sparticle_Transverse("CombSplitSq_ISR", "Minimize M_{T}_{ISR}^{2} + M_{T}_{S}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_ISR_Sparticle_Transverse("CombSplitSq_S_V", "Minimize M_{T}_{Sa}^{2} + M_{T}_{Sb}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_L_ISR_Sparticle_Transverse("CombSplitSq_S_L", "Minimize M_{T}_{Sa}^{2} + M_{T}_{Sb}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_S_J_ISR_Sparticle_Transverse("CombSplitSq_S_J", "Minimize M_{T}_{Sa}^{2} + M_{T}_{Sb}^{2}",2,2);
  MinMassesSqInvJigsaw InvSplit_ISR_Sparticle("InvSplit_ISR_Sparticle", "INV -> I_{a} + I_{b}", 2);

  if(BDT_MET || S_MET){
    INV_ISR_Sparticle.AddJigsaw(MinMassJigsaw_BDT_ISR_Sparticle);
    INV_ISR_Sparticle.AddJigsaw(RapidityJigsaw_BDT_ISR_Sparticle);
    RapidityJigsaw_BDT_ISR_Sparticle.AddVisibleFrames(S_ISR_Sparticle.GetListVisibleFrames());
  }

  if(!HEM){
    if(S_leps){
      COMB_L_ISR_Sparticle.AddFrame(L_S_ISR_Sparticle);
      COMB_L_ISR_Sparticle.SetNElementsForFrame(L_S_ISR_Sparticle, 0);
      COMB_J_ISR_Sparticle.AddFrame(V_ISR_ISR_Sparticle);
      COMB_J_ISR_Sparticle.SetNElementsForFrame(V_ISR_ISR_Sparticle, 1);
      COMB_J_ISR_Sparticle.AddFrame(J_S_ISR_Sparticle);
      COMB_J_ISR_Sparticle.SetNElementsForFrame(J_S_ISR_Sparticle, 0);
      if(SqJIG){
        if(transverseJIG){
          COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_ISR_ISR_Sparticle_Transverse);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.SetTransverse();
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_S_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
        }
        else{
          COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_ISR_ISR_Sparticle);
          CombSplitSq_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle.AddCombFrame(J_S_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
        }
      } 
      else if(transverseJIG){
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle_Transverse);
        CombSplit_ISR_ISR_Sparticle.SetTransverse();
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_S_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      else{
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(J_S_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
      }
    }
    else{
      VIS_ISR_Sparticle.AddFrame(V_ISR_ISR_Sparticle);
      VIS_ISR_Sparticle.SetNElementsForFrame(V_ISR_ISR_Sparticle, 1);
      VIS_ISR_Sparticle.AddFrame(V_S_ISR_Sparticle);
      VIS_ISR_Sparticle.SetNElementsForFrame(V_S_ISR_Sparticle, 0);
      if(transverseJIG){
        VIS_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle_Transverse);
        CombSplit_ISR_ISR_Sparticle_Transverse.SetTransverse();
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_S_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      else{
        VIS_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_S_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
      }
    }
  } // if(!HEM)
  if(HEM){
    INV_ISR_Sparticle.AddJigsaw(InvSplit_ISR_Sparticle);
    InvSplit_ISR_Sparticle.AddInvisibleFrame(Ia_ISR_Sparticle, 0);
    InvSplit_ISR_Sparticle.AddInvisibleFrame(Ib_ISR_Sparticle, 1);
    if(S_leps){
      InvSplit_ISR_Sparticle.AddVisibleFrame(J_Sa_ISR_Sparticle, 0);
      InvSplit_ISR_Sparticle.AddVisibleFrame(L_Sa_ISR_Sparticle, 0);
      InvSplit_ISR_Sparticle.AddVisibleFrame(J_Sb_ISR_Sparticle, 1);
      InvSplit_ISR_Sparticle.AddVisibleFrame(L_Sb_ISR_Sparticle, 1);
    }
    else{
      InvSplit_ISR_Sparticle.AddVisibleFrame(V_Sa_ISR_Sparticle, 0);
      InvSplit_ISR_Sparticle.AddVisibleFrame(V_Sb_ISR_Sparticle, 1);
    }
    if(S_leps){
      // jigsaw(s) for splitting ISR & S system
      COMB_L_ISR_Sparticle.AddFrame(L_Sa_ISR_Sparticle);
      COMB_L_ISR_Sparticle.SetNElementsForFrame(L_Sa_ISR_Sparticle, 0);
      COMB_L_ISR_Sparticle.AddFrame(L_Sb_ISR_Sparticle);
      COMB_L_ISR_Sparticle.SetNElementsForFrame(L_Sb_ISR_Sparticle, 0);
      COMB_J_ISR_Sparticle.AddFrame(J_Sa_ISR_Sparticle);
      COMB_J_ISR_Sparticle.SetNElementsForFrame(J_Sa_ISR_Sparticle, 0);
      COMB_J_ISR_Sparticle.AddFrame(J_Sb_ISR_Sparticle);
      COMB_J_ISR_Sparticle.SetNElementsForFrame(J_Sb_ISR_Sparticle, 0);
      COMB_J_ISR_Sparticle.AddFrame(V_ISR_ISR_Sparticle);
      COMB_J_ISR_Sparticle.SetNElementsForFrame(V_ISR_ISR_Sparticle, 1);
      if(SqJIG){
        if(S_MET){
          COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_ISR_ISR_Sparticle_Transverse);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.SetTransverse();
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_Sa_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_Sb_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
        }
        else{
          COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_ISR_ISR_Sparticle);
          CombSplitSq_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle.AddCombFrame(J_Sa_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle.AddCombFrame(J_Sb_ISR_Sparticle, 1);
          CombSplitSq_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
          CombSplitSq_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
        }
      } 
      else if(transverseJIG){
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle_Transverse);
        CombSplit_ISR_ISR_Sparticle.SetTransverse();
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_Sa_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(J_Sb_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      else{
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(J_Sa_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(J_Sb_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      // jigsaw(s) for splitting S system
      if(transverseJIG){
        COMB_L_ISR_Sparticle.AddJigsaw(CombSplitSq_S_L_ISR_Sparticle_Transverse);
        CombSplitSq_S_L_ISR_Sparticle_Transverse.SetTransverse();
        CombSplitSq_S_L_ISR_Sparticle_Transverse.AddCombFrame(L_Sa_ISR_Sparticle, 0);
        CombSplitSq_S_L_ISR_Sparticle_Transverse.AddCombFrame(L_Sb_ISR_Sparticle, 1);
        CombSplitSq_S_L_ISR_Sparticle_Transverse.AddObjectFrame(Sa_ISR_Sparticle, 0);
        CombSplitSq_S_L_ISR_Sparticle_Transverse.AddObjectFrame(Sb_ISR_Sparticle, 1);
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_S_J_ISR_Sparticle_Transverse);
        CombSplitSq_S_J_ISR_Sparticle_Transverse.SetTransverse();
        CombSplitSq_S_J_ISR_Sparticle_Transverse.AddCombFrame(J_Sa_ISR_Sparticle, 0);
        CombSplitSq_S_J_ISR_Sparticle_Transverse.AddCombFrame(J_Sb_ISR_Sparticle, 1);
        CombSplitSq_S_J_ISR_Sparticle_Transverse.AddObjectFrame(Sa_ISR_Sparticle, 0);
        CombSplitSq_S_J_ISR_Sparticle_Transverse.AddObjectFrame(Sb_ISR_Sparticle, 1);

      }
      else{
        COMB_L_ISR_Sparticle.AddJigsaw(CombSplitSq_S_L_ISR_Sparticle);
        CombSplitSq_S_L_ISR_Sparticle.AddCombFrame(L_Sa_ISR_Sparticle, 0);
        CombSplitSq_S_L_ISR_Sparticle.AddCombFrame(L_Sb_ISR_Sparticle, 1);
        CombSplitSq_S_L_ISR_Sparticle.AddObjectFrame(Sa_ISR_Sparticle, 0);
        CombSplitSq_S_L_ISR_Sparticle.AddObjectFrame(Sb_ISR_Sparticle, 1);
        COMB_J_ISR_Sparticle.AddJigsaw(CombSplitSq_S_J_ISR_Sparticle);
        CombSplitSq_S_J_ISR_Sparticle.AddCombFrame(J_Sa_ISR_Sparticle, 0);
        CombSplitSq_S_J_ISR_Sparticle.AddCombFrame(J_Sb_ISR_Sparticle, 1);
        CombSplitSq_S_J_ISR_Sparticle.AddObjectFrame(Sa_ISR_Sparticle, 0);
        CombSplitSq_S_J_ISR_Sparticle.AddObjectFrame(Sb_ISR_Sparticle, 1);
      }
    } // if(S_leps)
    else{
      // jigsaw(s) for splitting ISR & S systems
      VIS_ISR_Sparticle.AddFrame(V_ISR_ISR_Sparticle);
      VIS_ISR_Sparticle.SetNElementsForFrame(V_ISR_ISR_Sparticle, 1);
      VIS_ISR_Sparticle.AddFrame(V_Sa_ISR_Sparticle);
      VIS_ISR_Sparticle.SetNElementsForFrame(V_Sa_ISR_Sparticle, 0);
      VIS_ISR_Sparticle.AddFrame(V_Sb_ISR_Sparticle);
      VIS_ISR_Sparticle.SetNElementsForFrame(V_Sb_ISR_Sparticle, 0);
      if(transverseJIG){
        VIS_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle_Transverse);
        CombSplit_ISR_ISR_Sparticle_Transverse.SetTransverse();
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_Sa_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddCombFrame(V_Sb_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle_Transverse.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      else{
        VIS_ISR_Sparticle.AddJigsaw(CombSplit_ISR_ISR_Sparticle);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_Sa_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddCombFrame(V_Sb_ISR_Sparticle, 1);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(BDT_ISR_ISR_Sparticle, 0);
        CombSplit_ISR_ISR_Sparticle.AddObjectFrame(S_ISR_Sparticle, 1);
      }
      // jigsaw(s) for splitting S system
      VIS_ISR_Sparticle.AddJigsaw(CombSplitSq_S_ISR_Sparticle);
      CombSplitSq_S_ISR_Sparticle.AddCombFrame(V_Sa_ISR_Sparticle, 0);
      CombSplitSq_S_ISR_Sparticle.AddCombFrame(V_Sb_ISR_Sparticle, 1);
      CombSplitSq_S_ISR_Sparticle.AddObjectFrames(Sa_ISR_Sparticle.GetListVisibleFrames(), 0);
      CombSplitSq_S_ISR_Sparticle.AddObjectFrames(Sb_ISR_Sparticle.GetListVisibleFrames(), 1);
    } // not S_leps
  }

  if(!LAB_ISR_Sparticle.InitializeAnalysis()) cout << "Failed to initialize analysis of LAB_ISR_Sparticle!" << endl;
  if(treeplot){
    tree_plot.SetTree(LAB_ISR_Sparticle);
    tree_plot.Draw("ISR_Sparticle_tree", "Analysis Tree", treeplot_INV);
    if(S_MET){
      tree_plot.SetTree(INV_ISR_Sparticle);
      tree_plot.Draw("ISR_Sparticle_inv", "Invisible Jigsaws", treeplot_INV);
    }
    if(S_leps){
      tree_plot.SetTree(COMB_J_ISR_Sparticle);
      tree_plot.Draw("ISR_Sparticle_jets", "Jet Jigsaws", treeplot_INV);
      tree_plot.SetTree(COMB_L_ISR_Sparticle);
      tree_plot.Draw("ISR_Sparticle_leps", "Lep Jigsaws", treeplot_INV);
    }
    else{
      tree_plot.SetTree(VIS_ISR_Sparticle);
      tree_plot.Draw("ISR_Sparticle_vis", "Vis Jigsaws", treeplot_INV);
    }
  }

  // BDT ISR Sparticle Splitting #2
  LabRecoFrame LAB_ISR_Sparticle2("LAB_ISR_Sparticle2","lab");
  DecayRecoFrame CM_ISR_Sparticle2("CM_ISR_Sparticle2","CM");
  DecayRecoFrame S_ISR_Sparticle2("S_ISR_Sparticle2","S");
  DecayRecoFrame X2a_ISR_Sparticle2("X2a_ISR_Sparticle2","X2a");
  DecayRecoFrame X2b_ISR_Sparticle2("X2b_ISR_Sparticle2","X2b");
  VisibleRecoFrame J_X2a_ISR_Sparticle2("J_X2a_ISR_Sparticle2","J_{a}");
  VisibleRecoFrame L_X2a_ISR_Sparticle2("L_X2a_ISR_Sparticle2","L_{a}");
  VisibleRecoFrame L_X2b_ISR_Sparticle2("L_X2b_ISR_Sparticle2","L_{b}");
  VisibleRecoFrame ISR_ISR_Sparticle2("ISR_ISR_Sparticle2","ISR");
  InvisibleRecoFrame Ia_ISR_Sparticle2("Ia_ISR_Sparticle2","Ia");
  InvisibleRecoFrame Ib_ISR_Sparticle2("Ib_ISR_Sparticle2","Ib");

  LAB_ISR_Sparticle2.SetChildFrame(CM_ISR_Sparticle2);
  CM_ISR_Sparticle2.AddChildFrame(S_ISR_Sparticle2);
  CM_ISR_Sparticle2.AddChildFrame(ISR_ISR_Sparticle2);
  S_ISR_Sparticle2.AddChildFrame(X2a_ISR_Sparticle2);
  S_ISR_Sparticle2.AddChildFrame(X2b_ISR_Sparticle2);
  X2a_ISR_Sparticle2.AddChildFrame(J_X2a_ISR_Sparticle2);
  X2a_ISR_Sparticle2.AddChildFrame(L_X2a_ISR_Sparticle2);
  X2a_ISR_Sparticle2.AddChildFrame(Ia_ISR_Sparticle2);
  X2b_ISR_Sparticle2.AddChildFrame(L_X2b_ISR_Sparticle2);
  X2b_ISR_Sparticle2.AddChildFrame(Ib_ISR_Sparticle2);

  LAB_ISR_Sparticle2.InitializeTree();

  CombinatoricGroup COMB_J_ISR_Sparticle2("COMB_J_ISR_Sparticle2", "Combinatoric System of Jets");
  MinMassesSqCombJigsaw CombSplitSq_J_ISR_Sparticle2_Transverse("CombSplitSq_ISR", "Minimize M_{T}_{ISR}^{2} + M_{T}_{S}^{2}",2,2);
  MinMassesSqCombJigsaw CombSplitSq_JS_ISR_Sparticle2_Transverse("CombSplitSq_S", "Minimize M_{T}_{Va}^{2} + M_{T}_{Vb}^{2}",2,2);
  
  CombinatoricGroup COMB_L_ISR_Sparticle2("COMB_L_ISR_Sparticle2", "Combinatoric System of Leps");
  MinMassesSqCombJigsaw CombSplitSq_L_ISR_Sparticle2_Transverse("CombSplitSq_L_ISR_Sparticle2", "Minimize M_{T}_{Va}^{2} + M_{T}_{Vb}^{2}",2,2);

  COMB_L_ISR_Sparticle2.AddFrame(L_X2a_ISR_Sparticle2);
  COMB_L_ISR_Sparticle2.SetNElementsForFrame(L_X2a_ISR_Sparticle2, 0);
  COMB_L_ISR_Sparticle2.AddFrame(L_X2b_ISR_Sparticle2);
  COMB_L_ISR_Sparticle2.SetNElementsForFrame(L_X2b_ISR_Sparticle2, 0);
          
  COMB_L_ISR_Sparticle2.AddJigsaw(CombSplitSq_L_ISR_Sparticle2_Transverse);
  CombSplitSq_L_ISR_Sparticle2_Transverse.SetTransverse();
  CombSplitSq_L_ISR_Sparticle2_Transverse.AddCombFrame(L_X2a_ISR_Sparticle2, 0);
  CombSplitSq_L_ISR_Sparticle2_Transverse.AddCombFrame(L_X2b_ISR_Sparticle2, 1);
  CombSplitSq_L_ISR_Sparticle2_Transverse.AddObjectFrame(X2a_ISR_Sparticle2, 0);
  CombSplitSq_L_ISR_Sparticle2_Transverse.AddObjectFrame(X2b_ISR_Sparticle2, 1);
        
  COMB_J_ISR_Sparticle2.AddFrame(ISR_ISR_Sparticle2);
  COMB_J_ISR_Sparticle2.SetNElementsForFrame(ISR_ISR_Sparticle2, 1);
  COMB_J_ISR_Sparticle2.AddFrame(J_X2a_ISR_Sparticle2);
  COMB_J_ISR_Sparticle2.SetNElementsForFrame(J_X2a_ISR_Sparticle2, 0);
      
  COMB_J_ISR_Sparticle2.AddJigsaw(CombSplitSq_J_ISR_Sparticle2_Transverse);
  CombSplitSq_J_ISR_Sparticle2_Transverse.SetTransverse();
  CombSplitSq_J_ISR_Sparticle2_Transverse.AddCombFrame(ISR_ISR_Sparticle2, 0);
  CombSplitSq_J_ISR_Sparticle2_Transverse.AddCombFrame(J_X2a_ISR_Sparticle2, 1);
  CombSplitSq_J_ISR_Sparticle2_Transverse.AddObjectFrame(ISR_ISR_Sparticle2, 0);
  CombSplitSq_J_ISR_Sparticle2_Transverse.AddObjectFrame(S_ISR_Sparticle2, 1);

  InvisibleGroup INV_ISR_Sparticle2("INV_ISR_Sparticle2","Invisible System");
  INV_ISR_Sparticle2.AddFrame(Ia_ISR_Sparticle2);
  INV_ISR_Sparticle2.AddFrame(Ib_ISR_Sparticle2);
  
  SetMassInvJigsaw InvM_ISR_Sparticle2("InvM_ISR_Sparticle2", "Set inv. system mass");
  INV_ISR_Sparticle2.AddJigsaw(InvM_ISR_Sparticle2);
  
  SetRapidityInvJigsaw InvEta_ISR_Sparticle2("InvEta_ISR_Sparticle2", "Set inv. system rapidity");
  INV_ISR_Sparticle2.AddJigsaw(InvEta_ISR_Sparticle2);
  InvEta_ISR_Sparticle2.AddVisibleFrames(S_ISR_Sparticle2.GetListVisibleFrames());
  
  MinMassesSqInvJigsaw InvSplit_ISR_Sparticle2("InvSplit_ISR_Sparticle2", "INV -> I_{a} + I_{b}", 2);
  INV_ISR_Sparticle2.AddJigsaw(InvSplit_ISR_Sparticle2);
  InvSplit_ISR_Sparticle2.AddVisibleFrame(J_X2a_ISR_Sparticle2, 0);
  InvSplit_ISR_Sparticle2.AddVisibleFrame(L_X2a_ISR_Sparticle2, 0);
  InvSplit_ISR_Sparticle2.AddVisibleFrame(L_X2b_ISR_Sparticle2, 1);
  InvSplit_ISR_Sparticle2.AddInvisibleFrame(Ia_ISR_Sparticle2, 0);
  InvSplit_ISR_Sparticle2.AddInvisibleFrame(Ib_ISR_Sparticle2, 1);

  LAB_ISR_Sparticle2.InitializeAnalysis();
  if(treeplot){
    tree_plot.SetTree(LAB_ISR_Sparticle2);
    tree_plot.Draw("ISR_Sparticle2_tree", "Analysis Tree", treeplot_INV);
    tree_plot.SetTree(INV_ISR_Sparticle2);
    tree_plot.Draw("ISR_Sparticle2_inv", "Invisible Jigsaws", treeplot_INV);
    tree_plot.SetTree(COMB_J_ISR_Sparticle2);
    tree_plot.Draw("ISR_Sparticle2_jets", "Jet Jigsaws", treeplot_INV);
    tree_plot.SetTree(COMB_L_ISR_Sparticle2);
    tree_plot.Draw("ISR_Sparticle2_leps", "Lep Jigsaws", treeplot_INV);
  }

  // ISR after BDT tree (default)

  LabRecoFrame LAB("LAB","lab");
  DecayRecoFrame CM("CM","CM");
  VisibleRecoFrame ISR("ISR","ISR");
  DecayRecoFrame S("S","S");
  DecayRecoFrame X3a("X3a","X3a");
  DecayRecoFrame X3b("X3b","X3b");
  DecayRecoFrame X2a("X2a","X2a");
  DecayRecoFrame X2b("X2b","X2b");
  SelfAssemblingRecoFrame saJa("saJa","saJa");
  SelfAssemblingRecoFrame saJb("saJb","saJb");
  SelfAssemblingRecoFrame saLa("saLa","saLa");
  SelfAssemblingRecoFrame saLb("saLb","saLb");
  VisibleRecoFrame Ja("Ja","Ja");
  VisibleRecoFrame Jb("Jb","Jb");
  VisibleRecoFrame La("La","La");
  VisibleRecoFrame Lb("Lb","Lb");
  InvisibleRecoFrame Ia("Ia","Ia");
  InvisibleRecoFrame Ib("Ib","Ib");

  LAB.SetChildFrame(CM);
  CM.AddChildFrame(ISR);
  CM.AddChildFrame(S);
  S.AddChildFrame(X3a);
  S.AddChildFrame(X3b);
  X3a.AddChildFrame(X2a);
  X3b.AddChildFrame(X2b);
  X3a.AddChildFrame(saJa);
  X2a.AddChildFrame(saLa);
  X2a.AddChildFrame(Ia);
  X3b.AddChildFrame(saJb);
  X2b.AddChildFrame(saLb);
  X2b.AddChildFrame(Ib);
  saJa.AddChildFrame(Ja);
  saJb.AddChildFrame(Jb);
  saLa.AddChildFrame(La);
  saLb.AddChildFrame(Lb);

  LAB.InitializeTree();

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
  InvSplit.AddVisibleFrame(Ja, 0);
  InvSplit.AddVisibleFrame(Jb, 1);
  InvSplit.AddVisibleFrame(La, 0);
  InvSplit.AddVisibleFrame(Lb, 1);
  InvSplit.AddInvisibleFrame(Ia, 0);
  InvSplit.AddInvisibleFrame(Ib, 1);
  
  CombinatoricGroup COMB_J("COMB_J", "Combinatoric System of Jets");
  MinMassesCombJigsaw CombSplit_ISR("CombSplit_ISR", "Minimize M_{T}^{ISR} and M_{T}^{S}");
  MinMassesSqCombJigsaw CombSplit_J("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
  
  CombinatoricGroup COMB_L("COMB_L", "Combinatoric System of Leps");
  MinMassesSqCombJigsaw CombSplit_L("CombSplit_L", "Minimize M_{Va}^{2} + M_{Vb}^{2}",2,2);
  
  COMB_J.AddFrame(ISR);
  COMB_J.SetNElementsForFrame(ISR, 1);
  COMB_J.AddFrame(Ja);
  COMB_J.SetNElementsForFrame(Ja, 0);
  COMB_J.AddFrame(Jb);
  COMB_J.SetNElementsForFrame(Jb, 0);
  
  COMB_J.AddJigsaw(CombSplit_ISR);
  CombSplit_ISR.SetTransverse();
  CombSplit_ISR.AddCombFrame(ISR, 0);
  CombSplit_ISR.AddCombFrame(Ja, 1);
  CombSplit_ISR.AddCombFrame(Jb, 1);
  CombSplit_ISR.AddObjectFrame(ISR, 0);
  CombSplit_ISR.AddObjectFrame(S, 1);
    
  COMB_J.AddJigsaw(CombSplit_J);
  CombSplit_J.AddCombFrame(Ja, 0);
  CombSplit_J.AddCombFrame(Jb, 1);
  CombSplit_J.AddObjectFrames(X3a.GetListVisibleFrames(), 0);
  CombSplit_J.AddObjectFrames(X3b.GetListVisibleFrames(), 1);
    
  COMB_L.AddFrame(La);
  COMB_L.SetNElementsForFrame(La, 0);
  COMB_L.AddFrame(Lb);
  COMB_L.SetNElementsForFrame(Lb, 0);
    
  COMB_L.AddJigsaw(CombSplit_L);
  CombSplit_L.AddCombFrame(La, 0);
  CombSplit_L.AddCombFrame(Lb, 1);
  CombSplit_L.AddObjectFrames(X3a.GetListVisibleFrames(), 0);
  CombSplit_L.AddObjectFrames(X3b.GetListVisibleFrames(), 1);

  LAB.InitializeAnalysis();

  // new tree that runs BDT at end

  LabRecoFrame LAB_ISR_BDT("LAB_ISR_BDT","lab_ISR_BDT");
  DecayRecoFrame CM_ISR_BDT("CM_ISR_BDT","CM_ISR_BDT");
  VisibleRecoFrame ISR_ISR_BDT("ISR_ISR_BDT","ISR_ISR_BDT");
  DecayRecoFrame S_ISR_BDT("S_ISR_BDT","S_ISR_BDT");
  DecayRecoFrame X2a_ISR_BDT("X2a_ISR_BDT","X2a_ISR_BDT");
  DecayRecoFrame X2b_ISR_BDT("X2b_ISR_BDT","X2b_ISR_BDT");
  SelfAssemblingRecoFrame saISR_ISR_BDT("saISR_ISR_BDT","saISR_ISR_BDT");
  SelfAssemblingRecoFrame saVa_ISR_BDT("saVa_ISR_BDT","saVa_ISR_BDT");
  SelfAssemblingRecoFrame saVb_ISR_BDT("saVb_ISR_BDT","saVb_ISR_BDT");
  VisibleRecoFrame Ja_ISR_BDT("Ja_ISR_BDT","Ja_ISR_BDT");
  VisibleRecoFrame Jb_ISR_BDT("Jb_ISR_BDT","Jb_ISR_BDT");
  VisibleRecoFrame La_ISR_BDT("La_ISR_BDT","La_ISR_BDT");
  VisibleRecoFrame Lb_ISR_BDT("Lb_ISR_BDT","Lb_ISR_BDT");
  InvisibleRecoFrame Ia_ISR_BDT("Ia_ISR_BDT","Ia_ISR_BDT");
  InvisibleRecoFrame Ib_ISR_BDT("Ib_ISR_BDT","Ib_ISR_BDT");

  LAB_ISR_BDT.SetChildFrame(CM_ISR_BDT);
  CM_ISR_BDT.AddChildFrame(saISR_ISR_BDT);
  saISR_ISR_BDT.AddChildFrame(ISR_ISR_BDT);
  CM_ISR_BDT.AddChildFrame(S_ISR_BDT);
  S_ISR_BDT.AddChildFrame(X2a_ISR_BDT);
  S_ISR_BDT.AddChildFrame(X2b_ISR_BDT);
  X2a_ISR_BDT.AddChildFrame(Ia_ISR_BDT);
  X2b_ISR_BDT.AddChildFrame(Ib_ISR_BDT);
  X2a_ISR_BDT.AddChildFrame(saVa_ISR_BDT);
  X2b_ISR_BDT.AddChildFrame(saVb_ISR_BDT);
  saVa_ISR_BDT.AddChildFrame(Ja_ISR_BDT);
  saVa_ISR_BDT.AddChildFrame(La_ISR_BDT);
  saVb_ISR_BDT.AddChildFrame(Jb_ISR_BDT);
  saVb_ISR_BDT.AddChildFrame(Lb_ISR_BDT);
  LAB_ISR_BDT.InitializeTree();

  InvisibleGroup INV_ISR_BDT("INV","Invisible System");
  INV_ISR_BDT.AddFrame(Ia_ISR_BDT);
  INV_ISR_BDT.AddFrame(Ib_ISR_BDT);
  
  SetMassInvJigsaw InvM_ISR_BDT("InvM", "Set inv. system mass");
  INV_ISR_BDT.AddJigsaw(InvM_ISR_BDT);
  
  SetRapidityInvJigsaw InvEta_ISR_BDT("InvEta", "Set inv. system rapidity");
  INV_ISR_BDT.AddJigsaw(InvEta_ISR_BDT);
  InvEta_ISR_BDT.AddVisibleFrames(S_ISR_BDT.GetListVisibleFrames());
  
  MinMassesSqInvJigsaw InvSplit_ISR_BDT("InvSplit", "INV -> I_{a} + I_{b}", 2);
  INV_ISR_BDT.AddJigsaw(InvSplit_ISR_BDT);
  InvSplit_ISR_BDT.AddVisibleFrame(Ja_ISR_BDT, 0);
  InvSplit_ISR_BDT.AddVisibleFrame(Jb_ISR_BDT, 1);
  InvSplit_ISR_BDT.AddVisibleFrame(La_ISR_BDT, 0);
  InvSplit_ISR_BDT.AddVisibleFrame(Lb_ISR_BDT, 1);
  InvSplit_ISR_BDT.AddInvisibleFrame(Ia_ISR_BDT, 0);
  InvSplit_ISR_BDT.AddInvisibleFrame(Ib_ISR_BDT, 1);
  
  CombinatoricGroup COMB_J_ISR_BDT("COMB_J", "Combinatoric System of Jets");
  MinMassesCombJigsaw CombSplit_ISR_ISR_BDT("CombSplit_ISR", "Minimize M_{T}^{ISR} and M_{T}^{S}");
  MinMassesSqCombJigsaw CombSplit_J_ISR_BDT("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
  
  CombinatoricGroup COMB_L_ISR_BDT("COMB_L", "Combinatoric System of Leps");
  MinMassesSqCombJigsaw CombSplit_L_ISR_BDT("CombSplit_L", "Minimize M_{Va}^{2} + M_{Vb}^{2}",2,2);
  
  COMB_J_ISR_BDT.AddFrame(ISR_ISR_BDT);
  COMB_J_ISR_BDT.SetNElementsForFrame(ISR_ISR_BDT, 1);
  COMB_J_ISR_BDT.AddFrame(Ja_ISR_BDT);
  COMB_J_ISR_BDT.SetNElementsForFrame(Ja_ISR_BDT, 0);
  COMB_J_ISR_BDT.AddFrame(Jb_ISR_BDT);
  COMB_J_ISR_BDT.SetNElementsForFrame(Jb_ISR_BDT, 0);
  
  COMB_J_ISR_BDT.AddJigsaw(CombSplit_ISR_ISR_BDT);
  CombSplit_ISR_ISR_BDT.SetTransverse();
  CombSplit_ISR_ISR_BDT.AddCombFrame(ISR_ISR_BDT, 0);
  CombSplit_ISR_ISR_BDT.AddCombFrame(Ja_ISR_BDT, 1);
  CombSplit_ISR_ISR_BDT.AddCombFrame(Jb_ISR_BDT, 1);
  CombSplit_ISR_ISR_BDT.AddObjectFrame(ISR_ISR_BDT, 0);
  CombSplit_ISR_ISR_BDT.AddObjectFrame(S_ISR_BDT, 1);
    
  COMB_J_ISR_BDT.AddJigsaw(CombSplit_J_ISR_BDT);
  CombSplit_J_ISR_BDT.AddCombFrame(Ja_ISR_BDT, 0);
  CombSplit_J_ISR_BDT.AddCombFrame(Jb_ISR_BDT, 1);
  CombSplit_J_ISR_BDT.AddObjectFrames(X2a_ISR_BDT.GetListVisibleFrames(), 0);
  CombSplit_J_ISR_BDT.AddObjectFrames(X2b_ISR_BDT.GetListVisibleFrames(), 1);

  COMB_L_ISR_BDT.AddFrame(La_ISR_BDT);
  COMB_L_ISR_BDT.SetNElementsForFrame(La_ISR_BDT, 0);
  COMB_L_ISR_BDT.AddFrame(Lb_ISR_BDT);
  COMB_L_ISR_BDT.SetNElementsForFrame(Lb_ISR_BDT, 0);
    
  COMB_L_ISR_BDT.AddJigsaw(CombSplit_L_ISR_BDT);
  CombSplit_L_ISR_BDT.AddCombFrame(La_ISR_BDT, 0);
  CombSplit_L_ISR_BDT.AddCombFrame(Lb_ISR_BDT, 1);
  CombSplit_L_ISR_BDT.AddObjectFrames(X2a_ISR_BDT.GetListVisibleFrames(), 0);
  CombSplit_L_ISR_BDT.AddObjectFrames(X2b_ISR_BDT.GetListVisibleFrames(), 1);

  LAB_ISR_BDT.InitializeAnalysis();

  // tree that takes BDT and assigns jets associated with leptons to ISR
  // then we take BDT and look for candidates
  // BDT_ISR

  LabRecoFrame LAB_BDT_ISR("LAB_BDT_ISR","lab_BDT_ISR");
  DecayRecoFrame CM_BDT_ISR("CM_BDT_ISR","CM_BDT_ISR");
  VisibleRecoFrame ISR_BDT_ISR("ISR_BDT_ISR","ISR_BDT_ISR");
  DecayRecoFrame S_BDT_ISR("S_BDT_ISR","S_BDT_ISR");
  DecayRecoFrame X2a_BDT_ISR("X2a_BDT_ISR","X2a_BDT_ISR");
  DecayRecoFrame X2b_BDT_ISR("X2b_BDT_ISR","X2b_BDT_ISR");
  SelfAssemblingRecoFrame saVa_BDT_ISR("saVa_BDT_ISR","saVa_BDT_ISR");
  SelfAssemblingRecoFrame saVb_BDT_ISR("saVb_BDT_ISR","saVb_BDT_ISR");
  VisibleRecoFrame Ja_BDT_ISR("Ja_BDT_ISR","Ja_BDT_ISR");
  VisibleRecoFrame Jb_BDT_ISR("Jb_BDT_ISR","Jb_BDT_ISR");
  VisibleRecoFrame La_BDT_ISR("La_BDT_ISR","La_BDT_ISR");
  VisibleRecoFrame Lb_BDT_ISR("Lb_BDT_ISR","Lb_BDT_ISR");
  InvisibleRecoFrame Ia_BDT_ISR("Ia_BDT_ISR","Ia_BDT_ISR");
  InvisibleRecoFrame Ib_BDT_ISR("Ib_BDT_ISR","Ib_BDT_ISR");

  LAB_BDT_ISR.SetChildFrame(CM_BDT_ISR);
  CM_BDT_ISR.AddChildFrame(ISR_BDT_ISR);
  CM_BDT_ISR.AddChildFrame(S_BDT_ISR);
  S_BDT_ISR.AddChildFrame(X2a_BDT_ISR);
  S_BDT_ISR.AddChildFrame(X2b_BDT_ISR);
  X2a_BDT_ISR.AddChildFrame(Ia_BDT_ISR);
  X2b_BDT_ISR.AddChildFrame(Ib_BDT_ISR);
  X2a_BDT_ISR.AddChildFrame(saVa_BDT_ISR);
  X2b_BDT_ISR.AddChildFrame(saVb_BDT_ISR);
  saVa_BDT_ISR.AddChildFrame(Ja_BDT_ISR);
  saVa_BDT_ISR.AddChildFrame(La_BDT_ISR);
  saVb_BDT_ISR.AddChildFrame(Jb_BDT_ISR);
  saVb_BDT_ISR.AddChildFrame(Lb_BDT_ISR);
  LAB_BDT_ISR.InitializeTree();

  InvisibleGroup INV_BDT_ISR("INV","Invisible System");
  INV_BDT_ISR.AddFrame(Ia_BDT_ISR);
  INV_BDT_ISR.AddFrame(Ib_BDT_ISR);
  
  SetMassInvJigsaw InvM_BDT_ISR("InvM", "Set inv. system mass");
  INV_BDT_ISR.AddJigsaw(InvM_BDT_ISR);
  
  SetRapidityInvJigsaw InvEta_BDT_ISR("InvEta", "Set inv. system rapidity");
  INV_BDT_ISR.AddJigsaw(InvEta_BDT_ISR);
  InvEta_BDT_ISR.AddVisibleFrames(S_BDT_ISR.GetListVisibleFrames());
  
  MinMassesSqInvJigsaw InvSplit_BDT_ISR("InvSplit", "INV -> I_{a} + I_{b}", 2);
  INV_BDT_ISR.AddJigsaw(InvSplit_BDT_ISR);
  InvSplit_BDT_ISR.AddVisibleFrame(Ja_BDT_ISR, 0);
  InvSplit_BDT_ISR.AddVisibleFrame(Jb_BDT_ISR, 1);
  InvSplit_BDT_ISR.AddVisibleFrame(La_BDT_ISR, 0);
  InvSplit_BDT_ISR.AddVisibleFrame(Lb_BDT_ISR, 1);
  InvSplit_BDT_ISR.AddInvisibleFrame(Ia_BDT_ISR, 0);
  InvSplit_BDT_ISR.AddInvisibleFrame(Ib_BDT_ISR, 1);
  
  CombinatoricGroup COMB_J_BDT_ISR("COMB_J", "Combinatoric System of Jets");
  MinMassesSqCombJigsaw CombSplit_J_BDT_ISR("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
  
  CombinatoricGroup COMB_L_BDT_ISR("COMB_L", "Combinatoric System of Leps");
  MinMassesSqCombJigsaw CombSplit_L_BDT_ISR("CombSplit_L", "Minimize M_{Va}^{2} + M_{Vb}^{2}",2,2);
  
  COMB_J_BDT_ISR.AddFrame(Ja_BDT_ISR);
  COMB_J_BDT_ISR.SetNElementsForFrame(Ja_BDT_ISR, 0);
  COMB_J_BDT_ISR.AddFrame(Jb_BDT_ISR);
  COMB_J_BDT_ISR.SetNElementsForFrame(Jb_BDT_ISR, 0);
    
  COMB_J_BDT_ISR.AddJigsaw(CombSplit_J_BDT_ISR);
  CombSplit_J_BDT_ISR.AddCombFrame(Ja_BDT_ISR, 0);
  CombSplit_J_BDT_ISR.AddCombFrame(Jb_BDT_ISR, 1);
  CombSplit_J_BDT_ISR.AddObjectFrames(X2a_BDT_ISR.GetListVisibleFrames(), 0);
  CombSplit_J_BDT_ISR.AddObjectFrames(X2b_BDT_ISR.GetListVisibleFrames(), 1);

  COMB_L_BDT_ISR.AddFrame(La_BDT_ISR);
  COMB_L_BDT_ISR.SetNElementsForFrame(La_BDT_ISR, 0);
  COMB_L_BDT_ISR.AddFrame(Lb_BDT_ISR);
  COMB_L_BDT_ISR.SetNElementsForFrame(Lb_BDT_ISR, 0);
    
  COMB_L_BDT_ISR.AddJigsaw(CombSplit_L_BDT_ISR);
  CombSplit_L_BDT_ISR.AddCombFrame(La_BDT_ISR, 0);
  CombSplit_L_BDT_ISR.AddCombFrame(Lb_BDT_ISR, 1);
  CombSplit_L_BDT_ISR.AddObjectFrames(X2a_BDT_ISR.GetListVisibleFrames(), 0);
  CombSplit_L_BDT_ISR.AddObjectFrames(X2b_BDT_ISR.GetListVisibleFrames(), 1);

  LAB_BDT_ISR.InitializeAnalysis();

  if(treeplot){
    tree_plot.WriteOutput("trees.root");
    std::cout << "Writing trees to trees.root" << endl;
    return 0;
  }

  int quark_cut = 0;
  int tot_quark = 0;

  int prong2 = 0;
  int prong3 = 0;
  int auntprong = 0;
  int lepprong = 0;
  int totcands = 0;
  int genHBosons = 0;
  int matchedAside = 0;
  int matchedBside = 0;
  int unmatchedAside = 0;
  int unmatchedBside = 0;

  int prong2_ISR_Sparticle = 0;
  int prong3_ISR_Sparticle = 0;
  int auntprong_ISR_Sparticle = 0;
  int lepprong_ISR_Sparticle = 0;
  int totcands_ISR_Sparticle = 0;

  int prong2_ISR_BDT = 0;
  int prong3_ISR_BDT = 0;
  int auntprong_ISR_BDT = 0;
  int lepprong_ISR_BDT = 0;
  int totcands_ISR_BDT = 0;
  int genHBosons_ISR_BDT = 0;
  int matchedAside_ISR_BDT = 0;
  int matchedBside_ISR_BDT = 0;
  int unmatchedAside_ISR_BDT = 0;
  int unmatchedBside_ISR_BDT = 0;

  int prong2_BDT_ISR_lep = 0;
  int prong3_BDT_ISR_lep = 0;
  int auntprong_BDT_ISR_lep = 0;
  int lepprong_BDT_ISR_lep = 0;
  int totcands_BDT_ISR_lep = 0;
  int genHBosons_BDT_ISR_lep = 0;
  int matchedAside_BDT_ISR_lep = 0;
  int matchedBside_BDT_ISR_lep = 0;
  int unmatchedAside_BDT_ISR_lep = 0;
  int unmatchedBside_BDT_ISR_lep = 0;

  int prong2_BDT_ISR_singlet = 0;
  int prong3_BDT_ISR_singlet = 0;
  int auntprong_BDT_ISR_singlet = 0;
  int lepprong_BDT_ISR_singlet = 0;
  int totcands_BDT_ISR_singlet = 0;
  int genHBosons_BDT_ISR_singlet = 0;
  int matchedAside_BDT_ISR_singlet = 0;
  int matchedBside_BDT_ISR_singlet = 0;
  int unmatchedAside_BDT_ISR_singlet = 0;
  int unmatchedBside_BDT_ISR_singlet = 0;

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
  //hists1.push_back(hist_NQuarks);

  TH1D* hist_V_had_cand = new TH1D((title+"_V_had_cand").c_str(), (title+"_V_had_cand;V_had_cand;").c_str(), 10, 0, 10);
  if(Cand) hists1.push_back(hist_V_had_cand);
  TH2D* hist_Njets_V_had_cand = new TH2D((title+"_Njets_V_had_cand").c_str(), (title+"_Njets_V_had_cand;Njets;Total cands").c_str(), 20, 0, 20, 10, 0, 10);
  if(Cand) hists2.push_back(hist_Njets_V_had_cand);
  TH2D* hist_Njets_matched_V_had_cand = new TH2D((title+"_Njets_matched_V_had_cand").c_str(), (title+"_Njets_matched_V_had_cand;Njets;Total Matched Cands").c_str(), 20, 0, 20, 5, 0, 5);
  if(Cand) hists2.push_back(hist_Njets_matched_V_had_cand);
  TH2D* hist_Njets_unmatched_V_had_cand = new TH2D((title+"_Njets_unmatched_V_had_cand").c_str(), (title+"_Njets_unmatched_V_had_cand;Njets;Total Unmatched Cands").c_str(), 20, 0, 20, 10, 0, 10);
  if(Cand) hists2.push_back(hist_Njets_unmatched_V_had_cand);
  TH2D* hist_Njets_partial_V_had_cand = new TH2D((title+"_Njets_partial_V_had_cand").c_str(), (title+"_Njets_partial_V_had_cand;Njets;Total Partial Cands").c_str(), 20, 0, 20, 5, 0, 5);
  if(Cand) hists2.push_back(hist_Njets_partial_V_had_cand);

  TH1D* hist_MatchedVJetsPt = new TH1D((title+"_MatchedVJetsPt").c_str(), (title+"_MatchedVJetsPt;MatchedVJetsPt").c_str(), g_NX, 0., 200.);
  if(Cand) hists1.push_back(hist_MatchedVJetsPt);
  TH1D* hist_UnmatchedVJetsPt = new TH1D((title+"_UnmatchedVJetsPt").c_str(), (title+"_UnmatchedVJetsPt;UnmatchedVJetsPt").c_str(), g_NX, 0., 200.);
  if(Cand) hists1.push_back(hist_UnmatchedVJetsPt);
  TH1D* hist_MatchedVJetsEta = new TH1D((title+"_MatchedVJetsEta").c_str(), (title+"_MatchedVJetsEta;MatchedVJetsEta").c_str(), g_NX, -2.5, 2.5);
  if(Cand) hists1.push_back(hist_MatchedVJetsEta);
  TH1D* hist_UnmatchedVJetsEta = new TH1D((title+"_UnmatchedVJetsEta").c_str(), (title+"_UnmatchedVJetsEta;UnmatchedVJetsEta").c_str(), g_NX, -2.5, 2.5);
  if(Cand) hists1.push_back(hist_UnmatchedVJetsEta);
  TH1D* hist_MatchedVJetsMass = new TH1D((title+"_MatchedVJetsMass").c_str(), (title+"_MatchedVJetsMass;MatchedVJetsMass").c_str(), g_NX, 0, 50.);
  if(Cand) hists1.push_back(hist_MatchedVJetsMass);
  TH1D* hist_UnmatchedVJetsMass = new TH1D((title+"_UnmatchedVJetsMass").c_str(), (title+"_UnmatchedVJetsMass;UnmatchedVJetsMass").c_str(), g_NX, 0., 50.);
  if(Cand) hists1.push_back(hist_UnmatchedVJetsMass);

  TH1D* hist_MatchedVPt = new TH1D((title+"_MatchedVPt").c_str(), (title+"_MatchedVPt;MatchedVPt").c_str(), g_NX, 0., 200.);
  if(Cand) hists1.push_back(hist_MatchedVPt);
  TH1D* hist_UnmatchedVPt = new TH1D((title+"_UnmatchedVPt").c_str(), (title+"_UnmatchedVPt;UnmatchedVPt").c_str(), g_NX, 0., 200.);
  if(Cand) hists1.push_back(hist_UnmatchedVPt);
  TH1D* hist_MatchedVEta = new TH1D((title+"_MatchedVEta").c_str(), (title+"_MatchedVEta;MatchedVEta").c_str(), g_NX, -2.5, 2.5);
  if(Cand) hists1.push_back(hist_MatchedVEta);
  TH1D* hist_UnmatchedVEta = new TH1D((title+"_UnmatchedVEta").c_str(), (title+"_UnmatchedVEta;UnmatchedVEta").c_str(), g_NX, -2.5, 2.5);
  if(Cand) hists1.push_back(hist_UnmatchedVEta);
  TH1D* hist_MatchedVMass = new TH1D((title+"_MatchedVMass").c_str(), (title+"_MatchedVMass;MatchedVMass").c_str(), g_NX, 0, 110.);
  if(Cand) hists1.push_back(hist_MatchedVMass);
  TH1D* hist_UnmatchedVMass = new TH1D((title+"_UnmatchedVMass").c_str(), (title+"_UnmatchedVMass;UnmatchedVMass").c_str(), g_NX, 0., 110.);
  if(Cand) hists1.push_back(hist_UnmatchedVMass);
  TH1D* hist_PartialmatchedVMass = new TH1D((title+"_PartialmatchedVMass").c_str(), (title+"_PartialmatchedVMass;PartialmatchedVMass").c_str(), g_NX, 0, 110.);
  if(Cand) hists1.push_back(hist_PartialmatchedVMass);
  TH1D* hist_MatchedVCosDecayAngle = new TH1D((title+"_MatchedVCosDecayAngle").c_str(), (title+"_MatchedVCosDecayAngle;MatchedVCosDecayAngle").c_str(), g_NX, -1., 1.);
  if(Cand) hists1.push_back(hist_MatchedVCosDecayAngle);
  TH1D* hist_UnmatchedVCosDecayAngle = new TH1D((title+"_UnmatchedVCosDecayAngle").c_str(), (title+"_UnmatchedVCosDecayAngle;UnmatchedVCosDecayAngle").c_str(), g_NX, -1., 1.);
  if(Cand) hists1.push_back(hist_UnmatchedVCosDecayAngle);

  TH1D* hist_MatchedVProngDeltaPhi = new TH1D((title+"_MatchedVProngDeltaPhi").c_str(), (title+"_MatchedVProngDeltaPhi;MatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
  if(Cand) hists1.push_back(hist_MatchedVProngDeltaPhi);
  TH1D* hist_PartialmatchedVProngDeltaPhi = new TH1D((title+"_PartialmatchedVProngDeltaPhi").c_str(), (title+"_PartialmatchedVProngDeltaPhi;PartialmatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
  if(Cand) hists1.push_back(hist_PartialmatchedVProngDeltaPhi);
  TH1D* hist_UnmatchedVProngDeltaPhi = new TH1D((title+"_UnmatchedVProngDeltaPhi").c_str(), (title+"_UnmatchedVProngDeltaPhi;UnmatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
  if(Cand) hists1.push_back(hist_UnmatchedVProngDeltaPhi);
  TH2D* hist_MatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_MatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_MatchedV_Mass_ProngDeltaPhi;matched Mass;matched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_ProngDeltaPhi);
  TH2D* hist_PartialmatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_PartialmatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_PartialmatchedV_Mass_ProngDeltaPhi;partial matched Mass;partial matched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_ProngDeltaPhi);
  TH2D* hist_UnmatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_UnmatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_UnmatchedV_Mass_ProngDeltaPhi;unmatched Mass;unmatched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_ProngDeltaPhi);
  TH1D* hist_MatchedVProngDeltaEta = new TH1D((title+"_MatchedVProngDeltaEta").c_str(), (title+"_MatchedVProngDeltaEta;MatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
  if(Cand) hists1.push_back(hist_MatchedVProngDeltaEta);
  TH1D* hist_PartialmatchedVProngDeltaEta = new TH1D((title+"_PartialmatchedVProngDeltaEta").c_str(), (title+"_PartialmatchedVProngDeltaEta;PartialmatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
  if(Cand) hists1.push_back(hist_PartialmatchedVProngDeltaEta);
  TH1D* hist_UnmatchedVProngDeltaEta = new TH1D((title+"_UnmatchedVProngDeltaEta").c_str(), (title+"_UnmatchedVProngDeltaEta;UnmatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
  if(Cand) hists1.push_back(hist_UnmatchedVProngDeltaEta);
  TH2D* hist_MatchedV_Mass_ProngDeltaEta = new TH2D((title+"_MatchedV_Mass_ProngDeltaEta").c_str(), (title+"_MatchedV_Mass_ProngDeltaEta;matched Mass;matched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_ProngDeltaEta);
  TH2D* hist_PartialmatchedV_Mass_ProngDeltaEta = new TH2D((title+"_PartialmatchedV_Mass_ProngDeltaEta").c_str(), (title+"_PartialmatchedV_Mass_ProngDeltaEta;partial matched Mass;partial matched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_ProngDeltaEta);
  TH2D* hist_UnmatchedV_Mass_ProngDeltaEta = new TH2D((title+"_UnmatchedV_Mass_ProngDeltaEta").c_str(), (title+"_UnmatchedV_Mass_ProngDeltaEta;unmatched Mass;unmatched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_ProngDeltaEta);
  TH1D* hist_MatchedVProngAbsDeltaEta = new TH1D((title+"_MatchedVProngAbsDeltaEta").c_str(), (title+"_MatchedVProngAbsDeltaEta;MatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
  if(Cand) hists1.push_back(hist_MatchedVProngAbsDeltaEta);
  TH1D* hist_PartialmatchedVProngAbsDeltaEta = new TH1D((title+"_PartialmatchedVProngAbsDeltaEta").c_str(), (title+"_PartialmatchedVProngAbsDeltaEta;PartialmatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
  if(Cand) hists1.push_back(hist_PartialmatchedVProngAbsDeltaEta);
  TH1D* hist_UnmatchedVProngAbsDeltaEta = new TH1D((title+"_UnmatchedVProngAbsDeltaEta").c_str(), (title+"_UnmatchedVProngAbsDeltaEta;UnmatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
  if(Cand) hists1.push_back(hist_UnmatchedVProngAbsDeltaEta);
  TH2D* hist_MatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_MatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_MatchedV_Mass_ProngAbsDeltaEta;matched Mass;matched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_ProngAbsDeltaEta);
  TH2D* hist_PartialmatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_PartialmatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_PartialmatchedV_Mass_ProngAbsDeltaEta;partial matched Mass;partial matched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_ProngAbsDeltaEta);
  TH2D* hist_UnmatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_UnmatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_UnmatchedV_Mass_ProngAbsDeltaEta;unmatched Mass;unmatched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_ProngAbsDeltaEta);

  TH2D* hist_MatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngDeltaPhi;matched CosDecayAngle_CM;matched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngDeltaPhi);
  TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi;partial matched CosDecayAngle_CM;partial matched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi);
  TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi;unmatched CosDecayAngle_CM;unmatched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
  if(Cand) hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi);
  TH2D* hist_MatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngDeltaEta;matched CosDecayAngle_CM;matched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngDeltaEta);
  TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta;partial matched CosDecayAngle_CM;partial matched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta);
  TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta;unmatched CosDecayAngle_CM;unmatched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
  if(Cand) hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta);
  TH2D* hist_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;matched CosDecayAngle_CM;matched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);
  TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;partial matched CosDecayAngle_CM;partial matched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);
  TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;unmatched CosDecayAngle_CM;unmatched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
  if(Cand) hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);

  TH2D* hist_MatchedV_Mass_CosDecayAngle = new TH2D((title+"_MatchedV_Mass_CosDecayAngle").c_str(), (title+"_MatchedV_Mass_CosDecayAngle;matched V mass; matched V CosDecayAngle;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_CosDecayAngle);
  TH2D* hist_UnmatchedV_Mass_CosDecayAngle = new TH2D((title+"_UnmatchedV_Mass_CosDecayAngle").c_str(), (title+"_UnmatchedV_Mass_CosDecayAngle;unmatched V mass; unmatched V CosDecayAngle;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_CosDecayAngle);

  TH2D* hist_MatchedV_PCM_PCand = new TH2D((title+"_MatchedV_PCM_PCand").c_str(), (title+"_MatchedV_PCM_PCand;PCM;matched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
  if(Cand) hists2.push_back(hist_MatchedV_PCM_PCand);
  TH2D* hist_UnmatchedV_PCM_PCand = new TH2D((title+"_UnmatchedV_PCM_PCand").c_str(), (title+"_UnmatchedV_PCM_PCand;PCM;unmatched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
  if(Cand) hists2.push_back(hist_UnmatchedV_PCM_PCand);
  TH2D* hist_PartialmatchedV_PCM_PCand = new TH2D((title+"_PartialmatchedV_PCM_PCand").c_str(), (title+"_PartialmatchedV_PCM_PCand;PCM;Partialmatched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_PCM_PCand);

  TH1D* hist_UnmatchedV_CosDecayAngle_CM = new TH1D((title+"_UnmatchedV_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM;unmatchedV_CosDecayAngle_CM").c_str(), g_NX/2., -1., 1.);
  if(Cand) hists1.push_back(hist_UnmatchedV_CosDecayAngle_CM);
  TH1D* hist_UnmatchedV_DeltaPhiDecayAngle_CM = new TH1D((title+"_UnmatchedV_DeltaPhiDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayAngle_CM;unmatchedV_DeltaPhiDecayAngle_CM").c_str(), g_NX/2., 0., 3.5);
  if(Cand) hists1.push_back(hist_UnmatchedV_DeltaPhiDecayAngle_CM);
  TH1D* hist_UnmatchedV_DeltaPhiDecayPlanes_CM = new TH1D((title+"_UnmatchedV_DeltaPhiDecayPlanes_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayPlanes_CM;unmatchedV_DeltaPhiDecayPlanes_CM").c_str(), g_NX/2., 0., 6.5);
  if(Cand) hists1.push_back(hist_UnmatchedV_DeltaPhiDecayPlanes_CM);
  TH2D* hist_UnmatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_Mass_CosDecayAngle_CM;unmatched Mass;unmatched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_UnmatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_UnmatchedV_Mass_DeltaPhiDecayAngle_CM;unmatched Mass;unmatched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_DeltaPhiDecayAngle_CM);
  TH2D* hist_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM;unmatched Mass;unmatched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
  if(Cand) hists2.push_back(hist_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM);

  TH2D* hist_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM = new TH2D((title+"_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM;DeltaPhiDecayAngle_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 3.5, g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM;DeltaPhiDecayAngle_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 3.5, g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM);
  TH2D* hist_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM = new TH2D((title+"_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM;DeltaPhiDecayPlanes_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 6.5, g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM;DeltaPhiDecayPlanes_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 6.5, g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM);

  TH1D* hist_MatchedV_CosDecayAngle_CM = new TH1D((title+"_MatchedV_CosDecayAngle_CM").c_str(), (title+"_MatchedV_CosDecayAngle_CM;matchedV_CosDecayAngle_CM").c_str(), g_NX/2., -1., 1.);
  if(Cand) hists1.push_back(hist_MatchedV_CosDecayAngle_CM);
  TH1D* hist_MatchedV_DeltaPhiDecayAngle_CM = new TH1D((title+"_MatchedV_DeltaPhiDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayAngle_CM;matchedV_DeltaPhiDecayAngle_CM").c_str(), g_NX/2., 0., 3.5);
  if(Cand) hists1.push_back(hist_MatchedV_DeltaPhiDecayAngle_CM);
  TH1D* hist_MatchedV_DeltaPhiDecayPlanes_CM = new TH1D((title+"_MatchedV_DeltaPhiDecayPlanes_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayPlanes_CM;matchedV_DeltaPhiDecayPlanes_CM").c_str(), g_NX/2., 0., 6.5);
  if(Cand) hists1.push_back(hist_MatchedV_DeltaPhiDecayPlanes_CM);
  TH2D* hist_MatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_MatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_MatchedV_Mass_CosDecayAngle_CM;matched Mass;matched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_CosDecayAngle_CM);
  TH2D* hist_MatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_MatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_MatchedV_Mass_DeltaPhiDecayAngle_CM;matched Mass;matched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_DeltaPhiDecayAngle_CM);
  TH2D* hist_MatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_MatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_MatchedV_Mass_DeltaPhiDecayPlanes_CM;matched Mass;matched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
  if(Cand) hists2.push_back(hist_MatchedV_Mass_DeltaPhiDecayPlanes_CM);
  
  TH2D* hist_PartialmatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_Mass_CosDecayAngle_CM;partial matched Mass;partial matched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_CosDecayAngle_CM);
  TH2D* hist_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM;partial matched Mass;partial matched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM);
  TH2D* hist_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM;partial matched Mass;partial matched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
  if(Cand) hists2.push_back(hist_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM);

  TH2D* hist_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM;partial matched CosDecayAngle;partial matched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM);
  TH2D* hist_MatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_MatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_MatchedV_CosDecayAngle_CosDecayAngle_CM;matched CosDecayAngle;matched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_CosDecayAngle_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_CosDecayAngle_CosDecayAngle_CM;unmatched CosDecayAngle;unmatched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_CosDecayAngle_CosDecayAngle_CM);

  TH1D* hist_Matched_CosDecayAngle_Hemi = new TH1D((title+"_Matched_CosDecayAngle_Hemi").c_str(), (title+"_Matched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
  if(Cand) hists1.push_back(hist_Matched_CosDecayAngle_Hemi);
  TH1D* hist_Partialmatched_CosDecayAngle_Hemi = new TH1D((title+"_Partialmatched_CosDecayAngle_Hemi").c_str(), (title+"_Partialmatched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
  if(Cand) hists1.push_back(hist_Partialmatched_CosDecayAngle_Hemi);
  TH1D* hist_Unmatched_CosDecayAngle_Hemi = new TH1D((title+"_Unmatched_CosDecayAngle_Hemi").c_str(), (title+"_Unmatched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
  if(Cand) hists1.push_back(hist_Unmatched_CosDecayAngle_Hemi);
  TH2D* hist_Matched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Matched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Matched_Mass_CosDecayAngle_Hemi;matched Mass;matched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_Matched_Mass_CosDecayAngle_Hemi);
  TH2D* hist_Unmatched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Unmatched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Unmatched_Mass_CosDecayAngle_Hemi;unmatched Mass;unmatched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_Unmatched_Mass_CosDecayAngle_Hemi);
  TH2D* hist_Partialmatched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Partialmatched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Partialmatched_Mass_CosDecayAngle_Hemi;partial matched Mass;partial matched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_Partialmatched_Mass_CosDecayAngle_Hemi);

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

  TH1D* hist_Count = new TH1D((title+"_Count").c_str(), (title+"_Count;").c_str(), 20, 0, 20); // generic counting
  TH1D* hist_CandCount = new TH1D((title+"_CandCount").c_str(), (title+"_CandCount;").c_str(), 5, 0, 5);
  TH1D* hist_CandCount_Combo = new TH1D((title+"_CandCount_Combo").c_str(), (title+"_CandCount_Combo;").c_str(), 5, 0, 5);
  TH1D* hist_CandCount_Chi2 = new TH1D((title+"_CandCount_Chi2").c_str(), (title+"_CandCount_Chi2;").c_str(), 5, 0, 5);
  TH1D* hist_Count_ISR_BDT = new TH1D((title+"_Count_ISR_BDT").c_str(), (title+"_Count_ISR_BDT;").c_str(), 20, 0, 20); // generic counting ISR_BDT
  TH1D* hist_CandCount_ISR_BDT = new TH1D((title+"_CandCount_ISR_BDT").c_str(), (title+"_CandCount_ISR_BDT;").c_str(), 5, 0, 5);
  TH1D* hist_Count_BDT_ISR_lep = new TH1D((title+"_Count_BDT_ISR_lep").c_str(), (title+"_Count_BDT_ISR_lep;").c_str(), 20, 0, 20); // generic counting BDT_ISR_lep
  TH1D* hist_CandCount_BDT_ISR_lep = new TH1D((title+"_CandCount_BDT_ISR_lep").c_str(), (title+"_CandCount_BDT_ISR_lep;").c_str(), 5, 0, 5);
  TH1D* hist_Count_BDT_ISR_singlet = new TH1D((title+"_Count_BDT_ISR_singlet").c_str(), (title+"_Count_BDT_ISR_singlet;").c_str(), 20, 0, 20); // generic counting BDT_ISR_singlet
  TH1D* hist_CandCount_BDT_ISR_singlet = new TH1D((title+"_CandCount_BDT_ISR_singlet").c_str(), (title+"_CandCount_BDT_ISR_singlet;").c_str(), 5, 0, 5);

  TH1D* hist_Count_ISR_Sparticle = new TH1D((title+"_Count_ISR_Sparticle").c_str(), (title+"_Count_ISR_Sparticle;").c_str(), 20, 0, 20); // generic counting

  TH1D* hist_N_V_had_ISR_Sparticle = new TH1D((title+"_N_V_had_ISR_Sparticle").c_str(), (title+"_N_V_had_ISR_Sparticle;Total Cands;").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_N_V_had_ISR_Sparticle);
  TH1D* hist_N_V_had_S_ISR_Sparticle = new TH1D((title+"_N_V_had_S_ISR_Sparticle").c_str(), (title+"_N_V_had_S_ISR_Sparticle;Total Cands In S System").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_N_V_had_S_ISR_Sparticle);
  TH1D* hist_N_V_had_ISR_ISR_Sparticle = new TH1D((title+"_N_V_had_ISR_ISR_Sparticle").c_str(), (title+"_N_V_had_ISR_ISR_Sparticle;Total Cands In ISR;").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_N_V_had_ISR_ISR_Sparticle);

  TH1D* hist_matched_S_ISR_Sparticle = new TH1D((title+"_matched_S_ISR_Sparticle").c_str(), (title+"_matched_S_ISR_Sparticle;matched cands S").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_matched_S_ISR_Sparticle);
  TH1D* hist_unmatched_S_ISR_Sparticle = new TH1D((title+"_unmatched_S_ISR_Sparticle").c_str(), (title+"_unmatched_S_ISR_Sparticle;unmatched cands S").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_unmatched_S_ISR_Sparticle);
  TH1D* hist_matched_ISR_ISR_Sparticle = new TH1D((title+"_matched_ISR_ISR_Sparticle").c_str(), (title+"_matched_ISR_ISR_Sparticle;matched cands ISR").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_matched_ISR_ISR_Sparticle);
  TH1D* hist_unmatched_ISR_ISR_Sparticle = new TH1D((title+"_unmatched_ISR_ISR_Sparticle").c_str(), (title+"_unmatched_ISR_ISR_Sparticle;unmatched cands ISR").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_unmatched_ISR_ISR_Sparticle);
  TH1D* hist_partial_S_ISR_Sparticle = new TH1D((title+"_partial_S_ISR_Sparticle").c_str(), (title+"_partial_S_ISR_Sparticle;partial cands S").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_partial_S_ISR_Sparticle);
  TH1D* hist_partial_ISR_ISR_Sparticle = new TH1D((title+"_partial_ISR_ISR_Sparticle").c_str(), (title+"_partial_ISR_ISR_Sparticle;partial cands ISR").c_str(), 5, 0, 5);
  if(Sparticle1) hists1.push_back(hist_partial_ISR_ISR_Sparticle);

  TH2D* hist_Njets_N_V_had_ISR_Sparticle = new TH2D((title+"_Njets_N_V_had_ISR_Sparticle").c_str(), (title+"_Njets_N_V_had_ISR_Sparticle;Njets;Total Cands").c_str(), 15, 0, 15, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_N_V_had_ISR_Sparticle);
  TH2D* hist_Njets_N_V_had_S_ISR_Sparticle = new TH2D((title+"_Njets_N_V_had_S_ISR_Sparticle").c_str(), (title+"_Njets_N_V_had_S_ISR_Sparticle;Njets;Total Cands In S System").c_str(), 15, 0, 15, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_N_V_had_S_ISR_Sparticle);
  TH2D* hist_Njets_N_V_had_ISR_ISR_Sparticle = new TH2D((title+"_Njets_N_V_had_ISR_ISR_Sparticle").c_str(), (title+"_Njets_N_V_had_ISR_ISR_Sparticle;Njets;Total Cands In ISR").c_str(), 15, 0, 15, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_N_V_had_ISR_ISR_Sparticle);

  TH2D* hist_Njets_S_matched_S_ISR_Sparticle = new TH2D((title+"_Njets_S_matched_S_ISR_Sparticle").c_str(), (title+"_Njets_S_matched_S_ISR_Sparticle;Njets S;matched cands S").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_S_matched_S_ISR_Sparticle);
  TH2D* hist_Njets_S_unmatched_S_ISR_Sparticle = new TH2D((title+"_Njets_S_unmatched_S_ISR_Sparticle").c_str(), (title+"_Njets_S_unmatched_S_ISR_Sparticle;Njets S;unmatched cands S").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_S_unmatched_S_ISR_Sparticle);
  TH2D* hist_Njets_ISR_matched_ISR_ISR_Sparticle = new TH2D((title+"_Njets_ISR_matched_ISR_ISR_Sparticle").c_str(), (title+"_Njets_ISR_matched_ISR_ISR_Sparticle;Njets ISR;matched cands ISR").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_ISR_matched_ISR_ISR_Sparticle);
  TH2D* hist_Njets_ISR_unmatched_ISR_ISR_Sparticle = new TH2D((title+"_Njets_ISR_unmatched_ISR_ISR_Sparticle").c_str(), (title+"_Njets_ISR_unmatched_ISR_ISR_Sparticle;Njets ISR;unmatched cands ISR").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_ISR_unmatched_ISR_ISR_Sparticle);
  TH2D* hist_Njets_S_partial_S_ISR_Sparticle = new TH2D((title+"_Njets_S_partial_S_ISR_Sparticle").c_str(), (title+"_Njets_S_partial_S_ISR_Sparticle;Njets S;partial cands S").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_S_partial_S_ISR_Sparticle);
  TH2D* hist_Njets_ISR_partial_ISR_ISR_Sparticle = new TH2D((title+"_Njets_ISR_partial_ISR_ISR_Sparticle").c_str(), (title+"_Njets_ISR_partial_ISR_ISR_Sparticle;Njets ISR;partial cands ISR").c_str(), 10, 0, 10, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_Njets_ISR_partial_ISR_ISR_Sparticle);

  TH2D* hist_matched_Njets_S_matched_S_ISR_Sparticle = new TH2D((title+"_matched_Njets_S_matched_S_ISR_Sparticle").c_str(), (title+"_matched_Njets_S_matched_S_ISR_Sparticle;matched Njets S;matched cands S").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_Njets_S_matched_S_ISR_Sparticle);
  TH2D* hist_unmatched_Njets_S_unmatched_S_ISR_Sparticle = new TH2D((title+"_unmatched_Njets_S_unmatched_S_ISR_Sparticle").c_str(), (title+"_unmatched_Njets_S_unmatched_S_ISR_Sparticle;unmatched Njets S;unmatched cands S").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_unmatched_Njets_S_unmatched_S_ISR_Sparticle);
  TH2D* hist_matched_Njets_ISR_matched_ISR_ISR_Sparticle = new TH2D((title+"_matched_Njets_ISR_matched_ISR_ISR_Sparticle").c_str(), (title+"_matched_Njets_ISR_matched_ISR_ISR_Sparticle;matched Njets ISR;matched cands ISR").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_Njets_ISR_matched_ISR_ISR_Sparticle);
  TH2D* hist_unmatched_Njets_ISR_unmatched_ISR_ISR_Sparticle = new TH2D((title+"_unmatched_Njets_ISR_unmatched_ISR_ISR_Sparticle").c_str(), (title+"_unmatched_Njets_ISR_unmatched_ISR_ISR_Sparticle;unmatched Njets ISR;unmatched cands ISR").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_unmatched_Njets_ISR_unmatched_ISR_ISR_Sparticle);
  TH2D* hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle = new TH2D((title+"_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle").c_str(), (title+"_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle;matched Njets S;unmatched Njets S").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle);
  TH2D* hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle = new TH2D((title+"_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle").c_str(), (title+"_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle;matched Njets ISR;unmatched Njets ISR").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle);
  TH2D* hist_matched_S_matched_Sa_ISR_Sparticle = new TH2D((title+"_matched_S_matched_Sa_ISR_Sparticle").c_str(), (title+"_matched_S_matched_Sa_ISR_Sparticle;matched cands S;matched cands Sa").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_S_matched_Sa_ISR_Sparticle);
  TH2D* hist_matched_S_matched_Sb_ISR_Sparticle = new TH2D((title+"_matched_S_matched_Sb_ISR_Sparticle").c_str(), (title+"_matched_S_matched_Sb_ISR_Sparticle;matched cands S;matched cands Sb").c_str(), 5, 0, 5, 5, 0, 5);
  if(Sparticle1) hists2.push_back(hist_matched_S_matched_Sb_ISR_Sparticle);
  
  TH2D* hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle2 = new TH2D((title+"_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle2").c_str(), (title+"_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle2;matched Njets S;unmatched Njets S").c_str(), 5, 0, 5, 5, 0, 5);
  //if(Sparticle2) hists2.push_back(hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle2);
  TH2D* hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle2 = new TH2D((title+"_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle2").c_str(), (title+"_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle2;matched Njets ISR;unmatched Njets ISR").c_str(), 5, 0, 5, 5, 0, 5);
  //if(Sparticle2) hists2.push_back(hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle2);
  TH2D* hist_matched_Njets_Sa_matched_Njets_Sb_ISR_Sparticle2 = new TH2D((title+"_matched_Njets_Sa_matched_Njets_Sb_ISR_Sparticle2").c_str(), (title+"_matched_Njets_Sa_matched_Njets_Sb_ISR_Sparticle2;matched Njets Sa;matched Njets Sb").c_str(), 5, 0, 5, 5, 0, 5);
  //if(Sparticle2) hists2.push_back(hist_matched_Njets_Sa_matched_Njets_Sb_ISR_Sparticle2);
  TH2D* hist_unmatched_Njets_Sa_unmatched_Njets_Sb_ISR_Sparticle2 = new TH2D((title+"_unmatched_Njets_Sa_unmatched_Njets_Sb_ISR_Sparticle2").c_str(), (title+"_unmatched_Njets_Sa_unmatched_Njets_Sb_ISR_Sparticle2;unmatched Njets Sa;unmatched Njets Sb").c_str(), 5, 0, 5, 5, 0, 5);
  //if(Sparticle2) hists2.push_back(hist_unmatched_Njets_Sa_unmatched_Njets_Sb_ISR_Sparticle2);

  TH2D* hist_Matched_DeltaR_P = new TH2D((title+"_Matched_DeltaR_P").c_str(), (title+"_Matched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
  if(Cand) hists2.push_back(hist_Matched_DeltaR_P);
  TH2D* hist_Unmatched_DeltaR_P = new TH2D((title+"_Unmatched_DeltaR_P").c_str(), (title+"_Unmatched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
  if(Cand) hists2.push_back(hist_Unmatched_DeltaR_P);

  TH2D* hist_combo_Matched_DeltaR_P = new TH2D((title+"_combo_Matched_DeltaR_P").c_str(), (title+"_combo_Matched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
  if(Cand) hists2.push_back(hist_combo_Matched_DeltaR_P);
  TH2D* hist_combo_Unmatched_DeltaR_P = new TH2D((title+"_combo_Unmatched_DeltaR_P").c_str(), (title+"_combo_Unmatched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
  if(Cand) hists2.push_back(hist_combo_Unmatched_DeltaR_P);

  TH2D* hist_Matched_lnDeltaR_lnP = new TH2D((title+"_Matched_lnDeltaR_lnP").c_str(), (title+"_Matched_lnDeltaR_lnP;lnDeltaR;lnP;").c_str(), g_NX/2., 0., 1.5, g_NX/2., 0., 7.);
  if(Cand) hists2.push_back(hist_Matched_lnDeltaR_lnP);
  TH2D* hist_Unmatched_lnDeltaR_lnP = new TH2D((title+"_Unmatched_lnDeltaR_lnP").c_str(), (title+"_Unmatched_lnDeltaR_lnP;lnDeltaR;lnP;").c_str(), g_NX/2., 0., 1.5, g_NX/2., 0., 7.);
  if(Cand) hists2.push_back(hist_Unmatched_lnDeltaR_lnP);

  TH2D* hist_Matched_Mass_P_Multi_DeltaR = new TH2D((title+"_Matched_Mass_P_Multi_DeltaR").c_str(), (title+"_Matched_Mass_P_Multi_DeltaR;Mass;P*DeltaR;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 1000.);
  if(Cand) hists2.push_back(hist_Matched_Mass_P_Multi_DeltaR);
  TH2D* hist_Unmatched_Mass_P_Multi_DeltaR = new TH2D((title+"_Unmatched_Mass_P_Multi_DeltaR").c_str(), (title+"_Unmatched_Mass_P_Multi_DeltaR;Mass;P*DeltaR;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 1000.);
  if(Cand) hists2.push_back(hist_Unmatched_Mass_P_Multi_DeltaR);

  TH2D* hist_Matched_CosDecayAngle_CM_P_Multi_DeltaR = new TH2D((title+"_Matched_CosDecayAngle_CM_P_Multi_DeltaR").c_str(), (title+"_Matched_CosDecayAngle_CM_P_Multi_DeltaR;CosDecayAngle_CM;P*DeltaR;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 1000.);
  if(Cand) hists2.push_back(hist_Matched_CosDecayAngle_CM_P_Multi_DeltaR);
  TH2D* hist_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR = new TH2D((title+"_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR").c_str(), (title+"_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR;CosDecayAngle_CM;P*DeltaR;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 1000.);
  if(Cand) hists2.push_back(hist_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR);

  TH1D* hist_BaconNum = new TH1D((title+"_BaconNum").c_str(), (title+"_BaconNum;BaconNum").c_str(), 10, 0, 10);
  if(Cand) hists1.push_back(hist_BaconNum);

  TH1D* hist_MatchedlepCandBaconNum = new TH1D((title+"_MatchedlepCandBaconNum").c_str(), (title+"_MatchedlepCandBaconNum;MatchedlepCandBaconNum").c_str(), 10, 0, 10);
  if(Cand) hists1.push_back(hist_MatchedlepCandBaconNum);
  TH1D* hist_PartialMatchedlepCandBaconNum = new TH1D((title+"_PartialMatchedlepCandBaconNum").c_str(), (title+"_PartialMatchedlepCandBaconNum;PartialMatchedlepCandBaconNum").c_str(), 10, 0, 10);
  if(Cand) hists1.push_back(hist_PartialMatchedlepCandBaconNum);
  TH1D* hist_UnmatchedlepCandBaconNum = new TH1D((title+"_UnmatchedlepCandBaconNum").c_str(), (title+"_UnmatchedlepCandBaconNum;UnmatchedlepCandBaconNum").c_str(), 10, 0, 10);
  if(Cand) hists1.push_back(hist_UnmatchedlepCandBaconNum);

  TH2D* hist_MatchedlepCandBaconNum_Njets = new TH2D((title+"_MatchedlepCandBaconNum_Njets").c_str(), (title+"_MatchedlepCandBaconNum;MatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
  if(Cand) hists2.push_back(hist_MatchedlepCandBaconNum_Njets);
  TH2D* hist_UnmatchedlepCandBaconNum_Njets = new TH2D((title+"_UnmatchedlepCandBaconNum_Njets").c_str(), (title+"_UnmatchedlepCandBaconNum;UnmatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
  if(Cand) hists2.push_back(hist_UnmatchedlepCandBaconNum_Njets);
  TH2D* hist_PartialMatchedlepCandBaconNum_Njets = new TH2D((title+"_PartialMatchedlepCandBaconNum_Njets").c_str(), (title+"_PartialMatchedlepCandBaconNum;PartialMatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
  if(Cand) hists2.push_back(hist_PartialMatchedlepCandBaconNum_Njets);

  TH1D* hist_chi2WMass = new TH1D((title+"_chi2WMass").c_str(), (title+"_chi2WMass;chi2WMass").c_str(), g_NX, 30., 120.);
  if(Cand) hists1.push_back(hist_chi2WMass);

  TH2D* hist_MatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_MatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_MatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
  if(Cand) hists2.push_back(hist_MatchedCand_JetA_Mass_JetB_Mass);
  TH2D* hist_UnmatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_UnmatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_UnmatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
  if(Cand) hists2.push_back(hist_UnmatchedCand_JetA_Mass_JetB_Mass);
  TH2D* hist_PartialMatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_PartialMatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_PartialMatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
  if(Cand) hists2.push_back(hist_PartialMatchedCand_JetA_Mass_JetB_Mass);

  // Prong Mass Ratio
  TH2D* hist_MatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_MatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_MatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedCand_PMR_CosDecayAngle_CM);
  TH2D* hist_PartialMatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_PartialMatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_PartialMatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialMatchedCand_PMR_CosDecayAngle_CM);
  TH2D* hist_UnmatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedCand_PMR_CosDecayAngle_CM);

  TH2D* hist_MatchedCand_PMR_Mass = new TH2D((title+"_MatchedCand_PMR_Mass").c_str(), (title+"_MatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
  if(Cand) hists2.push_back(hist_MatchedCand_PMR_Mass);
  TH2D* hist_PartialMatchedCand_PMR_Mass = new TH2D((title+"_PartialMatchedCand_PMR_Mass").c_str(), (title+"_PartialMatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
  if(Cand) hists2.push_back(hist_PartialMatchedCand_PMR_Mass);
  TH2D* hist_UnmatchedCand_PMR_Mass = new TH2D((title+"_UnmatchedCand_PMR_Mass").c_str(), (title+"_UnmatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
  if(Cand) hists2.push_back(hist_UnmatchedCand_PMR_Mass);


  TEfficiency* eff_quarkMatchedGenJet_eta = new TEfficiency((title+"_eff_quarkMatchedGenJet_eta").c_str(),"Efficiency of Quark Getting Matched to GenJet;quarkEta;Quark Matching Efficiency", g_NX, -5., 5.);
  if(Cand) effs.push_back(eff_quarkMatchedGenJet_eta);
  TEfficiency* eff_GenJetMatchedAll_V_cand_eta = new TEfficiency((title+"_eff_GenJetMatchedAll_V_cand_eta").c_str(),"Efficiency of V cand having all Matched GenJet;V_candEta;Efficiency", g_NX, -5., 5.);
  if(Cand) effs.push_back(eff_GenJetMatchedAll_V_cand_eta);
  TEfficiency* eff_GenJetMatchedAll_V_cand_mass = new TEfficiency((title+"_eff_GenJetMatchedAll_V_cand_mass").c_str(),"Efficiency of V cand having all Matched GenJet;V_candMass;Efficiency", g_NX, 0., 110.);
  if(Cand) effs.push_back(eff_GenJetMatchedAll_V_cand_mass);

  TEfficiency* eff_bosonMatched_MatchedV_genboson_Mass = new TEfficiency((title+"_eff_bosonMatched_MatchedV_genboson_Mass").c_str(),"Efficiency of V Cand Matching to Gen Boson;genBosonMass;Efficiency", g_NX, 50., 110.);
  if(Cand) effs.push_back(eff_bosonMatched_MatchedV_genboson_Mass);
  TEfficiency* eff_bosonMatched_MatchedV_genboson_Eta = new TEfficiency((title+"_eff_bosonMatched_MatchedV_genboson_Eta").c_str(),"Efficiency of V Cand Matching to Gen Boson;genBosonEta;Efficiency", g_NX, -5., 5.);
  if(Cand) effs.push_back(eff_bosonMatched_MatchedV_genboson_Eta);
  TEfficiency* eff_bosonMatched_UnmatchedV_genboson_Mass = new TEfficiency((title+"_eff_bosonMatched_UnmatchedV_genboson_Mass").c_str(),"Efficiency of V Cand Not Matching to Gen Boson;genBosonMass;Efficiency", g_NX, 50., 110.);
  if(Cand) effs.push_back(eff_bosonMatched_UnmatchedV_genboson_Mass);
  TEfficiency* eff_bosonMatched_UnmatchedV_genboson_Eta = new TEfficiency((title+"_eff_bosonMatched_UnmatchedV_genboson_Eta").c_str(),"Efficiency of V Cand Not Matching to Gen Boson;genBosonEta;Efficiency", g_NX, -5., 5.);
  if(Cand) effs.push_back(eff_bosonMatched_UnmatchedV_genboson_Eta);

  TH2D* hist_MatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_EMR_CosDecayAngle_CM);
  TH2D* hist_PartialmatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_EMR_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_EMR_CosDecayAngle_CM);
  TH2D* hist_MatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_EPR_CosDecayAngle_CM);
  TH2D* hist_PartialmatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_EPR_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_EPR_CosDecayAngle_CM);
  TH2D* hist_MatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_EMR_other_CosDecayAngle_CM);
  TH2D* hist_PartialmatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_EMR_other_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_EMR_other_CosDecayAngle_CM);
  TH2D* hist_MatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_MatchedV_EPR_other_CosDecayAngle_CM);
  TH2D* hist_PartialmatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_PartialmatchedV_EPR_other_CosDecayAngle_CM);
  TH2D* hist_UnmatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
  if(Cand) hists2.push_back(hist_UnmatchedV_EPR_other_CosDecayAngle_CM);

  TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Cand) hists2.push_back(hist_RISR_PTISR);
  TH2D* hist_RISR_PTISR_Cands = new TH2D((title+"_RISR_PTISR_Cands").c_str(), (title+"_RISR_PTISR_Cands;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Cand) hists2.push_back(hist_RISR_PTISR_Cands);
  TH2D* hist_RISR_PTISR_ISR_BDT = new TH2D((title+"_RISR_PTISR_ISR_BDT").c_str(), (title+"_RISR_PTISR_ISR_BDT;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Cand) hists2.push_back(hist_RISR_PTISR_ISR_BDT);
  TH2D* hist_RISR_PTISR_BDT_ISR_lep = new TH2D((title+"_RISR_PTISR_BDT_ISR_lep").c_str(), (title+"_RISR_PTISR_BDT_ISR_lep;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Cand) hists2.push_back(hist_RISR_PTISR_BDT_ISR_lep);
  TH2D* hist_RISR_PTISR_BDT_ISR_singlet = new TH2D((title+"_RISR_PTISR_BDT_ISR_singlet").c_str(), (title+"_RISR_PTISR_BDT_ISR_singlet;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Cand) hists2.push_back(hist_RISR_PTISR_BDT_ISR_singlet);
  TH2D* hist_RISR_PTISR_ISR_Sparticle = new TH2D((title+"_RISR_PTISR_ISR_Sparticle").c_str(), (title+"_RISR_PTISR_ISR_Sparticle;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Sparticle1) hists2.push_back(hist_RISR_PTISR_ISR_Sparticle);

  // Sparticle2 KIN Plots
  TH2D* hist_RISR_PTISR_ISR_Sparticle2 = new TH2D((title+"_RISR_PTISR_ISR_Sparticle2").c_str(), (title+"_RISR_PTISR_ISR_Sparticle2;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
  if(Sparticle2) hists2.push_back(hist_RISR_PTISR_ISR_Sparticle2);
  TH2D* hist_dphiCMI_PTCM_ISR_Sparticle2 = new TH2D((title+"_dphiCMI_PTCM_ISR_Sparticle2").c_str(), (title+"_dphiCMI_PTCM_ISR_Sparticle2;dphiCMI;PTCM").c_str(), g_NX/2., 0., 3.15, g_NX/2., 0., 500.);
  if(Sparticle2) hists2.push_back(hist_dphiCMI_PTCM_ISR_Sparticle2);
  TH2D* hist_dphiMETV_PTISR_ISR_Sparticle2 = new TH2D((title+"_dphiMETV_PTISR_ISR_Sparticle2").c_str(), (title+"_dphiMETV_PTISR_ISR_Sparticle2;dphiMETV;PTISR").c_str(), g_NX/2., 0., 3.15, g_NX/2., 200., 800.);
  if(Sparticle2) hists2.push_back(hist_dphiMETV_PTISR_ISR_Sparticle2);
  TH2D* hist_gammaPerp_RISR_ISR_Sparticle2 = new TH2D((title+"_gammaPerp_RISR_ISR_Sparticle2").c_str(), (title+"_gammaPerp_RISR_ISR_Sparticle2;#gamma_{#perp};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
  if(Sparticle2) hists2.push_back(hist_gammaPerp_RISR_ISR_Sparticle2);
  TH2D* hist_Mperp_RISR_ISR_Sparticle2 = new TH2D((title+"_Mperp_RISR_ISR_Sparticle2").c_str(), (title+"_Mperp_RISR_ISR_Sparticle2;M_{#perp};RISR").c_str(), g_NX/2., 0., 100., g_NX/2., 0.5, 1.);
  if(Sparticle2) hists2.push_back(hist_Mperp_RISR_ISR_Sparticle2);
  TH2D* hist_RISR_mll_ISR_Sparticle2 = new TH2D((title+"_RISR_mll_ISR_Sparticle2").c_str(), (title+"_RISR_mll_ISR_Sparticle2;RISR;m_{ll}").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
  if(Sparticle2) hists2.push_back(hist_RISR_mll_ISR_Sparticle2);
  TH2D* hist_RISR_mllHEM_ISR_Sparticle2 = new TH2D((title+"_RISR_mllHEM_ISR_Sparticle2").c_str(), (title+"_RISR_mllHEM_ISR_Sparticle2;RISR;m_{ll}HEM").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.); // invariant mass of leps in same hemisphere
  if(Sparticle2) hists2.push_back(hist_RISR_mllHEM_ISR_Sparticle2);
  TH2D* hist_RISR_mSJ_ISR_Sparticle2 = new TH2D((title+"_RISR_mSJ_ISR_Sparticle2").c_str(), (title+"_RISR_mSJ_ISR_Sparticle2;RISR;m_{SJ}").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
  if(Sparticle2) hists2.push_back(hist_RISR_mSJ_ISR_Sparticle2);
  TH2D* hist_CosThetaCM_mSJ_ISR_Sparticle2 = new TH2D((title+"_CosThetaCM_mSJ_ISR_Sparticle2").c_str(), (title+"_CosThetaCM_mSJ_ISR_Sparticle2;Cos{#theta}_{CM};m_{SJ}").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 150.);
  if(Sparticle2) hists2.push_back(hist_CosThetaCM_mSJ_ISR_Sparticle2);
  TH2D* hist_RISR_mllLEAD_ISR_Sparticle2 = new TH2D((title+"_RISR_mllLEAD_ISR_Sparticle2").c_str(), (title+"_RISR_mllLEAD_ISR_Sparticle2;RISR;m_{ll}LEAD").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
  if(Sparticle2) hists2.push_back(hist_RISR_mllLEAD_ISR_Sparticle2);
  TH2D* hist_RISR_mL_ISR_Sparticle2 = new TH2D((title+"_RISR_mL_ISR_Sparticle2").c_str(), (title+"_RISR_mL_ISR_Sparticle2;RISR;mL").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.); // mass of leptonic system
  if(Sparticle2) hists2.push_back(hist_RISR_mL_ISR_Sparticle2);
  TH2D* hist_Ma_Mb_ISR_Sparticle2 = new TH2D((title+"_Ma_Mb_ISR_Sparticle2").c_str(), (title+"_Ma_Mb_ISR_Sparticle2;Ma;Mb").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 350.);
  if(Sparticle2) hists2.push_back(hist_Ma_Mb_ISR_Sparticle2);
  TH2D* hist_MLa_MLb_ISR_Sparticle2 = new TH2D((title+"_MLa_MLb_ISR_Sparticle2").c_str(), (title+"_MLa_MLb_ISR_Sparticle2;MLa;MLb").c_str(), g_NX/2., 0., 200., g_NX/2., 0., 200.);
  if(Sparticle2) hists2.push_back(hist_MLa_MLb_ISR_Sparticle2);
  TH2D* hist_RISR_DiffML_ISR_Sparticle2 = new TH2D((title+"_RISR_DiffML_ISR_Sparticle2").c_str(), (title+"_RISR_DiffML_ISR_Sparticle2;RISR;MLa-MLb").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
  if(Sparticle2) hists2.push_back(hist_RISR_DiffML_ISR_Sparticle2);
  TH2D* hist_RISR_MVis_ISR_Sparticle2 = new TH2D((title+"_RISR_MVis_ISR_Sparticle2").c_str(), (title+"_RISR_MVis_ISR_Sparticle2;RISR;MVis").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 400.);
  if(Sparticle2) hists2.push_back(hist_RISR_MVis_ISR_Sparticle2);
  TH2D* hist_RISR_MS_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MS_S_ISR_Sparticle2").c_str(), (title+"_RISR_MS_S_ISR_Sparticle2;RISR;MS^{S}").c_str(), g_NX/2., 0.5, 1.,g_NX/2., 0., 500.); // Mass of S system using four vectors evaluated in S frame
  if(Sparticle2) hists2.push_back(hist_RISR_MS_S_ISR_Sparticle2);
  TH2D* hist_RISR_MSa_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MSa_S_ISR_Sparticle2").c_str(), (title+"_RISR_MSa_S_ISR_Sparticle2;RISR;MSa^{S}").c_str(), g_NX/2., 0.5, 1.,g_NX/2., 0., 500.); // Mass of Sa system using four vectors evaluated in S frame
  if(Sparticle2) hists2.push_back(hist_RISR_MSa_S_ISR_Sparticle2);
  TH2D* hist_RISR_MSb_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MSb_S_ISR_Sparticle2").c_str(), (title+"_RISR_MSb_S_ISR_Sparticle2;RISR;MSb^{S}").c_str(), g_NX/2., 0.5, 1.,g_NX/2., 0., 250.); // Mass of Sb system using four vectors evaluated in S frame
  if(Sparticle2) hists2.push_back(hist_RISR_MSb_S_ISR_Sparticle2);
  TH2D* hist_MSa_S_MSb_S_ISR_Sparticle2 = new TH2D((title+"_MSa_S_MSb_S_ISR_Sparticle2").c_str(), (title+"_MSa_S_MSb_S_ISR_Sparticle2;MSa^{S};MSb^{S}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 250.); // Mass of Sb system using four vectors evaluated in S frame
  if(Sparticle2) hists2.push_back(hist_MSa_S_MSb_S_ISR_Sparticle2);
  TH2D* hist_RISR_MVisS_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MVisS_S_ISR_Sparticle2").c_str(), (title+"_RISR_MVisS_S_ISR_Sparticle2;RISR;MVisS^{S}").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 500.); // Mass of visible part of S system using four vectors evaluated in S frame ("subtracting out invisible part")
  if(Sparticle2) hists2.push_back(hist_RISR_MVisS_S_ISR_Sparticle2);
  TH2D* hist_RISR_MVisSa_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MVisSa_S_ISR_Sparticle2").c_str(), (title+"_RISR_MVisSa_S_ISR_Sparticle2;RISR;MVisSa^{S}").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 500.); // Mass of visible part of Sa system using four vectors evaluated in S frame ("subtracting out invisible part")
  if(Sparticle2) hists2.push_back(hist_RISR_MVisSa_S_ISR_Sparticle2);
  TH2D* hist_RISR_MVisSb_S_ISR_Sparticle2 = new TH2D((title+"_RISR_MVisSb_S_ISR_Sparticle2").c_str(), (title+"_RISR_MVisSb_S_ISR_Sparticle2;RISR;MVisSb^{S}").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 250.); // Mass of visible part of Sb system using four vectors evaluated in S frame ("subtracting out invisible part")
  if(Sparticle2) hists2.push_back(hist_RISR_MVisSb_S_ISR_Sparticle2);
  TH2D* hist_MVisSa_S_MVisSb_S_ISR_Sparticle2 = new TH2D((title+"_MVisSa_S_MVisSb_S_ISR_Sparticle2").c_str(), (title+"_MVisSa_S_MVisSb_S_ISR_Sparticle2;MVisSa^{S};MVisSb^{S}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 250.); // Mass of visible part of Sb system using four vectors evaluated in S frame ("subtracting out invisible part")
  if(Sparticle2) hists2.push_back(hist_MVisSa_S_MVisSb_S_ISR_Sparticle2);
  TH2D* hist_MSperp_RISR_ISR_Sparticle2 = new TH2D((title+"_MSperp_RISR_Sparticle2").c_str(), (title+"_MSperp_RISR_ISR_Sparticle2;MSperp;RISR").c_str(), g_NX/2., 0., 200., g_NX/2., 0., 1.);
  if(Sparticle2) hists2.push_back(hist_MSperp_RISR_ISR_Sparticle2);

  
  

  if(debug) cout << "Loading TChains" << endl;
  TChain* chain = new TChain("Events");
  chain->Add(ifile.c_str());

  if(debug) cout << "Making base class" << endl;
  SUSYNANOBase* base = new SUSYNANOBase(chain, era == "Run3"); // Run2
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
    
    //baseRun3->GetEntry(e);
    base->GetEntry(e);
    int mymod = (N1-N0)/10;
    if(mymod < 1)
      mymod = 1;
    if(e%mymod == 0)
      cout << " event = " << e << " : [" << N0 << " , " << N1 << "]" << endl;
    //if(debug) cout << "event: " << e << endl;
    double hweight = base->genWeight*weight*double(SKIP)/1.e7; // need 1.e7 factor to normalize from input!
    int N = base->nGenPart;
    int PDGID;
    int MP = 0;
    int MC = 0;
    if(proc.find("TChi") != std::string::npos){
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

    if(DM > 0.)
      if(MP != 0 && MC != 0 && ((MP - MC) > (DM + 1.) || (MP - MC) < (DM - 1.))) continue;

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
           if(base->Electron_mvaNoIso_WP80[i]){ // run3
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
           if(base->Electron_mvaFall17V2noIso_WP80[i]){ // run2 2017
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
      //if(!((abs(PDGID%1000) > 100 && abs(PDGID%1000) < 600)
      //        || (abs(PDGID%1000) > 1000 && abs(PDGID%1000) < 6000)
      //        || (abs(PDGID) > 0 && abs(PDGID) < 6)
      //        || PDGID == 21)) continue;
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
      //if(abs(momID) == 24 && base->GenPart_status[mom] == 22){
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
    if(Nleps < reco_lep_cut) continue;
    int Ngen_jets = int(gen_jets.size());
    int Ngen_leps = int(gen_leptons.size());

    int Ngen_sig_leps = 0; // leps coming from something interesting
    for(int i = 0; i < Ngen_leps; i++){
      if(gen_leptons[i].SourceID() == kSignal)
        Ngen_sig_leps++;
    }
    //if(gen_lepton_cut && Ngen_sig_leps != 1) continue; // only keep events with exactly one lepton coming from W
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

    // binary decay tree
    std::vector<V_Cand> V_had_cands;
    std::vector<V_Cand> V_lep_cands;
    LAB_BDT.ClearEvent();
    vector<int>   jet_BDT_singlet; // jets next to DecayFrame in BDT
    vector<int>   jet_BDT_nonsinglet; // jets NOT next to DecayFrame in BDT
    vector<int>   jet_BDT_lep; // jets next to lepton in BDT
    vector<int>   jet_BDT_nonlep; // jets NOT next to lepton in BDT
    vector<int>   jet_BDT_noncand; // jets not part of cand in BDT
    vector<RFKey> jetID_BDT;
    vector<RFKey> lepID_BDT;
    vector<RFKey> objID_BDT;
    if(BDT_MET)
      INV.SetLabFrameThreeVector(MET);

    // use reco jets for BDT
    if(!use_gen_jets){
      for(int i = 0; i < Njets; i++){
        if(Zero_RJR_BDT_JetMass){
          RFKey key = VIS.AddLabFrameFourVector(TLorentzVector(jets[i].Pt(),jets[i].Eta(),jets[i].Phi(),0.));
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
        else if(transverseOBJ){
          Particle temp = jets[i];
          temp.SetZ(0.);
          RFKey key = VIS.AddLabFrameFourVector(temp);
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
        else{
          RFKey key = VIS.AddLabFrameFourVector(jets[i]);
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
      }
    }

    // use gen jets for BDT
    if(use_gen_jets){
      for(int i = 0; i < Ngen_jets; i++){
        if(Zero_RJR_BDT_JetMass){
          RFKey key = VIS.AddLabFrameFourVector(TLorentzVector(gen_jets[i].Pt(),gen_jets[i].Eta(),gen_jets[i].Phi(),0.));
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
        else if(transverseOBJ){
          Particle temp = gen_jets[i];
          temp.SetZ(0.);
          RFKey key = VIS.AddLabFrameFourVector(temp);
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
        else{
          RFKey key = VIS.AddLabFrameFourVector(gen_jets[i]);
          jetID_BDT.push_back(key);
          objID_BDT.push_back(key);
        }
      }
    }


    for(int i = 0; i < Nleps; i++){
      if(transverseOBJ){
        Particle temp = leptons[i];
        temp.SetZ(0.);
        RFKey key = VIS.AddLabFrameFourVector(temp);
        lepID_BDT.push_back(key);
        objID_BDT.push_back(key);
      }
      else{
        RFKey key = VIS.AddLabFrameFourVector(leptons[i]);
        lepID_BDT.push_back(key);
        objID_BDT.push_back(key);
      }
    }
    if(!LAB_BDT.AnalyzeEvent()) cout << "Problem with RF BDT Analyze Event \n";

    int N_V_had = 0;
    //if(Cand){
      double ECM = 0.; // sum of energy of all objects in CM frame
      double ECM_LAB_BDT = 0.; // sum of energy of all objects in CM frame evalutated in LAB_BDT frame
      for(int i = 0; i < Nobjs; i++){
        ECM += BDT_CM.GetFrame(objID_BDT[i]).GetFourVector(BDT_CM).E();
        ECM_LAB_BDT += BDT_CM.GetFrame(objID_BDT[i]).GetFourVector().E();
      }
  
      if(use_gen_jets)
        V_had_cands = cand_list(gen_jets,leptons,VIS,BDT_CM,V1,BDT_CM,jetID_BDT,lepID_BDT);
      if(!use_gen_jets)
        V_had_cands = cand_list(jets,leptons,VIS,BDT_CM,V1,BDT_CM,jetID_BDT,lepID_BDT);
  
      N_V_had = V_had_cands.size();
      totcands += N_V_had;
      //if(N_V_had < 1) continue;
      for(int i = 0; i < N_V_had; i++){
        if(V_had_cands[i].Type() == kSib)
          prong2++;
        else if(V_had_cands[i].Type() == kAunt)
          auntprong++;
        if(V_had_cands[i].Type() == kLep)
          lepprong++;
      }
  
      if(!use_gen_jets){
      
        for(int i = 0; i < Njets; i++){
          const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[i]);
          const RestFrame& sibling = frame.GetSiblingFrame();
          jet_BDT_nonsinglet.push_back(i);
          jet_BDT_nonlep.push_back(i);
          if(sibling.IsDecayFrame()){
            jet_BDT_singlet.push_back(i);
            jet_BDT_nonsinglet.pop_back();
          }
          if(int(lepID_BDT.size()) == Nleps){
            for(int j = 0; j < Nleps; j++){
              if(sibling == BDT_CM.GetFrame(lepID_BDT[j])){
                jet_BDT_lep.push_back(i);
                jet_BDT_nonlep.pop_back();
              } // if(sibling == BDT_CM.GetFrame(lepID_BDT[j]))
            } // for(int j = 0; j < Nleps; j++)
          } // if(int(lepID_BDT.size()) == Nleps)
        bool noncand = true;
        for(int j = 0; j < N_V_had; j++){
          if(V_had_cands[j].IsProng(frame)){
            noncand = false;
            break;
          }
        }
        if(noncand)
          jet_BDT_noncand.push_back(i);
        } // for(int i = 0; i < Njets; i++)
      }  // if(!use_gen_jets)
  
      cand_side(V_had_cands, BDT_CM);
      cand_matching(V_had_cands);
      hist_CandCount->SetBinContent(0,hist_CandCount->GetBinContent(0)+Nhadbosons);
      for(int i = 0; i < N_V_had; i++){
        bool matched = false;
        if(V_had_cands[i].Match() == kMatched){
          hist_CandCount->SetBinContent(1,hist_CandCount->GetBinContent(1)+1);
          matched = true;
        }
        else{
          if(V_had_cands[i].Match() == kW){
            hist_CandCount->SetBinContent(2,hist_CandCount->GetBinContent(2)+1);
          } // kW
          else if(V_had_cands[i].Match() == kZ){
            hist_CandCount->SetBinContent(3,hist_CandCount->GetBinContent(3)+1);
          } // kZ
          else if(V_had_cands[i].Match() == kB){
            hist_CandCount->SetBinContent(4,hist_CandCount->GetBinContent(4)+1);
          } // kB
          else{
            hist_CandCount->SetBinContent(5,hist_CandCount->GetBinContent(5)+1);
          } // else
        } // else
        eff_GenJetMatchedAll_V_cand_eta->Fill(matched, V_had_cands[i].Eta(), hweight);
        eff_GenJetMatchedAll_V_cand_mass->Fill(matched, V_had_cands[i].M(), hweight);
      } // for(int i = 0; i < N_V_had; i++)
      
      for(int i = 0; i < Nhadbosons; i++){
        bool matched = false;
        for(int j = 0; j < N_V_had; j++){
          if(V_had_cands[j].Match() == kMatched){
            if(gen_hadbosons[i].GenIndex() == V_had_cands[j].PL()[0].GenMomIndex()){
              matched = true;
              break;
            }
          }
          else continue;
        }
        eff_bosonMatched_MatchedV_genboson_Mass->Fill(matched, gen_hadbosons[i].M(), hweight);
        eff_bosonMatched_MatchedV_genboson_Eta->Fill(matched, gen_hadbosons[i].Eta(), hweight);
        eff_bosonMatched_UnmatchedV_genboson_Mass->Fill(!matched, gen_hadbosons[i].M(), hweight);
        eff_bosonMatched_UnmatchedV_genboson_Eta->Fill(!matched, gen_hadbosons[i].Eta(), hweight);
      }
  
      std::map<int,int> jet_pairs;
      if(use_gen_jets){
        for(int i = 0; i < Ngen_jets; i++){
          for(int j = i+1; j < Ngen_jets; j++){
            //if(abs(gen_jets[i].MomPDGID()) == 24 && gen_jets[i].GenMomIndex() == gen_jets[j].GenMomIndex()){
            if((abs(gen_jets[i].MomPDGID()) == 23 || abs(gen_jets[i].MomPDGID()) == 24) && gen_jets[i].GenMomIndex() == gen_jets[j].GenMomIndex()){
              jet_pairs.insert(std::make_pair(i,j));
            }
          }
        }
      }
      if(!use_gen_jets){
        for(int i = 0; i < Njets; i++){
          for(int j = i+1; j < Njets; j++){
            if((abs(jets[i].MomPDGID()) == 23 || abs(jets[i].MomPDGID()) == 24) && jets[i].GenMomIndex() == jets[j].GenMomIndex()){
              jet_pairs.insert(std::make_pair(i,j));
            }
          }
        }
      }
  
      for(std::map<int,int>::iterator jet_pair = jet_pairs.begin(); jet_pair != jet_pairs.end(); ++jet_pair){
        const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[jet_pair->first]);
        const RestFrame& siblingFrame = BDT_CM.GetFrame(jetID_BDT[jet_pair->second]);
        int BaconNum = FrameDistance(frame, siblingFrame);
        hist_BaconNum->Fill(BaconNum, hweight);
      }
  
      int BDT_matched = 0;
      int BDT_unmatched = 0;
      int BDT_partial = 0;
      const RestFrame& Prod_Frame = BDT_CM;
      for(int i = 0; i < N_V_had; i++){
        if(V_had_cands[i].Match() == kMatched){
          BDT_matched++;
          if(V_had_cands[i].Side() == kAside)
            matchedAside++;
          else
            matchedBside++;
          for(int j = 0; j < Nleps; j++){
            hist_MatchedlepCandBaconNum->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),hweight);
            hist_MatchedlepCandBaconNum_Njets->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),Njets,hweight);
          }
          for(int j = 0; j < int(V_had_cands[i].size()); j++){
            hist_MatchedVJetsPt->Fill(V_had_cands[i][j].Pt(), hweight);
            hist_MatchedVJetsEta->Fill(V_had_cands[i][j].Eta(), hweight);
            hist_MatchedVJetsMass->Fill(V_had_cands[i][j].M(), hweight);
          }
          hist_MatchedVPt->Fill(V_had_cands[i].Pt(), hweight);
          hist_MatchedVEta->Fill(V_had_cands[i].Eta(), hweight);
          hist_MatchedVMass->Fill(V_had_cands[i].M(), hweight);
          TLorentzVector P_1_P = V_had_cands[i].RL().Get(0).GetFourVector(Prod_Frame);
          TLorentzVector P_2_P = V_had_cands[i].RL().Get(1).GetFourVector(Prod_Frame);
          TLorentzVector P_W_P = P_1_P + P_2_P;
          TVector3 Beta_W_P = P_W_P.BoostVector();
          P_1_P.Boost(-Beta_W_P);
          double CosDecayAngle = P_1_P.Vect().Unit().Dot(Beta_W_P.Unit());
          double CosDecayAngle_CM = BDT_CM.GetCosDecayAngle(V_had_cands[i].CandFrame());
          double DeltaPhiDecayAngle_CM = BDT_CM.GetDeltaPhiDecayAngle(RestFrame::GetAxis(), V_had_cands[i].CandFrame());
          double DeltaPhiDecayPlanes_CM = BDT_CM.GetDeltaPhiDecayPlanes(V_had_cands[i].CandFrame());
          double EMR = V_had_cands[i].E()/V_had_cands[i].M();
          double EPR = V_had_cands[i].E()/V_had_cands[i].P();
          double EMR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector() - V_had_cands[i].CandFrame().GetFourVector()).Mag();
          double EPR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector().Vect() - V_had_cands[i].CandFrame().GetFourVector().Vect()).Mag();
          hist_MatchedVCosDecayAngle->Fill(CosDecayAngle, hweight);
          hist_MatchedV_Mass_CosDecayAngle->Fill(V_had_cands[i].M(), CosDecayAngle, hweight);
          hist_MatchedV_CosDecayAngle_CM->Fill(CosDecayAngle_CM, hweight);
          hist_MatchedV_DeltaPhiDecayAngle_CM->Fill(DeltaPhiDecayAngle_CM, hweight);
          hist_MatchedV_DeltaPhiDecayPlanes_CM->Fill(DeltaPhiDecayPlanes_CM, hweight);
          hist_MatchedV_Mass_CosDecayAngle_CM->Fill(V_had_cands[i].M(), CosDecayAngle_CM, hweight);
          hist_MatchedV_Mass_DeltaPhiDecayAngle_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayAngle_CM, hweight);
          hist_MatchedV_Mass_DeltaPhiDecayPlanes_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayPlanes_CM, hweight);
          hist_MatchedVProngDeltaPhi->Fill(V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_MatchedV_Mass_ProngDeltaPhi->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_MatchedV_CosDecayAngle_CM_ProngDeltaPhi->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_MatchedVProngDeltaEta->Fill(V_had_cands[i].ProngDeltaEta(), hweight);
          hist_MatchedV_Mass_ProngDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaEta(), hweight);
          hist_MatchedV_CosDecayAngle_CM_ProngDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaEta(), hweight);
          hist_MatchedVProngAbsDeltaEta->Fill(V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_MatchedV_Mass_ProngAbsDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM->Fill(DeltaPhiDecayAngle_CM, CosDecayAngle_CM, hweight);
          hist_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM->Fill(DeltaPhiDecayPlanes_CM, CosDecayAngle_CM, hweight);
          hist_MatchedV_CosDecayAngle_CosDecayAngle_CM->Fill(CosDecayAngle, CosDecayAngle_CM, hweight);
          hist_Matched_DeltaR_P->Fill(V_had_cands[i].ProngDeltaR(), V_had_cands[i].P(), hweight);
          hist_Matched_lnDeltaR_lnP->Fill(log(V_had_cands[i].ProngDeltaR()), log(V_had_cands[i].P()), hweight);
          hist_Matched_Mass_P_Multi_DeltaR->Fill(V_had_cands[i].M(), V_had_cands[i].P()*V_had_cands[i].ProngDeltaR(), hweight);
          hist_Matched_CosDecayAngle_CM_P_Multi_DeltaR->Fill(CosDecayAngle_CM, V_had_cands[i].P()*V_had_cands[i].ProngDeltaR(), hweight);
          hist_MatchedV_PCM_PCand->Fill(BDT_CM.GetFourVector().P(), V_had_cands[i].P(), hweight);
          hist_MatchedCand_JetA_Mass_JetB_Mass->Fill(V_had_cands[i].PL()[0].M(), V_had_cands[i].PL()[1].M(), hweight);
          hist_MatchedCand_PMR_CosDecayAngle_CM->Fill(V_had_cands[i].ProngMassRatio(), CosDecayAngle_CM, hweight);
          hist_MatchedCand_PMR_Mass->Fill(V_had_cands[i].ProngMassRatio(), V_had_cands[i].M(), hweight);
          hist_MatchedV_EMR_CosDecayAngle_CM->Fill(EMR, CosDecayAngle_CM, hweight);
          hist_MatchedV_EPR_CosDecayAngle_CM->Fill(EPR, CosDecayAngle_CM, hweight);
          hist_MatchedV_EMR_other_CosDecayAngle_CM->Fill(EMR_other, CosDecayAngle_CM, hweight);
          hist_MatchedV_EPR_other_CosDecayAngle_CM->Fill(EPR_other, CosDecayAngle_CM, hweight);
          if(objID_BDT.size() > 5){
            if(V_had_cands[i].Side() == kAside){
              hist_Matched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Matched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
            else{
              hist_Matched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Matched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
          }
        }
        else if(V_had_cands[i].Match() == kUnmatched){
          BDT_unmatched++;
          if(V_had_cands[i].Side() == kAside)
            unmatchedAside++;
          else
            unmatchedBside++;
          for(int j = 0; j < Nleps; j++){
            hist_UnmatchedlepCandBaconNum->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),hweight);
            hist_UnmatchedlepCandBaconNum_Njets->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),Njets,hweight);
          }
          for(int j = 0; j < int(V_had_cands[i].size()); j++){
           hist_UnmatchedVJetsPt->Fill(V_had_cands[i][j].Pt(), hweight);
           hist_UnmatchedVJetsEta->Fill(V_had_cands[i][j].Eta(), hweight);
           hist_UnmatchedVJetsMass->Fill(V_had_cands[i][j].M(), hweight);
          }
          hist_UnmatchedVPt->Fill(V_had_cands[i].Pt(), hweight);
          hist_UnmatchedVEta->Fill(V_had_cands[i].Eta(), hweight);
          hist_UnmatchedVMass->Fill(V_had_cands[i].M(), hweight);
          TLorentzVector P_1_P = V_had_cands[i].RL().Get(0).GetFourVector(Prod_Frame);
          TLorentzVector P_2_P = V_had_cands[i].RL().Get(1).GetFourVector(Prod_Frame);
          TLorentzVector P_W_P = P_1_P + P_2_P;
          TVector3 Beta_W_P = P_W_P.BoostVector();
          P_1_P.Boost(-Beta_W_P);
          double CosDecayAngle = P_1_P.Vect().Unit().Dot(Beta_W_P.Unit());
          double CosDecayAngle_CM = BDT_CM.GetCosDecayAngle(V_had_cands[i].CandFrame());
          double DeltaPhiDecayAngle_CM = BDT_CM.GetDeltaPhiDecayAngle(RestFrame::GetAxis(), V_had_cands[i].CandFrame());
          double DeltaPhiDecayPlanes_CM = BDT_CM.GetDeltaPhiDecayPlanes(V_had_cands[i].CandFrame());
          double EMR = V_had_cands[i].E()/V_had_cands[i].M();
          double EPR = V_had_cands[i].E()/V_had_cands[i].P();
          double EMR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector() - V_had_cands[i].CandFrame().GetFourVector()).Mag();
          double EPR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector().Vect() - V_had_cands[i].CandFrame().GetFourVector().Vect()).Mag();
          hist_UnmatchedVCosDecayAngle->Fill(CosDecayAngle, hweight);
          hist_UnmatchedV_Mass_CosDecayAngle->Fill(V_had_cands[i].M(), CosDecayAngle, hweight);
          hist_UnmatchedV_CosDecayAngle_CM->Fill(CosDecayAngle_CM, hweight);
          hist_UnmatchedV_DeltaPhiDecayAngle_CM->Fill(DeltaPhiDecayAngle_CM, hweight);
          hist_UnmatchedV_DeltaPhiDecayPlanes_CM->Fill(DeltaPhiDecayPlanes_CM, hweight);
          hist_UnmatchedV_Mass_CosDecayAngle_CM->Fill(V_had_cands[i].M(), CosDecayAngle_CM, hweight);
          hist_UnmatchedV_Mass_DeltaPhiDecayAngle_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayAngle_CM, hweight);
          hist_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayPlanes_CM, hweight);
          hist_UnmatchedVProngDeltaPhi->Fill(V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_UnmatchedV_Mass_ProngDeltaPhi->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_UnmatchedVProngDeltaEta->Fill(V_had_cands[i].ProngDeltaEta(), hweight);
          hist_UnmatchedV_Mass_ProngDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaEta(), hweight);
          hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaEta(), hweight);
          hist_UnmatchedVProngAbsDeltaEta->Fill(V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_UnmatchedV_Mass_ProngAbsDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_UnmatchedV_CosDecayAngle_CosDecayAngle_CM->Fill(CosDecayAngle, CosDecayAngle_CM, hweight);
          hist_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM->Fill(DeltaPhiDecayAngle_CM, CosDecayAngle_CM, hweight);
          hist_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM->Fill(DeltaPhiDecayPlanes_CM, CosDecayAngle_CM, hweight);
          hist_Unmatched_DeltaR_P->Fill(V_had_cands[i].ProngDeltaR(), V_had_cands[i].P(), hweight);
          hist_Unmatched_lnDeltaR_lnP->Fill(log(V_had_cands[i].ProngDeltaR()), log(V_had_cands[i].P()), hweight);
          hist_Unmatched_Mass_P_Multi_DeltaR->Fill(V_had_cands[i].M(), V_had_cands[i].P()*V_had_cands[i].ProngDeltaR(), hweight);
          hist_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR->Fill(CosDecayAngle_CM, V_had_cands[i].P()*V_had_cands[i].ProngDeltaR(), hweight);
          hist_UnmatchedV_PCM_PCand->Fill(BDT_CM.GetMomentum(LAB_BDT), V_had_cands[i].P(), hweight);
          hist_UnmatchedCand_JetA_Mass_JetB_Mass->Fill(V_had_cands[i].PL()[0].M(), V_had_cands[i].PL()[1].M(), hweight);
          hist_UnmatchedCand_PMR_CosDecayAngle_CM->Fill(V_had_cands[i].ProngMassRatio(), CosDecayAngle_CM, hweight);
          hist_UnmatchedCand_PMR_Mass->Fill(V_had_cands[i].ProngMassRatio(), V_had_cands[i].M(), hweight);
          hist_UnmatchedV_EMR_CosDecayAngle_CM->Fill(EMR, CosDecayAngle_CM, hweight);
          hist_UnmatchedV_EPR_CosDecayAngle_CM->Fill(EPR, CosDecayAngle_CM, hweight);
          hist_UnmatchedV_EMR_other_CosDecayAngle_CM->Fill(EMR_other, CosDecayAngle_CM, hweight);
          hist_UnmatchedV_EPR_other_CosDecayAngle_CM->Fill(EPR_other, CosDecayAngle_CM, hweight);
          if(objID_BDT.size() > 5){
            if(V_had_cands[i].Side() == kAside){
              hist_Unmatched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Unmatched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
            else{
              hist_Unmatched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Unmatched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
          }
        }
        else{
          BDT_partial++;
          for(int j = 0; j < Nleps; j++){
            hist_PartialMatchedlepCandBaconNum->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),hweight);
            hist_PartialMatchedlepCandBaconNum_Njets->Fill(FrameDistance(V_had_cands[i].CandFrame(),BDT_CM.GetFrame(lepID_BDT[j])),Njets,hweight);
          }
          hist_PartialmatchedVMass->Fill(V_had_cands[i].M(), hweight);
          TLorentzVector P_1_P = V_had_cands[i].RL().Get(0).GetFourVector(Prod_Frame);
          TLorentzVector P_2_P = V_had_cands[i].RL().Get(1).GetFourVector(Prod_Frame);
          TLorentzVector P_W_P = P_1_P + P_2_P;
          TVector3 Beta_W_P = P_W_P.BoostVector();
          P_1_P.Boost(-Beta_W_P);
          double CosDecayAngle = P_1_P.Vect().Unit().Dot(Beta_W_P.Unit());
          double CosDecayAngle_CM = BDT_CM.GetCosDecayAngle(V_had_cands[i].CandFrame());
          double DeltaPhiDecayAngle_CM = BDT_CM.GetDeltaPhiDecayAngle(RestFrame::GetAxis(), V_had_cands[i].CandFrame());
          double DeltaPhiDecayPlanes_CM = BDT_CM.GetDeltaPhiDecayPlanes(V_had_cands[i].CandFrame());
          double EMR = V_had_cands[i].E()/V_had_cands[i].M();
          double EPR = V_had_cands[i].E()/V_had_cands[i].P();
          double EMR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector() - V_had_cands[i].CandFrame().GetFourVector()).Mag();
          double EPR_other = (ECM_LAB_BDT - V_had_cands[i].E()) / (BDT_CM.GetFourVector().Vect() - V_had_cands[i].CandFrame().GetFourVector().Vect()).Mag();
          hist_PartialmatchedV_Mass_CosDecayAngle_CM->Fill(V_had_cands[i].M(), CosDecayAngle_CM, hweight);
          hist_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayAngle_CM, hweight);
          hist_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM->Fill(V_had_cands[i].M(), DeltaPhiDecayPlanes_CM, hweight);
          hist_PartialmatchedVProngDeltaPhi->Fill(V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_PartialmatchedV_Mass_ProngDeltaPhi->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaPhi(), hweight);
          hist_PartialmatchedVProngDeltaEta->Fill(V_had_cands[i].ProngDeltaEta(), hweight);
          hist_PartialmatchedV_Mass_ProngDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngDeltaEta(), hweight);
          hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngDeltaEta(), hweight);
          hist_PartialmatchedVProngAbsDeltaEta->Fill(V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_PartialmatchedV_Mass_ProngAbsDeltaEta->Fill(V_had_cands[i].M(), V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta->Fill(CosDecayAngle_CM, V_had_cands[i].ProngAbsDeltaEta(), hweight);
          hist_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM->Fill(CosDecayAngle, CosDecayAngle_CM, hweight);
          hist_PartialmatchedV_PCM_PCand->Fill(BDT_CM.GetMomentum(LAB_BDT), V_had_cands[i].P(), hweight);
          hist_PartialMatchedCand_JetA_Mass_JetB_Mass->Fill(V_had_cands[i].PL()[0].M(), V_had_cands[i].PL()[1].M(), hweight);
          hist_PartialMatchedCand_PMR_CosDecayAngle_CM->Fill(V_had_cands[i].ProngMassRatio(), CosDecayAngle_CM, hweight);
          hist_PartialMatchedCand_PMR_Mass->Fill(V_had_cands[i].ProngMassRatio(), V_had_cands[i].M(), hweight);
          hist_PartialmatchedV_EMR_CosDecayAngle_CM->Fill(EMR, CosDecayAngle_CM, hweight);
          hist_PartialmatchedV_EPR_CosDecayAngle_CM->Fill(EPR, CosDecayAngle_CM, hweight);
          hist_PartialmatchedV_EMR_other_CosDecayAngle_CM->Fill(EMR_other, CosDecayAngle_CM, hweight);
          hist_PartialmatchedV_EPR_other_CosDecayAngle_CM->Fill(EPR_other, CosDecayAngle_CM, hweight);
          if(objID_BDT.size() > 5){
            if(V_had_cands[i].Side() == kAside){
              hist_Partialmatched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Partialmatched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(0).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
            else{
              hist_Partialmatched_CosDecayAngle_Hemi->Fill(BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
              hist_Partialmatched_Mass_CosDecayAngle_Hemi->Fill(V_had_cands[i].M(), BDT_CM.GetChildFrame(1).GetCosDecayAngle(V_had_cands[i].CandFrame()), hweight);
            }
          }
        }
      }
  
      hist_V_had_cand->Fill(N_V_had, hweight);
      hist_Njets_V_had_cand->Fill(Njets, N_V_had, hweight);
      hist_Njets_matched_V_had_cand->Fill(Njets, BDT_matched); 
      hist_Njets_unmatched_V_had_cand->Fill(Njets, BDT_unmatched); 
      hist_Njets_partial_V_had_cand->Fill(Njets, BDT_partial); 
  
      //-/-/-/-/-/-/-/-/-//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
      // chi2 min bs
      // gen jets
      if(use_gen_jets){
        hist_CandCount_Chi2->SetBinContent(0,hist_CandCount_Chi2->GetBinContent(0)+Nhadbosons);
        if(Ngen_jets > 3){
          int X2_j1 = -1;
          int X2_j2 = -1;
          getJetsMinX2(X2_j1, X2_j2, gen_jets);
          if(X2_j1 != X2_j2 && X2_j1 >= 0){
            TLorentzVector j1 = gen_jets[X2_j1];
            TLorentzVector j2 = gen_jets[X2_j2];
            double mJJ = (j1 + j2).M();
            hist_chi2WMass->Fill(mJJ, hweight);
            //if(abs(gen_jets[X2_j1].MomPDGID()) == 24 && gen_jets[X2_j1].GenMomIndex() == gen_jets[X2_j2].GenMomIndex()){
            if((abs(gen_jets[X2_j1].MomPDGID()) == 23 || abs(gen_jets[X2_j1].MomPDGID()) == 24) && gen_jets[X2_j1].GenMomIndex() == gen_jets[X2_j2].GenMomIndex()){
              hist_CandCount_Chi2->SetBinContent(1,hist_CandCount_Chi2->GetBinContent(1)+1);
            }
            else{
              if(abs(gen_jets[X2_j1].MomPDGID()) == 24 || abs(gen_jets[X2_j2].MomPDGID()) == 24){
                hist_CandCount_Chi2->SetBinContent(2,hist_CandCount_Chi2->GetBinContent(2)+1);
              }
              else if(abs(gen_jets[X2_j1].MomPDGID()) == 23 || abs(gen_jets[X2_j2].MomPDGID()) == 23){
                hist_CandCount_Chi2->SetBinContent(3,hist_CandCount_Chi2->GetBinContent(3)+1);
              }
              else if(abs(gen_jets[X2_j1].MomPDGID()) == 6 || abs(gen_jets[X2_j2].MomPDGID()) == 6){
                hist_CandCount_Chi2->SetBinContent(4,hist_CandCount_Chi2->GetBinContent(4)+1);
              }
              else{
                hist_CandCount_Chi2->SetBinContent(5,hist_CandCount_Chi2->GetBinContent(5)+1);
              }
            }
          }
        }
      }
      // end gen jets
      // reco jets
      if(!use_gen_jets){
        hist_CandCount_Chi2->SetBinContent(0,hist_CandCount_Chi2->GetBinContent(0)+Nhadbosons);
        if(Njets > 3){
          int X2_j1 = -1;
          int X2_j2 = -1;
          getJetsMinX2(X2_j1, X2_j2, jets);
          if(X2_j1 != X2_j2 && X2_j1 >= 0){
            TLorentzVector j1 = jets[X2_j1];
            TLorentzVector j2 = jets[X2_j2];
            double mJJ = (j1 + j2).M();
            hist_chi2WMass->Fill(mJJ, hweight);
            //if(abs(jets[X2_j1].MomPDGID()) == 24 && jets[X2_j1].GenMomIndex() == jets[X2_j2].GenMomIndex()){
            if((abs(jets[X2_j1].MomPDGID()) == 23 || abs(jets[X2_j1].MomPDGID()) == 24) && jets[X2_j1].GenMomIndex() == jets[X2_j2].GenMomIndex()){
              hist_CandCount_Chi2->SetBinContent(1,hist_CandCount_Chi2->GetBinContent(1)+1);
            }
            else{
              if(abs(jets[X2_j1].MomPDGID()) == 24 || abs(jets[X2_j2].MomPDGID()) == 24){
                hist_CandCount_Chi2->SetBinContent(2,hist_CandCount_Chi2->GetBinContent(2)+1);
              }
              else if(abs(jets[X2_j1].MomPDGID()) == 23 || abs(jets[X2_j2].MomPDGID()) == 23){
                hist_CandCount_Chi2->SetBinContent(3,hist_CandCount_Chi2->GetBinContent(3)+1);
              }
              else if(abs(jets[X2_j1].MomPDGID()) == 6 || abs(jets[X2_j2].MomPDGID()) == 6){
                hist_CandCount_Chi2->SetBinContent(4,hist_CandCount_Chi2->GetBinContent(4)+1);
              }
              else{
                hist_CandCount_Chi2->SetBinContent(5,hist_CandCount_Chi2->GetBinContent(5)+1);
              }
            }
          }
        }
      }
      // end reco jets
      //-/-/-/-/-/-/-/-/-//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
  
      std::vector<V_Cand> V_had_cands_combo; // pairwise combine all jet combinations
  
      if(use_gen_jets){
        for(int i = 0; i < Ngen_jets; i++){
          break; // speed up code by skipping this
          for(int j = i+1; j < Ngen_jets; j++){
            ParticleList V_cand_Part;
            ConstRestFrameList V_cand_RF;
            const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[i]);
            const RestFrame& sibling = BDT_CM.GetFrame(jetID_BDT[j]);
            V_cand_Part.push_back(gen_jets[i]);
            V_cand_Part.push_back(gen_jets[j]); // 2-prong "sibling"
            V_cand_RF.Add(frame);
            V_cand_RF.Add(sibling);
            V_had_cands_combo.push_back(V_Cand(V_cand_Part,V_cand_RF));
          }
        }
      }
      if(!use_gen_jets){
        for(int i = 0; i < Njets; i++){
          break; // speed up code by skipping this
          for(int j = i+1; j < Njets; j++){
            ParticleList V_cand_Part;
            ConstRestFrameList V_cand_RF;
            const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[i]);
            const RestFrame& sibling = BDT_CM.GetFrame(jetID_BDT[j]);
            V_cand_Part.push_back(jets[i]);
            V_cand_Part.push_back(jets[j]); // 2-prong "sibling"
            V_cand_RF.Add(frame);
            V_cand_RF.Add(sibling);
            V_had_cands_combo.push_back(V_Cand(V_cand_Part,V_cand_RF));
          }
        }
      }
  
      int N_V_had_combo = V_had_cands_combo.size();
      genHBosons+=Nhadbosons;
      hist_CandCount_Combo->SetBinContent(0,hist_CandCount_Combo->GetBinContent(0)+Nhadbosons);
      cand_matching(V_had_cands_combo);
      for(int i = 0; i < N_V_had_combo; i++){
        if(V_had_cands_combo[i].Match() == kMatched){
          hist_CandCount_Combo->SetBinContent(1,hist_CandCount_Combo->GetBinContent(1)+1);
        }
        else{
          if(V_had_cands_combo[i].Match() == kW){
            hist_CandCount_Combo->SetBinContent(2,hist_CandCount_Combo->GetBinContent(2)+1);
          } // kW
          else if(V_had_cands_combo[i].Match() == kZ){
            hist_CandCount_Combo->SetBinContent(3,hist_CandCount_Combo->GetBinContent(3)+1);
          } // kZ
          else if(V_had_cands_combo[i].Match() == kB){
            hist_CandCount_Combo->SetBinContent(4,hist_CandCount_Combo->GetBinContent(4)+1);
          } // kB
          else{
            hist_CandCount_Combo->SetBinContent(5,hist_CandCount_Combo->GetBinContent(5)+1);
          } // else
        } // else
      } // for(int i = 0; i < N_V_had_combo; i++)
  
      for(int i = 0; i < N_V_had_combo; i++){
        if(V_had_cands_combo[i].Match() == kMatched){
          hist_combo_Matched_DeltaR_P->Fill(V_had_cands_combo[i].ProngDeltaR(), V_had_cands_combo[i].P(), hweight);
        }
        else if(V_had_cands_combo[i].Match() == kUnmatched){
          hist_combo_Unmatched_DeltaR_P->Fill(V_had_cands_combo[i].ProngDeltaR(), V_had_cands_combo[i].P(), hweight);
        }
      }
  
    //} // if(Cand)

    // Original RJR Tree
    LAB.ClearEvent();
    INV.SetLabFrameThreeVector(MET);

    vector<RFKey> jetID;
    vector<RFKey> lepID;
    vector<RFKey> objID;
    // use reco jets for RJR
    if(!use_gen_jets){
      for(int i = 0; i < Njets; i++){
        RFKey key = COMB_J.AddLabFrameFourVector(jets[i]);
        jetID.push_back(key);
        objID.push_back(key);
      }
    }
    // use gen jets for RJR
    if(use_gen_jets){
      for(int i = 0; i < Ngen_jets; i++){
        RFKey key = COMB_J.AddLabFrameFourVector(gen_jets[i]);
        jetID.push_back(key);
        objID.push_back(key);
      }
    }
    for(int i = 0; i < Nleps; i++){
      RFKey key = COMB_L.AddLabFrameFourVector(leptons[i]);
      lepID.push_back(key);
      objID.push_back(key);
    }
    int Njet_S = 0;
    TVector3 vPISR(0.,0.,0.);
    double PTISR = 0.;
    TVector3 vPINV(0.,0.,0.);
    double RISR = 0.;
    //if(Cand){
      if(!LAB.AnalyzeEvent()) cout << "Problem with RF Analyze Event! Njets: " << jetID.size() << " Nleps: " << lepID.size() << " MET: " << MET.Mag() << endl;
      for(int i = 0; i < Njets; i++)
        if(COMB_J.GetFrame(jetID[i]) == Ja || COMB_J.GetFrame(jetID[i]) == Jb)
          Njet_S++;
      if(Nleps + Njet_S < 1){
        vPISR = S.GetFourVector(CM).Vect();
        PTISR = S.GetTransverseFourVector(CM).Vect().Mag();
        vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();
        RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
        hist_RISR_PTISR->Fill(RISR, PTISR, hweight);
      }
    //} // if(Cand)

    // binary decay tree with S & ISR splitting
    LAB_ISR_Sparticle.ClearEvent();
    vector<RFKey> jetID_BDT_ISR_Sparticle;
    vector<RFKey> lepID_BDT_ISR_Sparticle;
    vector<RFKey> objID_BDT_ISR_Sparticle;
    if(BDT_MET || S_MET)
      INV_ISR_Sparticle.SetLabFrameThreeVector(MET);

    // use reco jets for BDT
    if(!use_gen_jets){
      for(int i = 0; i < Njets; i++){
        if(Zero_RJR_BDT_JetMass){
          RFKey key = VIS_ISR_Sparticle.AddLabFrameFourVector(TLorentzVector(jets[i].Pt(),jets[i].Eta(),jets[i].Phi(),0.));
          jetID_BDT_ISR_Sparticle.push_back(key);
          objID_BDT_ISR_Sparticle.push_back(key);
        }
        Particle temp = jets[i];
        temp.SetZ(0.);
        if(S_leps){
          if(transverseOBJ){
            RFKey key = COMB_J_ISR_Sparticle.AddLabFrameFourVector(temp);
            jetID_BDT_ISR_Sparticle.push_back(key);
            objID_BDT_ISR_Sparticle.push_back(key);
          }
          else{
            RFKey key = COMB_J_ISR_Sparticle.AddLabFrameFourVector(jets[i]);
            jetID_BDT_ISR_Sparticle.push_back(key);
            objID_BDT_ISR_Sparticle.push_back(key);
          }
        }
        else{
          if(transverseOBJ){
            RFKey key = VIS_ISR_Sparticle.AddLabFrameFourVector(temp);
            jetID_BDT_ISR_Sparticle.push_back(key);
            objID_BDT_ISR_Sparticle.push_back(key);
          }
          else{
            RFKey key = VIS_ISR_Sparticle.AddLabFrameFourVector(jets[i]);
            jetID_BDT_ISR_Sparticle.push_back(key);
            objID_BDT_ISR_Sparticle.push_back(key);
          }
        }
      }
    }

    for(int i = 0; i < Nleps; i++){
      Particle temp = leptons[i];
      temp.SetZ(0.);
      if(S_leps){
        if(transverseOBJ){
          RFKey key = COMB_L_ISR_Sparticle.AddLabFrameFourVector(temp);
          lepID_BDT_ISR_Sparticle.push_back(key);
          objID_BDT_ISR_Sparticle.push_back(key);
        }
        else{
          RFKey key = COMB_L_ISR_Sparticle.AddLabFrameFourVector(leptons[i]);
          lepID_BDT_ISR_Sparticle.push_back(key);
          objID_BDT_ISR_Sparticle.push_back(key);
        }
      }
      else{
        if(transverseOBJ){
          RFKey key = VIS_ISR_Sparticle.AddLabFrameFourVector(temp);
          lepID_BDT_ISR_Sparticle.push_back(key);
          objID_BDT_ISR_Sparticle.push_back(key);
        }
        else{
          RFKey key = VIS_ISR_Sparticle.AddLabFrameFourVector(leptons[i]);
          lepID_BDT_ISR_Sparticle.push_back(key);
          objID_BDT_ISR_Sparticle.push_back(key);
        }
      }
    }
    int Njets_S = 0;
    int Njets_ISR = 0;
    int Njets_S_matched = 0;
    int Njets_S_unmatched = 0;
    int Njets_ISR_matched = 0;
    int Njets_ISR_unmatched = 0;
    if(!LAB_ISR_Sparticle.AnalyzeEvent()) cout << "Problem with RF ISR Sparticle Analyze Event \n";
    //if(Cand){
      vPISR = S_ISR_Sparticle.GetFourVector(CM_ISR_Sparticle).Vect();
      PTISR = S_ISR_Sparticle.GetTransverseFourVector(CM_ISR_Sparticle).Vect().Mag();
      if(BDT_MET || S_MET){
        if(HEM)
          vPINV = (Ia_ISR_Sparticle.GetFourVector(CM_ISR_Sparticle)+Ib_ISR_Sparticle.GetFourVector(CM_ISR_Sparticle)).Vect();
        else
          vPINV = I_ISR_Sparticle.GetFourVector(CM_ISR_Sparticle).Vect();
      }
      else
        vPINV = MET;
      RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
      //if(KIN && (RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut)) continue;
        hist_RISR_PTISR_ISR_Sparticle->Fill(RISR, PTISR, hweight);

        std::vector<V_Cand> V_had_cands_S_ISR_Sparticle;
        std::vector<V_Cand> V_had_cands_Sa_ISR_Sparticle;
        std::vector<V_Cand> V_had_cands_Sb_ISR_Sparticle;
        if(HEM){
          if(S_leps){
            V_had_cands_Sa_ISR_Sparticle = cand_list(jets,leptons,COMB_J_ISR_Sparticle,BDT_Sa_ISR_Sparticle,J_Sa_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
            V_had_cands_Sb_ISR_Sparticle = cand_list(jets,leptons,COMB_J_ISR_Sparticle,BDT_Sb_ISR_Sparticle,J_Sb_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
          }
          else{
            V_had_cands_Sa_ISR_Sparticle = cand_list(jets,leptons,VIS_ISR_Sparticle,BDT_Sa_ISR_Sparticle,V_S_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
            V_had_cands_Sb_ISR_Sparticle = cand_list(jets,leptons,VIS_ISR_Sparticle,BDT_Sb_ISR_Sparticle,V_S_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
          }
        } // if(HEM)
        else{
          if(S_leps){
            V_had_cands_S_ISR_Sparticle = cand_list(jets,leptons,COMB_J_ISR_Sparticle,BDT_S_ISR_Sparticle,J_S_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
          }
          else{
            V_had_cands_S_ISR_Sparticle = cand_list(jets,leptons,VIS_ISR_Sparticle,BDT_S_ISR_Sparticle,V_S_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
          }
        } // not HEM

        std::vector<V_Cand> V_had_cands_ISR_ISR_Sparticle;
        if(S_leps)
          V_had_cands_ISR_ISR_Sparticle = cand_list(jets,leptons,COMB_J_ISR_Sparticle,BDT_ISR_ISR_Sparticle,V_ISR_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);
        else
          V_had_cands_ISR_ISR_Sparticle = cand_list(jets,leptons,VIS_ISR_Sparticle,BDT_ISR_ISR_Sparticle,V_ISR_ISR_Sparticle,CM_ISR_Sparticle,jetID_BDT_ISR_Sparticle,lepID_BDT_ISR_Sparticle);

        if(HEM){ // merge cand lists from Sa and Sb sides
          for(int i = 0; i < int(V_had_cands_Sa_ISR_Sparticle.size()); i++)
            V_had_cands_S_ISR_Sparticle.push_back(V_had_cands_Sa_ISR_Sparticle[i]);
          for(int i = 0; i < int(V_had_cands_Sb_ISR_Sparticle.size()); i++)
            V_had_cands_S_ISR_Sparticle.push_back(V_had_cands_Sb_ISR_Sparticle[i]);
        }

        for(int i = 0; i < Njets; i++){
          if(VIS_ISR_Sparticle.GetFrame(jetID_BDT_ISR_Sparticle[i]) == V_S_ISR_Sparticle || COMB_J_ISR_Sparticle.GetFrame(jetID_BDT_ISR_Sparticle[i]) == J_S_ISR_Sparticle
          || COMB_J_ISR_Sparticle.GetFrame(jetID_BDT_ISR_Sparticle[i]) == J_Sa_ISR_Sparticle || COMB_J_ISR_Sparticle.GetFrame(jetID_BDT_ISR_Sparticle[i]) == J_Sb_ISR_Sparticle){
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


        int N_V_had_Sa_ISR_Sparticle = V_had_cands_Sa_ISR_Sparticle.size();
        int N_V_had_Sb_ISR_Sparticle = V_had_cands_Sb_ISR_Sparticle.size();
        int N_V_had_ISR_Sparticle = V_had_cands_S_ISR_Sparticle.size()+V_had_cands_ISR_ISR_Sparticle.size();
        totcands_ISR_Sparticle += N_V_had_ISR_Sparticle;
        int N_V_had_S_ISR_Sparticle = V_had_cands_S_ISR_Sparticle.size();
        int N_V_had_ISR_ISR_Sparticle = V_had_cands_ISR_ISR_Sparticle.size();
        hist_N_V_had_ISR_Sparticle->Fill(N_V_had_ISR_Sparticle, hweight);
        hist_N_V_had_S_ISR_Sparticle->Fill(N_V_had_S_ISR_Sparticle, hweight);
        hist_N_V_had_ISR_ISR_Sparticle->Fill(N_V_had_ISR_ISR_Sparticle, hweight);
        hist_Njets_N_V_had_ISR_Sparticle->Fill(Njets, N_V_had_ISR_Sparticle, hweight);
        hist_Njets_N_V_had_S_ISR_Sparticle->Fill(Njets, N_V_had_S_ISR_Sparticle, hweight);
        hist_Njets_N_V_had_ISR_ISR_Sparticle->Fill(Njets, N_V_had_ISR_ISR_Sparticle, hweight);

        cand_side(V_had_cands_S_ISR_Sparticle, BDT_S_ISR_Sparticle);
        cand_side(V_had_cands_ISR_ISR_Sparticle, BDT_ISR_ISR_Sparticle);
        cand_matching(V_had_cands_S_ISR_Sparticle);
        cand_matching(V_had_cands_Sa_ISR_Sparticle);
        cand_matching(V_had_cands_Sb_ISR_Sparticle);
        cand_matching(V_had_cands_ISR_ISR_Sparticle);

        int matched_S_ISR_Sparticle = 0; // matched cands in the Sparticle system of the ISR_Sparticle tree
        int unmatched_S_ISR_Sparticle = 0; // unmatched cands in the Sparticle system of the ISR_Sparticle tree
        int matched_Sa_ISR_Sparticle = 0; // matched cands in the A side of Sparticle system of the ISR_Sparticle tree
        int unmatched_Sa_ISR_Sparticle = 0; // unmatched cands in A side of the Sparticle system of the ISR_Sparticle tree
        int matched_Sb_ISR_Sparticle = 0; // matched cands in the B side of Sparticle system of the ISR_Sparticle tree
        int unmatched_Sb_ISR_Sparticle = 0; // unmatched cands in B side of the Sparticle system of the ISR_Sparticle tree
        int matched_ISR_ISR_Sparticle = 0; // matched cands in the ISR system of the ISR_Sparticle tree
        int unmatched_ISR_ISR_Sparticle = 0; // unmatched cands in the ISR system of the ISR_Sparticle tree
        int partial_S_ISR_Sparticle = 0; // partial cands in the Sparticle system of the ISR_Sparticle tree
        int partial_ISR_ISR_Sparticle = 0; // partial cands in the ISR system of the ISR_Sparticle tree

        for(int i = 0; i < N_V_had_S_ISR_Sparticle; i++)
          if(V_had_cands_S_ISR_Sparticle[i].Match() == kMatched)
            matched_S_ISR_Sparticle++;
          else if(V_had_cands_S_ISR_Sparticle[i].Match() == kUnmatched)
            unmatched_S_ISR_Sparticle++;
          else
            partial_S_ISR_Sparticle++;
        for(int i = 0; i < N_V_had_ISR_ISR_Sparticle; i++)
          if(V_had_cands_ISR_ISR_Sparticle[i].Match() == kMatched)
            matched_ISR_ISR_Sparticle++;
          else if(V_had_cands_ISR_ISR_Sparticle[i].Match() == kUnmatched)
            unmatched_ISR_ISR_Sparticle++;
          else
            partial_ISR_ISR_Sparticle++;

        for(int i = 0; i < N_V_had_Sa_ISR_Sparticle; i++)
          if(V_had_cands_Sa_ISR_Sparticle[i].Match() == kMatched)
            matched_Sa_ISR_Sparticle++;
          else if(V_had_cands_S_ISR_Sparticle[i].Match() == kUnmatched)
            unmatched_Sa_ISR_Sparticle++;
        for(int i = 0; i < N_V_had_Sb_ISR_Sparticle; i++)
          if(V_had_cands_Sb_ISR_Sparticle[i].Match() == kMatched)
            matched_Sb_ISR_Sparticle++;
          else if(V_had_cands_S_ISR_Sparticle[i].Match() == kUnmatched)
            unmatched_Sb_ISR_Sparticle++;

        hist_matched_S_ISR_Sparticle->Fill(matched_S_ISR_Sparticle, hweight);
        hist_unmatched_S_ISR_Sparticle->Fill(unmatched_S_ISR_Sparticle, hweight);
        hist_matched_ISR_ISR_Sparticle->Fill(matched_ISR_ISR_Sparticle, hweight);
        hist_unmatched_ISR_ISR_Sparticle->Fill(unmatched_ISR_ISR_Sparticle, hweight);
        hist_partial_S_ISR_Sparticle->Fill(partial_S_ISR_Sparticle, hweight);
        hist_partial_ISR_ISR_Sparticle->Fill(partial_ISR_ISR_Sparticle, hweight);

        hist_Njets_S_matched_S_ISR_Sparticle->Fill(Njets_S, matched_S_ISR_Sparticle, hweight);
        hist_Njets_S_unmatched_S_ISR_Sparticle->Fill(Njets_S, unmatched_S_ISR_Sparticle, hweight);
        hist_Njets_S_partial_S_ISR_Sparticle->Fill(Njets_S, partial_S_ISR_Sparticle, hweight);
        hist_Njets_ISR_matched_ISR_ISR_Sparticle->Fill(Njets_ISR, matched_ISR_ISR_Sparticle, hweight);
        hist_Njets_ISR_unmatched_ISR_ISR_Sparticle->Fill(Njets_ISR, unmatched_ISR_ISR_Sparticle, hweight);
        hist_Njets_ISR_partial_ISR_ISR_Sparticle->Fill(Njets_ISR, partial_ISR_ISR_Sparticle, hweight);

        hist_matched_Njets_S_matched_S_ISR_Sparticle->Fill(Njets_S_matched, matched_S_ISR_Sparticle, hweight);
        hist_unmatched_Njets_S_unmatched_S_ISR_Sparticle->Fill(Njets_S_unmatched, unmatched_S_ISR_Sparticle, hweight);
        hist_matched_Njets_ISR_matched_ISR_ISR_Sparticle->Fill(Njets_ISR_matched, matched_ISR_ISR_Sparticle, hweight);
        hist_unmatched_Njets_ISR_unmatched_ISR_ISR_Sparticle->Fill(Njets_ISR_unmatched, unmatched_ISR_ISR_Sparticle, hweight);
        hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle->Fill(Njets_S_matched, Njets_S_unmatched, hweight);
        hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle->Fill(Njets_ISR_matched, Njets_ISR_unmatched, hweight);
        hist_matched_S_matched_Sa_ISR_Sparticle->Fill(matched_S_ISR_Sparticle, matched_Sa_ISR_Sparticle, hweight);
        hist_matched_S_matched_Sb_ISR_Sparticle->Fill(matched_S_ISR_Sparticle, matched_Sb_ISR_Sparticle, hweight);
    //} // if(Cand)

    // binary decay tree with S & ISR splitting #2
    LAB_ISR_Sparticle2.ClearEvent();
    vector<RFKey> jetID_BDT_ISR_Sparticle2;
    vector<RFKey> lepID_BDT_ISR_Sparticle2;
    vector<RFKey> objID_BDT_ISR_Sparticle2;
    INV_ISR_Sparticle2.SetLabFrameThreeVector(MET);

    // use reco jets for BDT
    if(!use_gen_jets){
      for(int i = 0; i < Njets; i++){
        if(Zero_RJR_BDT_JetMass){
          RFKey key = COMB_J_ISR_Sparticle2.AddLabFrameFourVector(TLorentzVector(jets[i].Pt(),jets[i].Eta(),jets[i].Phi(),0.));
          jetID_BDT_ISR_Sparticle2.push_back(key);
          objID_BDT_ISR_Sparticle2.push_back(key);
        }
        else if(transverseOBJ){
          Particle temp = jets[i];
          temp.SetZ(0.);
          RFKey key = COMB_J_ISR_Sparticle2.AddLabFrameFourVector(temp);
          jetID_BDT_ISR_Sparticle2.push_back(key);
          objID_BDT_ISR_Sparticle2.push_back(key);
        }
        else{
          RFKey key = COMB_J_ISR_Sparticle2.AddLabFrameFourVector(jets[i]);
          jetID_BDT_ISR_Sparticle2.push_back(key);
          objID_BDT_ISR_Sparticle2.push_back(key);
        }
      }
    }

    for(int i = 0; i < Nleps; i++){
      if(transverseOBJ){
        Particle temp = leptons[i];
        temp.SetZ(0.);
        RFKey key = COMB_L_ISR_Sparticle2.AddLabFrameFourVector(temp);
        lepID_BDT_ISR_Sparticle2.push_back(key);
        objID_BDT_ISR_Sparticle2.push_back(key);
      }
      else{
        RFKey key = COMB_L_ISR_Sparticle2.AddLabFrameFourVector(leptons[i]);
        lepID_BDT_ISR_Sparticle2.push_back(key);
        objID_BDT_ISR_Sparticle2.push_back(key);
      }
    }
    if(!LAB_ISR_Sparticle2.AnalyzeEvent()) cout << "Problem with RF ISR Sparticle Analyze Event #2 \n";
    vPISR = S_ISR_Sparticle2.GetFourVector(CM_ISR_Sparticle2).Vect();
    PTISR = S_ISR_Sparticle2.GetTransverseFourVector(CM_ISR_Sparticle2).Vect().Mag();
    vPINV = (Ia_ISR_Sparticle2.GetFourVector(CM_ISR_Sparticle2)+Ib_ISR_Sparticle2.GetFourVector(CM_ISR_Sparticle2)).Vect();
    RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
    hist_RISR_PTISR_ISR_Sparticle2->Fill(RISR, PTISR, hweight);
    if(RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut) continue;

    int Nleps_A = 0;
    int Nleps_B = 0;
    Njets_S = 0;
    Njets_ISR = 0;
    Njets_S_matched = 0;
    Njets_S_unmatched = 0;
    int Njets_Sa_matched = 0;
    int Njets_Sa_unmatched = 0;
    int Njets_Sb_matched = 0;
    int Njets_Sb_unmatched = 0;
    Njets_ISR_matched = 0;
    Njets_ISR_unmatched = 0;
    for(int i = 0; i < Njets; i++){
      if(COMB_J_ISR_Sparticle2.GetFrame(jetID_BDT_ISR_Sparticle2[i]) == J_X2a_ISR_Sparticle2){
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
      if(COMB_L_ISR_Sparticle2.GetFrame(lepID_BDT_ISR_Sparticle2[i]) == L_X2a_ISR_Sparticle2)
        Nleps_A++;
      else
        Nleps_B++;
    }
    hist_matched_Njets_S_unmatched_Njets_S_ISR_Sparticle2->Fill(Njets_S_matched, Njets_S_unmatched, hweight);
    hist_matched_Njets_ISR_unmatched_Njets_ISR_ISR_Sparticle2->Fill(Njets_ISR_matched, Njets_ISR_unmatched, hweight);
    hist_matched_Njets_Sa_matched_Njets_Sb_ISR_Sparticle2->Fill(Njets_Sa_matched, Njets_Sb_matched, hweight);
    hist_unmatched_Njets_Sa_unmatched_Njets_Sb_ISR_Sparticle2->Fill(Njets_Sa_unmatched, Njets_Sb_unmatched, hweight);
    if(Nleps >= 2){
      TLorentzVector l1l2 = leptons[0]+leptons[1];
      double mllLEAD = l1l2.M();
      hist_RISR_mllLEAD_ISR_Sparticle2->Fill(RISR, mllLEAD, hweight);
    }
    if(Nleps > 2) hist_RISR_mL_ISR_Sparticle2->Fill(RISR, (L_X2a_ISR_Sparticle2.GetFourVector() + L_X2b_ISR_Sparticle2.GetFourVector()).M(), hweight);
    for(int i = 0; i < Nleps; i++){
      for(int j = i+1; j < Nleps; j++){
        TLorentzVector l1l2 = leptons[i]+leptons[j];
        double mll = l1l2.M();
        hist_RISR_mll_ISR_Sparticle2->Fill(RISR,mll,hweight);
      }
    }
    if(Nleps_A >= 2){
      for(int i = 0; i < Nleps; i++){
        if(COMB_L_ISR_Sparticle2.GetFrame(lepID_BDT_ISR_Sparticle2[i]) != L_X2a_ISR_Sparticle2) continue;
        for(int j = i+1; j < Nleps; j++){
          if(COMB_L_ISR_Sparticle2.GetFrame(lepID_BDT_ISR_Sparticle2[j]) != L_X2a_ISR_Sparticle2) continue;
          TLorentzVector l1l2 = leptons[i]+leptons[j];
          double mll = l1l2.M();
          hist_RISR_mllHEM_ISR_Sparticle2->Fill(RISR,mll,hweight);
        }
      }
    }
    if(Nleps_B >= 2){
      for(int i = 0; i < Nleps; i++){
        if(COMB_L_ISR_Sparticle2.GetFrame(lepID_BDT_ISR_Sparticle2[i]) == L_X2a_ISR_Sparticle2) continue;
        for(int j = i+1; j < Nleps; j++){
          if(COMB_L_ISR_Sparticle2.GetFrame(lepID_BDT_ISR_Sparticle2[j]) == L_X2a_ISR_Sparticle2) continue;
          TLorentzVector l1l2 = leptons[i]+leptons[j];
          double mll = l1l2.M();
          hist_RISR_mllHEM_ISR_Sparticle2->Fill(RISR,mll,hweight);
        }
      }
    }
    if(Njets_S > 0){
      hist_RISR_mSJ_ISR_Sparticle2->Fill(RISR, J_X2a_ISR_Sparticle2.GetFourVector().M(), hweight);
      hist_CosThetaCM_mSJ_ISR_Sparticle2->Fill(CM_ISR_Sparticle2.GetCosDecayAngle(J_X2a_ISR_Sparticle2), J_X2a_ISR_Sparticle2.GetFourVector().M(), hweight);
    }
    hist_dphiCMI_PTCM_ISR_Sparticle2->Fill(CM_ISR_Sparticle2.GetDeltaPhiBoostVisible(), CM_ISR_Sparticle2.GetFourVector().Pt(), hweight);
    hist_dphiMETV_PTISR_ISR_Sparticle2->Fill(S_ISR_Sparticle2.GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(MET), PTISR, hweight);
    hist_Ma_Mb_ISR_Sparticle2->Fill(X2a_ISR_Sparticle2.GetFourVector().M(), X2b_ISR_Sparticle2.GetFourVector().M(), hweight);
    hist_MLa_MLb_ISR_Sparticle2->Fill(L_X2a_ISR_Sparticle2.GetFourVector().M(), L_X2b_ISR_Sparticle2.GetFourVector().M(), hweight);
    hist_RISR_DiffML_ISR_Sparticle2->Fill(RISR, L_X2a_ISR_Sparticle2.GetFourVector().M()-L_X2b_ISR_Sparticle2.GetFourVector().M(), hweight);
    double MVis = (J_X2a_ISR_Sparticle2.GetFourVector() + L_X2a_ISR_Sparticle2.GetFourVector() + L_X2b_ISR_Sparticle2.GetFourVector()).M();
    hist_RISR_MVis_ISR_Sparticle2->Fill(RISR, MVis, hweight);

    TLorentzVector vP_S_CM = S_ISR_Sparticle2.GetFourVector(CM_ISR_Sparticle2);
    TLorentzVector vP_Ja_S = J_X2a_ISR_Sparticle2.GetFourVector(S_ISR_Sparticle2);
    TLorentzVector vP_La_S = L_X2a_ISR_Sparticle2.GetFourVector(S_ISR_Sparticle2);
    TLorentzVector vP_Lb_S = L_X2b_ISR_Sparticle2.GetFourVector(S_ISR_Sparticle2);
    TLorentzVector vP_Ia_S = Ia_ISR_Sparticle2.GetFourVector(S_ISR_Sparticle2);
    TLorentzVector vP_Ib_S = Ib_ISR_Sparticle2.GetFourVector(S_ISR_Sparticle2);

    double MS_S = (vP_Ja_S+vP_La_S+vP_Lb_S+vP_Ia_S+vP_Ib_S).M();
    hist_RISR_MS_S_ISR_Sparticle2->Fill(RISR, MS_S, hweight);
    double MSa_S = (vP_Ja_S+vP_La_S+vP_Ia_S).M();
    double MSb_S = (vP_Lb_S+vP_Ib_S).M();
    hist_RISR_MSa_S_ISR_Sparticle2->Fill(RISR, MSa_S, hweight);
    hist_RISR_MSb_S_ISR_Sparticle2->Fill(RISR, MSb_S, hweight);
    hist_MSa_S_MSb_S_ISR_Sparticle2->Fill(MSa_S, MSb_S, hweight);

    double MVisS_S = (vP_Ja_S+vP_La_S+vP_Lb_S).M();
    hist_RISR_MVisS_S_ISR_Sparticle2->Fill(RISR, MVisS_S, hweight);
    double MVisSa_S = (vP_Ja_S+vP_La_S).M();
    hist_RISR_MVisSa_S_ISR_Sparticle2->Fill(RISR, MVisSa_S, hweight);
    double MVisSb_S = (vP_Lb_S).M();
    hist_RISR_MVisSb_S_ISR_Sparticle2->Fill(RISR, MVisSb_S, hweight);
    hist_MVisSa_S_MVisSb_S_ISR_Sparticle2->Fill(MVisSa_S, MVisSb_S, hweight);

    TVector3 boostVis = (vP_Ja_S+vP_La_S+vP_Lb_S).BoostVector();
    TVector3 boostInv = (vP_Ia_S+vP_Ib_S).BoostVector();
    TVector3 daBoost = vP_S_CM.Vect().Unit();
    boostVis = (boostVis.Dot(daBoost))*daBoost;
    boostInv = (boostInv.Dot(daBoost))*daBoost;
    if((!std::isnan(boostVis.Mag())) &&
     (boostVis.Mag() < 1)){
    vP_Ja_S.Boost(-boostVis);
    vP_La_S.Boost(-boostVis);
    vP_Lb_S.Boost(-boostVis);
    } else {
      vP_Ja_S.SetVectM(TVector3(0.,0.,0.),std::max(0.,vP_Ja_S.M()));
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
    double PX3_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).P();
    double MX3a_BoostT = (vP_Ja_S+vP_La_S+vP_Ia_S).M();
    double MX3b_BoostT = (vP_Lb_S+vP_Ib_S).M();
    double Mperp = sqrt(MX3a_BoostT*MX3a_BoostT+MX3b_BoostT*MX3b_BoostT)/sqrt(2.);
    double gammaPerp = 2*Mperp/(sqrt(MX3a_BoostT*MX3a_BoostT+PX3_BoostT*PX3_BoostT) +
		                sqrt(MX3b_BoostT*MX3b_BoostT+PX3_BoostT*PX3_BoostT));
    hist_Mperp_RISR_ISR_Sparticle2->Fill(Mperp, RISR, hweight);
    hist_gammaPerp_RISR_ISR_Sparticle2->Fill(gammaPerp, RISR, hweight);

    TVector3 vP_Ja_S0;
    TVector3 vP_La_S0;
    TVector3 vP_Lb_S0;
    TVector3 vP_Ia_S0;
    TVector3 vP_Ib_S0;
    vP_Ja_S0.SetPtEtaPhi(vP_Ja_S0.Pt(),vP_Ja_S0.Eta(),vP_Ja_S0.Phi());
    vP_La_S0.SetPtEtaPhi(vP_La_S0.Pt(),vP_La_S0.Eta(),vP_La_S0.Phi());
    vP_Lb_S0.SetPtEtaPhi(vP_Lb_S0.Pt(),vP_Lb_S0.Eta(),vP_Lb_S0.Phi());
    vP_Ia_S0.SetPtEtaPhi(vP_Ia_S0.Pt(),vP_Ia_S0.Eta(),vP_Ia_S0.Phi());
    vP_Ib_S0.SetPtEtaPhi(vP_Ib_S0.Pt(),vP_Ib_S0.Eta(),vP_Ib_S0.Phi());

    vP_Ja_S0 = vP_Ja_S0 - vP_Ja_S0.Dot(vP_S_CM.Vect())*daBoost;
    vP_La_S0 = vP_La_S0 - vP_La_S0.Dot(vP_S_CM.Vect())*daBoost;
    vP_Lb_S0 = vP_Lb_S0 - vP_Lb_S0.Dot(vP_S_CM.Vect())*daBoost;
    vP_Ia_S0 = vP_Ia_S0 - vP_Ia_S0.Dot(vP_S_CM.Vect())*daBoost;
    vP_Ib_S0 = vP_Ib_S0 - vP_Ib_S0.Dot(vP_S_CM.Vect())*daBoost;

    double MSperp = (vP_Ja_S0 + vP_La_S0 + vP_Lb_S0 + vP_Ia_S0 + vP_Ib_S0).Mag();
    hist_MSperp_RISR_ISR_Sparticle2->Fill(MSperp, RISR, hweight);

    // fill original RJR Tree but use candidates and other jets as separate inputs
    // BDT_ISR approach
    LAB.ClearEvent();
    vector<RFKey> candID;
    jetID.clear();
    lepID.clear();
    objID.clear();
    INV.SetLabFrameThreeVector(MET);
    for(int i = 0; i < N_V_had; i++){
      RFKey key = COMB_J.AddLabFrameFourVector(V_had_cands[i].CandFrame().GetFourVector());
      candID.push_back(key);
      objID.push_back(key);
    }
    for(int i = 0; i < int(jet_BDT_noncand.size()); i++){
      RFKey key = COMB_J.AddLabFrameFourVector(jets[jet_BDT_noncand[i]]);
      jetID.push_back(key);
      objID.push_back(key);
    }
    for(int i = 0; i < Nleps; i++){
      RFKey key = COMB_L.AddLabFrameFourVector(leptons[i]);
      lepID.push_back(key);
      objID.push_back(key);
    }
    if(!LAB.AnalyzeEvent()) cout << "Problem with RF Cand Analyze Event \n";
    //if(Cand){
      Njet_S = 0;
      for(int i = 0; i < Njets; i++)
        if(COMB_J.GetFrame(jetID[i]) == Ja || COMB_J.GetFrame(jetID[i]) == Jb)
          Njet_S++;
      vPISR = S.GetFourVector(CM).Vect();
      PTISR = S.GetTransverseFourVector(CM).Vect().Mag();
      vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();
      RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
      hist_RISR_PTISR_Cands->Fill(RISR, PTISR, hweight);
    //} // if(Cand)

      // Use BDT after running ISR tree
      LAB_ISR_BDT.ClearEvent();
      vector<RFKey> jetID_ISR_BDT;
      vector<RFKey> lepID_ISR_BDT;
      vector<RFKey> objID_ISR_BDT;
      INV_ISR_BDT.SetLabFrameThreeVector(MET);
      // use reco jets for RJR
      if(!use_gen_jets){
        for(int i = 0; i < Njets; i++){
          RFKey key = COMB_J_ISR_BDT.AddLabFrameFourVector(jets[i]);
          jetID_ISR_BDT.push_back(key);
          objID_ISR_BDT.push_back(key);
        }
      }
      // use gen jets for RJR
      if(use_gen_jets){
        for(int i = 0; i < Ngen_jets; i++){
          RFKey key = COMB_J_ISR_BDT.AddLabFrameFourVector(gen_jets[i]);
          jetID_ISR_BDT.push_back(key);
          objID_ISR_BDT.push_back(key);
        }
      }
      for(int i = 0; i < Nleps; i++){
        RFKey key = COMB_L_ISR_BDT.AddLabFrameFourVector(leptons[i]);
        lepID_ISR_BDT.push_back(key);
        objID_ISR_BDT.push_back(key);
      }
    int Njet_S_ISR_BDT = 0;
    if(!LAB_ISR_BDT.AnalyzeEvent()) cout << "Problem with RF ISR BDT Analyze Event \n";
    //if(Cand){
      vPISR = S_ISR_BDT.GetFourVector(CM_ISR_BDT).Vect();
      PTISR = S_ISR_BDT.GetTransverseFourVector(CM_ISR_BDT).Vect().Mag();
      vPINV = (Ia_ISR_BDT.GetFourVector(CM_ISR_BDT)+Ib_ISR_BDT.GetFourVector(CM_ISR_BDT)).Vect();
      RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
      //if(KIN && (RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut)) continue;
        hist_RISR_PTISR_ISR_BDT->Fill(RISR, PTISR, hweight);

        for(int i = 0; i < Njets; i++)
          if(COMB_J_ISR_BDT.GetFrame(jetID_ISR_BDT[i]) == Ja_ISR_BDT || COMB_J_ISR_BDT.GetFrame(jetID_ISR_BDT[i]) == Jb_ISR_BDT)
            Njet_S_ISR_BDT++;

        std::vector<V_Cand> V_had_cands_ISR_BDT_Sa;
        if(!use_gen_jets)
          V_had_cands_ISR_BDT_Sa = cand_list(jets,leptons,COMB_J_ISR_BDT,saVa_ISR_BDT,Ja_ISR_BDT,CM_ISR_BDT,jetID_ISR_BDT,lepID_ISR_BDT);
        else if(use_gen_jets)
          V_had_cands_ISR_BDT_Sa = cand_list(gen_jets,leptons,COMB_J_ISR_BDT,saVa_ISR_BDT,Ja_ISR_BDT,CM_ISR_BDT,jetID_ISR_BDT,lepID_ISR_BDT);
        std::vector<V_Cand> V_had_cands_ISR_BDT_Sb;
        if(!use_gen_jets)
          V_had_cands_ISR_BDT_Sb = cand_list(jets,leptons,COMB_J_ISR_BDT,saVb_ISR_BDT,Jb_ISR_BDT,CM_ISR_BDT,jetID_ISR_BDT,lepID_ISR_BDT);
        else if(use_gen_jets)
          V_had_cands_ISR_BDT_Sb = cand_list(gen_jets,leptons,COMB_J_ISR_BDT,saVb_ISR_BDT,Jb_ISR_BDT,CM_ISR_BDT,jetID_ISR_BDT,lepID_ISR_BDT);

        int N_V_had_cands_ISR_BDT_Sa = V_had_cands_ISR_BDT_Sa.size();
        int N_V_had_cands_ISR_BDT_Sb = V_had_cands_ISR_BDT_Sb.size();

        cand_matching(V_had_cands_ISR_BDT_Sa);
        cand_matching(V_had_cands_ISR_BDT_Sb);
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sa; i++)
          V_had_cands_ISR_BDT_Sa[i].SetSide(kAside);
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sb; i++)
          V_had_cands_ISR_BDT_Sb[i].SetSide(kBside);

        int N_V_had_ISR_BDT = N_V_had_cands_ISR_BDT_Sa+N_V_had_cands_ISR_BDT_Sb;
        totcands_ISR_BDT += N_V_had_ISR_BDT;
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sa; i++){
          if(V_had_cands_ISR_BDT_Sa[i].Type() == kSib)
            prong2_ISR_BDT++;
          else if(V_had_cands_ISR_BDT_Sa[i].Type() == kAunt)
            auntprong_ISR_BDT++;
          if(V_had_cands_ISR_BDT_Sa[i].Type() == kLep)
            lepprong_ISR_BDT++;
        }
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sb; i++){
          if(V_had_cands_ISR_BDT_Sb[i].Type() == kSib)
            prong2_ISR_BDT++;
          else if(V_had_cands_ISR_BDT_Sb[i].Type() == kAunt)
            auntprong_ISR_BDT++;
          if(V_had_cands_ISR_BDT_Sb[i].Type() == kLep)
            lepprong_ISR_BDT++;
        }

        hist_CandCount_ISR_BDT->SetBinContent(0,hist_CandCount_ISR_BDT->GetBinContent(0)+Nhadbosons);
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sa; i++){
          bool matched = false;
          if(V_had_cands_ISR_BDT_Sa[i].Match() == kMatched){
            matchedAside_ISR_BDT++;
            hist_CandCount_ISR_BDT->SetBinContent(1,hist_CandCount_ISR_BDT->GetBinContent(1)+1);
            matched = true;
          }
          else{
            if(V_had_cands_ISR_BDT_Sa[i].Match() == kW){
              hist_CandCount_ISR_BDT->SetBinContent(2,hist_CandCount_ISR_BDT->GetBinContent(2)+1);
            } // kW
            else if(V_had_cands_ISR_BDT_Sa[i].Match() == kZ){
              hist_CandCount_ISR_BDT->SetBinContent(3,hist_CandCount_ISR_BDT->GetBinContent(3)+1);
            } // kZ
            else if(V_had_cands_ISR_BDT_Sa[i].Match() == kB){
              hist_CandCount_ISR_BDT->SetBinContent(4,hist_CandCount_ISR_BDT->GetBinContent(4)+1);
            } // kB
            else{
              hist_CandCount_ISR_BDT->SetBinContent(5,hist_CandCount_ISR_BDT->GetBinContent(5)+1);
              unmatchedAside_ISR_BDT++;
            } // unmatched
          } // else
          eff_GenJetMatchedAll_V_cand_eta->Fill(matched, V_had_cands_ISR_BDT_Sa[i].Eta(), hweight);
          eff_GenJetMatchedAll_V_cand_mass->Fill(matched, V_had_cands_ISR_BDT_Sa[i].M(), hweight);
        } // for(int i = 0; i < N_V_had; i++)
        for(int i = 0; i < N_V_had_cands_ISR_BDT_Sb; i++){
          bool matched = false;
          if(V_had_cands_ISR_BDT_Sb[i].Match() == kMatched){
            matchedBside_ISR_BDT++;
            hist_CandCount_ISR_BDT->SetBinContent(1,hist_CandCount_ISR_BDT->GetBinContent(1)+1);
            matched = true;
          }
          else{
            if(V_had_cands_ISR_BDT_Sb[i].Match() == kW){
              hist_CandCount_ISR_BDT->SetBinContent(2,hist_CandCount_ISR_BDT->GetBinContent(2)+1);
            } // kW
            else if(V_had_cands_ISR_BDT_Sb[i].Match() == kZ){
              hist_CandCount_ISR_BDT->SetBinContent(3,hist_CandCount_ISR_BDT->GetBinContent(3)+1);
            } // kZ
            else if(V_had_cands_ISR_BDT_Sb[i].Match() == kB){
              hist_CandCount_ISR_BDT->SetBinContent(4,hist_CandCount_ISR_BDT->GetBinContent(4)+1);
            } // kB
            else{
              hist_CandCount_ISR_BDT->SetBinContent(5,hist_CandCount_ISR_BDT->GetBinContent(5)+1);
              unmatchedBside_ISR_BDT++;
            } // unmatched
          } // else
          eff_GenJetMatchedAll_V_cand_eta->Fill(matched, V_had_cands_ISR_BDT_Sb[i].Eta(), hweight);
          eff_GenJetMatchedAll_V_cand_mass->Fill(matched, V_had_cands_ISR_BDT_Sb[i].M(), hweight);
        } // for(int i = 0; i < N_V_had; i++)
    //} // if(Cand)

    // Use BDT after running ISR tree but this time don't add jets associated with leptons
    if(jet_BDT_lep.size() + jet_BDT_nonlep.size() != Njets) cout << "WTF: " << jet_BDT_lep.size() << " " << jet_BDT_nonlep.size() << endl;
    LAB_BDT_ISR.ClearEvent();
    vector<RFKey> jetID_BDT_ISR;
    vector<RFKey> lepID_BDT_ISR;
    vector<RFKey> objID_BDT_ISR;
    INV_BDT_ISR.SetLabFrameThreeVector(MET);
    // use reco jets for RJR
    if(!use_gen_jets){
      if(jet_BDT_lep.size() > 0){
        Particle Jet_ISR = jets[jet_BDT_lep[0]];
        for(int i = 1; i < int(jet_BDT_lep.size()); i++){
          Jet_ISR.Merge(jets[jet_BDT_lep[i]]);
        }
        ISR_BDT_ISR.SetLabFrameFourVector(Jet_ISR);
      }
      for(int i = 0; i < int(jet_BDT_nonlep.size()); i++){
        RFKey key = COMB_J_BDT_ISR.AddLabFrameFourVector(jets[jet_BDT_nonlep[i]]);
        jetID_BDT_ISR.push_back(key);
        objID_BDT_ISR.push_back(key);
      }
    }
    for(int i = 0; i < Nleps; i++){
      RFKey key = COMB_L_BDT_ISR.AddLabFrameFourVector(leptons[i]);
      lepID_BDT_ISR.push_back(key);
      objID_BDT_ISR.push_back(key);
    }
    if(!LAB_BDT_ISR.AnalyzeEvent()) cout << "Problem with RF ISR BDT Analyze Event \n";
    //if(Cand){
      vPISR = S_BDT_ISR.GetFourVector(CM_BDT_ISR).Vect();
      PTISR = S_BDT_ISR.GetTransverseFourVector(CM_BDT_ISR).Vect().Mag();
      vPINV = (Ia_BDT_ISR.GetFourVector(CM_BDT_ISR)+Ib_BDT_ISR.GetFourVector(CM_BDT_ISR)).Vect();
      RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
      //if(KIN && (RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut)) continue;
        Njet_S_ISR_BDT = 0;
        for(int i = 0; i < Njets; i++)
          if(COMB_J_ISR_BDT.GetFrame(jetID_ISR_BDT[i]) == Ja_ISR_BDT || COMB_J_ISR_BDT.GetFrame(jetID_ISR_BDT[i]) == Jb_ISR_BDT)
            Njet_S_ISR_BDT++;

        std::vector<V_Cand> V_had_cands_BDT_ISR_lep_Sa;
        if(!use_gen_jets)
          V_had_cands_BDT_ISR_lep_Sa = cand_list(jets,leptons,COMB_J_BDT_ISR,saVa_BDT_ISR,Ja_BDT_ISR,CM_BDT_ISR,jetID_BDT_ISR,lepID_BDT_ISR);
        else if(use_gen_jets)
          V_had_cands_BDT_ISR_lep_Sa = cand_list(gen_jets,leptons,COMB_J_BDT_ISR,saVa_BDT_ISR,Ja_BDT_ISR,CM_BDT_ISR,jetID_BDT_ISR,lepID_BDT_ISR);
        std::vector<V_Cand> V_had_cands_BDT_ISR_lep_Sb;
        if(!use_gen_jets)
          V_had_cands_BDT_ISR_lep_Sb = cand_list(jets,leptons,COMB_J_BDT_ISR,saVb_BDT_ISR,Jb_BDT_ISR,CM_BDT_ISR,jetID_BDT_ISR,lepID_BDT_ISR);
        else if(use_gen_jets)
          V_had_cands_BDT_ISR_lep_Sb = cand_list(gen_jets,leptons,COMB_J_BDT_ISR,saVb_BDT_ISR,Jb_BDT_ISR,CM_BDT_ISR,jetID_BDT_ISR,lepID_BDT_ISR);

        int N_V_had_cands_BDT_ISR_lep_Sa = V_had_cands_BDT_ISR_lep_Sa.size();
        int N_V_had_cands_BDT_ISR_lep_Sb = V_had_cands_BDT_ISR_lep_Sb.size();

        cand_matching(V_had_cands_BDT_ISR_lep_Sa);
        cand_matching(V_had_cands_BDT_ISR_lep_Sb);
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sa; i++)
          V_had_cands_BDT_ISR_lep_Sa[i].SetSide(kAside);
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sb; i++)
          V_had_cands_BDT_ISR_lep_Sb[i].SetSide(kBside);

        int N_V_had_BDT_ISR_lep = N_V_had_cands_BDT_ISR_lep_Sa+N_V_had_cands_BDT_ISR_lep_Sb;
        totcands_BDT_ISR_lep += N_V_had_BDT_ISR_lep;
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sa; i++){
          if(V_had_cands_BDT_ISR_lep_Sa[i].Type() == kSib)
            prong2_BDT_ISR_lep++;
          else if(V_had_cands_BDT_ISR_lep_Sa[i].Type() == kAunt)
            auntprong_BDT_ISR_lep++;
          if(V_had_cands_BDT_ISR_lep_Sa[i].Type() == kLep)
            lepprong_BDT_ISR_lep++;
        }
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sb; i++){
          if(V_had_cands_BDT_ISR_lep_Sb[i].Type() == kSib)
            prong2_BDT_ISR_lep++;
          else if(V_had_cands_BDT_ISR_lep_Sb[i].Type() == kAunt)
            auntprong_BDT_ISR_lep++;
          if(V_had_cands_BDT_ISR_lep_Sb[i].Type() == kLep)
            lepprong_BDT_ISR_lep++;
        }

        hist_CandCount_BDT_ISR_lep->SetBinContent(0,hist_CandCount_BDT_ISR_lep->GetBinContent(0)+Nhadbosons);
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sa; i++){
          bool matched = false;
          if(V_had_cands_BDT_ISR_lep_Sa[i].Match() == kMatched){
            matchedAside_BDT_ISR_lep++;
            hist_CandCount_BDT_ISR_lep->SetBinContent(1,hist_CandCount_BDT_ISR_lep->GetBinContent(1)+1);
            matched = true;
          }
          else{
            if(V_had_cands_BDT_ISR_lep_Sa[i].Match() == kW){
              hist_CandCount_BDT_ISR_lep->SetBinContent(2,hist_CandCount_BDT_ISR_lep->GetBinContent(2)+1);
            } // kW
            else if(V_had_cands_BDT_ISR_lep_Sa[i].Match() == kZ){
              hist_CandCount_BDT_ISR_lep->SetBinContent(3,hist_CandCount_BDT_ISR_lep->GetBinContent(3)+1);
            } // kZ
            else if(V_had_cands_BDT_ISR_lep_Sa[i].Match() == kB){
              hist_CandCount_BDT_ISR_lep->SetBinContent(4,hist_CandCount_BDT_ISR_lep->GetBinContent(4)+1);
            } // kB
            else{
              hist_CandCount_BDT_ISR_lep->SetBinContent(5,hist_CandCount_BDT_ISR_lep->GetBinContent(5)+1);
              unmatchedAside_BDT_ISR_lep++;
            } // unmatched
          } // else
          eff_GenJetMatchedAll_V_cand_eta->Fill(matched, V_had_cands_BDT_ISR_lep_Sa[i].Eta(), hweight);
          eff_GenJetMatchedAll_V_cand_mass->Fill(matched, V_had_cands_BDT_ISR_lep_Sa[i].M(), hweight);
        } // for(int i = 0; i < N_V_had; i++)
        for(int i = 0; i < N_V_had_cands_BDT_ISR_lep_Sb; i++){
          bool matched = false;
          if(V_had_cands_BDT_ISR_lep_Sb[i].Match() == kMatched){
            matchedBside_BDT_ISR_lep++;
            hist_CandCount_BDT_ISR_lep->SetBinContent(1,hist_CandCount_BDT_ISR_lep->GetBinContent(1)+1);
            matched = true;
          }
          else{
            if(V_had_cands_BDT_ISR_lep_Sb[i].Match() == kW){
              hist_CandCount_BDT_ISR_lep->SetBinContent(2,hist_CandCount_BDT_ISR_lep->GetBinContent(2)+1);
            } // kW
            else if(V_had_cands_BDT_ISR_lep_Sb[i].Match() == kZ){
              hist_CandCount_BDT_ISR_lep->SetBinContent(3,hist_CandCount_BDT_ISR_lep->GetBinContent(3)+1);
            } // kZ
            else if(V_had_cands_BDT_ISR_lep_Sb[i].Match() == kB){
              hist_CandCount_BDT_ISR_lep->SetBinContent(4,hist_CandCount_BDT_ISR_lep->GetBinContent(4)+1);
            } // kB
            else{
              hist_CandCount_BDT_ISR_lep->SetBinContent(5,hist_CandCount_BDT_ISR_lep->GetBinContent(5)+1);
              unmatchedBside_BDT_ISR_lep++;
            } // unmatched
          } // else
          eff_GenJetMatchedAll_V_cand_eta->Fill(matched, V_had_cands_BDT_ISR_lep_Sb[i].Eta(), hweight);
          eff_GenJetMatchedAll_V_cand_mass->Fill(matched, V_had_cands_BDT_ISR_lep_Sb[i].M(), hweight);
        } // for(int i = 0; i < N_V_had; i++)
    //} // if(Cand)
    // Use BDT after running ISR tree but this time force singlet jets into ISR system
    LAB_BDT_ISR.ClearEvent();
    jetID_BDT_ISR.clear();
    lepID_BDT_ISR.clear();
    objID_BDT_ISR.clear();
    INV_BDT_ISR.SetLabFrameThreeVector(MET);
    // use reco jets for RJR
    if(!use_gen_jets){
      if(jet_BDT_singlet.size() > 0){
        Particle Jet_ISR = jets[jet_BDT_singlet[0]];
        for(int i = 1; i < int(jet_BDT_singlet.size()); i++){
          Jet_ISR.Merge(jets[jet_BDT_singlet[i]]);
        }
        ISR_BDT_ISR.SetLabFrameFourVector(Jet_ISR);
      }
      for(int i = 0; i < int(jet_BDT_nonsinglet.size()); i++){
        RFKey key = COMB_J_BDT_ISR.AddLabFrameFourVector(jets[jet_BDT_nonsinglet[i]]);
        jetID_BDT_ISR.push_back(key);
        objID_BDT_ISR.push_back(key);
      }
    }
    for(int i = 0; i < Nleps; i++){
      RFKey key = COMB_L_BDT_ISR.AddLabFrameFourVector(leptons[i]);
      lepID_BDT_ISR.push_back(key);
      objID_BDT_ISR.push_back(key);
    }
    if(!LAB_BDT_ISR.AnalyzeEvent()) cout << "Problem with RF ISR BDT Analyze Event \n";
    //if(Cand){
      vPISR = S_BDT_ISR.GetFourVector(CM_BDT_ISR).Vect();
      PTISR = S_BDT_ISR.GetTransverseFourVector(CM_BDT_ISR).Vect().Mag();
      vPINV = (Ia_BDT_ISR.GetFourVector(CM_BDT_ISR)+Ib_BDT_ISR.GetFourVector(CM_BDT_ISR)).Vect();
      RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
      //if(KIN && (RISR < RISR_cut || PTISR < PTISR_cut || MET.Mag() < MET_cut)) continue;
        hist_RISR_PTISR_BDT_ISR_singlet->Fill(RISR, PTISR, hweight);

        std::vector<V_Cand> V_had_cands_BDT_ISR_singlet;
        hist_CandCount_BDT_ISR_singlet->SetBinContent(0,hist_CandCount_BDT_ISR_singlet->GetBinContent(0)+Nhadbosons);

        Njet_S_ISR_BDT = 0;
        for(int i = 0; i < Njets; i++)
          if(COMB_J_BDT_ISR.GetFrame(jetID_BDT_ISR[i]) == Ja_BDT_ISR || COMB_J_BDT_ISR.GetFrame(jetID_BDT_ISR[i]) == Jb_BDT_ISR)
            Njet_S_ISR_BDT++;
        
        if(!use_gen_jets){
          for(int i = 0; i < Njets; i++){
            ParticleList V_cand_Part;
            ConstRestFrameList V_cand_RF;
            if(COMB_J_BDT_ISR.GetFrame(jetID_BDT_ISR[i]) != Ja_BDT_ISR) continue;
            const RestFrame& frame = saVa_BDT_ISR.GetFrame(jetID_BDT_ISR[i]);
            const RestFrame& sibling = frame.GetSiblingFrame();
            const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
              for(int j = 0; j < Nleps; j++){
                if(sibling == saVa_BDT_ISR.GetFrame(lepID_BDT_ISR[j])){
                  for(int k = i+1; k < Njets; k++){
                    if(aunt == saVa_BDT_ISR.GetFrame(jetID_BDT_ISR[k])){
                      if(frame.GetProductionFrame().GetMass() - sibling.GetFourVector(frame.GetProductionFrame()).M() < 100.){
                        V_cand_Part.push_back(jets[i]);
                        V_cand_Part.push_back(jets[k]);
                        V_cand_RF.Add(frame);
                        V_cand_RF.Add(aunt);
                        V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                        cand.SetSide(kAside);
                        if(use_lep_prong){
                          //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                            V_had_cands_BDT_ISR_singlet.push_back(cand);
                            lepprong_BDT_ISR_singlet++;
                          //}
                        }
                      }
                    }
                  }
                }
              }
            for(int j = i+1; j < Njets; j++){
              if(sibling == saVa_BDT_ISR.GetFrame(jetID_BDT_ISR[j])){
                if(frame.GetProductionFrame().GetMass() < 100.){
                  V_cand_Part.clear();
                  V_cand_RF.Clear();
                  V_cand_Part.push_back(jets[i]);
                  V_cand_Part.push_back(jets[j]); // 2-prong "sibling"
                  V_cand_RF.Add(frame);
                  V_cand_RF.Add(sibling);
                  V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                  cand.SetSide(kAside);
                  if(use_prong2){
                    //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                      V_had_cands_BDT_ISR_singlet.push_back(cand);
                      prong2_BDT_ISR_singlet++;
                    //}
                  }
                }
              }
              if(aunt == saVa_BDT_ISR.GetFrame(jetID_BDT_ISR[j])){
                bool lep_sibling = false;
                if(int(lepID_BDT_ISR.size()) == Nleps){
                  for(int k = 0; k < Nleps; k++){
                    if(sibling == saVa_BDT_ISR.GetFrame(lepID_BDT_ISR[k])){
                      lep_sibling = true;
                    }
                  }
                }
                if(lep_sibling) continue;
                if((frame.GetMass() + aunt.GetMass()) < 100.){
                  V_cand_Part.clear();
                  V_cand_RF.Clear();
                  V_cand_Part.push_back(jets[i]);
                  V_cand_Part.push_back(jets[j]); // 2-prong "aunt"
                  V_cand_RF.Add(frame);
                  V_cand_RF.Add(aunt);
                  V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                  cand.SetSide(kAside);
                  if(use_aunt_prong){
                    //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                      V_had_cands_BDT_ISR_singlet.push_back(cand);
                      auntprong_BDT_ISR_singlet++;
                    //}
                  }
                }
              }
            }
          }
        }
        if(!use_gen_jets){
          for(int i = 0; i < Njets; i++){
            ParticleList V_cand_Part;
            ConstRestFrameList V_cand_RF;
            if(COMB_J_BDT_ISR.GetFrame(jetID_BDT_ISR[i]) != Jb_BDT_ISR) continue;
            const RestFrame& frame = saVb_BDT_ISR.GetFrame(jetID_BDT_ISR[i]);
            const RestFrame& sibling = frame.GetSiblingFrame();
            const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
              for(int j = 0; j < Nleps; j++){
                if(sibling == saVb_BDT_ISR.GetFrame(lepID_BDT_ISR[j])){
                  for(int k = i+1; k < Njets; k++){
                    if(aunt == saVb_BDT_ISR.GetFrame(jetID_BDT_ISR[k])){
                      if(frame.GetProductionFrame().GetMass() - sibling.GetFourVector(frame.GetProductionFrame()).M() < 100.){
                        V_cand_Part.push_back(jets[i]);
                        V_cand_Part.push_back(jets[k]);
                        V_cand_RF.Add(frame);
                        V_cand_RF.Add(aunt);
                        V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                        cand.SetSide(kBside);
                        if(use_lep_prong){
                          //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                            V_had_cands_BDT_ISR_singlet.push_back(cand);
                            lepprong_BDT_ISR_singlet++;
                          //}
                        }
                      }
                    }
                  }
                }
              }
            for(int j = i+1; j < Njets; j++){
              if(sibling == saVb_BDT_ISR.GetFrame(jetID_BDT_ISR[j])){
                if(frame.GetProductionFrame().GetMass() < 100.){
                  V_cand_Part.clear();
                  V_cand_RF.Clear();
                  V_cand_Part.push_back(jets[i]);
                  V_cand_Part.push_back(jets[j]); // 2-prong "sibling"
                  V_cand_RF.Add(frame);
                  V_cand_RF.Add(sibling);
                  V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                  cand.SetSide(kBside);
                  if(use_prong2){
                    //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                      V_had_cands_BDT_ISR_singlet.push_back(cand);
                      prong2_BDT_ISR_singlet++;
                    //}
                  }
                }
              }
              if(aunt == saVb_BDT_ISR.GetFrame(jetID_BDT_ISR[j])){
                bool lep_sibling = false;
                if(int(lepID_BDT_ISR.size()) == Nleps){
                  for(int k = 0; k < Nleps; k++){
                    if(sibling == saVb_BDT_ISR.GetFrame(lepID_BDT_ISR[k])){
                      lep_sibling = true;
                    }
                  }
                }
                if(lep_sibling) continue;
                if((frame.GetMass() + aunt.GetMass()) < 100.){
                  V_cand_Part.clear();
                  V_cand_RF.Clear();
                  V_cand_Part.push_back(jets[i]);
                  V_cand_Part.push_back(jets[j]); // 2-prong "aunt"
                  V_cand_RF.Add(frame);
                  V_cand_RF.Add(aunt);
                  V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                  cand.SetSide(kBside);
                  if(use_aunt_prong){
                    //if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                      V_had_cands_BDT_ISR_singlet.push_back(cand);
                      auntprong_BDT_ISR_singlet++;
                    //}
                  }
                }
              }
            }
          }
        }

        int N_V_had_BDT_ISR_singlet = V_had_cands_BDT_ISR_singlet.size();
        totcands_BDT_ISR_singlet += N_V_had_BDT_ISR_singlet;
        for(int i = 0; i < N_V_had_BDT_ISR_singlet; i++){
          bool matched = false; // perfect match
          bool unmatched = true; // both jets are radiative
          bool found24 = false;
          bool found23 = false;
          bool found6 = false;
          for(int j = 0; j < int(V_had_cands_BDT_ISR_singlet[i].size()); j++){
            for(int k = j+1; k < int(V_had_cands_BDT_ISR_singlet[i].size()); k++){
              //if(abs(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID()) == 24 && V_had_cands_BDT_ISR_singlet[i][j].GenMomIndex() == V_had_cands_BDT_ISR_singlet[i][k].GenMomIndex()){
              if((abs(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID()) == 23 || abs(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID()) == 24) && V_had_cands_BDT_ISR_singlet[i][j].GenMomIndex() == V_had_cands_BDT_ISR_singlet[i][k].GenMomIndex()){
                V_had_cands_BDT_ISR_singlet[i].SetMatch(kMatched);
                hist_CandCount_BDT_ISR_singlet->SetBinContent(1,hist_CandCount_BDT_ISR_singlet->GetBinContent(1)+1);
                matched = true;
                unmatched = false;
                break;
              }
            }
            if(!matched){
              if(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID() == 24){
                if(found24)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(2,hist_CandCount_BDT_ISR_singlet->GetBinContent(2)-1);
                if(found23)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(3,hist_CandCount_BDT_ISR_singlet->GetBinContent(3)-1);
                if(found6)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(4,hist_CandCount_BDT_ISR_singlet->GetBinContent(4)-1);
                V_had_cands_BDT_ISR_singlet[i].SetMatch(kW);
                hist_CandCount_BDT_ISR_singlet->SetBinContent(2,hist_CandCount_BDT_ISR_singlet->GetBinContent(2)+1);
                found24 = true;
                unmatched = false;
              }
              else if(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID() == 23){
                if(found24)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(2,hist_CandCount_BDT_ISR_singlet->GetBinContent(2)-1);
                if(found23)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(3,hist_CandCount_BDT_ISR_singlet->GetBinContent(3)-1);
                if(found6)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(4,hist_CandCount_BDT_ISR_singlet->GetBinContent(4)-1);
                V_had_cands_BDT_ISR_singlet[i].SetMatch(kZ);
                hist_CandCount_BDT_ISR_singlet->SetBinContent(3,hist_CandCount_BDT_ISR_singlet->GetBinContent(3)+1);
                found23 = true;
                unmatched = false;
              }
              else if(V_had_cands_BDT_ISR_singlet[i][j].MomPDGID() == 6){
                if(found24)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(2,hist_CandCount_BDT_ISR_singlet->GetBinContent(2)-1);
                if(found23)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(3,hist_CandCount_BDT_ISR_singlet->GetBinContent(3)-1);
                if(found6)
                  hist_CandCount_BDT_ISR_singlet->SetBinContent(4,hist_CandCount_BDT_ISR_singlet->GetBinContent(4)-1);
                V_had_cands_BDT_ISR_singlet[i].SetMatch(kB);
                hist_CandCount_BDT_ISR_singlet->SetBinContent(4,hist_CandCount_BDT_ISR_singlet->GetBinContent(4)+1);
                found6 = true;
                unmatched = false;
              }
            } // if(!matched)
          } // for(int j = 0; j < int(V_had_cands_BDT_ISR_singlet[i].size()); j++)
          if(unmatched){
            V_had_cands_BDT_ISR_singlet[i].SetMatch(kUnmatched);
            hist_CandCount_BDT_ISR_singlet->SetBinContent(5,hist_CandCount_BDT_ISR_singlet->GetBinContent(5)+1);
          }
        } // for(int i = 0; i < N_V_had; i++)

        //const RestFrame& Prod_Frame_BDT_ISR_singlet = saVa_BDT_ISR_singlet;
        for(int i = 0; i < N_V_had_BDT_ISR_singlet; i++){
          //if(V_had_cands_BDT_ISR_singlet[i].Side() == kBside)
          //Prod_Frame_BDT_ISR_singlet = saVb_BDT_ISR_singlet;
          if(V_had_cands_BDT_ISR_singlet[i].Match() == kMatched){
            if(V_had_cands_BDT_ISR_singlet[i].Side() == kAside)
              matchedAside_BDT_ISR_singlet++;
            else
              matchedBside_BDT_ISR_singlet++;
          }
          else if(V_had_cands_BDT_ISR_singlet[i].Match() == kUnmatched){
            if(V_had_cands_BDT_ISR_singlet[i].Side() == kAside)
              unmatchedAside_BDT_ISR_singlet++;
            else
              unmatchedBside_BDT_ISR_singlet++;
          }
        }
    //} // if(Cand)


    kept_file_events++;
    hist_MET->Fill(base->MET_pt, hweight);
  }
  
  hist_Count->SetBinContent(0,processed_events);
  hist_Count->SetBinContent(1,kept_file_events);
  hist_Count->SetBinContent(2,quark_cut); // quark eff
  hist_Count->SetBinContent(3,tot_quark); // quark eff
  hist_Count->SetBinContent(4,genHBosons);
  hist_Count->SetBinContent(5,totcands);
  hist_Count->SetBinContent(6,prong2);
  hist_Count->SetBinContent(7,auntprong);
  hist_Count->SetBinContent(8,lepprong);
  hist_Count->SetBinContent(9,hist_CandCount->GetBinContent(1)); // matched
  hist_Count->SetBinContent(10,hist_CandCount->GetBinContent(2)); // W
  hist_Count->SetBinContent(11,hist_CandCount->GetBinContent(3)); // Z
  hist_Count->SetBinContent(12,hist_CandCount->GetBinContent(4)); // B
  hist_Count->SetBinContent(13,hist_CandCount->GetBinContent(5)); // unmatched
  hist_Count->SetBinContent(14,matchedAside);
  hist_Count->SetBinContent(15,matchedBside);
  hist_Count->SetBinContent(16,unmatchedAside);
  hist_Count->SetBinContent(17,unmatchedBside);

  hist_Count_ISR_BDT->SetBinContent(0,processed_events);
  hist_Count_ISR_BDT->SetBinContent(1,kept_file_events);
  hist_Count_ISR_BDT->SetBinContent(2,quark_cut); // quark eff
  hist_Count_ISR_BDT->SetBinContent(3,tot_quark); // quark eff
  hist_Count_ISR_BDT->SetBinContent(4,genHBosons);
  hist_Count_ISR_BDT->SetBinContent(5,totcands_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(6,prong2_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(7,auntprong_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(8,lepprong_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(9,hist_CandCount_ISR_BDT->GetBinContent(1)); // matched
  hist_Count_ISR_BDT->SetBinContent(10,hist_CandCount_ISR_BDT->GetBinContent(2)); // W
  hist_Count_ISR_BDT->SetBinContent(11,hist_CandCount_ISR_BDT->GetBinContent(3)); // Z
  hist_Count_ISR_BDT->SetBinContent(12,hist_CandCount_ISR_BDT->GetBinContent(4)); // B
  hist_Count_ISR_BDT->SetBinContent(13,hist_CandCount_ISR_BDT->GetBinContent(5)); // unmatched
  hist_Count_ISR_BDT->SetBinContent(14,matchedAside_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(15,matchedBside_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(16,unmatchedAside_ISR_BDT);
  hist_Count_ISR_BDT->SetBinContent(17,unmatchedBside_ISR_BDT);

  hist_Count_BDT_ISR_lep->SetBinContent(0,processed_events);
  hist_Count_BDT_ISR_lep->SetBinContent(1,kept_file_events);
  hist_Count_BDT_ISR_lep->SetBinContent(2,quark_cut); // quark eff
  hist_Count_BDT_ISR_lep->SetBinContent(3,tot_quark); // quark eff
  hist_Count_BDT_ISR_lep->SetBinContent(4,genHBosons);
  hist_Count_BDT_ISR_lep->SetBinContent(5,totcands_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(6,prong2_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(7,auntprong_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(8,lepprong_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(9,hist_CandCount_BDT_ISR_lep->GetBinContent(1)); // matched
  hist_Count_BDT_ISR_lep->SetBinContent(10,hist_CandCount_BDT_ISR_lep->GetBinContent(2)); // W
  hist_Count_BDT_ISR_lep->SetBinContent(11,hist_CandCount_BDT_ISR_lep->GetBinContent(3)); // Z
  hist_Count_BDT_ISR_lep->SetBinContent(12,hist_CandCount_BDT_ISR_lep->GetBinContent(4)); // B
  hist_Count_BDT_ISR_lep->SetBinContent(13,hist_CandCount_BDT_ISR_lep->GetBinContent(5)); // unmatched
  hist_Count_BDT_ISR_lep->SetBinContent(14,matchedAside_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(15,matchedBside_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(16,unmatchedAside_BDT_ISR_lep);
  hist_Count_BDT_ISR_lep->SetBinContent(17,unmatchedBside_BDT_ISR_lep);

  hist_Count_BDT_ISR_singlet->SetBinContent(0,processed_events);
  hist_Count_BDT_ISR_singlet->SetBinContent(1,kept_file_events);
  hist_Count_BDT_ISR_singlet->SetBinContent(2,quark_cut); // quark eff
  hist_Count_BDT_ISR_singlet->SetBinContent(3,tot_quark); // quark eff
  hist_Count_BDT_ISR_singlet->SetBinContent(4,genHBosons);
  hist_Count_BDT_ISR_singlet->SetBinContent(5,totcands_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(6,prong2_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(7,auntprong_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(8,lepprong_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(9,hist_CandCount_BDT_ISR_singlet->GetBinContent(1)); // matched
  hist_Count_BDT_ISR_singlet->SetBinContent(10,hist_CandCount_BDT_ISR_singlet->GetBinContent(2)); // W
  hist_Count_BDT_ISR_singlet->SetBinContent(11,hist_CandCount_BDT_ISR_singlet->GetBinContent(3)); // Z
  hist_Count_BDT_ISR_singlet->SetBinContent(12,hist_CandCount_BDT_ISR_singlet->GetBinContent(4)); // B
  hist_Count_BDT_ISR_singlet->SetBinContent(13,hist_CandCount_BDT_ISR_singlet->GetBinContent(5)); // unmatched
  hist_Count_BDT_ISR_singlet->SetBinContent(14,matchedAside_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(15,matchedBside_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(16,unmatchedAside_BDT_ISR_singlet);
  hist_Count_BDT_ISR_singlet->SetBinContent(17,unmatchedBside_BDT_ISR_singlet);
  
  string output_root_filename = "output_Plot_1D_NANO_"+ofile;
  TFile* output_file = new TFile((output_root_filename).c_str(),"RECREATE");
  output_file->mkdir(plot_folder.c_str());
  output_file->cd(plot_folder.c_str());
  for(int hist1 = 0; hist1 < int(hists1.size()); hist1++) hists1[hist1]->Write();
  for(int hist2 = 0; hist2 < int(hists2.size()); hist2++) hists2[hist2]->Write();
  for(int eff = 0; eff < int(effs.size()); eff++) effs[eff]->Write();
  hist_Count->Write();
  if(Cand) hist_CandCount->Write();
  if(Cand) hist_CandCount_Combo->Write();
  if(Cand) hist_CandCount_Chi2->Write();
  if(Cand) hist_Count_ISR_BDT->Write();
  if(Cand) hist_CandCount_ISR_BDT->Write();
  if(Cand) hist_Count_BDT_ISR_lep->Write();
  if(Cand) hist_CandCount_BDT_ISR_lep->Write();
  if(Cand) hist_Count_BDT_ISR_singlet->Write();
  if(Cand) hist_CandCount_BDT_ISR_singlet->Write();
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

//-/-/-/-/-/-/-/-/-//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
// chi2 min bs
double getX2(const double& pMass, const double& mass, const double& width){
  return ((pMass - mass)*(pMass - mass))/(width*width);
}

void getJetsMinX2(int& j1_index, int& j2_index, ParticleList& jets){
return; // speed up code by skipping this
  int NJ = jets.size();
  if(NJ < 3) return; // not enough jets
  double mW = 82.5; // W mass
  double wW = 12.6; // W width
  double mT = 172.8; // top mass
  double wT = 18.9; // top width
  TLorentzVector j0 = jets[0];
  TLorentzVector j1 = jets[1];
  TLorentzVector j2 = jets[2];
  double mTJJ = (j0 + j1 + j2).M();
  double X2T = getX2(mT,mTJJ,wT);
  double minX2T = X2T;
  int top_index[3] = {0, 1, 2};
  for(int i = 0; i < NJ; i++){
    for(int j = i+1; j < NJ; j++){
      for(int k = j+1; k < NJ; k++){
        j0 = jets[i];
	j1 = jets[j];
	j2 = jets[k];
	X2T = getX2(mT,(j0 + j1 + j2).M(),wT);
	if(X2T < minX2T){
	  top_index[0] = i;
	  top_index[1] = j;
	  top_index[2] = k;
	}
      }
    }
  }
  j0 = jets[top_index[0]];
  j1 = jets[top_index[1]];
  j2 = jets[top_index[2]];
  double X2W01 = getX2(mW,(j0 + j1).M(),wW);
  double X2W02 = getX2(mW,(j0 + j2).M(),wW);
  double X2W12 = getX2(mW,(j1 + j2).M(),wW);
  if(X2W01 < X2W02 && X2W01 < X2W12){
    j1_index = top_index[0];
    j2_index = top_index[1];
  }
  else if(X2W02 < X2W01 && X2W02 < X2W12){
    j1_index = top_index[0];
    j2_index = top_index[2];
  }
  else if(X2W12 < X2W01 && X2W12 < X2W02){
    j1_index = top_index[1];
    j2_index = top_index[2];
  }
  else{
    j1_index = top_index[0];
    j2_index = top_index[1];
  }
}
//-/-/-/-/-/-/-/-/-//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-

int FrameDistance(const RestFrame& frame, const RestFrame& sibling){
  int Distance = -2;
  const RestFrame& CM = frame.GetLabFrame().GetChildFrame(0);

  int depth = CM.GetFrameDepth(frame);
  int siblingDepth = CM.GetFrameDepth(sibling);
  int minDepth = depth;
  if(siblingDepth < depth){
    minDepth = siblingDepth;
  }
  minDepth--; // depth of Least Common Ancestor
  bool useCM = false;
  const RestFrame& LCA = CM.GetFrameAtDepth(minDepth,frame); // Least Common Ancestor
  if(!LCA.GetListFrames().Contains(frame) || !LCA.GetListFrames().Contains(sibling))
    useCM = true;
  if(!useCM){
    depth = LCA.GetFrameDepth(frame);
    siblingDepth = LCA.GetFrameDepth(sibling);
  }
  Distance += depth;
  Distance += siblingDepth;
  return Distance;
}

void load_files_from_list(vector<string>& files, string fileListName)
{
 std::ifstream infile(fileListName);
 string line = "";
 while(getline(infile,line))
 {
  files.push_back(line);
 }
}

string double_to_string(double val){
  std::stringstream strstr_val;
  strstr_val << std::fixed << std::setprecision(3) << val;
  string str_val = strstr_val.str();
  std::replace(str_val.begin(), str_val.end(), '.', 'p');
  return str_val;
}

std::vector<V_Cand> cand_list(const ParticleList& jets, const ParticleList& leps,
  const CombinatoricGroup& Comb, const SelfAssemblingRecoFrame& BDT, const VisibleRecoFrame& VIS, const RestFrame& CM,
  const std::vector<RFKey>& jetID, const std::vector<RFKey>& lepID)
{
  std::vector<V_Cand> V_had_cands;
  int Njets = jets.size();
  int Nleps = leps.size();
  for(int i = 0; i < Njets; i++){
    if(Comb.GetFrame(jetID[i]) != VIS) continue;
    ParticleList V_cand_Part;
    ConstRestFrameList V_cand_RF;
    const RestFrame& frame = BDT.GetFrame(jetID[i]);
    const RestFrame& sibling = frame.GetSiblingFrame();
    const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
    if(int(lepID.size()) == Nleps){
      for(int j = 0; j < Nleps; j++){
        if(sibling == BDT.GetFrame(lepID[j])){
          for(int k = i+1; k < Njets; k++){
            if(aunt == BDT.GetFrame(jetID[k])){
              if(!BDT.GetListFrames().Contains(aunt)) continue;
              if(frame.GetProductionFrame().GetMass() - sibling.GetFourVector(frame.GetProductionFrame()).M() < 100.){
                V_cand_Part.push_back(jets[i]);
                V_cand_Part.push_back(jets[k]);
                V_cand_RF.Add(frame);
                V_cand_RF.Add(aunt);
                V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                if(use_lep_prong){
                  if(CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                    cand.SetType(kLep);
                    V_had_cands.push_back(cand);
                  } // if(cos_cut)
                } // if(use_lep_prong)
              } // if(mass_cut)
            } // if(aunt == BDT.GetFrame(jetID[k]))
          } // for(int k = i+1; k < Njets; k++)
        } // if(sibling == BDT.GetFrame(lepID[j]))
      } // for(int j = 0; j < Nleps; j++)
    } // if(int(lepID.size()) == Nleps)
    for(int j = i+1; j < Njets; j++){
      if(sibling == BDT.GetFrame(jetID[j])){
        if(!BDT.GetListFrames().Contains(sibling)) continue;
        if(frame.GetProductionFrame().GetMass() < 100.){
          V_cand_Part.clear();
          V_cand_RF.Clear();
          V_cand_Part.push_back(jets[i]);
          V_cand_Part.push_back(jets[j]); // 2-prong "sibling"
          V_cand_RF.Add(frame);
          V_cand_RF.Add(sibling);
          V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
          if(use_prong2){
            if(CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
              cand.SetType(kSib);
              V_had_cands.push_back(cand);
            } // if(cos_cut)
          } // if(use_prong2)
        } // if(mass_cut)
      } // if(sibling == BDT.GetFrame(jetID[j]))
      if(aunt == BDT.GetFrame(jetID[j])){
        if(!BDT.GetListFrames().Contains(aunt)) continue;
        bool lep_sibling = false;
        if(int(lepID.size()) == Nleps){
          for(int k = 0; k < Nleps; k++){
            if(sibling == BDT.GetFrame(lepID[k])){
              lep_sibling = true;
            }
          }
        }
        if(lep_sibling) continue;
        if((frame.GetMass() + aunt.GetMass()) < 100.){
          V_cand_Part.clear();
          V_cand_RF.Clear();
          V_cand_Part.push_back(jets[i]);
          V_cand_Part.push_back(jets[j]); // 2-prong "aunt"
          V_cand_RF.Add(frame);
          V_cand_RF.Add(aunt);
          V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
          if(use_aunt_prong){
            if(CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
              cand.SetType(kAunt);
              V_had_cands.push_back(cand);
            } // if(cos_cut)
          } // if(use_aunt_prong)
        } // if(mass_cut)
      } // if(aunt == BDT.GetFrame(jetID[j]))
    } // for(int j = i+1; j < Njets; j++)
  } // for(int i = 0; i < Njets; i++)
  return V_had_cands;
} // std::vector<V_Cand> cand_list()

void cand_matching(std::vector<V_Cand>& cand_list){ 
  int N_cands = cand_list.size();
    for(int i = 0; i < N_cands; i++){
      bool unmatched = true; // both jets are radiative
      bool matched = false; // both jets come from same boson
      for(int j = 0; j < int(cand_list[i].size()); j++){
        for(int k = j+1; k < int(cand_list[i].size()); k++){
          if((abs(cand_list[i][j].MomPDGID()) == 23 || abs(cand_list[i][j].MomPDGID()) == 24) && cand_list[i][j].GenMomIndex() == cand_list[i][k].GenMomIndex()){
            cand_list[i].SetMatch(kMatched);
            unmatched = false;
            matched = true;
            break;
          }
        }
        if(!matched){
          if(cand_list[i][j].MomPDGID() == 24){
            cand_list[i].SetMatch(kW);
            unmatched = false;
          }
          else if(cand_list[i][j].MomPDGID() == 23){
            cand_list[i].SetMatch(kZ);
            unmatched = false;
          }
          else if(cand_list[i][j].MomPDGID() == 6){
            cand_list[i].SetMatch(kB);
            unmatched = false;
          }
        } // if(!matched)
      } // for(int j = 0; j < int(cand_list[i].size()); j++)
      if(unmatched){
        cand_list[i].SetMatch(kUnmatched);
      }
    } // for(int i = 0; i < N_V_had; i++)
} // cand_matching(const std::vector<V_Cand>& cand_list)

void cand_side(std::vector<V_Cand>& cand_list, const SelfAssemblingRecoFrame& BDT){
  for(int i = 0; i < int(cand_list.size()); i++){
      const RestFrame& BDT_ChildA = BDT.GetChildFrame(0);
      const RestFrame& BDT_ChildB = BDT.GetChildFrame(1);
      ConstRestFrameList BDT_ChildA_Children = BDT_ChildA.GetListFrames();
      ConstRestFrameList BDT_ChildB_Children = BDT_ChildB.GetListFrames();
      for(int r = 0; r < BDT_ChildA_Children.GetN(); r++)
        if(cand_list[i].CandFrame() == BDT_ChildA_Children[r])
          cand_list[i].SetSide(kAside);
      for(int r = 0; r < BDT_ChildB_Children.GetN(); r++)
        if(cand_list[i].CandFrame() == BDT_ChildB_Children[r])
          cand_list[i].SetSide(kBside);
  }

}
