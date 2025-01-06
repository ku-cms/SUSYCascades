// things to check
// check efficiency of truth matched candidates/generated W bosons
// check where events in 'horns' of cosdecayangle plot are ending up in other distributions
// consider other 'production' frames like the child frame of the CM

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
#include <TSystem.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TColorWheel.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TError.h>
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include "../include/SUSYNANOBase.hh"
#include "../include/NANORun3.hh"
#include "../include/SampleTool.hh"
#include "../include/CategoryTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/V_Cand.hh"
#include "../include/XsecTool.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
void Plot_Hist(TH1D* h, bool Scale=true);
void Plot_Hist(TH2D* h);
void Plot_Eff(TEfficiency* e);

string g_PlotTitle;
string g_Label;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;

bool make_plots = true;
bool use_gen_jets = false;
bool Zero_RJR_BDT_JetMass = false; // set mass of jets to zero before adding into RJR tree
bool use_prong2 = true;
bool use_aunt_prong = false;
bool use_lep_prong = false;
bool boson_acceptance_cut = true;
bool gen_lepton_cut = false;
// boson_acceptance_cut && gen_lepton_cut would be to target events like WV -> l nu, q qbar
int num_files = 3; // number of nano input root files to run over
double min_Cand_CosDecayAngle_CM = -1.; // default: -1
double max_Cand_CosDecayAngle_CM = 1.; // default: 1.

string plot_folder = "";
string output_root_file = "output.root";

//-/-/-/-/-/-/-/-/-//-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
// chi2 min bs
double getX2(const double& pMass, const double& mass, const double& width){
  return ((pMass - mass)*(pMass - mass))/(width*width);
}

void getJetsMinX2(int& j1_index, int& j2_index, ParticleList& jets){
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

using namespace RestFrames;

std::string get_str_between_two_str(const std::string &s, const std::string &start_delim, const std::string &stop_delim)
{
 unsigned first_delim_pos = s.find(start_delim);
 unsigned end_pos_of_first_delim = first_delim_pos + start_delim.length();
 unsigned last_delim_pos = s.find_first_of(stop_delim, end_pos_of_first_delim);
 return s.substr(end_pos_of_first_delim,last_delim_pos - end_pos_of_first_delim);
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

template <typename V>
bool inVec(const std::vector<V>& vect, const V& value){
  return std::find(vect.begin(), vect.end(), value) != vect.end();
}

double getWeightPerEvent(string input_dataset, string input_filetag, string EventCountFile);

string double_to_string(double val);

void Plot_NANO_Cands(){
  plot_folder += "plots";
  if(use_gen_jets) plot_folder += "_genjets";
  else plot_folder += "_recojets";
  if(Zero_RJR_BDT_JetMass) plot_folder += "_ZeroJetMass";
  if(use_prong2) plot_folder += "_prong2";
  if(use_aunt_prong) plot_folder += "_auntprong";
  if(use_lep_prong) plot_folder += "_lepprong";
  if(boson_acceptance_cut) plot_folder += "_qqbar";
  else plot_folder += "_XX";
  if(gen_lepton_cut) plot_folder += "_lnu";
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
  plot_folder += "_numfiles"+std::to_string(num_files);  
  plot_folder += "/"; // need '/' at end
  string num_info_filename = plot_folder+"num_info.txt";
  string file_event_info_filename = plot_folder+"file_event_info.txt";
  //if(!gSystem->AccessPathName(num_info_filename.c_str()))
  //  gSystem->Exec(("rm "+num_info_filename).c_str());
  //if(!gSystem->AccessPathName(file_event_info_filename.c_str()))
  //  gSystem->Exec(("rm "+file_event_info_filename).c_str());
  output_root_file = plot_folder+"output_Plot_NANO_Cands.root";
  gSystem->Exec(("mkdir -p "+plot_folder).c_str());
  std::cout << "Outputting to: " << plot_folder << std::endl;

  RestFrames::SetStyle();
  g_Xmin = 0.;
  g_Xmax = 500.; 
  g_NX = 128;
  g_PlotTitle = "";
  string path_to_lists = "/home/zflowers/CMSSW_13_3_1/src/SUSYCascades/samples/NANO/";

  LabRecoFrame LAB_BDT("LAB_BDT","lab");
  SelfAssemblingRecoFrame BDT_CM("BDT_CM","CM");
  VisibleRecoFrame V1("V1","#vec{p}");
  CombinatoricGroup VIS("VIS","Visible Object Jigsaws");
  VIS.AddFrame(V1);
  VIS.SetNElementsForFrame(V1,1,false);
  LAB_BDT.SetChildFrame(BDT_CM);
  BDT_CM.AddChildFrame(V1);
  LAB_BDT.InitializeTree();
  LAB_BDT.InitializeAnalysis();

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

/*
   vector<string> datasets_list_2017 = {
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
   };
   filetag = "Fall17_102X";
   eventcount = "root/EventCount/EventCount_NANO_Fall17_102X.root";
   for(int i = 0; i < int(datasets_list_2017.size()); i++) { cout << "Weight Per Event: " << filetag << " " << datasets_list_2017[i] << " " << getWeightPerEvent(datasets_list_2017[i],filetag,eventcount) << endl; }
*/

   //std::map<string,double> WJets_dataset_list_2017 = {
   // {"Fall17_102X/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.00443 * 57.48 * 1.21 * 1.009},
   //};
   
   std::map<string,double> TTJets_dataset_list_2023 = {
    {"Summer23_130X/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.txt",1.},
   };

   std::map<string,double> TTJets_dataset_list_2017 = {
    {"Fall17_102X/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
    //{"Fall17_102X/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
    //{"Fall17_102X/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
    //{"Fall17_102X/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8.txt",1.},
   };

   std::map<string,double> TTZToQQ_dataset_list_2017 = {
    {"Fall17_102X/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8.txt",1.},
   };

   std::map<string,double> WZ_dataset_list_2017 = {
    {"Fall17_102X/WZ_TuneCP5_13TeV-pythia8.txt",1.},
   };

   std::map<string,double> WWZ_dataset_list_2017 = {
    {"Fall17_102X/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",1.},
   };

  std::map<string,double> Cascade_dataset_list_2023BPix = {
    //{"Summer23BPix_130X_SMS/local_cascades_40.txt",1.},
    {"Summer23BPix_130X_SMS/local_cascades_10.txt",1.},
  };

  std::map<string,double> TChiWZ_3_50_dataset_list_2017 = {
    {"Fall17_102X_SMS/TChiWZ_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,double> TChiWZ_60_90_dataset_list_2017 = {
    {"Fall17_102X_SMS/SMS-TChiWZ_dM-60to90_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8.txt",1.},
  };
  std::map<string,double> TChiWZ_100_Inf_dataset_list_2017 = {
    {"Fall17_102X_SMS/SMS-TChiWZ_TuneCP2_13TeV-madgraphMLM-pythia8.txt",1.},
  };

  std::map<string,std::map<string,double>> datasets = {
    {"ttbar",TTJets_dataset_list_2017},
    //{"TTZToQQ",TTZToQQ_dataset_list_2017},
    {"DB_WZ",WZ_dataset_list_2017},
    {"TB_WWZ",WWZ_dataset_list_2017},
    {"TChiWZ_20",TChiWZ_3_50_dataset_list_2017},
    {"TChiWZ_50",TChiWZ_3_50_dataset_list_2017},
    //{"TChiWZ_90",TChiWZ_60_90_dataset_list_2017},
    //{"TChiWZ_100",TChiWZ_100_Inf_dataset_list_2017},
    ////{"cascades",Cascade_dataset_list_2023BPix},
    //{"ttbar_2023",TTJets_dataset_list_2023},
  };
  int SKIP = 1;

  for(std::map<string,std::map<string,double>>::iterator iter1 = datasets.begin(); iter1 != datasets.end(); ++iter1)
  {
    int kept_sample_events = 0;

    int quark_cut = 0;
    int tot_quark = 0;

    double prong2 = 0;
    double prong3 = 0;
    double auntprong = 0;
    double lepprong = 0;
    double totcands = 0;
    double genHBosons = 0;
    int matchedAside = 0;
    int matchedBside = 0;
    int unmatchedAside = 0;
    int unmatchedBside = 0;

    string proc_Name = iter1->first;
    std::map<string,double> dataset = iter1->second;
    vector<string> files;
    vector<TH1D*> hists1;
    vector<TH2D*> hists2;
    vector<TEfficiency*> effs;
    string title = proc_Name;
    g_PlotTitle = title;

    double DM = 0.;
    if(proc_Name.find("TChi") != std::string::npos)
      DM = std::stod(proc_Name.substr(proc_Name.find("_") + 1));

    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET;").c_str(), g_NX, g_Xmin, g_Xmax);
    //hists1.push_back(hist_MET);
    TH1D* hist_V_had_cand = new TH1D((title+"_V_had_cand").c_str(), (title+"_V_had_cand;V_had_cand;").c_str(), 10, 0, 10);
    //hists1.push_back(hist_V_had_cand);

    TH2D* hist_Njets_v_V_had_cand = new TH2D((title+"_Njets_v_V_had_cand").c_str(), (title+"_Njets_v_V_had_cand;Njets;V_had_cand").c_str(), 20, 0, 20, 10, 0, 10);
    hists2.push_back(hist_Njets_v_V_had_cand);
    TH2D* hist_Njets_v_V_lep_cand = new TH2D((title+"_Njets_v_V_lep_cand").c_str(), (title+"_Njets_v_V_lep_cand;Njets;V_lep_cand").c_str(), 20, 0, 20, 5, 0, 5);
    //hists2.push_back(hist_Njets_v_V_lep_cand);
    TH2D* hist_V_had_cand_v_V_lep_cand = new TH2D((title+"_V_had_cand_v_V_lep_cand").c_str(), (title+"_V_had_cand_v_V_lep_cand;V_had_cand;V_lep_cand").c_str(), 10, 0, 10, 5, 0, 5);
    //hists2.push_back(hist_V_had_cand_v_V_lep_cand);
    TH2D* hist_Nleps_v_V_lep_cand = new TH2D((title+"_Nleps_v_V_lep_cand").c_str(), (title+"_Nleps_v_V_lep_cand;Nleps;V_lep_cand").c_str(), 5, 0, 5, 5, 0, 5);
    //hists2.push_back(hist_Nleps_v_V_lep_cand);

    TH1D* hist_NQuarks = new TH1D((title+"_NQuarks").c_str(), (title+"_NQuarks;NQuarks").c_str(), 6, 0, 6);
    //hists1.push_back(hist_NQuarks);

    TH1D* hist_MatchedVJetsPt = new TH1D((title+"_MatchedVJetsPt").c_str(), (title+"_MatchedVJetsPt;MatchedVJetsPt").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_MatchedVJetsPt);
    TH1D* hist_UnmatchedVJetsPt = new TH1D((title+"_UnmatchedVJetsPt").c_str(), (title+"_UnmatchedVJetsPt;UnmatchedVJetsPt").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_UnmatchedVJetsPt);
    TH1D* hist_MatchedVJetsEta = new TH1D((title+"_MatchedVJetsEta").c_str(), (title+"_MatchedVJetsEta;MatchedVJetsEta").c_str(), g_NX, -2.5, 2.5);
    hists1.push_back(hist_MatchedVJetsEta);
    TH1D* hist_UnmatchedVJetsEta = new TH1D((title+"_UnmatchedVJetsEta").c_str(), (title+"_UnmatchedVJetsEta;UnmatchedVJetsEta").c_str(), g_NX, -2.5, 2.5);
    hists1.push_back(hist_UnmatchedVJetsEta);
    TH1D* hist_MatchedVJetsMass = new TH1D((title+"_MatchedVJetsMass").c_str(), (title+"_MatchedVJetsMass;MatchedVJetsMass").c_str(), g_NX, 0, 50.);
    hists1.push_back(hist_MatchedVJetsMass);
    TH1D* hist_UnmatchedVJetsMass = new TH1D((title+"_UnmatchedVJetsMass").c_str(), (title+"_UnmatchedVJetsMass;UnmatchedVJetsMass").c_str(), g_NX, 0., 50.);
    hists1.push_back(hist_UnmatchedVJetsMass);

    TH1D* hist_MatchedVPt = new TH1D((title+"_MatchedVPt").c_str(), (title+"_MatchedVPt;MatchedVPt").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_MatchedVPt);
    TH1D* hist_UnmatchedVPt = new TH1D((title+"_UnmatchedVPt").c_str(), (title+"_UnmatchedVPt;UnmatchedVPt").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_UnmatchedVPt);
    TH1D* hist_MatchedVEta = new TH1D((title+"_MatchedVEta").c_str(), (title+"_MatchedVEta;MatchedVEta").c_str(), g_NX, -2.5, 2.5);
    hists1.push_back(hist_MatchedVEta);
    TH1D* hist_UnmatchedVEta = new TH1D((title+"_UnmatchedVEta").c_str(), (title+"_UnmatchedVEta;UnmatchedVEta").c_str(), g_NX, -2.5, 2.5);
    hists1.push_back(hist_UnmatchedVEta);
    TH1D* hist_MatchedVMass = new TH1D((title+"_MatchedVMass").c_str(), (title+"_MatchedVMass;MatchedVMass").c_str(), g_NX, 0, 110.);
    hists1.push_back(hist_MatchedVMass);
    TH1D* hist_UnmatchedVMass = new TH1D((title+"_UnmatchedVMass").c_str(), (title+"_UnmatchedVMass;UnmatchedVMass").c_str(), g_NX, 0., 110.);
    hists1.push_back(hist_UnmatchedVMass);
    TH1D* hist_PartialmatchedVMass = new TH1D((title+"_PartialmatchedVMass").c_str(), (title+"_PartialmatchedVMass;PartialmatchedVMass").c_str(), g_NX, 0, 110.);
    hists1.push_back(hist_PartialmatchedVMass);
    TH1D* hist_MatchedVCosDecayAngle = new TH1D((title+"_MatchedVCosDecayAngle").c_str(), (title+"_MatchedVCosDecayAngle;MatchedVCosDecayAngle").c_str(), g_NX, -1., 1.);
    hists1.push_back(hist_MatchedVCosDecayAngle);
    TH1D* hist_UnmatchedVCosDecayAngle = new TH1D((title+"_UnmatchedVCosDecayAngle").c_str(), (title+"_UnmatchedVCosDecayAngle;UnmatchedVCosDecayAngle").c_str(), g_NX, -1., 1.);
    hists1.push_back(hist_UnmatchedVCosDecayAngle);

    TH1D* hist_MatchedVProngDeltaPhi = new TH1D((title+"_MatchedVProngDeltaPhi").c_str(), (title+"_MatchedVProngDeltaPhi;MatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
    hists1.push_back(hist_MatchedVProngDeltaPhi);
    TH1D* hist_PartialmatchedVProngDeltaPhi = new TH1D((title+"_PartialmatchedVProngDeltaPhi").c_str(), (title+"_PartialmatchedVProngDeltaPhi;PartialmatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
    hists1.push_back(hist_PartialmatchedVProngDeltaPhi);
    TH1D* hist_UnmatchedVProngDeltaPhi = new TH1D((title+"_UnmatchedVProngDeltaPhi").c_str(), (title+"_UnmatchedVProngDeltaPhi;UnmatchedVProngDeltaPhi").c_str(), g_NX, -3.2, 3.2);
    hists1.push_back(hist_UnmatchedVProngDeltaPhi);
    TH2D* hist_MatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_MatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_MatchedV_Mass_ProngDeltaPhi;matched Mass;matched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_MatchedV_Mass_ProngDeltaPhi);
    TH2D* hist_PartialmatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_PartialmatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_PartialmatchedV_Mass_ProngDeltaPhi;partial matched Mass;partial matched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_PartialmatchedV_Mass_ProngDeltaPhi);
    TH2D* hist_UnmatchedV_Mass_ProngDeltaPhi = new TH2D((title+"_UnmatchedV_Mass_ProngDeltaPhi").c_str(), (title+"_UnmatchedV_Mass_ProngDeltaPhi;unmatched Mass;unmatched ProngDeltaPhi;").c_str(), g_NX/2., -3.2, 110., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_UnmatchedV_Mass_ProngDeltaPhi);
    TH1D* hist_MatchedVProngDeltaEta = new TH1D((title+"_MatchedVProngDeltaEta").c_str(), (title+"_MatchedVProngDeltaEta;MatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
    hists1.push_back(hist_MatchedVProngDeltaEta);
    TH1D* hist_PartialmatchedVProngDeltaEta = new TH1D((title+"_PartialmatchedVProngDeltaEta").c_str(), (title+"_PartialmatchedVProngDeltaEta;PartialmatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
    hists1.push_back(hist_PartialmatchedVProngDeltaEta);
    TH1D* hist_UnmatchedVProngDeltaEta = new TH1D((title+"_UnmatchedVProngDeltaEta").c_str(), (title+"_UnmatchedVProngDeltaEta;UnmatchedVProngDeltaEta").c_str(), g_NX, -3., 3.);
    hists1.push_back(hist_UnmatchedVProngDeltaEta);
    TH2D* hist_MatchedV_Mass_ProngDeltaEta = new TH2D((title+"_MatchedV_Mass_ProngDeltaEta").c_str(), (title+"_MatchedV_Mass_ProngDeltaEta;matched Mass;matched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
    hists2.push_back(hist_MatchedV_Mass_ProngDeltaEta);
    TH2D* hist_PartialmatchedV_Mass_ProngDeltaEta = new TH2D((title+"_PartialmatchedV_Mass_ProngDeltaEta").c_str(), (title+"_PartialmatchedV_Mass_ProngDeltaEta;partial matched Mass;partial matched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
    hists2.push_back(hist_PartialmatchedV_Mass_ProngDeltaEta);
    TH2D* hist_UnmatchedV_Mass_ProngDeltaEta = new TH2D((title+"_UnmatchedV_Mass_ProngDeltaEta").c_str(), (title+"_UnmatchedV_Mass_ProngDeltaEta;unmatched Mass;unmatched ProngDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., -3., 3.);
    hists2.push_back(hist_UnmatchedV_Mass_ProngDeltaEta);
    TH1D* hist_MatchedVProngAbsDeltaEta = new TH1D((title+"_MatchedVProngAbsDeltaEta").c_str(), (title+"_MatchedVProngAbsDeltaEta;MatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
    hists1.push_back(hist_MatchedVProngAbsDeltaEta);
    TH1D* hist_PartialmatchedVProngAbsDeltaEta = new TH1D((title+"_PartialmatchedVProngAbsDeltaEta").c_str(), (title+"_PartialmatchedVProngAbsDeltaEta;PartialmatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
    hists1.push_back(hist_PartialmatchedVProngAbsDeltaEta);
    TH1D* hist_UnmatchedVProngAbsDeltaEta = new TH1D((title+"_UnmatchedVProngAbsDeltaEta").c_str(), (title+"_UnmatchedVProngAbsDeltaEta;UnmatchedVProngAbsDeltaEta").c_str(), g_NX, 0., 3.);
    hists1.push_back(hist_UnmatchedVProngAbsDeltaEta);
    TH2D* hist_MatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_MatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_MatchedV_Mass_ProngAbsDeltaEta;matched Mass;matched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
    hists2.push_back(hist_MatchedV_Mass_ProngAbsDeltaEta);
    TH2D* hist_PartialmatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_PartialmatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_PartialmatchedV_Mass_ProngAbsDeltaEta;partial matched Mass;partial matched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
    hists2.push_back(hist_PartialmatchedV_Mass_ProngAbsDeltaEta);
    TH2D* hist_UnmatchedV_Mass_ProngAbsDeltaEta = new TH2D((title+"_UnmatchedV_Mass_ProngAbsDeltaEta").c_str(), (title+"_UnmatchedV_Mass_ProngAbsDeltaEta;unmatched Mass;unmatched ProngAbsDeltaEta;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.);
    hists2.push_back(hist_UnmatchedV_Mass_ProngAbsDeltaEta);

    TH2D* hist_MatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngDeltaPhi;matched CosDecayAngle_CM;matched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngDeltaPhi);
    TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi;partial matched CosDecayAngle_CM;partial matched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaPhi);
    TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi;unmatched CosDecayAngle_CM;unmatched ProngDeltaPhi;").c_str(), g_NX/2., -1., 1., g_NX/2., -3.2, 3.2);
    hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaPhi);
    TH2D* hist_MatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngDeltaEta;matched CosDecayAngle_CM;matched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
    hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngDeltaEta);
    TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta;partial matched CosDecayAngle_CM;partial matched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
    hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngDeltaEta);
    TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta;unmatched CosDecayAngle_CM;unmatched ProngDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., -3., 3.);
    hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngDeltaEta);
    TH2D* hist_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;matched CosDecayAngle_CM;matched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
    hists2.push_back(hist_MatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);
    TH2D* hist_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;partial matched CosDecayAngle_CM;partial matched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
    hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);
    TH2D* hist_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta = new TH2D((title+"_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta;unmatched CosDecayAngle_CM;unmatched ProngAbsDeltaEta;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 3.);
    hists2.push_back(hist_UnmatchedV_CosDecayAngle_CM_ProngAbsDeltaEta);

    TH2D* hist_MatchedV_Mass_CosDecayAngle = new TH2D((title+"_MatchedV_Mass_CosDecayAngle").c_str(), (title+"_MatchedV_Mass_CosDecayAngle;matched V mass; matched V CosDecayAngle;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_Mass_CosDecayAngle);
    TH2D* hist_UnmatchedV_Mass_CosDecayAngle = new TH2D((title+"_UnmatchedV_Mass_CosDecayAngle").c_str(), (title+"_UnmatchedV_Mass_CosDecayAngle;unmatched V mass; unmatched V CosDecayAngle;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_Mass_CosDecayAngle);

    TH2D* hist_MatchedV_PCM_PCand = new TH2D((title+"_MatchedV_PCM_PCand").c_str(), (title+"_MatchedV_PCM_PCand;PCM;matched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
    hists2.push_back(hist_MatchedV_PCM_PCand);
    TH2D* hist_UnmatchedV_PCM_PCand = new TH2D((title+"_UnmatchedV_PCM_PCand").c_str(), (title+"_UnmatchedV_PCM_PCand;PCM;unmatched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
    hists2.push_back(hist_UnmatchedV_PCM_PCand);
    TH2D* hist_PartialmatchedV_PCM_PCand = new TH2D((title+"_PartialmatchedV_PCM_PCand").c_str(), (title+"_PartialmatchedV_PCM_PCand;PCM;Partialmatched PCand;").c_str(), g_NX/2., 0., 3000., g_NX/2., 0., 2000.);
    hists2.push_back(hist_PartialmatchedV_PCM_PCand);

    TH1D* hist_UnmatchedV_CosDecayAngle_CM = new TH1D((title+"_UnmatchedV_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_CosDecayAngle_CM;unmatchedV_CosDecayAngle_CM").c_str(), g_NX/2., -1., 1.);
    //hists1.push_back(hist_UnmatchedV_CosDecayAngle_CM);
    TH1D* hist_UnmatchedV_DeltaPhiDecayAngle_CM = new TH1D((title+"_UnmatchedV_DeltaPhiDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayAngle_CM;unmatchedV_DeltaPhiDecayAngle_CM").c_str(), g_NX/2., 0., 3.5);
    //hists1.push_back(hist_UnmatchedV_DeltaPhiDecayAngle_CM);
    TH1D* hist_UnmatchedV_DeltaPhiDecayPlanes_CM = new TH1D((title+"_UnmatchedV_DeltaPhiDecayPlanes_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayPlanes_CM;unmatchedV_DeltaPhiDecayPlanes_CM").c_str(), g_NX/2., 0., 6.5);
    //hists1.push_back(hist_UnmatchedV_DeltaPhiDecayPlanes_CM);
    TH2D* hist_UnmatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_Mass_CosDecayAngle_CM;unmatched Mass;unmatched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_Mass_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_UnmatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_UnmatchedV_Mass_DeltaPhiDecayAngle_CM;unmatched Mass;unmatched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
    hists2.push_back(hist_UnmatchedV_Mass_DeltaPhiDecayAngle_CM);
    TH2D* hist_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM;unmatched Mass;unmatched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
    hists2.push_back(hist_UnmatchedV_Mass_DeltaPhiDecayPlanes_CM);

    TH2D* hist_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM = new TH2D((title+"_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM;DeltaPhiDecayAngle_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 3.5, g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM;DeltaPhiDecayAngle_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 3.5, g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_DeltaPhiDecayAngle_CM_CosDecayAngle_CM);
    TH2D* hist_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM = new TH2D((title+"_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM;DeltaPhiDecayPlanes_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 6.5, g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM;DeltaPhiDecayPlanes_CM;CosDecayAngle_CM").c_str(), g_NX/2., 0., 6.5, g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_DeltaPhiDecayPlanes_CM_CosDecayAngle_CM);

    TH1D* hist_MatchedV_CosDecayAngle_CM = new TH1D((title+"_MatchedV_CosDecayAngle_CM").c_str(), (title+"_MatchedV_CosDecayAngle_CM;matchedV_CosDecayAngle_CM").c_str(), g_NX/2., -1., 1.);
    //hists1.push_back(hist_MatchedV_CosDecayAngle_CM);
    TH1D* hist_MatchedV_DeltaPhiDecayAngle_CM = new TH1D((title+"_MatchedV_DeltaPhiDecayAngle_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayAngle_CM;matchedV_DeltaPhiDecayAngle_CM").c_str(), g_NX/2., 0., 3.5);
    //hists1.push_back(hist_MatchedV_DeltaPhiDecayAngle_CM);
    TH1D* hist_MatchedV_DeltaPhiDecayPlanes_CM = new TH1D((title+"_MatchedV_DeltaPhiDecayPlanes_CM").c_str(), (title+"_MatchedV_DeltaPhiDecayPlanes_CM;matchedV_DeltaPhiDecayPlanes_CM").c_str(), g_NX/2., 0., 6.5);
    //hists1.push_back(hist_MatchedV_DeltaPhiDecayPlanes_CM);
    TH2D* hist_MatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_MatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_MatchedV_Mass_CosDecayAngle_CM;matched Mass;matched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_Mass_CosDecayAngle_CM);
    TH2D* hist_MatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_MatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_MatchedV_Mass_DeltaPhiDecayAngle_CM;matched Mass;matched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
    hists2.push_back(hist_MatchedV_Mass_DeltaPhiDecayAngle_CM);
    TH2D* hist_MatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_MatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_MatchedV_Mass_DeltaPhiDecayPlanes_CM;matched Mass;matched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
    hists2.push_back(hist_MatchedV_Mass_DeltaPhiDecayPlanes_CM);
    
    TH2D* hist_PartialmatchedV_Mass_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_Mass_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_Mass_CosDecayAngle_CM;partial matched Mass;partial matched CosDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_Mass_CosDecayAngle_CM);
    TH2D* hist_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM = new TH2D((title+"_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM").c_str(), (title+"_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM;partial matched Mass;partial matched DeltaPhiDecayAngle_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 3.5);
    hists2.push_back(hist_PartialmatchedV_Mass_DeltaPhiDecayAngle_CM);
    TH2D* hist_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM = new TH2D((title+"_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM").c_str(), (title+"_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM;partial matched Mass;partial matched DeltaPhiDecayPlanes_CM;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 6.5);
    hists2.push_back(hist_PartialmatchedV_Mass_DeltaPhiDecayPlanes_CM);

    TH2D* hist_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM;partial matched CosDecayAngle;partial matched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_CosDecayAngle_CosDecayAngle_CM);
    TH2D* hist_MatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_MatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_MatchedV_CosDecayAngle_CosDecayAngle_CM;matched CosDecayAngle;matched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_CosDecayAngle_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_CosDecayAngle_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_CosDecayAngle_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_CosDecayAngle_CosDecayAngle_CM;unmatched CosDecayAngle;unmatched CosDecayAngle_CM;").c_str(), g_NX/2., -1., 1., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_CosDecayAngle_CosDecayAngle_CM);

    TH1D* hist_Matched_CosDecayAngle_Hemi = new TH1D((title+"_Matched_CosDecayAngle_Hemi").c_str(), (title+"_Matched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
    hists1.push_back(hist_Matched_CosDecayAngle_Hemi);
    TH1D* hist_Partialmatched_CosDecayAngle_Hemi = new TH1D((title+"_Partialmatched_CosDecayAngle_Hemi").c_str(), (title+"_Partialmatched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
    hists1.push_back(hist_Partialmatched_CosDecayAngle_Hemi);
    TH1D* hist_Unmatched_CosDecayAngle_Hemi = new TH1D((title+"_Unmatched_CosDecayAngle_Hemi").c_str(), (title+"_Unmatched_CosDecayAngle_Hemi;CosDecayAngle Hemi;").c_str(), g_NX/2., -1., 1.);
    hists1.push_back(hist_Unmatched_CosDecayAngle_Hemi);
    TH2D* hist_Matched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Matched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Matched_Mass_CosDecayAngle_Hemi;matched Mass;matched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_Matched_Mass_CosDecayAngle_Hemi);
    TH2D* hist_Unmatched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Unmatched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Unmatched_Mass_CosDecayAngle_Hemi;unmatched Mass;unmatched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_Unmatched_Mass_CosDecayAngle_Hemi);
    TH2D* hist_Partialmatched_Mass_CosDecayAngle_Hemi = new TH2D((title+"_Partialmatched_Mass_CosDecayAngle_Hemi").c_str(), (title+"_Partialmatched_Mass_CosDecayAngle_Hemi;partial matched Mass;partial matched CosDecayAngle_Hemi;").c_str(), g_NX/2., 0., 110., g_NX/2., -1., 1.);
    hists2.push_back(hist_Partialmatched_Mass_CosDecayAngle_Hemi);

    TH1D* hist_quarkPt = new TH1D((title+"_quarkPt").c_str(), (title+"_quarkPt;quarkPt").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_quarkPt);
    TH1D* hist_quarkEta = new TH1D((title+"_quarkEta").c_str(), (title+"_quarkEta;quarkEta").c_str(), g_NX, -5., 5.);
    hists1.push_back(hist_quarkEta);

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

    TH1D* hist_CandCount = new TH1D((title+"_CandCount").c_str(), (title+"_CandCount;").c_str(), 5, 0, 5);
    TH1D* hist_CandCount_Combo = new TH1D((title+"_CandCount_Combo").c_str(), (title+"_CandCount_Combo;").c_str(), 5, 0, 5);
    TH1D* hist_CandCount_Chi2 = new TH1D((title+"_CandCount_Chi2").c_str(), (title+"_CandCount_Chi2;").c_str(), 5, 0, 5);

    TH2D* hist_Matched_DeltaR_P = new TH2D((title+"_Matched_DeltaR_P").c_str(), (title+"_Matched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
    //hists2.push_back(hist_Matched_DeltaR_P);
    TH2D* hist_Unmatched_DeltaR_P = new TH2D((title+"_Unmatched_DeltaR_P").c_str(), (title+"_Unmatched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
    //hists2.push_back(hist_Unmatched_DeltaR_P);

    TH2D* hist_combo_Matched_DeltaR_P = new TH2D((title+"_combo_Matched_DeltaR_P").c_str(), (title+"_combo_Matched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
    //hists2.push_back(hist_combo_Matched_DeltaR_P);
    TH2D* hist_combo_Unmatched_DeltaR_P = new TH2D((title+"_combo_Unmatched_DeltaR_P").c_str(), (title+"_combo_Unmatched_DeltaR_P;DeltaR;P;").c_str(), g_NX/2., 0., 4., g_NX/2., 0., 500.);
    //hists2.push_back(hist_combo_Unmatched_DeltaR_P);

    TH2D* hist_Matched_lnDeltaR_lnP = new TH2D((title+"_Matched_lnDeltaR_lnP").c_str(), (title+"_Matched_lnDeltaR_lnP;lnDeltaR;lnP;").c_str(), g_NX/2., 0., 1.5, g_NX/2., 0., 7.);
    //hists2.push_back(hist_Matched_lnDeltaR_lnP);
    TH2D* hist_Unmatched_lnDeltaR_lnP = new TH2D((title+"_Unmatched_lnDeltaR_lnP").c_str(), (title+"_Unmatched_lnDeltaR_lnP;lnDeltaR;lnP;").c_str(), g_NX/2., 0., 1.5, g_NX/2., 0., 7.);
    //hists2.push_back(hist_Unmatched_lnDeltaR_lnP);

    TH2D* hist_Matched_Mass_P_Multi_DeltaR = new TH2D((title+"_Matched_Mass_P_Multi_DeltaR").c_str(), (title+"_Matched_Mass_P_Multi_DeltaR;Mass;P*DeltaR;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 1000.);
    hists2.push_back(hist_Matched_Mass_P_Multi_DeltaR);
    TH2D* hist_Unmatched_Mass_P_Multi_DeltaR = new TH2D((title+"_Unmatched_Mass_P_Multi_DeltaR").c_str(), (title+"_Unmatched_Mass_P_Multi_DeltaR;Mass;P*DeltaR;").c_str(), g_NX/2., 0., 110., g_NX/2., 0., 1000.);
    hists2.push_back(hist_Unmatched_Mass_P_Multi_DeltaR);

    TH2D* hist_Matched_CosDecayAngle_CM_P_Multi_DeltaR = new TH2D((title+"_Matched_CosDecayAngle_CM_P_Multi_DeltaR").c_str(), (title+"_Matched_CosDecayAngle_CM_P_Multi_DeltaR;CosDecayAngle_CM;P*DeltaR;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 1000.);
    hists2.push_back(hist_Matched_CosDecayAngle_CM_P_Multi_DeltaR);
    TH2D* hist_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR = new TH2D((title+"_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR").c_str(), (title+"_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR;CosDecayAngle_CM;P*DeltaR;").c_str(), g_NX/2., -1., 1., g_NX/2., 0., 1000.);
    hists2.push_back(hist_Unmatched_CosDecayAngle_CM_P_Multi_DeltaR);

    TH1D* hist_BaconNum = new TH1D((title+"_BaconNum").c_str(), (title+"_BaconNum;BaconNum").c_str(), 10, 0, 10);
    hists1.push_back(hist_BaconNum);

    TH1D* hist_MatchedlepCandBaconNum = new TH1D((title+"_MatchedlepCandBaconNum").c_str(), (title+"_MatchedlepCandBaconNum;MatchedlepCandBaconNum").c_str(), 10, 0, 10);
    hists1.push_back(hist_MatchedlepCandBaconNum);
    TH1D* hist_PartialMatchedlepCandBaconNum = new TH1D((title+"_PartialMatchedlepCandBaconNum").c_str(), (title+"_PartialMatchedlepCandBaconNum;PartialMatchedlepCandBaconNum").c_str(), 10, 0, 10);
    hists1.push_back(hist_PartialMatchedlepCandBaconNum);
    TH1D* hist_UnmatchedlepCandBaconNum = new TH1D((title+"_UnmatchedlepCandBaconNum").c_str(), (title+"_UnmatchedlepCandBaconNum;UnmatchedlepCandBaconNum").c_str(), 10, 0, 10);
    hists1.push_back(hist_UnmatchedlepCandBaconNum);

    TH2D* hist_MatchedlepCandBaconNum_Njets = new TH2D((title+"_MatchedlepCandBaconNum_Njets").c_str(), (title+"_MatchedlepCandBaconNum;MatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
    hists2.push_back(hist_MatchedlepCandBaconNum_Njets);
    TH2D* hist_UnmatchedlepCandBaconNum_Njets = new TH2D((title+"_UnmatchedlepCandBaconNum_Njets").c_str(), (title+"_UnmatchedlepCandBaconNum;UnmatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
    hists2.push_back(hist_UnmatchedlepCandBaconNum_Njets);
    TH2D* hist_PartialMatchedlepCandBaconNum_Njets = new TH2D((title+"_PartialMatchedlepCandBaconNum_Njets").c_str(), (title+"_PartialMatchedlepCandBaconNum;PartialMatchedlepCandBaconNum;Njets").c_str(), 10, 0, 10, 15, 0, 15);
    hists2.push_back(hist_PartialMatchedlepCandBaconNum_Njets);

    TH1D* hist_chi2WMass = new TH1D((title+"_chi2WMass").c_str(), (title+"_chi2WMass;chi2WMass").c_str(), g_NX, 30., 120.);
    //hists1.push_back(hist_chi2WMass);

    TH2D* hist_MatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_MatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_MatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
    hists2.push_back(hist_MatchedCand_JetA_Mass_JetB_Mass);
    TH2D* hist_UnmatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_UnmatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_UnmatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
    hists2.push_back(hist_UnmatchedCand_JetA_Mass_JetB_Mass);
    TH2D* hist_PartialMatchedCand_JetA_Mass_JetB_Mass = new TH2D((title+"_PartialMatchedCand_JetA_Mass_JetB_Mass").c_str(), (title+"_PartialMatchedCand_JetA_Mass_JetB_Mass;JetA_Mass;JetB_Mass").c_str(), g_NX/2., 0., 60., g_NX/2., 0., 40.);
    hists2.push_back(hist_PartialMatchedCand_JetA_Mass_JetB_Mass);

    // Prong Mass Ratio
    TH2D* hist_MatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_MatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_MatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedCand_PMR_CosDecayAngle_CM);
    TH2D* hist_PartialMatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_PartialMatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_PartialMatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialMatchedCand_PMR_CosDecayAngle_CM);
    TH2D* hist_UnmatchedCand_PMR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedCand_PMR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedCand_PMR_CosDecayAngle_CM;PMR;CosDecayAngle_CM").c_str(), g_NX/2., 0., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedCand_PMR_CosDecayAngle_CM);

    TH2D* hist_MatchedCand_PMR_Mass = new TH2D((title+"_MatchedCand_PMR_Mass").c_str(), (title+"_MatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
    hists2.push_back(hist_MatchedCand_PMR_Mass);
    TH2D* hist_PartialMatchedCand_PMR_Mass = new TH2D((title+"_PartialMatchedCand_PMR_Mass").c_str(), (title+"_PartialMatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
    hists2.push_back(hist_PartialMatchedCand_PMR_Mass);
    TH2D* hist_UnmatchedCand_PMR_Mass = new TH2D((title+"_UnmatchedCand_PMR_Mass").c_str(), (title+"_UnmatchedCand_PMR_Mass;PMR;Mass").c_str(), g_NX/2., 0., 10., g_NX/2., 0., 110.);
    hists2.push_back(hist_UnmatchedCand_PMR_Mass);


    TEfficiency* eff_quarkMatchedGenJet_eta = new TEfficiency((title+"_eff_quarkMatchedGenJet_eta").c_str(),"Efficiency of Quark Getting Matched to GenJet;quarkEta;Quark Matching Efficiency", g_NX, -5., 5.);
    //effs.push_back(eff_quarkMatchedGenJet_eta);
    TEfficiency* eff_GenJetMatchedAll_V_cand_eta = new TEfficiency((title+"_eff_GenJetMatchedAll_V_cand_eta").c_str(),"Efficiency of V cand having all Matched GenJet;V_candEta;Efficiency", g_NX, -5., 5.);
    //effs.push_back(eff_GenJetMatchedAll_V_cand_eta);
    TEfficiency* eff_GenJetMatchedAll_V_cand_mass = new TEfficiency((title+"_eff_GenJetMatchedAll_V_cand_mass").c_str(),"Efficiency of V cand having all Matched GenJet;V_candMass;Efficiency", g_NX, 0., 110.);
    //effs.push_back(eff_GenJetMatchedAll_V_cand_mass);

    TEfficiency* eff_bosonMatched_MatchedV_genboson_Mass = new TEfficiency((title+"_eff_bosonMatched_MatchedV_genboson_Mass").c_str(),"Efficiency of V Cand Matching to Gen Boson;genBosonMass;Efficiency", g_NX, 50., 110.);
    //effs.push_back(eff_bosonMatched_MatchedV_genboson_Mass);
    TEfficiency* eff_bosonMatched_MatchedV_genboson_Eta = new TEfficiency((title+"_eff_bosonMatched_MatchedV_genboson_Eta").c_str(),"Efficiency of V Cand Matching to Gen Boson;genBosonEta;Efficiency", g_NX, -5., 5.);
    //effs.push_back(eff_bosonMatched_MatchedV_genboson_Eta);
    TEfficiency* eff_bosonMatched_UnmatchedV_genboson_Mass = new TEfficiency((title+"_eff_bosonMatched_UnmatchedV_genboson_Mass").c_str(),"Efficiency of V Cand Not Matching to Gen Boson;genBosonMass;Efficiency", g_NX, 50., 110.);
    //effs.push_back(eff_bosonMatched_UnmatchedV_genboson_Mass);
    TEfficiency* eff_bosonMatched_UnmatchedV_genboson_Eta = new TEfficiency((title+"_eff_bosonMatched_UnmatchedV_genboson_Eta").c_str(),"Efficiency of V Cand Not Matching to Gen Boson;genBosonEta;Efficiency", g_NX, -5., 5.);
    //effs.push_back(eff_bosonMatched_UnmatchedV_genboson_Eta);


    TH2D* hist_MatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_EMR_CosDecayAngle_CM);
    TH2D* hist_PartialmatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_EMR_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_EMR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EMR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EMR_CosDecayAngle_CM;EMR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_EMR_CosDecayAngle_CM);
    TH2D* hist_MatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_EPR_CosDecayAngle_CM);
    TH2D* hist_PartialmatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_EPR_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_EPR_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EPR_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EPR_CosDecayAngle_CM;EPR;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_EPR_CosDecayAngle_CM);
    TH2D* hist_MatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_EMR_other_CosDecayAngle_CM);
    TH2D* hist_PartialmatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_EMR_other_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_EMR_other_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EMR_other_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EMR_other_CosDecayAngle_CM;EMR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_EMR_other_CosDecayAngle_CM);
    TH2D* hist_MatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_MatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_MatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_MatchedV_EPR_other_CosDecayAngle_CM);
    TH2D* hist_PartialmatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_PartialmatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_PartialmatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_PartialmatchedV_EPR_other_CosDecayAngle_CM);
    TH2D* hist_UnmatchedV_EPR_other_CosDecayAngle_CM = new TH2D((title+"_UnmatchedV_EPR_other_CosDecayAngle_CM").c_str(), (title+"_UnmatchedV_EPR_other_CosDecayAngle_CM;EPR_other;CosDecayAngle_CM;").c_str(), g_NX/2., -1., 10., g_NX/2., -1., 1.);
    hists2.push_back(hist_UnmatchedV_EPR_other_CosDecayAngle_CM);


    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
    hists2.push_back(hist_RISR_PTISR);
    TH2D* hist_RISR_PTISR_Cands = new TH2D((title+"_RISR_PTISR_Cands").c_str(), (title+"_RISR_PTISR;RISR;PTISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
    hists2.push_back(hist_RISR_PTISR_Cands);

    for(std::map<string,double>::iterator iter2 = dataset.begin(); iter2 != dataset.end(); ++iter2)
    {
      files.clear();
      string dataset_name = iter2->first;
      double weight = iter2->second;
      load_files_from_list(files,path_to_lists+dataset_name);
      dataset_name = get_str_between_two_str(dataset_name,"/",".");
      for(int f = 0; f < int(files.size()); f++) 
      {
        //if(proc_Name.find("TChi") != std::string::npos){
        //  if(f>=num_files) break;
        //}
        //else if(f>=num_files/2 && f>0) break;
        if(f>=num_files) break;
        TChain* chain = new TChain("Events");
        chain->Add(files[f].c_str());
        
        SUSYNANOBase* base = new SUSYNANOBase(chain); //Run2
        //NANORun3* base = new NANORun3(chain); // Run3
        //AnalysisBase<NANORun3>* base = new AnalysisBase<NANORun3>(chain); // do not use
        gErrorIgnoreLevel = kFatal;
        int Nentry = base->fChain->GetEntries();
        gErrorIgnoreLevel = 0;
        int kept_file_events = 0;
        // event loop
        for(int e = 0; e < Nentry; e += SKIP){
          base->GetEntry(e);
//std::cout << "NEW EVENT" << std::endl;
          if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
            std::cout << "      event " << e << " | " << Nentry << std::endl;
          double hweight = base->genWeight*weight*double(SKIP);
          int N = base->nGenPart;
          int PDGID;
          int MP = 0;
          int MC = 0;
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

          if(MP != 0 && MC != 0 && ((MP - MC) > (DM + 1.) || (MP - MC) < (DM - 1.))) continue;

          int year = 2017;
          // loop over jets 
          ParticleList jets;
          bool passID = true;
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
            if(base->Jet_jetId[i] < 3)
              continue;
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

	      if(base->Electron_mvaFall17V2noIso_WPL[i]){ // run2 2017
	      //if(base->Electron_mvaNoIso[i]){ // run3
                //if(lep.ParticleID() < kMedium || lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
                if(lep.MiniIso()*lep.Pt() >= 4. || lep.RelIso()*lep.Pt() >= 4.)
                       lep.SetLepQual(kBronze);
                     else if(lep.SIP3D() > 2.)
                       lep.SetLepQual(kSilver);
                     else
                       lep.SetLepQual(kGold);
                 electrons.push_back(lep);
               }
             }
           }
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
               int mom = base->GenPart_genPartIdxMother[i];
               if(mom >= 0 && mom < N){
         
         	int momID = base->GenPart_pdgId[mom];
         	int momStatus = base->GenPart_status[mom];
         
         	while(abs(momID) == 11){
         
         	  if(momStatus == 23){
         	    lep.SetMomPDGID(PDGID);
         	    lep.SetSourceID(GetLepSource(PDGID, PDGID, PDGID));
         	    break;
         	  }
         
         	  mom = base->GenPart_genPartIdxMother[mom];
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
               int mom = base->GenPart_genPartIdxMother[i];
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
         	  mom = base->GenPart_genPartIdxMother[mom];
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
               int mom = base->GenPart_genPartIdxMother[i];
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
            int mom = base->GenPart_genPartIdxMother[i];
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
//if(boson_acceptance_cut && (base->GenPart_pt[i] < 20. || abs(base->GenPart_eta[i]) > 2.4)) break;
if(boson_acceptance_cut && (base->GenPart_pt[i] < 5. || abs(base->GenPart_eta[i]) > 2.7)) break;
              //if(abs(momID) == 24 && base->GenPart_status[mom] == 22){
              if((abs(momID) == 24 || abs(momID) == 23) && base->GenPart_status[mom] == 22){
                int grandmom = base->GenPart_genPartIdxMother[mom];
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
                  grandmom = base->GenPart_genPartIdxMother[grandmom];
                }
              }
              mom = base->GenPart_genPartIdxMother[mom];
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
          gen_electrons = gen_electrons.PtEtaCut(1.,2.7);
          gen_muons = gen_muons.PtEtaCut(1.,2.7);
          ParticleList gen_leptons = gen_electrons+gen_muons;
          gen_leptons.SortByPt();

          int Njets = int(jets.size());
          int Nleps = int(leptons.size());
          int Nobjs = Njets + Nleps;
          int Ngen_jets = int(gen_jets.size());
          int Ngen_leps = int(gen_leptons.size());
          //if(Njets < 3) continue;

          int Ngen_boson_leps = 0;
          for(int i = 0; i < Ngen_leps; i++){
            //if(abs(gen_leptons[i].MomPDGID()) == 23 || abs(gen_leptons[i].MomPDGID()) == 24){
            if(abs(gen_leptons[i].MomPDGID()) == 24){
              Ngen_boson_leps++;
            }
          }
//if(gen_lepton_cut && Ngen_boson_leps < 1) continue; // only keep events with exactly one lepton coming from W
if(gen_lepton_cut && Ngen_boson_leps != 1) continue; // only keep events with exactly one lepton coming from W

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


// matched tagging analysis
          std::vector<V_Cand> V_had_cands; // binary decay tree
          std::vector<ParticleList> V_lep_cands;
          LAB_BDT.ClearEvent();
          vector<int>   jet_BDT_singlet;
          vector<RFKey> jetID_BDT;
          vector<RFKey> lepID_BDT;
          vector<RFKey> objID_BDT;
          TVector3 MET(0.,0.,0.);
          MET.SetPtEtaPhi(base->MET_pt,0.,base->MET_phi);
          MET.SetZ(0.);

// use reco jets for RJR
          if(!use_gen_jets){
            for(int i = 0; i < Njets; i++){
              if(Zero_RJR_BDT_JetMass)
                jets[i].SetPtEtaPhiM(jets[i].Pt(),jets[i].Eta(),jets[i].Phi(),0.);
              RFKey key = VIS.AddLabFrameFourVector(jets[i]);
              jetID_BDT.push_back(key);
              objID_BDT.push_back(key);
            }
          }

// use gen jets for RJR
          if(use_gen_jets){
            for(int i = 0; i < Ngen_jets; i++){
              // use gen jets except b-jets from tops for RJR
              // if(abs(gen_jets[i].MomPDGID()) == 6) continue;
              if(Zero_RJR_BDT_JetMass)
                gen_jets[i].SetPtEtaPhiM(gen_jets[i].Pt(),gen_jets[i].Eta(),gen_jets[i].Phi(),0.);
              RFKey key = VIS.AddLabFrameFourVector(gen_jets[i]);
              jetID_BDT.push_back(key);
              objID_BDT.push_back(key);
            }
          }


          for(int i = 0; i < Nleps; i++){
            RFKey key = VIS.AddLabFrameFourVector(leptons[i]);
            lepID_BDT.push_back(key);
            objID_BDT.push_back(key);
          }
          if(!LAB_BDT.AnalyzeEvent()) cout << "Problem with RF BDT Analyze Event \n";

          //for(int i = 0; i < Nleps; i++){
          //  ParticleList V_cand;
          //  const RestFrame& frame = BDT_CM.GetFrame(lepID_BDT[i]);
          //  const RestFrame& sibling = frame.GetSiblingFrame();
          //  const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
          //  for(int j = i+1; j < Nleps; j++){
          //    if(sibling == BDT_CM.GetFrame(lepID_BDT[j]) && frame.GetProductionFrame().GetMass() < 100.){
          //    //if(sibling == BDT_CM.GetFrame(lepID_BDT[j]) && frame.GetProductionFrame().GetMass() < 100. && frame.GetProductionFrame().GetMass() > 70.){
          //      V_cand.push_back(leptons[i]);
          //      V_cand.push_back(leptons[j]);
          //      V_lep_cands.push_back(V_cand); // 2-prong "sibling"
          //    }
          //    else if(aunt == BDT_CM.GetFrame(lepID_BDT[j]) && aunt.GetProductionFrame().GetMass() < 100.){
          //    //else if(aunt == BDT_CM.GetFrame(lepID_BDT[j]) && aunt.GetProductionFrame().GetMass() < 100. && aunt.GetProductionFrame().GetMass() > 70.){
          //      V_cand.push_back(leptons[i]);
          //      V_cand.push_back(leptons[j]);
          //      //V_lep_cands.push_back(V_cand); // 3-prong "aunt"
          //    }
          //    else if(frame.GetMass() < 100.){
          //    //else if(frame.GetMass() < 100. && frame.GetMass() > 70.){
          //      V_cand.push_back(leptons[i]); // 1-prong "orphan"
          //      //V_lep_cands.push_back(V_cand);
          //    }
          //  }
          //}

          double ECM = 0.; // sum of energy of all objects in CM frame
          double ECM_LAB_BDT = 0.; // sum of energy of all objects in CM frame evalutated in LAB_BDT frame
          for(int i = 0; i < Nobjs; i++){
            ECM += BDT_CM.GetFrame(objID_BDT[i]).GetFourVector(BDT_CM).E();
            ECM_LAB_BDT += BDT_CM.GetFrame(objID_BDT[i]).GetFourVector().E();
          }
          

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
          for(int i = 0; i < N_V_had_combo; i++){
            bool matched = false; // perfect match
            bool unmatched = true; // both jets are radiative
            bool found24 = false;
            bool found23 = false;
            bool found6 = false;
            //if(V_had_cands_combo[i].size() == 2 && abs(V_had_cands_combo[i][0].MomPDGID()) == 24 && V_had_cands_combo[i][0].GenMomIndex() == V_had_cands_combo[i][1].GenMomIndex()){
            if(V_had_cands_combo[i].size() == 2 && (abs(V_had_cands_combo[i][0].MomPDGID()) == 23 || abs(V_had_cands_combo[i][0].MomPDGID()) == 24) && V_had_cands_combo[i][0].GenMomIndex() == V_had_cands_combo[i][1].GenMomIndex()){
              V_had_cands_combo[i].SetMatch(kMatched);
              hist_CandCount_Combo->SetBinContent(1,hist_CandCount_Combo->GetBinContent(1)+1);
              matched = true;
              unmatched = false;
            }
            else{
              for(int j = 0; j < int(V_had_cands_combo[i].size()); j++){
                if(V_had_cands_combo[i][j].MomPDGID() == 24){
                  if(found24)
                    hist_CandCount_Combo->SetBinContent(2,hist_CandCount_Combo->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount_Combo->SetBinContent(3,hist_CandCount_Combo->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount_Combo->SetBinContent(4,hist_CandCount_Combo->GetBinContent(4)-1);
                  V_had_cands_combo[i].SetMatch(kW);
                  hist_CandCount_Combo->SetBinContent(2,hist_CandCount_Combo->GetBinContent(2)+1);
                  found24 = true;
                  unmatched = false;
                }
                else if(V_had_cands_combo[i][j].MomPDGID() == 23){
                  if(found24)
                    hist_CandCount_Combo->SetBinContent(2,hist_CandCount_Combo->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount_Combo->SetBinContent(3,hist_CandCount_Combo->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount_Combo->SetBinContent(4,hist_CandCount_Combo->GetBinContent(4)-1);
                  V_had_cands_combo[i].SetMatch(kZ);
                  hist_CandCount_Combo->SetBinContent(3,hist_CandCount_Combo->GetBinContent(3)+1);
                  found23 = true;
                  unmatched = false;
                }
                else if(V_had_cands_combo[i][j].MomPDGID() == 6){
                  if(found24)
                    hist_CandCount_Combo->SetBinContent(2,hist_CandCount_Combo->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount_Combo->SetBinContent(3,hist_CandCount_Combo->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount_Combo->SetBinContent(4,hist_CandCount_Combo->GetBinContent(4)-1);
                  V_had_cands_combo[i].SetMatch(kB);
                  hist_CandCount_Combo->SetBinContent(4,hist_CandCount_Combo->GetBinContent(4)+1);
                  found6 = true;
                  unmatched = false;
                }
              }
            }
            if(unmatched){
              V_had_cands_combo[i].SetMatch(kUnmatched);
              hist_CandCount_Combo->SetBinContent(5,hist_CandCount_Combo->GetBinContent(5)+1);
            }
          }

          if(use_gen_jets){
            for(int i = 0; i < Ngen_jets; i++){
              ParticleList V_cand_Part;
              ConstRestFrameList V_cand_RF;
              const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[i]);
              const RestFrame& sibling = frame.GetSiblingFrame();
              const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
              if(int(lepID_BDT.size()) == Nleps){
                for(int j = 0; j < Nleps; j++){
                  if(sibling == BDT_CM.GetFrame(lepID_BDT[j])){
                    for(int k = i+1; k < Ngen_jets; k++){
                      if(aunt == BDT_CM.GetFrame(jetID_BDT[k])){
                        if(frame.GetProductionFrame().GetMass() - sibling.GetFourVector(frame.GetProductionFrame()).M() < 100.){
                          V_cand_Part.push_back(gen_jets[i]);
                          V_cand_Part.push_back(gen_jets[k]);
                          V_cand_RF.Add(frame);
                          V_cand_RF.Add(aunt);
                          V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                          if(use_lep_prong){
                            if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                              V_had_cands.push_back(cand);
                              lepprong++;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              for(int j = i+1; j < Ngen_jets; j++){
                if(sibling == BDT_CM.GetFrame(jetID_BDT[j])){
                  if(frame.GetProductionFrame().GetMass() < 100.){
                    V_cand_Part.clear();
                    V_cand_RF.Clear();
                    V_cand_Part.push_back(gen_jets[i]);
                    V_cand_Part.push_back(gen_jets[j]); // 2-prong "sibling"
                    V_cand_RF.Add(frame);
                    V_cand_RF.Add(sibling);
                    V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                    if(use_prong2){
                      if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                        V_had_cands.push_back(cand);
                        prong2++;
                      }
                    }
                  }
                  //else if(frame.GetMass() < 100.){
                  //    V_cand_Part.push_back(gen_jets[i]); // 1-prong "orphan"
                  //    V_cand_RF.Add(frame);
                  //  }
                  //for(int k = j+1; k < Ngen_jets; k++){
                    //if(aunt == BDT_CM.GetFrame(jetID_BDT[k]) && aunt.GetProductionFrame().GetMass() < 100.){
                    //    V_cand_Part.push_back(gen_jets[k]); // 3-prong "aunt"
                    //    V_cand_RF.Add(aunt);
                    //    prong3++;
                    //    prong2--;
                    //}
                  //}
                }
                if(aunt == BDT_CM.GetFrame(jetID_BDT[j])){
                  bool lep_sibling = false;
                  if(int(lepID_BDT.size()) == Nleps){
                    for(int k = 0; k < Nleps; k++){
                      if(sibling == BDT_CM.GetFrame(lepID_BDT[k])){
                        lep_sibling = true;
                      }
                    }
                  }
                  if(lep_sibling) continue;
                  if((frame.GetMass() + aunt.GetMass()) < 100.){
                    V_cand_Part.clear();
                    V_cand_RF.Clear();
                    V_cand_Part.push_back(gen_jets[i]);
                    V_cand_Part.push_back(gen_jets[j]); // 2-prong "aunt"
                    V_cand_RF.Add(frame);
                    V_cand_RF.Add(aunt);
                    V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                    if(use_aunt_prong){
                      if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                        V_had_cands.push_back(cand);
                        auntprong++;
                      }
                    }
                  }
                }
              }
            }
          }  // if(use_gen_jets)


          if(!use_gen_jets){
            for(int i = 0; i < Njets; i++){
              ParticleList V_cand_Part;
              ConstRestFrameList V_cand_RF;
              const RestFrame& frame = BDT_CM.GetFrame(jetID_BDT[i]);
              const RestFrame& sibling = frame.GetSiblingFrame();
              const RestFrame& aunt = frame.GetProductionFrame().GetSiblingFrame();
              if(sibling.IsDecayFrame()) jet_BDT_singlet.push_back(i);
              if(int(lepID_BDT.size()) == Nleps){
                for(int j = 0; j < Nleps; j++){
                  if(sibling == BDT_CM.GetFrame(lepID_BDT[j])){
                    for(int k = i+1; k < Njets; k++){
                      if(aunt == BDT_CM.GetFrame(jetID_BDT[k])){
                        if(frame.GetProductionFrame().GetMass() - sibling.GetFourVector(frame.GetProductionFrame()).M() < 100.){
                          V_cand_Part.push_back(jets[i]);
                          V_cand_Part.push_back(jets[k]);
                          V_cand_RF.Add(frame);
                          V_cand_RF.Add(aunt);
                          V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                          if(use_lep_prong){
                            if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                              V_had_cands.push_back(cand);
                              lepprong++;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              for(int j = i+1; j < Njets; j++){
                if(sibling == BDT_CM.GetFrame(jetID_BDT[j])){
                  if(frame.GetProductionFrame().GetMass() < 100.){
                    V_cand_Part.clear();
                    V_cand_RF.Clear();
                    V_cand_Part.push_back(jets[i]);
                    V_cand_Part.push_back(jets[j]); // 2-prong "sibling"
                    V_cand_RF.Add(frame);
                    V_cand_RF.Add(sibling);
                    V_Cand cand = V_Cand(V_cand_Part,V_cand_RF);
                    if(use_prong2){
                      if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                        V_had_cands.push_back(cand);
                        prong2++;
                      }
                    }
                  }
                  //else if(frame.GetMass() < 100.){
                  //    V_cand_Part.push_back(jets[i]); // 1-prong "orphan"
                  //    V_cand_RF.Add(frame);
                  //  }
                  //for(int k = j+1; k < Njets; k++){
                    //if(aunt == BDT_CM.GetFrame(jetID_BDT[k]) && aunt.GetProductionFrame().GetMass() < 100.){
                    //    V_cand_Part.push_back(jets[k]); // 3-prong "aunt"
                    //    V_cand_RF.Add(aunt);
                    //    prong3++;
                    //    prong2--;
                    //}
                  //}
                }
                if(aunt == BDT_CM.GetFrame(jetID_BDT[j])){
                  bool lep_sibling = false;
                  if(int(lepID_BDT.size()) == Nleps){
                    for(int k = 0; k < Nleps; k++){
                      if(sibling == BDT_CM.GetFrame(lepID_BDT[k])){
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
                      if(BDT_CM.GetCosDecayAngle(cand.CandFrame()) > min_Cand_CosDecayAngle_CM && BDT_CM.GetCosDecayAngle(cand.CandFrame()) < max_Cand_CosDecayAngle_CM){
                        V_had_cands.push_back(cand);
                        auntprong++;
                      }
                    }
                  }
                }
              }
            }
          }



          int N_V_had = V_had_cands.size();
totcands += N_V_had;
          if(N_V_had < 1) continue;
          //hist_CandCount->SetBinContent(0,hist_CandCount->GetBinContent(0)+N_V_had);

          hist_CandCount->SetBinContent(0,hist_CandCount->GetBinContent(0)+Nhadbosons);

          for(int i = 0; i < N_V_had; i++){
            if(objID_BDT.size() >= 2){
              const RestFrame& CM_ChildA = BDT_CM.GetChildFrame(0);
              const RestFrame& CM_ChildB = BDT_CM.GetChildFrame(1);
              ConstRestFrameList CM_ChildA_Children = CM_ChildA.GetListFrames();
              ConstRestFrameList CM_ChildB_Children = CM_ChildB.GetListFrames();
              for(int r = 0; r < CM_ChildA_Children.GetN(); r++)
                if(V_had_cands[i].CandFrame() == CM_ChildA_Children[r])
                  V_had_cands[i].SetSide(kAside);
              for(int r = 0; r < CM_ChildB_Children.GetN(); r++)
                if(V_had_cands[i].CandFrame() == CM_ChildB_Children[r])
                  V_had_cands[i].SetSide(kBside);
            }
            bool matched = false; // perfect match
            bool unmatched = true; // both jets are radiative
            bool found24 = false;
            bool found23 = false;
            bool found6 = false;
            for(int j = 0; j < int(V_had_cands[i].size()); j++){
              for(int k = j+1; k < int(V_had_cands[i].size()); k++){
                //if(abs(V_had_cands[i][j].MomPDGID()) == 24 && V_had_cands[i][j].GenMomIndex() == V_had_cands[i][k].GenMomIndex()){
                if((abs(V_had_cands[i][j].MomPDGID()) == 23 || abs(V_had_cands[i][j].MomPDGID()) == 24) && V_had_cands[i][j].GenMomIndex() == V_had_cands[i][k].GenMomIndex()){
                  V_had_cands[i].SetMatch(kMatched);
                  hist_CandCount->SetBinContent(1,hist_CandCount->GetBinContent(1)+1);
                  matched = true;
                  unmatched = false;
                  break;
                }
              }
              if(!matched){
                if(V_had_cands[i][j].MomPDGID() == 24){
                  if(found24)
                    hist_CandCount->SetBinContent(2,hist_CandCount->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount->SetBinContent(3,hist_CandCount->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount->SetBinContent(4,hist_CandCount->GetBinContent(4)-1);
                  V_had_cands[i].SetMatch(kW);
                  hist_CandCount->SetBinContent(2,hist_CandCount->GetBinContent(2)+1);
                  found24 = true;
                  unmatched = false;
                }
                else if(V_had_cands[i][j].MomPDGID() == 23){
                  if(found24)
                    hist_CandCount->SetBinContent(2,hist_CandCount->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount->SetBinContent(3,hist_CandCount->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount->SetBinContent(4,hist_CandCount->GetBinContent(4)-1);
                  V_had_cands[i].SetMatch(kZ);
                  hist_CandCount->SetBinContent(3,hist_CandCount->GetBinContent(3)+1);
                  found23 = true;
                  unmatched = false;
                }
                else if(V_had_cands[i][j].MomPDGID() == 6){
                  if(found24)
                    hist_CandCount->SetBinContent(2,hist_CandCount->GetBinContent(2)-1);
                  if(found23)
                    hist_CandCount->SetBinContent(3,hist_CandCount->GetBinContent(3)-1);
                  if(found6)
                    hist_CandCount->SetBinContent(4,hist_CandCount->GetBinContent(4)-1);
                  V_had_cands[i].SetMatch(kB);
                  hist_CandCount->SetBinContent(4,hist_CandCount->GetBinContent(4)+1);
                  found6 = true;
                  unmatched = false;
                }
              } // if(!matched)
            } // for(int j = 0; j < int(V_had_cands[i].size()); j++)
            if(unmatched){
              V_had_cands[i].SetMatch(kUnmatched);
              hist_CandCount->SetBinContent(5,hist_CandCount->GetBinContent(5)+1);
            }
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

          for(int i = 0; i < N_V_had_combo; i++){
            if(V_had_cands_combo[i].Match() == kMatched){
              hist_combo_Matched_DeltaR_P->Fill(V_had_cands_combo[i].ProngDeltaR(), V_had_cands_combo[i].P(), hweight);
            }
            else if(V_had_cands_combo[i].Match() == kUnmatched){
              hist_combo_Unmatched_DeltaR_P->Fill(V_had_cands_combo[i].ProngDeltaR(), V_had_cands_combo[i].P(), hweight);
            }
            //else if(V_had_cands[i].Match() == kW || V_had_cands[i].Match() == kB){
            //}
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
                //if(abs(jets[i].MomPDGID()) == 24 && jets[i].GenMomIndex() == jets[j].GenMomIndex()){
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

          const RestFrame& Prod_Frame = BDT_CM;
          for(int i = 0; i < N_V_had; i++){
            if(V_had_cands[i].Match() == kMatched){
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
          hist_Njets_v_V_had_cand->Fill(Njets, N_V_had, hweight);
          hist_Njets_v_V_lep_cand->Fill(Njets, V_lep_cands.size(), hweight);
          hist_Nleps_v_V_lep_cand->Fill(Nleps, V_lep_cands.size(), hweight);
          hist_V_had_cand_v_V_lep_cand->Fill(N_V_had, V_lep_cands.size(), hweight);

          hist_MET->Fill(base->MET_pt, hweight);

          // RJR Tree
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
          if(!LAB.AnalyzeEvent()) cout << "Problem with RF Analyze Event \n";
          TVector3 vPISR = S.GetFourVector(CM).Vect();
          double PTISR = S.GetTransverseFourVector(CM).Vect().Mag();
          TVector3 vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();
          double RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
          hist_RISR_PTISR->Fill(RISR, PTISR, hweight);


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


          for(int i = 0; i < int(jet_BDT_singlet.size()); i++){
            RFKey key = COMB_J.AddLabFrameFourVector(jets[jet_BDT_singlet[i]]);
            jetID.push_back(key);
            objID.push_back(key);
          }
          for(int i = 0; i < Nleps; i++){
            RFKey key = COMB_L.AddLabFrameFourVector(leptons[i]);
            lepID.push_back(key);
            objID.push_back(key);
          }
          if(!LAB.AnalyzeEvent()) cout << "Problem with RF Cand Analyze Event \n";
          vPISR = S.GetFourVector(CM).Vect();
          PTISR = S.GetTransverseFourVector(CM).Vect().Mag();
          vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();
          RISR = fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
          hist_RISR_PTISR_Cands->Fill(RISR, PTISR, hweight);


          kept_sample_events++;
          kept_file_events++;
        }
        
        fstream file_event_info;
        file_event_info.open(file_event_info_filename, ios::app);
        if(f == 0)
          file_event_info << "sample: " << proc_Name << endl;
        file_event_info << f << "    " << kept_file_events << endl;
        file_event_info.close();
        delete base;
        delete chain;
      }
    }
    //cout << "Total for: " << title << " " << hists1[0]->Integral() << endl;
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11111111);

    fstream num_info;
    num_info.open(num_info_filename, ios::app);
    num_info << "sample: " << proc_Name << ", events: " << kept_sample_events << endl; 
    num_info << "quark eff: " << 100.-100.*quark_cut/tot_quark << endl;
    num_info << "genBosons: " << genHBosons << " total cands: " << totcands << " 2-prong cands: " << prong2 << " auntprong cands: " << auntprong << " lep-prong cands: " << lepprong << endl;
    num_info << "matched: " << hist_CandCount->GetBinContent(1) << " W: " << hist_CandCount->GetBinContent(2) << " Z: " << hist_CandCount->GetBinContent(3) << " B: " << hist_CandCount->GetBinContent(4) << " unmatched: " << hist_CandCount->GetBinContent(5) << endl;
    num_info << "matched A side: " << matchedAside << " matched B side: " << matchedBside << " unmatchedAside: " << unmatchedAside << " unmatchedBside: " << unmatchedBside << endl;
    num_info << endl;
    num_info.close();
  
    if(make_plots){
      for(int hist1 = 0; hist1 < int(hists1.size()); hist1++) Plot_Hist(hists1[hist1]);
      for(int hist2 = 0; hist2 < int(hists2.size()); hist2++) Plot_Hist(hists2[hist2]);
      for(int eff = 0; eff < int(effs.size()); eff++) Plot_Eff(effs[eff]);
      hists1.clear();
      hists2.clear();
      effs.clear();
  
      for(int i = 1; i <= hist_CandCount->GetNbinsX(); i++)
        hist_CandCount->SetBinContent(i,hist_CandCount->GetBinContent(i)/hist_CandCount->GetBinContent(0));
      Plot_Hist(hist_CandCount,false); 
      for(int i = 1; i <= hist_CandCount_Combo->GetNbinsX(); i++)
        hist_CandCount_Combo->SetBinContent(i,hist_CandCount_Combo->GetBinContent(i)/hist_CandCount_Combo->GetBinContent(0));
      Plot_Hist(hist_CandCount_Combo,false); 
      for(int i = 1; i <= hist_CandCount_Chi2->GetNbinsX(); i++)
        hist_CandCount_Chi2->SetBinContent(i,hist_CandCount_Chi2->GetBinContent(i)/hist_CandCount_Chi2->GetBinContent(0));
      Plot_Hist(hist_CandCount_Chi2,false); 
    }

  }
  std::cout << "FINISHED!" << std::endl;
  gApplication->Terminate(0);
}

void Plot_Eff(TEfficiency* e){
  //h->Scale(1./h->Integral());
  string title = e->GetName();
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+title).c_str(),("can_"+title).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
  can->SetLogz();
  can->Draw();
  can->cd();
  e->Draw("AP");
  gPad->Update();
  e->GetPaintedGraph()->GetXaxis()->CenterTitle();
  e->GetPaintedGraph()->GetXaxis()->SetTitleFont(42);
  e->GetPaintedGraph()->GetXaxis()->SetTitleSize(0.06);
  e->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.06);
  e->GetPaintedGraph()->GetXaxis()->SetLabelFont(42);
  e->GetPaintedGraph()->GetXaxis()->SetLabelSize(0.05);
  double xmin = e->GetTotalHistogram()->GetXaxis()->GetXmin();
  double xmax = e->GetTotalHistogram()->GetXaxis()->GetXmax();
  if(xmin < 0) xmin = xmin*1.1;
  else xmin = xmin*0.9;
  if(xmax > 0) xmax = xmax*1.1;
  else xmax = xmax*0.9;
  e->GetPaintedGraph()->GetXaxis()->SetRangeUser(xmin,xmax);
  //e->GetPaintedGraph()->GetXaxis()->SetTitle(g_Xname.c_str());
  e->GetPaintedGraph()->GetYaxis()->CenterTitle();
  e->GetPaintedGraph()->GetYaxis()->SetTitleFont(42);
  e->GetPaintedGraph()->GetYaxis()->SetTitleSize(0.06);
  e->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.12);
  e->GetPaintedGraph()->GetYaxis()->SetLabelFont(42);
  e->GetPaintedGraph()->GetYaxis()->SetLabelSize(0.05);
  //e->GetPaintedGraph()->GetYaxis()->SetTitle("N_{events} / 137 fb^{-1}");
  e->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.05);
  //e->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.9*h->GetMinimum(0.0),1.1*h->GetMaximum());

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
  l.DrawLatex(0.71,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = can->GetTitle();
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write();
  file->Close();
  delete e;
  delete can;
}

void Plot_Hist(TH1D* h, bool Scale){
  if(Scale) h->Scale(1./h->Integral());
  string title = h->GetName();
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+title).c_str(),("can_"+title).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
  can->SetLogz();
  can->Draw();
  can->cd();
  h->Draw("HIST");
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  //h->GetXaxis()->SetTitle(g_Xname.c_str());
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.12);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  //h->GetYaxis()->SetTitle("N_{events} / 137 fb^{-1}");
  //h->GetYaxis()->SetRangeUser(0.9*h->GetMinimum(0.0),1.1*h->GetMaximum());
  h->GetYaxis()->SetRangeUser(0.0,1.1*h->GetMaximum());

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
  l.DrawLatex(0.71,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = can->GetTitle();
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write();
  file->Close();
  delete h;
  delete can;
}

void Plot_Hist(TH2D* h){
  h->Scale(1./h->Integral());
  string title = h->GetName();
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+title).c_str(),("can_"+title).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
  can->SetLogz();
  can->Draw();
  can->cd();
  h->Draw("COLZ");
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleFont(42);
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.06);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.05);
  //h->GetXaxis()->SetTitle(g_Xname.c_str());
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(1.12);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelSize(0.05);
  //h->GetYaxis()->SetTitle(g_Yname.c_str());
  h->GetZaxis()->CenterTitle();
  h->GetZaxis()->SetTitleFont(42);
  h->GetZaxis()->SetTitleSize(0.055);
  h->GetZaxis()->SetTitleOffset(1.05);
  h->GetZaxis()->SetLabelFont(42);
  h->GetZaxis()->SetLabelSize(0.05);
  //h->GetZaxis()->SetTitle("N_{events} / 137 fb^{-1}");
  h->GetZaxis()->SetRangeUser(0.9*h->GetMinimum(0.0),1.1*h->GetMaximum());

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
  l.DrawLatex(0.71,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = can->GetTitle();
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write();
  file->Close();
  delete h;
  delete can;
}

string double_to_string(double val){
  std::stringstream strstr_val;
  strstr_val << std::fixed << std::setprecision(3) << val;
  string str_val = strstr_val.str();
  std::replace(str_val.begin(), str_val.end(), '.', 'p');
  return str_val;
}

double getWeightPerEvent(string input_dataset, string input_filetag, string EventCountFile){
  TFile* fout = new TFile(EventCountFile.c_str(),"READ");
  TBranch* b_dataset = nullptr;
  TBranch* b_filetag = nullptr;
  TBranch* b_Nweight = nullptr;
  TBranch* b_Nevent = nullptr;
  TBranch* b_MP = nullptr;
  TBranch* b_MC = nullptr;
  string* dataset = nullptr;
  string* filetag = nullptr;
  double Nweight = 0.0;
  double Nevent = 0.0;
  int MP = 0;
  int MC = 0;
  TTree* tree = nullptr;
  tree = (TTree*) fout->Get("EventCount");
  
  tree->SetMakeClass(1);
  tree->SetBranchAddress("Nevent", &Nevent,&b_Nevent);
  tree->SetBranchAddress("Nweight", &Nweight,&b_Nweight);
  tree->SetBranchAddress("filetag", &filetag,&b_filetag);
  tree->SetBranchAddress("dataset", &dataset,&b_dataset);
  tree->SetBranchAddress("MP", &MP,&b_MP);
  tree->SetBranchAddress("MC", &MC,&b_MC);

  double tot_Nevent = 0;
  double tot_Nweight = 0;

  for(int i = 0; i < tree->GetEntries(); i++)
  {
   tree->GetEntry(i);
   if(dataset->find(input_dataset) != std::string::npos && filetag->find(input_filetag) != std::string::npos)
   {
    tot_Nevent += Nevent;
    tot_Nweight += Nweight;
   }
  }
  fout->Close();
  delete fout;
  return tot_Nevent/tot_Nweight;
}


