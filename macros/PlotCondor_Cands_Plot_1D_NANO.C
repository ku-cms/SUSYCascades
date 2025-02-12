#include <iostream>
#include <string>
#include <fstream>

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
#include <TLine.h>
#include <TLorentzVector.h>
#include <TKey.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TError.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include "../include/SUSYNANOBase_slim.hh"
#include "../include/SampleTool.hh"
#include "../include/CategoryTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/Hadronic.hh"
#include "../include/FitReader.hh"
#include "../include/XsecTool.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

string g_PlotTitle;
string g_Label = "";
string plot_folder = "";
string output_root_file = "output.root";
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;

void hCount_to_text(string proc, TH1* h, string num_info_filename);
vector<string> get_dirs(const char *fname);
vector<TH1*> list_hists1(const char *fname, const char *dir);
vector<TH2*> list_hists2(const char *fname, const char *dir);
vector<TEfficiency*> list_effs(const char *fname, const char *dir);
void Plot_Hist(TH1* h, bool Scale=true, double Scale_Val = 0);
void Plot_Hist(TH2* h, bool Scale=true, double Scale_Val = 0);
void Plot_Eff(TEfficiency* e);

void PlotCondor_Cands_Plot_1D_NANO(){
  cout << "RUNNING PLOTTER..." << endl;
  bool condor = false;
  bool hadd = true;
  if(condor){
    cout << "Attempting to hold any running jobs to help prevent file corruption" << endl;
    gSystem->Exec("condor_hold $USER");
  }
  RestFrames::SetStyle();
  g_Xname = "";
  g_Xmin = 0.;
  g_Xmax = 500.; 
  g_NX = 128;
  g_PlotTitle = "";
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);

  bool Scale = false;
  bool Norm = true; // scale hists by 1/Integral
  string base_output = "output_Plot_1D_NANO";

  vector<string> extra_tags = {
    "_Sparticle2_RecoLep2_ge1JS_MET100_PTISR200_RISR0p5_SplitSJet",
    "_Sparticle2_RecoLep2_0JS_MET100_PTISR200_RISR0p5_SplitSJet",
  };

  std::map<string,double> ttbar_procs_2017 = {
    {"ttbar",1.},
  };
  std::map<string,double> TTZToQQ_procs_2017 = {
    {"TTZToQQ",1.},
  };
  std::map<string,double> DB_WZ_procs_2017 = {
    {"DB_WZ",1.},
  };
  std::map<string,double> DB_WW_procs_2017 = {
   {"DB_WW",1.},
  };
  std::map<string,double> ZDY_procs_2017 = {
   {"ZDY",1.},
  };
  std::map<string,double> WJets_procs_2017 = {
   {"WJets",1.},
  };
  std::map<string,double> TB_WWZ_procs_2017 = {
    {"TB_WWZ",1.},
  };
  std::map<string,double> TChiWZ_lowDM_procs_2017 = {
    {"TChiWZ_20",1.},
    {"TChiWZ_50",1.},
  };
  std::map<string,double> TChiWZ_medDM_procs_2017 = {
    {"TChiWZ_90",1.},
  };
  std::map<string,double> Cascades_20_procs_2023BPix = {
    {"Cascades_20",1.},
  };
  std::map<string,double> Cascades_10_procs_2023BPix = {
    {"Cascades_10",1.},
  };
  std::map<string,double> Cascades_5_procs_2023BPix = {
    {"Cascades_5",1.},
  };
  std::map<string,double> Cascades_10_Justin_procs_2023BPix = {
    {"Cascades_10_Justin",1.},
  };
  
  std::map<string,std::map<string,double>> WZ_dataset_2023BPix = {
    {"WZ_TuneCP5_13p6TeV_pythia8",DB_WZ_procs_2017},
  };
  std::map<string,std::map<string,double>> ttbar_dataset_2023BPix = {
    {"TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",ttbar_procs_2017},
    // {"TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",ttbar_procs_2017},
  };
  std::map<string,std::map<string,double>> ttbar_dataset_2017 = {
    {"TTJets_TuneCP5_13TeV-madgraphMLM-pythia8",ttbar_procs_2017},
  };
  std::map<string,std::map<string,double>> TTZToQQ_dataset_2017 = {
    {"TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8",TTZToQQ_procs_2017},
  };
  std::map<string,std::map<string,double>> DB_WZ_dataset_2017 = {
    {"WZ_TuneCP5_13TeV-pythia8",DB_WZ_procs_2017},
  };
  std::map<string,std::map<string,double>> DB_WW_dataset_2017 = {
    {"WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8",DB_WW_procs_2017},
  };
  std::map<string,std::map<string,double>> ZDY_dataset_2017 = {
    {"DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8",ZDY_procs_2017},
  };
  std::map<string,std::map<string,double>> WJets_dataset_2017 = {
    {"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",WJets_procs_2017},
  };
  std::map<string,std::map<string,double>> TB_WWZ_dataset_2017 = {
    {"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8",TB_WWZ_procs_2017},
  };
  std::map<string,std::map<string,double>> TChiWZ_lowDM_dataset_2017 = {
    {"TChiWZ_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8",TChiWZ_lowDM_procs_2017},
  };
  std::map<string,std::map<string,double>> TChiWZ_medDM_dataset_2017 = {
    {"SMS-TChiWZ_dM-60to90_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8",TChiWZ_medDM_procs_2017},
  };
  std::map<string,std::map<string,double>> Cascades_20_dataset_2023BPix = {
    {"SlepSnuCascade_MN1-220_MN2-260_MC1-240_TuneCP5_13p6TeV_madgraphMLM-pythia8",Cascades_20_procs_2023BPix},
  };
  std::map<string,std::map<string,double>> Cascades_10_dataset_2023BPix = {
    {"SlepSnuCascade_MN1-260_MN2-280_MC1-270_TuneCP5_13p6TeV_madgraphMLM-pythia8",Cascades_10_procs_2023BPix},
  };
  std::map<string,std::map<string,double>> Cascades_10_Justin_dataset_2023BPix = {
    {"SlepSnuCascade_MN1-260_MN2-280_MC1-270_TuneCP5_13p6TeV_madgraphMLM-pythia8_Justin",Cascades_10_Justin_procs_2023BPix},
  };
  std::map<string,std::map<string,double>> Cascades_5_dataset_2023BPix = {
    {"SlepSnuCascade_MN1-270_MN2-280_MC1-275_TuneCP5_13p6TeV_madgraphMLM-pythia8",Cascades_5_procs_2023BPix},
  };

  vector<std::map<string,std::map<string,double>>> bkg_datasets_2017 = {
    ttbar_dataset_2017,
    //TTZToQQ_dataset_2017,
    //DB_WZ_dataset_2017,
    DB_WW_dataset_2017,
    WJets_dataset_2017,
    //ZDY_dataset_2017,
    //TB_WWZ_dataset_2017,
  };
  vector<std::map<string,std::map<string,double>>> sms_datasets_2017 = {
    //TChiWZ_lowDM_dataset_2017,
    //TChiWZ_medDM_dataset_2017,
  };
  vector<std::map<string,std::map<string,double>>> sms_datasets_2023_BPix = {
    //Cascades_10_Justin_dataset_2023BPix,
    Cascades_20_dataset_2023BPix,
    Cascades_10_dataset_2023BPix,
    Cascades_5_dataset_2023BPix,
  };
  vector<std::map<string,std::map<string,double>>> bkg_datasets_2023_BPix = {
    //ttbar_dataset_2023BPix,
    //WZ_dataset_2023BPix,
  };

  std::map<string,vector<std::map<string,std::map<string,double>>>> filetags = {
    //{"Fall17_102X_SMS",sms_datasets_2017},
    {"Fall17_102X",bkg_datasets_2017},
    {"Summer23BPix_130X_SMS",sms_datasets_2023_BPix},
    //{"Summer23BPix_130X",bkg_datasets_2023_BPix},
  };

  for(int tag_i = 0; tag_i < int(extra_tags.size()); tag_i++){
    string extra_tag = extra_tags[tag_i];
    cout << "Processing: " << extra_tag << endl;
    string path_to_inputfiles = "/home/zflowers/CMSSW_13_3_1/src/SUSYCascades/SubmitCondor_Plot_1D_NANO"+extra_tag+"/";
    output_root_file = base_output+extra_tag+".root";
    TH2D* hist_EventCount = new TH2D("EventCount", "EventCount;", 12, 0., 12., 20, 0., 20.);
    for(std::map<string,vector<std::map<string,std::map<string,double>>>>::iterator iter0 = filetags.begin(); iter0 != filetags.end(); ++iter0){
      string filetag = iter0->first;
      vector<std::map<string,std::map<string,double>>> datasets = iter0->second;
      for(int iter1 = 0; iter1 < int(datasets.size()); iter1++){
        for(std::map<string,std::map<string,double>>::iterator iter2 = datasets[iter1].begin(); iter2 != datasets[iter1].end(); ++iter2){
          string dataset_name = iter2->first;
          cout << "  Running on " << dataset_name << endl;
          std::map<string,double> procs = iter2->second;
          for(std::map<string,double>::iterator iter3 = procs.begin(); iter3 != procs.end(); ++ iter3){
            string proc = iter3->first;
            double SF = iter3->second;
            string hadd_cmd = "hadd -k -f -v 0 "+path_to_inputfiles+base_output+"_"+filetag+"_"+dataset_name+"_"+proc+".root "+path_to_inputfiles+"root/"+base_output+"_"+filetag+"_"+dataset_name+"_"+proc+"_*.root";
            if(hadd) gSystem->Exec(hadd_cmd.c_str());
            string inputRootFileName = path_to_inputfiles+base_output+"_"+filetag+"_"+dataset_name+"_"+proc+".root";
            vector<string> dirs = get_dirs(inputRootFileName.c_str());
            for(int d = 0; d < int(dirs.size()); d++){
              string dir = dirs[d];
              plot_folder = dir+"_"+extra_tag+"/";
              gSystem->Exec(("mkdir -p "+plot_folder).c_str());
              vector<TH1*> hists1 = list_hists1(inputRootFileName.c_str(), dir.c_str());
              vector<TH2*> hists2 = list_hists2(inputRootFileName.c_str(), dir.c_str());
              vector<TEfficiency*> effs = list_effs(inputRootFileName.c_str(), dir.c_str());
              TFile* inputRootFile = new TFile(inputRootFileName.c_str(), "READ");
              string count_hist = dir+"/"+proc+"_MET";
              TH1* hist_MET = (TH1*)inputRootFile->Get(count_hist.c_str())->Clone();
              double events = hist_MET->Integral();
              inputRootFile->Close();
              delete inputRootFile;
              g_PlotTitle = proc+"_"+filetag;
              //cout << "Plotting 1D hists... " << endl;
              string num_info_filename = plot_folder+"/num_info";
              for(int hist1 = 0; hist1 < int(hists1.size()); hist1++){
                string h_name = hists1[hist1]->GetName();
                if(h_name.find("EventCount") != std::string::npos) {hist_EventCount->Add(hists1[hist1]); continue;}
                if(h_name.find("_Count") != std::string::npos)
                  hCount_to_text(proc, hists1[hist1], num_info_filename);
                else
                  Plot_Hist(hists1[hist1], Norm, events);
              } // for(int hist1 = 0; hist1 < int(hists1.size()); hist1++){
              //cout << "Plotting 2D hists... " << endl;
              for(int hist2 = 0; hist2 < int(hists2.size()); hist2++) Plot_Hist(hists2[hist2], Norm, events);
              //cout << "Plotting effs... " << endl;
              for(int eff = 0; eff < int(effs.size()); eff++) Plot_Eff(effs[eff]);
            } // for(int d = 0; d < int(dirs.size()); d++){
          } // for(std::map<string,double>::iterator iter3 = procs.begin(); iter3 != procs.end(); ++ iter3){
        } // for(vector<std::map<string,std::map<string,std::map<string,double>>>>::iterator iter2 = datasets[iter1].begin(); iter2 != datasets[iter1].end(); ++iter2){
      } // for(int iter1 = 0; iter1 < int(datasets.size()); iter1++){
    } //for(std::map<string,vector<std::map<string,std::map<string,std::map<string,double>>>>>::iterator iter0 = filetags.begin(); iter0 != filetags.end(); ++iter0){

    int min_Nleps = 2;
    int max_Nleps = 4; // >=
    int min_Njets_S = 0;
    int max_Njets_S = 2; // >=
    int bin = 1;
    for(int l = min_Nleps; l <= max_Nleps; l++){
      for(int j = min_Njets_S; j <= max_Njets_S; j++){
        hist_EventCount->GetXaxis()->SetBinLabel(bin, (std::to_string(l)+"L "+std::to_string(j)+"J").c_str());
        bin++;
      }
    }
    hist_EventCount->GetYaxis()->SetBinLabel(1, "ttbar");
    hist_EventCount->GetYaxis()->SetBinLabel(2, "WW");
    hist_EventCount->GetYaxis()->SetBinLabel(3, "WJets");
    hist_EventCount->GetYaxis()->SetBinLabel(4, ""); // DY
    int DM = 20;
    int lastbkg = 4;
    for(int i = 0; i < 3; i++){
      hist_EventCount->GetYaxis()->SetBinLabel(lastbkg+1, ("Cascades #tilde{#nu} #tilde{#nu} "+std::to_string(DM)).c_str());
      hist_EventCount->GetYaxis()->SetBinLabel(lastbkg+2, ("Cascades #tilde{#it{l}} #tilde{#nu} "+std::to_string(DM)).c_str());
      hist_EventCount->GetYaxis()->SetBinLabel(lastbkg+3, ("Cascades #tilde{#it{l}} #tilde{#it{l}} "+std::to_string(DM)).c_str());
      DM /= 2;
      lastbkg += 3;
    }
    for(int i = 1; i < hist_EventCount->GetNbinsX(); i++)
      for(int j = 1; j < hist_EventCount->GetNbinsY(); j++)
        if(hist_EventCount->GetBinContent(i,j) > 0)
          hist_EventCount->SetBinContent(i,j, hist_EventCount->GetBinContent(i,j)/hist_EventCount->GetBinContent(0,j));
    Plot_Hist(hist_EventCount, false, 1.);

  } // for(int tag_i = 0; tag_i < int(extra_tags.size()); tag_i++)
  if(condor){
    cout << "Releasing previously held jobs" << endl;
    gSystem->Exec("condor_release $USER");
  }
  cout << "FINISHED!" << endl;
  gApplication->Terminate(0);
} // void PlotCondor_Cands_Plot_1D_NANO(){

vector<string> get_dirs(const char *fname)
{   
  vector<string> vect_dirnames;
  TFile *f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie())
  {   
    cout << "Unable to open " << fname << " for reading..." << endl;
    return vect_dirnames;
  }           
  TKey *key = nullptr;
  TIter next((TList *)f->GetListOfKeys());
  while((key = (TKey *)next()))
    vect_dirnames.push_back(key->GetName());
  return vect_dirnames;
}

vector<TH1*> list_hists1(const char *fname, const char *dir)
{   
  std::vector<TH1*> vect_hist;
  TKey *key;
  TFile *f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){   
    cout << "Unable to open " << fname << " for reading..." << endl;
    return vect_hist;
  }           
  TDirectoryFile* folder = nullptr;
  f->GetObject(dir,folder);
  TIter next((TList *)folder->GetListOfKeys());
  while((key = (TKey *)next())){
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if(cl->InheritsFrom("TH1")){
      string h_name = key->GetName();
      string str_dir = dir;
      h_name = str_dir+"/"+h_name;
      vect_hist.push_back((TH1*)f->Get(h_name.c_str()));
    }
  }
  return vect_hist;
}

vector<TH2*> list_hists2(const char *fname, const char *dir)
{   
  std::vector<TH2*> vect_hist;
  TKey *key;
  TFile *f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){   
    cout << "Unable to open " << fname << " for reading..." << endl;
    return vect_hist;
  }           
  TDirectoryFile* folder = nullptr;
  f->GetObject(dir,folder);
  TIter next((TList *)folder->GetListOfKeys());
  while((key = (TKey *)next())){
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if(cl->InheritsFrom("TH2")){
      string h_name = key->GetName();
      string str_dir = dir;
      h_name = str_dir+"/"+h_name;
      vect_hist.push_back((TH2*)f->Get(h_name.c_str()));
    }
  }
  return vect_hist;
}

vector<TEfficiency*> list_effs(const char *fname, const char *dir)
{   
  std::vector<TEfficiency*> vect_eff;
  TKey *key;
  TFile *f = TFile::Open(fname, "READ");
  if(!f || f->IsZombie()){   
    cout << "Unable to open " << fname << " for reading..." << endl;
    return vect_eff;
  }           
  TDirectoryFile* folder = nullptr;
  f->GetObject(dir,folder);
  TIter next((TList *)folder->GetListOfKeys());
  while((key = (TKey *)next())){
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if(cl->InheritsFrom("TEfficiency")){
      string e_name = key->GetName();
      string str_dir = dir;
      e_name = str_dir+"/"+e_name;
      vect_eff.push_back((TEfficiency*)f->Get(e_name.c_str()));
    }
  }
  return vect_eff;
}

void Plot_Eff(TEfficiency* e){
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
  l.DrawLatex(0.57,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  //l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete e;
  delete can;
}

void Plot_Hist(TH1* h, bool Scale, double Scale_Val){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  string title = h->GetName();
  if(title.find("CandCount") != std::string::npos){
    for(int i = 1; i <= h->GetNbinsX(); i++)
      h->SetBinContent(i,h->GetBinContent(i)/h->GetBinContent(0));
  }
  if(Scale && title.find("Count") == std::string::npos) h->Scale(1./Scale_Val);
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+title).c_str(),("can_"+title).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
  can->Draw();
  can->cd();
  if(title.find("EventCount") != std::string::npos) h->Draw("TEXT");
  else h->Draw("HIST");
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
  l.DrawLatex(0.57,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  //l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete h;
  delete can;
}

void Plot_Hist(TH2* h, bool Scale, double Scale_Val){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  h->Scale(1./Scale_Val);
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
  l.DrawLatex(0.57,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01,0.943,"#bf{CMS} Simulation Preliminary");

  l.SetTextSize(0.045);
  l.SetTextFont(42);
  //l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  can->SaveAs((plot_folder+pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete h;
  delete can;
}

void hCount_to_text(string proc, TH1* h, string num_info_filename){
  fstream num_info;
  string h_name = h->GetName();
  if(h_name.find("ISR_BDT") != std::string::npos)
    num_info_filename += "_ISR_BDT.txt";
  else if(h_name.find("BDT_ISR_lep") != std::string::npos)
    num_info_filename += "_BDT_ISR_lep.txt";
  else if(h_name.find("BDT_ISR_singlet") != std::string::npos)
    num_info_filename += "_BDT_ISR_singlet.txt";
  else if(h_name.find("BDT_ISR") != std::string::npos)
    num_info_filename += "_BDT_ISR.txt";
  else
    num_info_filename += ".txt";
  num_info.open(num_info_filename, ios::app);
  num_info << "process: " << proc << endl;
  num_info << std::fixed << "total events processed: " << int(h->GetBinContent(0)) << endl;
  num_info << std::fixed << "kept events: " << int(h->GetBinContent(1)) << endl;
  num_info << std::fixed << "quark eff: " << 100.-100.*int(h->GetBinContent(2))/int(h->GetBinContent(3)) << endl;
  num_info << std::fixed << "genBosons: " << int(h->GetBinContent(4)) << " total cands: " << int(h->GetBinContent(5)) << " 2-prong cands: " << int(h->GetBinContent(6)) << " auntprong cands: " << int(h->GetBinContent(7)) << " lep-prong cands: " << int(h->GetBinContent(8)) << endl;
  num_info << std::fixed << "matched: " << int(h->GetBinContent(9)) << " W: " << int(h->GetBinContent(10)) << " Z: " << int(h->GetBinContent(11)) << " B: " << int(h->GetBinContent(12)) << " unmatched: " << int(h->GetBinContent(13)) << endl;
  num_info << std::fixed << "matched A side: " << int(h->GetBinContent(14)) << " matched B side: " << int(h->GetBinContent(15)) << " unmatchedAside: " << int(h->GetBinContent(16)) << " unmatchedBside: " << int(h->GetBinContent(17)) << endl;
  num_info << "Ratios: " << endl;
  num_info << "Total cands/bosons: " << h->GetBinContent(5)/h->GetBinContent(4) << endl;
  num_info << "Matched/bosons: " << h->GetBinContent(9)/h->GetBinContent(4) << endl;
  num_info << "Unmatched/bosons: " << h->GetBinContent(13)/h->GetBinContent(4) << endl;
  num_info << "Matched/Unmatched: " << h->GetBinContent(9)/h->GetBinContent(13) << endl;
  num_info << "Matched A/Matched B: " << h->GetBinContent(14)/h->GetBinContent(15) << endl;
  num_info << "Unmatched A/Unmatched B: " << h->GetBinContent(16)/h->GetBinContent(17) << endl;
  num_info << "Matched A/Unmatched A: " << h->GetBinContent(14)/h->GetBinContent(16) << endl;
  num_info << "Matched B/Unmatched B: " << h->GetBinContent(15)/h->GetBinContent(17) << endl;
  num_info << endl;
  num_info.close();
}
