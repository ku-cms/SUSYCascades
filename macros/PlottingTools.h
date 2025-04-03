#include <iostream>
#include <string>

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
#include <TH1D.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TError.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TApplication.h>

#include "../include/ReducedBase_V2.hh"
#include "../include/SampleTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/FitPlotter.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

// declare global vars
string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX = 128.;
string g_Yname;
double g_Ymin;
double g_Ymax;
double g_NY;
string g_Label;
// root file to store output of plots
string output_root_file = "output_Plot_";
bool SavePDF = true; // whether or not to save pdfs of plots
string folder_name = "";

FitPlotter FP("","","");

// Function to sort histograms and reorder the map accordingly
void SortHistogramsAndProcesses(std::vector<TH1*>& histograms, std::vector<std::pair<std::string, ProcessList>>& vec_samples){
  // Step 1: Create a vector of indices for sorting
  std::vector<size_t> indices(histograms.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    indices[i] = i;
  }

  // Step 2: Sort indices based on histogram integral values (descending order)
  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
    return histograms[a]->Integral() > histograms[b]->Integral();
  });

  // Step 3: Create sorted versions of histograms and key-value pairs
  std::vector<TH1*> sorted_histograms(histograms.size());
  std::vector<std::pair<std::string, ProcessList>> sorted_vec_samples(vec_samples.size());

  for (size_t i = 0; i < indices.size(); ++i) {
    sorted_histograms[i] = histograms[indices[i]];
    sorted_vec_samples[i] = vec_samples[indices[i]]; // Keep key-value mapping consistent
  }

  // Step 4: Update original vectors with sorted versions
  histograms = std::move(sorted_histograms);
  vec_samples = std::move(sorted_vec_samples);
}

void GetMinMaxIntegral(const std::vector<TH1*>& histograms, double& minIntegral, double& maxIntegral){
  if (histograms.empty()) return; // Handle empty vector case

  auto [minHist, maxHist] = std::minmax_element(histograms.begin(), histograms.end(),
      [](TH1* a, TH1* b) {
        return a->Integral() < b->Integral();  // Compare by integral value
      });

  minIntegral = (*minHist)->GetMinimum();  // Minimum y-value of the lowest integral histogram
  maxIntegral = (*maxHist)->GetMaximum();  // Maximum y-value of the highest integral histogram
}

void SetMinimumBinContent(TH1* hist, double min_value) {
  if (!hist) return; // Safety check for null pointer

  int nBins = hist->GetNbinsX();
  for (int bin = 1; bin <= nBins; ++bin) // Loop over bins (excluding underflow/overflow)
    if (hist->GetBinContent(bin) < min_value)
      hist->SetBinContent(bin, min_value);
}


void Plot_Stack(vector<TH1*>& vect_h, std::vector<std::pair<std::string, ProcessList>> vec_samples, double signal_boost){
  TH1D* h_BKG = nullptr;
  TH1D* h_DATA = nullptr;
  bool isBKG = false;
  int Nsample = vec_samples.size();
  if (Nsample != int(vect_h.size())){
    std::cout << "hist stack size does not match sample stack size!" << std::endl;
    return;
  }

  SortHistogramsAndProcesses(vect_h, vec_samples);

  int index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    SetMinimumBinContent(vect_h[index],1.e-6);
    if(p->second[0].Type() == kBkg){
      if(!isBKG){
        h_BKG = (TH1D*) vect_h[index]->Clone("TOT_BKG");
        isBKG = true;
      } else {
        for(int k = 0; k < index; k++){
          vect_h[k]->Add(vect_h[index]);
        }
        h_BKG->Add(vect_h[index]);
      }
    } // if(p->second[0].Type() == kBkg)
    if(p->second[0].Type() == kData){
      h_DATA = (TH1D*) vect_h[index]->Clone("TOT_DATA");
    }
    index++;
  }
  
  double h_min, h_max;
  GetMinMaxIntegral(vect_h, h_min, h_max);
  if(h_min <= 0.) h_min = 1.e-1;
  double fmax = -1.;
  int imax = -1;
  for(int i = 0; i < Nsample; i++){
    if(vect_h[i]->GetMaximum() > fmax){
      fmax = vect_h[i]->GetMaximum();
      imax = i;
    }
  }
  float width = vect_h[0]->GetBinWidth(1);
  char *yaxis = new char[100];
  sprintf(yaxis,"Events / %.2f bin", width);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+g_PlotTitle).c_str(),("can_"+g_PlotTitle).c_str(),1200,600); // 1200, 700 is default

  double hlo = 0.09;
  double hhi = 0.22;
  double hbo = 0.27;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridx();
  can->SetGridy();
  can->SetLogy();
  can->Draw();
  can->cd();

  TPad* pad = new TPad("pad","pad",0,0.32,1.,1.);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetBottomMargin(0.02);
  if(h_BKG != nullptr && h_DATA != nullptr){
    if(h_BKG->GetEntries() > 0 && h_DATA->GetEntries() > 0){
      pad->Draw();
      pad->cd();
      pad->Update();
      can->Update();
    }
  }

  vect_h[imax]->Draw("hist");
  vect_h[imax]->GetXaxis()->CenterTitle();
  vect_h[imax]->GetXaxis()->SetTitleFont(42);
  vect_h[imax]->GetXaxis()->SetTitleSize(0.05);
  vect_h[imax]->GetXaxis()->SetTitleOffset(1.0);
  vect_h[imax]->GetXaxis()->SetLabelFont(42);
  vect_h[imax]->GetXaxis()->SetLabelSize(0.04);
  vect_h[imax]->GetXaxis()->SetTickSize(0.);
  vect_h[imax]->GetYaxis()->CenterTitle();
  vect_h[imax]->GetYaxis()->SetTitleFont(42);
  vect_h[imax]->GetYaxis()->SetTitleSize(0.04);
  vect_h[imax]->GetYaxis()->SetTitleOffset(0.9);
  vect_h[imax]->GetYaxis()->SetLabelFont(42);
  vect_h[imax]->GetYaxis()->SetLabelSize(0.035);
  vect_h[imax]->GetYaxis()->SetRangeUser(0.9*h_min,1.1*h_max);

  index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    if(p->second[0].Type() == kBkg){
      if(vect_h[index]->GetEntries() == 0) {index++; continue;}
      vect_h[index]->SetLineColor(kBlack);
      vect_h[index]->SetLineWidth(1.0);
      vect_h[index]->SetFillColor(FP.getColor(p->first));
      vect_h[index]->SetFillStyle(1001);
      vect_h[index]->Draw("SAME HIST");
    }
    index++;
  }

  if(isBKG){
    h_BKG->SetLineWidth(3.0);
    h_BKG->SetLineColor(kRed);
    h_BKG->SetMarkerSize(0);
    h_BKG->Draw("SAME HIST");
  }
  
  index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    if(p->second[0].Type() == kSig){
      vect_h[index]->SetLineWidth(3.0);
      vect_h[index]->SetMarkerSize(0.);
      vect_h[index]->SetMarkerColor(kBlack);
      vect_h[index]->SetLineStyle(7);
      vect_h[index]->SetLineColor(FP.getColor(p->first));
      vect_h[index]->Draw("SAME HIST");
    }
    if(p->second[0].Type() == kData){
      vect_h[index]->SetLineWidth(3.0);
      vect_h[index]->SetMarkerSize(0.);
      vect_h[index]->SetMarkerColor(kBlack);
      vect_h[index]->SetLineStyle(7);
      vect_h[index]->SetLineColor(FP.getColor(p->first));
      vect_h[index]->Draw("SAME");
    }
    index++;
  }

  TLegend* leg = new TLegend(1.-hhi+0.01, 1.- (vect_h.size()+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(132);
  leg->SetTextSize(0.042);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  if(isBKG) leg->AddEntry(h_BKG, "SM total");
  index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    if(vect_h[index]->GetEntries() == 0) {index++; continue;}
    if(p->second[0].Type() == kBkg)
      leg->AddEntry(vect_h[index],FP.getTitle(p->first).c_str(),"F");
    else if(p->second[0].Type() == kSig && signal_boost != 1.)
      leg->AddEntry(vect_h[index],(FP.getTitle(p->first)+" * "+std::to_string(int(signal_boost))).c_str(),"F");
    else
      leg->AddEntry(vect_h[index],FP.getTitle(p->first).c_str());
    index++;
  }
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->Draw("SAME");

  TLatex l;
  l.SetTextFont(132);
  l.SetNDC();
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  l.DrawLatex(0.55,0.943,g_Label.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.09,0.943,"#bf{#it{CMS}} Internal 13 TeV Simulation");
  l.SetTextSize(0.05);
  l.SetTextFont(132);
  //string s_lumi = "#scale[0.6]{#int} #it{L dt} = "+to_string(int(g_Lumi))+" fb^{-1}";
  //l.DrawLatex(0.43,0.79,s_lumi.c_str());	
  l.SetTextSize(0.045);
  l.SetTextFont(42);

  TPad *pad_res = new TPad("pad_res","pad_res",0.,0.03,1.,0.31);
  pad_res->SetGridx(); 
  pad_res->SetGridy();
  pad_res->SetTopMargin(0.05);
  pad_res->SetBottomMargin(0.2);
  if(h_BKG != nullptr && h_DATA != nullptr){
    if(h_BKG->GetEntries() > 0 && h_DATA->GetEntries() > 0){
      vect_h[imax]->GetXaxis()->SetLabelOffset(0.05);
      TH1D* h_res = (TH1D*) h_BKG->Clone("res");
      h_res->Divide(h_DATA);
      can->Update();
      can->cd();
      pad_res->Draw();
      pad_res->cd();
      pad_res->Update();
      can->Update();
      h_res->Draw("");
      h_res->GetYaxis()->SetTitle("TOT BKG / DATA");
      h_res->GetYaxis()->SetRangeUser(0.9*h_res->GetMinimum(),1.1*h_res->GetMaximum());
      //h_res->GetYaxis()->SetRangeUser(0., 8.);
      h_res->GetXaxis()->CenterTitle();
      h_res->GetXaxis()->SetTitleFont(132);
      h_res->GetXaxis()->SetTitleSize(0.08);
      h_res->GetXaxis()->SetTitleOffset(1.02);
      h_res->GetXaxis()->SetLabelFont(132);
      h_res->GetXaxis()->SetLabelSize(0.08);
      h_res->GetXaxis()->SetTitle(g_Xname.c_str());
      h_res->GetYaxis()->CenterTitle();
      h_res->GetYaxis()->SetTitleFont(132);
      h_res->GetYaxis()->SetTitleSize(0.06);
      h_res->GetYaxis()->SetTitleOffset(1.02);
      h_res->GetYaxis()->SetLabelFont(132);
      h_res->GetYaxis()->SetLabelSize(0.08);
      pad_res->Modified();
      pad_res->Update();
      can->Update();
    }
  }
  
  string pdf_title = folder_name+"/"+g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  if(SavePDF)
    can->SaveAs((pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  for(int i = 0; i < Nsample; i++){
    delete vect_h[i];
  }
  delete[] yaxis;
  delete leg;
  delete can;
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
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = folder_name+"/"+g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  if(SavePDF)
    can->SaveAs((pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete e;
  delete can;
}

void Plot_Hist(TH1* h, bool Scale, double Scale_Val, double signal_boost, bool IsSMS){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  if(IsSMS) Scale_Val *= signal_boost;
  string title = h->GetName();
  // do not scale data
  if(title.find("data") == std::string::npos) h->Scale(Scale_Val);
  if(IsSMS) Scale_Val /= signal_boost;
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+title).c_str(),("can_"+title).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
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
  h->GetYaxis()->SetTitle(("N_{events} / "+std::to_string(int(Scale_Val))+" fb^{-1}").c_str());
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
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = folder_name+"/"+g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  if(SavePDF)
    can->SaveAs((pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete can;
}

void Plot_Hist(TH2* h, bool Scale, double Scale_Val, double signal_boost, bool IsSMS){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  if(IsSMS) Scale_Val *= signal_boost;
  string title = h->GetName();
  // do not scale data
  if(title.find("data") == std::string::npos) h->Scale(Scale_Val);
  if(IsSMS) Scale_Val /= signal_boost;
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
  h->GetZaxis()->SetTitle(("N_{events} / "+std::to_string(int(Scale_Val))+" fb^{-1}").c_str());
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
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = folder_name+"/"+g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  if(SavePDF)
    can->SaveAs((pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete h;
  delete can;
}

void Plot_Ratio(TH1* h_num, TH1* h_den){
  h_num->Divide(h_den);  
  float width = h_num->GetBinWidth(1);
  char *yaxis = new char[100];
  sprintf(yaxis,"Events / %f", width);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",700.,600);

  can->SetLeftMargin(0.13);
  can->SetRightMargin(0.04);
  can->SetBottomMargin(0.16);
  can->SetTopMargin(0.085);
  can->SetGridx();
  can->SetGridy();
  can->Draw();
  can->cd();

  h_num->Draw("hist");
  h_num->GetXaxis()->CenterTitle();
  h_num->GetXaxis()->SetTitleFont(42);
  h_num->GetXaxis()->SetTitleSize(0.06);
  h_num->GetXaxis()->SetTitleOffset(1.06);
  h_num->GetXaxis()->SetLabelFont(42);
  h_num->GetXaxis()->SetLabelSize(0.05);
  h_num->GetXaxis()->SetTitle(g_Xname.c_str());
  h_num->GetYaxis()->CenterTitle();
  h_num->GetYaxis()->SetTitleFont(42);
  h_num->GetYaxis()->SetTitleSize(0.06);
  h_num->GetYaxis()->SetTitleOffset(1.1);
  h_num->GetYaxis()->SetLabelFont(42);
  h_num->GetYaxis()->SetLabelSize(0.05);
  h_num->GetYaxis()->SetTitle("a. u.");
  h_num->GetYaxis()->SetRangeUser(0., h_num->GetMaximum()*1.1);
  //hist->GetYaxis()->SetTitle(yaxis);
  //hist->GetYaxis()->SetTitle("N_{evt} / fb^{-1}");
  h_num->Draw("hist SAME");

  TLatex l;
  l.SetTextFont(132);
  l.SetNDC();
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.6,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.135,0.943,"#bf{CMS} Simulation Preliminary");

  string pdf_title = folder_name+"/"+g_PlotTitle+"_";
  pdf_title += can->GetTitle();
  gErrorIgnoreLevel = 1001;
  if(SavePDF)
    can->SaveAs((pdf_title+".pdf").c_str());
  gErrorIgnoreLevel = 0;
  TFile* file = new TFile(output_root_file.c_str(),"UPDATE");
  can->Write(0,TObject::kWriteDelete);
  file->Close();
  delete[] yaxis;
  delete can;
}
