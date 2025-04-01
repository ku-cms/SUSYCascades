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

#include "RestFrames/RestFrames.hh"
#include "../include/ReducedBase_V2.hh"
#include "../include/SampleSet.hh"

using namespace std;

using namespace RestFrames;

vector<SampleSet*> g_Samples;

string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;
string g_Yname;
double g_Ymin;
double g_Ymax;
double g_NY;

void Plot_1D(){
  RestFrames::SetStyle();

  g_PlotTitle = "di-lepton selection";

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";

  int BKG_SKIP = 1;
  
  SampleSet ttX;
  ttX.SetBkg(true);
  ttX.SetTitle("t #bar{t} + jets");
  ttX.SetColor(kAzure+1);
  ttX.AddFile(NtuplePath+"Summer23BPix_130X/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Summer23BPix_130X.root");
  ttX.SetSkip(BKG_SKIP);
  g_Samples.push_back(&ttX);
  
  //SampleSet SIG1;
  //SIG1.SetBkg(false);
  //SIG1.SetTitle("#scale[0.8]{Cascades; m_{#chi^{0}_{2}} = 280, m_{#chi^{0}_{1}} = 220 GeV");
  //SIG1.SetTreeName("SMS_300_220");
  //SIG1.SetColor(kMagenta);
  //SIG1.AddFile(NtuplePath+"Summer23BPix_130X/");
  //SIG1.SetSkip(1);
  //g_Samples.push_back(&SIG1);

  int Nsample = g_Samples.size();
 
  g_Xname = "RISR";
  g_Xmin = 0.;
  g_Xmax = 1.;
  g_NX = 128;

  TH1D* hist[Nsample];
  for(int i = 0; i < Nsample; i++)
    hist[i] = new TH1D(("h"+to_string(i)).c_str(),
		       ("h"+to_string(i)).c_str(),
		       g_NX,g_Xmin,g_Xmax);

  for(int s = 0; s < Nsample; s++){
    
    int Nfile = g_Samples[s]->GetNFile();
    cout << "Processing " << Nfile << " files for sample " << g_Samples[s]->GetTitle() << endl;
    for(int f = 0; f < Nfile; f++){
      cout << "   Processing file " << g_Samples[s]->GetFile(f) << " w/ tree " << g_Samples[s]->GetTreeName() << endl;
    
      TChain* chain = new TChain(g_Samples[s]->GetTreeName().c_str());
      chain->Add(g_Samples[s]->GetFile(f).c_str());

      ReducedBase_V2* base = new ReducedBase_V2(chain);

      int Nentry = base->fChain->GetEntries();

      int SKIP = g_Samples[s]->GetSkip();
     
      for(int e = 0; e < Nentry; e += SKIP){
	base->GetEntry(e);
	if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
	  cout << "      event " << e << " | " << Nentry << endl;
	
	///////////// 2 lepton
	if(base->Nlep != 2)
	  continue;

	if(base->ID_lep->at(0) < 3 ||
	   base->ID_lep->at(1) < 3)
	  continue;

	if(base->MiniIso_lep->at(0)*base->PT_lep->at(0) > 4. ||
	   base->MiniIso_lep->at(1)*base->PT_lep->at(1) > 4.)
	  continue;

	if(base->SIP3D_lep->at(0) > 8 ||
	   base->SIP3D_lep->at(1) > 8)
	continue;
	 
	if(base->MET < 100)
	  continue;
	
	if(base->PTISR < 200)
	  continue;

	hist[s]->Fill(base->RISR, base->weight*double(SKIP));
      }
      
      delete base;
      delete chain;
    }
  }
      

  
  double max = -1.;
  int imax = -1;
  for(int i = 0; i < Nsample; i++){
    cout << g_Samples[i]->GetTitle().c_str() << " " << hist[i]->Integral()*137 << " total" << endl;
    hist[i]->Scale(1./hist[i]->Integral());
    if(hist[i]->GetMaximum() > max){
      max = hist[i]->GetMaximum();
      imax = i;
    }
  }
  
  float width = hist[0]->GetBinWidth(1);
  char *yaxis = new char[100];
  sprintf(yaxis,"Events / %f", width);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",600.,500);

  can->SetLeftMargin(0.13);
  can->SetRightMargin(0.04);
  can->SetBottomMargin(0.16);
  can->SetTopMargin(0.085);
  can->SetGridx();
  can->SetGridy();
  can->Draw();
  can->cd();
  hist[imax]->Draw("hist");
  hist[imax]->GetXaxis()->CenterTitle();
  hist[imax]->GetXaxis()->SetTitleFont(42);
  hist[imax]->GetXaxis()->SetTitleSize(0.06);
  hist[imax]->GetXaxis()->SetTitleOffset(1.06);
  hist[imax]->GetXaxis()->SetLabelFont(42);
  hist[imax]->GetXaxis()->SetLabelSize(0.05);
  hist[imax]->GetXaxis()->SetTitle(g_Xname.c_str());
  hist[imax]->GetYaxis()->CenterTitle();
  hist[imax]->GetYaxis()->SetTitleFont(42);
  hist[imax]->GetYaxis()->SetTitleSize(0.06);
  hist[imax]->GetYaxis()->SetTitleOffset(1.1);
  hist[imax]->GetYaxis()->SetLabelFont(42);
  hist[imax]->GetYaxis()->SetLabelSize(0.05);
  hist[imax]->GetYaxis()->SetTitle("a. u.");
  hist[imax]->GetYaxis()->SetRangeUser(0., hist[imax]->GetMaximum()*1.1);
  int Ntype[3];

  int mycolor[8];
  mycolor[0] = kBlue+2;
  mycolor[1] = kGreen+3;
  mycolor[2] = kRed+1;
  mycolor[3] = kYellow+2;
  mycolor[4] = kMagenta+1;
  mycolor[5] = kMagenta+2;
  mycolor[6] = kCyan+2;
  mycolor[7] = kCyan+3;
  
  Ntype[0] = 0;
  for(int i = Nsample-1; i >= 0; i--){
    hist[i]->SetLineColor(mycolor[i]);
    hist[i]->SetLineWidth(3);
    if(i >= 4 && false){
      if(i%2 == 0)
	hist[i]->SetLineStyle(7);
      if(i%2 == 1)
	hist[i]->SetLineStyle(9);
      hist[i]->SetLineWidth(4);
    }
    hist[i]->SetMarkerColor(mycolor[i]);
    hist[i]->SetMarkerSize(0);
    hist[i]->SetFillColor(kWhite);
    Ntype[0]++;
    hist[i]->Draw("hist SAME");
  }

  TLegend* leg = new TLegend(0.688,0.22,0.93,0.42);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int i = 0; i < Nsample; i++)
    leg->AddEntry(hist[i],g_Samples[i]->GetTitle().c_str());
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(kWhite);
  leg->Draw("SAME");

  TLatex l;
  l.SetTextFont(132);
  l.SetNDC();
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  // l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
  l.DrawLatex(0.6,0.943,g_PlotTitle.c_str());
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.135,0.943,"#bf{CMS} Simulation Preliminary");
  l.SetTextSize(0.04);
  can->SaveAs("Plot_1D_Canvas.pdf");
}

