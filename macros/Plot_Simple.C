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

#include "../include/ReducedBase_V2.hh"
#include "../include/SampleTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/Hadronic.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

// declare global vars
string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;
string g_Yname;
double g_Ymin;
double g_Ymax;
double g_NY;
string g_Label;
// root file to store output of plots
string output_root_file = "output_Plot_Simple.root";
bool SavePDF = false; // whether or not to save pdfs of plots

// Plotting helper functions
void Plot_Hist(TH1* h, bool Scale=true, double Scale_Val = 1);
void Plot_Hist(TH2* h, bool Scale=true, double Scale_Val = 1);
void Plot_Eff(TEfficiency* e);

void Plot_Simple(){

  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  VS vsignals;
  vsignals.a("Cascades_2890220");
  vsignals.a("Cascades_3000220");

  ProcessList backgrounds = ST.Get(kBkg);

  ProcessList signals;
  for(auto s : vsignals)
    signals += ST.Get(s);

  ProcessList samples = signals;
  samples += backgrounds;
  
  g_Label = "2 lepton";
  g_NX = 128;

  vector<TH1*> hists1;
  vector<TH2*> hists2;
  vector<TEfficiency*> effs;

  int SKIP = 1;
  //int lumi = 138.; // Run 2
  int lumi = 138.+109+27+34; // Run 2&3
  bool Norm = true; // scale hists by lumi
   
  int Nsample = samples.GetN();
  
  for(int s = 0; s < Nsample; s++){
    Process proc = samples[s];
    g_PlotTitle = proc.Name();
    
    string title = proc.Name();

    bool is_data   = (proc.Type() == kData);
    bool is_bkg    = (proc.Type() == kBkg);
    bool is_signal = (proc.Type() == kSig);
    
    int Nfile = ST.NTrees(proc);

    // Declare hists here
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET;").c_str(), g_NX, 0., 500.);
    // push_back hists that you want to plot at the end (hists are filled regardless of whether or not you push_back)
    hists1.push_back(hist_MET);
    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
    hists2.push_back(hist_RISR_PTISR);   
    TEfficiency* eff_METtrig = new TEfficiency((title+"_eff_METtrig").c_str(), "Efficiency of MET trigger;Eff;MET [GeV]", g_NX, 0., 500.);
    effs.push_back(eff_METtrig);
 
    cout << "Processing " << Nfile << " files for process " << title << endl;
    for(int f = 0; f < Nfile; f++){
      string file = ST.FileName(proc, f);
      string tree = ST.TreeName(proc, f);
      
      bool is_FastSim = ST.IsFastSim(proc, f);
      bool do_FilterDilepton = ST.FilterDilepton(proc, f);
      double sample_weight = ST.GetSampleWeight(proc, f);
      
      if(is_signal)
        sample_weight *= SF.GetX20BRSF(file, tree);
      
      cout << "   Processing file " << file << " w/ tree " << tree << endl;
      cout << "      Sample weight is " << sample_weight << endl;
      if(is_FastSim)
	cout << "      Is FastSim" << endl;
      if(do_FilterDilepton)
	cout << "      Filter Out dilepton events" << endl;
      
      TChain* chain = ST.Tree(proc, f);
      
      ReducedBase_V2* base = new ReducedBase_V2(chain);
      
      int Nentry = base->fChain->GetEntries(); 
      
      // event loop
      for(int e = 0; e < Nentry; e += SKIP){
	base->GetEntry(e);
	
	if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
	  cout << "      event " << e << " | " << Nentry << endl;

        // Apply PreSelection
	
	if(do_FilterDilepton)
	  if(SF.DileptonEvent(base))
	    continue;
	
	// apply trigger to data and FullSim events
	//if(!base->METORtrigger && !is_FastSim)
	//  continue;
		
	double MET = base->MET;
	if(MET < 100)
	  continue;

	if(base->PTISR < 200.)
	  continue;

	// Cleaning cuts...
	// double x = fabs(base->dphiCMI);
	// 
	// if(base->PTCM > 200.)
	//   continue;
	// if(base->PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
	//    -2.777*x*x+1.388*x+0.8264 > 0.)
	//   continue;
	// if(base->PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
	//    -1.5625*x*x+7.8125*x-8.766 > 0.)
	//   continue;
	// End Cleaning cuts...
	  
	//if(fabs(base->dphiMET_V) > acos(-1.)/2.)
	//  continue;
	  
	if(base->RISR < 0. || base->RISR > 1.0)
	  continue;

        // Get Physics Objects

	int Nlep     = base->Nlep;
	int NjetS    = base->Njet_S;
	int NbjetS   = base->Nbjet_S;
	int NjetISR  = base->Njet_ISR;
	int NbjetISR = base->Nbjet_ISR;

	if(Nlep != 2)
	  continue;
	//if(NjetS != 1)
        //  continue;

	double minDR = 1000;
	double minMLL = 1000;
	
	int index_1, index_2;
	for(int i = 0; i < Nlep-1; i++){
	  TLorentzVector lep_1;
	  lep_1.SetPtEtaPhiM( base->PT_lep->at(i),
			      base->Eta_lep->at(i),
			      base->Phi_lep->at(i),
			      std::max(0.,base->M_lep->at(i)) );

	  
	  for(int j = i+1; j < Nlep; j++){
	    TLorentzVector lep_2;
	    lep_2.SetPtEtaPhiM( base->PT_lep->at(j),
				base->Eta_lep->at(j),
				base->Phi_lep->at(j),
				std::max(0.,base->M_lep->at(j)) );

	    if(lep_1.DeltaR(lep_2) < minDR)
	      minDR = lep_1.DeltaR(lep_2);
	    if(lep_1.DeltaR(lep_2) < minMLL)
	      minMLL = (lep_1 + lep_2).M();
	    
	  }	  
	}
	
	LepList list_a;
	LepList list_b;
	  
	int index;
	  
        // ID leps as gold, silver, or bronze
	for(int i = 0; i < base->Nlep_a; i++){
	  index = (*base->index_lep_a)[i];
	    
	  int PDGID = base->PDGID_lep->at(index);
	    
	  LepID id;
	  if(base->ID_lep->at(index) < 3 ||
	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)
	    id = kBronze;
	  else if(base->SIP3D_lep->at(index) > 2.)
	    id = kSilver;
	  else
	    id = kGold;
	  LepFlavor flavor;
	  if(abs(PDGID) == 11)
	    flavor = LepFlavor::kElectron;
	  else
	    flavor = LepFlavor::kMuon;
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  LepSource source = LepSource(base->SourceID_lep->at(index));
	    
	  list_a += Lep(flavor, charge, id, source);
	}
	for(int i = 0; i < base->Nlep_b; i++){
	  index = (*base->index_lep_b)[i];
	  
	  int PDGID = base->PDGID_lep->at(index);

	  LepID id;
	  if(base->ID_lep->at(index) < 3 ||
	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)
	    id = kBronze;
	  else if(base->SIP3D_lep->at(index) > 2.)
	    id = kSilver;
	  else
	    id = kGold;
	  LepFlavor flavor;
	  if(abs(PDGID) == 11)
	    flavor = LepFlavor::kElectron;
	  else
	    flavor = LepFlavor::kMuon;
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  LepSource source = LepSource(base->SourceID_lep->at(index));
	  
	  list_b += Lep(flavor, charge, id, source);
	}

        // get variables from root files using base class
	double gammaT = base->gammaT;
	double Mperp = base->Mperp;
	double RISR = base->RISR;
	double PTISR = base->PTISR;
	
	double weight = (base->weight != 0.) ? base->weight : 1.;
	weight *= double(SKIP);
	
        // Fill hists, effs, etc.
	hist_MET->Fill(MET, weight);
        hist_RISR_PTISR->Fill(RISR, PTISR, weight);
        eff_METtrig->Fill(base->METtrigger, MET, weight);
      }
      delete base;
      delete chain;
    } // for(int f = 0; f < Nfile; f++)

    // call plotting functions after looping over files in sample
    for(int hist1 = 0; hist1 < int(hists1.size()); hist1++)
      Plot_Hist(hists1[hist1], Norm, lumi);
    for(int hist2 = 0; hist2 < int(hists2.size()); hist2++)
      Plot_Hist(hists2[hist2], Norm, lumi);
    for(int eff = 0; eff < int(effs.size()); eff++)
      Plot_Eff(effs[eff]);
    hists1.clear();
    hists2.clear();
    effs.clear();

  } // for(int s = 0; s < Nsample; s++)

} // End of macro

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
  l.DrawLatex(0.7,0.04,g_Label.c_str());

  string pdf_title = g_PlotTitle+"_";
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

void Plot_Hist(TH1* h, bool Scale, double Scale_Val){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  h->Scale(Scale_Val);
  string title = h->GetName();
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
  h->GetYaxis()->SetTitle(("N_{events} / "+std::to_string(Scale_Val)+" fb^{-1}").c_str());
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

  string pdf_title = g_PlotTitle+"_";
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

void Plot_Hist(TH2* h, bool Scale, double Scale_Val){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  h->Scale(Scale_Val);
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
  h->GetZaxis()->SetTitle(("N_{events} / "+std::to_string(Scale_Val)+" fb^{-1}").c_str());
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

  string pdf_title = g_PlotTitle+"_";
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
