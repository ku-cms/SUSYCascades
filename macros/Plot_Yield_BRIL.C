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
#include <TApplication.h>
#include <TLorentzVector.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include "../include/ReducedBase_slim.hh"
#include "../include/SampleTool.hh"
#include "../include/CategoryTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/Hadronic.hh"
//#include "../include/FitReader.hh"
#include "../include/BRILTool.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;

using namespace RestFrames;

TCanvas* PlotGraphs(vector<TGraphErrors*>& grs, vector<string>& labels, string title, string label, string scan);

void Plot_Yield_BRIL(){
  RestFrames::SetStyle();

  int year = 2018;
  bool lumi_only = false;
  BRILTool bril;
  bril.BuildMap("json/BRIL/brilcalc_"+std::to_string(year)+".txt");

  string NtuplePath = "/home/t3-ku/z374f439/storage/crogan/";
  //string NtuplePath = "root://xrootd.unl.edu//store/user/zflowers/crogan/";
  //string NtuplePath = "/Users/christopherrogan/Dropbox/SAMPLES/EWKino/NANO/NEW_21_09_20/";


  cout << "Initializing sample maps from path " << NtuplePath << " for year " << year << endl;
  SampleTool ST(NtuplePath, year);
  
  ScaleFactorTool SF;

  std::map<int,double> Ncount;
  std::map<int,double> Ncount_0L[4];
  std::map<int,double> Ncount_1L[4];
  std::map<int,double> Ncount_2L[3];
  std::map<int,double> Ncount_3L;
  
  // g_Xname = "fill";
  // if(year == 2016){
  //   g_Xmin = 4915;
  //   g_Xmax = 5451; 
  // }
  // else if(year == 2017){
  //   g_Xmin = 5830;
  //   g_Xmax = 6380; 
  // }
  // else if(year == 2018){
  //   g_Xmin = 6615;
  //   g_Xmax = 7334; 
  // }
  // else cout << "check year! " << endl;

  // g_NX = g_Xmax - g_Xmin;

  CategoryTool CT;
  CategoryList Categories;
  Categories += CT.GetCategories(false);

  ProcessList samples = ST.Get(kData);
    
  int Nsample = samples.GetN();
  
  for(int s = 0; s < Nsample; s++){
    Process proc = samples[s];
    
    string title = proc.Name();
   
    bool is_data   = (proc.Type() == kData);
    bool is_bkg    = (proc.Type() == kBkg);
    bool is_signal = (proc.Type() == kSig);
    
    int Nfile = ST.NTrees(proc);
    
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
      
      ReducedBase* base = new ReducedBase(chain);
      
      int Nentry = base->fChain->GetEntries();
      
      int SKIP = 1;
      
      // event loop
      for(int e = 0; e < Nentry; e += SKIP){
	base->GetEntry(e);
	
	if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
	  cout << "      event " << e << " | " << Nentry << endl;
	
	if(do_FilterDilepton)
	  if(SF.DileptonEvent(base))
	    continue;

	if(base->runnum >= 319077 && is_data && year == 2018){
	  if(base->HEM_Veto){
	    continue;
	  }
	}
	
	// apply trigger to data and FullSim events
	if(!base->METORtrigger && !is_FastSim)
	  continue;
		
	if(base->MET < 150)
	  continue;
	
	if(base->PTISR < 250.)
	  continue;

	if(base->RISR < 0.45 || base->RISR > 1.0)
	  continue;

	double x = fabs(base->dphiCMI);
	
	if(base->PTCM > 200.)
	  continue;
	if(base->PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
	   -2.777*x*x+1.388*x+0.8264 > 0.)
	  continue;
	if(base->PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
	   -1.5625*x*x+7.8125*x-8.766 > 0.)
	  continue;
	  
	if(fabs(base->dphiMET_V) > acos(-1.)/2.)
	  continue;

	int Nlep     = base->Nlep;
	int NjetS    = base->Njet_S;
	int NbjetS   = base->Nbjet_S;
	int NjetISR  = base->Njet_ISR;
	int NbjetISR = base->Nbjet_ISR;
	int NSV      = base->NSV_S;

	if(Nlep + NjetS + NSV < 1)
	  continue;

	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	
	LepList list_a;
	LepList list_b;
	  
	int index;

	for(int i = 0; i < base->Nlep_a; i++){
	  index = (*base->index_lep_a)[i];
	    
	  int PDGID = base->PDGID_lep->at(index);
	    
	  LepID id;
	  if(base->ID_lep->at(index*2) < 3 ||
	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)
	    id = kBronze;
	  else if(base->SIP3D_lep->at(index) > 2.)
	    id = kSilver;
	  else
	    id = kGold;
	  LepFlavor flavor;
	  if(abs(PDGID) == 11)
	    flavor = (LepFlavor)kElec;
	  else
	    flavor = (LepFlavor)kMuon;
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  //LepSource source = LepSource(base->SourceID_lep->at(index));
	  LepSource source = LepSource(base->ID_lep->at(index*2+1));
	    
	  list_a += Lep(flavor, charge, id, source);
	}
	for(int i = 0; i < base->Nlep_b; i++){
	  index = (*base->index_lep_b)[i];
	  
	  int PDGID = base->PDGID_lep->at(index);

	  LepID id;
	  if(base->ID_lep->at(index*2) < 3 ||
	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)
	    id = kBronze;
	  else if(base->SIP3D_lep->at(index) > 2.)
	    id = kSilver;
	  else
	    id = kGold;
	  LepFlavor flavor;
	  if(abs(PDGID) == 11)
	    flavor = (LepFlavor)kElec;
	  else
	    flavor = (LepFlavor)kMuon;
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  //LepSource source = LepSource(base->SourceID_lep->at(index));
	  LepSource source = LepSource(base->ID_lep->at(index*2+1));
	  
	  list_b += Lep(flavor, charge, id, source);
	}

	// SV eta
	double SVmaxeta = 0.0001; // 1 is fine b/c less than 1.5 cutoff
	for(int ie = 0; ie < base->NSV_S; ie++)
	  if(fabs(base->Eta_SV->at(ie)) > SVmaxeta)
	    SVmaxeta = fabs(base->Eta_SV->at(ie));

	double sveta = SVmaxeta;

	// gammaT calc
	double MST =
	  sqrt(base->MX3a_BoostT*base->MX3a_BoostT+base->PX3_BoostT*base->PX3_BoostT) +
	  sqrt(base->MX3b_BoostT*base->MX3b_BoostT+base->PX3_BoostT*base->PX3_BoostT);

	double gammaT = 2.*base->Mperp / MST;
	double Mperp = base->Mperp;
	double RISR = base->RISR;
	double PTISR = base->PTISR;

	Category Event(Leptonic(list_a, list_b),
		       Hadronic(NjetS, NbjetS, NSV),
		       Hadronic(NjetISR, NbjetISR, base->NSV_ISR));
	Event.AddGenericVal(GenericVal(base->PTISR));
	Event.AddGenericVal(double(gammaT));
	Event.AddGenericVal(SVmaxeta);
	
	int eindex = Categories.Find(Event);
	if(eindex < 0){
	  continue;
	}	
       
	/////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////
	
        int fill = bril.GetFillFromRun(base->runnum);

	if(Ncount.count(fill) < 1)
	  Ncount[fill] = 0.;
	Ncount[fill]++;
	
	if(Nlep == 0){
	  if(Ncount_0L[0].count(fill) < 1)
	    Ncount_0L[0][fill] = 0.;
	  Ncount_0L[0][fill]++;

	  if(NjetS >= 4){
	    if(Ncount_0L[3].count(fill) < 1)
	      Ncount_0L[3][fill] = 0.;
	    Ncount_0L[3][fill]++;
	  } else if(NjetS >= 2){
	    if(Ncount_0L[2].count(fill) < 1)
	      Ncount_0L[2][fill] = 0.;
	    Ncount_0L[2][fill]++;
	  } else {
	     if(Ncount_0L[1].count(fill) < 1)
	       Ncount_0L[1][fill] = 0.;
	     Ncount_0L[1][fill]++;
	  }
	}

	if(Nlep == 1){
	  if(Ncount_1L[0].count(fill) < 1)
	    Ncount_1L[0][fill] = 0.;
	  Ncount_1L[0][fill]++;

	  if(NjetS >= 3){
	    if(Ncount_1L[3].count(fill) < 1)
	      Ncount_1L[3][fill] = 0.;
	    Ncount_1L[3][fill]++;
	  } else if(NjetS >= 1){
	    if(Ncount_1L[2].count(fill) < 1)
	      Ncount_1L[2][fill] = 0.;
	    Ncount_1L[2][fill]++;
	  } else {
	    if(Ncount_1L[1].count(fill) < 1)
	       Ncount_1L[1][fill] = 0.;
	     Ncount_1L[1][fill]++;
	  }
	}

	if(Nlep == 2){
	  if(Ncount_2L[0].count(fill) < 1)
	    Ncount_2L[0][fill] = 0.;
	  Ncount_2L[0][fill]++;

	  if(NjetS >= 1){
	    if(Ncount_2L[2].count(fill) < 1)
	      Ncount_2L[2][fill] = 0.;
	    Ncount_2L[2][fill]++;
	  } else {
	    if(Ncount_2L[1].count(fill) < 1)
	      Ncount_2L[1][fill] = 0.;
	    Ncount_2L[1][fill]++;
	  }
	}

	if(Nlep == 3){
	  if(Ncount_3L.count(fill) < 1)
	    Ncount_3L[fill] = 0.;
	  Ncount_3L[fill]++;
	}
      }
      
      delete base;
      delete chain;
    }
  }

  vector<double> X;
  vector<double> Y;
  vector<double> Y_err;
  vector<double> X_Nlep[4];
  vector<double> Y_Nlep[4];
  vector<double> Y_err_Nlep[4];
  vector<double> X_0L[3];
  vector<double> Y_0L[3];
  vector<double> Y_err_0L[3];
  vector<double> X_1L[3];
  vector<double> Y_1L[3];
  vector<double> Y_err_1L[3];
  vector<double> X_2L[2];
  vector<double> Y_2L[2];
  vector<double> Y_err_2L[2];

  for(auto c : Ncount){
    int fill = c.first;
    int N = c.second;
    double L = bril.GetIntegratedLumi(fill);
    if(L > 0){
      X.push_back(fill);
      Y.push_back(N / L);
      Y_err.push_back(sqrt(N) / L);
    }
  }

  for(auto c : Ncount_0L[0]){
    int fill = c.first;
    int N = c.second;
    double L = bril.GetIntegratedLumi(fill);
    if(L > 0){
      X_Nlep[0].push_back(fill);
      Y_Nlep[0].push_back(N / L);
      Y_err_Nlep[0].push_back(sqrt(N) / L);
    }
  }

  for(auto c : Ncount_1L[0]){
    int fill = c.first;
    int N = c.second;
    double L = bril.GetIntegratedLumi(fill);
    if(L > 0){
      X_Nlep[1].push_back(fill);
      Y_Nlep[1].push_back(N / L);
      Y_err_Nlep[1].push_back(sqrt(N) / L);
    }
  }

  for(auto c : Ncount_2L[0]){
    int fill = c.first;
    int N = c.second;
    double L = bril.GetIntegratedLumi(fill);
    if(L > 0){
      X_Nlep[2].push_back(fill);
      Y_Nlep[2].push_back(N / L);
      Y_err_Nlep[2].push_back(sqrt(N) / L);
    }
  }

  for(auto c : Ncount_3L){
    int fill = c.first;
    int N = c.second;
    double L = bril.GetIntegratedLumi(fill);
    if(L > 0){
      X_Nlep[3].push_back(fill);
      Y_Nlep[3].push_back(N / L);
      Y_err_Nlep[3].push_back(sqrt(N) / L);
    }
  }

  for(int i = 0; i < 3; i++){
    for(auto c : Ncount_0L[i+1]){
      int fill = c.first;
      int N = c.second;
      double L = bril.GetIntegratedLumi(fill);
      if(L > 0){
	X_0L[i].push_back(fill);
	Y_0L[i].push_back(N / L);
	Y_err_0L[i].push_back(sqrt(N) / L);
      }
    }
  }

  for(int i = 0; i < 3; i++){
    for(auto c : Ncount_1L[i+1]){
      int fill = c.first;
      int N = c.second;
      double L = bril.GetIntegratedLumi(fill);
      if(L > 0){
	X_1L[i].push_back(fill);
	Y_1L[i].push_back(N / L);
	Y_err_1L[i].push_back(sqrt(N) / L);
      }
    }
  }

  for(int i = 0; i < 2; i++){
    for(auto c : Ncount_2L[i+1]){
      int fill = c.first;
      int N = c.second;
      double L = bril.GetIntegratedLumi(fill);
      if(L > 0){
	X_2L[i].push_back(fill);
	Y_2L[i].push_back(N / L);
	Y_err_2L[i].push_back(sqrt(N) / L);
      }
    }
  }

  string title = std::to_string(year)+" Data";
  
  vector<TGraphErrors*> gr_incl;
  gr_incl.push_back(new TGraphErrors(X.size(), &X[0], &Y[0], 0, &Y_err[0]));
  vector<string> labels;
  
  TCanvas* can_incl = PlotGraphs(gr_incl, labels, title, "Preselection", "incl");
  
  vector<TGraphErrors*> gr_Nlep;
  gr_Nlep.push_back(new TGraphErrors(X_Nlep[0].size(), &X_Nlep[0][0], &Y_Nlep[0][0], 0, &Y_err_Nlep[0][0]));
  gr_Nlep.push_back(new TGraphErrors(X_Nlep[1].size(), &X_Nlep[1][0], &Y_Nlep[1][0], 0, &Y_err_Nlep[1][0]));
  gr_Nlep.push_back(new TGraphErrors(X_Nlep[2].size(), &X_Nlep[2][0], &Y_Nlep[2][0], 0, &Y_err_Nlep[2][0]));
  gr_Nlep.push_back(new TGraphErrors(X_Nlep[3].size(), &X_Nlep[3][0], &Y_Nlep[3][0], 0, &Y_err_Nlep[3][0]));
  labels.clear();
  labels.push_back("0L J X");
  labels.push_back("1L J X");
  labels.push_back("2L J X");
  labels.push_back("3L J X");

  TCanvas* can_Nlep = PlotGraphs(gr_Nlep, labels, title, "Preselection", "Nlep");
  
  vector<TGraphErrors*> gr_0L;
  gr_0L.push_back(new TGraphErrors(X_0L[0].size(), &X_0L[0][0], &Y_0L[0][0], 0, &Y_err_0L[0][0]));
  gr_0L.push_back(new TGraphErrors(X_0L[1].size(), &X_0L[1][0], &Y_0L[1][0], 0, &Y_err_0L[1][0]));
  gr_0L.push_back(new TGraphErrors(X_0L[2].size(), &X_0L[2][0], &Y_0L[2][0], 0, &Y_err_0L[2][0]));
  labels.clear();
  labels.push_back("0L 0-1J X");
  labels.push_back("0L 2-3J X");
  labels.push_back("0L 4-5J X");

  TCanvas* can_0L = PlotGraphs(gr_0L, labels, title, "Preselection + 0L", "0L");

  vector<TGraphErrors*> gr_1L;
  gr_1L.push_back(new TGraphErrors(X_1L[0].size(), &X_1L[0][0], &Y_1L[0][0], 0, &Y_err_1L[0][0]));
  gr_1L.push_back(new TGraphErrors(X_1L[1].size(), &X_1L[1][0], &Y_1L[1][0], 0, &Y_err_1L[1][0]));
  gr_1L.push_back(new TGraphErrors(X_1L[2].size(), &X_1L[2][0], &Y_1L[2][0], 0, &Y_err_1L[2][0]));
  labels.clear();
  labels.push_back("1L 0-1J X");
  labels.push_back("1L 2-3J X");
  labels.push_back("1L 4-5J X");

  TCanvas* can_1L = PlotGraphs(gr_1L, labels, title, "Preselection + 1L", "1L");

  vector<TGraphErrors*> gr_2L;
  gr_2L.push_back(new TGraphErrors(X_2L[0].size(), &X_2L[0][0], &Y_2L[0][0], 0, &Y_err_2L[0][0]));
  gr_2L.push_back(new TGraphErrors(X_2L[1].size(), &X_2L[1][0], &Y_2L[1][0], 0, &Y_err_2L[1][0]));
  labels.clear();
  labels.push_back("2L 0J X");
  labels.push_back("2L 1-2J X");

  TCanvas* can_2L = PlotGraphs(gr_2L, labels, title, "Preselection + 2L", "2L");

  string ofile = "Yields_"+std::to_string(year);
  TFile* fout = new TFile((ofile+".root").c_str(),"RECREATE");
  fout->cd();
  
  can_incl->SaveAs((ofile+"_incl.gif").c_str());
  can_incl->Write("",TObject::kOverwrite);
  can_Nlep->SaveAs((ofile+"_Nlep.gif").c_str());
  can_Nlep->Write("",TObject::kOverwrite);
  can_0L->SaveAs((ofile+"_0L.gif").c_str());
  can_0L->Write("",TObject::kOverwrite);
  can_1L->SaveAs((ofile+"_1L.gif").c_str());
  can_1L->Write("",TObject::kOverwrite);
  can_2L->SaveAs((ofile+"_2L.gif").c_str());
  can_2L->Write("",TObject::kOverwrite);
  
  fout->Close();
  delete fout;

}


TCanvas* PlotGraphs(vector<TGraphErrors*>& grs, vector<string>& labels, string title, string label, string scan){
  int Ngr = grs.size();
  if(Ngr < 1)
    return nullptr; 

  TCanvas* can = (TCanvas*) new TCanvas(("can_"+scan).c_str(),
					("can_"+scan).c_str(),700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.04);
  can->SetBottomMargin(0.15);
  can->SetGridy();
  can->SetLogy();
  can->Draw();
  can->cd();

  int mycolor[8];
  mycolor[0] = kBlue+2;
  mycolor[1] = kGreen+3;
  mycolor[2] = kRed+1;
  mycolor[3] = kYellow+2;
  mycolor[4] = kMagenta+1;
  mycolor[5] = kMagenta+2;
  mycolor[6] = kCyan+2;
  mycolor[7] = kCyan+3;

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < Ngr; i++){
    grs[i]->SetMarkerSize(1);
    grs[i]->SetLineWidth(2);
    if(Ngr > 1){
      grs[i]->SetMarkerColor(mycolor[i]);
      grs[i]->SetLineColor(mycolor[i]);
    }
    mg->Add(grs[i]);
  }
  
  mg->Draw("AP");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleFont(42);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->SetLabelFont(42);
  mg->GetXaxis()->SetLabelSize(0.045);
  mg->GetXaxis()->SetTitle("fill number");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleFont(42);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(1.14);
  mg->GetYaxis()->SetLabelFont(42);
  mg->GetYaxis()->SetLabelSize(0.045);
  mg->GetYaxis()->SetTitle("N_{event} / fill integrated lumi [fb^{-1}]");
  if(Ngr > 1)
    mg->GetYaxis()->SetRangeUser(mg->GetYaxis()->GetXmin(), 8.*mg->GetYaxis()->GetXmax());

  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.SetTextAlign(12);
  l.DrawLatex(0.15,0.943,"#bf{#it{CMS}} work-in-progress");
  l.SetTextAlign(32);
  l.DrawLatex(0.95,0.943,title.c_str());
  l.SetTextAlign(12);
  l.DrawLatex(0.18,0.86,label.c_str());

  if(Ngr < 2)
    return can;
  
  TLegend* leg = new TLegend(0.75,0.675,0.955,0.892);
  leg->SetTextFont(132);
  leg->SetTextSize(0.044);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  for(int i = 0; i < Ngr; i++)
    leg->AddEntry(grs[i],labels[i].c_str());
  leg->Draw("SAME");

  return can;
  
}
