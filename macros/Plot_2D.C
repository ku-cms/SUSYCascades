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

#include "../include/ReducedBase_slim.hh"
#include "../include/SampleTool.hh"
#include "../include/CategoryTool.hh"
#include "../include/ScaleFactorTool.hh"
#include "../include/Leptonic.hh"
#include "../include/Hadronic.hh"
#include "../include/FitReader.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;

string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;
string g_Yname;
double g_Ymin;
double g_Ymax;
double g_NY;

using namespace RestFrames;

void Plot_2D(){

  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";

  cout << "Initializing sample maps from path " << NtuplePath << " for year " << 2017 << endl;
  SampleTool ST(NtuplePath, 2017);
  
  ScaleFactorTool SF;
  CategoryTool CT;
  CategoryList Categories = CT.GetCategories();

  CategoryTreeTool CTTool;
  CategoryTree CT_0L = CTTool.GetCategories_0L();
  CategoryTree CT_1L = CTTool.GetCategories_1L();
  CategoryTree CT_2L = CTTool.GetCategories_2L();
  CategoryTree CT_3L = CTTool.GetCategories_3L();

  int depth0 = CT_1L.GetDepth();
  vector<const CategoryTree*> CTs;
  CT_1L.GetListDepth(CTs, depth0-2);

  int iCAT = 2;

  //const CategoryTree* myCT = CTs[iCAT];
  const CategoryTree* myCT = &CT_0L;
  //Categories = Categories.Filter(CT_2L);

  VS vsignals;
  //vsignals.a("TChiWZ_2500245").a("TChiWZ_2500240").a("TChiWZ_2500220").a("TChiWZ_2500200").a("TChiWZ_2500160");
  //vsignals.a("TSlepSlep_2500245").a("TSlepSlep_2500240").a("TSlepSlep_2500220").a("TSlepSlep_2500200").a("TSlepSlep_2500170");
  //vsignals.a("T2tt_5000490").a("T2bW_5000490").a("T2tt_5000480").a("T2bW_5000480");
  //vsignals.a("T2tt_5000480").a("T2tt_5000450").a("T2tt_5000400").a("T2tt_5000325");
  //vsignals.a("T2bW_5000480").a("T2bW_5000450").a("T2bW_5000420").a("T2bW_5000325");
  //vsignals.a("T2tt_5000450").a("T2bW_5000420").a("T2tt_5000375").a("T2bW_5000325");
  //vsignals.a("T2bb_5000475").a("T2bb_5000450").a("T2bb_5000400").a("T2bb_5000350");
  // VS all = signals;
  // all.a("ttbar").a("ST").a("DB").a("ZDY").a("TB").a("QCD").a("Wjets").a("Fakes");

  //vsignals.a("T2tt_6000590");
  vsignals.a("TChiWZ_2750245");
  //vsignals.a("ttbar");
  
  //ProcessList backgrounds = ST.Get(kBkg).Remove("QCD");
  ProcessList backgrounds = ST.Get(kBkg).Filter("Wjets");

  ProcessList signals;
  for(auto s : vsignals)
    signals += ST.Get(s);
  
  string g_Label = "No selection";
  //string g_Label = "1 lepton";


  //g_Xname = "E_{lep #scale[0.8]{#perp}}^{ S}  [GeV]";
  //g_Xname = "#bar{M}_{#tilde{#chi_{2}} #scale[0.8]{#perp}}  [GeV]";
  // g_Xname = "p_{T}^{lep} [GeV]";
  g_Yname = "M_{#perp}  [GeV]";
  //g_Xname = "#Delta R (#it{l} #it{l} )";
  //g_Xmin = 0.0;
  //g_Xmax = 1.; 
  //g_NX = 32;
  //g_Yname = "#Delta R (#it{l} #it{l} )";
  //g_Yname = "m(#it{l} #it{l} ) [GeV]";
  //g_Yname = "#gamma_{T}";
  //g_Yname = "p_{T}^{ISR}";
  g_Ymin = 0.;
  g_Ymax = 100.;
  g_NY = 32;

   g_Xname = "R_{ISR}";
   g_Xmin = 0.5;
   g_Xmax = 1.; 
   g_NX = 32;
  // g_Yname = "p_{T}^{ISR} [GeV]";
  // g_Ymin = -1.;
  // g_Ymax = 1.;
  // g_NY = 32;

  TH2D* hist = new TH2D("hist", "hist",
			g_NX, g_Xmin, g_Xmax,
			g_NY, g_Ymin, g_Ymax);
  
  int SKIP = 100;
  //ProcessList samples = signals;
  ProcessList samples = backgrounds;
  //samples += backgrounds;

   g_PlotTitle = samples[0].Name();
   //g_PlotTitle = "#tilde{#chi}_{1}^{#pm} #tilde{#chi}_{2}^{0} -> #tilde{#chi}_{1}^{0} #tilde{#chi}_{1}^{0}, m_{#tilde{#chi}_{1}^{#pm}} = m_{#tilde{#chi}_{2}^{0}} = 275, m_{#tilde{#chi}_{1}^{0}} = 245";
   
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
      
      // event loop
      for(int e = 0; e < Nentry; e += SKIP){
	base->GetEntry(e);
	
	if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
	  cout << "      event " << e << " | " << Nentry << endl;
	
	if(do_FilterDilepton)
	  if(SF.DileptonEvent(base))
	    continue;
	
	// apply trigger to data and FullSim events
	if(!base->METORtrigger && !is_FastSim)
	  continue;
		
	if(base->MET < 150)
	  continue;
	  
	// if(base->PTISR < 500.)
	//   continue;

		// if(base->MET < 200)
	//   continue;
	
	if(base->PTISR < 250.)
	  continue;

	// if(base->Nbjet_ISR->at(1) > 0)
	//   continue;

	// if(base->index_SV_b->size() < 2)
	//   continue;
	// if(base->index_SV_b->at(1).size() < 1)
	//   continue;
	// if(base->NSV <= base->index_SV_b->at(1)[0])
	//   continue;
	// double etaSV = base->Eta_SV->at(base->index_SV_b->at(1)[0]);
	// double phiSV = base->Phi_SV->at(base->index_SV_b->at(1)[0]);
	// TVector3 SV;
	// SV.SetPtEtaPhi(1,etaSV,phiSV);
	// double dphiSVM = fabs(SV.DeltaPhi(ETMiss));

	// if(dphiML > 1.)
	//   continue;

	// if(base->NSV_S->at(1) > 0)
	//   continue;

	// if(base->NSV_S->at(1) < 1)
	//   continue;

	// if(base->MS->at(1) > 80.)
	//   continue;
	
	// if(base->Njet_a->at(1)+base->NSV_a->at(1) != 0)
	//   continue;

	// if(base->Njet_b->at(1)+base->NSV_b->at(1) != 1)
	//   continue;

	double x = fabs(base->dphiCMI);
	
	if(base->PTCM > 200.)
	  continue;
	if(base->PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
	   -2.777*x*x+1.388*x+0.8264 > 0.)
	  continue;
	if(base->PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
	   -1.5625*x*x+7.8125*x-8.766 > 0.)
	  continue;
	  
	if(base->RISR < 0.5 || base->RISR > 1.0)
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

	g_Label = "2L 1J";
	if(Nlep != 2)
	  continue;
	if(NjetS != 1)
          continue;

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
	    flavor = kElectron;
	  else
	    flavor = kMuon;
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
	    flavor = kElectron;
	  else
	    flavor = kMuon;
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
	Event.AddGenericVal(GenericVal(500.));
	//Event.AddGenericVal(GenericVal(base->PTISR));
	//Event.AddGenericVal(double(gammaT));
	Event.AddGenericVal(GenericVal(1.));
	Event.AddGenericVal(SVmaxeta);
	
	int eindex = Categories.Find(Event);
	if(eindex < 0){
	  continue;
	}	

	/////////////////
	
	double weight = (base->weight != 0.) ? base->weight : 1.;

	// if(base->RISR->at(1) < 0.975)
	//   continue;

	
	hist->Fill(RISR, Mperp, weight*double(SKIP));
      }
      
      delete base;
      delete chain;
    }
  }

  cout << "Total " << hist->Integral()*137. << endl;
  
  hist->Scale(137.);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = (TCanvas*) new TCanvas("can","can",700.,600);

  can->SetLeftMargin(0.15);
  can->SetRightMargin(0.18);
  can->SetBottomMargin(0.15);
  can->SetGridx();
  can->SetGridy();
  can->SetLogz();
  can->Draw();
  can->cd();
  hist->Draw("COLZ");
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitle(g_Xname.c_str());
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitle(g_Yname.c_str());
  hist->GetZaxis()->CenterTitle();
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(0.055);
  hist->GetZaxis()->SetTitleOffset(1.05);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelSize(0.05);
  hist->GetZaxis()->SetTitle("N_{events} / 137 fb^{-1}");
  hist->GetZaxis()->SetRangeUser(0.9*hist->GetMinimum(0.0),1.1*hist->GetMaximum());

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

  can->SaveAs("can.pdf");
  TFile* file = new TFile("output_Plot_2D.root","RECREATE");
  can->Write();
  file->Close();
  //delete can;

}
