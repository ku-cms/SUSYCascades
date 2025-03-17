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
vector<int> colors = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-1,kBlack,kAzure+2,kGray};
vector<int> markers = {20,21,22,23,29,33,34,43,49};

// Plotting helper functions
void Plot_Hist(TH1* h, bool Scale=true, double Scale_Val = 1);
void Plot_Hist(TH2* h, bool Scale=true, double Scale_Val = 1);
void Plot_Eff(TEfficiency* e);
void Plot_Stack(vector<TH1*> vect_h, const std::map<std::string, ProcessList>& map_samples);
void Plot_Ratio(TH1* h_num, TH1* h_den);

void Plot_Simple(){

  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  g_Label = "2 lepton SR";
  //g_Label = "2 lepton ttbar CR";
  g_NX = 128;

  int SKIP = 1; // TODO: need to change this to only apply SKIP to BKG
  //double lumi = 138.; // Run 2
  //double lumi = 138.+109+27+34; // Run 2&3
  double lumi = 9.451; // Summer23BPix
  bool Norm = true; // scale hists by lumi

  std::map<std::string, ProcessList> map_samples;

  std::map<std::string, VS> map_vsignals;
  map_vsignals.insert(std::make_pair("Cascades_220", VS({"Cascades_*_220_*"})));
  map_vsignals.insert(std::make_pair("Cascades_260", VS({"Cascades_*_260_*"})));
  map_vsignals.insert(std::make_pair("Cascades_270", VS({"Cascades_*_270_*"})));
   
  // loop over backgrounds and add to map
  ProcessList backgrounds = ST.Get(kBkg);
  for(int s = 0; s < int(backgrounds.GetN()); s++){
    ProcessList bkg;
    bkg += backgrounds[s];
    map_samples.insert(std::make_pair(backgrounds[s].Name(), bkg));
  }

  // loop over data and add to map
  ProcessList data = ST.Get(kData);
  for(int s = 0; s < int(data.GetN()); s++){
    ProcessList datum;
    datum += data[s];
    //map_samples.insert(std::make_pair(data[s].Name(), datum));
  }

  // loop over signals and add to map
  for(auto p = map_vsignals.begin(); p != map_vsignals.end(); p++){
    ProcessList signals;
    for(auto s : p->second){
      signals += ST.Get(s);
    }
    map_samples.insert(std::make_pair(p->first, signals));
  }
  
  // vectors to hold 'generic' plotting objects
  vector<vector<TH1*>*> hist_stacks; // vector for holding hist stacks
  vector<TH1*> hist_stack_MET; // example of vector for stacking all MET hists
  hist_stacks.push_back(&hist_stack_MET);

  for(auto p = map_samples.begin(); p != map_samples.end(); p++){
    // vectors to hold 'generic' plotting objects
    vector<TH1*> hists1;
    vector<TH2*> hists2;
    vector<TEfficiency*> effs;

    g_PlotTitle = p->first; 
    int Nsample = p->second.GetN();
    string title = p->first;

    // Declare hists here
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET;").c_str(), g_NX/4, 100., 400.);
    // push_back hists that you want to plot at the end (hists are filled regardless of whether or not you push_back)
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot
    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
    hists2.push_back(hist_RISR_PTISR);   
    TEfficiency* eff_METtrig = new TEfficiency((title+"_eff_METtrig").c_str(), "Efficiency of MET trigger;Eff;MET [GeV]", g_NX, 0., 500.);
    effs.push_back(eff_METtrig);

    for(int s = 0; s < Nsample; s++){
      Process proc = p->second[s];
      
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
          if(!base->METORtrigger && !is_FastSim)
            continue;
          	
          double MET = base->MET;
          if(MET < 100)
            continue;

          if(base->PTISR < 200.)
          //if(base->PTISR < 400.) // SR
            continue;

          // Cleaning cuts...
          double x = fabs(base->dphiCMI);
          
          if(base->PTCM > 200.)
            continue;
          if(base->PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
             -2.777*x*x+1.388*x+0.8264 > 0.)
            continue;
          if(base->PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
             -1.5625*x*x+7.8125*x-8.766 > 0.)
            continue;
          // End Cleaning cuts...
            
          if(fabs(base->dphiMET_V) > acos(-1.)/2.)
            continue;
            
          //if(base->RISR < 0. || base->RISR > 1.0)
          //if(base->RISR < 0.4 || base->RISR > 0.7) // CR
          if(base->RISR < 0.9)
            continue;

          // Get Physics Objects
          int Nlep     = base->Nlep;
          int NjetS    = base->Njet_S;
          int NbjetS   = base->Nbjet_S;
          int NjetISR  = base->Njet_ISR;
          int NbjetISR = base->Nbjet_ISR;

          //if(NbjetISR + NbjetS != 2) continue; // CR
          if(NbjetISR + NbjetS > 1) continue; // SR

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
          LepList list_leps;  
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
            list_leps += Lep(flavor, charge, id, source);
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
            list_leps += Lep(flavor, charge, id, source);
          }

          bool skip = false;
          for(int i = 0; i < list_leps.GetN(); i++)
            if(list_leps[i].ID() != kGold) skip = true;
          if(skip) continue;

          // get variables from root files using base class
          double gammaT = base->gammaT;
          double Mperp = base->Mperp;
          double RISR = base->RISR;
          double PTISR = base->PTISR;
          
          double weight = (base->weight != 0.) ? base->weight : 1.;
          if(!is_data)
            weight *= double(SKIP);

          // Fill hists, effs, etc.
          hist_MET->Fill(MET, weight);
          hist_RISR_PTISR->Fill(RISR, PTISR, weight);
          eff_METtrig->Fill(base->METtrigger, MET, weight);
        }
        delete base;
        delete chain;
      } // for(int f = 0; f < Nfile; f++)

    } // for(int s = 0; s < Nsample; s++)

    // call plotting functions after looping over files in samples
    for(int hist1 = 0; hist1 < int(hists1.size()); hist1++)
      Plot_Hist(hists1[hist1], Norm, lumi);
    for(int hist2 = 0; hist2 < int(hists2.size()); hist2++)
      Plot_Hist(hists2[hist2], Norm, lumi);
    for(int eff = 0; eff < int(effs.size()); eff++)
      Plot_Eff(effs[eff]);
    hists1.clear();
    hists2.clear();
    effs.clear();

  } // for(auto p = map_samples.begin(); p != map_samples.end(); p++){
  for(int stack_h = 0; stack_h < int(hist_stacks.size()); stack_h++){
     g_PlotTitle = (*hist_stacks[stack_h])[0]->GetName();
     // loop through process names, if process name in title, remove that part from title and add '_stack'
     for(auto p = map_samples.begin(); p != map_samples.end(); p++){
       string name = p->first;
       // look for name of process in title
       size_t pos = g_PlotTitle.find(name);
       if(pos != std::string::npos){
         // remove name from title
         g_PlotTitle.erase(pos, name.length());
       }
     }
     // add stack to title
     g_PlotTitle = "stack"+g_PlotTitle;
     Plot_Stack(*hist_stacks[stack_h], map_samples);
  }
  gApplication->Terminate(0);
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
  string title = h->GetName();
  // do not scale data
  if(title.find("data") == std::string::npos) { h->Scale(Scale_Val); }
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

  string pdf_title = g_PlotTitle+"_";
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

void Plot_Hist(TH2* h, bool Scale, double Scale_Val){
  if(Scale_Val == 0) Scale_Val = h->Integral();
  string title = h->GetName();
  // do not scale data
  if(title.find("data") == std::string::npos) h->Scale(Scale_Val);
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

void Plot_Stack(vector<TH1*> vect_h, const std::map<std::string, ProcessList>& map_samples){
  TH1D* h_BKG = nullptr;
  TH1D* h_DATA = nullptr;
  bool isBKG = false;
  int Nsample = map_samples.size();
  if (Nsample != int(vect_h.size())){
    std::cout << "hist stack size does not match sample stack size!" << std::endl;
    return;
  }
  int index = 0;
  for(auto p = map_samples.begin(); p != map_samples.end(); p++){
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
  TCanvas* can = (TCanvas*) new TCanvas(("can_"+g_PlotTitle).c_str(),("can_"+g_PlotTitle).c_str(),1278,775);

  can->SetLeftMargin(0.13);
  can->SetRightMargin(0.04);
  can->SetBottomMargin(0.15);
  can->SetTopMargin(0.085);
  can->SetGridx();
  can->SetGridy();
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
  vect_h[imax]->GetXaxis()->SetTitleFont(132);
  vect_h[imax]->GetXaxis()->SetTitleSize(0.06);
  vect_h[imax]->GetXaxis()->SetTitleOffset(1.06);
  vect_h[imax]->GetXaxis()->SetLabelFont(132);
  vect_h[imax]->GetXaxis()->SetLabelSize(0.05);
  vect_h[imax]->GetXaxis()->SetTitle(g_Xname.c_str());
  vect_h[imax]->GetYaxis()->CenterTitle();
  vect_h[imax]->GetYaxis()->SetTitleFont(132);
  vect_h[imax]->GetYaxis()->SetTitleSize(0.06);
  vect_h[imax]->GetYaxis()->SetTitleOffset(1.06);
  vect_h[imax]->GetYaxis()->SetLabelFont(132);
  vect_h[imax]->GetYaxis()->SetLabelSize(0.05);
  //vect_h[imax]->GetYaxis()->SetTitle("a. u.");
  vect_h[imax]->GetYaxis()->SetTitle("N_{evt} / fb^{-1}");

  index = 0;
  for(auto p = map_samples.begin(); p != map_samples.end(); p++){
    if(p->second[0].Type() == kBkg){
      //if(vect_h[index]->GetEntries() == 0) continue;
      vect_h[index]->SetLineColor(kBlack);
      vect_h[index]->SetLineWidth(1.0);
      vect_h[index]->SetFillColor(colors[index]);
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
  for(auto p = map_samples.begin(); p != map_samples.end(); p++){
    if(p->second[0].Type() == kSig){
      vect_h[index]->SetLineWidth(3.0);
      vect_h[index]->SetMarkerSize(0.);
      vect_h[index]->SetMarkerColor(kBlack);
      vect_h[index]->SetLineStyle(7);
      vect_h[index]->SetLineColor(colors[index]);
      vect_h[index]->Draw("SAME HIST");
    }
    if(p->second[0].Type() == kData){
      vect_h[index]->SetLineWidth(3.0);
      vect_h[index]->SetMarkerSize(0.);
      vect_h[index]->SetMarkerColor(kBlack);
      vect_h[index]->SetLineStyle(7);
      vect_h[index]->SetLineColor(colors[index]);
      vect_h[index]->Draw("SAME");
    }
    index++;
  }

  TLegend* leg = new TLegend(0.65,0.7,0.8,0.85);
  leg->SetTextFont(132);
  leg->SetTextSize(0.045);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);
  if(isBKG) leg->AddEntry(h_BKG, "SM total");
  index = 0;
  for(auto p = map_samples.begin(); p != map_samples.end(); p++){
    if(p->second[0].Type() == kBkg)
      leg->AddEntry(vect_h[index],p->first.c_str(),"F");
    else
      leg->AddEntry(vect_h[index],p->first.c_str());
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
  l.DrawLatex(0.65,0.943,g_Label.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.15,0.943,"#bf{#it{CMS}} Internal 13 TeV Simulation");
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
      h_res->GetYaxis()->SetRangeUser(0., 8.);
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
  
  string pdf_title = g_PlotTitle+"_";
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

  string pdf_title = g_PlotTitle+"_";
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
