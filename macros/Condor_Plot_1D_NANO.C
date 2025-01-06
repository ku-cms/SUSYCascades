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
#include <TH1D.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TError.h>
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include "SUSYNANOBase.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

string g_PlotTitle;
string g_Xname;
double g_Xmin;
double g_Xmax;
double g_NX;

string load_file_from_list(string fileListName, int line_number);

void Condor_Plot_1D_NANO(string filetag = "", string dataset = "", int fileline = 0, string proc_Name = "", double weight = 0.0){
  
  RestFrames::SetStyle();

  string plot_folder = "";
  plot_folder += "plots/";
  string g_Label = "";
  g_Xmin = 0.;
  g_Xmax = 500.; 
  g_NX = 128;
  vector<TH1D*> hists1;
  vector<TH2D*> hists2;
  vector<TEfficiency*> effs;
  string title = proc_Name;
  g_PlotTitle = title;

  //TH1D* hist = new TH1D(("hist_"+filetag+"_"+dataset+"_"+std::to_string(fileline)).c_str(), ("hist_"+filetag+"_"+dataset+"_"+std::to_string(fileline)).c_str(), g_NX, g_Xmin, g_Xmax);
  TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET;").c_str(), g_NX, g_Xmin, g_Xmax);
  hists1.push_back(hist_MET);

  TChain* chain = new TChain("Events");
  string inputfile_name = load_file_from_list(dataset+".txt",fileline);
  chain->Add(inputfile_name.c_str());

  SUSYNANOBase* base = new SUSYNANOBase(chain);
  int Nentry = base->fChain->GetEntries();
  int SKIP = 1000;
  
  // event loop
  for(int e = 0; e < Nentry; e += SKIP){
    base->GetEntry(e);
    
    if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
      cout << "      event " << e << " | " << Nentry << endl;
    double hweight = base->genWeight*weight*double(SKIP)/1.e7; // need 1.e7 factor to normalize from input!
    hist_MET->Fill(base->MET_pt, hweight);
  }
  
  delete base;
  delete chain;
  
  TFile* output_file = new TFile(("output_Plot_1D_NANO_"+filetag+"_"+dataset+"_"+std::to_string(fileline)+".root").c_str(),"RECREATE");
  output_file->mkdir(plot_folder.c_str());
  output_file->cd(plot_folder.c_str());
  //output_file->SetDirectory(gDirectory);
  for(int hist1 = 0; hist1 < int(hists1.size()); hist1++) hists1[hist1]->Write();
  for(int hist2 = 0; hist2 < int(hists2.size()); hist2++) hists2[hist2]->Write();
  for(int eff = 0; eff < int(effs.size()); eff++) effs[eff]->Write();
  output_file->Close(); 
  gApplication->Terminate(0);
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

void load_files_from_list(vector<string>& files, string fileListName)
{
 std::ifstream infile(fileListName);
 string line = "";
 while(getline(infile,line))
 {
  files.push_back(line);
 }
}

