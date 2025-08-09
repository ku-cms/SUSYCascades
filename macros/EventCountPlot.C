#include "PlottingTools.h"

void format(TH2D* hist, string Xname, string Yname)
{
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitle(Xname.c_str());
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.12);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitle(Yname.c_str());
  hist->GetZaxis()->CenterTitle();
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleOffset(1.15);
  hist->GetZaxis()->SetLabelFont(42);
  hist->GetZaxis()->SetLabelSize(0.05);
  hist->GetZaxis()->SetTitle("N_{events}");
  hist->GetZaxis()->SetRangeUser(0.9 * hist->GetMinimum(0.0), 1.1 * hist->GetMaximum());
}

void drawLatex(string title, string label)
{
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  //l.DrawLatex(0.17, 0.855, title.c_str());
  l.DrawLatex(0.71, 0.943, title.c_str());
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.01, 0.943, "#bf{CMS} Simulation Preliminary");
  
  l.SetTextSize(0.045);
  l.SetTextFont(42);
  l.DrawLatex(0.7, 0.04, label.c_str());
}

void EventCountPlot()
{
  // speed up by not displaying canvas
  gROOT->SetBatch(kTRUE);
  
  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  std::vector<std::pair<std::string, ProcessList>> vec_samples;

  std::map<std::string, VS> map_vsignals;
  map_vsignals.insert(std::make_pair("T1bbbb", VS({"T1bbbb_*"})));
  // loop over signals and add to map
  for(auto p = map_vsignals.begin(); p != map_vsignals.end(); p++){
    ProcessList signals;
    for(auto s : p->second){
      signals += ST.Get(s);
    }
    vec_samples.push_back(std::make_pair(p->first, signals));
  }

  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    g_Label = p->first;
    g_Xname = "M_{NLSP} [GeV]";
    g_Xmin  = 500.;
    g_Xmax  = 3000.; 
    g_NX    = 60;
    
    g_Yname = "M_{LSP} [GeV]";
    g_Ymin  = 500.;
    g_Ymax  = 2000.;
    g_NY    = 40;

    // mass diff plot
    g_dm_Yname = "M_{NLSP} - M_{LSP} [GeV]";
    g_dm_Ymin  = 0.;
    g_dm_Ymax  = 1600.;
    g_dm_YN    = 20;
    
    // Change label for T2bW
    if (g_Label.find("T2bW") != string::npos)
    {
      g_Xname    = "M_{stop} [GeV]";
      g_dm_Yname = "M_{stop} - M_{LSP} [GeV]";
    }
    
    cout << "Processing " << g_Label << endl;

    // m_lsp vs m_nlsp
    TH2D* hist = new TH2D(
      (g_Label+"_EventCount").c_str(),
      (g_Label+"_EventCount").c_str(),
      g_NX, g_Xmin, g_Xmax,
      g_NY, g_Ymin, g_Ymax
    );
    // mass diff (dm) vs m_nlsp
    TH2D* hist_dm = new TH2D(
      (g_Label+"_dm_EventCount").c_str(),
      (g_Label+"_dm_EventCount").c_str(),
      g_NX,    g_Xmin,    g_Xmax,
      g_dm_YN, g_dm_Ymin, g_dm_Ymax
    );
    
    int Nsample = p->second.GetN();
    for(int s = 0; s < Nsample; s++){
      Process proc = p->second[s];
      int Nfile = ST.NTrees(proc);
      for(int f = 0; f < Nfile; f++){
        string file = ST.FileName(proc, f);
        string tree = ST.TreeName(proc, f);
        TChain* chain = ST.Tree(proc, f);
        ReducedBase_V2* base = new ReducedBase_V2(chain);
        gErrorIgnoreLevel = kFatal;
        int events = base->fChain->GetEntries(); 
        gErrorIgnoreLevel = 0;
        delete chain;
        string NLSP_Mass_string = tree.erase(0,4); //remove SMS_
        string LSP_Mass_string = NLSP_Mass_string;
        NLSP_Mass_string = NLSP_Mass_string.substr(0,NLSP_Mass_string.find("_"));
        double NLSP_Mass = std::stod(NLSP_Mass_string);
        size_t pos = LSP_Mass_string.find("_");
        LSP_Mass_string = LSP_Mass_string.erase(0,pos+1);
        double LSP_Mass = std::stod(LSP_Mass_string);
        hist->Fill(NLSP_Mass, LSP_Mass, events);
        hist_dm->Fill(NLSP_Mass, NLSP_Mass - LSP_Mass, events);
      }
    }

    cout << "Total " << hist->Integral() << endl;
    
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

    // format histograms
    format(hist,    g_Xname, g_Yname);
    format(hist_dm, g_Xname, g_dm_Yname);

    // output plot files
    string plot_dir     = "plots";
    string plot_name    = plot_dir + "/";
    string plot_name_dm = plot_dir + "/";
    plot_name    += hist->GetName();
    plot_name_dm += hist_dm->GetName();
    
    // create plot directory
    //boost::filesystem::create_directories(plot_dir);
    gSystem->mkdir(plot_dir.c_str());
    
    // m_lsp vs m_nlsp
    hist->Draw("COLZ");
    drawLatex(g_PlotTitle, g_Label);
    can->SaveAs((plot_name + ".pdf").c_str());
    
    // mass diff (dm) vs m_nlsp
    hist_dm->Draw("COLZ");
    drawLatex(g_PlotTitle, g_Label);
    can->SaveAs((plot_name_dm + ".pdf").c_str());
    
    // output root file
    TFile* output_file = new TFile("output_EventCountPlot.root", "UPDATE");
    can->Write();
    output_file->Close();
    
    delete can;
  }
  gApplication->Terminate(0);
}
