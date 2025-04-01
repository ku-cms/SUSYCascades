#include <iostream>
#include <map>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TKey.h>

#include "RestFrames/RestFrames.hh"

#include "PlotReader.hh"

///////////////////////////////////////////
////////// PlotReader class
///////////////////////////////////////////

PlotReader::PlotReader(const string& inputfile)
  : m_File(inputfile.c_str(), "READ") {
  
  Initialize();
}

PlotReader::~PlotReader(){
  m_File.Close();
}

void PlotReader::Initialize(){

  ReadHists();

  m_CMSLabel = "#bf{#it{CMS}} work-in-progress";
  
  m_Strings["Data"] = VS().a("data_obs").a("Data-j").a("Data-bj").a("Data-sv").a("Data-el").a("Data-mu");
  
  m_Strings["Bkg"] = VS().a("Bkg-j").a("Bkg-bj").a("Bkg-sv").a("Bkg-el").a("Bkg-mu");
  m_Title["Bkg"] = "SM backgrounds";
  m_Color["Bkg"] = 7001;
  
  m_Title["ttbar"] = "t #bar{t} + X";
  m_Color["ttbar"] = 7011;
  
  m_Title["ST"] = "single top";
  m_Color["ST"] = 7010;
  
  m_Title["DB"] = "di-bosons";
  m_Color["DB"] = 7051;

  m_Title["TB"] = "tri-bosons";
  m_Color["TB"] = 7050;

  m_Title["ZDY"] = "Z / #gamma* + jets";
  m_Color["ZDY"] = 7000;

  m_Title["Wjets"] = "W + jets";
  m_Color["Wjets"] = 7001;

  m_Title["QCD"] = "QCD multi-jets";
  m_Color["QCD"] = 7023;

  m_Title["MET"]  = "MET [GeV]";
  m_Title["RISR"] = "R_{ISR}";
  m_Title["Mperp"] = "M_{#perp } [GeV]";
  m_Title["PTISR"] = "P_{T}^{ISR} [GeV]";
  m_Title["Gperp"] = "#gamma_{#perp}";

  m_Title["NSjet"] = "N_{jet}^{S}";
  m_Title["Nbjet"] = "N_{bjet}";
  m_Title["NSbjet"] = "N_{bjet}^{S}";
  m_Title["NSV"] = "N_{SV}";
  m_Title["Nel"] = "N_{el}";
  m_Title["Nmu"] = "N_{mu}";

  m_Title["JetEta"] = "#eta_{jet}";
  m_Title["JetPhi"] = "#phi_{jet}";
  m_Title["JetPt"] = "p_{T}^{jet} [GeV]";
  m_Title["BjetEta"] = "#eta_{bjet}";
  m_Title["BjetPhi"] = "#phi_{bjet}";
  m_Title["BjetPt"] = "p_{T}^{bjet} [GeV]";
  m_Title["SVJetEta"] = "#eta_{SV}";
  m_Title["SVPhi"] = "#phi_{SV}";
  m_Title["SVPt"] = "p_{T}^{SV} [GeV]";
  m_Title["ElEta"] = "#eta_{el}";
  m_Title["ElPhi"] = "#phi_{el}";
  m_Title["ElPt"] = "p_{T}^{el} [GeV]";
  m_Title["MuEta"] = "#eta_{mu}";
  m_Title["MuPhi"] = "#phi_{mu}";
  m_Title["MuPt"] = "p_{T}^{mu} [GeV]";

  m_Title["p"] = "preselection";
  m_Title["incl"]  = "preselection";
  m_Title["0J"]  = "L 0J X";
  m_Title["1J"]  = "L 1J X";
  m_Title["2J"]  = "L 2J X";
  m_Title["3J"]  = "L 3J X";
  m_Title["4J"]  = "L 4J X";
  m_Title["5J"]  = "L 5J X";
  m_Title["0L"]  = "0L J X";
  m_Title["0L_0J"]  = "0L 0J X";
  m_Title["0L_1J"]  = "0L 1J X";
  m_Title["0L_2J"]  = "0L 2J X";
  m_Title["0L_3J"]  = "0L 3J X";
  m_Title["0L_4J"]  = "0L 4J X";
  m_Title["0L_5J"]  = "0L 5J X";

  m_Title["1L"]  = "1L J X";
  m_Title["1L_0J"]  = "1L 0J X";
  m_Title["1L_1J"]  = "1L 1J X";
  m_Title["1L_2J"]  = "1L 2J X";
  m_Title["1L_3J"]  = "1L 3J X";
  m_Title["1L_4J"]  = "1L 4J X";
  m_Title["gold1L"]  = "1L_{G} J X";
  m_Title["gold1L_0J"]  = "1L_{G} 0J X";
  m_Title["gold1L_1J"]  = "1L_{G} 1J X";
  m_Title["gold1L_2J"]  = "1L_{G} 2J X";
  m_Title["gold1L_3J"]  = "1L_{G} 3J X";
  m_Title["gold1L_4J"]  = "1L_{G} 4J X";
  m_Title["slvr1L"]  = "1L_{S} J X";
  m_Title["slvr1L_0J"]  = "1L_{S} 0J X";
  m_Title["slvr1L_1J"]  = "1L_{S} 1J X";
  m_Title["slvr1L_2J"]  = "1L_{S} 2J X";
  m_Title["slvr1L_3J"]  = "1L_{S} 3J X";
  m_Title["slvr1L_4J"]  = "1L_{S} 4J X";
  m_Title["bron1L"]  = "1L_{B} J X";
  m_Title["bron1L_0J"]  = "1L_{B} 0J X";
  m_Title["bron1L_1J"]  = "1L_{B} 1J X";
  m_Title["bron1L_2J"]  = "1L_{B} 2J X";
  m_Title["bron1L_3J"]  = "1L_{B} 3J X";
  m_Title["bron1L_4J"]  = "1L_{B} 4J X";
  m_Title["2L"]  = "2L J X";
  m_Title["2L_0J"]  = "2L 0J X";
  m_Title["2L_1J"]  = "2L 1J X";
  m_Title["2L_2J"]  = "2L 2J X";
  m_Title["gold2L"]  = "2L_{G} J X";
  m_Title["gold2L_0J"]  = "2L_{G} 0J X";
  m_Title["gold2L_1J"]  = "2L_{G} 1J X";
  m_Title["gold2L_2J"]  = "2L_{G} 2J X";
  m_Title["slvr2L"]  = "2L_{S} J X";
  m_Title["slvr2L_0J"]  = "2L_{S} 0J X";
  m_Title["slvr2L_1J"]  = "2L_{S} 1J X";
  m_Title["slvr2L_2J"]  = "2L_{S} 2J X";
  m_Title["bron2L"]  = "2L_{B} J X";
  m_Title["bron2L_0J"]  = "2L_{B} 0J X";
  m_Title["bron2L_1J"]  = "2L_{B} 1J X";
  m_Title["bron2L_2J"]  = "2L_{B} 2J X";
  m_Title["3L"]  = "3L J X";
  m_Title["3L_0J"]  = "3L 0J X";
  m_Title["3L_1J"]  = "3L 1J X";
  m_Title["gold3L"]  = "3L_{G} J X";
  m_Title["gold3L_0J"]  = "3L_{G} 0J X";
  m_Title["gold3L_1J"]  = "3L_{G} 1J X";
  m_Title["slvr3L"]  = "3L_{S} J X";
  m_Title["slvr3L_0J"]  = "3L_{S} 0J X";
  m_Title["slvr3L_1J"]  = "3L_{S} 1J X";
  m_Title["bron3L"]  = "3L_{B} J X";
  m_Title["bron3L_0J"]  = "3L_{B} 0J X";
  m_Title["bron3L_1J"]  = "3L_{B} 1J X";
    
}

  

const VS& PlotReader::GetHists() const {
  return m_Hists;
}

void PlotReader::ReadHists(){
  if(!m_File.IsOpen())
    return;

  m_File.cd("plot/hist");
  TIter list(((TDirectory*)m_File.Get("plot/hist"))->GetListOfKeys());
  TKey* key;
  string name;
  while((key = (TKey*)list.Next())){
    name = string(key->GetName());
    m_Hists += name;
  }
}

TH1D* PlotReader::GetHist(const string& name,
			  const string& var,
			  const string& proc,
			  const VS& regions) const {

  VS list = GetHists();
  VS labels;
  for(auto r : regions){
    labels += var+"_"+proc+"_"+r;
  }

  // VS list = GetHists().Filter(var);
  // list = list.Filter(proc);
  list = list.FilterOR(labels);
  
  if(list.size() < 1)
    return nullptr;

  TH1D* hist = nullptr;

  for(auto h : list){
    if(h.find("_v_") != std::string::npos)
      continue;
    TH1D* temp = (TH1D*) m_File.Get(("plot/hist/"+h).c_str());
    if(!temp)
      continue;
    if(!hist)
      hist = (TH1D*) temp->Clone(("hist_"+name).c_str());
    else
      hist->Add(temp);
  }

  return hist;
}


TH1D* PlotReader::GetHist(const string& name,
			  const string& var,
			  const string& proc,
			  const string& region) const {
  VS regions;
  regions += region;
  return GetHist(name, var, proc, regions);
}



string PlotReader::GetSignalTitle(const string& label){
  size_t p = label.rfind("_");
  if(p == std::string::npos)
    return label;

  string title = label.substr(0, p);
  int    mass  = stoi(label.substr(p+1,label.length()-p));
  return title+" "+std::to_string(mass/10000)+" "+std::to_string(mass%100000);
}

TCanvas* PlotReader::Plot1Dstack(const string& can_name,
				 const string& var,
				 const VS& procs,
				 const string& region,
				 bool norm, int rebin, double xmin, double xmax){
  RestFrames::SetStyle();

  bool b_ratio = false;

  int Nproc = procs.size();
  
  if(Nproc == 0)
    return nullptr;

  // Processes
  VS                  processes;
  vector<TH1D*>       hists;
  
  TH1D* hist_data = nullptr;
  double total_data = 0.;

  TH1D* hist_totbkg = nullptr;
  double total_totbkg = 0.;

  for(int i = 0; i < Nproc; i++){
    bool isbkg = true;
    if(procs[i] == "Data" ||
       procs[i] == "data_obs")
      isbkg = false;

    VS vproc;
    if(m_Strings.count(procs[i]) != 0)
      vproc = m_Strings[procs[i]];
    else
      vproc += procs[i];

    TH1D* hist_proc = nullptr;
    
    for(auto p : vproc){
      if(!hist_proc){
	hist_proc = (TH1D*) GetHist(can_name+"hist_"+p,
				    var, p, region);
      } else {
	TH1D* temp = (TH1D*) GetHist(can_name+"temp_"+p,
				     var, p, region);
	if(temp){
	  hist_proc->Add(temp);
	  delete temp;
	}
      }
    }

    if(!hist_proc)
      continue;

    hist_proc->Rebin(rebin);
    
    if(isbkg){
      hists.push_back(hist_proc);
      processes += procs[i];
      if(hist_totbkg)
	hist_totbkg->Add(hist_proc);
      else {
	hist_totbkg = (TH1D*) hist_proc->Clone((can_name+"_totalbkg"+procs[i]).c_str());
      }
    } else {
      hist_data = hist_proc;
    }
  }

  if(hists.size() < 1)
    return nullptr;

  double scale_bkg = 1.;
  if(hist_data && hist_totbkg){
    b_ratio = true;
    if(norm){
      scale_bkg = hist_data->Integral() / hist_totbkg->Integral();
      hist_totbkg->Scale(scale_bkg);
    }
  }

  VS             vprocs;
  vector<TH1D*>  vhists;
  vector<double> vtotal;

  string stemp;
  TH1D*  htemp;
  double ttemp;
  
  int Nhist = hists.size();

  for(int i = 0; i < Nhist; i++){
    vprocs += processes[i];
    int Nbin = hists[i]->GetNbinsX();
    for(int j = 0; j < Nbin; j++){
      if(hists[i]->GetBinContent(j+1) < 0.)
    	hists[i]->SetBinContent(j+1, 0.);
    }
    if(norm)
      hists[i]->Scale(scale_bkg);
    vhists.push_back(hists[i]);
    vtotal.push_back(hists[i]->Integral());
    for(int j = vtotal.size()-2; j >= 0; j--){
      if(vtotal[j] < vtotal[j+1]){
	stemp = vprocs[j+1];
	htemp = vhists[j+1];
	ttemp = vtotal[j+1];
	vprocs[j+1] = vprocs[j];
	vhists[j+1]  = vhists[j];
	vtotal[j+1]  = vtotal[j];
	vprocs[j] = stemp;
	vhists[j] = htemp;
	vtotal[j] = ttemp;
      } else {
	break;
      }
    }
  }

  for(int i = 0; i < Nhist-1; i++){
    for(int j = i+1; j < Nhist; j++)
      vhists[i]->Add(vhists[j]);
  }
    
  hists = vhists;

  VS          labels;
  vector<int> colors;

  for(auto p : vprocs){
    if(m_Title.count(p) > 0)
      labels += m_Title[p];
    else
      labels += p;
    if(m_Color.count(p) > 0)
      colors.push_back(m_Color[p]);
    else
      colors.push_back(kBlack);
  }
  
  TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
			     Form("can_%s", can_name.c_str()),
			     1000, 700);

  double hlo = 0.105;
  double hhi = 0.22;
  double hbo = 0.19;
  double hto = 0.07;

  double ratio_h = 0.18;
  
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridy();
  //can->SetLogy();
  can->Draw();
  can->cd();

  TPad* pad = new TPad(Form("pad_%s", can_name.c_str()),
  		       Form("pad_%s", can_name.c_str()),
		       0.0, 0.0, 1., 1.);

  if(b_ratio)
    pad->SetBottomMargin(hbo+ratio_h);
  else
    pad->SetBottomMargin(hbo);
  pad->SetLeftMargin(hlo);
  pad->SetRightMargin(hhi);
  pad->SetTopMargin(hto);
  pad->SetGridy();
  pad->SetGridx();
  //pad->SetLogy();
  pad->Draw();
  pad->cd();

  string varlabel = var;
  if(m_Title.count(var) > 0)
    varlabel = m_Title[var];
  
  double hmax = hists[0]->GetMaximum();
  double hmin = std::min(0.05, hists[0]->GetMinimum()/2.);
  if(hmin <= 0.)
    hmin = 0.05;
  
  hists[0]->Draw("hist");
  hists[0]->GetXaxis()->SetNdivisions(10, 0, 0);
  hists[0]->GetXaxis()->CenterTitle();
  hists[0]->GetXaxis()->SetTitleFont(42);
  hists[0]->GetXaxis()->SetTitleSize(0.045);
  hists[0]->GetXaxis()->SetTitleOffset(1.0);
  hists[0]->GetXaxis()->SetLabelFont(42);
  hists[0]->GetXaxis()->SetLabelSize(0.035);
  hists[0]->GetXaxis()->SetTitle(varlabel.c_str());
  if(b_ratio){
    hists[0]->GetXaxis()->SetLabelSize(0.0);
    hists[0]->GetXaxis()->SetTitle("");
  }
  hists[0]->GetXaxis()->SetTickSize(0.);
  hists[0]->GetYaxis()->CenterTitle();
  hists[0]->GetYaxis()->SetTitleFont(42);
  hists[0]->GetYaxis()->SetTitleSize(0.04);
  hists[0]->GetYaxis()->SetTitleOffset(1.1);
  hists[0]->GetYaxis()->SetLabelFont(42);
  hists[0]->GetYaxis()->SetLabelSize(0.035);
  hists[0]->GetYaxis()->SetTitle("number of events");
  hists[0]->GetYaxis()->SetRangeUser(hmin, hmax*5);

  if(xmin >= 0. && xmax > xmin)
    hists[0]->GetXaxis()->SetRangeUser(xmin, xmax);
  pad->SetLogy();
   
  for(int i = 0; i < Nhist; i++){
    hists[i]->SetLineColor(kBlack);
    hists[i]->SetLineWidth(1.0);
    hists[i]->SetFillColor(colors[i]);
    hists[i]->SetFillStyle(1001);
    hists[i]->Draw("SAME HIST");
  }

  TH1D* htot = nullptr;
  if(hist_totbkg)
    htot = hist_totbkg;
  else
    htot = hists[0];

  int Nbin = htot->GetNbinsX();
  
  vector<double> X;
  vector<double> Xerr;
  vector<double> Y;
  vector<double> Yerr;
  vector<double> Y_bkg_ratio;
  vector<double> Yerr_bkg_ratio;
  for(int i = 0; i < Nbin; i++){
    X.push_back(htot->GetXaxis()->GetBinCenter(i+1));
    Xerr.push_back(htot->GetXaxis()->GetBinWidth(i+1)/2.);
    Y.push_back(htot->GetBinContent(i+1));
    Yerr.push_back(htot->GetBinError(i+1));
    Y_bkg_ratio.push_back(1.);
    if(htot->GetBinContent(i+1) > 0.)
      Yerr_bkg_ratio.push_back(htot->GetBinError(i+1)/htot->GetBinContent(i+1));
    else
      Yerr_bkg_ratio.push_back(0.);
  }

  TGraphErrors* gr = (TGraphErrors*) new TGraphErrors(Nbin, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);  
  gr->SetMarkerSize(0);
  gr->SetLineColor(kBlack);
  gr->SetFillColor(kBlack);
  gr->SetFillStyle(3244);
  gr->Draw("same p2");

  TGraphErrors* gr_bkg_ratio = (TGraphErrors*) new TGraphErrors(Nbin, &X[0],
								&Y_bkg_ratio[0],  &Xerr[0], &Yerr_bkg_ratio[0]);  
  gr_bkg_ratio->SetMarkerSize(0);
  gr_bkg_ratio->SetLineColor(kBlack);
  gr_bkg_ratio->SetFillColor(kBlack);
  gr_bkg_ratio->SetFillStyle(3244);

  TH1D* hist_data_ratio = nullptr;

  if(hist_data){
    hist_data->SetLineColor(kBlack);
    hist_data->SetFillColor(kWhite);
    hist_data->SetMarkerStyle(8);
    hist_data->SetMarkerSize(1.);
    hist_data->SetLineWidth(2);
    hist_data->Draw("SAME ep");
    if(hist_data->GetMaximum() > hmax)
      hmax = hist_data->GetMaximum();

    hist_data_ratio = (TH1D*) hist_data->Clone(Form("data_ratio_%s", can_name.c_str()));
    for(int i = 0; i < Nbin; i++){
      if(htot->GetBinContent(i+1) > 0.) {
	hist_data_ratio->SetBinContent(i+1, hist_data_ratio->GetBinContent(i+1)/htot->GetBinContent(i+1));
	hist_data_ratio->SetBinError(i+1, hist_data_ratio->GetBinError(i+1)/htot->GetBinContent(i+1));
      } else {
	hist_data_ratio->SetBinContent(i+1, 0.);
	hist_data_ratio->SetBinError(i+1, 0.);
      }
    }
  }

  double eps = 0.0015;
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  l.SetTextSize(std::min(0.03, 1.5/double(Nbin)));
  l.SetTextFont(42);
  l.SetTextAlign(32);
  l.SetTextAngle(90);

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nhist+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  if(hist_data)
    leg->AddEntry(hist_data, "data");
  leg->AddEntry(gr, "MC uncertainty","F");
  for(int i = 0; i < Nhist; i++)
    leg->AddEntry(hists[i], labels[i].c_str(), "F");
  leg->Draw("SAME");
  
  TPad* pad_ratio = nullptr;
  if(b_ratio){
    can->cd();
    pad_ratio = new TPad(Form("pad_%s", can_name.c_str()),
			 Form("pad_%s", can_name.c_str()),
			 0.0, 0., 1.-hhi, hbo+ratio_h);

    pad_ratio->SetBottomMargin(0.4);
    pad_ratio->SetLeftMargin(hlo*1.26);
    pad_ratio->SetRightMargin(0.);
    pad_ratio->SetTopMargin(0.05);
    pad_ratio->SetGridy();
    pad_ratio->SetGridx();
    pad_ratio->Draw();
    pad_ratio->cd();

    hist_data_ratio->Draw("ep");
    hist_data_ratio->GetXaxis()->SetNdivisions(10, 0, 0);
    hist_data_ratio->GetXaxis()->CenterTitle();
    hist_data_ratio->GetXaxis()->SetTitleFont(42);
    hist_data_ratio->GetXaxis()->SetTitleSize(0.11);
    hist_data_ratio->GetXaxis()->SetTitleOffset(1.47);
    hist_data_ratio->GetXaxis()->SetLabelFont(42);
    hist_data_ratio->GetXaxis()->SetLabelSize(0.095);
    hist_data_ratio->GetXaxis()->SetTitle(varlabel.c_str());
    hist_data_ratio->GetYaxis()->SetNdivisions(6, 5, 0);
    hist_data_ratio->GetYaxis()->CenterTitle();
    hist_data_ratio->GetYaxis()->SetTitleFont(42);
    hist_data_ratio->GetYaxis()->SetTitleSize(0.11);
    hist_data_ratio->GetYaxis()->SetTitleOffset(0.45);
    hist_data_ratio->GetYaxis()->SetLabelFont(42);
    hist_data_ratio->GetYaxis()->SetLabelSize(0.095);
    hist_data_ratio->GetYaxis()->SetTitle("#frac{data}{MC}");

    if(hist_data_ratio->GetMaximum() > 2.)
      hist_data_ratio->GetYaxis()->SetRangeUser(0., 2.8);

    if(xmin >= 0. && xmax > xmin)
       hist_data_ratio->GetXaxis()->SetRangeUser(xmin, xmax);
    
    gr_bkg_ratio->Draw("same p2");
    hist_data_ratio->Draw("same ep");
    
  }

  pad->cd();

  string reglabel = region;
  if(m_Title.count(region) > 0)
    reglabel = m_Title[region];
  
  l.SetTextAngle(0);
  l.SetTextColor(kBlack);
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, m_PlotTitle.c_str());
  l.DrawLatex(1.-hhi-.05, 1.-hto-0.05, string("Regions "+reglabel).c_str());

  
  l.SetTextAlign(11);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, m_CMSLabel.c_str());
  
  //can->Update();
  return can;
}






