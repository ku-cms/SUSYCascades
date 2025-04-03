#include <iostream>
#include <numeric>
#include <map>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TMultiGraph.h>

#include "FitPlotter.hh"
#include "shapeComparison.hh"
#include "plotShapeComparison.hh"

#include "RestFrames/RestFrames.hh"

///////////////////////////////////////////
////////// FitPlotter class
///////////////////////////////////////////

FitPlotter::FitPlotter(const string& inputfile,
		       const string& otherfile,
		       const string& otherfold)
  : FitReader(inputfile, otherfile, otherfold) {
  InitializeRecipes();
}

//FitPlotter::FitPlotter(){
//  InitializeRecipes();
//}

FitPlotter::~FitPlotter(){
  
}

TGraphErrors* FitPlotter::GetTotalBackground(const CategoryList& cat){
  TH1D* hist = nullptr;
  int Ncat = cat.GetN();
  for(int i = 0; i < Ncat; i++){
    string shist = m_FileFold+cat[i].Label()+"_"+cat[i].GetLabel()+"/total_background";
    //cout << shist << endl;
    if(hist == nullptr)
      hist = (TH1D*) m_FilePtr->Get((m_FileFold+"/"+shist).c_str())->Clone((shist+"_total").c_str());
    else
      hist->Add((TH1D*) m_FilePtr->Get((m_FileFold+"/"+shist).c_str()));
  }

  int NB = hist->GetNbinsX();
  
  vector<double> X;
  vector<double> Xerr;
  vector<double> Y;
  vector<double> Yerr;
  for(int i = 0; i < NB; i++){
    cout << "bin #: " << i+1 << " total_background: " << hist->GetBinContent(i+1) << endl;
	//X.push_back(0.5 + i);
    X.push_back(hist->GetXaxis()->GetBinCenter(i+1));
    Xerr.push_back(0.5);
    Y.push_back(hist->GetBinContent(i+1));
    Yerr.push_back(hist->GetBinError(i+1));
  }
   
  TGraphErrors* gr = new TGraphErrors(NB, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);

  return gr;
   
}

string FitPlotter::GetSignalTitle(const string& label){
  size_t p = label.rfind("_");
  if(p == std::string::npos)
    return label;

  string title = label.substr(0, p);
  int    mass  = stoi(label.substr(p+1,label.length()-p));
  return title+" "+std::to_string(mass/10000)+" "+std::to_string(mass%100000);
}



TCanvas* FitPlotter::Plot1DShape(const string& proc,
        const VS& lep_cat,
        const VS& hadS_cat,
        const VS& hadI_cat,
	const string& name){
  RestFrames::SetStyle();

  
  // int Nproc = proc.size();
  int Nlep  = lep_cat.size();
  int NhadS = hadS_cat.size();
  int NhadI = hadI_cat.size();
  if(Nlep  == 0 ||
     NhadS == 0 ||
     NhadI == 0)
    return nullptr;

  CategoryList cat = GetCategories();
  // cat.Print();
  
  // Leptonic
  VS lep_labels;
  VS vlep;
  for(int i = 0; i < Nlep; i++){
    if(m_Title.count(lep_cat[i]) != 0)
      lep_labels.push_back(m_Title[lep_cat[i]]);
    else
      lep_labels.push_back(lep_cat[i]);

    if(m_Strings.count(lep_cat[i]) != 0){
      int N = m_Strings[lep_cat[i]].size();
      for(int j = 0; j < N; j++){
  vlep.push_back(m_Strings[lep_cat[i]][j]);
      }
    } else {
      vlep.push_back(lep_cat[i]);
    }
  }

  cat = cat.FilterOR(vlep);

  // Hadronic S
  VS hadS_labels;
  VS vhadS;
  for(int i = 0; i < NhadS; i++){
    if(m_Title.count(hadS_cat[i]) != 0)
      hadS_labels.push_back(m_Title[hadS_cat[i]]);
    else
      hadS_labels.push_back(hadS_cat[i]);

    if(m_Strings.count(hadS_cat[i]) != 0){
      int N = m_Strings[hadS_cat[i]].size();
      for(int j = 0; j < N; j++)
  vhadS.push_back(m_Strings[hadS_cat[i]][j]);
    } else {
      vhadS.push_back(hadS_cat[i]);
    }
  }

  cat = cat.FilterOR(vhadS);

  // Hadronic ISR
  VS hadI_labels;
  VS vhadI;
  for(int i = 0; i < NhadI; i++){
    if(m_Title.count(hadI_cat[i]) != 0)
      hadI_labels.push_back(m_Title[hadI_cat[i]]);
    else
      hadI_labels.push_back(hadI_cat[i]);

    if(m_Strings.count(hadI_cat[i]) != 0){
      int N = m_Strings[hadI_cat[i]].size();
      for(int j = 0; j < N; j++)
  vhadI.push_back(m_Strings[hadI_cat[i]][j]);
    } else {
      vhadI.push_back(hadI_cat[i]);
    }
  }

  cat = cat.FilterOR(vhadI);

  int Ncat = cat.GetN();
  
  if(Ncat < 1)
    return nullptr;
  
  // Processes
  VS            labels;
  vector<int>   colors;
  // TH1D* hist;
  // vector<TH1D*> hists;
  
  // VS             labels_sig;
  // vector<TH1D*>  hists_sig;

  // TH1D* hist_data = nullptr;
  
  // for(int i = 0; i < Nproc; i++){
     VS vproc;
    if(m_Strings.count(proc) != 0)
      vproc = m_Strings[proc];
    else
      vproc += proc;

    ProcessType type = kBkg;
    TH1D*       hist = nullptr;
    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
  continue;
      
      Process pp = GetProcesses()[index];

      if(pp.Type() == kSig)
  type = kSig;
      if(pp.Type() == kData)
  type = kData;
      
      for(int c = 0; c < Ncat; c++){
  //cout << cat[c].GetLabel() << " " << pp.Name() << endl;
  if(!IsFilled(cat[c], pp))
    continue;

  //cout << "filled " << cat[c].GetLabel() << " " << pp.Name() << endl;
  
  if(!hist){
    hist = (TH1D*) GetHistogram(cat[c], pp)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
  } else {
    hist->Add(GetHistogram(cat[c], pp));
  }
      }
    }

    if(hist == nullptr){
      cout << "null hist" << endl;
      return NULL;
    }
    
    // if(type == kData){
      // hist_data = hist;
    // }

    // if(type == kSig){
      // labels_sig.push_back(GetSignalTitle(proc[i]));
      
      // hists_sig.push_back(hist);
    // }

    if(type == kBkg){
      if(m_Title.count(proc) != 0)
  labels += m_Title[proc];
      else
  labels += proc;
      
      if(m_Color.count(proc) != 0)
  colors.push_back(m_Color[proc]);
      else
  colors.push_back(m_ColorDefault[0]);
      
      // hists.push_back(hist);
    } 
  // }
  //normalize histogram
  hist->Scale(1/hist->Integral());

  for(int i=0; i < hist->GetNbinsX(); i++){
	cout << "Bin #" << i+1 << ": " << hist->GetBinContent(i+1) << endl;
}

  const FitBin& bin = cat[0].GetFitBin();

  int NR = bin.NRBins();
  int NB = bin.NBins();
  VS blabels;
  for(int r = 0; r < NR; r++)
    blabels += bin[r].GetMBinLabels();

  int lmax = 0;
  for(int b = 0; b < NB; b++){
    int len = blabels[b].length();
    if(blabels[b].find("#infty") != std::string::npos)
      len -= 5;
    if(len > lmax)
      lmax = len;
  }
  string space = "";
  for(int l = 0; l < 1.2*lmax; l++)
    space += " ";
  
  for(int b = 0; b < NB; b++){
    if(b%2 == 1)
      hist->GetXaxis()->SetBinLabel(b+1, (blabels[b]+space).c_str());
    else
      hist->GetXaxis()->SetBinLabel(b+1, blabels[b].c_str());
  }
  
  blabels.clear();

  for(int r = 0; r < NR; r++)
    blabels += bin[r].GetRBinLabel();
  
  hist->LabelsOption("v","X");
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
           Form("can_%s", name.c_str()),
           1200, 700);
  double hlo = 0.09;
  double hhi = 0.22;
  double hbo = 0.27;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridy();
  can->Draw();
  can->cd();

  double hmax = hist->GetMaximum();
  
  hist->Draw("p");
  hist->SetMarkerSize(3);
  hist->SetMarkerStyle(8);
  hist->SetMarkerColor(kRed);
  hist->SetLineColor(kRed);
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTitle("");
  hist->GetXaxis()->SetTickSize(0.);
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleOffset(0.85);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitle("a.u.");
   
  // for(int i = 0; i < Nbkg; i++){
  //   hists[i]->SetLineColor(kBlack);
  //   hists[i]->SetLineWidth(1.0);
  //   hists[i]->SetFillColor(colors[i]);
  //   hists[i]->SetFillStyle(1001);
  //   hists[i]->Draw("SAME HIST");
  // }

  //TGraphErrors* gr = nullptr;
  //if(!m_FilePtr){
  //  vector<double> X;
  //  vector<double> Xerr;
  //  vector<double> Y;
  //  vector<double> Yerr;
  //  for(int i = 0; i < NB; i++){
  //    X.push_back(hist->GetXaxis()->GetBinCenter(i+1));
  //    Xerr.push_back(0.5);
  //    Y.push_back(hist->GetBinContent(i+1));
  //    Yerr.push_back(hist->GetBinError(i+1));
  //  }
  //  gr = (TGraphErrors*) new TGraphErrors(NB, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);
  //} else {
  //  cout << "here " << gr << endl;
  //  gr = (TGraphErrors*) GetTotalBackground(cat);
  //  cout << "here " << gr << endl;
  //}
  //  
  //gr->SetMarkerSize(0);
  //gr->SetLineColor(kBlack);
  //gr->SetFillColor(kBlack);
  //gr->SetFillStyle(3244);
  //gr->Draw("same p2");

  // for(int i = 0; i < Nsig; i++){
  //   hists_sig[i]->SetLineColor(7030+i*10);
  //   hists_sig[i]->SetLineWidth(5);
  //   hists_sig[i]->SetFillColor(kWhite);
  //   hists_sig[i]->Draw("SAME HIST");
  //   if(hists_sig[i]->GetMaximum() > hmax)
  //     hmax = hists_sig[i]->GetMaximum();
  // }

  // if(hist_data){
  //   hist_data->SetLineColor(kBlack);
  //   hist_data->SetFillColor(kWhite);
  //   hist_data->SetMarkerStyle(8);
  //   hist_data->SetMarkerSize(1.);
  //   hist_data->SetLineWidth(2);
  //   hist_data->Draw("SAME ep");
  //   if(hist_data->GetMaximum() > hmax)
  //     hmax = hist_data->GetMaximum();
  // }

  hist->GetYaxis()->SetRangeUser(0.05, 1.1*hmax);

//   TLegend* leg = new TLegend(1.-hhi+0.01, 1.- (2)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
//   leg->SetName("legend");
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.035);
//   leg->SetFillColor(kWhite);
//   leg->SetLineColor(kWhite);
//   leg->SetShadowColor(kWhite);

  // if(hist_data)
  //   leg->AddEntry(hist_data, "data");
  // leg->AddEntry(gr, "total uncertainty","F");
  // for(int i = 0; i < Nbkg; i++)
//     leg->AddEntry(hist, labels[0].c_str(), "F");
 //  for(int i = 0; i < Nsig; i++)
 //    leg->AddEntry(hists_sig[i], labels_sig[i].c_str(), "L");
  // leg->Draw("SAME");


//CANVAS FORMATTING
   double eps = 0.0015;
  
   TLatex l;
   l.SetTextFont(42);
   l.SetNDC();

   TLine* line = new TLine();
   line->SetLineWidth(2);
   line->SetLineColor(kBlack);

   // line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 
   l.SetTextSize(0.025);
   l.SetTextFont(42);
   l.SetTextAlign(23);
   line->SetLineWidth(1);
   double lo = hlo;
   double hi = hlo;
   double yline = hbo-0.024*lmax;
   int ib = 0;
   for(int r = 0; r < NR; r++){
     int NM = bin[r].NBins();
     lo = hi;
     hi = double(NM)/double(NB)*(1.-hhi-hlo) + lo;
  
     line->SetLineStyle(1);
     line->DrawLineNDC(lo + eps, yline,
             lo + eps, yline + 6*eps);
     line->DrawLineNDC(hi - eps, yline,
             hi - eps, yline + 6*eps);
     line->DrawLineNDC(lo + eps, yline,
             hi - eps, yline);
     line->SetLineStyle(5);
     line->DrawLineNDC(hi, hbo, hi, 1.-hto);
     line->SetLineStyle(3);
     for(int b = 0; b < NM; b++){
       if(ib%2 == 1)
     line->DrawLineNDC(lo + (hi-lo)*(b+0.5)/double(NM), hbo,
           lo + (hi-lo)*(b+0.5)/double(NM), (hbo+yline)/2.+eps);
       ib++;
     }
  
     l.DrawLatex((hi+lo)/2., yline - 8*eps, blabels[r].c_str());
   }
     
   l.SetTextAlign(32);
   l.SetTextSize(0.03);
   l.SetTextFont(42);
   l.DrawLatex(hlo, (hbo+yline)/2.+eps, "M_{#perp}   [GeV] #in");

   l.SetTextSize(0.03);
   l.SetTextFont(42);
   l.DrawLatex(hlo, yline - 15*eps, "#scale[1.15]{R_{ISR}} #in");
 
   l.SetTextAlign(31);
   l.SetTextSize(0.04);
   l.SetTextFont(42);
   l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, "2017 MC KUEWKino");
   l.SetTextAlign(11);
   l.SetTextSize(0.04);
   l.SetTextFont(42);
   l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
   l.SetTextSize(0.05);

   string plotlabel = "#scale[1.2]{"+labels[0]+"}"; 
   plotlabel += " #color[7014]{"+lep_labels[0]+"} + ";
   plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
   plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
   plotlabel += "p_{T}^{ISR} > 300 GeV";
  
   l.SetTextColor(kBlack);
   l.SetTextAlign(13);
   l.SetTextSize(0.035);
   l.SetTextFont(42);
   l.DrawLatex(hlo+0.02, 1-hto-0.012, plotlabel.c_str());
  
  return can;
  
}



TCanvas* FitPlotter::Plot1Dstack(const VS& proc,
				const VS& lep_cat,
				const VS& hadS_cat,
				const VS& hadI_cat,
				const string& name,
				const VS& extra){

  RestFrames::SetStyle();
  
  int Nproc = proc.size();
  int Nlep  = lep_cat.size();
  int NhadS = hadS_cat.size();
  int NhadI = hadI_cat.size();
  int Nextra = -999;

  if(!extra.empty()) 
    Nextra = extra.size();

  if(Nproc == 0 ||
     Nlep  == 0 ||
     NhadS == 0 ||
     NhadI == 0)
    return nullptr;
  
  CategoryList cat = GetCategories();

  // Leptonic
  VS lep_labels = FetchLabels(lep_cat);
  VS vlep = FetchCats(lep_cat);
  cat = cat.FilterOR(vlep);

  // Hadronic S
  VS hadS_labels = FetchLabels(hadS_cat);
  VS vhadS = FetchCats(hadS_cat);
  cat = cat.FilterOR(vhadS);
  
  // Hadronic ISR
  VS hadI_labels = FetchLabels(hadI_cat);
  VS vhadI = FetchCats(hadI_cat);
  cat = cat.FilterOR(vhadI);

  //Extra
  VS extra_labels = FetchLabels(extra);
  VS vextra = FetchCats(extra);
  cat = cat.FilterOR(vextra);

  int Ncat = cat.GetN();
  
  if(Ncat == 0){
    cout << "no categories found" << endl;
    return nullptr;
  }
  /*
  VS proc_bkg;
  VS proc_sig;
  VS proc_data;

  for(int i = 0; i < Nproc; i++){
    VS vproc;

    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    ProcessType type = kBkg;

    for(int p = 0; p < int(vproc.size()); p++){
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
        continue;

      Process pp = GetProcesses()[index];

      if(pp.Type() == kSig)
	proc_sig += vproc[p];

      if(pp.Type() == kBkg)
        proc_bkg += vproc[p];

      if(pp.Type() == kData)
        proc_data += vproc[p];

    }
  }

  cout << "signals: " << endl;
  for(int i = 0; i < proc_sig.size(); i++)
    cout << proc_sig[i] << endl;
  
  cout << "backgrounds: " << endl;
  for(int i = 0; i < proc_bkg.size(); i++)
    cout << proc_bkg[i] << endl;
  
  cout << "data: " << endl;
  for(int i = 0; i < proc_data.size(); i++)
    cout << proc_data[i] << endl;
  */
  // Processes
  VS            labels;
  vector<int>   colors;
  vector<TH1D*> hists;
  
  VS             labels_sig;
  vector<TH1D*>  hists_sig;

  TH1D* hist_data = nullptr;
  
  for(int i = 0; i < Nproc; i++){
    VS vproc;
    
    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    ProcessType type = kBkg;
    TH1D*       hist = nullptr;
    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;
      
      Process pp = GetProcesses()[index];
      
      if(pp.Type() == kSig){
	type = kSig;
	if(m_FileFold != nullptr){
	  string name = pp.Name();
	  name = name.substr(0,name.find("_")+1);
	  pp = Process(name, kSig);
	}
      }
      if(pp.Type() == kData)
	type = kData;
      
      for(int c = 0; c < Ncat; c++){
	if(!IsFilled(cat[c], pp)){
	  continue;
	}
	if(!hist){
	  hist = (TH1D*) GetHistogram(cat[c], pp)->Clone(Form("plothist_%d_%s", i, name.c_str()));
	} else {
	  hist->Add(GetHistogram(cat[c], pp));
	}
      }
    }

    if(hist == nullptr)
      continue;
    
    if(type == kData){
      hist_data = hist;
    }
    
    if(type == kSig){
      labels_sig += GetSignalTitle(proc[i]);      
      hists_sig.push_back(hist);
    }

    if(type == kBkg){
      if(m_Title.count(proc[i]) != 0)
	labels += m_Title[proc[i]];
      else
	labels += proc[i];
      
      if(m_Color.count(proc[i]) != 0)
	colors.push_back(m_Color[proc[i]]);
      else
	colors.push_back(m_ColorDefault[i]);
      
      hists.push_back(hist);
    } 
  }

  if(hists.size() < 1) return nullptr;
  int Nsig = hists_sig.size();
  
  // sort the histograms by integral (N^2/2 brute force)
  int Nbkg = hists.size();
  VS             vlabels;
  vector<int>    vcolors;
  vector<TH1D*>  vhists;
  string stemp;
  int    itemp;
  TH1D*  htemp;
  int nBins = hists[0]->GetNbinsX();
  
  for(int i = 0; i < Nbkg; i++){
    vlabels.push_back(labels[i]);
    vcolors.push_back(colors[i]);
    vhists.push_back(hists[i]);
    for(int j = vhists.size()-2; j >= 0; j--){
      if(vhists[j]->Integral() < vhists[j+1]->Integral()){
	stemp = vlabels[j+1];
	itemp = vcolors[j+1];
	htemp = vhists[j+1];
	vlabels[j+1] = vlabels[j];
	vcolors[j+1] = vcolors[j];
	vhists[j+1]  = vhists[j];
	vlabels[j] = stemp;
	vcolors[j] = itemp;
	vhists[j]  = htemp;
      } else {
	break;
      }
    }
  }
  
  //"stack" the histograms by adding
  for(int i = Nbkg-2; i >= 0; i--)
    vhists[i]->Add(vhists[i+1]);
  
  hists  = vhists;
  labels = vlabels;
  colors = vcolors;
  
  const FitBin& bin = cat[0].GetFitBin();
  
  int NR = bin.NRBins();
  int NB = bin.NBins();
  VS blabels;
  for(int r = 0; r < NR; r++)
    blabels += bin[r].GetMBinLabels();

  int lmax = 0;
  for(int b = 0; b < NB; b++){
    int len = blabels[b].length();
    if(blabels[b].find("#infty") != std::string::npos)
      len -= 5;
    if(len > lmax)
      lmax = len;
  }
  
  string space = "";
  for(int l = 0; l < 1.6*lmax; l++)
    space += " ";
  
  for(int b = 0; b < NB; b++){
    if(b%2 == 1)
      hists[0]->GetXaxis()->SetBinLabel(b+1, (blabels[b]+space).c_str());
    else
      hists[0]->GetXaxis()->SetBinLabel(b+1, blabels[b].c_str());
  }
  
  blabels.clear();

  for(int r = 0; r < NR; r++)
    blabels += bin[r].GetRBinLabel();
  
  hists[0]->LabelsOption("v","X");
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
			     Form("can_%s", name.c_str()),
			     1200, 700);
  double hlo = 0.09;
  double hhi = 0.22;
  double hbo = 0.27;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridy();
  can->Draw();
  can->cd();
  
  double hmax = hists[0]->GetMaximum();
  
  hists[0]->SetTitle("");
  hists[0]->Draw("hist");
  hists[0]->GetXaxis()->CenterTitle();
  hists[0]->GetXaxis()->SetTitleFont(42);
  hists[0]->GetXaxis()->SetTitleSize(0.05);
  hists[0]->GetXaxis()->SetTitleOffset(1.0);
  hists[0]->GetXaxis()->SetLabelFont(42);
  hists[0]->GetXaxis()->SetLabelSize(0.04);
  hists[0]->GetXaxis()->SetTitle("");
  hists[0]->GetXaxis()->SetTickSize(0.);
  hists[0]->GetYaxis()->CenterTitle();
  hists[0]->GetYaxis()->SetTitleFont(42);
  hists[0]->GetYaxis()->SetTitleSize(0.04);
  hists[0]->GetYaxis()->SetTitleOffset(0.85);
  hists[0]->GetYaxis()->SetLabelFont(42);
  hists[0]->GetYaxis()->SetLabelSize(0.035);
  hists[0]->GetYaxis()->SetTitle("number of events");
  for(int i = 0; i < Nbkg; i++){
    hists[i]->SetLineColor(kBlack);
    hists[i]->SetLineWidth(1.0);
    hists[i]->SetFillColor(colors[i]);
    hists[i]->SetFillStyle(1001);
    hists[i]->Draw("SAME HIST");
  }

  TGraphErrors* gr = nullptr;
  if(!m_FilePtr){
    vector<double> X;
    vector<double> Xerr;
    vector<double> Y;
    vector<double> Yerr;
    for(int i = 0; i < NB; i++){
      X.push_back(hists[0]->GetXaxis()->GetBinCenter(i+1));
      Xerr.push_back(0.5);
      Y.push_back(hists[0]->GetBinContent(i+1));
      Yerr.push_back(hists[0]->GetBinError(i+1));
    }
    gr = (TGraphErrors*) new TGraphErrors(NB, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);
  } 
  else 
    gr = (TGraphErrors*) GetTotalBackground(cat);
    
  gr->SetMarkerSize(0);
  gr->SetLineColor(kBlack);
  gr->SetFillColor(kBlack);
  gr->SetFillStyle(3244);
  gr->Draw("same p2");

  for(int i = 0; i < Nsig; i++){
    hists_sig[i]->SetLineColor(7030+i*10);
    hists_sig[i]->SetLineWidth(5);
    hists_sig[i]->SetFillColor(kWhite);
    hists_sig[i]->Draw("SAME HIST");
    if(hists_sig[i]->GetMaximum() > hmax)
      hmax = hists_sig[i]->GetMaximum();
  }

  if(hist_data){
    hist_data->SetLineColor(kBlack);
    hist_data->SetFillColor(kWhite);
    hist_data->SetMarkerStyle(8);
    hist_data->SetMarkerSize(1.);
    hist_data->SetLineWidth(2);
    hist_data->Draw("SAME ep");
    if(hist_data->GetMaximum() > hmax)
      hmax = hist_data->GetMaximum();
  }

  hists[0]->GetYaxis()->SetRangeUser(0.05, 1.1*hmax);

  TLegend* leg = new TLegend(1.-hhi+0.01, 1.- (Nbkg+Nsig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetName("legend");
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  if(hist_data)
    leg->AddEntry(hist_data, "data");
  leg->AddEntry(gr, "total uncertainty","F");
  for(int i = 0; i < Nbkg; i++)
    leg->AddEntry(hists[i], labels[i].c_str(), "F");
  for(int i = 0; i < Nsig; i++)
    leg->AddEntry(hists_sig[i], labels_sig[i].c_str(), "L");
  leg->Draw("SAME");

  double eps = 0.0015;
  
  TLatex l;
  l.SetName("latex");
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);

  // line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 
  l.SetTextSize(0.025);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(1);
  double lo = hlo;
  double hi = hlo;
  double yline = hbo-0.024*lmax;
  int ib = 0;
  for(int r = 0; r < NR; r++){
    int NM = bin[r].NBins();
    lo = hi;
    hi = double(NM)/double(NB)*(1.-hhi-hlo) + lo;
    
    line->SetLineStyle(1);
    line->DrawLineNDC(lo + eps, yline,
  		      lo + eps, yline + 6*eps);
    line->DrawLineNDC(hi - eps, yline,
  		      hi - eps, yline + 6*eps);
    line->DrawLineNDC(lo + eps, yline,
  		      hi - eps, yline);
    line->SetLineStyle(5);
    line->DrawLineNDC(hi, hbo, hi, 1.-hto);
    line->SetLineStyle(3);
    for(int b = 0; b < NM; b++){
      if(ib%2 == 1)
  	line->DrawLineNDC(lo + (hi-lo)*(b+0.5)/double(NM), hbo,
  			  lo + (hi-lo)*(b+0.5)/double(NM), (hbo+yline)/2.+eps);
      ib++;
    }
    
    l.DrawLatex((hi+lo)/2., yline - 8*eps, blabels[r].c_str());
  }
       
  l.SetTextAlign(32);
  l.SetTextSize(0.03);
  l.SetTextFont(42);
  l.DrawLatex(hlo, (hbo+yline)/2.+eps, "M_{#perp}   [GeV] #in");

  l.SetTextSize(0.03);
  l.SetTextFont(42);
  l.DrawLatex(hlo, yline - 15*eps, "#scale[1.15]{R_{ISR}} #in");
 
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, "2017 MC KUEWKino");
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
  l.SetTextSize(0.05);

  string plotlabel = "#color[7014]{"+lep_labels[0]+"} + ";
  plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
  plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
  if(extra_labels[0] != "") plotlabel += "#color[7024]{"+extra_labels[0]+"} + ";
  plotlabel += "p_{T}^{ISR} > 300 GeV";
  
  l.SetTextColor(kBlack);
  l.SetTextAlign(13);
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(hlo+0.02, 1-hto-0.012, plotlabel.c_str());
  
  return can;
  
}

TCanvas* FitPlotter::Plot2D(const VS& proc,
			   const VS& lep_cat,
			   const VS& hadS_cat,
			   const VS& hadI_cat,
			   const string& name,
			   const string& extra){
  RestFrames::SetStyle();

  
  int Nproc = proc.size();
  int Nlep  = lep_cat.size();
  int NhadS = hadS_cat.size();
  int NhadI = hadI_cat.size();
  if(Nproc == 0 ||
     Nlep  == 0 ||
     NhadS == 0 ||
     NhadI == 0)
    return nullptr;

  CategoryList cat = GetCategories();
  //cat.Print();
  
  // Leptonic
  VS lep_labels;
  VS vlep;
  for(int i = 0; i < Nlep; i++){
    if(m_Title.count(lep_cat[i]) != 0)
      lep_labels.push_back(m_Title[lep_cat[i]]);
    else
      lep_labels.push_back(lep_cat[i]);

    if(m_Strings.count(lep_cat[i]) != 0){
      int N = m_Strings[lep_cat[i]].size();
      for(int j = 0; j < N; j++){
	vlep.push_back(m_Strings[lep_cat[i]][j]);
      }
    } else {
      vlep.push_back(lep_cat[i]);
    }
  }

  cat = cat.FilterOR(vlep);

  // Hadronic S
  VS hadS_labels;
  VS vhadS;
  for(int i = 0; i < NhadS; i++){
    if(m_Title.count(hadS_cat[i]) != 0)
      hadS_labels.push_back(m_Title[hadS_cat[i]]);
    else
      hadS_labels.push_back(hadS_cat[i]);

    if(m_Strings.count(hadS_cat[i]) != 0){
      int N = m_Strings[hadS_cat[i]].size();
      for(int j = 0; j < N; j++)
	vhadS.push_back(m_Strings[hadS_cat[i]][j]);
    } else {
      vhadS.push_back(hadS_cat[i]);
    }
  }

  cat = cat.FilterOR(vhadS);

  // Hadronic ISR
  VS hadI_labels;
  VS vhadI;
  for(int i = 0; i < NhadI; i++){
    if(m_Title.count(hadI_cat[i]) != 0)
      hadI_labels.push_back(m_Title[hadI_cat[i]]);
    else
      hadI_labels.push_back(hadI_cat[i]);

    if(m_Strings.count(hadI_cat[i]) != 0){
      int N = m_Strings[hadI_cat[i]].size();
      for(int j = 0; j < N; j++)
	vhadI.push_back(m_Strings[hadI_cat[i]][j]);
    } else {
      vhadI.push_back(hadI_cat[i]);
    }
  }

  cat = cat.FilterOR(vhadI);

  if(extra != "")
    cat = cat.Filter(extra);

  int Ncat = cat.GetN();
  
  if(Ncat < 1)
    return nullptr;
  
  // Processes
  string label;
  TH2D* hist = nullptr;
  
  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    ProcessType type = kBkg;
    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;
      
      Process pp = GetProcesses()[index];

      if(pp.Type() == kSig)
	type = kSig;
      if(pp.Type() == kData)
	type = kData;
      
      for(int c = 0; c < Ncat; c++){
	//cout << cat[c].GetLabel() << " " << pp.Name() << endl;
	if(!IsFilled(cat[c], pp))
	  continue;

	cout << "filled " << cat[c].GetLabel() << " " << pp.Name() << endl;
	
	if(!hist){
	  hist = (TH2D*) GetHistogram2D(cat[c], pp)->Clone(Form("plothist_%dhis_%s_2D", i, name.c_str()));
	} else {
	  hist->Add(GetHistogram2D(cat[c], pp));
	}
      }
    }

    if(hist == nullptr)
      continue;
    
    if(type == kData){
      label = "Data";
    }

    if(type == kSig){
      label = GetSignalTitle(proc[i]);
    }

    if(type == kBkg){
      if(m_Title.count(proc[i]) != 0)
	label = m_Title[proc[i]];
      else
	label = proc[i];
    } 
  }

  const FitBin& bin = cat[0].GetFitBin();

  int NR = bin.NRBins();
  int NB = bin.NBins();
  // VS blabels;
  // for(int r = 0; r < NR; r++)
  //   blabels += bin[r].GetMBinLabels();

  // int lmax = 0;
  // for(int b = 0; b < NB; b++){
  //   int len = blabels[b].length();
  //   if(blabels[b].find("#infty") != std::string::npos)
  //     len -= 5;
  //   if(len > lmax)
  //     lmax = len;
  // }
  // string space = "";
  // for(int l = 0; l < 1.6*lmax; l++)
  //   space += " ";
  
  // for(int b = 0; b < NB; b++){
  //   if(b%2 == 1)
  //     hists[0]->GetXaxis()->SetBinLabel(b+1, (blabels[b]+space).c_str());
  //   else
  //     hists[0]->GetXaxis()->SetBinLabel(b+1, blabels[b].c_str());
  // }
  
  // blabels.clear();

  // for(int r = 0; r < NR; r++)
  //   blabels += bin[r].GetRBinLabel();
  
  // hists[0]->LabelsOption("v","X");
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
			     Form("can_%s", name.c_str()),
			     700, 700);
  double eps = 0.0015;

  double hlo = 0.14;
  double hhi = 0.18;
  double hbo = 0.13;
  double hto = 0.17;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->Draw();
  can->cd();
  
  hist->Draw("colz");
  hist->GetXaxis()->SetTitle("M_{#perp}   [GeV]");
  hist->GetYaxis()->SetTitle("R_{ISR}");
  hist->GetZaxis()->SetTitle("number of events");
  
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTickSize(0.);
  
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTickSize(0.);
  
  hist->GetZaxis()->CenterTitle();
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleOffset(0.9);

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(7024);

  // line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 
  l.SetTextSize(0.025);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(3);
  line->SetLineColor(7024);
  line->SetLineStyle(7);
 
  for(int r = 0; r < NR; r++){
    line->DrawLine(hist->GetXaxis()->GetXmin(), bin[r].Rlow(),
		   hist->GetXaxis()->GetXmax(), bin[r].Rlow());
    line->DrawLine(hist->GetXaxis()->GetXmin(), bin[r].Rhigh(),
		   hist->GetXaxis()->GetXmax(), bin[r].Rhigh());

    for(int b = 0; b <  bin[r].BinEdges().size()-1; b++){
      line->DrawLine(bin[r].BinEdges()[b], bin[r].Rlow(),
		     bin[r].BinEdges()[b], bin[r].Rhigh());
    }
    
    // int NM = bin[r].NBins();
    // lo = hi;
    // hi = double(NM)/double(NB)*(1.-hhi-hlo) + lo;
    
    // line->SetLineStyle(1);
    // line->DrawLineNDC(lo + eps, yline,
    // 		      lo + eps, yline + 6*eps);
    // line->DrawLineNDC(hi - eps, yline,
    // 		      hi - eps, yline + 6*eps);
    // line->DrawLineNDC(lo + eps, yline,
    // 		      hi - eps, yline);
    // line->SetLineStyle(5);
    // line->DrawLineNDC(hi, hbo, hi, 1.-hto);
    // line->SetLineStyle(3);
    // for(int b = 1; b < NM; b++)
    //   if(b%2 == 1)
    // 	line->DrawLineNDC(lo + (hi-lo)*(b+0.5)/double(NM), hbo,
    // 			  lo + (hi-lo)*(b+0.5)/double(NM), (hbo+yline)/2.+eps);
    
    // l.DrawLatex((hi+lo)/2., yline - 8*eps, blabels[r].c_str());
  }
 
 
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  //l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.12, "2017 MC KUEWKino");
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, label.c_str());
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
  l.SetTextSize(0.05);

  string plotlabel = "#color[7014]{"+lep_labels[0]+"} + ";
  plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
  plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
  plotlabel += "p_{T}^{ISR} > 300 GeV";
  
  l.SetTextColor(kBlack);
  l.SetTextAlign(13);
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(hlo+0.02, 1-hto+0.15, plotlabel.c_str());
  
  
  return can;
  
}

TCanvas* FitPlotter::PlotYields(const string& can_name,
			       const VS& proc,
			       const CategoryTree& CT){

  RestFrames::SetStyle();

  int Nproc = proc.size();
  if(Nproc == 0)
    return nullptr;

  int Nvis = CT.GetNVisible();
  if(Nvis < 1)
    return nullptr;

  vector<const CategoryTree*> CatTrees;
  CT.GetListVisible(CatTrees);
  
  CategoryList CatList = GetCategories();

  if(CatList.GetN() < 1)
    return nullptr;

  // Processes
  VS                  labels;
  vector<int>         colors;
  vector<TH1D*>       hists[Nvis];
  vector<double>      total;
  
  VS                  labels_sig;
  vector<TH1D*>       hists_sig[Nvis];

  TH1D* hist_data[Nvis];
  double total_data = 0.;

  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    double itot = 0.;
    
    ProcessList procs;
    
    ProcessType type = kBkg;
    TH1D* hist[Nvis];
    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;
      
      Process pp = GetProcesses()[index];


      if(pp.Type() == kSig){
	type = kSig;
	if(m_FileFold != nullptr){
		string name = pp.Name();
		name = name.substr(0,name.find("_")+1);
		pp = Process(name,kSig);
	}
	}
      
      procs += pp;
      if(pp.Type() == kData){
	type = kData;
      }
    }

    for(int v = 0; v < Nvis; v++){
      CategoryList cat = CatList.Filter(*CatTrees[v]);
      TH1D* h = GetIntegralHist(Form("plothist_%d_%d_%s", i, v, can_name.c_str()), cat, procs);
      if(h){
	itot += h->Integral();
      }
	hist[v] = h;
    }
    
    if(itot <= 1e-4)
      continue;
    
    if(type == kData){
      for(int v = 0; v < Nvis; v++)
	hist_data[v] = hist[v];
      total_data = itot;
    }
    
    if(type == kSig){
      labels_sig.push_back(GetSignalTitle(proc[i]));
      for(int v = 0; v < Nvis; v++)
	hists_sig[v].push_back(hist[v]);
    }
    
    if(type == kBkg){
      if(m_Title.count(proc[i]) != 0)
	labels += m_Title[proc[i]];
      else
	labels += proc[i];
      
      if(m_Color.count(proc[i]) != 0)
	colors.push_back(m_Color[proc[i]]);
      else
	colors.push_back(m_ColorDefault[i]);
      
      for(int v = 0; v < Nvis; v++)
	hists[v].push_back(hist[v]);
      total.push_back(itot);
    } 
  }
  
  int Nsig = hists_sig[0].size();
  
  // sort the histograms by integral (N^2/2 brute force)
  int Nbkg = total.size();
  VS             vlabels;
  vector<int>    vcolors;
  vector<TH1D*>  vhists[Nvis];
  vector<double> vtotal;
  string stemp;
  int    itemp;
  TH1D*  htemp[Nvis];
  double ttemp;
  
  for(int i = 0; i < Nbkg; i++){
    vlabels.push_back(labels[i]);
    vcolors.push_back(colors[i]);
    for(int v = 0; v < Nvis; v++)
      vhists[v].push_back(hists[v][i]);
    vtotal.push_back(total[i]);
    for(int j = vtotal.size()-2; j >= 0; j--){
      if(vtotal[j] < vtotal[j+1]){
	stemp = vlabels[j+1];
	itemp = vcolors[j+1];
	for(int v = 0; v < Nvis; v++)
	  htemp[v] = vhists[v][j+1];
	ttemp = vtotal[j+1];
	vlabels[j+1] = vlabels[j];
	vcolors[j+1] = vcolors[j];
	for(int v = 0; v < Nvis; v++)
	  vhists[v][j+1]  = vhists[v][j];
	vtotal[j+1]  = vtotal[j];
	vlabels[j] = stemp;
	vcolors[j] = itemp;
	for(int v = 0; v < Nvis; v++)
	  vhists[v][j]  = htemp[v];
	vtotal[j] = ttemp;
      } else {
	break;
      }
    }
  }
  
  vector<TH1D*> fhists;
  vector<TH1D*> fhists_sig;
  TH1D*         fhist_data = nullptr;

  for(int i = 0; i < Nsig; i++){
    fhists_sig.push_back(new TH1D(Form("fhistsig_%d_%s", i, can_name.c_str()),
				  Form("fhistsig_%d_%s", i, can_name.c_str()),
				  Nvis, 0., Nvis));
    for(int v = 0; v < Nvis; v++){
      if(hists_sig[v][i]){
	fhists_sig[i]->SetBinContent(v+1, hists_sig[v][i]->GetBinContent(1));
	fhists_sig[i]->SetBinError(v+1, hists_sig[v][i]->GetBinError(1));
      }
    }
  }
 
  if(total_data > 0.){
    fhist_data = new TH1D(Form("fhistdata_%s", can_name.c_str()),
			  Form("fhistdata_%s", can_name.c_str()),
			  Nvis, 0., Nvis);
    for(int v = 0; v < Nvis; v++){
      if(hist_data[v]){
	fhist_data->SetBinContent(v+1, hist_data[v]->GetBinContent(1));
	fhist_data->SetBinError(v+1, hist_data[v]->GetBinError(1));
      }
    }
  }

  for(int i = 0; i < Nbkg; i++){
    for(int v = 0; v < Nvis; v++){
      if(vhists[v][i]){
	for(int j = i+1; j < Nbkg; j++)
	  if(vhists[v][j])
	    vhists[v][i]->Add(vhists[v][j]);
      }
    }
  }

  for(int i = 0; i < Nbkg; i++){
    fhists.push_back(new TH1D(Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Nvis, 0., Nvis));
    for(int v = 0; v < Nvis; v++){
      int j = i;
      TH1D* hptr = nullptr;
      while(j < Nbkg && hptr == nullptr){
	hptr = vhists[v][j];
	j++;
      }
      if(hptr){
	fhists[i]->SetBinContent(v+1, hptr->GetBinContent(1));
	fhists[i]->SetBinError(v+1, hptr->GetBinError(1));
      }
    }
  }

  labels = vlabels;
  colors = vcolors;

  for(int b = 0; b < Nvis; b++){
    fhists[0]->GetXaxis()->SetBinLabel(b+1, CatTrees[b]->GetSpectroscopicLabel().c_str());
  }
  
  double hlo = 0.09;
  double hhi = 0.2;
  double hbo = 0.19;
  double hto = 0.07;

  fhists[0]->LabelsOption("v","X");
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
			     Form("can_%s", can_name.c_str()),
			     1200, 700); 
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridy();
  can->SetLogy();
  can->Draw();
  can->cd();

  double hmax = fhists[0]->GetMaximum();
  double hmin = std::max(0.1, fhists[Nbkg-1]->GetMinimum());
  
  fhists[0]->Draw("hist");
  fhists[0]->GetXaxis()->CenterTitle();
  fhists[0]->GetXaxis()->SetTitleFont(42);
  fhists[0]->GetXaxis()->SetTitleSize(0.05);
  fhists[0]->GetXaxis()->SetTitleOffset(1.0);
  fhists[0]->GetXaxis()->SetLabelFont(42);
  fhists[0]->GetXaxis()->SetLabelSize(std::min(2./double(Nvis), 0.04));
  fhists[0]->GetXaxis()->SetTitle("");
  fhists[0]->GetXaxis()->SetTickSize(0.);
  fhists[0]->GetYaxis()->CenterTitle();
  fhists[0]->GetYaxis()->SetTitleFont(42);
  fhists[0]->GetYaxis()->SetTitleSize(0.04);
  fhists[0]->GetYaxis()->SetTitleOffset(0.85);
  fhists[0]->GetYaxis()->SetLabelFont(42);
  fhists[0]->GetYaxis()->SetLabelSize(0.035);
  fhists[0]->GetYaxis()->SetTitle("number of events");
   
  for(int i = 0; i < Nbkg; i++){
    fhists[i]->SetLineColor(kBlack);
    fhists[i]->SetLineWidth(1.0);
    fhists[i]->SetFillColor(colors[i]);
    fhists[i]->SetFillStyle(1001);
    fhists[i]->Draw("SAME HIST");
  }

  TGraphErrors* gr = nullptr;
  if(true){// !m_FilePtr
    vector<double> X;
    vector<double> Xerr;
    vector<double> Y;
    vector<double> Yerr;
    for(int i = 0; i < Nvis; i++){
      X.push_back(fhists[0]->GetXaxis()->GetBinCenter(i+1));
      Xerr.push_back(0.5);
      Y.push_back(fhists[0]->GetBinContent(i+1));
      Yerr.push_back(fhists[0]->GetBinError(i+1));
    }
    gr = (TGraphErrors*) new TGraphErrors(Nvis, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);
  } else {
    // cout << "here " << gr << endl;
    // gr = (TGraphErrors*) GetTotalBackground(cat);
    // cout << "here " << gr << endl;
  }
    
  gr->SetMarkerSize(0);
  gr->SetLineColor(kBlack);
  gr->SetFillColor(kBlack);
  gr->SetFillStyle(3244);
  gr->Draw("same p2");
  
  for(int i = 0; i < Nsig; i++){
    fhists_sig[i]->SetLineColor(m_SignalColor[i]);
    fhists_sig[i]->SetLineWidth(8);
    fhists_sig[i]->SetFillColor(kWhite);
    fhists_sig[i]->Draw("SAME HIST");
    if(fhists_sig[i]->GetMaximum() > hmax)
      hmax = fhists_sig[i]->GetMaximum();
  }

  if(fhist_data){
    fhist_data->SetLineColor(kBlack);
    fhist_data->SetFillColor(kWhite);
    fhist_data->SetMarkerStyle(8);
    fhist_data->SetMarkerSize(1.);
    fhist_data->SetLineWidth(2);
    fhist_data->Draw("SAME ep");
    if(fhist_data->GetMaximum() > hmax)
      hmax = fhist_data->GetMaximum();
  }

  double logrange = log(hmax) - log(hmin);
  hmax = exp(log(hmax) + logrange*0.3);
  fhists[0]->GetYaxis()->SetRangeUser(hmin/1.5, hmax);

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nbkg+Nsig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  if(fhist_data)
    leg->AddEntry(fhist_data, "data");
  leg->AddEntry(gr, "total uncertainty","F");
  for(int i = 0; i < Nbkg; i++)
    leg->AddEntry(fhists[i], labels[i].c_str(), "F");
  for(int i = 0; i < Nsig; i++)
    leg->AddEntry(fhists_sig[i], labels_sig[i].c_str(), "L");
  leg->Draw("SAME");

  DrawCatTree(CT, can);
  
  
  double eps = 0.0015;
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);

  // line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 
  l.SetTextSize(0.025);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(1);
  double lo = hlo;
  double hi = hlo;
 
 
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, string("Regions "+CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, m_CMSLabel.c_str());
  l.SetTextSize(0.05);
  
  return can;

}

void FitPlotter::DrawCatTree(const CategoryTree& CT, TCanvas* can){
  if(!can)
    return;

  double hlo = can->GetLeftMargin();
  double hhi = 1. - can->GetRightMargin();
  double hbo = can->GetBottomMargin();
  double hto = 1. - can->GetTopMargin();

  double frac = 0.15;
  
  can->cd();
  //const CategoryTree* GetParent()
  //GetBareLabel()
  int Depth = CT.GetDepth();
  int Nvis = CT.GetNVisible();
  int ibin = 0;
  int irun = 0;
  int brun = 0;

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.03);
  l.SetTextAlign(33);
  l.SetTextAngle(90);

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(7024);
  line->SetLineStyle(7);
  
  for(int i = 0; i < Depth; i++){
    int icolor = 7004 + 10*((i)%8);
    //int icolor = 7024;
    l.SetTextColor(icolor);
    line->SetLineColor(icolor);
    vector<const CategoryTree*> cats;
    CT.GetListDepth(cats, i);

    const CategoryTree* CurPar = cats[0]->GetVisibleParent();
    
    int Nc = cats.size();
    for(int c = 0; c < Nc; c++){
      irun += 1;
      ibin += 1;
      if(c+1 >= Nc ||
	 cats[c+1]->GetVisibleParent() != CurPar){
	if(irun > 1 && CurPar){
	  line->DrawLineNDC(hlo + (hhi-hlo)*double(brun)/double(Nvis), hto - double(i)/double(Depth)*frac*(hto-hbo),
	  		    hlo + (hhi-hlo)*double(brun)/double(Nvis), hbo);
	  l.DrawLatex(hlo + (hhi-hlo)*double(brun)/double(Nvis) + 0.0055,
		      hto - double(i)/double(Depth)*frac*(hto-hbo),
		      CurPar->GetSpectroscopicLabel().c_str());
	}
	irun = 0;
	brun = ibin;
	if(c+1 < Nc)
	  CurPar = cats[c+1]->GetVisibleParent();
      }
    }
  }
}
//////////////////////////////////////////////////
//modify
TCanvas* FitPlotter::Plot1Dstack(const string& can_name,
				 const VS& proc,
				 const CategoryTree& CT,
				 PlotType pType,
				 bool do_ratio){
  RestFrames::SetStyle();
  
  bool b_ratio = false;

  int Nproc = proc.size();
  if(Nproc == 0)
    return nullptr;
  
  vector<const CategoryTree*> CatTrees;
  CT.GetListDeepest(CatTrees);
  int Nvis = CatTrees.size();
  int Nbin;
  vector<int> Nbins;
  if(Nvis < 1)
    return nullptr;
  CategoryList CatList = GetCategories();
  CategoryList dumcat = CatList.Filter(*CatTrees[0]);
  const FitBin& fitbin = dumcat[0].GetFitBin();
  vector<CategoryList> dumcats;  
  if(pType == kFull)
    Nbin = fitbin.NBins();
  if(pType == kRISR || pType == kInv)
    Nbin = fitbin.NRBins();
  for(int v = 0; v < Nvis; v++){
    dumcats.push_back(CatList.Filter(*CatTrees[v]));
    CategoryList dumcat_tmp = CatList.Filter(*CatTrees[v]); 
    if(pType == kFull) Nbins.push_back(dumcats[v][0].GetFitBin().NBins());
    if(pType == kRISR || pType == kInv) Nbins.push_back(dumcats[v][0].GetFitBin().NRBins());
  } 
  
  if(CatList.GetN() < 1)
    return nullptr;
  
  // Processes
  VS                  labels;
  vector<int>         colors;
  vector<TH1D*>       hists[Nvis];
  vector<double>      total;
  
  VS                  labels_sig;
  vector<TH1D*>       hists_sig[Nvis];
  
  TH1D* hist_data[Nvis];
  double total_data = 0.;

  TH1D* hist_totbkg[Nvis];
  double total_totbkg = 0.;
  
  //int dumCatIdx = -999;
  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    double itot = 0.;
    
    ProcessList procs;
    
    ProcessType type = kBkg;
    TH1D* hist[Nvis];
    vector<double> hist_val;  

    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;
      
      Process pp = GetProcesses()[index];

      procs += pp;

      if(pp.Type() == kSig){
	type = kSig;
      }
      
      if(pp.Type() == kData){
	type = kData;
	if(do_ratio)
	  b_ratio = true;
      }
    }
    
    for(int v = 0; v < Nvis; v++){
      CategoryList cat = CatList.Filter(*CatTrees[v]);
      TH1D* h = nullptr;

      if(pType == kFull)
	h = GetAddedHist(Form("plothist_%d_%d_%s", i, v, can_name.c_str()), cat, procs);
      if(pType == kRISR || pType == kInv)
	h = IntegrateMperp(Form("plotintegratedhist_%d_%d_%s", i, v, can_name.c_str()),cat[0].GetFitBin(),
			   GetAddedHist(Form("plothist_%d_%d_%s", i, v, can_name.c_str()), cat, procs));

      if(h) itot += h->Integral();
	/*
	if(type == kBkg && proc[i] == "QCD") {
	  //cout << proc[i] << endl;
	  //cout << v << " " << cat[0].GetLabel() << " " << h->GetNbinsX() << endl;
	  
	  for(int bin = 1; bin < h->GetNbinsX()+1; bin++) {
	    cout << h->GetBinContent(bin) << endl;
	  }
	  cout << "integral: " << itot << endl;
	  
	  }*/
	//h->Draw(); 
      
      hist[v] = h;
    }

    if(type == kData && itot <= 1e-4) 
      b_ratio = false;
    if(itot <= 1e-4){// && type != kSig){
      //cout << itot << endl;
      continue;
    }
    
    if(type == kData){
      for(int v = 0; v < Nvis; v++)
	hist_data[v] = hist[v];
      total_data = itot;
    }

    if(type == kSig){
 
      labels_sig.push_back(GetSignalTitle(proc[i]));
      //cout << "here" << endl;
      //cout << GetSignalTitle(proc[i]) << endl;
      for(int v = 0; v < Nvis; v++){
	//CategoryList cat = CatList.Filter(*CatTrees[v]);
	//cout << "here" << v << endl;
	hists_sig[v].push_back(hist[v]);
	//cout << "here" << endl;
	/*
	for(int bin = 1; bin < hists_sig[v][i]->GetNbinsX()+1; bin++){
	  //cout << v << " " << hists_sig[v][i]->GetNbinsX() << endl;
	  //cout << hists_sig[v][i]->GetBinContent(bin) << endl;
	}
	*/
      }
    }
    //cout << "here" << endl;
    if(type == kBkg){
      //cout << proc[i] << endl;
      if(m_Title.count(proc[i]) != 0)
	labels += m_Title[proc[i]];
      else
	labels += proc[i];
      
      if(m_Color.count(proc[i]) != 0)
	colors.push_back(m_Color[proc[i]]);
      else
	colors.push_back(m_ColorDefault[i]);
      
      for(int v = 0; v < Nvis; v++){
	hists[v].push_back(hist[v]);
      }
      total.push_back(itot);
    }
  }

  if(m_FilePtr){
    Process totbkg("total_background", kBkg);
    ProcessList totbkgs;
    totbkgs += totbkg;
    for(int v = 0; v < Nvis; v++){
      CategoryList cat = CatList.Filter(*CatTrees[v]);
      TH1D* h = nullptr;

      if(pType == kFull)
        h = GetAddedHist(Form("plothist_%d_%s",  v, can_name.c_str()), cat, totbkgs);
      if(pType == kRISR || pType == kInv)
        h = IntegrateMperp(Form("plotintegratedhist_%d_%s", v, can_name.c_str()),cat[0].GetFitBin(),
			   GetAddedHist(Form("plothist_tot_%d_%s", v, can_name.c_str()), cat, totbkgs));

      hist_totbkg[v] = h;
      if(h) total_totbkg += h->Integral();
    }
  }

  if(total.size() < 1) return nullptr; 
  int Nsig = hists_sig[0].size();

  // sort the histograms by integral (N^2/2 brute force)
  int Nbkg = total.size();
  VS             vlabels;
  vector<int>    vcolors;
  vector<TH1D*>  vhists[Nvis];
  vector<double> vtotal;
  string stemp;
  int    itemp;
  TH1D*  htemp[Nvis];
  double ttemp;

  for(int i = 0; i < Nbkg; i++){
    vlabels.push_back(labels[i]);
    vcolors.push_back(colors[i]);
    for(int v = 0; v < Nvis; v++){
      vhists[v].push_back(hists[v][i]);
    }
    vtotal.push_back(total[i]);
    for(int j = vtotal.size()-2; j >= 0; j--){
      if(vtotal[j] < vtotal[j+1]){
	stemp = vlabels[j+1];
	itemp = vcolors[j+1];
	for(int v = 0; v < Nvis; v++)
	  htemp[v] = vhists[v][j+1];
	ttemp = vtotal[j+1];
	vlabels[j+1] = vlabels[j];
	vcolors[j+1] = vcolors[j];
	for(int v = 0; v < Nvis; v++)
	  vhists[v][j+1]  = vhists[v][j];
	vtotal[j+1]  = vtotal[j];
	vlabels[j] = stemp;
	vcolors[j] = itemp;
	for(int v = 0; v < Nvis; v++)
	  vhists[v][j]  = htemp[v];
	vtotal[j] = ttemp;
      } else {
	break;
      }
    }
  }

  vector<TH1D*> fhists;
  vector<TH1D*> fhists_sig;
  TH1D*         fhist_data   = nullptr;
  TH1D*         fhist_totbkg = nullptr;
  int		Nbins_total = std::accumulate(Nbins.begin(),Nbins.end(),0);
  int		Nbins_max = *std::max_element(Nbins.begin(),Nbins.end());
  /*
  for(int i = 0; i < Nsig; i++){
    fhists_sig.push_back(new TH1D(Form("fhistsig_%d_%s", i, can_name.c_str()),
				  Form("fhistsig_%d_%s", i, can_name.c_str()),
				  Nbin*Nvis, 0., Nbin*Nvis));
    for(int b = 0; b < Nbin; b++){
      for(int v = 0; v < Nvis; v++){
	if(hists_sig[v][i]){
	  int index = -1;
	  if(pType == kInv)
	    index = v*Nbin+b+1;
	  else 
	    index = b*Nvis+v+1;
	  cout << labels_sig[i] << " content: " << hists_sig[v][i]->GetBinContent(b+1) << endl;
	  cout << index << endl;
	  fhists_sig[i]->SetBinContent(index, hists_sig[v][i]->GetBinContent(b+1));
	  fhists_sig[i]->SetBinError(index, hists_sig[v][i]->GetBinError(b+1));
	}
      }
    }
  }
  */
  for(int i = 0; i < Nsig; i++){
    fhists_sig.push_back(new TH1D(Form("fhistsig_%d_%s", i, can_name.c_str()),
                                  Form("fhistsig_%d_%s", i, can_name.c_str()),
                                  Nbins_total, 0., Nbins_total));
    //cout << labels_sig[i] << endl;
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
        if(hists_sig[v][i]){
          int index = -1;
	  if(b+1 > Nbins[v]){
            offset += Nbins[v];
            continue;
          };

	  if(pType == kInv)
            index = offset+b+1;
          else
            index = b*Nvis+v+1;

          //cout << labels_sig[i] << " content: " << hists_sig[v][i]->GetBinContent(b+1) << endl;
          //cout << index << endl;
	  /*
	  cout << " total signal;  cat " << v
               << " has nbins " << Nbins[v]
               << ", Nbins_max: " << Nbins_max
               << ", b " << b
               << ", index " << index
               << " with " << hists_sig[v][i]->GetBinContent(b+1) << endl;
	  */
          fhists_sig[i]->SetBinContent(index, hists_sig[v][i]->GetBinContent(b+1));
          fhists_sig[i]->SetBinError(index, hists_sig[v][i]->GetBinError(b+1));
        }
	offset += Nbins[v];
      }
    }
  }

  if(total_data > 0.){
    fhist_data = new TH1D(Form("fhistdata_%s", can_name.c_str()),
			  Form("fhistdata_%s", can_name.c_str()),
			  Nbins_total, 0., Nbins_total);
    //int offset = 0; 
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	
	if(hist_data[v]){
	  int index;
	  if(b+1 > Nbins[v]){ 
	    offset += Nbins[v];
	    continue;
	  };
	  if(pType == kInv)
            //index = v*Nbins[v]+b+1-offset;
	    index = offset+b+1;
          else
            index = b*Nvis+v+1;
	  /*
	  cout << " total data;  cat " << v
               << " has nbins " << Nbins[v]
	       << " Nbins_max: " << Nbins_max
               << " b " << b
               << " index " << index
	       << " with " << hist_data[v]->GetBinContent(b+1) << endl;
	  */
	  fhist_data->SetBinContent(index, hist_data[v]->GetBinContent(b+1));
	  fhist_data->SetBinError(index, hist_data[v]->GetBinError(b+1));
	  offset += Nbins[v];
	  //if(b+2 > Nbins[v] && b+2 <= Nbins_max) offset++; 
	}
      }
    }
  }

  if(total_totbkg > 0.){
    fhist_totbkg = new TH1D(Form("fhisttotbkg_%s", can_name.c_str()),
			  Form("fhisttotbkg_%s", can_name.c_str()),
			  Nbins_total, 0., Nbins_total);
    //int offset = 0; 
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	if(hist_totbkg[v]){
	  int index;
	  if(b+1 > Nbins[v]){
            offset += Nbins[v];
            continue;
          };
	 if(pType == kInv)
	   //index = v*Nbins[v]+b+1-offset;
	   index = offset+b+1;
	  else
	    index = b*Nvis+v+1;

	  fhist_totbkg->SetBinContent(index, hist_totbkg[v]->GetBinContent(b+1));
	  fhist_totbkg->SetBinError(index, hist_totbkg[v]->GetBinError(b+1));
	  offset += Nbins[v];
	  //if(b+2 > Nbins[v] && b+2 <= Nbins_max) offset++; 
	  //cout <<  "total bkg;  cat " << v << " has nbins " << Nbins[v] << " b " << b << " index " << index << " with " << hist_totbkg[v]->GetBinContent(b+1) << " offset: " << offset << endl;	  
	}
      }
    }
  }
  for(int i = 0; i < Nbkg; i++){
    for(int v = 0; v < Nvis; v++){
      if(vhists[v][i]){
	for(int j = i+1; j < Nbkg; j++)
	  if(vhists[v][j])
	    vhists[v][i]->Add(vhists[v][j]);
      }
    }
  }
 for(int i = 0; i < Nbkg; i++){
    fhists.push_back(new TH1D(Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Nbins_total, 0., Nbins_total));
    //int offset = 0; 
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	int j = i;
	TH1D* hptr = nullptr;
	while(j < Nbkg && hptr == nullptr){
	hptr = vhists[v][j];
	  //	if(hptr != nullptr) cout << "Nvis: " << v << " Nbkg: " << j << " true Nbins: " << hptr->GetNbinsX() << " Nbin: " << Nbins[v] << endl;
	  j++;
	}
	//if(hptr == nullptr) 
	//cout << "nullptr" << endl;
	//continue;	
	if(b+1 > Nbins[v]){
	  offset += Nbins[v];
	  continue;
	};
	double content = 0;
	if(hptr) content = hptr->GetBinContent(b+1);

	int index;
	if(pType == kInv)
	  index =  offset+b+1;//v*Nbins[v]+b+1-offset;
	else
	  index = b*Nvis+v+1;
	fhists[i]->SetBinContent(index, content);//hptr->GetBinContent(b+1));
	fhists[i]->SetBinError(index, content);//hptr->GetBinError(b+1));
	offset += Nbins[v];
	//cout << labels_bkg[i] << " content: " << hists_sig[v][i]->GetBinContent(b+1) << endl; 
	//if(b+2 > Nbins[v] && b+2 <= Nbins_max) offset++; 
	//cout <<  proc[i] << " cat " << v << " has nbins " << Nbins[v] << " b " << b << " index " << index << " with " << content << " offset: " << offset << endl;	  
	//}
      }
    }
 }
 labels = vlabels;
 colors = vcolors;
 for(int b = 0; b < Nbins_total; b++)
   fhists[0]->GetXaxis()->SetBinLabel(b+1, "");
 
 fhists[0]->LabelsOption("v","X");
 
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11111111);
 TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
			    Form("can_%s", can_name.c_str()),
			    1200, 700);
 
 double yline = can->GetBottomMargin() - 0.02;
 double hlo = 0.105;
 double hhi = 0.2;
 double hbo = 0.25;//0.19;
 double hto = 0.07;
 
 double ratio_h = 0.18;
 
 can->SetLeftMargin(hlo);
 can->SetRightMargin(hhi);
 can->SetBottomMargin(hbo);
 can->SetTopMargin(hto);
 can->SetGridy();
 can->SetLogy();
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
 pad->SetLogy();
 pad->Draw();
 pad->cd();
 
 double hmax = fhists[0]->GetMaximum();
 double hmin = std::min(0.5, fhists[0]->GetMinimum());
 
 fhists[0]->Draw("hist");
 fhists[0]->GetXaxis()->SetNdivisions(30, 0, 0);
 fhists[0]->GetXaxis()->CenterTitle();
 fhists[0]->GetXaxis()->SetTitleFont(42);
 fhists[0]->GetXaxis()->SetTitleSize(0.045);
 fhists[0]->GetXaxis()->SetTitleOffset(1.0);
 fhists[0]->GetXaxis()->SetLabelFont(42);
 fhists[0]->GetXaxis()->SetLabelSize(0.05);
 fhists[0]->GetXaxis()->SetTitle("");
 fhists[0]->GetXaxis()->SetTickSize(0.);
 fhists[0]->GetYaxis()->CenterTitle();
 fhists[0]->GetYaxis()->SetTitleFont(42);
 fhists[0]->GetYaxis()->SetTitleSize(0.04);
 fhists[0]->GetYaxis()->SetTitleOffset(0.85);
 fhists[0]->GetYaxis()->SetLabelFont(42);
 fhists[0]->GetYaxis()->SetLabelSize(0.035);
 fhists[0]->GetYaxis()->SetTitle("number of events");
 
 for(int i = 0; i < Nbkg; i++){
   fhists[i]->SetLineColor(kBlack);
   fhists[i]->SetLineWidth(1.0);
   fhists[i]->SetFillColor(colors[i]);
   fhists[i]->SetFillStyle(1001);
   fhists[i]->Draw("SAME HIST");
 }
 
 TH1D* htot = nullptr;
 if(fhist_totbkg)
   htot = fhist_totbkg;
 else
   htot = fhists[0];
 
 vector<double> X;
 vector<double> Xerr;
 vector<double> Y;
 vector<double> Yerr;
 vector<double> Y_bkg_ratio;
 vector<double> Yerr_bkg_ratio;

 for(int i = 0; i < Nbins_total; i++){
   X.push_back(htot->GetXaxis()->GetBinCenter(i+1));
   Xerr.push_back(0.5);
   Y.push_back(htot->GetBinContent(i+1));
   Yerr.push_back(htot->GetBinError(i+1));
   Y_bkg_ratio.push_back(1.);
   Yerr_bkg_ratio.push_back(htot->GetBinError(i+1)/htot->GetBinContent(i+1));
 }

 TGraphErrors* gr = (TGraphErrors*) new TGraphErrors(Nbins_total, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);  
 gr->SetMarkerSize(0);
 gr->SetLineColor(kBlack);
 gr->SetFillColor(kBlack);
 gr->SetFillStyle(3244);
 gr->Draw("same p2");
 
 TGraphErrors* gr_bkg_ratio = (TGraphErrors*) new TGraphErrors(Nbins_total, &X[0],
							       &Y_bkg_ratio[0],  &Xerr[0], &Yerr_bkg_ratio[0]);  
 gr_bkg_ratio->SetMarkerSize(0);
 gr_bkg_ratio->SetLineColor(kBlack);
 gr_bkg_ratio->SetFillColor(kBlack);
 gr_bkg_ratio->SetFillStyle(3244);

 for(int i = 0; i < Nsig; i++){
   fhists_sig[i]->SetLineColor(m_SignalColor[i]);
   fhists_sig[i]->SetMarkerColor(m_SignalColor[i]);
   fhists_sig[i]->SetLineWidth(5);
   fhists_sig[i]->SetMarkerStyle(22);
   fhists_sig[i]->SetFillColor(kWhite);
   fhists_sig[i]->Draw("SAME hist");
   if(fhists_sig[i]->GetMaximum() > hmax)
     hmax = fhists_sig[i]->GetMaximum();
 }

 TH1D* fhist_data_ratio = nullptr;
 
 if(fhist_data){
   fhist_data->SetLineColor(kBlack);
   fhist_data->SetFillColor(kWhite);
   fhist_data->SetMarkerStyle(8);
   fhist_data->SetMarkerSize(1.);
   fhist_data->SetLineWidth(2);
   fhist_data->Draw("SAME ep");
   if(fhist_data->GetMaximum() > hmax)
     hmax = fhist_data->GetMaximum();
   
   fhist_data_ratio = (TH1D*) fhist_data->Clone(Form("data_ratio_%s", can_name.c_str()));
   for(int i = 0; i < Nbins_total; i++){
     fhist_data_ratio->SetBinContent(i+1, fhist_data_ratio->GetBinContent(i+1)/htot->GetBinContent(i+1));
     if(fhist_data_ratio->GetBinContent(i+1) > 0.)
       fhist_data_ratio->SetBinError(i+1, fhist_data_ratio->GetBinError(i+1)/htot->GetBinContent(i+1));
     else
       fhist_data_ratio->SetBinError(i+1, 0.);
   }
 }

 fhists[0]->GetYaxis()->SetRangeUser(hmin, hmax*10.);
 
 double eps = 0.0015;
 
 TLatex l;
 l.SetTextFont(42);
 l.SetNDC();
 
 l.SetTextSize(std::min(0.03, 1.5/double(Nbins_total)));
 l.SetTextFont(42);
 l.SetTextAlign(32);
 l.SetTextAngle(90);
 
 int Depth = CT.GetDepth();
 VS rlabels;
 VS vLabels;
 
 for(int v = 0; v < Nvis; v++){
   CategoryList tmpcat = CatList.Filter(*CatTrees[v]);
   for(int r = 0; r < Nbins[v]; r++){
     rlabels += tmpcat[0].GetFitBin()[r].GetRBinLabel();
   }
 }

 for(int v = 0; v < Nvis; v++)
   vLabels += CatTrees[v]->GetPlainLabel(Depth);
 
 for(int b = 0; b < Nbins_total; b++){
   string label;
   double xpos, ypos;
   if(pType == kInv){
      label = rlabels[b];
      xpos = hlo+(1.-hhi-hlo)/double(Nbins_total)*(0.5+b);
      ypos = yline+0.085;//0.0015;
   }
   else{
      label = CatTrees[b%Nvis]->GetPlainLabel(Depth);
      xpos = hlo+(1.-hhi-hlo)/double(Nbins_total)*(0.5+b);
      ypos = 1. - hto - 4*eps;
      l.SetTextColor(7004 + 10*((b%Nvis)%8));
    }

    l.DrawLatex(xpos,ypos,label.c_str());
  }
  /*
  for(int b = 0; b < Nvis*Nbin; b++){
    string label;
    if(pType == kInv)
      label = rlabels[b%Nbin];
    else 
      label = vLabels[b%Nvis];//CatTrees[b%Nvis]->GetPlainLabel(Depth);
    // if(colors[Nbkg-1]%1000 >=3)
    //   l.SetTextColor(7000 + 10*((b%Nvis)%8));
    // else
    l.SetTextColor(7004 + 10*((b%Nvis)%8));
    // l.DrawLatex(hlo+(1.-hhi-hlo)/double(Nvis*Nbin)*(0.5+b), hbo + 4*eps,
    // 		CatTrees[b%Nvis]->GetPlainLabel(Depth).c_str());
    l.DrawLatex(hlo+(1.-hhi-hlo)/double(Nvis*Nbin)*(0.5+b), 1.-hto - 4*eps,
		label.c_str());
  }
  */

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nbkg+Nsig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  can->Update();

  if(fhist_data)
    leg->AddEntry(fhist_data, "data");
  leg->AddEntry(gr, "total uncertainty","F");
  for(int i = 0; i < Nbkg; i++)
    leg->AddEntry(fhists[i], labels[i].c_str(), "F");
  for(int i = 0; i < Nsig; i++)
    leg->AddEntry(fhists_sig[i], labels_sig[i].c_str(), "L");
  leg->Draw("SAME");

  TPad* pad_ratio = nullptr;
  if(b_ratio){
    can->cd();
    pad_ratio = new TPad(Form("pad_%s", can_name.c_str()),
			 Form("pad_%s", can_name.c_str()),
			 0.0, hbo, 1.-hhi, hbo+ratio_h);

    pad_ratio->SetBottomMargin(0.);
    pad_ratio->SetLeftMargin(hlo*1.25);
    pad_ratio->SetRightMargin(0.);
    pad_ratio->SetTopMargin(0.);
    pad_ratio->SetGridy();
    pad_ratio->SetGridx();
    pad_ratio->Draw();
    pad_ratio->cd();


    fhist_data_ratio->Draw("ep");
    fhist_data_ratio->GetXaxis()->SetNdivisions(30, 0, 0);
    fhist_data_ratio->GetYaxis()->SetNdivisions(6, 5, 0);
    fhist_data_ratio->GetYaxis()->CenterTitle();
    fhist_data_ratio->GetYaxis()->SetTitleFont(42);
    fhist_data_ratio->GetYaxis()->SetTitleSize(0.21);
    fhist_data_ratio->GetYaxis()->SetTitleOffset(0.22);
    fhist_data_ratio->GetYaxis()->SetLabelFont(42);
    fhist_data_ratio->GetYaxis()->SetLabelSize(0.2);
    fhist_data_ratio->GetYaxis()->SetTitle("#frac{data}{bkg model}");

    gr_bkg_ratio->Draw("same p2");
    fhist_data_ratio->Draw("same ep");
  }

  if(pType == kFull)
    DrawMR(fitbin, can, pad, pad_ratio);
  if(pType == kRISR)
    DrawRM(fitbin, can, pad, pad_ratio);
  if(pType == kInv)
    //DrawRM(fitbin,can,can,vLabels);
    DrawRM(dumcats,can,can,vLabels);

  l.SetTextAngle(0);
  l.SetTextColor(kBlack);
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, string("Regions "+CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02,  "#bf{CMS} 138 fb^{-1} (13 TeV)");// m_CMSLabel.c_str());

  //cout << 1.-hhi-eps*4 << " " <<  1.-hto+0.02 << " " <<  string("Regions "+CT.GetSpectroscopicLabel()).c_str() << endl;
  //cout << hlo+eps*4 << " " << 1.-hto+0.02 << " " << m_CMSLabel.c_str() << endl;

  can->Update();

  // make y-axis in ratio plot be symmetric and cover all data and mc points
  // fhist_data_ratio doesn't return the correct value for GetYaxis()->GetXmax/min() so we brute force it
  // start with getting the min and max deviation from the first bin (0 bin is underflow)
  double ratio_axis_min = fhist_data_ratio->GetBinContent(1)-fhist_data_ratio->GetBinError(1);
  double ratio_axis_max = fhist_data_ratio->GetBinContent(1)+fhist_data_ratio->GetBinError(1);

  // loop over remaining bins to find max deviation up and down for the data (b = nbins is last bin)
  for (int b = 2; b <= fhist_data_ratio->GetNbinsX(); b++) {
      if(fhist_data_ratio->GetBinContent(b)-fhist_data_ratio->GetBinError(b) < ratio_axis_min)
        ratio_axis_min = fhist_data_ratio->GetBinContent(b)-fhist_data_ratio->GetBinError(b);
      if(fhist_data_ratio->GetBinContent(b)+fhist_data_ratio->GetBinError(b) > ratio_axis_max)
        ratio_axis_max = fhist_data_ratio->GetBinContent(b)+fhist_data_ratio->GetBinError(b);
    }

  // check if the bkg ratio has larger deviation than the data
  if(ratio_axis_min > gr_bkg_ratio->GetYaxis()->GetXmin())
    ratio_axis_min = gr_bkg_ratio->GetYaxis()->GetXmin();
  if(ratio_axis_max < gr_bkg_ratio->GetYaxis()->GetXmax())
    ratio_axis_max = gr_bkg_ratio->GetYaxis()->GetXmax();

  // round to nearest three decimal places
  ratio_axis_min = std::floor(ratio_axis_min*1000.)/1000.;
  ratio_axis_max = std::ceil(ratio_axis_max*1000.)/1000.;
  
  // make it so we have the same range from 1 for min and max
  //   example: 0.9, 1.15 becomes 0.85, 1.15
  if(fabs(1.-ratio_axis_min) > fabs(ratio_axis_max-1.)){
    ratio_axis_max = 1.-ratio_axis_min+1.;
  }
  else if(fabs(ratio_axis_max-1.) > fabs(1.-ratio_axis_min)){
    ratio_axis_min = 1.+1.-ratio_axis_max;
  }

  // add offset so there is a bit of white space between histograms and axis
  if(ratio_axis_max > 1.5){
    ratio_axis_min -= 0.11;
    ratio_axis_max += 0.11;
  }
  else if(ratio_axis_max > 1.25){
    ratio_axis_min -= 0.06;
    ratio_axis_max += 0.06;
  }
  else{
    ratio_axis_min -= 0.02;
    ratio_axis_max += 0.02;
  }
  
  // minimum deviation from 1 is 7% so we can clearly see 5% on the axis label (5% is the default minimum deviation to show)
  if(ratio_axis_min > 0.93 && ratio_axis_max < 1.07){
    ratio_axis_min = 0.93;
    ratio_axis_max = 1.07;
  }
  fhist_data_ratio->GetYaxis()->SetRangeUser(ratio_axis_min, ratio_axis_max);
  
  can->Update();
  can->SaveAs("plots/"+TString(can_name)+".pdf");
  can->SaveAs("plots/"+TString(can_name)+".gif");

  return can;
}

TCanvas* FitPlotter::Plot1DstackPullRatio(const string& can_name,
					  const VS& proc,
					  const CategoryTree& CT,
					  PlotType pType,
					  bool do_ratio){
  RestFrames::SetStyle();
  
  bool b_ratio = false;

  int Nproc = proc.size();
  if(Nproc == 0)
    return nullptr;
  
  vector<const CategoryTree*> CatTrees;
  CT.GetListDeepest(CatTrees);
  int Nvis = CatTrees.size();
  int Nbin;
  vector<int> Nbins;
  if(Nvis < 1)
    return nullptr;
  CategoryList CatList = GetCategories();
  CategoryList dumcat = CatList.Filter(*CatTrees[0]);
  const FitBin& fitbin = dumcat[0].GetFitBin();
  vector<CategoryList> dumcats;  
  if(pType == kFull)
    Nbin = fitbin.NBins();
  if(pType == kRISR || pType == kInv)
    Nbin = fitbin.NRBins();
  for(int v = 0; v < Nvis; v++){
    dumcats.push_back(CatList.Filter(*CatTrees[v]));
    CategoryList dumcat_tmp = CatList.Filter(*CatTrees[v]); 
    if(pType == kFull) Nbins.push_back(dumcats[v][0].GetFitBin().NBins());
    if(pType == kRISR || pType == kInv) Nbins.push_back(dumcats[v][0].GetFitBin().NRBins());
  } 
  
  if(CatList.GetN() < 1)
    return nullptr;
  
  // Processes
  VS                  labels;
  vector<int>         colors;
  vector<TH1D*>       hists[Nvis];
  vector<double>      total;
  
  VS                  labels_sig;
  vector<TH1D*>       hists_sig[Nvis];
  
  TH1D* hist_data[Nvis];
  double total_data = 0.;

  TH1D* hist_totbkg[Nvis];
  double total_totbkg = 0.;
  
  //int dumCatIdx = -999;
  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc[i]) != 0)
      vproc = m_Strings[proc[i]];
    else
      vproc += proc[i];

    double itot = 0.;
    
    ProcessList procs;
    
    ProcessType type = kBkg;
    TH1D* hist[Nvis];
    vector<double> hist_val;  

    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
      continue;
      
      Process pp = GetProcesses()[index];

      procs += pp;

      if(pp.Type() == kSig){
      type = kSig;
      }
      
      if(pp.Type() == kData){
      type = kData;
      if(do_ratio)
        b_ratio = true;
      }
    }
    
    for(int v = 0; v < Nvis; v++){
      CategoryList cat = CatList.Filter(*CatTrees[v]);
      TH1D* h = nullptr;

      if(pType == kFull)
      h = GetAddedHist(Form("plothist_%d_%d_%s", i, v, can_name.c_str()), cat, procs);
      if(pType == kRISR || pType == kInv)
      h = IntegrateMperp(Form("plotintegratedhist_%d_%d_%s", i, v, can_name.c_str()),cat[0].GetFitBin(),
      	     GetAddedHist(Form("plothist_%d_%d_%s", i, v, can_name.c_str()), cat, procs));

      if(h) itot += h->Integral();      
      hist[v] = h;
    }

    if(type == kData && itot <= 1e-4) 
      b_ratio = false;
    if(itot <= 1e-4)
      continue;
    
    if(type == kData){
      for(int v = 0; v < Nvis; v++)
      hist_data[v] = hist[v];
      total_data = itot;
    }

    if(type == kSig) {
      labels_sig.push_back(GetSignalTitle(proc[i]));

      for(int v = 0; v < Nvis; v++)
	hists_sig[v].push_back(hist[v]);
      
    }

    if(type == kBkg){
      if(m_Title.count(proc[i]) != 0)
	labels += m_Title[proc[i]];
      else
	labels += proc[i];
      
      if(m_Color.count(proc[i]) != 0)
	colors.push_back(m_Color[proc[i]]);
      else
	colors.push_back(m_ColorDefault[i]);
      
      for(int v = 0; v < Nvis; v++){
	hists[v].push_back(hist[v]);
      }
      total.push_back(itot);
    }
  }

  if(m_FilePtr){
    Process totbkg("total_background", kBkg);
    ProcessList totbkgs;
    totbkgs += totbkg;
    for(int v = 0; v < Nvis; v++){
      CategoryList cat = CatList.Filter(*CatTrees[v]);
      TH1D* h = nullptr;

      if(pType == kFull)
        h = GetAddedHist(Form("plothist_%d_%s",  v, can_name.c_str()), cat, totbkgs);
      if(pType == kRISR || pType == kInv)
        h = IntegrateMperp(Form("plotintegratedhist_%d_%s", v, can_name.c_str()),cat[0].GetFitBin(),
	       GetAddedHist(Form("plothist_tot_%d_%s", v, can_name.c_str()), cat, totbkgs));

      hist_totbkg[v] = h;
      if(h) total_totbkg += h->Integral();
    }
  }

  if(total.size() < 1) return nullptr; 
  int Nsig = hists_sig[0].size();

  // sort the histograms by integral (N^2/2 brute force)
  int Nbkg = total.size();
  VS             vlabels;
  vector<int>    vcolors;
  vector<TH1D*>  vhists[Nvis];
  vector<double> vtotal;
  string stemp;
  int    itemp;
  TH1D*  htemp[Nvis];
  double ttemp;

  for(int i = 0; i < Nbkg; i++){
    vlabels.push_back(labels[i]);
    vcolors.push_back(colors[i]);
    for(int v = 0; v < Nvis; v++){
      vhists[v].push_back(hists[v][i]);
    }
    vtotal.push_back(total[i]);
    for(int j = vtotal.size()-2; j >= 0; j--){
      if(vtotal[j] < vtotal[j+1]){
      stemp = vlabels[j+1];
      itemp = vcolors[j+1];
      for(int v = 0; v < Nvis; v++)
        htemp[v] = vhists[v][j+1];
	ttemp = vtotal[j+1];
	vlabels[j+1] = vlabels[j];
	vcolors[j+1] = vcolors[j];
	for(int v = 0; v < Nvis; v++)
	  vhists[v][j+1]  = vhists[v][j];
	  vtotal[j+1]  = vtotal[j];
	  vlabels[j] = stemp;
	  vcolors[j] = itemp;
	  for(int v = 0; v < Nvis; v++)
	    vhists[v][j]  = htemp[v];
	    vtotal[j] = ttemp;
      } else {
      break;
      }
    }
  }

  vector<TH1D*> fhists;
  vector<TH1D*> fhists_sig;
  TH1D*         fhist_data   = nullptr;
  TH1D*         fhist_totbkg = nullptr;
  int			     Nbins_total = std::accumulate(Nbins.begin(),Nbins.end(),0);
  int			     		 Nbins_max = *std::max_element(Nbins.begin(),Nbins.end());

  for(int i = 0; i < Nsig; i++){
    fhists_sig.push_back(new TH1D(Form("fhistsig_%d_%s", i, can_name.c_str()),
                                  Form("fhistsig_%d_%s", i, can_name.c_str()),
                                  Nbins_total, 0., Nbins_total));

    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
        if(hists_sig[v][i]){
          int index = -1;
	  if(b+1 > Nbins[v]){
	    offset += Nbins[v];
            continue;
	  }
	  
	  if(pType == kInv)
            index = offset+b+1;
          else
            index = b*Nvis+v+1;
	  
          fhists_sig[i]->SetBinContent(index, hists_sig[v][i]->GetBinContent(b+1));
          fhists_sig[i]->SetBinError(index, hists_sig[v][i]->GetBinError(b+1));
        }
	offset += Nbins[v];
      }
    }
  }
  
  if(total_data > 0.){
    fhist_data = new TH1D(Form("fhistdata_%s", can_name.c_str()),
			  Form("fhistdata_%s", can_name.c_str()),
			  Nbins_total, 0., Nbins_total);
    
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	
	if(hist_data[v]){
	  int index;
	  if(b+1 > Nbins[v]){ 
	    offset += Nbins[v];
	  continue;
	  };
	  if(pType == kInv)
	    index = offset+b+1;
	  else
	    index = b*Nvis+v+1;
	  fhist_data->SetBinContent(index, hist_data[v]->GetBinContent(b+1));
	  fhist_data->SetBinError(index, hist_data[v]->GetBinError(b+1));
	  offset += Nbins[v];
	}
      }
    }
  }

  if(total_totbkg > 0.){
    fhist_totbkg = new TH1D(Form("fhisttotbkg_%s", can_name.c_str()),
			    Form("fhisttotbkg_%s", can_name.c_str()),
			    Nbins_total, 0., Nbins_total);

    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	if(hist_totbkg[v]){
	  int index;
	  if(b+1 > Nbins[v]){
            offset += Nbins[v];
            continue;
          }
	  if(pType == kInv)
	    index = offset+b+1;
	  else
	    index = b*Nvis+v+1;
	  
	  fhist_totbkg->SetBinContent(index, hist_totbkg[v]->GetBinContent(b+1));
	  fhist_totbkg->SetBinError(index, hist_totbkg[v]->GetBinError(b+1));
	  offset += Nbins[v];
	}
      }
    }
  }

  for(int i = 0; i < Nbkg; i++){
    for(int v = 0; v < Nvis; v++){
      if(vhists[v][i]){
      for(int j = i+1; j < Nbkg; j++)
        if(vhists[v][j])
	    vhists[v][i]->Add(vhists[v][j]);
      }
    }
  }
  for(int i = 0; i < Nbkg; i++){
    fhists.push_back(new TH1D(Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Form("fhistsbkg_%d_%s", i, can_name.c_str()),
			      Nbins_total, 0., Nbins_total));
    
    for(int b = 0; b < Nbins_max; b++){
      int offset = 0;
      for(int v = 0; v < Nvis; v++){
	int j = i;
	TH1D* hptr = nullptr;
	while(j < Nbkg && hptr == nullptr){
	  hptr = vhists[v][j];
	  j++;
	}
	if(b+1 > Nbins[v]){
	  offset += Nbins[v];
	  continue;
	}
	double content = 0;
	if(hptr) content = hptr->GetBinContent(b+1);
	
	int index;
	if(pType == kInv)
	  index =  offset+b+1;//v*Nbins[v]+b+1-offset;
	else
	  index = b*Nvis+v+1;
	fhists[i]->SetBinContent(index, content);//hptr->GetBinContent(b+1));
	fhists[i]->SetBinError(index, content);//hptr->GetBinError(b+1));
	offset += Nbins[v];
      }
    }
 }

 labels = vlabels;
 colors = vcolors;
 for(int b = 0; b < Nbins_total; b++)
   fhists[0]->GetXaxis()->SetBinLabel(b+1, "");
 
 fhists[0]->LabelsOption("v","X");
 
 gStyle->SetOptTitle(0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11111111);
 TCanvas* can = new TCanvas(Form("can_pullRatio_%s", can_name.c_str()),
			    Form("can_pullRatio_%s", can_name.c_str()),
			    1200, 700);
 
 double yline = can->GetBottomMargin() - 0.02;
 double hlo = 0.105;
 double hhi = 0.2;
 double hbo = 0.25;//0.19;
 double hto = 0.07;
 
 double ratio_h = 0.18;
 
 can->SetLeftMargin(hlo);
 can->SetRightMargin(hhi);
 can->SetBottomMargin(hbo);
 can->SetTopMargin(hto);
 can->SetGridy();
 can->SetLogy();
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
 pad->SetLogy();
 pad->Draw();
 pad->cd();
 
 double hmax = fhists[0]->GetMaximum();
 double hmin = std::min(0.5, fhists[0]->GetMinimum());
 
 fhists[0]->Draw("hist");
 fhists[0]->GetXaxis()->SetNdivisions(30, 0, 0);
 fhists[0]->GetXaxis()->CenterTitle();
 fhists[0]->GetXaxis()->SetTitleFont(42);
 fhists[0]->GetXaxis()->SetTitleSize(0.045);
 fhists[0]->GetXaxis()->SetTitleOffset(1.0);
 fhists[0]->GetXaxis()->SetLabelFont(42);
 fhists[0]->GetXaxis()->SetLabelSize(0.05);
 fhists[0]->GetXaxis()->SetTitle("");
 fhists[0]->GetXaxis()->SetTickSize(0.);
 fhists[0]->GetYaxis()->CenterTitle();
 fhists[0]->GetYaxis()->SetTitleFont(42);
 fhists[0]->GetYaxis()->SetTitleSize(0.04);
 fhists[0]->GetYaxis()->SetTitleOffset(0.85);
 fhists[0]->GetYaxis()->SetLabelFont(42);
 fhists[0]->GetYaxis()->SetLabelSize(0.035);
 fhists[0]->GetYaxis()->SetTitle("number of events");
 
 for(int i = 0; i < Nbkg; i++){
   fhists[i]->SetLineColor(kBlack);
   fhists[i]->SetLineWidth(1.0);
   fhists[i]->SetFillColor(colors[i]);
   fhists[i]->SetFillStyle(1001);
   fhists[i]->Draw("SAME HIST");
 }
 
 TH1D* htot = nullptr;
 if(fhist_totbkg)
   htot = fhist_totbkg;
 else
   htot = fhists[0];
 
 vector<double> X;
 vector<double> Xerr;
 vector<double> Y;
 vector<double> Yerr;
 vector<double> Y_bkg_ratio;
 vector<double> Yerr_bkg_ratio;

 for(int i = 0; i < Nbins_total; i++){
   X.push_back(htot->GetXaxis()->GetBinCenter(i+1));
   Xerr.push_back(0.5);
   Y.push_back(htot->GetBinContent(i+1));
   Yerr.push_back(htot->GetBinError(i+1));
   Y_bkg_ratio.push_back(0.);
   Yerr_bkg_ratio.push_back(1.);
 }

 TGraphErrors* gr = (TGraphErrors*) new TGraphErrors(Nbins_total, &X[0], &Y[0],  &Xerr[0], &Yerr[0]);  
 gr->SetMarkerSize(0);
 gr->SetLineColor(kBlack);
 gr->SetFillColor(kBlack);
 gr->SetFillStyle(3244);
 gr->Draw("same p2");
 
 TGraphErrors* gr_bkg_ratio = (TGraphErrors*) new TGraphErrors(Nbins_total, &X[0],
 	       		      		      	  			           &Y_bkg_ratio[0],  &Xerr[0], &Yerr_bkg_ratio[0]);  
 gr_bkg_ratio->SetMarkerSize(0);
 gr_bkg_ratio->SetLineColor(kBlack);
 gr_bkg_ratio->SetFillColor(kBlack);
 gr_bkg_ratio->SetFillStyle(3244);

 for(int i = 0; i < Nsig; i++){
   fhists_sig[i]->SetLineColor(m_SignalColor[i]);
   fhists_sig[i]->SetMarkerColor(m_SignalColor[i]);
   fhists_sig[i]->SetLineWidth(5);
   fhists_sig[i]->SetMarkerStyle(22);
   fhists_sig[i]->SetFillColor(kWhite);
   fhists_sig[i]->Draw("SAME hist");
   if(fhists_sig[i]->GetMaximum() > hmax)
     hmax = fhists_sig[i]->GetMaximum();
 }
 
 TH1D* fhist_data_ratio = nullptr;
 
 if(fhist_data){
   fhist_data->SetLineColor(kBlack);
   fhist_data->SetFillColor(kWhite);
   fhist_data->SetMarkerStyle(8);
   fhist_data->SetMarkerSize(1.);
   fhist_data->SetLineWidth(2);
   fhist_data->Draw("SAME ep");
   if(fhist_data->GetMaximum() > hmax)
     hmax = fhist_data->GetMaximum();
   
   fhist_data_ratio = (TH1D*) fhist_data->Clone(Form("data_ratio_%s", can_name.c_str()));
   for(int i = 0; i < Nbins_total; i++){
     const double data_yield = fhist_data_ratio->GetBinContent(i+1);
     const double bkg_yield = htot->GetBinContent(i+1);
     const double bkg_error = htot->GetBinError(i+1);
     fhist_data_ratio->SetBinContent(i+1, (data_yield - bkg_yield) / bkg_error);//sqrt(bkg_yield + bkg_error*bkg_error)); 
     fhist_data_ratio->SetBinError(i+1, sqrt(data_yield / (bkg_error*bkg_error)) );//(bkg_yield + bkg_error*bkg_error)) );
   }
 }

 fhists[0]->GetYaxis()->SetRangeUser(hmin, hmax*10.);
 
 double eps = 0.0015;
 
 TLatex l;
 l.SetTextFont(42);
 l.SetNDC();
 
 l.SetTextSize(std::min(0.03, 1.5/double(Nbins_total)));
 l.SetTextFont(42);
 l.SetTextAlign(32);
 l.SetTextAngle(90);
 
 int Depth = CT.GetDepth();
 VS rlabels;
 VS vLabels;
 
 for(int v = 0; v < Nvis; v++){
   CategoryList tmpcat = CatList.Filter(*CatTrees[v]);
   for(int r = 0; r < Nbins[v]; r++){
     rlabels += tmpcat[0].GetFitBin()[r].GetRBinLabel();
   }
 }

 for(int v = 0; v < Nvis; v++)
   vLabels += CatTrees[v]->GetPlainLabel(Depth);
 
 for(int b = 0; b < Nbins_total; b++){
   string label;
   double xpos, ypos;
   if(pType == kInv){
      label = rlabels[b];
      xpos = hlo+(1.-hhi-hlo)/double(Nbins_total)*(0.5+b);
      ypos = yline+0.085;//0.0015;
   }
   else{
      label = CatTrees[b%Nvis]->GetPlainLabel(Depth);
      xpos = hlo+(1.-hhi-hlo)/double(Nbins_total)*(0.5+b);
      ypos = 1. - hto - 4*eps;
      l.SetTextColor(7004 + 10*((b%Nvis)%8));
   }

    l.DrawLatex(xpos,ypos,label.c_str());
  }

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nbkg+Nsig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  can->Update();

  if(fhist_data)
    leg->AddEntry(fhist_data, "data");
  leg->AddEntry(gr, "total uncertainty","F");
  for(int i = 0; i < Nbkg; i++)
    leg->AddEntry(fhists[i], labels[i].c_str(), "F");
  for(int i = 0; i < Nsig; i++)
    leg->AddEntry(fhists_sig[i], labels_sig[i].c_str(), "L");
  leg->Draw("SAME");

  TPad* pad_ratio = nullptr;
  if(b_ratio){
    can->cd();
    pad_ratio = new TPad(Form("pad_%s", can_name.c_str()),
			 Form("pad_%s", can_name.c_str()),
			 0.0, hbo, 1.-hhi, hbo+ratio_h);

    pad_ratio->SetBottomMargin(0.);
    pad_ratio->SetLeftMargin(hlo*1.25);
    pad_ratio->SetRightMargin(0.);
    pad_ratio->SetTopMargin(0.);
    pad_ratio->SetGridy();
    pad_ratio->SetGridx();
    pad_ratio->Draw();
    pad_ratio->cd();

    fhist_data_ratio->Draw("ep");
    fhist_data_ratio->GetXaxis()->SetNdivisions(30, 0, 0);
    fhist_data_ratio->GetYaxis()->SetNdivisions(6, 5, 0);
    fhist_data_ratio->GetYaxis()->CenterTitle();
    fhist_data_ratio->GetYaxis()->SetTitleFont(42);
    fhist_data_ratio->GetYaxis()->SetTitleSize(0.21);
    fhist_data_ratio->GetYaxis()->SetTitleOffset(0.22);
    fhist_data_ratio->GetYaxis()->SetLabelFont(42);
    fhist_data_ratio->GetYaxis()->SetLabelSize(0.2);
    fhist_data_ratio->GetYaxis()->SetTitle("#frac{data-model}{#sigma_{post-fit}}");

    gr_bkg_ratio->Draw("same p2");
    fhist_data_ratio->Draw("same ep");
  }

  if(pType == kFull)
    DrawMR(fitbin, can, pad, pad_ratio);
  if(pType == kRISR)
    DrawRM(fitbin, can, pad, pad_ratio);
  if(pType == kInv)
    DrawRM(dumcats,can,can,vLabels);

  l.SetTextAngle(0);
  l.SetTextColor(kBlack);
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, string("Regions "+CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, "#bf{CMS} #it{Preliminary} 138 fb^{-1} (13 TeV)");//m_CMSLabel.c_str());

  can->Update();

  // make y-axis in ratio plot be symmetric and cover all data and mc points
  // fhist_data_ratio doesn't return the correct value for GetYaxis()->GetXmax/min() so we brute force it
  // start with getting the min and max deviation from the first bin (0 bin is underflow)
  double ratio_axis_min = fhist_data_ratio->GetBinContent(1)-fhist_data_ratio->GetBinError(1);
  double ratio_axis_max = fhist_data_ratio->GetBinContent(1)+fhist_data_ratio->GetBinError(1);
  
  // loop over remaining bins to find max deviation up and down for the data (b = nbins is last bin)
  for (int b = 2; b <= fhist_data_ratio->GetNbinsX(); b++) {
      if(fhist_data_ratio->GetBinContent(b)-fhist_data_ratio->GetBinError(b) < ratio_axis_min)
        ratio_axis_min = fhist_data_ratio->GetBinContent(b)-fhist_data_ratio->GetBinError(b);
      if(fhist_data_ratio->GetBinContent(b)+fhist_data_ratio->GetBinError(b) > ratio_axis_max)
        ratio_axis_max = fhist_data_ratio->GetBinContent(b)+fhist_data_ratio->GetBinError(b);
    }

  // check if the bkg ratio has larger deviation than the data
  if(ratio_axis_min > gr_bkg_ratio->GetYaxis()->GetXmin())
    ratio_axis_min = gr_bkg_ratio->GetYaxis()->GetXmin();
  if(ratio_axis_max < gr_bkg_ratio->GetYaxis()->GetXmax())
    ratio_axis_max = gr_bkg_ratio->GetYaxis()->GetXmax();

  // round to nearest three decimal places
  ratio_axis_min = std::floor(ratio_axis_min*1000.)/1000.;
  ratio_axis_max = std::ceil(ratio_axis_max*1000.)/1000.;
  
  // make it so we have the same range from 1 for min and max
  //   example: 0.9, 1.15 becomes 0.85, 1.15
  /*
  if(fabs(1.-ratio_axis_min) > fabs(ratio_axis_max-1.)){
    ratio_axis_max = 1.-ratio_axis_min+1.;
  }
  else if(fabs(ratio_axis_max-1.) > fabs(1.-ratio_axis_min)){
    ratio_axis_min = 1.+1.-ratio_axis_max;
  }
  */
  // add offset so there is a bit of white space between histograms and axis
  if(ratio_axis_max > 1.5){
    ratio_axis_min -= 0.11;
    ratio_axis_max += 0.11;
  }
  else if(ratio_axis_max > 1.25){
    ratio_axis_min -= 0.06;
    ratio_axis_max += 0.06;
  }
  else{
    ratio_axis_min -= 0.02;
    ratio_axis_max += 0.02;
  }
  
  // minimum deviation from 1 is 7% so we can clearly see 5% on the axis label (5% is the default minimum deviation to show)
  //if(ratio_axis_min > 0.93 && ratio_axis_max < 1.07){
  ratio_axis_min = -5.4;
  ratio_axis_max = 5.4;
    //}
    
  fhist_data_ratio->GetYaxis()->SetRangeUser(ratio_axis_min, ratio_axis_max);
  
  can->Update();
  can->SaveAs("plots/"+TString(can_name)+"_pullRatio.pdf");
  can->SaveAs("plots/"+TString(can_name)+"_pullRatio.gif");

  return can;
}

void FitPlotter::DrawMR(const FitBin& fitbin, TCanvas* can, TPad* pad, TPad* pad_ratio){
  pad->cd();

  double eps = 0.0015;

  double hlo = can->GetLeftMargin();
  double hhi = can->GetRightMargin();
  double hbo = can->GetBottomMargin();
  double hto = can->GetTopMargin();

  double hlo_rat = 0.;
  double hhi_rat = 0.;
  if(pad_ratio){
    hlo_rat = pad_ratio->GetLeftMargin();
    hhi_rat = pad_ratio->GetRightMargin();
  }
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  int NR = fitbin.NRBins();
  int NB = fitbin.NBins();
  
  VS mlabels;
  VS rlabels;
  for(int r = 0; r < NR; r++)
    mlabels += fitbin[r].GetMBinLabels();

  int lmax = 0;
  for(int b = 0; b < NB; b++){
    int len = mlabels[b].length();
    if(mlabels[b].find("#infty") != std::string::npos)
      len -= 5;
    if(len > lmax)
      lmax = len;
  }
  
  string space = "";
  for(int l = 0; l < 1.6*lmax; l++)
    space += " ";

  for(int r = 0; r < NR; r++)
    rlabels += fitbin[r].GetRBinLabel();
  
  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  
  l.SetTextSize(0.033);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(2);
  double lo = hlo;
  double hi = hlo;
  double lo_ratio = hlo_rat;
  double hi_ratio = hlo_rat;
  double yline = hbo-0.018*lmax;
  int ib = 0;
  for(int r = 0; r < NR; r++){
    int NM = fitbin[r].NBins();
    lo = hi;
    lo_ratio = hi_ratio;
    hi = double(NM)/double(NB)*(1.-hhi-hlo) + lo;
    hi_ratio = double(NM)/double(NB)*(1.-hhi_rat-hlo_rat) + lo_ratio;
    
    line->SetLineStyle(1);
    line->DrawLineNDC(lo + eps, yline,
  		      lo + eps, yline + 6*eps);
    line->DrawLineNDC(hi - eps, yline,
  		      hi - eps, yline + 6*eps);
   
    line->DrawLineNDC(lo + eps, yline,
  		      hi - eps, yline);
  
    line->SetLineStyle(1);
    if(r < NR-1){
      line->DrawLineNDC(hi, yline + 10*eps , hi, 1.-hto);
      if(pad_ratio){
	pad_ratio->cd();
	line->DrawLineNDC(hi_ratio, 0., hi_ratio, 1.);
	pad->cd();
      }
    }
    l.DrawLatex((hi+lo)/2., yline - 8*eps, rlabels[r].c_str());
  }
  
  line->SetLineStyle(1);
  l.SetTextAngle(90);
  l.SetTextAlign(32);
  for(int b = 0; b < NB; b++){
    l.DrawLatex(hlo + (1.-hhi-hlo)*(0.5+b)/double(NB), hbo - 4*eps, mlabels[b].c_str());
    if(b > 0){
      line->DrawLineNDC(hlo + (1.-hhi-hlo)/double(NB)*b , hbo , hlo + (1.-hhi-hlo)/double(NB)*b, 1.-hto);
      if(pad_ratio){
	pad_ratio->cd();
	line->DrawLineNDC(hlo_rat + (1.-hhi_rat-hlo_rat)/double(NB)*b, 0.,
			  hlo_rat + (1.-hhi_rat-hlo_rat)/double(NB)*b, 1.);
	pad->cd();
      }
    }
  }
  
  l.SetTextAngle(0);
  l.SetTextAlign(32);
  l.SetTextSize(0.035);
  l.DrawLatex(hlo, (hbo+yline)/2.+eps, "M_{#perp}   [GeV] #in");
  l.DrawLatex(hlo, yline - 16*eps, "#scale[1.15]{R_{ISR}} #in");
}

void FitPlotter::DrawRM(const FitBin& fitbin, TCanvas* can, TPad* pad, TPad* pad_ratio){
  pad->cd();

  double eps = 0.0015;

  double hlo = can->GetLeftMargin();
  double hhi = can->GetRightMargin();
  double hbo = can->GetBottomMargin();
  double hto = can->GetTopMargin();

  double hlo_rat = 0.;
  double hhi_rat = 0.;

  if(pad_ratio){
    hlo_rat = pad_ratio->GetLeftMargin();
    hhi_rat = pad_ratio->GetRightMargin();
  }
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  int NR = fitbin.NBins();
  int NB = fitbin.NRBins();

  VS rlabels;
  for(int r = 0; r < NR; r++)
    rlabels += fitbin[r].GetRBinLabel();
  
  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  
  l.SetTextSize(0.033);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(2);

  double lo = hlo;
  double hi = hlo;
  double lo_ratio = hlo_rat;
  double hi_ratio = hlo_rat;
  double yline = hbo-0.02;

  for(int r = 0; r < NB; r++){
    lo = hi;
    lo_ratio = hi_ratio;
    hi = 1./double(NB)*(1.-hhi-hlo) + lo;
    hi_ratio = 1./double(NB)*(1.-hhi_rat-hlo_rat) + lo_ratio;
    /*
    line->SetLineStyle(1);
    line->DrawLineNDC(lo + eps, yline,
                      lo + eps, yline + 6*eps);
    line->DrawLineNDC(hi - eps, yline,
                      hi - eps, yline + 6*eps);

    line->DrawLineNDC(lo + eps, yline,
                      hi - eps, yline);
    */
    line->SetLineStyle(1);
    if(r < NR-1){
      line->DrawLineNDC(hi, hbo, hi, 1.-hto);
      if(pad_ratio){
	pad_ratio->cd();
	line->DrawLineNDC(hi_ratio, 0., hi_ratio, 1.);
	pad->cd();
      }
    }
    l.DrawLatex(hlo + (1.-hhi-hlo)*(0.5+r)/double(NB), yline, rlabels[r].c_str());
  }
  
  l.SetTextAngle(0);
  l.SetTextAlign(32);
  l.SetTextSize(0.035);
  l.DrawLatex((hi-hlo)/2.+hlo+0.02, hbo-0.1, "#scale[1.15]{R_{ISR}}");
}

void FitPlotter::DrawRM(vector<CategoryList> cats, TCanvas* can, TPad* pad, VS labels, TPad* pad_ratio){
  pad->cd();

  double eps = 0.0015;

  double hlo = can->GetLeftMargin();
  double hhi = can->GetRightMargin();
  double hbo = can->GetBottomMargin();
  double hto = can->GetTopMargin();

  double hlo_rat = 0.;
  double hhi_rat = 0.;
  if(pad_ratio){
    hlo_rat = pad_ratio->GetLeftMargin();
    hhi_rat = pad_ratio->GetRightMargin();
  }
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();


  vector<int> NR;// fitbin.NBins();
  vector<int> NB; //fitbin.NRBins();
  int Ncats = cats.size();
  for(int v = 0; v < Ncats; v++){
   NR.push_back(cats[v][0].GetFitBin().NBins());
   NB.push_back(cats[v][0].GetFitBin().NRBins());
}
  vector<VS> rlabels;
  int maxSize = 0;
  for(int v = 0; v < Ncats; v++){
     rlabels.push_back(VS());
     for(int r = 0; r < NR[v]; r++){
      string label = cats[v][0].GetFitBin()[r].GetRBinLabel();
      rlabels[v] += label;
      if(maxSize < labels.size())
        maxSize = label.size();
      }
  }
  
  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  
  l.SetTextSize(0.033);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(2);

  double lo = hlo;
  double hi = hlo;
  double lo_ratio = hlo_rat;
  double hi_ratio = hlo_rat;
  double yline = hbo-0.02;
  int Nl = labels.size();
  int Nbins_total = std::accumulate(NB.begin(),NB.end(),0);
  double hi_next;
  double hi_last;
  for(int r = 0; r < Ncats; r++){
    lo = hi;
    hi_last = hi;
    lo_ratio = hi_ratio;
    hi = NB[r]/double(Nbins_total)*(1.-hhi-hlo) + lo;
    hi_ratio = 1./double(Nl)*(1.-hhi_rat-hlo_rat) + lo_ratio;

    line->SetLineStyle(1);
    if(r < Nl-1){
      line->DrawLineNDC(hi, hbo-0.011*(maxSize+1), hi, 1.-hto);
      if(pad_ratio){
	pad_ratio->cd();
	line->DrawLineNDC(hi_ratio, 0., hi_ratio, 1.);
	pad->cd();
      }
    }
//this if is a hack
if (labels[r].find("0L 0J #geq 1SV") != std::string::npos)
 {
  std::cout << "Fixing label for: " << labels[r] << std::endl;
  l.SetTextAngle(90);
  l.DrawLatex((hi + hi_last)/2.-0.0075, 1 - hto - 4*eps - 0.1, labels[r].c_str());
  l.SetTextAngle(0);
 }
else
  l.DrawLatex((hi + hi_last)/2., 1 - hto - 4*eps - 0.015, labels[r].c_str());
  
  }

  l.SetTextAngle(0);
  l.SetTextAlign(32);
  l.SetTextSize(0.035);
  l.DrawLatex((hi-hlo)/2.+hlo+0.02, hbo-0.015*(maxSize+1), "#scale[1.15]{R_{ISR}}");
}
void FitPlotter::DrawRM(const FitBin& fitbin, TCanvas* can, TPad* pad, VS labels, TPad* pad_ratio){
  pad->cd();

  double eps = 0.0015;

  double hlo = can->GetLeftMargin();
  double hhi = can->GetRightMargin();
  double hbo = can->GetBottomMargin();
  double hto = can->GetTopMargin();

  double hlo_rat = 0.;
  double hhi_rat = 0.;
  if(pad_ratio){
    hlo_rat = pad_ratio->GetLeftMargin();
    hhi_rat = pad_ratio->GetRightMargin();
  }
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  int NR = fitbin.NBins();
  int NB = fitbin.NRBins();

  VS rlabels;
  int maxSize = 0;
  for(int r = 0; r < NR; r++){
    string label = fitbin[r].GetRBinLabel();
    rlabels += label;
    if(maxSize < label.size())
      maxSize = label.size();
  }
  
  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
  
  l.SetTextSize(0.033);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(2);

  double lo = hlo;
  double hi = hlo;
  double lo_ratio = hlo_rat;
  double hi_ratio = hlo_rat;
  double yline = hbo-0.02;

  int Nl = labels.size();

  for(int r = 0; r < Nl; r++){
    lo = hi;
    lo_ratio = hi_ratio;
    hi = 1./double(Nl)*(1.-hhi-hlo) + lo;
    hi_ratio = 1./double(Nl)*(1.-hhi_rat-hlo_rat) + lo_ratio;

    line->SetLineStyle(1);
    if(r < Nl-1){
      line->DrawLineNDC(hi, hbo-0.011*(maxSize+1), hi, 1.-hto);
      if(pad_ratio){
	pad_ratio->cd();
	line->DrawLineNDC(hi_ratio, 0., hi_ratio, 1.);
	pad->cd();
      }
    }

    l.DrawLatex(hlo + (1.-hhi-hlo)*(0.5+r)/double(Nl), 1 - hto - 4*eps, labels[r].c_str());
  }
  
  l.SetTextAngle(0);
  l.SetTextAlign(32);
  l.SetTextSize(0.035);
  l.DrawLatex((hi-hlo)/2.+hlo+0.02, hbo-0.015*(maxSize+1), "#scale[1.15]{R_{ISR}}");
}

TCanvas* FitPlotter::Plot2D(const string& can_name,
			   const VS& proc,
			   const CategoryTree& CT){
  RestFrames::SetStyle();

  int Nproc = proc.size();
  if(Nproc == 0)
    return nullptr;

 // cout << "Nproc: " << Nproc << endl;
  
  vector<const CategoryTree*> CatTrees;
  CT.GetListDeepest(CatTrees);
  
  int Nvis = CatTrees.size();
  
  //cout << "Nvis: " << Nvis << endl;

  if(Nvis < 1)
    return nullptr;
  
  CategoryList CatList = GetCategories();
  
  if(CatList.GetN() < 1)
    return nullptr;
  
  // Processes
  string label;
  TH2D* hist = nullptr;

  CategoryList cat = CatList.Filter(CT);
  int Ncat = cat.GetN();
  
  //cout << "Ncat: " << Ncat << endl;
  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc[i]) != 0){
      vproc = m_Strings[proc[i]];
      //cout << vproc[i] << endl;
    }
    else{
      vproc += proc[i];
      //cout << proc[i] << endl;
    }
    ProcessList procs;
    
    ProcessType type = kBkg;
    
    for(int p = 0; p < int(vproc.size()); p++){
      
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;
      
      Process pp = GetProcesses()[index];

      procs += pp;

      if(pp.Type() == kSig){
	type = kSig;
	//cout << "here" << endl;
      }
      
      if(pp.Type() == kData){
	type = kData;
      }

      for(int c = 0; c < Ncat; c++){
	if(GetHistogram2D(cat[c],pp)){
	  //cout << "here" << endl;
	  if(!hist)
	    hist = (TH2D*) GetHistogram2D(cat[c],pp)
	      ->Clone(Form("plothist2D_%s", can_name.c_str()));
	  else
	    hist->Add(GetHistogram2D(cat[c], pp));
	}
      }
    }
    //cout << "hist ptr: " << hist << endl;
    if(hist == nullptr)
      continue;
    
    if(type == kData){
      label = "Data";
    }

    if(type == kSig){
      label = GetSignalTitle(proc[i]);
    }

    if(type == kBkg){
      if(m_Title.count(proc[i]) != 0)
	label = m_Title[proc[i]];
      else
	label = proc[i];
    }
  }

  if(!hist)
    return nullptr;
  
  const FitBin& bin = cat[0].GetFitBin();

  int NR = bin.NRBins();
  int NB = bin.NBins();
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
			     Form("can_%s", can_name.c_str()),
			     700, 700);
  double eps = 0.0015;

  double hlo = 0.14;
  double hhi = 0.18;
  double hbo = 0.13;
  double hto = 0.17;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->Draw();
  can->cd();
  
  hist->Draw("colz");
  hist->GetXaxis()->SetTitle("M_{#perp}   [GeV]");
  hist->GetYaxis()->SetTitle("R_{ISR}");
  hist->GetZaxis()->SetTitle("number of events");
  
  hist->GetXaxis()->CenterTitle();
  hist->GetXaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.15);
  hist->GetXaxis()->SetLabelFont(42);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetTickSize(0.);
  
  hist->GetYaxis()->CenterTitle();
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.15);
  hist->GetYaxis()->SetLabelFont(42);
  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTickSize(0.);
  
  hist->GetZaxis()->CenterTitle();
  hist->GetZaxis()->SetTitleFont(42);
  hist->GetZaxis()->SetTitleSize(0.05);
  hist->GetZaxis()->SetTitleOffset(0.9);

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(7024);
 
  l.SetTextSize(0.025);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(3);
  line->SetLineColor(7024);
  line->SetLineStyle(7);
 
  for(int r = 0; r < NR; r++){
    line->DrawLine(hist->GetXaxis()->GetXmin(), bin[r].Rlow(),
		   hist->GetXaxis()->GetXmax(), bin[r].Rlow());
    line->DrawLine(hist->GetXaxis()->GetXmin(), bin[r].Rhigh(),
		   hist->GetXaxis()->GetXmax(), bin[r].Rhigh());

    for(int b = 0; b <  bin[r].BinEdges().size()-1; b++){
      line->DrawLine(bin[r].BinEdges()[b], bin[r].Rlow(),
		     bin[r].BinEdges()[b], bin[r].Rhigh());
    }
  }

  l.SetTextAngle(0);
  l.SetTextColor(kBlack);
  l.SetTextSize(0.04);
  l.SetTextAlign(31);
  l.DrawLatex(1.-hhi+eps*10, 1.-hto+0.02, string(CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, m_CMSLabel.c_str());
  l.SetTextAlign(33);
  l.DrawLatex(1.-hhi-eps*6, 1.-hto-0.02, label.c_str());
  
  return can;
}

TCanvas* FitPlotter::Plot1DratioSyst(const VS& proc,
 	   const Systematic& syst,
           const VS& lep_cat,
           const VS& hadS_cat,
           const VS& hadI_cat,
           const string& name,
           const VS& extra){

  //RestFrames::SetStyle();


  int Nlep  = lep_cat.size();
  int NhadS = hadS_cat.size();
  int NhadI = hadI_cat.size();
  int Nextra = -999;
  if(!extra.empty()) Nextra = extra.size();
  if(Nlep  == 0 ||
   NhadS == 0 ||
   NhadI == 0)
  return nullptr;
int catTest = 0;

  CategoryList cat = GetCategories();
  vector<CategoryList> cats;

  int Ncats;
  if(Nlep > 1) Ncats = Nlep;
  else if(NhadS > 1) Ncats = NhadS;
  else if(NhadI > 1) Ncats = NhadI;
  else if(Nextra > 1) Ncats = Nextra;
  if(Nlep  < 2 &&
   NhadS < 2 &&
   NhadI < 2  ){
    cout << "Need multiple categories for either process, lepton ID, hadS, hadI, extra argument, or systematics to create ratios" << endl;
    return nullptr;
  }

  for(int i = 0; i < Ncats; i++){
    cats.push_back(GetCategories());
  }
 // cat.Print();

cout << "total # of cats: " << cats[catTest].GetN() << endl;

  // Leptonic
  VS lep_labels;
  vector<VS> vleps;
  VS vlep;
  for(int i = 0; i < Nlep; i++){
    vleps.push_back(VS());
  if(m_Title.count(lep_cat[i]) != 0)
    lep_labels.push_back(m_Title[lep_cat[i]]);
  else
    lep_labels.push_back(lep_cat[i]);

  if(m_Strings.count(lep_cat[i]) != 0){
    int N = m_Strings[lep_cat[i]].size();
     cout << "number of strings in lep tag: " << N << endl;
  for(int j = 0; j < N; j++){
    if(Nlep > 1){
      vleps[i].push_back(m_Strings[lep_cat[i]][j]);
  }else
    vlep.push_back(m_Strings[lep_cat[i]][j]);
     cout << m_Strings[lep_cat[i]][j] << endl;
  }
  } else {
    if(Nlep > 1)
    vleps[i].push_back(lep_cat[i]);
  else 
    vlep.push_back(lep_cat[i]);
  }
  if(Nlep > 1){
    cats[i] = cats[i].FilterOR(vleps[i]);
  }
  }
if(Nlep == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vlep);
cout << cats[catTest].GetN() << endl;
cout << "leptonic cuts passed" << endl;


// Hadronic S
  VS hadS_labels;
  vector<VS> vhadSs;
  VS vhadS;
  for(int i = 0; i < NhadS; i++){
    vhadSs.push_back(VS());
  if(m_Title.count(hadS_cat[i]) != 0)
    hadS_labels.push_back(m_Title[hadS_cat[i]]);
  else
    hadS_labels.push_back(hadS_cat[i]);
cout << hadS_labels[i] << endl;
  if(m_Strings.count(hadS_cat[i]) != 0){
    int N = m_Strings[hadS_cat[i]].size();
    for(int j = 0; j < N; j++){
      if(NhadS > 1)
  vhadSs[i].push_back(m_Strings[hadS_cat[i]][j]);
else vhadS.push_back(m_Strings[hadS_cat[i]][j]);
     cout << m_Strings[hadS_cat[i]][j] << endl;
  } 
  }
  else {
    if(NhadS > 1)
    vhadSs[i].push_back(hadS_cat[i]);
  else
    vhadS.push_back(hadS_cat[i]);
  }
  if(NhadS > 1){
    cats[i] = cats[i].FilterOR(vhadSs[i]);
  }
  }
  if(NhadS == 1){
  for(int i = 0; i < Ncats; i++){
  if(Nlep == 1)
    cats[i] = cats[i].FilterOR(vhadS);
  else{
    if(!strstr(m_Title[lep_cat[i]].c_str(),"gold")){
      cats[i] = cats[i].FilterOR(vhadS); //if lep is !gold, filter with string given to macro for hadS_cat
      for(int i = 0; i < vhadS.size(); i++) cout << vhadS[i] << endl;
    }
    else{
      cats[i] = cats[i].FilterOR(m_Strings[hadS_cat[0]]); //if lep is gold, filter with VS() from m_Strings
      cout << hadS_cat[0] << endl;
    }

  }
  }
  
}

cout << cats[catTest].GetN() << endl;
 
cout << "hadronic S cuts passed" << endl;
 // Hadronic ISR
  VS hadI_labels;
  vector<VS> vhadIs;
  VS vhadI;
  for(int i = 0; i < NhadI; i++){
    vhadIs.push_back(VS());
  if(m_Title.count(hadI_cat[i]) != 0)
    hadI_labels.push_back(m_Title[hadI_cat[i]]);
  else
    hadI_labels.push_back(hadI_cat[i]);

  if(m_Strings.count(hadI_cat[i]) != 0){
    int N = m_Strings[hadI_cat[i]].size();
    for(int j = 0; j < N; j++){
      if(NhadI > 1)
  vhadIs[i].push_back(m_Strings[hadI_cat[i]][j]);
else vhadI.push_back(m_Strings[hadI_cat[i]][j]);
  } 
} else {
  if(NhadI > 1)
    vhadIs[i].push_back(hadI_cat[i]);
  else
    vhadI.push_back(hadI_cat[i]);
  }
  if(NhadI > 1){
    cats[i] = cats[i].FilterOR(vhadIs[i]);
  }
  }
  if(NhadI == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vhadI);

  if(NhadI == 1){
  for(int i = 0; i < Ncats; i++){
  if(Nlep == 1)
    cats[i] = cats[i].FilterOR(vhadI);
  else{
    if(!strstr(m_Title[lep_cat[i]].c_str(),"gold")){
      cats[i] = cats[i].FilterOR(vhadI); //if lep is !gold, filter with string given to macro for hadS_cat
    }
    else{
      cats[i] = cats[i].FilterOR(m_Strings[hadI_cat[0]]); //if lep is gold, filter with VS() from m_Strings
    }

  }
  }
  
}
  //cat = cat.FilterOR(vhadI);
cout << cats[catTest].GetN() << endl;
cout << "hadronic ISR cuts passed" << endl;
  // extra (PTISR, gammaT)
  VS extra_labels;
  vector<VS> vextras;
  VS vextra;
if(Nextra != -999){
  for(int i = 0; i < Nextra; i++){
    vextras.push_back(VS());
  if(m_Title.count(extra[i]) != 0)
    extra_labels.push_back(m_Title[extra[i]]);
  else
    extra_labels.push_back(extra[i]);

  if(m_Strings.count(extra[i]) != 0){
    int N = m_Strings[extra[i]].size();
    for(int j = 0; j < N; j++){
      if(Nextra > 1)
  vextras[i].push_back(m_Strings[extra[i]][j]);
else vextra.push_back(m_Strings[extra[i]][j]);
  } 
  }
  else {
    if(Nextra > 1)
    vextras[i].push_back(extra[i]);
  else 
    vextra.push_back(extra[i]);
  }
  if(Nextra > 1){
    cats[i] = cats[i].FilterOR(vextras[i]);
  }
  }
  if(Nextra == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vextra);
 // cat = cat.FilterOR(vextra);
}
cout << cats[catTest].GetN() << endl;
  cout << "extra cuts passed" << endl;


  for(int c = 0; c < Ncats; c++){
    int Ncat = cats[c].GetN();
    if(Ncat < 1){
      cout << "no categories found with specified lepton tags for categoryList #" << c << endl;
      return nullptr;
    }
  }
  // Processes
  VS            labels;
  vector<int>   colors;
  vector<TH1D*> hists;


int colorIdx = 0;
  

  for(int s = 0; s < 1; s++){
  	if(!(!syst)){
  		if(syst.IsUp()){
  			syst.Down();
  			s--;
  		}
  		else syst.Up();
  	}
  CategoryList cat1;
  for(int cc = 0; cc < Ncats; cc++){
    cat1 = cats[cc];
    int Ncat = cat1.GetN();
VS vproc;
//sum over fakes
for(int i = 0; i < proc.size(); i++){
	if(m_Strings.count(proc[i]) != 0)
	  vproc += m_Strings[proc[i]];
	else
	  vproc += proc[i];
}
    TH1D* histNom = nullptr;
    TH1D* hist = nullptr;
    for(int p = 0; p < int(vproc.size()); p++){
      cout << "vprocess: " << vproc[p] << endl;
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
    continue;
    
    Process pp = GetProcesses()[index];

    for(int c = 0; c < Ncat; c++){
      // cout << cat1[c].GetLabel() << " " << pp.Name() << endl;
      if(!IsFilled(cat1[c], pp,syst))
        continue;

  cout << "filled " << cat1[c].GetLabel() << " " << pp.Name() << endl;
  if(!histNom)
    histNom = (TH1D*) GetHistogram(cat1[c], pp)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
   else 
    histNom->Add(GetHistogram(cat1[c], pp));  
  if(!syst.IsDefault()){
  	if(!hist)
    hist = (TH1D*) GetHistogram(cat1[c], pp,syst)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
 else 
    hist->Add(GetHistogram(cat1[c], pp,syst));
	}
    }
  }

  if(hist == nullptr){
    cout << "hist not found for systematic: "<< syst.Label() << endl;
    continue;
  }
  if(histNom == nullptr){
    cout << "nominal hist not found" << endl;
    continue;
  }

if(syst.IsUp()){
  if(Nlep > 1)
    labels += "#splitline{"+m_Title[lep_cat[cc]]+"}{"+syst.Label()+"Up}";
  else if(NhadS > 1)
    labels += "#splitline{"+hadS_cat[cc]+"}{"+syst.Label()+"Up}";
  else if(NhadI > 1)
    labels += "#splitline{"+hadI_cat[cc]+"}{"+syst.Label()+"Up}";
  else if(Nextra > 1)
    labels += "#splitline{"+extra[cc]+"}{"+syst.Label()+"Up}";
}
else{
   if(Nlep > 1)
    labels += "#splitline{"+m_Title[lep_cat[cc]]+"}{"+syst.Label()+"Down}";
  else if(NhadS > 1)
    labels += "#splitline{"+hadS_cat[cc]+"}{"+syst.Label()+"Down}";
  else if(NhadI > 1)
    labels += "#splitline{"+hadI_cat[cc]+"}{"+syst.Label()+"Down}";
  else if(Nextra > 1)
    labels += "#splitline{"+extra[cc]+"}{"+syst.Label()+"Down}";
}
		// else labels += "nominal";
	//do ratio
	float norm = histNom->Integral()/hist->Integral();
	hist->Divide(histNom);
	hist->Scale(norm);
    
  colors.push_back(m_ColorDefault[colorIdx]);
  colorIdx++;
  hists.push_back(hist); 
 } 
}
if(hists.size() < 1) return nullptr;
int Nhist = hists.size();
int nBins;
int gBin; 

  const FitBin& bin = cats[0][0].GetFitBin();
  int NR = bin.NRBins();
  int NB = bin.NBins();
  // plotShapeComparison* psc = new plotShapeComparison(hists,colors,labels);
  //TCanvas* can = psc->formatPlots(name);
// vector<TGraphAsymmErrors*> graphs = psc->m_graphs;
// vector<TGraphAsymmErrors*> gLHs = psc->m_gLHs;
// vector<double> pvals = psc->m_pvals;

// TAxis* ax = gLHs[0]->GetHistogram()->GetXaxis();
// double ax_x1 = ax->GetBinLowEdge(1);
// double ax_x2 = ax->GetBinUpEdge(ax->GetNbins());
// gLHs[0]->GetHistogram()->GetXaxis()->Set(NB,ax_x1,ax_x2);

  VS blabels;
  for(int r = 0; r < NR; r++)
  blabels += bin[r].GetMBinLabels();


  int lmax = 0;
  for(int b = 0; b < NB; b++){
  int len = blabels[b].length();
  if(blabels[b].find("#infty") != std::string::npos)
    len -= 5;
  if(len > lmax)
    lmax = len;
  }
  string space = "";
  for(int l = 0; l < 1.2*lmax; l++)
  space += " ";

  for(int b = 0; b < NB; b++){
    if(b%2 == 1)
      hists[0]->GetXaxis()->SetBinLabel(b+1, (blabels[b]+space).c_str());
    else
      hists[0]->GetXaxis()->SetBinLabel(b+1, blabels[b].c_str());
  }
	blabels.clear();
  for(int r = 0; r < NR; r++)
  blabels += bin[r].GetRBinLabel();

  hists[0]->LabelsOption("v","X");

// double x; double y;

//  gStyle->SetPadTopMargin(0.09);
//    gStyle->SetPadRightMargin(0.25);
//    gStyle->SetPadBottomMargin(0.18);
//    gStyle->SetPadLeftMargin(0.15);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(11111111);
   TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
            Form("can_%s", name.c_str()),
            1200,700);
   // double hlo = 0.15;
   // double hhi = 0.22;
   // double hbo = 0.27;
   // double hto = 0.07;
	double hlo = 0.09;
	double hhi = 0.22;
	double hbo = 0.27;
	double hto = 0.07;
   can->SetLeftMargin(hlo);
   can->SetRightMargin(hhi);
   can->SetBottomMargin(hbo);
   can->SetTopMargin(hto);
   can->SetGridy();
   can->Draw();
   can->cd();



double hmax = hists[0]->GetMaximum();
hists[0]->Draw("hist");
hists[0]->GetXaxis()->CenterTitle();
hists[0]->GetXaxis()->SetTitleFont(42);
hists[0]->GetXaxis()->SetTitleSize(0.05);
hists[0]->GetXaxis()->SetTitleOffset(1.0);
hists[0]->GetXaxis()->SetLabelFont(42);
hists[0]->GetXaxis()->SetLabelSize(0.04);
hists[0]->GetXaxis()->SetTitle("");
hists[0]->GetXaxis()->SetTickSize(0.);
hists[0]->GetYaxis()->CenterTitle();
hists[0]->GetYaxis()->SetTitleFont(42);
hists[0]->GetYaxis()->SetTitleSize(0.04);
hists[0]->GetYaxis()->SetTitleOffset(0.85);
hists[0]->GetYaxis()->SetLabelFont(42);
hists[0]->GetYaxis()->SetTickLength(0.02);
hists[0]->GetYaxis()->SetLabelSize(0.03);
hists[0]->GetYaxis()->SetTitle("Ratio to nominal");
hists[0]->GetYaxis()->SetRangeUser(0.8, 1.1*hmax);
//hists[0]->GetXaxis()->SetLabelOffset(999);
//hists[0]->GetXaxis()->SetLabelSize(0);
//hists[0]->GetXaxis()->SetTickLength(0);
//gStyle->SetTickLength(0.);
//gStyle->SetLabelOffset(999);
//gStyle->SetLabelSize(0.);
for(int i = 0; i < Nhist; i++){
	if(i < Nhist/2){
	hists[i]->SetMarkerStyle(20);
	hists[i]->SetLineStyle(1);
	}
	else{
	hists[i]->SetMarkerStyle(21);
	hists[i]->SetLineStyle(2);
	}
	hists[i]->SetLineColor(colors[i%(Nhist/2)]);
	hists[i]->SetMarkerColor(colors[i%(Nhist/2)]);
	hists[i]->SetLineWidth(2);
	hists[i]->Draw("SAME HIST");
}

TLegend* leg = new TLegend(1.-hhi+0.01, 1.- (Nhist+1)*(1.-0.49)/9.-0.15, 0.98, 1.-hto);
leg->SetTextFont(132);
leg->SetTextSize(0.035);
leg->SetFillColor(kWhite);
leg->SetLineColor(kWhite);
leg->SetShadowColor(kWhite);

for(int i = 0; i < Nhist; i++)
leg->AddEntry(hists[i], (labels[i]).c_str(), "L");
leg->Draw("SAME");
double eps = 0.0015;

can->Modified();
can->Update();


TLatex l;
l.SetTextFont(42);
l.SetNDC();

TLine* line = new TLine();
line->SetLineWidth(2);
line->SetLineColor(kBlack);
 
// line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 l.SetTextSize(0.025);
 l.SetTextFont(42);
 l.SetTextAlign(23);
 line->SetLineWidth(1);
 double lo = hlo;
 double hi = hlo;
 double yline = hbo-0.024*lmax;
 int ib = 0;
 for(int r = 0; r < NR; r++){
   int NM = bin[r].NBins();
   lo = hi;
   hi = double(NM)/double(NB)*(1.-hhi-hlo) + lo;
   
   line->SetLineStyle(1);
   line->DrawLineNDC(lo + eps, yline,
		     lo + eps, yline + 6*eps);
   line->DrawLineNDC(hi - eps, yline,
		     hi - eps, yline + 6*eps);
   line->DrawLineNDC(lo + eps, yline,
		     hi - eps, yline);
   line->SetLineStyle(5);
   line->DrawLineNDC(hi, hbo, hi, 1.-hto);
   line->SetLineStyle(3);
   for(int b = 0; b < NM; b++){
     if(ib%2 == 1)
       line->DrawLineNDC(lo + (hi-lo)*(b+0.5)/double(NM), hbo,
			 lo + (hi-lo)*(b+0.5)/double(NM), (hbo+yline)/2.+eps);
     ib++;
   }
   
   l.DrawLatex((hi+lo)/2., yline - 8*eps, blabels[r].c_str());
 }      
 l.SetTextAlign(32);
 l.SetTextSize(0.03);
 l.SetTextFont(42);
 l.DrawLatex(hlo, (hbo+yline)/2.+eps, "M_{#perp}   [GeV] #in");
 
 l.SetTextSize(0.03);
 l.SetTextFont(42);
 l.DrawLatex(hlo, yline - 15*eps, "#scale[1.15]{R_{ISR}} #in");
 
 l.SetTextAlign(31);
 l.SetTextSize(0.04);
 l.SetTextFont(42);
 l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, "2017 MC KUEWKino");
 l.SetTextAlign(11);
 l.SetTextSize(0.04);
 l.SetTextFont(42);
 l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
 l.SetTextSize(0.05);
 
 
 string plotlabel;
 if(Nlep == 1)
   plotlabel += "#color[7014]{"+lep_labels[0]+"} + ";
 if(NhadS == 1)
   plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
 if(NhadI == 1)
   plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
 if(Nextra == 1)
   plotlabel += "#color[7024]{"+extra_labels[0]+"} + ";
 plotlabel += "p_{T}^{ISR} > 300 GeV, "+m_Title[proc[0]];

 l.SetTextColor(kBlack);
 l.SetTextAlign(13);
 l.SetTextSize(0.035);
 l.SetTextFont(42);
 l.DrawLatex(hlo+0.02, 1-hto-0.012, plotlabel.c_str());
 
 // l.DrawLatex(0.177796, 0.890436, plotlabel.c_str());
 
 return can;
 
 
}
///////


TCanvas* FitPlotter::PlotRatioSystDist(const VS& proc,
		const Systematic& syst,
		const VS& lep_cat,
		const VS& matchString,
		const string& name){
	CategoryList cat = GetCategories();
	vector<CategoryList> cats;
	int Ncats = 0;
	vector<vector<CategoryList>> cats2D;
	for(int i = 0; i < lep_cat.size(); i++){
	cout << "i: " << i << endl;
		vector<CategoryList> cat_tmp;
		for(int j = 0; j < matchString.size(); j++){ 
			cout << "j :" << j << endl; 
			cat_tmp.push_back(GetCategories());
			Ncats++;
		}
		cats2D.push_back(cat_tmp);
	}
	VS labels;
	vector<VS> filterStrings;
	for(int i = 0; i < matchString.size(); i++){
		filterStrings.push_back(VS());
		string lbl = matchString[i].substr(0,2)+"S";
		if(count(labels.begin(),labels.end(), lbl) == 0) labels.push_back(lbl);
		if(m_Strings.count(matchString[i]) != 0){
			int N = m_Strings[matchString[i]].size();
			for(int j = 0; j < N; j++){
				filterStrings[i].push_back(m_Strings[matchString[i]][j]);
			}
		}
		else{
			filterStrings[i].push_back(matchString[i]);
		}
		
	}
int Nlep = lep_cat.size();
 VS lep_labels;
  vector<VS> vleps;
  VS vlep;
  for(int i = 0; i < Nlep; i++){
    vleps.push_back(VS());
  if(m_Title.count(lep_cat[i]) != 0)
    lep_labels.push_back(m_Title[lep_cat[i]]);
  else
    lep_labels.push_back(lep_cat[i]);

  if(m_Strings.count(lep_cat[i]) != 0){
    int N = m_Strings[lep_cat[i]].size();
     cout << "number of strings in lep tag: " << N << endl;
  for(int j = 0; j < N; j++){
    if(Nlep > 1){
      vleps[i].push_back(m_Strings[lep_cat[i]][j]);
  }else
    vlep.push_back(m_Strings[lep_cat[i]][j]);
     cout << m_Strings[lep_cat[i]][j] << endl;
  }
  } else {
    if(Nlep > 1)
    vleps[i].push_back(lep_cat[i]);
  else 
    vlep.push_back(lep_cat[i]);
  }
  if(Nlep > 1){
    for(int j = 0; j < cats2D[i].size(); j++) cats2D[i][j] = cats2D[i][j].FilterOR(vleps[i]);
  }
  }
if(Nlep == 1)
  for(int i = 0; i < cats2D.size(); i++)
    for(int j = 0; j < cats2D[i].size(); j++) cats2D[i][j] = cats2D[i][j].FilterOR(vlep);
cout << "leptonic cuts passed" << endl;





for(int i = 0; i < cats2D.size(); i++){
	for(int j = 0; j < cats2D[i].size(); j++){
		cats2D[i][j] = cats2D[i][j].FilterOR(filterStrings[j]);
	}
}

//cout << "cats2D size: " << cats2D.size() << endl;
//cout << "flatten vector" << endl;
for(int i = 0; i < cats2D.size(); i++){
	for(int j = 0; j < cats2D[i].size(); j++){
		cats.push_back(cats2D[i][j]);
	}
}
//cout << "a" << endl;
//	//int Ncats = cats.size();
//cout << cats.size() << endl;
//cout << "b" << endl;
	double hmax = 0.;
	double hmin = 1.;
	vector<TH1D*> hists;
for(int s = 0; s < 1; s++){
if(!(!syst)){
	if(syst.IsUp()){
		syst.Down();
		s--;
	}
	else syst.Up();
}
CategoryList cat1;
for(int cc = 0; cc < Ncats; cc++){
cat1 = cats[cc];
int Ncat = cat1.GetN();
VS vproc;
//sum over fakes
for(int i = 0; i < proc.size(); i++){
if(m_Strings.count(proc[i]) != 0)
vproc += m_Strings[proc[i]];
else
vproc += proc[i];
}
TH1D* histNom = nullptr;
TH1D* hist = nullptr;
for(int p = 0; p < int(vproc.size()); p++){
cout << "vprocess: " << vproc[p] << endl;
int index = GetProcesses().Find(vproc[p]);
if(index < 0)
continue;
Process pp = GetProcesses()[index];
for(int c = 0; c < Ncat; c++){
if(!IsFilled(cat1[c], pp,syst))
continue;
cout << "filled " << cat1[c].GetLabel() << " " << pp.Name() << endl;
if(!histNom)
histNom = (TH1D*) GetHistogram(cat1[c], pp)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
else histNom->Add(GetHistogram(cat1[c], pp));  
if(!syst.IsDefault()){
if(!hist)
hist = (TH1D*) GetHistogram(cat1[c], pp,syst)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
else 
hist->Add(GetHistogram(cat1[c], pp,syst));
}
}
}
if(hist == nullptr){
	cout << "hist not found for systematic: "<< syst.Label() << endl;
	continue;
}
if(histNom == nullptr){
	cout << "nominal hist not found" << endl;
	continue;
}
	//do ratio
	float norm = histNom->Integral()/hist->Integral();
	hist->Divide(histNom);
	hist->Scale(norm);
			    
hists.push_back(hist); 
if(hist->GetMaximum() > hmax) hmax = hist->GetMaximum();
if(hist->GetMinimum() < hmin) hmin = hist->GetMinimum();
} 
}
TH1D* dist = new TH1D(syst.Label().c_str(),syst.Label().c_str(),80,hmin-0.05,hmax+0.05);
int Nhist = hists.size();
for(int h = 0; h < Nhist; h++){
	for(int b = 0; b < hists[h]->GetNbinsX(); b++){
		dist->Fill(hists[h]->GetBinContent(b+1));
	}


}
double eps = 0.0015;


   gStyle->SetOptTitle(0);
   //gStyle->SetOptStat(0);
   gStyle->SetOptFit(11111111);
   TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
            Form("can_%s", name.c_str()),
            1200,700);
double hlo = 0.09;
double hhi = 0.22;
double hbo = 0.17;
double hto = 0.07;
//can->SetLeftMargin(hlo);
//can->SetRightMargin(hhi);
can->SetBottomMargin(hbo);
//can->SetTopMargin(hto);
can->SetGridy();
can->Draw();
can->cd();
dist->GetXaxis()->CenterTitle();
dist->GetXaxis()->SetTitleFont(42);
dist->GetXaxis()->SetTitleSize(0.05);
dist->GetXaxis()->SetTitleOffset(1.0);
dist->GetXaxis()->SetLabelFont(42);
dist->GetXaxis()->SetLabelSize(0.04);
dist->GetXaxis()->SetTitle("Ratio to nominal");
dist->Draw("hist");

TLatex l;
l.SetTextAlign(31);
l.SetTextSize(0.04);
l.SetTextFont(42);
double x = hmin;//1.-hhi-eps*4;
double y = dist->GetMaximum();//1.-hto+0.02;
cout << "x: " << x << " y: " << y << endl;
l.DrawLatex(x+0.03, y*1.1, "2017 MC KUEWKino");
l.SetTextAlign(11);
l.SetTextSize(0.04);
l.SetTextFont(42);
l.DrawLatex(x+eps+0.03, y*1.1,"#bf{#it{CMS}} work-in-progress");
l.SetTextSize(0.05);

 string plotlabel;
  for(int i = 0; i < lep_cat.size(); i++)
    plotlabel += "#color[7014]{"+lep_labels[i]+"} + ";
  for(int i = 0; i < labels.size(); i++)
    plotlabel += "#color[7004]{"+labels[i]+"} + ";
  plotlabel += "p_{T}^{ISR} > 300 GeV, "+m_Title[proc[0]];
l.SetTextColor(kBlack);
	l.SetTextAlign(11);
	l.SetTextSize(0.035);
	l.SetTextFont(42);
	l.DrawLatex(x-0.05, y, plotlabel.c_str());
return can;
}



////////////////////////////////
////////////////////////////////
TCanvas* FitPlotter::Plot1Dratio(const string& proc,
           const VS& lep_cat,
           const VS& hadS_cat,
           const VS& hadI_cat,
           const string& name,
           const VS& extra){

  //RestFrames::SetStyle();


  int Nlep  = lep_cat.size();
  int NhadS = hadS_cat.size();
  int NhadI = hadI_cat.size();
  int Nextra = -999;
  if(!extra.empty()) Nextra = extra.size();
  if(Nlep  == 0 ||
   NhadS == 0 ||
   NhadI == 0)
  return nullptr;

int catTest = 0;

  CategoryList cat = GetCategories();
  vector<CategoryList> cats;

  int Ncats;
  if(Nlep > 1) Ncats = Nlep;
  else if(NhadS > 1) Ncats = NhadS;
  else if(NhadI > 1) Ncats = NhadI;
  else if(Nextra > 1) Ncats = Nextra;
  if(Nlep  < 2 &&
   NhadS < 2 &&
   NhadI < 2){
    cout << "Need multiple categories for either process, lepton ID, hadS, hadI, or extra argument to create ratios" << endl;
    return nullptr;
  }

  for(int i = 0; i < Ncats; i++){
    cats.push_back(GetCategories());
  }
 // cat.Print();

cout << "total # of cats: " << cats[catTest].GetN() << endl;

  // Leptonic
  VS lep_labels;
  vector<VS> vleps;
  VS vlep;
  for(int i = 0; i < Nlep; i++){
    vleps.push_back(VS());
  if(m_Title.count(lep_cat[i]) != 0)
    lep_labels.push_back(m_Title[lep_cat[i]]);
  else
    lep_labels.push_back(lep_cat[i]);

  if(m_Strings.count(lep_cat[i]) != 0){
    int N = m_Strings[lep_cat[i]].size();
     cout << "number of strings in lep tag: " << N << endl;
  for(int j = 0; j < N; j++){
    if(Nlep > 1){
      vleps[i].push_back(m_Strings[lep_cat[i]][j]);
  }else
    vlep.push_back(m_Strings[lep_cat[i]][j]);
     cout << m_Strings[lep_cat[i]][j] << endl;
  }
  } else {
    if(Nlep > 1)
    vleps[i].push_back(lep_cat[i]);
  else 
    vlep.push_back(lep_cat[i]);
  }
  if(Nlep > 1){
    cats[i] = cats[i].FilterOR(vleps[i]);
  }
  }
if(Nlep == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vlep);
cout << cats[catTest].GetN() << endl;
cout << "leptonic cuts passed" << endl;
// Hadronic S
  VS hadS_labels;
  vector<VS> vhadSs;
  VS vhadS;
  for(int i = 0; i < NhadS; i++){
    vhadSs.push_back(VS());
  if(m_Title.count(hadS_cat[i]) != 0)
    hadS_labels.push_back(m_Title[hadS_cat[i]]);
  else
    hadS_labels.push_back(hadS_cat[i]);

  if(m_Strings.count(hadS_cat[i]) != 0){
    int N = m_Strings[hadS_cat[i]].size();
    for(int j = 0; j < N; j++){
      if(NhadS > 1)
  vhadSs[i].push_back(m_Strings[hadS_cat[i]][j]);
else vhadS.push_back(m_Strings[hadS_cat[i]][j]);
     cout << m_Strings[hadS_cat[i]][j] << endl;
  } 
  }
  else {
    if(NhadS > 1)
    vhadSs[i].push_back(hadS_cat[i]);
  else
    vhadS.push_back(hadS_cat[i]);
  }
  if(NhadS > 1){
    cats[i] = cats[i].FilterOR(vhadSs[i]);
  }
  }
  if(NhadS == 1){
  for(int i = 0; i < Ncats; i++){
  if(Nlep == 1)
    cats[i] = cats[i].FilterOR(vhadS);
  else{
    if(!strstr(m_Title[lep_cat[i]].c_str(),"gold")){
      cats[i] = cats[i].FilterOR(vhadS); //if lep is !gold, filter with string given to macro for hadS_cat
      for(int i = 0; i < vhadS.size(); i++) cout << vhadS[i] << endl;
    }
    else{
      cats[i] = cats[i].FilterOR(m_Strings[hadS_cat[0]]); //if lep is gold, filter with VS() from m_Strings
      cout << hadS_cat[0] << endl;
    }

  }
  }
  
}
//cat = cat.FilterOR(vhadS);
cout << cats[catTest].GetN() << endl;
 
cout << "hadronic S cuts passed" << endl;
 // Hadronic ISR
  VS hadI_labels;
  vector<VS> vhadIs;
  VS vhadI;
  for(int i = 0; i < NhadI; i++){
    vhadIs.push_back(VS());
  if(m_Title.count(hadI_cat[i]) != 0)
    hadI_labels.push_back(m_Title[hadI_cat[i]]);
  else
    hadI_labels.push_back(hadI_cat[i]);

  if(m_Strings.count(hadI_cat[i]) != 0){
    int N = m_Strings[hadI_cat[i]].size();
    for(int j = 0; j < N; j++){
      if(NhadI > 1)
  vhadIs[i].push_back(m_Strings[hadI_cat[i]][j]);
else vhadI.push_back(m_Strings[hadI_cat[i]][j]);
  } 
} else {
  if(NhadI > 1)
    vhadIs[i].push_back(hadI_cat[i]);
  else
    vhadI.push_back(hadI_cat[i]);
  }
  if(NhadI > 1){
    cats[i] = cats[i].FilterOR(vhadIs[i]);
  }
  }
  if(NhadI == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vhadI);

  if(NhadI == 1){
  for(int i = 0; i < Ncats; i++){
  if(Nlep == 1)
    cats[i] = cats[i].FilterOR(vhadI);
  else{
    if(!strstr(m_Title[lep_cat[i]].c_str(),"gold")){
      cats[i] = cats[i].FilterOR(vhadI); //if lep is !gold, filter with string given to macro for hadS_cat
    }
    else{
      cats[i] = cats[i].FilterOR(m_Strings[hadI_cat[0]]); //if lep is gold, filter with VS() from m_Strings
    }

  }
  }
  
}
  //cat = cat.FilterOR(vhadI);
cout << cats[catTest].GetN() << endl;
cout << "hadronic ISR cuts passed" << endl;
  // extra (PTISR, gammaT)
  VS extra_labels;
  vector<VS> vextras;
  VS vextra;
if(Nextra != -999){
  for(int i = 0; i < Nextra; i++){
    vextras.push_back(VS());
  if(m_Title.count(extra[i]) != 0)
    extra_labels.push_back(m_Title[extra[i]]);
  else
    extra_labels.push_back(extra[i]);

  if(m_Strings.count(extra[i]) != 0){
    int N = m_Strings[extra[i]].size();
    for(int j = 0; j < N; j++){
      if(Nextra > 1)
  vextras[i].push_back(m_Strings[extra[i]][j]);
else vextra.push_back(m_Strings[extra[i]][j]);
  } 
  }
  else {
    if(Nextra > 1)
    vextras[i].push_back(extra[i]);
  else 
    vextra.push_back(extra[i]);
  }
  if(Nextra > 1){
    cats[i] = cats[i].FilterOR(vextras[i]);
  }
  }
  if(Nextra == 1)
  for(int i = 0; i < Ncats; i++)
    cats[i] = cats[i].FilterOR(vextra);
 // cat = cat.FilterOR(vextra);
}
cout << cats[catTest].GetN() << endl;
  cout << "extra cuts passed" << endl;


  for(int c = 0; c < Ncats; c++){
    int Ncat = cats[c].GetN();
    if(Ncat < 1){
      cout << "no categories found with specified lepton tags for categoryList #" << c << endl;
      return nullptr;
    }
  }
  // Processes
  VS            labels;
  vector<int>   colors;
  vector<TH1D*> hists;


  CategoryList cat1;
  for(int cc = 0; cc < Ncats; cc++){
    cat1 = cats[cc];
    int Ncat = cat1.GetN();
VS vproc;
    if(m_Strings.count(proc) != 0)
      vproc = m_Strings[proc];
    else
      vproc += proc;
    TH1D* hist = nullptr;
    for(int p = 0; p < int(vproc.size()); p++){
      cout << "vprocess: " << vproc[p] << endl;
      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
    continue;
    
    Process pp = GetProcesses()[index];

    for(int c = 0; c < Ncat; c++){
      // cout << cat1[c].GetLabel() << " " << pp.Name() << endl;
      if(!IsFilled(cat1[c], pp))
        continue;

  cout << "filled " << cat1[c].GetLabel() << " " << pp.Name() << endl;

  if(!hist){
    hist = (TH1D*) GetHistogram(cat1[c], pp)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
  //  if(hist->GetBinContent(0) != 0) cout << hist->GetTitle() << " has bin #0 filled with " << hist->GetBinContent(0) << endl;
  } else {
    hist->Add(GetHistogram(cat1[c], pp));
    //if(hist->GetBinContent(0) != 0) cout << GetHistogram(cat1[c], pp)->GetTitle() << " has bin #0 filled with " << GetHistogram(cat1[c], pp)->GetBinContent(0) << endl;
  }
    }
  }

  if(hist == nullptr){
    cout << "hist not found" << endl;
    continue;
  }


  else if(Nlep > 1)
    labels += m_Title[lep_cat[cc]];
  else if(NhadS > 1)
    labels += hadS_cat[cc];
  else if(NhadI > 1)
    labels += hadI_cat[cc];
  else if(Nextra > 1)
    labels += extra[cc];
    
  colors.push_back(m_ColorDefault[cc]);
  hists.push_back(hist); //one hist per cat group
 } 
int Nhist = hists.size();

int nBins;
int gBin; 

  const FitBin& bin = cats[0][0].GetFitBin();
  int NR = bin.NRBins();
  int NB = bin.NBins();
  plotShapeComparison* psc = new plotShapeComparison(hists,colors,labels);
  //TCanvas* can = psc->formatPlots(name);
vector<TGraphAsymmErrors*> graphs = psc->m_graphs;
vector<TGraphAsymmErrors*> gLHs = psc->m_gLHs;
vector<double> pvals = psc->m_pvals;

TAxis* ax = gLHs[0]->GetHistogram()->GetXaxis();
double ax_x1 = ax->GetBinLowEdge(1);
double ax_x2 = ax->GetBinUpEdge(ax->GetNbins());
gLHs[0]->GetHistogram()->GetXaxis()->Set(NB,ax_x1,ax_x2);

  VS blabels_RISR;
  VS blabels_Mperp;
  for(int r = 0; r < NR; r++)
  blabels_Mperp += bin[r].GetMBinLabels();


  int lmax = 0;
  for(int b = 0; b < NB; b++){
  int len = blabels_Mperp[b].length();
  if(blabels_Mperp[b].find("#infty") != std::string::npos)
    len -= 5;
  if(len > lmax)
    lmax = len;
  }
  string space = "";
  for(int l = 0; l < 1.2*lmax; l++)
  space += " ";

  for(int r = 0; r < NR; r++)
  blabels_RISR += bin[r].GetRBinLabel();

  //gLHs[0]->GetHistogram()->LabelsOption("v","X");

double x; double y;
//for(int i = 0; i < graphs[0]->GetN(); i++) { graphs[0]->GetPoint(i,x,y); cout << "point #" << i << "| x: " << x << " y: " << y << endl;}

  gStyle->SetPadTopMargin(0.09);
    gStyle->SetPadRightMargin(0.25);
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadLeftMargin(0.15);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(11111111);
   TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
            Form("can_%s", name.c_str()),
            0,45,1337, 738);
   double hlo = 0.15;
   double hhi = 0.22;
   double hbo = 0.27;
   double hto = 0.07;
   can->SetLeftMargin(hlo);
   can->SetRightMargin(hhi);
   can->SetBottomMargin(hbo);
   can->SetTopMargin(hto);
   //can->SetGridy();
   can->Draw();
   can->cd();
 TPad* pRatio = new TPad("pRatio","pRatio",0,0.3,1.,1.);
  pRatio->SetGridy();
  pRatio->SetGridx();
  pRatio->SetBottomMargin(0.02);
  pRatio->Draw();
  pRatio->cd();
  can->Update();


  pRatio->cd();


   double hmax = graphs[0]->GetMaximum();
TMultiGraph* mg = new TMultiGraph("ratios","ratios");
   mg->GetYaxis()->SetRangeUser(0.0, 1.1*hmax);
   //mg->Draw("AP");
   mg->GetXaxis()->CenterTitle();
   mg->GetXaxis()->SetTitleFont(42);
   mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetTitleOffset(1.0);
   mg->GetXaxis()->SetLabelFont(42);
   mg->GetXaxis()->SetLabelSize(0.04);
   mg->GetXaxis()->SetTitle("");
   mg->GetXaxis()->SetTickSize(0.);
   mg->GetYaxis()->CenterTitle();
   mg->GetYaxis()->SetTitleFont(42);
   mg->GetYaxis()->SetTitleSize(0.04);
   mg->GetYaxis()->SetTitleOffset(0.85);
   mg->GetYaxis()->SetLabelFont(42);
   mg->GetYaxis()->SetTickLength(0.02);
   mg->GetYaxis()->SetLabelSize(0.03);
   mg->GetYaxis()->SetTitle("Ratio to total");
   mg->GetXaxis()->SetLabelOffset(999);
   mg->GetXaxis()->SetLabelSize(0);
   mg->GetXaxis()->SetTickLength(0);
   gStyle->SetTickLength(0.);
   gStyle->SetLabelOffset(999);
   gStyle->SetLabelSize(0.);
   for(int i = 0; i < Nhist; i++){
   graphs[i]->SetLineColor(colors[i]);
   graphs[i]->SetMarkerColor(colors[i]);
   graphs[i]->SetLineWidth(1.0);
   graphs[i]->SetMarkerStyle(20+i);
   graphs[i]->SetLineStyle(i+1);
   mg->Add(graphs[i]);
   }
mg->Draw("AP");
can->Modified();
can->Update();
can->cd();
  TPad* pLH = new TPad("pLikelihood","pLikelihood",0,0,1,0.3);
  pLH->SetGridy();
  pLH->SetGridx();
  pLH->SetTopMargin(0.);
  pLH->Draw();
  pLH->cd();
  pLH->Update();
  can->Update();
TMultiGraph* mgLH = new TMultiGraph("LH","LH");
//gLHs[0]->Draw("AP");
   mgLH->GetXaxis()->CenterTitle();
//   gLHs[0]->GetYaxis()->SetRangeUser(0.0, 1.1*gLHs[0]->GetMaximum());
   mgLH->GetXaxis()->SetTitleFont(42);
   mgLH->GetXaxis()->SetTitleSize(0.3);
   mgLH->GetXaxis()->SetTitleOffset(1.0);
   mgLH->GetXaxis()->SetLabelFont(42);
   mgLH->GetXaxis()->SetLabelSize(0.3);
   mgLH->GetXaxis()->SetTitle("Bin Number");
   mgLH->GetXaxis()->SetTickSize(0.03);
   mgLH->GetYaxis()->SetTickSize(0.02);
   mgLH->GetYaxis()->CenterTitle();
   mgLH->GetYaxis()->SetTitleFont(42);
   mgLH->GetYaxis()->SetTitleSize(0.3);
   mgLH->GetYaxis()->SetTitleOffset(0.37);
   mgLH->GetYaxis()->SetLabelFont(42);
   mgLH->GetYaxis()->SetLabelSize(0.07);
   mgLH->GetYaxis()->SetTitle("Sigma");
   gStyle->SetTickLength(0.02,"X");
   gStyle->SetTickLength(0.02,"Y");
   gStyle->SetTitleOffset(0.5,"Y");
   gStyle->SetTitleOffset(0.9,"X");
   gStyle->SetLabelOffset(0.01,"X");
   gStyle->SetLabelSize(0.07,"X");
   gStyle->SetLabelSize(0.07,"Y");
   gStyle->SetTitleSize(0.07,"X");
   gStyle->SetTitleSize(0.07,"Y");
for(int i = 0; i < gLHs.size(); i++){
   gLHs[i]->SetLineColor(colors[i]);
   gLHs[i]->SetMarkerColor(colors[i]);
   gLHs[i]->SetLineWidth(1.0);
   gLHs[i]->SetMarkerStyle(20+i);
   gLHs[i]->SetLineStyle(i+1);
   mgLH->Add(gLHs[i]);
 //gLHs[i]->Draw("AP same");
}
mgLH->Draw("AP");
can->Modified();
can->Update();
can->cd();
pRatio->cd();
   TLegend* leg = new TLegend(1.-hhi+0.03, 1.- (Nhist+1)*(1.-0.49)/9.-0.15, 0.98, 1.-hto);
   leg->SetTextFont(132);
   leg->SetTextSize(0.035);
   leg->SetFillColor(kWhite);
   leg->SetLineColor(kWhite);
   leg->SetShadowColor(kWhite);
   // if(hist_data)
   // leg->AddEntry(hist_data, "data");
   // leg->AddEntry(gr, "total uncertainty","F");
   for(int i = 0; i < Nhist; i++)
   leg->AddEntry(graphs[i], ("#splitline{"+labels[i]+"}{Integral: "+std::to_string(hists[i]->Integral(1,hists[i]->GetNbinsX()))+"}").c_str(), "LP");
  leg->Draw("SAME");
  double eps = 0.0015;

can->Modified();
can->Update();
can->cd();
pLH->cd();
   TLegend* legLH = new TLegend(1.-hhi-0.01, 1.- (Nhist+1)*(1.-0.49)/9.-0.4, 0.98, 1.-hto);
   legLH->SetTextFont(132);
   legLH->SetTextSize(0.06);
   legLH->SetFillColor(kWhite);
   legLH->SetLineColor(kWhite);
   legLH->SetShadowColor(kWhite);
   for(int i = 0; i < gLHs.size(); i++)
   legLH->AddEntry(gLHs[i], psc->m_LHlabels[i].c_str(), "LP");
  legLH->Draw("same");

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
can->Modified();
can->Update();
can->cd();
pLH->cd();  
// line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
  l.SetTextSize(0.03);
  l.SetTextFont(132);
  l.SetTextAlign(13);
//  line->SetLineWidth(1);
  double lo = hlo;
  double hi = hlo;
  double yline = hbo-0.024*lmax;

can->Modified();
can->Update();
can->cd();
pRatio->cd();
int ib = 0;
for(int r = 0; r < NR; r++){
int NM = bin[r].NBins();
  for(int m = 0; m < NM; m++){
  string binLabel = "Bin #"+std::to_string(ib)+" | RISR: "+blabels_RISR[r]+" | M_{#perp}  : "+bin[r].GetMBinLabels()[m]; 
  //string binLabel = "Bin #"+std::to_string(ib)+" | RISR: "+blabels_RISR[r]+" | M_{#perp}  : "+blabels_Mperp[m]; 
  l.DrawLatex(1.-hhi+0.01, 1.- (Nhist+1)*(1.-0.49)/9. - 0.23- 0.05*ib,binLabel.c_str());
   //TLegend* leg = new TLegend(1.-hhi+0.03, 1.- (Nhist+1)*(1.-0.49)/9.-0.2, 0.98, 1.-hto);
  //1.-hhi+0.01, 1.- (Nhist+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005
  ib++;
  } 
}      
can->Modified();
can->Update();
can->cd();
pLH->cd();


can->Modified();
can->Update();
can->cd();
pRatio->cd();

  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, "2017 MC KUEWKino");
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
  l.SetTextSize(0.05);


  string plotlabel;
  if(Nlep == 1)
    plotlabel += "#color[7014]{"+lep_labels[0]+"} + ";
  if(NhadS == 1)
    plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
  if(NhadI == 1)
    plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
  if(Nextra == 1)
    plotlabel += "#color[7024]{"+extra_labels[0]+"} + ";
  plotlabel += "p_{T}^{ISR} > 300 GeV, "+m_Title[proc];

  l.SetTextColor(kBlack);
  l.SetTextAlign(13);
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(0.177796, 0.890436, plotlabel.c_str());

  return can;
  

  }
////////////
TCanvas* FitPlotter::Plot1DratioProc(const VS& proc,
           const string& lep_cat,
           const string& hadS_cat,
           const string& hadI_cat,
           const string& name,
           const string& extra){

 // RestFrames::SetStyle();


  int Nproc = proc.size();
  if(Nproc  == 0)
  return nullptr;

  


  CategoryList cat = GetCategories();
  // vector<CategoryList> cats;

  int Ncats = Nproc;

  // for(int i = 0; i < Ncats; i++){
  //   cats.push_back(GetCategories());
  // }
  //cat.Print();

  // Leptonic
  VS lep_labels;
  // vector<VS> vleps;
  VS vlep;
cout << "N original cats: " << cat.GetN() << endl;
  if(m_Title.count(lep_cat) != 0)
    lep_labels.push_back(m_Title[lep_cat]);
  else
    lep_labels.push_back(lep_cat);

  if(m_Strings.count(lep_cat) != 0){
    int N = m_Strings[lep_cat].size();
    // cout << "number of strings in lep tag: " << N << endl;
  for(int j = 0; j < N; j++){
    vlep.push_back(m_Strings[lep_cat][j]);
    // cout << m_Strings[lep_cat[i]][j] << endl;
  }
  } else {
    vlep.push_back(lep_cat);
  }
  // for(int i = 0; i < Ncats; i++)
  //   cats[i] = cats[i].FilterOR(vlep);
  cat = cat.FilterOR(vlep);
for(int i = 0; i < int(vlep.size()); i++) cout << vlep[i] << endl;
cout << "leptonic cuts passed, Ncats: " << cat.GetN() << endl;
// Hadronic S
  VS hadS_labels;
  // vector<VS> vhadSs;
  VS vhadS;
  if(m_Title.count(hadS_cat) != 0)
    hadS_labels.push_back(m_Title[hadS_cat]);
  else
    hadS_labels.push_back(hadS_cat);

  if(m_Strings.count(hadS_cat) != 0){
    int N = m_Strings[hadS_cat].size();
    for(int j = 0; j < N; j++){
      vhadS.push_back(m_Strings[hadS_cat][j]);
  } 
  }
  else {
    vhadS.push_back(hadS_cat);
  }
  
 //   for(int i = 0; i < Ncats; i++){
  //  cats[i] = cats[i].FilterOR(vhadS);
  // }
  
// }
cat = cat.FilterOR(vhadS);
 
cout << "hadronic S cuts passed, Ncats: " << cat.GetN() << endl;
 // Hadronic ISR
  VS hadI_labels;
  // vector<VS> vhadIs;
  VS vhadI;

  if(m_Title.count(hadI_cat) != 0)
    hadI_labels.push_back(m_Title[hadI_cat]);
  else
    hadI_labels.push_back(hadI_cat);

  if(m_Strings.count(hadI_cat) != 0){
    int N = m_Strings[hadI_cat].size();
    for(int j = 0; j < N; j++){
      vhadI.push_back(m_Strings[hadI_cat][j]);
  } 
} else {
    vhadI.push_back(hadI_cat);
  }
  
  
  // // if(NhadI == 1)
  // for(int i = 0; i < Ncats; i++)
  //   cats[i] = cats[i].FilterOR(vhadI);
cat.FilterOR(vhadI);
  
  
// }
  //cat = cat.FilterOR(vhadI);
cout << "hadronic ISR cuts passed, Ncats: " << cat.GetN() << endl;
  // extra (PTISR, gammaT)
  VS extra_labels;
  // vector<VS> vextras;
  VS vextra;

  if(m_Title.count(extra) != 0)
    extra_labels.push_back(m_Title[extra]);
  else
    extra_labels.push_back(extra);

  if(m_Strings.count(extra) != 0){
    int N = m_Strings[extra].size();
    for(int j = 0; j < N; j++){
    vextra.push_back(m_Strings[extra][j]);
  } 
  }
  else {
    vextra.push_back(extra);
  }
  
  // for(int i = 0; i < Ncats; i++)
  //   cats[i] = cats[i].FilterOR(vextra);
cat.FilterOR(vextra);

  cout << "extra cuts passed, Ncats: " << cat.GetN() << endl;


    int Ncat = cat.GetN();
    if(Ncat < 1){
      cout << "no categories found"  << endl;
      return nullptr;
    }
  
  // Processes
  VS            labels;
  vector<int>   colors;
  vector<TH1D*> hists;


  for(int i = 0; i < Nproc; i++){
    TH1D* hist = nullptr;
    VS vproc;
cout << "process: " << proc[i] << endl;
      if(m_Strings.count(proc[i]) != 0)
        vproc = m_Strings[proc[i]];
      else
        vproc += proc[i];
    for(int p = 0; p < int(vproc.size()); p++){
      // TH1D* hist = nullptr;

        cout << "vprocess: " << vproc[p] << endl;
        int index = GetProcesses().Find(vproc[p]);
        if(index < 0)
      continue;

      Process pp = GetProcesses()[index];

      for(int c = 0; c < Ncat; c++){
          // cout << cat1[c].GetLabel() << " " << pp.Name() << endl;
          if(!IsFilled(cat[c], pp))
            continue;

        cout << "filled " << cat[c].GetLabel() << " " << pp.Name() << endl;

        if(!hist){
    cout << "first histogram: " << cat[c].GetLabel() << " " << pp.Name() << endl;
          hist = (TH1D*) GetHistogram(cat[c], pp)->Clone(Form("plothist_%d_%s", 0, name.c_str()));
        } else {
          hist->Add(GetHistogram(cat[c], pp));
      cout << "Added histogram: " << cat[c].GetLabel() << " " << pp.Name() << endl;
          }
      }
    }

    if(hist == nullptr){
      cout << "hist not found" << endl;
      continue;
    }

    labels += proc[i];


    colors.push_back(m_ColorDefault[i]);
    hist->Sumw2(); 
    hists.push_back(hist); //one hist per process
  } 




  // int Nsig = hists_sig.size();

  //hists[0]->Scale(1/hists[0]->Integral());


  // sort the histograms by integral (N^2/2 brute force)
  int Nhist = hists.size();

  

  int nBins;
int gBin; 

  const FitBin& bin = cat[0].GetFitBin();

  int NR = bin.NRBins();
  int NB = bin.NBins();

  plotShapeComparison* psc = new plotShapeComparison(hists,colors,labels);
  //TCanvas* can = psc->formatPlots(name);
vector<TGraphAsymmErrors*> graphs = psc->m_graphs;
vector<TGraphAsymmErrors*> gLHs = psc->m_gLHs;
vector<double> pvals = psc->m_pvals;

TAxis* ax = gLHs[0]->GetHistogram()->GetXaxis();
double ax_x1 = ax->GetBinLowEdge(1);
double ax_x2 = ax->GetBinUpEdge(ax->GetNbins());
gLHs[0]->GetHistogram()->GetXaxis()->Set(NB,ax_x1,ax_x2);

  VS blabels_RISR;
  VS blabels_Mperp;
  for(int r = 0; r < NR; r++)
  blabels_Mperp += bin[r].GetMBinLabels();

for(int i = 0; i < blabels_Mperp.size(); i++) cout << "Mperp bin #: " << i << " bin label: " << blabels_Mperp[i] << endl;

  int lmax = 0;
  for(int b = 0; b < NB; b++){
  int len = blabels_Mperp[b].length();
  if(blabels_Mperp[b].find("#infty") != std::string::npos)
    len -= 5;
  if(len > lmax)
    lmax = len;
  }
  string space = "";
  for(int l = 0; l < 1.2*lmax; l++)
  space += " ";

  for(int r = 0; r < NR; r++)
  blabels_RISR += bin[r].GetRBinLabel();

  //gLHs[0]->GetHistogram()->LabelsOption("v","X");

double x; double y;
//for(int i = 0; i < graphs[0]->GetN(); i++) { graphs[0]->GetPoint(i,x,y); cout << "point #" << i << "| x: " << x << " y: " << y << endl;}

  gStyle->SetPadTopMargin(0.09);
    gStyle->SetPadRightMargin(0.25);
    gStyle->SetPadBottomMargin(0.18);
    gStyle->SetPadLeftMargin(0.15);
   gStyle->SetOptTitle(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(11111111);
   TCanvas* can = new TCanvas(Form("can_%s", name.c_str()),
            Form("can_%s", name.c_str()),
            0,45,1337, 738);
   double hlo = 0.15;
   double hhi = 0.22;
   double hbo = 0.27;
   double hto = 0.07;
   can->SetLeftMargin(hlo);
   can->SetRightMargin(hhi);
   can->SetBottomMargin(hbo);
   can->SetTopMargin(hto);
   //can->SetGridy();
   can->Draw();
   can->cd();
 TPad* pRatio = new TPad("pRatio","pRatio",0,0.3,1.,1.);
  pRatio->SetGridy();
  pRatio->SetGridx();
  pRatio->SetBottomMargin(0.02);
  pRatio->Draw();
  pRatio->cd();
  can->Update();


  pRatio->cd();


   double hmax = graphs[0]->GetMaximum();
TMultiGraph* mg = new TMultiGraph("ratios","ratios");
   mg->GetYaxis()->SetRangeUser(0.0, 1.1*hmax);
   //mg->Draw("AP");
   mg->GetXaxis()->CenterTitle();
   mg->GetXaxis()->SetTitleFont(42);
   mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetXaxis()->SetTitleOffset(1.0);
   mg->GetXaxis()->SetLabelFont(42);
   mg->GetXaxis()->SetLabelSize(0.04);
   mg->GetXaxis()->SetTitle("");
   mg->GetXaxis()->SetTickSize(0.);
   mg->GetYaxis()->CenterTitle();
   mg->GetYaxis()->SetTitleFont(42);
   mg->GetYaxis()->SetTitleSize(0.04);
   mg->GetYaxis()->SetTitleOffset(0.85);
   mg->GetYaxis()->SetLabelFont(42);
   mg->GetYaxis()->SetTickLength(0.02);
   mg->GetYaxis()->SetLabelSize(0.03);
   mg->GetYaxis()->SetTitle("Ratio to total");
   mg->GetXaxis()->SetLabelOffset(999);
   mg->GetXaxis()->SetLabelSize(0);
   mg->GetXaxis()->SetTickLength(0);
   gStyle->SetTickLength(0.);
   gStyle->SetLabelOffset(999);
   gStyle->SetLabelSize(0.);
   for(int i = 0; i < Nhist; i++){
   graphs[i]->SetLineColor(colors[i]);
   graphs[i]->SetMarkerColor(colors[i]);
   graphs[i]->SetLineWidth(1.0);
   graphs[i]->SetMarkerStyle(20);
   graphs[i]->SetLineStyle(i+1);
   mg->Add(graphs[i]);
   }
mg->Draw("AP");
can->Modified();
can->Update();
can->cd();
  TPad* pLH = new TPad("pLikelihood","pLikelihood",0,0,1,0.3);
  pLH->SetGridy();
  pLH->SetGridx();
  pLH->SetTopMargin(0.);
  pLH->Draw();
  pLH->cd();
  pLH->Update();
  can->Update();
TMultiGraph* mgLH = new TMultiGraph("LH","LH");
//gLHs[0]->Draw("AP");
   mgLH->GetXaxis()->CenterTitle();
//   gLHs[0]->GetYaxis()->SetRangeUser(0.0, 1.1*gLHs[0]->GetMaximum());
   mgLH->GetXaxis()->SetTitleFont(42);
   mgLH->GetXaxis()->SetTitleSize(0.3);
   mgLH->GetXaxis()->SetTitleOffset(1.0);
   mgLH->GetXaxis()->SetLabelFont(42);
   mgLH->GetXaxis()->SetLabelSize(0.3);
   mgLH->GetXaxis()->SetTitle("Bin Number");
   mgLH->GetXaxis()->SetTickSize(0.03);
   mgLH->GetYaxis()->SetTickSize(0.02);
   mgLH->GetYaxis()->CenterTitle();
   mgLH->GetYaxis()->SetTitleFont(42);
   mgLH->GetYaxis()->SetTitleSize(0.3);
   mgLH->GetYaxis()->SetTitleOffset(0.37);
   mgLH->GetYaxis()->SetLabelFont(42);
   mgLH->GetYaxis()->SetLabelSize(0.07);
   mgLH->GetYaxis()->SetTitle("Sigma");
   gStyle->SetTickLength(0.02,"X");
   gStyle->SetTickLength(0.02,"Y");
   gStyle->SetTitleOffset(0.5,"Y");
   gStyle->SetTitleOffset(0.9,"X");
   gStyle->SetLabelOffset(0.01,"X");
   gStyle->SetLabelSize(0.07,"X");
   gStyle->SetLabelSize(0.07,"Y");
   gStyle->SetTitleSize(0.07,"X");
   gStyle->SetTitleSize(0.07,"Y");
for(int i = 0; i < gLHs.size(); i++){
   gLHs[i]->SetLineColor(colors[i]);
   gLHs[i]->SetMarkerColor(colors[i]);
   gLHs[i]->SetLineWidth(1.0);
   gLHs[i]->SetMarkerStyle(20);
   gLHs[i]->SetLineStyle(i+1);
   mgLH->Add(gLHs[i]);
 //gLHs[i]->Draw("AP same");
}
mgLH->Draw("AP");
can->Modified();
can->Update();
can->cd();
pRatio->cd();
   TLegend* leg = new TLegend(1.-hhi+0.03, 1.- (Nhist+1)*(1.-0.49)/9.-0.15, 0.98, 1.-hto);
   leg->SetTextFont(132);
   leg->SetTextSize(0.035);
   leg->SetFillColor(kWhite);
   leg->SetLineColor(kWhite);
   leg->SetShadowColor(kWhite);
   // if(hist_data)
   // leg->AddEntry(hist_data, "data");
   // leg->AddEntry(gr, "total uncertainty","F");
   for(int i = 0; i < Nhist; i++)
   leg->AddEntry(graphs[i], ("#splitline{"+labels[i]+"}{Integral: "+std::to_string(hists[i]->Integral(1,hists[i]->GetNbinsX()))+"}").c_str(), "LP");
  leg->Draw("SAME");
  double eps = 0.0015;

can->Modified();
can->Update();
can->cd();
pLH->cd();
   TLegend* legLH = new TLegend(1.-hhi-0.01, 1.- (Nhist+1)*(1.-0.49)/9.-0.4, 0.98, 1.-hto);
   legLH->SetTextFont(132);
   legLH->SetTextSize(0.06);
   legLH->SetFillColor(kWhite);
   legLH->SetLineColor(kWhite);
   legLH->SetShadowColor(kWhite);
   for(int i = 0; i < gLHs.size(); i++)
   legLH->AddEntry(gLHs[i], psc->m_LHlabels[i].c_str(), "LP");
  legLH->Draw("same");

  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);
can->Modified();
can->Update();
can->cd();
pLH->cd();  
// line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
  l.SetTextSize(0.03);
  l.SetTextFont(132);
  l.SetTextAlign(13);
//  line->SetLineWidth(1);
  double lo = hlo;
  double hi = hlo;
  double yline = hbo-0.024*lmax;

can->Modified();
can->Update();
can->cd();
pRatio->cd();
int ib = 0;
for(int r = 0; r < NR; r++){
int NM = bin[r].NBins();
  for(int m = 0; m < bin[r].GetMBinLabels().size(); m++){
  string binLabel = "Bin #"+std::to_string(ib)+" | RISR: "+blabels_RISR[r]+" | M_{#perp}  : "+bin[r].GetMBinLabels()[m]; 
  l.DrawLatex(1.-hhi+0.01, 1.- (Nhist+1)*(1.-0.49)/9. - 0.23- 0.05*ib,binLabel.c_str());
   //TLegend* leg = new TLegend(1.-hhi+0.03, 1.- (Nhist+1)*(1.-0.49)/9.-0.2, 0.98, 1.-hto);
  //1.-hhi+0.01, 1.- (Nhist+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005
  ib++;
  } 
}      
can->Modified();
can->Update();
can->cd();
pLH->cd();


can->Modified();
can->Update();
can->cd();
pRatio->cd();

  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, "2017 MC KUEWKino");
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02,"#bf{#it{CMS}} work-in-progress");
  l.SetTextSize(0.05);


  string plotlabel;
  //if(Nlep == 1)
    plotlabel += "#color[7014]{"+lep_labels[0]+"} + ";
  //if(NhadS == 1)
    plotlabel += "#color[7004]{"+hadS_labels[0]+"} + ";
  //if(NhadI == 1)
    plotlabel += "#color[7024]{"+hadI_labels[0]+"} + ";
  //if(Nextra == 1)
    plotlabel += "#color[7024]{"+extra_labels[0]+"} + ";
  plotlabel += "p_{T}^{ISR} > 300 GeV";
  l.SetTextColor(kBlack);
  l.SetTextAlign(13);
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(0.177796, 0.890436, plotlabel.c_str());

  return can;
  

}

void FitPlotter::InitializeRecipes(){
  // Processes

  m_Strings["Data"] = VS().a("data_obs");
  
  m_Title["ttbar"] = "t #bar{t} + X";
  //m_Color["ttbar"] = 7011;
  m_Color["ttbar"] = 8003;

  m_Title["ST"] = "single top";
  //m_Color["ST"] = 7010;
  m_Color["ST"] = 8004;

  m_Title["DB"] = "dibosons";
  //m_Color["DB"] = 7051;
  m_Color["DB"] = 8002;

  m_Title["TB"] = "tribosons";
  //m_Color["TB"] = 7050;
  m_Color["TB"] = 8006;

  m_Title["ZDY"] = "Z / #gamma* + jets";
  //m_Color["ZDY"] = 7000;
  m_Color["ZDY"] = 8005;

  m_Title["Wjets"] = "W + jets";
  //m_Color["Wjets"] = 7001;
  m_Color["Wjets"] = 8001; 

  m_Title["ZInv"] = "Z Inv";
  m_Color["ZInv"] = 8005;
  
  m_Title["DY"] = "DY";
  m_Color["DY"] = 8006;
  
  m_Title["Cascades_220"] = "Cascades 220";
  m_Color["Cascades_220"] = 8007;
  m_Title["Cascades_260"] = "Cascades 260";
  m_Color["Cascades_260"] = 8008;
  m_Title["Cascades_270"] = "Cascades 270";
  m_Color["Cascades_270"] = 8009;

  VS s_Fakes_HF;
  s_Fakes_HF.a("Fakes_elf1G").a("Fakes_muf1G").a("Fakes_elf1S")
    .a("Fakes_muf1S").a("Fakes_elf1B").a("Fakes_muf1B")
    .a("Fakes_elf1Gelf1G").a("Fakes_muf1Gmuf1G").a("Fakes_elf1Gelf1S").a("Fakes_muf1Gmuf1S")
    .a("Fakes_elf1Gelf1B").a("Fakes_muf1Gmuf1B").a("Fakes_elf1Self1S").a("Fakes_muf1Smuf1S")
    .a("Fakes_elf1Self1B").a("Fakes_muf1Smuf1B")
    .a("Fakes_elf0Gelf1G").a("Fakes_muf0Gmuf1G").a("Fakes_elf0Gelf1S").a("Fakes_muf0Gmuf1S")
    .a("Fakes_elf1Gelf0S").a("Fakes_muf1Gmuf0S").a("Fakes_elf0Gelf1B").a("Fakes_muf0Gmuf1B")
    .a("Fakes_elf1Gelf0B").a("Fakes_muf1Gmuf0B").a("Fakes_elf0Self1S").a("Fakes_muf0Smuf1S")
    .a("Fakes_elf0Self1B").a("Fakes_muf0Smuf1B").a("Fakes_elf1Self0B").a("Fakes_muf1Smuf0B")
    .a("Fakes_elf1Gmuf1G").a("Fakes_elf1Gmuf1S").a("Fakes_elf1Smuf1G").a("Fakes_elf1Gmuf1B")
    .a("Fakes_elf1Bmuf1G").a("Fakes_elf1Smuf1S").a("Fakes_elf1Smuf1B").a("Fakes_elf1Bmuf1S")
    .a("Fakes_elf1Gmuf0G").a("Fakes_elf1Gmuf0S").a("Fakes_elf1Smuf0G").a("Fakes_elf1Gmuf0B")
    .a("Fakes_elf1Bmuf0G").a("Fakes_elf1Smuf0S").a("Fakes_elf1Smuf0B").a("Fakes_elf1Bmuf0S")
    .a("Fakes_elf0Gmuf1G").a("Fakes_elf0Gmuf1S").a("Fakes_elf0Smuf1G").a("Fakes_elf0Gmuf1B")
    .a("Fakes_elf0Bmuf1G").a("Fakes_elf0Smuf1S").a("Fakes_elf0Smuf1B").a("Fakes_elf0Bmuf1S");
  
  VS s_Fakes_LF;
  s_Fakes_LF.a("Fakes_elf0G").a("Fakes_muf0G").a("Fakes_elf0S")
    .a("Fakes_muf0S").a("Fakes_elf0B").a("Fakes_muf0B")
    .a("Fakes_elf0Gelf0G").a("Fakes_muf0Gmuf0G").a("Fakes_elf0Gelf0S").a("Fakes_muf0Gmuf0S")
    .a("Fakes_elf0Gelf0B").a("Fakes_muf0Gmuf0B").a("Fakes_elf0Self0S").a("Fakes_muf0Smuf0S")
    .a("Fakes_elf0Self0B").a("Fakes_muf0Smuf0B")
    .a("Fakes_elf0Gmuf0G").a("Fakes_elf0Gmuf0S").a("Fakes_elf0Smuf0G").a("Fakes_elf0Gmuf0B")
    .a("Fakes_elf0Bmuf0G").a("Fakes_elf0Smuf0S").a("Fakes_elf0Smuf0B").a("Fakes_elf0Bmuf0S");
  
  m_Title["QCD"] = "QCD multijets";
  //m_Color["QCD"] = 7023;
  m_Color["QCD"] = 8007;
  m_Strings["QCD"] = VS().a("QCD");
  m_Strings["QCD"] += AddPrefix("QCD", s_Fakes_HF);
  m_Strings["QCD"] += AddPrefix("QCD", s_Fakes_LF);

  m_Title["HF_Fakes"] = "HF leptons";
  //m_Color["HF_Fakes"] = 7022;
  m_Color["HF_Fakes"] = 8008;
  m_Strings["HF_Fakes"] = VS();
  m_Strings["HF_Fakes"] += AddPrefix("ttbar", s_Fakes_HF);
  m_Strings["HF_Fakes"] += AddPrefix("Wjets", s_Fakes_HF);
  m_Strings["HF_Fakes"] += AddPrefix("ZDY", s_Fakes_HF);
  m_Strings["HF_Fakes"] += AddPrefix("DB", s_Fakes_HF);
  m_Strings["HF_Fakes"] += AddPrefix("ST", s_Fakes_HF);
  m_Strings["HF_Fakes"] += AddPrefix("TB", s_Fakes_HF);

  m_Title["LF_Fakes"] = "LF/fake leptons";
  //m_Color["LF_Fakes"] = 7021;
  m_Color["LF_Fakes"] = 8009;
  m_Strings["LF_Fakes"] = VS();
  m_Strings["LF_Fakes"] += AddPrefix("ttbar", s_Fakes_LF);
  m_Strings["LF_Fakes"] += AddPrefix("Wjets", s_Fakes_LF);
  m_Strings["LF_Fakes"] += AddPrefix("ZDY", s_Fakes_LF);
  m_Strings["LF_Fakes"] += AddPrefix("DB", s_Fakes_LF);
  m_Strings["LF_Fakes"] += AddPrefix("ST", s_Fakes_LF);
  m_Strings["LF_Fakes"] += AddPrefix("TB", s_Fakes_LF);
  
  m_Title["Fakes"] = "fake leptons";
  //m_Color["Fakes"] = 7021;
  m_Color["Fakes"] = 8010;
  m_Strings["Fakes"] = VS().a("ttbar_Fakes_elf0").a("ttbar_Fakes_elf1").a("ttbar_Fakes_muf0").a("ttbar_Fakes_muf1").a("Wjets_Fakes_elf0").a("Wjets_Fakes_elf1").a("Wjets_Fakes_muf0").a("Wjets_Fakes_muf1");

//VS().a("Fakes_elf0").a("Fakes_elf1").a("Fakes_elf2").a("Fakes_muf0").a("Fakes_muf1").a("Fakes_muf2");
  

  m_Title["HF"] = "heavy flavor";
  m_Color["HF"] = 7022;
  m_Strings["HF"] = VS().a("Fakes_elf0").a("Fakes_muf0");

  m_Title["LF"] = "light flavor";
  m_Color["LF"] = 7021;
  m_Strings["LF"] = VS().a("Fakes_elf1").a("Fakes_muf1");

  m_Title["ttbar_Fakes"] = "t #bar{t} fakes";
  m_Strings["ttbar_Fakes"] = VS();
  m_Strings["ttbar_Fakes"] += AddPrefix("ttbar",s_Fakes_LF);
  m_Strings["ttbar_Fakes"] += AddPrefix("ttbar",s_Fakes_HF);
  m_Color["ttbar_Fakes"] = 7020;

  m_Title["Wjets_Fakes"] = "W+jets fakes";
  m_Strings["Wjets_Fakes"] = VS();
  m_Strings["Wjets_Fakes"] += AddPrefix("Wjets",s_Fakes_LF);
  m_Strings["Wjets_Fakes"] += AddPrefix("Wjets",s_Fakes_HF);
  m_Color["Wjets_Fakes"] = 7023;

  m_Title["ST_Fakes"] = "single top fakes";
  m_Strings["ST_Fakes"] = VS();
  m_Strings["ST_Fakes"] += AddPrefix("ST",s_Fakes_LF);
  m_Strings["ST_Fakes"] += AddPrefix("ST",s_Fakes_HF);
  m_Color["ST_Fakes"] = 7024;



  m_Title["DB_Fakes"] = "di-boson fakes";
  m_Strings["DB_Fakes"] = VS();
  m_Strings["DB_Fakes"] += AddPrefix("DB",s_Fakes_LF);
  m_Strings["DB_Fakes"] += AddPrefix("DB",s_Fakes_HF);
  m_Color["DB_Fakes"] = 7012;

  m_Title["TB_Fakes"] = "tri-boson fakes";
  m_Strings["TB_Fakes"] = VS();
  m_Strings["TB_Fakes"] += AddPrefix("TB",s_Fakes_LF);
  m_Strings["TB_Fakes"] += AddPrefix("TB",s_Fakes_HF);
  m_Color["TB_Fakes"] = 7013;

  m_Title["ZDY_Fakes"] = "Z / #gamma* + jets fakes";
  m_Strings["ZDY_Fakes"] = VS();
  m_Strings["ZDY_Fakes"] += AddPrefix("ZDY",s_Fakes_LF);
  m_Strings["ZDY_Fakes"] += AddPrefix("ZDY",s_Fakes_HF);
  m_Color["ZDY_Fakes"] = 7014;

  m_Title["ttbar_Fakes_elf0"] = "HF el t #bar{t} + jets";
  m_Title["ttbar_Fakes_elf1"] = "LF+unm. el t #bar{t} + jets";
  m_Title["Wjets_Fakes_elf0"] = "HF el W + jets";
  m_Title["Wjets_Fakes_elf1"] = "LF+unm. el  W + jets";
  m_Title["ttbar_Fakes_muf0"] = "HF mu t #bar{t} + jets";
  m_Title["ttbar_Fakes_muf1"] = "LF+unm. mu t #bar{t} + jets";
  m_Title["Wjets_Fakes_muf0"] = "HF mu W + jets";
  m_Title["Wjets_Fakes_muf1"] = "LF+unm. mu  W + jets";

  m_Title["Total"] = "total background";
  m_Color["Total"] = 7000;
  m_Strings["Total"] = VS().a("ttbar").a("ST").a("DB").a("ZDY").a("Wjets").a("Fakes_elf0").a("Fakes_elf1").a("Fakes_elf2").a("Fakes_muf0").a("Fakes_muf1").a("Fakes_muf2");
  m_Strings["Fakes"] = VS();
  m_Strings["Fakes"] += m_Strings["HF_Fakes"];
  m_Strings["Fakes"] += m_Strings["LF_Fakes"];

  m_Title["ttbar_all"] = "t #bar{t} + jets";
  m_Color["ttbar_all"] = 7011;
  m_Strings["ttbar_all"] = VS().a("ttbar");
  m_Strings["ttbar_all"] += AddPrefix("ttbar", s_Fakes_HF);
  m_Strings["ttbar_all"] += AddPrefix("ttbar", s_Fakes_LF);

  m_Title["ST_all"] = "single top";
  m_Color["ST_all"] = 7010;
  m_Strings["ST_all"] = VS().a("ST");
  m_Strings["ST_all"] += AddPrefix("ST", s_Fakes_HF);
  m_Strings["ST_all"] += AddPrefix("ST", s_Fakes_LF);

  m_Title["DB_all"] = "di-bosons";
  m_Color["DB_all"] = 7051;
  m_Strings["DB_all"] = VS().a("DB");
  m_Strings["DB_all"] += AddPrefix("DB", s_Fakes_HF);
  m_Strings["DB_all"] += AddPrefix("DB", s_Fakes_LF);

  m_Title["TB_all"] = "tri-bosons";
  m_Color["TB_all"] = 7050;
  m_Strings["TB_all"] = VS().a("TB");
  m_Strings["TB_all"] += AddPrefix("TB", s_Fakes_HF);
  m_Strings["TB_all"] += AddPrefix("TB", s_Fakes_LF);

  m_Title["ZDY_all"] = "Z / #gamma* + jets";
  m_Color["ZDY_all"] = 7000;
  m_Strings["ZDY_all"] = VS().a("ZDY");
  m_Strings["ZDY_all"] += AddPrefix("ZDY", s_Fakes_HF);
  m_Strings["ZDY_all"] += AddPrefix("ZDY", s_Fakes_LF);

  m_Title["Wjets_all"] = "W + jets";
  m_Color["Wjets_all"] = 7001;
  m_Strings["Wjets_all"] = VS().a("Wjets");
  m_Strings["Wjets_all"] += AddPrefix("Wjets", s_Fakes_HF);
  m_Strings["Wjets_all"] += AddPrefix("Wjets", s_Fakes_LF);

  // m_Title["Total"] = "total background";
  // m_Color["Total"] = 7000;
  // m_Strings["Total"] = VS().a("ttbar").a("ST").a("DB").a("ZDY").a("Wjets").a("Fakes_elf0").a("Fakes_elf1").
  //   a("Fakes_elf2").a("Fakes_muf0").a("Fakes_muf1").a("Fakes_muf2");
  
  // leptonic categories
  m_Title["1L"] = "#scale[1.2]{single #it{l}}";
  m_Strings["1L"] = VS().a("Ch1L");//VS().a("1L_elm_elG").a("1L_elp_elG").a("1L_elpm_elG").a("1L_mupm_muG").a("1L_mup_muG").a("1L_mum_muG").a("1L_elm_elS").a("1L_elp_elS").a("1L_mup_muS").a("1L_mum_muS").a("1L_elp_elB").a("1L_elm_elB").a("1L_mup_muB").a("1L_mum_muB");
  
  m_Title["1Lgold"] = "#scale[1.2]{single #it{l}}";
  m_Strings["1Lgold"] = VS().a("1L_elm_elG").a("1L_elp_elG").a("1L_elpm_elG").a("1L_mupm_muG").a("1L_mup_muG").a("1L_mum_muG");
  
  m_Title["1Lel"] = "#scale[1.2]{single e}";
  m_Strings["1Lel"] = VS().a("1L_elp_elG").a("1L_elm_elG").a("1L_elpm_elG").a("1L_elp_elS").a("1L_elm_elS").a("1L_elpm_elS").a("1L_elp_elB").a("1L_elm_elB").a("1L_elpm_elB");

  m_Title["1Lmu"] = "#scale[1.2]{single #mu}";
  m_Strings["1Lmu"] = VS().a("1L_mup_muG").a("1L_mum_muG").a("1L_mupm_muG").a("1L_mup_muS").a("1L_mum_muS").a("1L_mupm_muS").a("1L_mup_muB").a("1L_mum_muB").a("1L_mupm_muB");
  
  m_Title["1Lmugold"] = "#scale[1.2]{single gold #mu}";
  m_Strings["1Lmugold"] = VS().a("1L_mup_muG").a("1L_mum_mG");

  m_Title["1Lelp"] = "#scale[1.2]{single e^{+}}";
  m_Strings["1Lelp"] = VS().a("1L_elp_elG");

  m_Title["1Lelm"] = "#scale[1.2]{single e^{-}}";
  m_Strings["1Lelm"] = VS().a("1L_elm_elG");

  m_Title["1Lmup"] = "#scale[1.2]{single #mu^{+}}";
  m_Strings["1Lmup"] = VS().a("1L_mup_muG");

  m_Title["1Lmum"] = "#scale[1.2]{single #mu^{-}}";
  m_Strings["1Lmum"] = VS().a("1L_mum_muG");

  m_Title["1Lp"] = "#scale[1.2]{single #it{l}^{+}}";
  m_Strings["1Lp"] = VS().a("1L_elp_elG").a("1L_mup_muG");

  m_Title["1Lm"] = "#scale[1.2]{single #it{l}^{-}}";
  m_Strings["1Lm"] = VS().a("1L_elm_elG").a("1L_mum_muG");

  m_Title["1Lsilver"] = "#scale[1.2]{single silver #it{l}}";
  m_Strings["1Lsilver"] = VS().a("elpm_slvr").a("mupm_slvr");//a("1L_elp_elS").a("1L_elm_elS").a("1L_mup_muS").a("1L_mum_muS");
  
  m_Title["1Lelsilver"] = "#scale[1.2]{single silver e}";
  m_Strings["1Lelsilver"] = VS().a("1L_elp_elS").a("1L_elm_elS").a("1L_elpm_elS").a("1L_elp-elS").a("1L_elm-elS").a("1L_elpm-elS");
  
  m_Title["1Lmusilver"] = "#scale[1.2]{single silver #mu}";
  m_Strings["1Lmusilver"] = VS().a("1L_mup_muS").a("1L_mum_muS");

  m_Title["1Lbronze"] = "#scale[1.2]{single bronze #it{l}}";
  m_Strings["1Lbronze"] = VS().a("elpm_bron").a("mupm_bron");//a("1L_elp_elB").a("1L_elm_elB").a("1L_mup_muB").a("1L_mum_muB");
  
  m_Title["1Lelbronze"] = "#scale[1.2]{single bronze e}";
  m_Strings["1Lelbronze"] = VS().a("1L_elp_elB").a("1L_elm_elB").a("1L_elpm_elB").a("1L_elp-elB").a("1L_elm-elB").a("1L_elpm-elB");
  
  m_Title["1Lmubronze"] = "#scale[1.2]{single bronze #mu}";
  m_Strings["1Lmubonze"] = VS().a("1L_mup_muB").a("1L_mum_muB");

  m_Title["2LOSSF"] = "#scale[1.2]{e^{#pm} e^{#mp} or #mu^{#pm} #mu^{#mp}}";
  m_Strings["2LOSSF"] = VS().a("2LOS_el^el_elGelG").a("2LOS_mu^mu_muGmuG").a("2LOS_elel^0_elGelG").a("2LOS_mumu^0_muGmuG");
  
  m_Title["2LOSOF"] = "#scale[1.2]{e^{#pm} #mu^{#mp}}";
  m_Strings["2LOSOF"] = VS().a("2LOS_el^mu_elGmuG").a("2LOS_elmu^0_elGmuG");

  m_Title["2LSSSF"] = "#scale[1.2]{e^{#pm} e^{#pm} or #mu^{#pm} #mu^{#pm}}";
  m_Strings["2LSSSF"] = VS().a("2LOS_el^el_elGelG").a("2LOS_mu^mu_muGmuG").a("2LOS_elel^0_elGelG").a("2LOS_mumu^0_muGmuG");
  
  m_Title["2LSSOF"] = "#scale[1.2]{e^{#pm} #mu^{#pm}}";
  m_Strings["2LSSOF"] = VS().a("2LSS_el^mu_elGmuG").a("2LSS_elmu^0_elGmuG");

  m_Title["2LOSSFsilver"] = "#scale[1.2]{e^{#pm} e^{#mp} or #mu^{#pm} #mu^{#mp}, #geq 1 silver #it{l}}";
  m_Strings["2LOSSFsilver"] = VS().a("2LOS_el^el_el0el1").a("2LOS_mu^mu_mu0mu1").a("2LOS_elel^0_el0el1").a("2LOS_mumu^0_mu0mu1")
    .a("2LOS_el^el_el1el1").a("2LOS_mu^mu_mu1mu1").a("2LOS_elel^0_el1el1").a("2LOS_mumu^0_mu1mu1")
    .a("2LOS_el^el_el1el2").a("2LOS_mu^mu_mu1mu2").a("2LOS_elel^0_el1el2").a("2LOS_mumu^0_mu1mu2");
  
  m_Title["2LOSOFsilver"] = "#scale[1.2]{e^{#pm} #mu^{#mp}, #geq 1 silver #it{l}}";
  m_Strings["2LOSOFsilver"] = VS().a("2LOS_el^mu_el0mu1").a("2LOS_elmu^0_el0mu1").a("2LOS_el^mu_mu0el1").a("2LOS_elmu^0_mu0el1")
    .a("2LOS_el^mu_el1mu1").a("2LOS_elmu^0_el1mu1").a("2LOS_el^mu_mu1el2").a("2LOS_elmu^0_mu1el2")
    .a("2LOS_el^mu_el1mu2").a("2LOS_elmu^1_el1mu2");

  m_Title["2LSSSFsilver"] = "#scale[1.2]{e^{#pm} e^{#pm} or #mu^{#pm} #mu^{#pm}, #geq 1 silver #it{l}}";
  m_Strings["2LSSSFsilver"] = VS().a("2LSS_el^el_el0el1").a("2LSS_mu^mu_mu0mu1").a("2LSS_elel^0_el0el1").a("2LSS_mumu^0_mu0mu1")
    .a("2LSS_el^el_el1el1").a("2LSS_mu^mu_mu1mu1").a("2LSS_elel^0_el1el1").a("2LSS_mumu^0_mu1mu1")
    .a("2LSS_el^el_el1el2").a("2LSS_mu^mu_mu1mu2").a("2LSS_elel^0_el1el2").a("2LSS_mumu^0_mu1mu2");
  
  m_Title["2LSSOFsilver"] = "#scale[1.2]{e^{#pm} #mu^{#pm}, #geq 1 silver #it{l}}";
  m_Strings["2LSSOFsilver"] = VS().a("2LSS_el^mu_el0mu1").a("2LSS_elmu^0_el0mu1").a("2LSS_el^mu_mu0el1").a("2LSS_elmu^0_mu0el1")
    .a("2LSS_el^mu_el1mu1").a("2LSS_elmu^0_el1mu1").a("2LSS_el^mu_mu1el2").a("2LSS_elmu^0_mu1el2")
    .a("2LSS_el^mu_el1mu2").a("2LSS_elmu^1_el1mu2");

  m_Title["2LOSSFbronze"] = "#scale[1.2]{e^{#pm} e^{#mp} or #mu^{#pm} #mu^{#mp}, #geq 1 bronze #it{l}}";
  m_Strings["2LOSSFbronze"] = VS().a("2LOS_el^el_el0el2").a("2LOS_mu^mu_mu0mu2").a("2LOS_elel^0_el0el2").a("2LOS_mumu^0_mu0mu2")
    .a("2LOS_el^el_el1el2").a("2LOS_mu^mu_mu1mu2").a("2LOS_elel^0_el1el2").a("2LOS_mumu^0_mu1mu2")
    .a("2LOS_el^el_el2el2").a("2LOS_mu^mu_mu2mu2").a("2LOS_elel^0_el2el2").a("2LOS_mumu^0_mu2mu2");
  
  m_Title["2LOSOFbronze"] = "#scale[1.2]{e^{#pm} #mu^{#mp}, #geq 1 bronze #it{l}}";
  m_Strings["2LOSOFbronze"] = VS().a("2LOS_el^mu_el0mu2").a("2LOS_elmu^0_el0mu2").a("2LOS_el^mu_mu0el2").a("2LOS_elmu^0_mu0el2")
    .a("2LOS_el^mu_el1mu2").a("2LOS_elmu^0_el1mu2").a("2LOS_el^mu_mu1el2").a("2LOS_elmu^0_mu1el2")
    .a("2LOS_el^mu_el2mu2").a("2LOS_elmu^1_el2mu2");

  m_Title["2LSSSFbronze"] = "#scale[1.2]{e^{#pm} e^{#pm} or #mu^{#pm} #mu^{#pm}, #geq 1 bronze #it{l}}";
  m_Strings["2LSSSFbronze"] = VS().a("2LSS_el^el_el0el2").a("2LSS_mu^mu_mu0mu2").a("2LSS_elel^0_el0el2").a("2LSS_mumu^0_mu0mu2")
    .a("2LSS_el^el_el1el2").a("2LSS_mu^mu_mu1mu2").a("2LSS_elel^0_el1el2").a("2LSS_mumu^0_mu1mu2")
    .a("2LSS_el^el_el2el2").a("2LSS_mu^mu_mu2mu2").a("2LSS_elel^0_el2el2").a("2LSS_mumu^0_mu2mu2");
  
  m_Title["2LSSOFbronze"] = "#scale[1.2]{e^{#pm} #mu^{#pm}, #geq 1 bronze #it{l}}";
  m_Strings["2LSSOFbronze"] = VS().a("2LSS_el^mu_el0mu2").a("2LSS_elmu^0_el0mu2").a("2LSS_el^mu_mu0el2").a("2LSS_elmu^0_mu0el2")
    .a("2LSS_el^mu_el1mu2").a("2LSS_elmu^0_el1mu2").a("2LSS_el^mu_mu1el2").a("2LSS_elmu^0_mu1el2")
    .a("2LSS_el^mu_el2mu2").a("2LSS_elmu^1_el2mu2");

  // hadronic categories
  m_Title["0j0svS"] = "#splitline{0 jets}{0 SV-tags} #scale[1.2]{#in S}";

  m_Title["0j1svS"] = "#splitline{0 jets}{1 SV-tag} #scale[1.2]{#in S}";

  m_Title["0jge1svS"] = "#splitline{0 jets}{#geq 1 SV-tag} #scale[1.2]{#in S}";

  m_Title["0jge2svS"] = "#splitline{0 jets}{#geq 2 SV-tags} #scale[1.2]{#in S}";

  m_Title["0jS"] = "#splitline{0 jets}{incl. SV-tags} #scale[1.2]{#in S}";
  m_Strings["0jS"] = VS().a("0j0svS").a("0jge1svS").a("0jS");

  m_Title["1j0bge1svS"] = "#splitline{1 jet, 0 b-tags}{#geq 1 SV-tag} #scale[1.2]{#in S}";
  
  m_Title["1j0b0svS"] = "#splitline{1 jet, 0 b-tags}{0 SV-tag} #scale[1.2]{#in S}";

  m_Title["1jge1svS"] = "#splitline{1 jet, incl. b-tags}{#geq 1 SV-tag} #scale[1.2]{#in S}";

  m_Title["1j0svS"] = "#splitline{1 jet, incl. b-tags}{0 SV-tags} #scale[1.2]{#in S}";

  m_Title["1jS"] = "#splitline{1 jet}{incl. b-tags} #scale[1.2]{#in S}";
  m_Strings["1jS"] = VS().a("1jS").a("1j0bge1svS").a("1j0svS").a("1jge1svS").a("1j0b0svS").a("1j1b0svS").a("1j1bge1svS");

  m_Title["2jS"] = "#splitline{2 jets}{incl. b-tags} #scale[1.2]{#in S}";
  m_Strings["2jS"] = VS().a("2jS").a("2j0bS").a("2j1bS").a("2j2bS");

  m_Title["2j0bS"] = "#splitline{2 jets}{0 b-tags} #scale[1.2]{#in S}";

  m_Title["2j1bS"] = "#splitline{2 jets}{1 b-tags} #scale[1.2]{#in S}";

  m_Title["2j2bS"] = "#splitline{2 jets}{2 b-tags} #scale[1.2]{#in S}";

  m_Title["3j0bS"] = "#splitline{3 jets}{0 b-tags} #scale[1.2]{#in S}";

  m_Title["3j1bS"] = "#splitline{3 jets}{1 b-tags} #scale[1.2]{#in S}";

  m_Title["3jge2bS"] = "#splitline{3 jets}{#geq 2 b-tags} #scale[1.2]{#in S}";

  m_Title["4j0bS"] = "#splitline{4 jets}{0 b-tags} #scale[1.2]{#in S}";

  m_Title["4j1bS"] = "#splitline{4 jets}{1 b-tags} #scale[1.2]{#in S}";

  m_Title["4jge2bS"] = "#splitline{4 jets}{#geq 2 b-tags} #scale[1.2]{#in S}";

  m_Title["5j0bS"] = "#splitline{5 jets}{0 b-tags} #scale[1.2]{#in S}";

  m_Title["5j1bS"] = "#splitline{5 jets}{1 b-tags} #scale[1.2]{#in S}";

  m_Title["5jge2bS"] = "#splitline{5 jets}{#geq 2 b-tags} #scale[1.2]{#in S}";
  m_Title["3jS"] = "#splitline{3 jets}{incl. b-tags} #scale[1.2]{#in S}";
  m_Strings["3jS"] = VS().a("3j0bS").a("3j1bS").a("3ge2bS");
  
  m_Title["4jS"] = "#splitline{4 jets}{incl. b-tags} #scale[1.2]{#in S}";
  m_Strings["4jS"] = VS().a("4j0bS").a("4j1bS").a("4ge2bS");
  
  m_Title["5jS"] = "#splitline{5 jets}{incl. b-tags} #scale[1.2]{#in S}";
  m_Strings["5jS"] = VS().a("ge5j0bS").a("ge5j1bS").a("ge5ge2bS");
  
  m_Title["ge1j0bISR"] = "#splitline{#geq 1 jet}{0 b-tags} #scale[1.2]{#in ISR}";

  m_Title["ge1jge1bISR"] = "#splitline{#geq 1 jet}{#geq 1 b-tags} #scale[1.2]{#in ISR}";

  m_Title["ge1jISR"] = "#splitline{#geq 1 jet}{incl. b-tags} #scale[1.2]{#in ISR}";
  m_Strings["ge1jISR"] = VS().a("ge1jISR").a("ge1j0bISR").a("ge1jge1bISR");


  m_Title["ge1j0bS"] = "#splitline{#geq 1 jet}{0 b-tags} #scale[1.2]{#in S}";

  m_Title["ge1j1bS"] = "#splitline{#geq 1 jet}{1 b-tag} #scale[1.2]{#in S}";

  m_Title["ge2bjetS"] = "#splitline{#geq 1 jet}{#geq 2 b-tags} #scale[1.2]{#in S}";

  m_ColorDefault.clear();
  for(int i = 0; i < 8; i++)
    m_ColorDefault.push_back(7002+i*10);
  for(int i = 0; i < 8; i++)
    m_ColorDefault.push_back(7000+i*10);
  for(int i = 0; i < 8; i++)
    m_ColorDefault.push_back(7004+i*10);

  m_SignalColor.clear();
  m_SignalColor.push_back(7043);
  m_SignalColor.push_back(7071);
  m_SignalColor.push_back(7041);
  m_SignalColor.push_back(7061);
  m_SignalColor.push_back(7040);
  m_SignalColor.push_back(7072);

  m_CMSLabel = "#bf{#it{CMS}} work-in-progress";
}

TCanvas* FitPlotter::PlotCatSignificance(const string& can_name,
					 const VS& proc_bkg,
					 const VS& proc_sig,
					 const CategoryTree& CT,
					 SignificanceType sType){

  RestFrames::SetStyle();

  ProcessList bkgs;
  ProcessList sigs;

  int Nproc_bkg = proc_bkg.size();
  if(Nproc_bkg == 0)
    return nullptr;

  int Nproc_sig = proc_sig.size();
  if(Nproc_sig == 0)
    return nullptr;

  int Nvis = CT.GetNVisible();
  if(Nvis < 1)
    return nullptr;

  vector<const CategoryTree*> CatTrees;
  CT.GetListVisible(CatTrees);
  CategoryList CatList = GetCategories();

  for(int b = 0; b < Nproc_bkg; b++){

    int index = GetProcesses().Find(proc_bkg[b]);
    if(index < 0)
      continue;

    Process pp = GetProcesses()[index];
    bkgs += pp;
  }

  for(int s = 0; s < Nproc_sig; s++){

    int index = GetProcesses().Find(proc_sig[s]);
    if(index < 0)
      continue;

    Process pp = GetProcesses()[index];
    sigs +=pp;
  }

  int Nbkg = bkgs.GetN();
  int Nsig = sigs.GetN();
  int Ncat = 0;

  double nBkgEvts[Nvis]; // vector over groups of categories
  double nSigEvts[Nsig][Nvis];
  //vector<double> X;

  for(int v = 0; v < Nvis; v++){
    CategoryList cat = CatList.Filter(*CatTrees[v]);

    Ncat = cat.GetN();

    //X.push_back(v);

    nBkgEvts[v] = 0.;
    for(int s = 0; s < Nsig; s++)
      nSigEvts[s][v] = 0.;

    for(int c = 0; c < Ncat; c++){
      // bkg loop
      for(int b = 0; b < Nbkg; b++){
        nBkgEvts[v] += Integral(cat[c], bkgs[b]);
      }

      // sig loop
      for(int s = 0; s < Nsig; s++){
        nSigEvts[s][v] += Integral(cat[c], sigs[s]);
      }
    }
  }

  vector<double> sOverB;
  vector<double> maxSB;
  double tempSB = 0.;
  double max = 0.;

  for(int s = 0; s < Nsig; s++){
    tempSB = 0.;
    for(int v = 0; v < Nvis; v++){

      if(sType == kSB && tempSB < nSigEvts[s][v]/nBkgEvts[v])
	  tempSB = nSigEvts[s][v]/nBkgEvts[v];
      
      if(sType == kSrootB && tempSB < nSigEvts[s][v]/TMath::Sqrt(nBkgEvts[v]))
	tempSB = nSigEvts[s][v]/TMath::Sqrt(nBkgEvts[v]);
      
      if(sType == kZbin && tempSB < CalculateZbi(nSigEvts[s][v],nBkgEvts[v]))
	tempSB = CalculateZbi(nSigEvts[s][v],nBkgEvts[v]);
    }
    if(max < tempSB)
      max = tempSB; 

    maxSB.push_back(tempSB);
  }

  //Plotting
  double ratioType = -999.;
  vector<int> clrs = {kBlue+2,kGreen+2,kRed+2,kOrange-3,kMagenta-1,kBlack,kAzure+2,kGray};
  vector<int> markers = {20,21,22,23,29,33,34,43,49};
  vector<TH1D*> fhists;

  for(int s = 0; s < Nsig; s++){
    fhists.push_back(new TH1D(Form("fhist_%d_%s", s, can_name.c_str()),
				  Form("fhist_%d_%s", s, can_name.c_str()),
				  Nvis, 0., Nvis));

    for(int v = 0; v < Nvis; v++){
      //if(nSigEvts[s][v]/nBkgEvts[v]){

      if(sType == kSB)
	ratioType = nSigEvts[s][v]/nBkgEvts[v];

      if(sType == kSrootB)
	ratioType = nSigEvts[s][v]/TMath::Sqrt(nBkgEvts[v]);

      if(sType == kZbin)
        ratioType = CalculateZbi(nSigEvts[s][v],nBkgEvts[v]);

      fhists[s]->SetBinContent(v+1, ratioType);
      fhists[s]->SetBinError(v+1, ratioType);
	//}
    }
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
                             Form("can_%s", can_name.c_str()),
                             1200, 700);
  double hlo = 0.09;
  double hhi = 0.2;
  double hbo = 0.19;
  double hto = 0.07;
  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->Draw();
  can->cd();

  TString yLabel = "";
  if(sType == kSB)
    yLabel = "S/B";
  if(sType == kSrootB)
    yLabel = "S/#sqrt{B}";
  if(sType == kZbin)
    yLabel = "Z_{Bi}";

  fhists[0]->SetMarkerSize(2);
  fhists[0]->SetMarkerColor(m_SignalColor[0]);
  fhists[0]->SetMarkerStyle(markers[0]);
  fhists[0]->Draw("hist p");
  fhists[0]->GetXaxis()->CenterTitle();
  fhists[0]->GetXaxis()->SetTitleFont(42);
  fhists[0]->GetXaxis()->SetTitleSize(0.05);
  fhists[0]->GetXaxis()->SetTitleOffset(1.0);
  fhists[0]->GetXaxis()->SetLabelFont(42);
  fhists[0]->GetXaxis()->SetLabelSize(std::min(2./double(Nvis), 0.04));
  fhists[0]->GetXaxis()->SetTitle("");
  fhists[0]->GetXaxis()->SetTickSize(0.);
  fhists[0]->GetYaxis()->CenterTitle();
  fhists[0]->GetYaxis()->SetTitleFont(42);
  fhists[0]->GetYaxis()->SetTitleSize(0.04);
  fhists[0]->GetYaxis()->SetTitleOffset(0.85);
  fhists[0]->GetYaxis()->SetLabelFont(42);
  fhists[0]->GetYaxis()->SetLabelSize(0.035);
  fhists[0]->GetYaxis()->SetRangeUser(0.,max*1.2);
  fhists[0]->GetYaxis()->SetTitle(yLabel);

  for(int v = 0; v < Nvis; v++){
    fhists[0]->GetXaxis()->SetBinLabel(v+1, CatTrees[v]->GetSpectroscopicLabel().c_str());
  }
  fhists[0]->LabelsOption("v","X");

  for(int s = 1; s < Nsig; s++){
    fhists[s]->SetMarkerSize(2);
    fhists[s]->SetMarkerColor(m_SignalColor[s]);
    fhists[s]->SetMarkerStyle(markers[s]);  
    fhists[s]->Draw("same hist p");
  }

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nproc_sig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  for(int s = 0; s < Nsig; s++){
    if(sigs[s].Type() == kSig)
      leg->AddEntry(fhists[s],TString(GetSignalTitle(proc_sig[s])), "P");
    if(sigs[s].Type() == kBkg)
      leg->AddEntry(fhists[s],TString(m_Title[proc_sig[s].c_str()]), "P");
  }
  leg->Draw("SAME");

  double eps = 0.0015;
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  TLine* line = new TLine();
  line->SetLineWidth(2);
  line->SetLineColor(kBlack);

  //line->DrawLineNDC(hlo, hbo-0.024*lmax, 1-hhi, hbo-0.0235*lmax);
 
  l.SetTextSize(0.025);
  l.SetTextFont(42);
  l.SetTextAlign(23);
  line->SetLineWidth(1);
  double lo = hlo;
  double hi = hlo;
 
 
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, string("Regions "+CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, m_CMSLabel.c_str());
  l.SetTextSize(0.05);

  can->SetGridy();
  can->SetGridx();
  can->Update(); 
  can->SaveAs("plots/"+TString(can_name)+".pdf");
  can->SaveAs("plots/"+TString(can_name)+".gif");
  return can;

}

TCanvas* FitPlotter::PlotRegionSignificance(const string& can_name,
					    const VS& proc_bkg,
					    const VS& proc_sig,
					    const CategoryTree& CT,
					    PlotType pType,
					    SignificanceType sType){

  RestFrames::SetStyle();

  CategoryList CatList = GetCategories();

  vector<const CategoryTree*> CatTrees;
  CT.GetListDeepest(CatTrees);

  int Nbkg = proc_bkg.size();
  int Nsig = proc_sig.size();
  int Nvis = CatTrees.size();
  int Ncat = CatList.GetN();
  
  if(Nsig == 0 || 
     Nbkg == 0 || 
     Nvis == 0 || 
     Ncat == 0)
    return nullptr;

  CategoryList dumcat = CatList.Filter(*CatTrees[0]);
  const FitBin& fitbin = dumcat[0].GetFitBin();

  int Nbin; 
  if(pType == kFull)
    Nbin = fitbin.NBins();
  if(pType == kRISR || pType == kInv)
    Nbin = fitbin.NRBins();

  ProcessList bkgs = FetchProcs(proc_bkg);
  ProcessList sigs = FetchProcs(proc_sig);

  //store background yields
  double nBkgEvts[Nvis][Nbin];

  for(int v = 0; v < Nvis; v++){

    CategoryList cat = CatList.Filter(*CatTrees[v]);

    vector<double> hbkgVec;
    if(pType == kFull)
      hbkgVec = GetAddedHistValues(cat, bkgs);
    if(pType == kRISR || pType == kInv)
      hbkgVec = IntegrateMperp(fitbin, GetAddedHistValues(cat, bkgs));

    for(int b = 0; b < Nbin; b++)
      nBkgEvts[v][b] = hbkgVec[b];
  }

  //store signal yields
  double nSigEvts[Nsig][Nvis][Nbin];

  for(int s = 0; s < Nsig; s++){
    ProcessList list;
    list += sigs[s];

    for(int v = 0; v < Nvis; v++){

      CategoryList cat = CatList.Filter(*CatTrees[v]);

      vector<double> hsigVec;
      if(pType == kFull)
	hsigVec = GetAddedHistValues(cat, list);
      if(pType == kRISR || pType == kInv)
	hsigVec = IntegrateMperp(fitbin, GetAddedHistValues(cat, list));

      for(int b = 0; b < Nbin; b++)
        nSigEvts[s][v][b] = hsigVec[b];
    }
  }
  
  //plot
  double max = 0.;
  double ratioType = -999.;
  vector<TH1D*> fhists;
  vector<int> markers = {20,21,22,23,29,33,34,43,49};

  for(int s = 0; s < Nsig; s++){
    fhists.push_back(new TH1D(Form("fhist_%d_%s", s, can_name.c_str()),
			      Form("fhist_%d_%s", s, can_name.c_str()),
			      Nvis*Nbin, 0., Nvis*Nbin));

    for(int v = 0; v < Nvis; v++){
      for(int b = 0; b < Nbin; b++){
	
	if(sType == kSB)
	  ratioType = nSigEvts[s][v][b]/nBkgEvts[v][b];

	if(sType == kSrootB)
          ratioType = nSigEvts[s][v][b]/TMath::Sqrt(nBkgEvts[v][b]);

	if(sType == kZbin)
          ratioType = CalculateZbi(nSigEvts[s][v][b],nBkgEvts[v][b]);

	int index;
	if(pType == kInv)
	  index = v*Nbin+b+1;
	else
	  index = b*Nvis+v+1;
	
	fhists[s]->SetBinContent(index, ratioType);
	fhists[s]->SetBinError(index, ratioType);
	
	//Get max plot value
	if (max < ratioType)
	  max = ratioType;
      }                                                                          
    }
  }
  
  TString yLabel = "";
  if(sType == kSB)
    yLabel = "S/B";
  if(sType == kSrootB)
    yLabel = "S/#sqrt{B}";
  if(sType == kZbin)
    yLabel = "Z_{Bi}";
    
  for(int b = 0; b < Nbin*Nvis; b++)
    fhists[0]->GetXaxis()->SetBinLabel(b+1, "");
  
  fhists[0]->LabelsOption("v","X");

  int Depth = CT.GetDepth();
  VS rlabels;
  VS vlabels;

  int maxSize = 0;
  for(int r = 0; r < Nbin; r++){
    string label = fitbin[r].GetRBinLabel();
    rlabels += label;
    if(maxSize < label.size())
      maxSize = label.size();
  }
    
  for(int v = 0; v < Nvis; v++)
    vlabels += CatTrees[v]->GetPlainLabel(Depth);

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  TCanvas* can = new TCanvas(Form("can_%s", can_name.c_str()),
                             Form("can_%s", can_name.c_str()),
                             1200, 700);

  double yline = can->GetBottomMargin() - 0.02;
  double hlo = 0.1;
  double hhi = 0.25;
  double hto = 0.07;
  double hbo = 0.19;

  can->SetLeftMargin(hlo);
  can->SetRightMargin(hhi);
  can->SetBottomMargin(hbo);
  can->SetTopMargin(hto);
  can->SetGridx();
  can->SetGridy();
  can->Draw();
  can->cd();

  fhists[0]->SetMarkerSize(2);
  fhists[0]->SetMarkerColor(m_SignalColor[0]);
  fhists[0]->SetMarkerStyle(markers[0]);
  fhists[0]->Draw("hist p");
  fhists[0]->GetXaxis()->CenterTitle();
  fhists[0]->GetXaxis()->SetTitleFont(42);
  fhists[0]->GetXaxis()->SetTitleSize(0.045);
  fhists[0]->GetXaxis()->SetTitleOffset(1.0);
  fhists[0]->GetXaxis()->SetLabelFont(42);
  fhists[0]->GetXaxis()->SetLabelSize(0.05);
  fhists[0]->GetXaxis()->SetTitle("");
  fhists[0]->GetXaxis()->SetTickSize(0.);
  fhists[0]->GetYaxis()->CenterTitle();
  fhists[0]->GetYaxis()->SetTitleFont(42);
  fhists[0]->GetYaxis()->SetTitleSize(0.04);
  fhists[0]->GetYaxis()->SetTitleOffset(0.85);
  fhists[0]->GetYaxis()->SetLabelFont(42);
  fhists[0]->GetYaxis()->SetLabelSize(0.035);
  fhists[0]->GetYaxis()->SetRangeUser(0.,max*1.35);
  fhists[0]->GetYaxis()->SetTitle(yLabel);
  
  for(int s = 1; s < Nsig; s++){
    fhists[s]->SetMarkerSize(2);
    fhists[s]->SetMarkerColor(m_SignalColor[s]);
    fhists[s]->SetMarkerStyle(markers[s]);
    fhists[s]->Draw("same hist p");
  }

  TLegend* leg = new TLegend(1.-hhi+0.007, 1.- (Nsig+1)*(1.-0.49)/9., 0.98, 1.-hto-0.005);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetFillColor(kWhite);
  leg->SetLineColor(kWhite);
  leg->SetShadowColor(kWhite);

  for(int s = 0; s < Nsig; s++){
    if(sigs[s].Type() == kSig)
      leg->AddEntry(fhists[s],TString(GetSignalTitle(proc_sig[s])), "P");
    if(sigs[s].Type() == kBkg)
      leg->AddEntry(fhists[s],TString(m_Title[proc_sig[s].c_str()]), "P");
  }
  leg->Draw("SAME");

  if(pType == kFull)
    DrawMR(fitbin,can,can);
  if(pType == kRISR)
    DrawRM(fitbin,can,can);
  if(pType == kInv)
    DrawRM(fitbin,can,can,vlabels);
    
  double eps = 0.0015;
  
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();

  l.SetTextSize(std::min(0.03, 1.5/double(Nvis*Nbin)));
  l.SetTextFont(42);
  l.SetTextAlign(32);
  l.SetTextAngle(90);

  for(int b = 0; b < Nvis*Nbin; b++){
    string label;
    double xpos, ypos;
    if(pType == kInv){
      label = rlabels[b%Nbin];
      xpos = hlo+(1.-hhi-hlo)/double(Nvis*Nbin)*(0.5+b);
      ypos = yline+0.015;
    }
    else{
      label = CatTrees[b%Nvis]->GetPlainLabel(Depth);
      xpos = hlo+(1.-hhi-hlo)/double(Nvis*Nbin)*(0.5+b);
      ypos = 1 - hto - 4*eps;
      l.SetTextColor(7004 + 10*((b%Nvis)%8));
    }
	
    l.DrawLatex(xpos,ypos,label.c_str());
  }

  l.SetTextAngle(0);
  l.SetTextColor(kBlack);
  l.SetTextAlign(31);
  l.SetTextSize(0.04);
  l.DrawLatex(1.-hhi-eps*4, 1.-hto+0.02, string("Regions "+CT.GetSpectroscopicLabel()).c_str());
  l.SetTextAlign(11);
  l.DrawLatex(hlo+eps*4, 1.-hto+0.02, m_CMSLabel.c_str());

  can->Update();
  can->SaveAs("plots/"+TString(can_name)+".pdf");
  can->SaveAs("plots/"+TString(can_name)+".gif");
  return can;
}

void FitPlotter::FindBkgZeros(const VS& proc_bkg){

  RestFrames::SetStyle();

  double thresh = 1e-5;

  ProcessList bkgs;
  int Nbkg = proc_bkg.size();

  CategoryList CatList = GetCategories();

  for(int b = 0; b < Nbkg; b++){

    VS vproc;

    if(m_Strings.count(proc_bkg[b]) != 0)
      vproc = m_Strings[proc_bkg[b]];
    else
      vproc += proc_bkg[b];

    for(int v = 0; v < vproc.size(); v++){
      
      int index = GetProcesses().Find(vproc[v]);
      if(index < 0)
	continue;
      
      double itot = 0.;
    
      Process pp = GetProcesses()[index];
      bkgs += pp;
    }
  }

  int count = 0;
  int totCat = 0;
  string tempStr = "";
  vector<double> bkgVec;
  VS catZeros;

  for(int c = 0; c < CatList.GetN(); c++){
    Category cat = CatList[c];
    CategoryList catlist;
    catlist += cat;
    const FitBin& fitbin = CatList[c].GetFitBin();
    
    bkgVec = GetAddedHistValues(catlist, bkgs);

    int ibin = 0;
    int NR = fitbin.NRBins();

    for(int r = 0; r < NR; r++){
      int NM =  fitbin.NMBins(r);
      //VS mlabels;
      VS mlabels = fitbin[r].GetMBinLabels();
      //VS rlabels;

      string rlabels = fitbin[r].GetRBinLabel();

      for(int m = 0; m < NM; m++){

	double numBkg = bkgVec[ibin];
	if(numBkg < thresh){

	  count++;
          cout << endl;
          cout << "zero number " << count << ": " << endl;
	  cout << "nBkg: " << numBkg << endl;
	  cout << "bkg bin: " << ibin << endl;
	  cout << "Mperp bin: " << mlabels[m] << endl;
          cout << "RISR bin: " << rlabels << endl;
          cout << cat.FullLabel() << endl;
	  /*
	  for(int b = 0; b < bkgs.GetN(); b++){

	    //cout << endl;
	    cout << "  Processing " << bkgs[b].Name() << " background: " << endl;

	    ProcessList bkgsZero;
	    vector<double> tempZero;

	    bkgsZero += bkgs[b];
	    //cout << "here" << endl;
	    tempZero = GetAddedHistValues(catlist, bkgsZero);
	    //cout << tempZero.size() << endl;
	    for(int i = 0; i < tempZero.size(); i++){
	      cout << "   Background bin " << i << ": " << endl;
	      cout << "    Value: " << tempZero[i] << endl;
	    }
	    tempZero.clear();

	  }
	  */
	  if(tempStr != cat.FullLabel()){
	    totCat++;
	    catZeros += cat.FullLabel();
	  }
	  tempStr = cat.FullLabel();
	}
	ibin++;
      }
    }
    bkgVec.clear();
  }
  cout << "Total zeros: " << count << endl;
  cout << "Total categories: " << totCat << endl;

  for(int i = 0; i < catZeros.size(); i++){
    cout << catZeros[i] << endl;
  }
}

void FitPlotter::FindBkgRare(){

  VS proc_bkg;
  proc_bkg.a("Non-Rare").a("ST_all").a("TB_all");

  int Nbkg = proc_bkg.size();

  CategoryList CatList = GetCategories();

  vector<double> bkgVec;

  VS rareOnlyCat;
  string tmpStr = "";

  for(int c = 0; c < CatList.GetN(); c++){

    bool foundZero = false;
    vector<int> binVec;

    Category cat = CatList[c];

    CategoryList catlist;
    catlist += cat;
    const FitBin& fitbin = CatList[c].GetFitBin();

    for(int b = 0; b < Nbkg; b++){

      ProcessList bkgs;
      VS vproc;

      if(m_Strings.count(proc_bkg[b]) != 0)
	vproc = m_Strings[proc_bkg[b]];
      else
	vproc += proc_bkg[b];

      for(int v = 0; v < vproc.size(); v++){

	int index = GetProcesses().Find(vproc[v]);
	if(index < 0)
	  continue;

	double itot = 0.;

	Process pp = GetProcesses()[index];
	bkgs += pp;
      }

      bkgVec = GetAddedHistValues(catlist, bkgs); 

      for(int i = 0; i < bkgVec.size(); i++){
	if(proc_bkg[b] == "Non-Rare" && bkgVec[i] < 1e-5){
	  foundZero = true;
	  binVec.push_back(i);
	}
      }

      if(foundZero){
	for(int j = 0; j < binVec.size(); j++){
	  if((proc_bkg[b] == "ST_all" || proc_bkg[b] == "TB_all") && bkgVec[binVec[j]] > 1e-5 && bkgVec.size()>1){
	    cout << endl;
	    cout << cat.FullLabel() << endl;
	    cout << "  " << proc_bkg[b] << ": " << endl;
	    cout << "      bin " << binVec[j] <<endl;           
	    cout << "       value: " << bkgVec[binVec[j]] << endl;
	    if(tmpStr != cat.FullLabel()){
	      rareOnlyCat +=  cat.FullLabel();
	      tmpStr = cat.FullLabel();
	    }
	  }
	}
      }    
      
    }
    foundZero = false;
  }

  cout << "Total categories with only rare backgrounds: " << rareOnlyCat.size() << endl;
  for(int i = 0; i < rareOnlyCat.size(); i++)
    cout << rareOnlyCat[i] << endl; 

}

void FitPlotter::YieldPerBkg(const VS& proc_bkg){

  int Nbkg = proc_bkg.size();

  CategoryList CatList = GetCategories();

  vector<double> bkgVec;

  VS rareOnlyCat;
  string tmpStr = "";

  for(int c = 0; c < CatList.GetN(); c++){

    bool foundZero = false;
    vector<int> binVec;

    Category cat = CatList[c];

    CategoryList catlist;
    catlist += cat;
    const FitBin& fitbin = CatList[c].GetFitBin();

    cout << endl;
    cout << cat.FullLabel() << endl;

    for(int b = 0; b < Nbkg; b++){

      ProcessList bkgs;
      VS vproc;

      if(m_Strings.count(proc_bkg[b]) != 0)
        vproc = m_Strings[proc_bkg[b]];
      else
        vproc += proc_bkg[b];

      for(int v = 0; v < vproc.size(); v++){

        int index = GetProcesses().Find(vproc[v]);
        if(index < 0)
          continue;

        double itot = 0.;

        Process pp = GetProcesses()[index];
        bkgs += pp;
      }

      bkgVec = GetAddedHistValues(catlist, bkgs);

      cout << "  " << proc_bkg[b] << ": " << endl;
      for(int i = 0; i < bkgVec.size(); i++){
	cout << "       value: " << bkgVec[i] << endl;
      }
    
    }
    cout << endl;
    foundZero = false;
  }
}

void FitPlotter::zeroBkgTest(){

  RestFrames::SetStyle();

  VS proc_bkg;
  proc_bkg.a("ttbar").a("ST").a("DB").a("ZDY").a("TB").a("QCD").a("Wjets").a("HF_Fakes").a("LF_Fakes");

  //ProcessList bkgs;
  int Nbkg = proc_bkg.size();

  CategoryList CatList = GetCategories();

  Category cat = CatList[0];
  CategoryList catlist;
  catlist += cat;

  vector<double> bkgVec;
  VS catZeros;


  for(int b = 0; b < Nbkg; b++){

    cout << endl;
    cout << "Processing " << proc_bkg[b] << " background: " << endl;

    ProcessList bkgs;

    int index = GetProcesses().Find(proc_bkg[b]);
    if(index < 0)
      continue;

    double itot = 0.;

    Process pp = GetProcesses()[index];
    bkgs += pp;

    bkgVec = GetAddedHistValues(catlist, bkgs);

    for(int i = 0; i < bkgVec.size(); i++){
      cout << " Background bin " << i << ": " << endl;
      cout << " Value: " << bkgVec[i] << endl;
    }
    bkgVec.clear();
  }
}

ProcessList FitPlotter::FetchProcs(VS proc_list){

  int Nproc = proc_list.size();
  ProcessList procs;

  for(int i = 0; i < Nproc; i++){
    VS vproc;
    if(m_Strings.count(proc_list[i]) != 0)
      vproc = m_Strings[proc_list[i]];
    else
      vproc += proc_list[i];

    for(int p = 0; p < int(vproc.size()); p++){

      int index = GetProcesses().Find(vproc[p]);
      if(index < 0)
	continue;

      Process pp = GetProcesses()[index];
      procs += pp;
    }
  }
  return procs;
}

VS FitPlotter::FetchLabels(VS categories){

  int Nentries = categories.size();
  VS labels;

  for(int i = 0; i < Nentries; i++){
    if(m_Title.count(categories[i]) != 0)
      labels += m_Title[categories[i]];
    else
      labels += categories[i];
  }

  return labels;
}

VS FitPlotter::FetchCats(VS categories){

  int Nentries = categories.size();
  VS cats;

  for(int i = 0; i < Nentries; i++){
    if(m_Strings.count(categories[i]) != 0){

      int N = m_Strings[categories[i]].size();      

      for(int j = 0; j < N; j++)
	cats += m_Strings[categories[i]][j];
    } 
    else
      cats += categories[i];
  }

  return cats;
}

VS FitPlotter::AddPrefix(const string& pre, const VS& post) const {
  VS ret;

  for(auto s : post)
    ret += pre+"_"+s;

  return ret;
}

int FitPlotter::getColor(const string& str){
  return m_Color[str];
}

string FitPlotter::getTitle(const string& str){
  return m_Title[str];
}
