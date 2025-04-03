#include <iostream>
#include <map>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TSystem.h>
#include <TMultiGraph.h>
#include <TEfficiency.h>

#include "FitReader.hh"
#include "shapeComparison.hh"
#include "plotShapeComparison.hh"
#include "RestFrames/RestFrames.hh"

///////////////////////////////////////////
////////// FitReader class
///////////////////////////////////////////

FitReader::FitReader(const string& inputfile,
		     const string& otherfile,
		     const string& otherfold)
  : m_File(inputfile.c_str(), "READ"),
    cms1(8001, 0.24705882352941178, 0.5647058823529412, 0.8549019607843137),  //blue
    cms2(8002, 1., 0.6627450980392157, 0.054901960784313725),                 // gold
    cms3(8003, 0.7411764705882353, 0.12156862745098039, 0.00392156862745098), //red
    cms4(8004, 0.5803921568627451, 0.6431372549019608, 0.6352941176470588),   // grey
    cms5(8005, 0.5137254901960784, 0.17647058823529413, 0.7137254901960784),  // purple
    cms6(8006, 0.6627450980392157, 0.4196078431372549, 0.34901960784313724),  // brown
    cms7(8007, 0.9058823529411765, 0.38823529411764707, 0.),                  // orange
    cms8(8008, 0.7254901960784313, 0.6745098039215687, 0.4392156862745098),   // tan
    cms9(8009, 0.44313725490196076, 0.4588235294117647, 0.5058823529411764),  // charcoal
    cms10(8010, 0.5725490196078431, 0.8549019607843137, 0.8666666666666667) { // light blue

  
  if(otherfile != ""){
    m_FilePtr = new TFile(otherfile.c_str(), "READ");
    if(!m_FilePtr || !m_FilePtr->IsOpen())
      m_FilePtr = nullptr;
    m_FileFold = otherfold+"/";
  } else {
    m_FilePtr = nullptr;
  }
  ReadProcesses();
  
  ReadCategories();
}

FitReader::~FitReader(){
  m_File.Close();
}

void FitReader::ReadProcesses(){
  if(!m_File.IsOpen())
    return;
  
  TTree* tree = (TTree*)m_File.Get("Process");
  if(!tree)
    return;

  m_ProcBranch.InitGet(tree);

  ProcessList ProcSys;

  int N = tree->GetEntries();
  for(int i = 0; i < N; i++){
    tree->GetEntry(i);
    Process p = m_ProcBranch.GetProcess();

  //  std::cout<<"assessing subproc "<<p.Name()<<" ";


    //std::cout << p.Name() << endl;

    if((p.Name().find("Up") != std::string::npos) ||
       (p.Name().find("Down") != std::string::npos)){
      ProcSys += p;
    //  std::cout<<"+ procsys\n";
    }
    else
      m_Proc += p;
  }
  if(m_FilePtr) m_Proc += Process("total_background",kBkg); 
  
  delete tree;
  int Nproc = m_Proc.GetN();
  int Nsys  = ProcSys.GetN();
  for(int p = 0; p < Nproc; p++){
    Systematics sys;
    string proc = m_Proc[p].Name();
    for(int s = 0; s < Nsys; s++){
      string label = ProcSys[s].Name();
      if((proc.find("Fakes") == std::string::npos) &&
	 (label.find("Fakes") != std::string::npos))
	continue;
      if(label.find(proc+"_") == std::string::npos) continue;
      if(label.find(proc) == 0){
	label.replace(0,proc.length()+1,"");
	if(label.find("Up") != std::string::npos)
	  label.replace(label.find("Up"),2,"");
	if(label.find("Down") != std::string::npos)
	  label.replace(label.find("Down"),4,"");
//	std::cout<<"in Nsys loop, adding systematic with label: "<<label<<"\n";
	sys += Systematic(label);
      }
    }
    if(sys.GetN() > 0)
      m_ProcSys[m_Proc[p]] = sys;
    m_Sys += sys;
  }
  //check sig procs sys
  /*if(m_FilePtr) m_ProcSys = Process("total_signal",kSig);
  m_ProcSys.Remove("METUncer_GenMET");
  Nproc = m_ProcSys.GetN();
  for(int p=0; p< Nproc; p++){
	string proc = m_ProcSys	
  }*/
  //hack to check signal systematics only
/*  std::vector<std::string> sysLabels{"JESUncer_Total", "JERUncer_Total",  "METUncer_UnClust", "METUncer_GenMET" };
  Systematics sysSig;
  for(int i=0; i<sysLabels.size(); i++){
    sysSig+= Systematic( sysLabels[i] );
  }
  m_Sys += sysSig;
*/
}

void FitReader::ReadCategories(){
  if(!m_File.IsOpen())
    return;
  
  TTree* tree = (TTree*)m_File.Get("Category");
  if(!tree)
    return;
  
  m_CatBranch.InitGet(tree);

  int N = tree->GetEntries();
  for(int i = 0; i < N; i++){
    tree->GetEntry(i);
    Category cat = m_CatBranch.GetCategory();
    if(m_CatLabel.count(cat.GetLabel()) == 0)
      m_CatLabel[cat.GetLabel()] = true;
    else
      continue;
    m_Cat += cat;
    if(m_Chan.count(cat.Label()) == 0)
      m_Chan[cat.Label()] = CategoryList();
    m_Chan[cat.Label()] += cat;
  }

  delete tree;
}

bool FitReader::HasSystematic(const Process& proc, const Systematic& sys) const {
  if(m_ProcSys.count(proc) == 0)
    return false;

  return m_ProcSys[proc].Contains(sys);
}

double FitReader::Integral(const Category&   cat,
			   const Process&    proc,
			   const Systematic& sys) const {
//cout << "Integral" << endl;
  const TH1D* hist = GetHistogram(cat, proc, sys);
  if(!hist)
    return 0.;
  //cout << "Integral end " << hist->Integral() << endl;
  return hist->Integral();
}


TH1D* FitReader::GetIntegralHist(const string& name,
				 const CategoryList& cats,
				 const ProcessList&  procs,
				 const Systematic&   sys) const {
  int Np = procs.GetN();
  int Nc = cats.GetN();

  TH1D* hist = nullptr;
  TH1D* histd = nullptr;

  for(int p = 0; p < Np; p++){
    for(int c = 0; c < Nc; c++){
      if(!IsFilled(cats[c], procs[p], sys))
	continue;
    
      if(!hist){
	histd = (TH1D*) GetHistogram(cats[c], procs[p], sys)->Clone("dum");
	histd->Rebin(histd->GetNbinsX());
	hist = (TH1D*) new TH1D(name.c_str(),name.c_str(), 1, 0., 1.);
	hist->SetBinContent(1, histd->GetBinContent(1));
	hist->SetBinError(1, histd->GetBinError(1));
	delete histd;
      } else {
	TH1D* dum = (TH1D*) GetHistogram(cats[c], procs[p], sys)->Clone("dum");
	dum->Rebin(dum->GetNbinsX());
	hist->SetBinContent(1, hist->GetBinContent(1)+dum->GetBinContent(1));
	hist->SetBinError(1, sqrt(hist->GetBinError(1)*hist->GetBinError(1)+dum->GetBinError(1)*dum->GetBinError(1)));
	delete dum;
      }
    }
  }
   
  return hist;
}

TH1D* FitReader::GetAddedHist(const string&       name,
			      const CategoryList& cats,
			      const ProcessList&  procs,
			      const Systematic&   sys) const {
  int Np = procs.GetN();
  int Nc = cats.GetN();
  TH1D* hist = nullptr;

  for(int p = 0; p < Np; p++){
    for(int c = 0; c < Nc; c++){
      if(!IsFilled(cats[c], procs[p], sys))
	continue;
      //if(procs[p].Name() != "total_background" && name.find("plothist_0_") != std::string::npos && p == 0) 
      //if(procs[p].Name() == "QCD" && name.find("plothist_0_") != std::string::npos && p == 0)
      //cout << cats[c].Label()+"_"+cats[c].GetLabel() << " " << procs[p].Name() << endl;  
  if(!hist){
	hist = (TH1D*) GetHistogram(cats[c], procs[p], sys)->Clone(name.c_str());
      } else {
	hist->Add(GetHistogram(cats[c], procs[p], sys));
      }
    }
  }
  //if(hist)
  //cout << "Added hist bins: " << hist->GetNbinsX() << endl;
   
  return hist;
}

vector<double> FitReader::GetAddedHistValues(const CategoryList& cats,
					     const ProcessList&  procs,
					     const Systematic&   sys) const {

  const FitBin& fitbin = cats[0].GetFitBin();
  int NBins = fitbin.NBins();

  vector<double> histValues;

  int Np = procs.GetN();
  int Nc = cats.GetN();

  TH1D* hist = nullptr;
  TString name = "dummy";

  for(int p = 0; p < Np; p++){
    for(int c = 0; c < Nc; c++){
      if(!IsThere(cats[c], procs[p], sys))
        continue;

      if(!hist){
        hist = (TH1D*) GetHistogram(cats[c], procs[p], sys)->Clone(name);
      } else {
        hist->Add(GetHistogram(cats[c], procs[p], sys));
      }
    }
  }

  if(hist){
    int bins = hist->GetNbinsX();
    
    for(int b = 0; b < bins; b++){
      histValues.push_back(hist->GetBinContent(b+1));
    } 

    delete hist;
  }
  else{
    for (int i = 0; i < NBins; i++)
      histValues.push_back(0.);
  }

  return histValues;
}

const TH1D* FitReader::GetHistogram(const Category&   cat,
				    const Process&    proc,
				    const Systematic& sys) const {
  if(!IsFilled(cat, proc, sys)){
    return nullptr;
  }
  if(!sys) {
    //cout << proc.Name() << endl;
    return m_ProcHist[proc][cat];
  } else {
    //cout << "GetHistogram end" << endl;
    return (sys.IsUp() ? m_ProcHistSys[proc][sys][cat].first :
	    m_ProcHistSys[proc][sys][cat].second);
  }
}

const TH2D* FitReader::GetHistogram2D(const Category&   cat,
				      const Process&    proc,
				      const Systematic& sys) const {
  //cout << "Getting histogram!" << endl;
  if(!IsFilled2D(cat, proc, sys)){
    //cout << "hist not filled!" << endl;
    return nullptr;
  }
  //cout << cat.GetLabel() << " " << proc.Name() << " hist integral: " << m_ProcHist_2D[proc][cat]->Integral() << endl;   
  if(!sys){
    //cout << "no sys, hist filled" << endl;
    return m_ProcHist_2D[proc][cat];
  } else {
    cout << "sys" << endl;
    return (sys.IsUp() ? m_ProcHistSys_2D[proc][sys][cat].first :
	    m_ProcHistSys_2D[proc][sys][cat].second);
  }
}

vector<double> FitReader::IntegrateMperp(const FitBin fitBin, const vector<double> histVec) const{
  int NR = fitBin.NRBins();
  int index = 0;
  vector<double> results;

  for(int r = 0; r < NR; r++){
    int NM = fitBin.NMBins(r);
    double sum = 0;
    for(int m = 0; m < NM; m++){
      sum += histVec[index+m];
    }
    index += NM;
    results.push_back(sum);
  }

  return results;
}

TH1D* FitReader::IntegrateMperp(const TString& name, const FitBin fitBin, TH1D* hist) const{

  if(!hist){
    //cout << "no histogram" << endl;
    return nullptr;

  }

  int NR = fitBin.NRBins();
  int index = 0;
  int bins = hist->GetNbinsX();

  TH1D* results = new TH1D(name, name, NR,0.,1.);

  for(int r = 0; r < NR; r++){
    int NM = fitBin.NMBins(r);
    double sum = 0;
    for(int m = 0; m < NM; m++){
      sum += hist->GetBinContent(index+m+1);
    }
    index += NM;
    results->SetBinContent(r+1,sum);
  }

  delete hist;
  return results;
}

double FitReader::CalculateZbi(double Nsig, double Nbkg, double deltaNbkg){
  double Nobs = Nsig+Nbkg;
  double tau = 1./Nbkg/(deltaNbkg*deltaNbkg);
  double aux = Nbkg*tau;
  double Pvalue = TMath::BetaIncomplete(1./(1.+tau),Nobs,aux+1.);
  double sigma = sqrt(2.)*TMath::ErfcInverse(Pvalue*2);

  return (isnan(sigma))?0:sigma;

}

bool FitReader::IsFilled(const Category&   cat,
			 const Process&    proc,
			 const Systematic& sys) const {
  //cout << "IsFilled" << endl;
  if(!sys){
    if(m_ProcHist.count(proc) == 0){
      m_ProcHist[proc] = map<Category,TH1D*>();
    }
    if(m_ProcHist[proc].count(cat) == 0){
      string shist = cat.Label()+"_"+cat.GetLabel()+"/"+proc.Name();
      if(!m_FilePtr){ //if there is no file pointer, just take histogram from original file
	m_ProcHist[proc][cat] = (TH1D*) m_File.Get(shist.c_str());
      } else { 
	if(proc.Type() == kData){ //if there is another file, and the histogram is data, get the info from the TGraphAsymmErrors
	  TGraphAsymmErrors* gr = (TGraphAsymmErrors*)m_FilePtr->Get((m_FileFold+cat.Label()+"_"+cat.GetLabel()+"/data").c_str());
	  if(gr == nullptr){ 
	    m_ProcHist[proc][cat] = (TH1D*) m_FilePtr->Get((m_FileFold+shist).c_str());
	  } else { //else if TGraphAsymmErrors is found
	    double x, y;
	    m_ProcHist[proc][cat] = new TH1D((m_FileFold+cat.Label()+"_"+cat.GetLabel()+"/data").c_str(),"data",(int)gr->GetN(),0.,(double)gr->GetN());
	    for(int i = 0; i < gr->GetN(); i++){
	      gr->GetPoint(i,x,y);
	      m_ProcHist[proc][cat]->SetBinContent(i+1,y);
	      m_ProcHist[proc][cat]->SetBinError(i+1,gr->GetErrorY(i));
	    }
	  }
	} else if (proc.Type() == kSig){
	  m_ProcHist[proc][cat] = (TH1D*) m_File.Get(shist.c_str());
	} else { //else if there is another file but it's not data, just get the histogram	  
	  m_ProcHist[proc][cat] = (TH1D*) m_FilePtr->Get((m_FileFold+shist).c_str());
	}
      }
    }
    return m_ProcHist[proc][cat];
    
  //else if there is a systematic in the function args
  } else {
    //cout << "do IsFilled w sys" << endl;
    if(m_ProcHistSys.count(proc) == 0)
      m_ProcHistSys[proc] = map<Systematic,map<Category,pair<TH1D*,TH1D*> > >();
    //cout << "proc count: " << m_ProcHistSys.count(proc) << endl;   
    if(m_ProcHistSys[proc].count(sys) == 0)
      m_ProcHistSys[proc][sys] = map<Category,pair<TH1D*,TH1D*> >();
    //cout << "sys count: " << m_ProcHistSys[proc].count(sys) << endl;   
    //cout << "cat count: " << m_ProcHistSys[proc][sys].count(cat) << endl;   
    if(m_ProcHistSys[proc][sys].count(cat) == 0){
      m_ProcHistSys[proc][sys][cat] = pair<TH1D*,TH1D*>(nullptr,nullptr);
      
      string label = cat.Label()+"_"+cat.GetLabel();
      string shistUp   = label+"/"+proc.Name()+"_"+sys.Label()+"Up";
      string shistDown = label+"/"+proc.Name()+"_"+sys.Label()+"Down";
      m_ProcHistSys[proc][sys][cat].first  = (TH1D*) m_File.Get(shistUp.c_str());
      m_ProcHistSys[proc][sys][cat].second = (TH1D*) m_File.Get(shistDown.c_str());
    }

//   cout << "up hist: " << m_ProcHistSys[proc][sys][cat].first->GetName() << endl;  
//   cout << "down hist: " << m_ProcHistSys[proc][sys][cat].second->GetName() << endl;  
//   cout << "IsFilled end" << endl; 
    return (sys.IsUp() ? m_ProcHistSys[proc][sys][cat].first :
	    m_ProcHistSys[proc][sys][cat].second);
  }    
}

bool FitReader::IsThere(const Category&   cat,
			const Process&    proc,
			const Systematic& sys) const {
  TH1D* hist = nullptr;
  TH1D* hist2 = nullptr;
  if(!sys){
   // if(m_ProcHist.count(proc) == 0)
     // m_ProcHist[proc] = map<Category,TH1D*>();
   // if(m_ProcHist[proc].count(cat) == 0){
      string shist = cat.Label()+"_"+cat.GetLabel()+"/"+proc.Name();
      if(proc.Type() == kData || !m_FilePtr)
	hist = (TH1D*) m_File.Get(shist.c_str());
      else
	hist = (TH1D*) m_FilePtr->Get((m_FileFold+shist).c_str());

      if(hist){
	delete hist;
	return true;
      } else
	return false;
 //   }
    return m_ProcHist[proc][cat];
    
  } else {
   // if(m_ProcHistSys.count(proc) == 0)
   //   m_ProcHistSys[proc] = map<Systematic,map<Category,pair<TH1D*,TH1D*> > >();
   // if(m_ProcHistSys[proc].count(sys) == 0)
   //   m_ProcHistSys[proc][sys] = map<Category,pair<TH1D*,TH1D*> >();
   // if(m_ProcHistSys[proc][sys].count(cat) == 0){
   //   m_ProcHistSys[proc][sys][cat] = pair<TH1D*,TH1D*>(nullptr,nullptr);
       
      string label = cat.Label()+"_"+cat.GetLabel();
      string shistUp   = label+"/"+proc.Name()+"_"+sys.Label()+"Up";
      string shistDown = label+"/"+proc.Name()+"_"+sys.Label()+"Down";
  //    std::cout<<"isTHere, trying to get labels: "<<label<<" "<<shistUp<<" "<<shistDown<<"\n"; 
//      hist  = m_ProcHistSys[proc][sys][cat].first  = (TH1D*) m_File.Get(shistUp.c_str());
//      hist2 = m_ProcHistSys[proc][sys][cat].second = (TH1D*) m_File.Get(shistDown.c_str());
      hist  = (TH1D*) m_File.Get(shistUp.c_str());
      hist2 = (TH1D*) m_File.Get(shistDown.c_str());

      bool there = true;
      if(hist)
	delete hist;
      else
	there = false;
      if(hist2)
	delete hist2;
      else
	there = false;
      return there;
    }
     
    return (sys.IsUp() ? m_ProcHistSys[proc][sys][cat].first :
	    m_ProcHistSys[proc][sys][cat].second);
  }    


bool FitReader::IsFilled2D(const Category&   cat,
			   const Process&    proc,
			   const Systematic& sys) const {
  
  if(!sys){
    if(m_ProcHist_2D.count(proc) == 0)
      m_ProcHist_2D[proc] = map<Category,TH2D*>();
    if(m_ProcHist_2D[proc].count(cat) == 0){
      string shist = cat.Label()+"_"+cat.GetLabel()+"/"+proc.Name()+"_2D";
      if(proc.Type() == kData || !m_FilePtr)
	m_ProcHist_2D[proc][cat] = (TH2D*) m_File.Get(shist.c_str());
      else
	m_ProcHist_2D[proc][cat] = (TH2D*) m_FilePtr->Get((m_FileFold+shist).c_str());
    }
    
    return m_ProcHist_2D[proc][cat];
    
  } else {
    if(m_ProcHistSys_2D.count(proc) == 0)
      m_ProcHistSys_2D[proc] = map<Systematic,map<Category,pair<TH2D*,TH2D*> > >();
    if(m_ProcHistSys_2D[proc].count(sys) == 0)
      m_ProcHistSys_2D[proc][sys] = map<Category,pair<TH2D*,TH2D*> >();
    if(m_ProcHistSys_2D[proc][sys].count(cat) == 0){
      m_ProcHistSys_2D[proc][sys][cat] = pair<TH2D*,TH2D*>(nullptr,nullptr);
       
      string label = cat.Label()+"_"+cat.GetLabel();
      string shistUp   = label+"/"+proc.Name()+"_"+sys.Label()+"Up_2D";
      string shistDown = label+"/"+proc.Name()+"_"+sys.Label()+"Down_2D";
       
      m_ProcHistSys_2D[proc][sys][cat].first  = (TH2D*) m_File.Get(shistUp.c_str());
      m_ProcHistSys_2D[proc][sys][cat].second = (TH2D*) m_File.Get(shistDown.c_str());
    }
     
    return (sys.IsUp() ? m_ProcHistSys[proc][sys][cat].first :
	    m_ProcHistSys[proc][sys][cat].second);
  }    
}

void FitReader::PrintCategories(bool verbose){
  cout << "*** Fit Categories ***" << endl;
  int N = m_Cat.GetN();
  for(int i = 0; i < N; i++)
    cout << m_Cat[i].Label()+"_"+m_Cat[i].GetLabel() << endl;
  cout << endl;
}

void FitReader::PrintProcesses(bool verbose){
  cout << "*** Fit Processes ***" << endl;
  int N = m_Proc.GetN();
  for(int i = 0; i < N; i++){
    cout << m_Proc[i].Name() << endl;
    if(m_ProcSys.count(m_Proc[i]) == 0 || !verbose)
      continue;
    cout << "  Sys:  " << m_ProcSys[m_Proc[i]][0].Label() << endl;
    int Nsys = m_ProcSys[m_Proc[i]].GetN();
    for(int s = 1; s < Nsys; s++)
      cout << "        " << m_ProcSys[m_Proc[i]][s].Label() << endl;
  }
  cout << endl;
}

const ProcessList&  FitReader::GetProcesses() const {
  return m_Proc;
}

const CategoryList& FitReader::GetCategories(const string& channel) const {
  if(m_Chan.count(channel) > 0)
    return m_Chan[channel];
  
  return m_Cat;
}

VS FitReader::GetChannels() const {
  VS list;
  for(auto c : m_Chan)
    list += c.first;
  
  return list;
}

const Systematics& FitReader::GetSystematics() const {
  return m_Sys;
}

string FitReader::GetSignalTitle(const string& label){
  size_t p = label.rfind("_");
  if(p == std::string::npos)
    return label;

  string title = label.substr(0, p);
  int    mass  = stoi(label.substr(p+1,label.length()-p));
  return title+" "+std::to_string(mass/10000)+" "+std::to_string(mass%100000);
}

