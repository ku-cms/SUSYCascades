// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <dirent.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>
//#include <TIter.h>
#include <TKey.h>

#include "ReducedBase.hh"
#include "FitInputBuilder.hh"
#include "Systematics.hh"
#include "SampleTool.hh"
#include "CategoryTool.hh"
#include "ScaleFactorTool.hh"
#include "Leptonic.hh"
#include "Hadronic.hh"

#include "RestFrames/RestFrames.hh"

using namespace std;
using namespace RestFrames;

int main(int argc, char* argv[]) {
int event_count = 0;
  int ifile = -1;
  //string NtuplePath = "/home/t3-ku/z374f439/storage/crogan/";
  string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v0/";
  string OutFile    = "BuildPlotInput_output.root";

  bool doSigFile = false;
  string SigFile = "";

  bool bprint = false;
  int  year   = 2017;
  bool addBkg  = false;
  bool addSig  = false;
  bool addData = false;
  bool extrahist = false;
  
  vector<string> proc_to_add;

  CategoryTool CT;
  CategoryList Categories;

  bool cat0L = false;
  bool cat1L = false;
  bool cat2L = false;
  bool cat3L = false;

  bool setLumi = false;
  double lumi;

  bool doSys = false;

  bool maskSR = false;
  
  for(int i = 0; i < argc; i++){
    if(strncmp(argv[i],"--help", 6) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-h", 2) == 0){
      bprint = true;
    }
    if(strncmp(argv[i],"-path", 5) == 0){
      i++;
      NtuplePath = string(argv[i]);
    }
    if(strncmp(argv[i],"-ifile", 6) == 0){
      i++;
      ifile = std::atoi(argv[i]);
    }
    if(strncmp(argv[i],"-o", 2) == 0){
      i++;
      OutFile = string(argv[i]);
    }
    if(strncmp(argv[i],"--output", 8) == 0){
      i++;
      OutFile = string(argv[i]);
    }
    if(strncmp(argv[i],"-year", 5) == 0){
      i++;
      year = std::atoi(argv[i]);
    }
    if(strncmp(argv[i],"+proc", 5) == 0){
      i++;
      proc_to_add.push_back(argv[i]);
    }
    if(strncmp(argv[i],"++bkg", 5) == 0){
      addBkg = true;
    }
    if(strncmp(argv[i],"++sig", 5) == 0){
      addSig = true;
    }
    if(strncmp(argv[i],"++data", 6) == 0){
      addData = true;
    }
    if(strncmp(argv[i],"++all", 5) == 0){
      addBkg  = true;
      addSig  = true;
      addData = true;
    }
    if(strncmp(argv[i],"-lumi", 5) == 0){
      i++;
      setLumi = true;
      lumi = std::stof(argv[i]);
    }
    if(strncmp(argv[i],"-maskSR", 7) == 0){
      maskSR = true;
    }
  }
      
  if((proc_to_add.size() == 0) &&
     (addBkg  == false) &&
     (addSig  == false) &&
     (addData == false))
    bprint = true;

  if(bprint){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)          print options" << endl;
    cout << "   -path [dest]        path to input ntuples" << endl;
    cout << "   --ouput(-o) [file]  output root file" << endl;
    cout << "   -year [year]        year to process" << endl;
    cout << "   +proc [label]       add processes matching label (can have >1 of these)" << endl;
    cout << "   ++bkg               add all background samples for year" << endl;
    cout << "   ++sig               add all signal samples" << endl;
    cout << "   ++data              add all background samples" << endl;
    cout << "   ++all               add all samples" << endl;
    cout << "   ++sys               turn on available systematics" << endl;
    cout << "   -lumi [lumi]        set luminosity to lumi" << endl;
    cout << "   -maskSR             mask high RISR bins" << endl;
    cout << "Example: ./BuildPlotInput.x ++bkg +proc T2tt" << endl;
   
    return 0;
  }
  
  cout << "Initializing sample maps from path " << NtuplePath << " for year " << year << endl;
  
  SampleTool ST(NtuplePath, year);

  ScaleFactorTool SF;

  ProcessList samples;
  if(addBkg){
    cout << "Adding all background processes" << endl;
    samples += ST.Get(kBkg);
  }
  if(addSig){
    cout << "Adding all signal processes" << endl;
    samples += ST.Get(kSig);
  }
  if(addData){
    cout << "Adding all data for year " << year << endl;
    samples += ST.Get(kData);
  }
  for(int p = 0; p < int(proc_to_add.size()); p++){
    cout << "Adding processes that match \"" << proc_to_add[p] << "\"" << endl;
    samples += ST.Get(proc_to_add[p]);
  }

  Categories      += CT.GetCategories(maskSR);
  
  // Declare observables for histogram booking
  HistPlot* hist_plot = new HistPlot("plot","");

  // variables
  const HistPlotVar& MET   = hist_plot->GetNewVar("MET", "MET", 150., 2000., "[GeV]");
  const HistPlotVar& RISR  = hist_plot->GetNewVar("RISR", "RISR", 0.45, 1.);
  const HistPlotVar& PTISR = hist_plot->GetNewVar("PTISR", "PTISR", 250., 2500.);
  const HistPlotVar& Gperp = hist_plot->GetNewVar("Gperp", "Gperp", 0., 1.);
  const HistPlotVar& Mperp = hist_plot->GetNewVar("Mperp","Mperp", 0., 800.);

  // object counting
  const HistPlotVar& NSjet  = hist_plot->GetNewVar("NSjet","NSjet",-0.5,15.5);
  const HistPlotVar& Nbjet  = hist_plot->GetNewVar("Nbjet","Nbjet",-0.5,7.5);
  const HistPlotVar& NSbjet = hist_plot->GetNewVar("NSbjet","NSbjet",-0.5,7.5);
  const HistPlotVar& NSV    = hist_plot->GetNewVar("NSV","SV",-0.5,7.5);
  const HistPlotVar& Nel    = hist_plot->GetNewVar("Nel","Nel",-0.5,7.5);
  const HistPlotVar& Nmu    = hist_plot->GetNewVar("Nmu","Nmu",-0.5,7.5);

  // object kin
  const HistPlotVar& BjetEta = hist_plot->GetNewVar("BjetEta","BjetEta",-2.4,2.4);
  const HistPlotVar& BjetPhi = hist_plot->GetNewVar("BjetPhi","BjetPhi",-acos(-1.),acos(-1.));
  const HistPlotVar& BjetPt  = hist_plot->GetNewVar("BjetPt","BjetPt",20., 1000.);
  const HistPlotVar& JetEta = hist_plot->GetNewVar("JetEta","JetEta",-2.4,2.4);
  const HistPlotVar& JetPhi = hist_plot->GetNewVar("JetPhi","JetPhi",-acos(-1.),acos(-1.));
  const HistPlotVar& JetPt = hist_plot->GetNewVar("JetPt","JetPt",20., 1000.);
  const HistPlotVar& ElEta = hist_plot->GetNewVar("ElEta","ElEta",-2.4,2.4);
  const HistPlotVar& ElPhi = hist_plot->GetNewVar("ElPhi","ElPhi",-acos(-1.),acos(-1.));
  const HistPlotVar& ElPt = hist_plot->GetNewVar("ElPt","ElPt",5., 200.);
  const HistPlotVar& MuEta = hist_plot->GetNewVar("MuEta","MuEta",-2.4,2.4);
  const HistPlotVar& MuPhi = hist_plot->GetNewVar("MuPhi","MuPhi",-acos(-1.),acos(-1.));
  const HistPlotVar& MuPt = hist_plot->GetNewVar("MuPt","MuPt",3., 200.);
  const HistPlotVar& SVEta = hist_plot->GetNewVar("SVEta","SVEta",-2.4,2.4);
  const HistPlotVar& SVPhi = hist_plot->GetNewVar("SVPhi","SVPhi",-acos(-1.),acos(-1.));
  const HistPlotVar& SVPt = hist_plot->GetNewVar("SVPt","SVPt",2., 20.);
  
  const HistPlotCategory& cat_data_bjet = hist_plot->GetNewCategory("Data-bj", "Data-bj");
  const HistPlotCategory& cat_data_jet = hist_plot->GetNewCategory("Data-j", "Data-j");
  const HistPlotCategory& cat_data_el = hist_plot->GetNewCategory("Data-el", "Data-el");
  const HistPlotCategory& cat_data_mu = hist_plot->GetNewCategory("Data-mu", "Data-mu");
  const HistPlotCategory& cat_data_sv = hist_plot->GetNewCategory("Data-mu", "Data-sv");

  const HistPlotCategory& cat_bkg_bjet = hist_plot->GetNewCategory("Bkg-bj", "Bkg-bj");
  const HistPlotCategory& cat_bkg_jet = hist_plot->GetNewCategory("Bkg-j", "Bkg-j");
  const HistPlotCategory& cat_bkg_el = hist_plot->GetNewCategory("Bkg-el", "Bkg-el");
  const HistPlotCategory& cat_bkg_mu = hist_plot->GetNewCategory("Bkg-mu", "Bkg-mu");
  const HistPlotCategory& cat_bkg_sv = hist_plot->GetNewCategory("Bkg-sv", "Bkg-sv");

  VS scat_incl;
  VS scat_Nlep[4];
  VS scat_gsb;

  scat_incl.a("incl").a("0J").a("1J").a("2J").a("3J").a("4J").a("5J");
  scat_Nlep[0].a("0L").a("0L_0J").a("0L_1J").a("0L_2J").a("0L_3J").a("0L_4J").a("0L_5J");
  scat_Nlep[1].a("1L").a("1L_0J").a("1L_1J").a("1L_2J").a("1L_3J").a("1L_4J");
  scat_Nlep[2].a("2L").a("2L_0J").a("2L_1J").a("2L_2J");
  scat_Nlep[3].a("3L").a("3L_0J").a("3L_1J");
  scat_gsb.a("").a("gold").a("slvr").a("bron");

  RFList<const HistPlotCategory> cat_all;

  map<Process,RFList<const HistPlotCategory> > cat_incl;
  map<Process,RFList<const HistPlotCategory> > cat_0L;
  map<Process,RFList<const HistPlotCategory> > cat_NL[3][4];
  
  // dummy in case there is no data requested
  Process data_obs("data_obs", kData);

  // to check RISR/Mperp bin
  FitInputBuilder FITBuilder;

  bool some_data = false;
  bool some_bkg  = false; 
  
  // sample (process) loop
  int Nsample = samples.GetN();
  for(int s = 0; s < Nsample; s++){
    Process proc = samples[s];

    bool is_data   = (proc.Type() == kData);
    bool is_bkg    = (proc.Type() == kBkg);
    bool is_signal = (proc.Type() == kSig);

    for(auto s : scat_incl){
      string label = proc.Name()+"_"+s;
      const HistPlotCategory& cat = hist_plot->GetNewCategory(label, label);
      if(cat_incl.count(proc) < 1)
	cat_incl[proc] = RFList<const HistPlotCategory>();
      
      cat_incl[proc] += cat;
      cat_all += cat;
    }
    for(auto s : scat_Nlep[0]){
      string label = proc.Name()+"_"+s;
      const HistPlotCategory& cat = hist_plot->GetNewCategory(label, label);
      if(cat_0L.count(proc) < 1)
	cat_0L[proc] = RFList<const HistPlotCategory>();
      
      cat_0L[proc] += cat;
      cat_all += cat;
    }
    for(int l = 1; l < 4; l++){
      for(int q = 0; q < 4; q++){
	for(auto s : scat_Nlep[l]){
	  string label = proc.Name()+"_"+scat_gsb[q]+s;
	  const HistPlotCategory& cat = hist_plot->GetNewCategory(label, label);
	  if(cat_NL[l-1][q].count(proc) < 1)
	    cat_NL[l-1][q][proc] = RFList<const HistPlotCategory>();

	  cat_NL[l-1][q][proc] += cat;
	  cat_all += cat;
	}
      }
    }
    
    if(is_data)
      some_data = true;
    if(is_bkg)
      some_bkg = true;
    
    cout << "processing sample " << proc.Name() << endl;
  }
  
  hist_plot->AddPlot(MET, cat_all);
  hist_plot->AddPlot(RISR, cat_all);
  hist_plot->AddPlot(PTISR, cat_all);
  hist_plot->AddPlot(Gperp, cat_all);
  hist_plot->AddPlot(Mperp, cat_all);
  hist_plot->AddPlot(NSjet, cat_all);
  hist_plot->AddPlot(Nbjet, cat_all);
  hist_plot->AddPlot(NSbjet, cat_all);
  hist_plot->AddPlot(NSV, cat_all);
  hist_plot->AddPlot(Nel, cat_all);
  hist_plot->AddPlot(Nmu, cat_all);
  
  if(some_data){
    hist_plot->AddPlot(BjetEta, cat_data_bjet);
    hist_plot->AddPlot(BjetPhi, cat_data_bjet);
    hist_plot->AddPlot(BjetPt,  cat_data_bjet);
    hist_plot->AddPlot(BjetEta, BjetPhi, cat_data_bjet);
    hist_plot->AddPlot(JetEta, cat_data_jet);
    hist_plot->AddPlot(JetPhi, cat_data_jet);
    hist_plot->AddPlot(JetPt,  cat_data_jet);
    hist_plot->AddPlot(JetEta, JetPhi, cat_data_jet);
    hist_plot->AddPlot(SVEta, cat_data_sv);
    hist_plot->AddPlot(SVPhi, cat_data_sv);
    hist_plot->AddPlot(SVPt,  cat_data_sv);
    hist_plot->AddPlot(SVEta, SVPhi, cat_data_sv);
    hist_plot->AddPlot(ElEta, cat_data_el);
    hist_plot->AddPlot(ElPhi, cat_data_el);
    hist_plot->AddPlot(ElPt,  cat_data_el);
    hist_plot->AddPlot(ElEta, ElPhi, cat_data_el);
    hist_plot->AddPlot(MuEta, cat_data_mu);
    hist_plot->AddPlot(MuPhi, cat_data_mu);
    hist_plot->AddPlot(MuPt,  cat_data_mu);
    hist_plot->AddPlot(MuEta, MuPhi, cat_data_mu);
  }
  if(some_bkg){
    hist_plot->AddPlot(BjetEta, cat_bkg_bjet);
    hist_plot->AddPlot(BjetPhi, cat_bkg_bjet);
    hist_plot->AddPlot(BjetPt,  cat_bkg_bjet);
    hist_plot->AddPlot(BjetEta, BjetPhi, cat_bkg_bjet);
    hist_plot->AddPlot(JetEta, cat_bkg_jet);
    hist_plot->AddPlot(JetPhi, cat_bkg_jet);
    hist_plot->AddPlot(JetPt,  cat_bkg_jet);
    hist_plot->AddPlot(JetEta, JetPhi, cat_bkg_jet);
    hist_plot->AddPlot(SVEta, cat_bkg_sv);
    hist_plot->AddPlot(SVPhi, cat_bkg_sv);
    hist_plot->AddPlot(SVPt,  cat_bkg_sv);
    hist_plot->AddPlot(SVEta, SVPhi, cat_bkg_sv);
    hist_plot->AddPlot(ElEta, cat_bkg_el);
    hist_plot->AddPlot(ElPhi, cat_bkg_el);
    hist_plot->AddPlot(ElPt,  cat_bkg_el);
    hist_plot->AddPlot(ElEta, ElPhi, cat_bkg_el);
    hist_plot->AddPlot(MuEta, cat_bkg_mu);
    hist_plot->AddPlot(MuPhi, cat_bkg_mu);
    hist_plot->AddPlot(MuPt,  cat_bkg_mu);
    hist_plot->AddPlot(MuEta, MuPhi, cat_bkg_mu);
  }

  for(int s = 0; s < Nsample; s++){
    Process proc = samples[s];

    if(doSigFile && proc.Type() == kSig){
      bool keep = false;
      int Nfile = ST.NTrees(proc);
      for(int f = 0; f < Nfile; f++){
	string file = ST.FileName(proc, f);
	if(file.find(SigFile) != string::npos)
	  keep = true;
      }
      if(!keep)
	continue;
    }
    
    string title = proc.Name();

    bool is_data   = (proc.Type() == kData);
    bool is_bkg    = (proc.Type() == kBkg);
    bool is_signal = (proc.Type() == kSig);
    
    int Nfile = ST.NTrees(proc);

    cout << "Processing " << Nfile << " files for process " << title << endl;

    for(int f = 0; f < Nfile; f++){
      if(ifile != -1)
        f = ifile;
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

	if(!base->EventFilter)
	  continue;

	if(base->runnum >= 319077 && is_data && year == 2018)
	  if(base->HEM_Veto)
	    continue;
	
	if(do_FilterDilepton)
	  if(SF.DileptonEvent(base))
	    continue;
	
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

	if(Nlep + NjetS + NSV < 1)
	  continue;

	NSjet = base->Njet_S;
	Nbjet = base->Nbjet_S+base->Nbjet_ISR;
	NSbjet = base->Nbjet_S;
	NSV = base->NSV_S;
	Nel = 0;
	Nmu = 0;
	
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
	  if(abs(PDGID) == 11){
	    flavor = kElec;
	    Nel += 1;
	  } else {
	    flavor = kMuon;
	    Nmu += 1;
	  }
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  LepSource source = LepSource(base->SourceID_lep->at(index));
	  //LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version
	    
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
	  if(abs(PDGID) == 11){
	    flavor = kElec;
	    Nel += 1;
	  } else {
	    flavor = kMuon;
	    Nmu += 1;
	  }
	  LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
	  LepSource source = LepSource(base->SourceID_lep->at(index));
	  //LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version
	  
	  list_b += Lep(flavor, charge, id, source);
	}

	// SV eta
	double SVmaxeta = 1.; // 1 is fine b/c less than 1.5 cutoff
	for(int ie = 0; ie < base->NSV_S; ie++)
	  if(fabs(base->Eta_SV->at(ie)) > SVmaxeta)
	    SVmaxeta = fabs(base->Eta_SV->at(ie));

	// gammaT calc
	double MST =
	  sqrt(base->MX3a_BoostT*base->MX3a_BoostT+base->PX3_BoostT*base->PX3_BoostT) +
	  sqrt(base->MX3b_BoostT*base->MX3b_BoostT+base->PX3_BoostT*base->PX3_BoostT);
	double gammaT = 2.*base->Mperp / MST;



	
	Category Event(Leptonic(list_a, list_b),
		       Hadronic(NjetS, NbjetS, NSV),
		       Hadronic(NjetISR, NbjetISR, base->NSV_ISR));
	Event.AddGenericVal(GenericVal(base->PTISR));
	Event.AddGenericVal(gammaT);
	Event.AddGenericVal(SVmaxeta);
	
	int eindex = Categories.Find(Event);

	if(eindex < 0){
	  continue;
	}
	
	string cname = Categories[eindex].FullLabel();
	
	MET = base->MET;
	RISR = base->RISR;
	PTISR = base->PTISR;
	Gperp = gammaT;
	Mperp = base->Mperp;
	if((1 == Nlep) && (NjetS == 0) && (NSV < 1))
	  Mperp = 2.*base->EL_BoostT;
	if((Nlep < 1) && (NjetS == 0) && (1 == NSV))
	  Mperp = 2.*base->EJ_BoostT;
	if((Nlep < 1) && (NjetS == 1) && (NSV < 1))
	  Mperp = 2.*base->EJ_BoostT;

	bool found = FITBuilder.AddEvent(1., Mperp, RISR,
					 Categories[eindex], data_obs);
	if(!found)
	  continue;
      
      
	double weight = 1.;
	if(!is_data)
	  weight = (setLumi ? lumi : ST.Lumi())*base->weight*sample_weight;

event_count++;    
	// region fills
	hist_plot->Fill(cat_incl[proc][0], weight);

	int ijet = NjetS+1;
	if(ijet > 6) ijet = 6;
	hist_plot->Fill(cat_incl[proc][ijet], weight);

	if(0 == Nlep){
	  hist_plot->Fill(cat_0L[proc][0], weight);
	  hist_plot->Fill(cat_0L[proc][ijet], weight);
	} else {
	  if(1 == Nlep)
	    if(ijet > 5) ijet = 5;
	  if(2 == Nlep)
	    if(ijet > 3) ijet = 3;
	  if(3 == Nlep)
	    if(ijet > 2) ijet = 2;
	  hist_plot->Fill(cat_NL[int(Nlep)-1][0][proc][0], weight);
	  hist_plot->Fill(cat_NL[int(Nlep)-1][0][proc][ijet], weight);

	  int igsb = 1;
	  if(cname.find("slvr") != std::string::npos)
	    igsb = 2;
	  if(cname.find("bron") != std::string::npos)
	    igsb = 3;

	  hist_plot->Fill(cat_NL[int(Nlep)-1][igsb][proc][0], weight);
	  hist_plot->Fill(cat_NL[int(Nlep)-1][igsb][proc][ijet], weight);
	}
	
	// object fills
	for(int j = 0; j < base->Njet; j++){
	  JetPt  = base->PT_jet->at(j);
	  JetEta = base->Eta_jet->at(j);
	  JetPhi = base->Phi_jet->at(j);

	  if(is_data)
	    hist_plot->Fill(cat_data_jet, weight);
	  if(is_bkg)
	    hist_plot->Fill(cat_bkg_jet, weight);

	  if(base->BtagID_jet->at(j) < 3)
	    continue;

	  BjetPt  = base->PT_jet->at(j);
	  BjetEta = base->Eta_jet->at(j);
	  BjetPhi = base->Phi_jet->at(j);

	  if(is_data)
	    hist_plot->Fill(cat_data_bjet, weight);
	  if(is_bkg)
	    hist_plot->Fill(cat_bkg_bjet, weight);
	}

	for(int s = 0; s < NSV; s++){
	  SVPt  = base->PT_SV->at(s);
	  SVEta = base->Eta_SV->at(s);
	  SVPhi = base->Phi_SV->at(s);

	 if(is_data)
	    hist_plot->Fill(cat_data_sv, weight);
	  if(is_bkg)
	    hist_plot->Fill(cat_bkg_sv, weight);
	}

	for(int l = 0; l < Nlep; l++){
	  if(abs(int(base->PDGID_lep->at(l))) == 11){
	    ElPt  = base->PT_lep->at(l);
	    ElEta = base->Eta_lep->at(l);
	    ElPhi = base->Phi_lep->at(l);

	    if(is_data)
	      hist_plot->Fill(cat_data_el, weight);
	    if(is_bkg)
	      hist_plot->Fill(cat_bkg_el, weight);
	  } else {
	    MuPt  = base->PT_lep->at(l);
	    MuEta = base->Eta_lep->at(l);
	    MuPhi = base->Phi_lep->at(l);

	    if(is_data)
	      hist_plot->Fill(cat_data_mu, weight);
	    if(is_bkg)
	      hist_plot->Fill(cat_bkg_mu, weight);
	  }
	}
	
      }
      delete base;
      delete chain;
      if(ifile != -1)
	break;
    }
  }

   hist_plot->WriteHist(OutFile);
  //FITBuilder.WriteFit(OutFile);
std::cout << "Number of Events: " << event_count << endl;
  
}
