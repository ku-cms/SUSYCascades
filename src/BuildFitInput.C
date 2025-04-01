// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <dirent.h>
#include <vector>
#include <cmath>

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

#include "ReducedBase_NTUPLES_v1.hh"
#include "FitInputBuilder.hh"
#include "Systematics.hh"
#include "SampleTool.hh"
#include "CategoryTool.hh"
#include "ScaleFactorTool.hh"
#include "Leptonic.hh"
#include "Hadronic.hh"
#include "METTriggerTool.hh"

using namespace std;

int main(int argc, char* argv[]) {
  int ifile = -1;
  //string NtuplePath = "root://xrootd.unl.edu//store/user/zflowers/crogan/";
  string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";
  string OutFile    = "BuildFitInput_output.root";

  bool doSigFile = false;
  string SigFile = "";

  bool bprint = false;
  int  year   = 2017;
  bool addBkg  = false;
  bool addSig  = false;
  bool addData = false;
  bool extrahist = false;
  bool sigProc = false; 
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
  bool doTreeSys = false;
  bool maskSR = false;
 
  bool debugVerbosity = false;
  bool debugData =false;

  string treeLoad = "";
  string treeSysName ="";
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
    if(strncmp(argv[i],"+hist", 5) == 0){
      extrahist  = true;
    }
    if(strncmp(argv[i],"+cat0L", 6) == 0){
      cat0L = true;
    }
    if(strncmp(argv[i],"+cat1L", 6) == 0){
      cat1L = true;
    }
    if(strncmp(argv[i],"+cat2L", 6) == 0){
      cat2L = true;
    }
    if(strncmp(argv[i],"+cat3L", 6) == 0){
      cat3L = true;
    }
    if(strncmp(argv[i],"++sys", 5) == 0){
      doSys = true;
    }
    if(strncmp(argv[i],"++treesys",9) ==0){
      doTreeSys = true;
    }
    if(strncmp(argv[i],"-lumi", 5) == 0){
      i++;
      setLumi = true;
      lumi = std::stof(argv[i]);
    }
    if(strncmp(argv[i],"-sigfile", 8) == 0){
      i++;
      doSigFile = true;
      SigFile = argv[i];
    }
    if(strncmp(argv[i],"-maskSR", 7) == 0){
      maskSR = true;
    }
    if(strncmp(argv[i],"-treeName", 9) == 0){
 	i++;
	treeLoad = std::string(argv[i]);
	
    } 
    if(strncmp(argv[i],"-treeSysName", 12) ==0){
	 i++;
	  treeSysName = std::string(argv[i]);	
	}
    
    if(strncmp(argv[i],"+SMS", 4)==0){
 	sigProc = true;     
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
    cout << "   -year [year]        index of year to process" << endl;
    cout << "   +proc [label]       add processes matching label (can have >1 of these)" << endl;
    cout << "   ++bkg               add all background samples for year" << endl;
    cout << "   ++sig               add all signal samples" << endl;
    cout << "   ++data              add all background samples" << endl;
    cout << "   ++all               add all samples" << endl;
    cout << "   +cat0L              add 0L categories" << endl;
    cout << "   +cat1L              add 1L categories" << endl;
    cout << "   +cat2L              add 2L categories" << endl;
    cout << "   +cat3L              add 3L categories" << endl;
    cout << "   ++sys               turn on available systematics" << endl;
    cout << "   +hist               book 2D histograms also" << endl;
    cout << "   -lumi [lumi]        set luminosity to lumi" << endl;
    cout << "   -sigfile            signal filename must match this string to be included" << endl;
    cout << "   -maskSR             mask high RISR bins" << endl;
    cout << "Example: ./BuildFitInput.x ++bkg +proc T2tt +cat1L ++sys" << endl;
   
    return 0;
  }

  if(cat0L)
    Categories += CT.GetCategories_0L(maskSR);
  if(cat1L)
    Categories += CT.GetCategories_1L(maskSR);
  if(cat2L)
    Categories += CT.GetCategories_2L(maskSR);
  if(cat3L)
    Categories += CT.GetCategories_3L(maskSR);
  
  cout << "Initializing sample maps from path " << NtuplePath << " for year " << year << endl;
  
  SampleTool ST(NtuplePath, year);

  ProcessList samples;
  if(addBkg){
    cout << "Adding all background processes" << endl;
    samples += ST.Get(kBkg);
  }
  if(addSig){
    cout << "Adding all signal processes" << endl;
    samples += ST.Get(kSig);
    //samples += ST.GetStrictSignalMatch( kSig);

  }
  if(addData){
    cout << "Adding all data for year " << year << endl;
    samples += ST.Get(kData);
  }
  for(int p = 0; p < int(proc_to_add.size()); p++){
    cout << "Adding processes that match \"" << proc_to_add[p] << "\"" << endl;
    //samples += ST.Get(proc_to_add[p]);
    //if its signal do strict match
    if(sigProc){
     samples += ST.GetStrictSignalMatch(proc_to_add[p]);
    }else{
    //otherwise do the normal thing 
     samples += ST.Get(proc_to_add[p]); 
    }
  } 

  if(Categories.GetN() == 0)
    Categories += CT.GetCategories(maskSR);
  
  //if a lepton region wasn't specified, turn them all on 
  if(!cat0L && !cat1L && !cat2L && !cat3L){
    cat0L = true;
    cat1L = true;
    cat2L = true;
    cat3L = true;
  }
  
  // cout << "Categories:" << endl;
  // Categories.Print();

  SystematicsTool SYS;
  
  METTriggerTool m_METTriggerTool;
  //if you run interactively, use csv/METTrigger path
  //for batch runs the csv file will be copied in the base directory
  m_METTriggerTool.BuildMap("Parameters.csv");
  //m_METTriggerTool.BuildMap("csv/METTrigger/Parameters.csv");

  ScaleFactorTool SF;
  SF.AddBtagFolder("./BtagSF");
  SF.init_bLambda();

  Systematics systematics(1);
 // Systematics nominalTreeSys(1);
  //Systematics treeSysList;//associated sys which contain Up or Down
  //this is the tree loading for bkg, ill do a separate one for signals
   Systematic treeSys(treeSysName); 

   if(doTreeSys){
	size_t found = treeLoad.find("Down");
	if(found != std::string::npos){
		std::cout<<"set "<<treeSysName<<" to Down\n";
		treeSys.Down();
	}
	else{
		std::cout<<"set "<<treeSysName<<" to Up\n";
		treeSys.Up();
	}
	size_t found_default = treeLoad.find("default");
	if(found_default != std::string::npos){
		treeLoad = "KUAnalysis";
		Systematics inputSys(0);
		inputSys += treeSys;
		systematics = inputSys;
	}
   }
   

//std::vector<std::string> TreesToProcess;//list of tree names to process
  //TreesToProcess.push_back("KUAnalysis");
  //TreesToProcess.push_back(treeLoad);
//  if(doSys)
 //   systematics += SYS.GetWeightSystematics();

 // if(doTreeSys){
  //  treeSysList = SYS.GetTreeSystematics();
  //  std::vector<std::string> tmp;
  //  tmp = SYS.GetTreeSystematicsNames();
  //  TreesToProcess.insert(TreesToProcess.end(),tmp.begin(),tmp.end());
   
 // }
  std::cout<<"Processing these trees:\n";
//  for(int t=0; t<TreesToProcess.size(); t++){
	//std::cout<<TreesToProcess[t]<<" ";
	std::cout<<treeLoad<<" ";
      
 // }
 // std::cout<<"checking systematic scope\n";
 // for(int t=0; t<treeSysList.GetN(); t++){
	/// std::cout<<" tree systematic "<<treeSysList[t].Label() <<" is: up="<<treeSysList[t].IsUp()<<" down="<<treeSysList[t].IsDown()<<"\n";
 // }

 std::cout<<"\n";
  	

  FitInputBuilder FITBuilder(extrahist);

  // dummy in case there is no data requested
  Process data_obs("data_obs", kData);

  // sample (process) loop
  int Nsample = samples.GetN();
 /* for(int s = 0; s < Nsample; s++){
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
    cout << "processing sample " << proc.Name() << endl;
  }*/
  //debugging loop
  if(sigProc){
    cout<<"Samples for processing: \n";
    for(int s = 0; s< Nsample; s++){
	cout<< samples[s].Name() << " ";
    }
    cout<<"\n";
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
    //cout << "kSig ="<<kSig<<" \n";
    cout << "processing sample " << proc.Name() << endl;
    string title = proc.Name();

    bool is_data   = (proc.Type() == kData);
    bool is_bkg    = (proc.Type() == kBkg);
    bool is_signal = (proc.Type() == kSig);

    int Nsys = (is_data ? 1 : systematics.GetN()); //weight systematics
    
    int Nfile = ST.NTrees(proc);

    cout << "Processing " << Nfile << " files for process " << title << endl;
    cout <<" The proc is: "<< proc.Name() << endl;
    for(int f =0; f<Nfile; f++){
     cout<< ST.FileName(proc,f) << " " << ST.TreeName(proc,f) << endl;
    }

    for(int f = 0; f < Nfile; f++){
      if(ifile != -1)
        f = ifile;
      string file = ST.FileName(proc, f);
      string tree = ST.TreeName(proc, f);

     //init Bmap here for per file processing
//      SF.AddBtagFolder("./BtagSF", file, year );
     double blambda0;
     double blambda12;
     std::vector<double> blams;
     if(proc.Name() != "ttbar" && !sigProc){
     blams = SF.GetbLambdas(proc.Name(), file, year);
     blambda0 = blams[0];
     blambda12 = blams[1];
     }
     else{
	std::cout<<"signal or ttbar detected, setting blambda = 1\n";
	blambda0=1.;
	blambda12 =1.;
     }
     
	
      bool is_FastSim = ST.IsFastSim(proc, f);
      bool do_FilterDilepton = ST.FilterDilepton(proc, f);
      double sample_weight = ST.GetSampleWeight(proc, f);
      SleptonFlavor do_FilterSleptons = ST.FilterSleptons(proc, f);

      if(is_signal)
	sample_weight *= SF.GetX20BRSF(file, tree);
      
      cout << "   Processing file " << file << endl; //<< " w/ tree " << tree << endl;
      cout << "      Sample weight is " << sample_weight << endl;
      if(is_FastSim)
	cout << "      Is FastSim" << endl;
      if(do_FilterDilepton)
	cout << "      Filter Out dilepton events" << endl;
      if(do_FilterSleptons == kSmu)
        cout << "      Filter out events containing smuons" << endl;
      if(do_FilterSleptons == kSel)
        cout << "      Filter out events containing selectrons" << endl;

  //  int treeSysCtr=0;
  //  for(int itree=0; itree<TreesToProcess.size(); itree++){ ///begin tree loop   
	 
//	tree = TreesToProcess[itree]; //rename this guy
//	std::cout<<"Processing Tree: "<<TreesToProcess[itree]<<"\n";  
	
//    TChain* chain = ST.Tree(proc, f);
	//debug dev testing
//	TChain* chain = ST.Tree(proc, f, "KUAnalysis_JESUncer_TotalUp");
	
	//TChain* chain = ST.Tree(proc, f, TreesToProcess[itree]);
	TChain* chain;
	if( doTreeSys ){
		chain = ST.Tree(proc, f, treeLoad);
	}
	else{
		chain = ST.Tree(proc, f);
	}
        cout<<" Chain Name "<< chain->GetName() << endl;
	//set up the correct systematic to use
	//always process the nominal kuanalysis tree first
	//replace with nominal systematics if not nominal tree
//	if(doTreeSys && (TreesToProcess[itree] != "KUAnalysis")){
//		systematics = nominalTreeSys;
//		std::cout<<"replacing syslist with "<<nominalTreeSys[0].Label()<<" is: up="<<nominalTreeSys[0].IsUp()<<" down="<<nominalTreeSys[0].IsDown()<<"\n";
//		
//	}
//	if( itree >= 1 ){
//		if( itree %2 != 0 ) treeSysList[treeSysCtr].Up();
//		if( itree %2 == 0 ) treeSysList[treeSysCtr].Down();
//		 std::cout<<" this tree systematic "<<treeSysList[treeSysCtr].Label() <<" is: up="<<treeSysList[treeSysCtr].IsUp()<<" down="<<treeSysList[treeSysCtr].IsDown()<<"\n";
//	}
	

      ReducedBase* base = new ReducedBase(chain);

      int Nentry = base->fChain->GetEntries();
      
      int SKIP = 1;

      // event loop
      for(int e = 0; e < Nentry; e += SKIP){
//	cout<<"in event loop with Nentry "<<e<<", "<<Nentry<<"\n";
	base->GetEntry(e);

	if((e/SKIP)%(std::max(1, int(Nentry/SKIP/10))) == 0)
	  cout << "      event " << e << " | " << Nentry << endl;

        if(base->Charge_lep->size() != base->Nlep ) {cout << "branches be fucked..." << endl; return -1;}

	if(!base->EventFilter)
	  continue;

        if(base->runnum >= 319077 && is_data && year == 2018)
          if(base->HEM_Veto)
            continue;
	
	if(do_FilterDilepton)
	  if(SF.DileptonEvent(base))
	    continue;

	//std::cout << "BEFORE genNsusy: " << base->genNsusy << std::endl;
        if(do_FilterSleptons == kSmu && SF.SleptonEvent(*(base->genPDGID_susy)) == kSmu)
         continue;

	if(do_FilterSleptons == kSel && SF.SleptonEvent(*(base->genPDGID_susy)) == kSel)
         continue;
	
	// apply trigger to data and FullSim events
	if(!base->METORtrigger && !is_FastSim)
	  continue;
		
	if(base->MET < 150)
	  continue;
	  
	if(base->PTISR < 250.)
	  continue;

	if(base->RISR < 0.45)
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
	  
	if(base->RISR < 0.45 || base->RISR > 1.0)
	  continue;

	if(fabs(base->dphiMET_V) > acos(-1.)/2.)
	  continue;
	
	int Nlep     = base->Nlep;
	int NjetS    = base->Njet_S;
	int Nbjet    = base->Nbjet;
	int NbjetS   = base->Nbjet_S;
	int NjetISR  = base->Njet_ISR;
	int NbjetISR = base->Nbjet_ISR;
	int NSV      = 0.; //base->NSV_S;

	if(Nlep + NjetS + NSV < 1)
	  continue;
	  
	LepList list_a;
	LepList list_b;

        std::vector<TLorentzVector> tlv_a;
        std::vector<TLorentzVector> tlv_b;
        double lep_pt,lep_eta,lep_phi,lep_m;
	  
	int index;

	  
	for(int i = 0; i < base->Nlep_a; i++){
	  index = (*base->index_lep_a)[i];
//	  std::cout<<"processing a lepton of index, Nlep_a(i): "<<i<<", "<<index<<"\n";   
	  int PDGID = base->PDGID_lep->at(index);
	    
	  LepID id;

          // "precalc"
//	  if(base->LepQual_lep->at(index) == 2)
	//debugging assume index problem has been fixed.. this will not work with older ntuples (applied in A and B)
//	  if(base->ID_lep->at(index*2) < 3 ||
	//    if(base->ID_lep->at(index) < 3 ||
//	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
//	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)
//	    id = kBronze;
//	  else if(base->LepQual_lep->at(index) == 1)
//	    id = kSilver;
//	  else
//	    id = kGold;
          // "on the fly"
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
	  //LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for old ntuple version
	    

	//  LepSource source = LepSource(base->SourceID_lep->at(index));
	//  LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version (this is the one turned on in master)
	  
	

	  list_a += Lep(flavor, charge, id, source);

	  lep_pt = base->PT_lep->at(index);
	  lep_eta = base->Eta_lep->at(index);
	  lep_phi = base->Phi_lep->at(index);
	  lep_m = base->M_lep->at(index);

          TLorentzVector tlv;
	  tlv.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,lep_m);
	  tlv_a.push_back(tlv);

	}
	for(int i = 0; i < base->Nlep_b; i++){
	  index = (*base->index_lep_b)[i];
  //	  std::cout<<"processing a lepton of index, Nlep_b(i): "<<i<<", "<<index<<"\n";
	  
	  int PDGID = base->PDGID_lep->at(index);

	  LepID id;

          // "precalc"
	 // if(base->LepQual_lep->at(index) == 2)
//	  if(base->ID_lep->at(index*2) < 3 || //index fixed for newly produced ntuples
	 //   if(base->ID_lep->at(index) < 3 ||
//	     base->MiniIso_lep->at(index)*base->PT_lep->at(index) >= 4. ||
//	     base->RelIso_lep->at(index)*base->PT_lep->at(index) >= 4.)

	 //   id = kBronze;
	//  else if(base->LepQual_lep->at(index) == 1)
	//    id = kSilver;
//	  else
//	    id = kGold;
          // "on the fly"
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
	  //LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for old ntuple version
	  
	  list_b += Lep(flavor, charge, id, source);


	//  LepSource source = LepSource(base->SourceID_lep->at(index));
	//  LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version
	//  list_b += Lep(flavor, charge, id, source);
	

	  lep_pt = base->PT_lep->at(index);
          lep_eta = base->Eta_lep->at(index);
          lep_phi = base->Phi_lep->at(index);
          lep_m = base->M_lep->at(index);

          TLorentzVector tlv;
          tlv.SetPtEtaPhiM(lep_pt,lep_eta,lep_phi,lep_m);
          tlv_b.push_back(tlv);

	}

	//loop over both lep lists and form all OSSF pairs
	//calculate all combinations mass. veto any event in j/psi window and break
	bool jpsi=false;
	bool upsilon=false;
	
	//need at least 2 Leps to try this
	if( base->Nlep >= 2 ){
	for( int i=0; i<list_a.GetN(); i++){
		for( int j=0; j<list_b.GetN(); j++){
			if(  (list_a[i].Flavor() == list_b[j].Flavor()) && (list_a[i].Charge() != list_b[j].Charge()) ){
				//OSSF pair calculate mass indexed by list_a(b)
				TLorentzVector tlv_ab = tlv_a[i] + tlv_b[j];
				//std::cout<<"MASS= "<<tlv_ab.M()<<"\n";	
				if( tlv_ab.M() < 3.2 && tlv_ab.M() > 3.0 ){
					//jpsi is present, veto event
					jpsi=true;
				// 	std::cout<<"tlv_ab M: "<<tlv_ab.M()<<" ";
				//	std::cout<<"found jpsi \n";	
				}
				//if( tlv_ab.M() < 10.5 && tlv_ab.M() > 9.0){
				//	upsilon=true;
				//	std::cout<<"tlv_ab M: "<<tlv_ab.M()<<" ";
				//	std::cout<<"found upsilon \n";
				//}
			}	
			if(jpsi || upsilon) break;
		}
		if(jpsi || upsilon) break;
	}
	}//end 2L check
	//veto event if flag flipped
	if(jpsi || upsilon) continue;

	// SV eta
	// double SVmaxeta = 1.; // 1 is fine b/c less than 1.5 cutoff
	// for(int ie = 0; ie < base->NSV_S; ie++)
	//   if(fabs(base->Eta_SV->at(ie)) > SVmaxeta)
	//     SVmaxeta = fabs(base->Eta_SV->at(ie));

	// gammaT calc
	//double MST =
	//  sqrt(base->MX3a_BoostT*base->MX3a_BoostT+base->PX3_BoostT*base->PX3_BoostT) +
	//  sqrt(base->MX3b_BoostT*base->MX3b_BoostT+base->PX3_BoostT*base->PX3_BoostT);
	//double gammaT = 2.*base->Mperp / MST;
	double gammaT = base->gammaT;


	int NSV_ISR = 0;
	double SVmaxeta = 1.;
	Category Event(Leptonic(list_a, list_b),
		       Hadronic(NjetS, NbjetS, NSV),
		       Hadronic(NjetISR, NbjetISR, NSV_ISR));
	Event.AddGenericVal(GenericVal(base->PTISR));
	Event.AddGenericVal(gammaT);
	Event.AddGenericVal(SVmaxeta);
	
	int eindex = Categories.Find(Event);

	if(eindex < 0){
	  continue;
	}

	double PTISR = base->PTISR;
        double HT = base->HT_eta5;
	double weight = 1.;
	double btag_weight = 1.;
	double PU_weight = 1.;
	double trig_weight = 1.;
        double PDF_weight = 1.;
        double MuR_weight = 1.;
        double MuF_weight = 1.;

	//lep sf
	double elID_weight = 1.;
        double elIso_weight = 1.;
        double elSIP_weight = 1.;
	double elVL_weight = 1.;
	
	double muID_weight = 1.;
        double muIso_weight = 1.;
        double muSIP_weight = 1.;
        double muVL_weight = 1.;
        double SF_weight = 1.;
	double prefire_weight = 1.;
	double blambda0_weight = 1.;
	double blambda1_weight = 1.;
	double evt_debug_weight =1.;
	
	if(!is_data){
	  if(year == 2018){//adjust lumi in 2018
		//if event is subject to veto, give preHEM lumi
		// if event is not vetoed by HEM, give full lumi
		weight = (setLumi ? lumi : ST.HEMLumi( base->HEM_Veto ))*base->weight*sample_weight;
	  }//end year check
	  else{// else do what we used to do
		weight = (setLumi ? lumi : ST.Lumi())*base->weight*sample_weight;
	  }//end else other years
	}//end if is data check
	evt_debug_weight = weight;
	// systematics loop
	// do down sys first
	string correct_sys = "";
	
	//for(int is = 0; is < Nsys; is++){
	//  Systematic& sys = systematics[is];
	//  Systematic sys = treeSys;
	   Systematic sys;
	   if(doTreeSys){
		sys = treeSys;
          }
	 else{//do default
		sys = systematics[0];
	 }
//	  std::cout<<"systemtatics list "<< systematics[is].Label()<<" isUp:"<<systematics[is].IsUp()<<" isDown:"<<systematics[is].IsDown()<<"\n";
//  	            std::cout<<"systemtatics list "<< sys.Label()<<" isUp:"<<sys.IsUp()<<" isDown:"<<sys.IsDown()<<"\n";

	/* 
	 if(!(!sys)  ){
	    if(sys.IsUp()){
	      sys.Down();
	      is--;
	    } else {
	      sys.Up();
	    }
	  }
	 */

	  btag_weight = 1.;
	  PU_weight = 1.;
	  trig_weight = 1.;
          PDF_weight = 1.;
          MuR_weight = 1.;
          MuF_weight = 1.;
	  prefire_weight = 1.;
          blambda0_weight = 1.;
          blambda1_weight = 1.;

	// if(!(!sys) && is_data) continue;      
//	if(is_data) continue;no more loop  dont continue
/*
	    trig_weight = m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0);
            if(is_FastSim)
	      trig_weight = m_METTriggerTool.Get_EFF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0)*
            	  m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0);
           //trig ntuples
           // trig_weight = base->MetTrigSFweight; 
            if(sys.Label().find("MET_TRIG") != std::string::npos && proc.Name() != "QCD")
            {
              //trig on the fly
            	        if(sys.IsUp())
            	          if(is_FastSim)
            	            trig_weight = m_METTriggerTool.Get_EFF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1)*
            	              m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1);
            	          else
            	            trig_weight = m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1);
            	        else
            	          if(is_FastSim)
            	            trig_weight = m_METTriggerTool.Get_EFF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1)*
            	              m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1);
            	          else
            	            trig_weight = m_METTriggerTool.Get_SF(base->MET, HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1);
             //trig ntuples
           // 	if(sys.IsUp())
           //      trig_weight = base->MetTrigSFweight_up;
           //     else
           //      trig_weight = base->MetTrigSFweight_down;
            
            }
*/	      
	  //
	  // BTAG systematics from the ntuples
	  //
	    if(sys == Systematic("BTAGHF_SF"))
	      if(sys.IsUp())
	        btag_weight *= base->BtagHFSFweight_up;
	      else
	        btag_weight *= base->BtagHFSFweight_down;
	    else 
	      btag_weight *= base->BtagHFSFweight;

	    if(sys == Systematic("BTAGLF_SF"))
	      if(sys.IsUp())
	        btag_weight *= base->BtagLFSFweight_up;
	      else
	        btag_weight *= base->BtagLFSFweight_down;
	    else 
	      btag_weight *= base->BtagLFSFweight;


	 

	//btag checking
	double otfn,otfu,otfd;
	
	if( e < 50 ){
	/* std::cout<<e<<" njet for this event: "<< base->Njet<<"\n";

	 std::cout<<e<<" base btagHF_SF (n,u,d) "<<base->BtagHFSFweight<<" "<< base->BtagHFSFweight_up <<" "<< base->BtagHFSFweight_down <<"\n";
	 std::cout<<e<<" base btagLF_SF (n,u,d) "<<base->BtagLFSFweight<<" "<< base->BtagLFSFweight_up <<" "<< base->BtagLFSFweight_down <<"\n";
//	 otfn = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, true, 0, kMedium,base->BtagHFSFweight);
//         otfu = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, true, 1, kMedium,base->BtagHFSFweight_up );
//	 otfn = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, true, -1,kMedium,base->BtagHFSFweight_down );
 	 std::cout<<e<<" otf  btagHF_SF (n,u,d) "<<otfn<<" "<< otfu<<" "<< otfd <<"\n";	      
//	 otfn = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, false, 0,kMedium,base->BtagLFSFweight );
//	 otfu = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, false, 1,kMedium,base->BtagLFSFweight );
//	 otfd = SF.GetBtagSFWeight(base->PT_jet, year, is_FastSim, false, -1,kMedium,base->BtagLFSFweight );
	 std::cout<<e<<" otf  btagLF_SF (n,u,d) "<<otfn<<" "<< otfu<<" "<< otfd<<"\n";
  	*/
	}

	  //
	  // BTAG systematics on the fly (needs jet collection in reduced ntuples)
	  //
	 /* 
	  if(sys == Systematic("BTAGHF_SF")){
	    if(sys.IsUp())
	      btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, true, 1);
	    else
	      btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, true, -1);
	  } else {
	    btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, true, 0);
	  }
	  
	  if(sys == Systematic("BTAGLF_SF")){
	    if(sys.IsUp())
	      btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, false, 1);
	    else
	      btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, false, -1);
	  } else {
	    btag_weight *= SF.GetBtagSFWeight(base, year, is_FastSim, false, 0);
	  }
          */

	  // turn off PU systematics for now
	   if(sys == Systematic("PU_SF"))
	     if(sys.IsUp())
	   	PU_weight = base->PUweight_up;
	     else
	   	PU_weight = base->PUweight_down;
	   else
	     PU_weight = base->PUweight;

	  //
	  // PDF systematics
	  //
	    if(sys == Systematic("PDF_SF"))
	      if(sys.IsUp())
	        PDF_weight = base->PDFweight_up;
	      else
	        PDF_weight = base->PDFweight_down;
	    else 
	      PDF_weight = base->PDFweight;

	    if(sys == Systematic("MuR_SF"))
	      if(sys.IsUp())
	        MuR_weight = base->MuRweight_up;
	      else
	        MuR_weight = base->MuRweight_down;
	    else 
	      MuR_weight = base->MuRweight;

	    if(sys == Systematic("MuF_SF"))
	      if(sys.IsUp())
	        MuF_weight = base->MuFweight_up;
	      else
	        MuF_weight = base->MuFweight_down;
	    else 
	      MuF_weight = base->MuFweight;
	 
	   if(sys == Systematic("elID_SF"))
		if(sys.IsUp())
		  elID_weight =base->elIDSFweight_up;
		else
		  elID_weight =base->elIDSFweight_down;
	   else
	     elID_weight =base->elIDSFweight;
	  
           if(sys == Systematic("elIso_SF"))
                if(sys.IsUp())
                  elIso_weight =base->elISOSFweight_up;
                else
                  elIso_weight =base->elISOSFweight_down;
           else
             elIso_weight =base->elISOSFweight;

	   if(sys == Systematic("elSIP_SF"))
                if(sys.IsUp())
                  elSIP_weight =base->elSIPSFweight_up;
                else
                  elSIP_weight =base->elSIPSFweight_down;
           else
             elSIP_weight =base->elSIPSFweight;
		 
	  if(sys == Systematic("elVL_SF"))
                if(sys.IsUp())
                  elVL_weight =base->elVLSFweight_up;
                else
                  elVL_weight =base->elVLSFweight_down;
           else
             elVL_weight =base->elVLSFweight;


	  if(sys == Systematic("muID_SF"))
                if(sys.IsUp())
                  muID_weight =base->muIDSFweight_up;
                else
                  muID_weight =base->muIDSFweight_down;
           else
             muID_weight =base->muIDSFweight;

           if(sys == Systematic("muIso_SF"))
                if(sys.IsUp())
                  muIso_weight =base->muISOSFweight_up;
                else
                  muIso_weight =base->muISOSFweight_down;
           else
             muIso_weight =base->muISOSFweight;

           if(sys == Systematic("muSIP_SF"))
                if(sys.IsUp())
                  muSIP_weight =base->muSIPSFweight_up;
                else
                  muSIP_weight =base->muSIPSFweight_down;
           else
             muSIP_weight =base->muSIPSFweight;

          if(sys == Systematic("muVL_SF"))
                if(sys.IsUp())
                  muVL_weight =base->muVLSFweight_up;
                else
                  muVL_weight =base->muVLSFweight_down;
           else
             muVL_weight =base->muVLSFweight;

	if(sys == Systematic("prefire_SF"))
             if(sys.IsUp()){
                prefire_weight = base->PrefireWeight_up;
	//	std::cout<<"apply prefire up:"<< prefire_weight<<"\n";
		}
             else{
                prefire_weight = base->PrefireWeight_down;
	//	std::cout<<"apply prefire down:"<< prefire_weight<<"\n";
		}
           else{
             prefire_weight = base->PrefireWeight;
	 //    std::cout<<"apply prefire nom:"<< prefire_weight<<"\n";
	   }



	//  weight *= btag_weight*PU_weight;
	  //if(MuR_weight == 0.) MuR_weight = 1.;
	  //if(MuF_weight == 0.) MuF_weight = 1.;
	  //weight *= btag_weight*PU_weight*PDF_weight*MuR_weight*MuF_weight;
	  //weight *= btag_weight*PU_weight*trig_weight;
	 // if(btag_weight*PU_weight*trig_weight*PDF_weight*MuR_weight*MuF_weight*elID_weight*elIso_weight*elSIP_weight*elVL_weight*muID_weight*muIso_weight*muSIP_weight*muVL_weight)
	 // weight = 1.;

	//check for SF nans
	if( isnan( btag_weight ) ) std::cout<<"NaN SF 1!!\n";
 	if( isnan( PU_weight  ) ) std::cout<<"NaN SF 2!!\n";
	if( isnan( trig_weight	 ) ) std::cout<<"NaN SF 3!!\n";
	if( isnan( PDF_weight  ) ) std::cout<<"NaN SF 4!!\n";
	if( isnan( MuR_weight  ) ) std::cout<<"NaN SF 5!!\n";
	if( isnan( MuF_weight  ) ) std::cout<<"NaN SF 6!!\n";
	if( isnan( elID_weight  ) ) std::cout<<"NaN SF 7!!\n";
	if( isnan( elIso_weight  ) ) std::cout<<"NaN SF 8!!\n";
	if( isnan( elSIP_weight  ) ) std::cout<<"NaN SF 9!!\n";
	if( isnan( elVL_weight  ) ) std::cout<<"NaN SF 10!!\n";
	if( isnan( muID_weight  ) ) std::cout<<"NaN SF 11!!\n";
	if( isnan( muIso_weight  ) ) std::cout<<"NaN SF 12!!\n";
	if( isnan( muSIP_weight  ) ) std::cout<<"NaN SF 13!!\n";
	if( isnan( muVL_weight  ) ) std::cout<<"NaN SF 14!!\n";
	if( isnan( prefire_weight) ) std::cout<<"Nan SF 15!!\n";
	//hack for build 110 - remove METtrig SF
	//build 115 everything but mettriiger, with 0 suppression
	//trig_weight=1.;	
	//hack PU weight to be off
	PU_weight=1.;
	SF_weight *= btag_weight*PU_weight*trig_weight*PDF_weight*MuR_weight*MuF_weight*elID_weight*elIso_weight*elSIP_weight*elVL_weight*muID_weight*muIso_weight*muSIP_weight*muVL_weight*prefire_weight;
//	std::cout<<"Prefire weight= "<<prefire_weight<<" SF_weight="<< SF_weight<<"\n";
	if( SF_weight<0.) SF_weight = 0.;
	if( isnan( SF_weight )) SF_weight = 0.;	
 
	//hack for build 111 - remove all nominal SF weight
//	SF_weight=1.;	
//	//hack for 112 theoretical sf
//	SF_weight *=PDF_weight*MuR_weight*MuF_weight;
//	// B114
//	if( SF_weight<0.) SF_weight = 0.;

	//hack for 113 lep sf
	//SF_weight *= elID_weight*elIso_weight*elSIP_weight*elVL_weight*muID_weight*muIso_weight*muSIP_weight*muVL_weight;

	//build 115 everything but mettriiger, with 0 suppression
	
	// get the by processi
	double Njet_nob = NjetS+NjetISR - Nbjet;
	blambda0_weight= pow(blambda0,Nbjet);
	blambda1_weight= pow(blambda12,Njet_nob);
	//if(Nbjet == 0) blambda0_weight=1.;
	//if(Njet_nob == 0) blambda1_weight=1.;
	
	double blambda_combined_weight = blambda0_weight* blambda1_weight;
	//std::cout<<"btag correction: "<<blambda_combined_weight<<" "<<Njet_nob<<" "<<Nbjet<<"\n";

	weight *= SF_weight;
	weight *= blambda_combined_weight;
	 
	 //if( weight < 0 ) weight = 0.;
	  //if(weight > 3 ){
	//	std::cout<<"Event "<<e<<" Large weight! "<<weight<<"\n";
	//	std::cout<<"btag:"<<btag_weight<<" PU:"<<PU_weight<<" TRG:"<<trig_weight<<" PDF:"<<PDF_weight<<" MuR:"<<MuR_weight<<" MuF:"<<MuF_weight<<" eID:"<<elID_weight<<" eIso:"<<elIso_weight<<" eSip:"<<elSIP_weight<<" eVL:"<<elVL_weight<<" muID:"<<muID_weight<<" muIso:"<<muIso_weight<<" muSip:"<<muSIP_weight<<" muVL:"<<muVL_weight<<std::endl;
	//	}
	  if(is_data) weight = 1.;
	  
	  LepList Fakes  = list_a.GetFakes();
	  Fakes         += list_b.GetFakes();
	  
	  double Mperp = base->Mperp;
	  
	  // use Eperp
	  if((Nlep == 1) && (NjetS == 0) && (NSV == 0))
	    Mperp = 2.*base->EL_BoostT;
	  if((Nlep == 0) && (NjetS == 0) && (NSV == 1))
	    Mperp = 2.*base->EJ_BoostT;
	  if((Nlep == 0) && (NjetS == 1) && (NSV == 0))
	    Mperp = 2.*base->EJ_BoostT;
	  
	
	  double RISR  = base->RISR;
        //weight fixing for debug samples
        //weight = 1.;		

	  if(Fakes.GetN() > 0 && is_bkg){
	    VS flabels = Fakes.GetFakeLabels(2); // processes w/ up to 2 "fake" leps
	    int Nf = flabels.size();
	  
	    for(int fl = 0; fl < Nf; fl++){
	      // if(title.find("QCD") == string::npos)
	      // 	FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR,
	      // 			    Categories[eindex], FITBuilder.FakeProcess(flabels[fl]), sys);
	      
             // if(TreesToProcess[itree] == "KUAnalysis"){
	     // FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR,
	//			  Categories[eindex], proc.FakeProcess(flabels[fl]), sys);
          //     }
	     // if(doTreeSys && (TreesToProcess[itree] != "KUAnalysis")){
	        if(doTreeSys){
		//	FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR,
               //                   Categories[eindex], proc.FakeProcess(flabels[fl]), treeSysList[treeSysCtr]);
			FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR, Categories[eindex], proc.FakeProcess(flabels[fl]), treeSys);
		}
		else{
			
			FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR, Categories[eindex], proc.FakeProcess(flabels[fl]), sys);
		 }

	//	}
		if(debugVerbosity){
			std::cout<<"Adding fakes event:"<<e<<" weight: "<<weight/double(Nf)<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].Label()<<"  flabel:"<<flabels[fl]<<"\n";
		}
	    }
	  } else {
		//std::cout<<"adding event "<< weight <<" "<< Mperp <<" "<< RISR <<" "<<Categories[eindex].Label()<<" "<<proc.Name()<<" "<<sys.Label()<<"\n";
	     // if(TreesToProcess[itree] == "KUAnalysis"){
	//	FITBuilder.AddEvent(weight, Mperp, RISR,
	//			Categories[eindex], proc, sys);
	//	}
	     // if(doTreeSys && (TreesToProcess[itree] != "KUAnalysis")){
	        if(doTreeSys){
	//	 FITBuilder.AddEvent(weight, Mperp, RISR,
          //                      Categories[eindex], proc,treeSysList[treeSysCtr]);
          	 FITBuilder.AddEvent(weight, Mperp, RISR, Categories[eindex], proc, treeSys);
		}
		else{
			  FITBuilder.AddEvent(weight, Mperp, RISR,Categories[eindex], proc, sys);
 		 }

		if(debugVerbosity){
			std::cout<<"Adding event:"<<e<<" weight: "<<weight<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].Label()<<" sysLabel:"<<sys.Label()<<"\n";
			std::cout<<"evt weight:"<<evt_debug_weight<<" SFweight: "<<SF_weight<<" lambda weight: "<<blambda_combined_weight<<"\n";
		}
		if(debugData){
			//std::cout<< Categories[eindex].FullLabel() << "\n";
			if( (RISR > 0.95) && (Categories[eindex].FullLabel() == "Ch1L_elpm_bron_1jge1svS_ge1jISR_PTISR0_gamT0_SVeta1"))
			std::cout<<"Adding event:"<<e<<" "<<base->luminum<<" "<<base->runnum <<" weight: "<<weight<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].FullLabel()<<" sysLabel:"<<sys.Label()<<"\n";
			if( base->luminum == 468 && base->runnum ==317640 && PTISR==1012.03){
			std::cout<<"FOUND MIGRATED EVENT!\n";
			std::cout<<"Adding event:"<<e<<" "<<base->luminum<<" "<<base->runnum <<" weight: "<<weight<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].FullLabel()<<" sysLabel:"<<sys.Label()<<"\n";
			}
		}
	  }
	  
	  // dummy data
	  // if(!addData && is_bkg && (title.find("QCD") == string::npos) && !sys)
	  //   FITBuilder.AddEvent(weight, Mperp, RISR,
	  // 			Categories[eindex], data_obs, sys);

	   //if(itree>0) break; //jump out of loop, dont duplicate entries on tree sys
//	   if(doTreeSys) break;//removed nsys loop so no break
//	}//end Nsys loop
      }
	//manage tree SYS ctr
//	if( itree>=1 && itree%2 ==0 ) treeSysCtr++;
      delete base;
      delete chain;
      if(ifile != -1)
        break;
//	}//end itree loop 
    }
  }

  FITBuilder.WriteFit(OutFile);

}
