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
#include "METTriggerTool.hh"

using namespace std;

int main(int argc, char* argv[]) {
  int ifile = -1;
  //string NtuplePath = "root://xrootd.unl.edu//store/user/zflowers/crogan/";
  string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v0/";
  string OutFile    = "BuildFitInput_output.root";

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
 
  bool debugVerbosity = false;
 
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
  }
  if(addData){
    cout << "Adding all data for year " << year << endl;
    samples += ST.Get(kData);
  }
  for(int p = 0; p < int(proc_to_add.size()); p++){
    cout << "Adding processes that match \"" << proc_to_add[p] << "\"" << endl;
    samples += ST.Get(proc_to_add[p]);
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
  m_METTriggerTool.BuildMap("Parameters.csv");
  //m_METTriggerTool.BuildMap("csv/METTrigger/Parameters.csv");

  ScaleFactorTool SF;
  SF.AddBtagFolder("./BtagSF");

  Systematics systematics(1);
  if(doSys)
    systematics += SYS.GetWeightSystematics();

  FitInputBuilder FITBuilder(extrahist);

  // dummy in case there is no data requested
  Process data_obs("data_obs", kData);

  // sample (process) loop
  int Nsample = samples.GetN();
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
    cout << "processing sample " << proc.Name() << endl;
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

    int Nsys = (is_data ? 1 : systematics.GetN());
    
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
	int NSV      = base->NSV_S;

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
	//debugging assume index problem has been fixed.. this will not work with older ntuples (applied in A and B)
	  if(base->ID_lep->at(index*2) < 3 ||
	//    if(base->ID_lep->at(index) < 3 ||
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
	//  LepSource source = LepSource(base->SourceID_lep->at(index));
	  LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version (this is the one turned on in master)
	  
	
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
	  if(base->ID_lep->at(index*2) < 3 || //index fixed for newly produced ntuples
	 //   if(base->ID_lep->at(index) < 3 ||
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
	//  LepSource source = LepSource(base->SourceID_lep->at(index));
	  LepSource source = LepSource(base->ID_lep->at(index*2+1)); // fix for current ntuple version
	  list_b += Lep(flavor, charge, id, source);
	
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

	double PTISR = base->PTISR;
	double PTISR_to_HT = PTISR*2.-150.;
	double weight = 1.;
	double btag_weight = 1.;
	double PU_weight = 1.;
	double trig_weight = 1.;

	if(!is_data){
	  weight = (setLumi ? lumi : ST.Lumi())*base->weight*sample_weight;
	}
	// systematics loop
	// do down sys first
	string correct_sys = "";
	for(int is = 0; is < Nsys; is++){
	  Systematic& sys = systematics[is];
//	  std::cout<<"systemtatics list "<< systematics[is].Label()<<"\n";
	  if(!(!sys)){
	    if(sys.IsUp()){
	      sys.Down();
	      is--;
	    } else {
	      sys.Up();
	    }
	  }
	  

	  btag_weight = 1.;
	  PU_weight = 1.;
	  trig_weight = 1.;
          if(!(!sys) && is_data) continue;      

		
	    trig_weight = m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0);
            if(is_FastSim)
	      trig_weight = m_METTriggerTool.Get_EFF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0)*
		m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 0);



            if(sys.Label().find("MET_TRIG") != std::string::npos && proc.Name() != "QCD")
            {
             /* 
              if(sys.Label() != correct_sys && correct_sys != "") continue;
              string scat = Categories[eindex].FullLabel();
            
              if(sys.Label() == "MET_TRIG_el")
            	if((scat.find("0L") == std::string::npos) && 
                       (scat.find("1L") == std::string::npos || (scat.find("_el") == std::string::npos && scat.find("_lp") == std::string::npos && scat.find("_lm") == std::string::npos)) &&
                       ((scat.find("2L") == std::string::npos && scat.find("3L") == std::string::npos) || (scat.find("elel") == std::string::npos && scat.find("_ll") == std::string::npos && 
                       scat.find("_noZ") == std::string::npos && scat.find("_Zstar") == std::string::npos && scat.find("_SS") == std::string::npos)))
            		continue;
            
              if(sys.Label() == "MET_TRIG_mu")
            	if((scat.find("0L") != std::string::npos) || 
                       (scat.find("1L") == std::string::npos || (scat.find("_mu") == std::string::npos && scat.find("_lp") == std::string::npos && scat.find("_lm") == std::string::npos)) &&
                       ((scat.find("2L") == std::string::npos && scat.find("3L") == std::string::npos) || (scat.find("mumu") == std::string::npos && scat.find("elmu") == std::string::npos && scat.find("_ll") == std::string::npos && 
                       scat.find("_noZ") == std::string::npos && scat.find("_Zstar") == std::string::npos && scat.find("_SS") == std::string::npos)))
            		continue;
              */
            	        if(sys.IsUp())
            	          if(is_FastSim)
            	            trig_weight = m_METTriggerTool.Get_EFF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1)*
            	              m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1);
            	          else
            	            trig_weight = m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, 1);
            	        else
            	          if(is_FastSim)
            	            trig_weight = m_METTriggerTool.Get_EFF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1)*
            	              m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1);
            	          else
            	            trig_weight = m_METTriggerTool.Get_SF(base->MET, PTISR_to_HT, year, (base->Nele > 0), (base->Nmu > 0), false, -1);
            
            }

	      
	  //
	  // BTAG systematics from the ntuples
	  //
	  /*
	    if(sys == Systematic("BTAGHF_SF"))
	    if(sys.IsUp())
	    btag_weight *= base->BtagSFweight_up;
	    else
	    btag_weight *= base->BtagSFweight_down;
	    else 
	    btag_weight *= base->BtagSFweight;

	    if(sys == Systematic("BTAGLF_SF"))
	    if(sys.IsUp())
	    btag_weight *= base->BtagSFweight_up;
	    else
	    btag_weight *= base->BtagSFweight_down;
	    else 
	    btag_weight *= base->BtagSFweight;
	  */

	  //
	  // BTAG systematics on the fly (needs jet collection in reduced ntuples)
	  //
	  
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


	  // turn off PU systematics for now
	  // if(sys == Systematic("PU_SF"))
	  //   if(sys.IsUp())
	  // 	PU_weight = base->PUweight_up;
	  //   else
	  // 	PU_weight = base->PUweight_down;
	  // else
	  //   PU_weight = base->PUweight;





	  weight *= btag_weight*PU_weight;
	  //weight *= btag_weight*PU_weight*trig_weight;
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
	      
	      FITBuilder.AddEvent(weight/double(Nf), Mperp, RISR,
				  Categories[eindex], proc.FakeProcess(flabels[fl]), sys);
		if(debugVerbosity){
			std::cout<<"Adding fakes event:"<<e<<" weight: "<<weight/double(Nf)<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].Label()<<"  flabel:"<<flabels[fl]<<"\n";
		}
	    }
	  } else {
		//std::cout<<"adding event "<< weight <<" "<< Mperp <<" "<< RISR <<" "<<Categories[eindex].Label()<<" "<<proc.Name()<<" "<<sys.Label()<<"\n";
	    FITBuilder.AddEvent(weight, Mperp, RISR,
				Categories[eindex], proc, sys);
		if(debugVerbosity){
			std::cout<<"Adding event:"<<e<<" weight: "<<weight<<" Mperp:"<<Mperp<<" RISR:"<<RISR<<" gammaT:"<<gammaT<<" PTISR:"<<PTISR<<" Cat:"<<Categories[eindex].Label()<<" sysLabel:"<<sys.Label()<<"\n";
		}
	  }
	  
	  // dummy data
	  // if(!addData && is_bkg && (title.find("QCD") == string::npos) && !sys)
	  //   FITBuilder.AddEvent(weight, Mperp, RISR,
	  // 			Categories[eindex], data_obs, sys);
	}
      }
      delete base;
      delete chain;
      if(ifile != -1)
        break;
    }
  }

  FITBuilder.WriteFit(OutFile);

}
