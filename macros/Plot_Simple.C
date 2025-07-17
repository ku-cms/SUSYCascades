#include "PlottingTools.h"

void Plot_Simple(){

  Long64_t start = gSystem->Now();
  RestFrames::SetStyle();

  string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_Cascades_v2/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Simple_";

  g_Label = "PreSelection";
  
  output_root_file += g_Label;
  SanitizeString(output_root_file);
  folder_name = output_root_file;
  output_root_file += ".root";
  if(SavePDF){
    std::cout << "making dir for plots: plot_outputs/" << folder_name << std::endl;
    gSystem->Exec(("mkdir -p plot_outputs/"+folder_name).c_str());
    gSystem->Exec(("cp macros/Plot_Simple.C plot_outputs/"+folder_name+"/").c_str());
  }
  std::cout << "Saving plots to: " << output_root_file << std::endl;

  int SKIP = 1; // note that this only applies to BKG
  //double lumi = 138.; // Run 2
  //double lumi = 138.+109+27+34; // Run 2&3
  //double lumi = 9.451; // Summer23BPix
  double lumi = 400.; // estimated Run2+Run3
  double signal_boost = 1.;
  bool Norm = true; // scale hists by lumi

  std::vector<std::pair<std::string, ProcessList>> vec_samples;

  std::map<std::string, VS> map_vsignals;
  map_vsignals.insert(std::make_pair("Cascades_220", VS({"Cascades_*_220_*"})));
  map_vsignals.insert(std::make_pair("Cascades_260", VS({"Cascades_*_260_*"})));
  map_vsignals.insert(std::make_pair("Cascades_270", VS({"Cascades_*_270_*"})));
    
  // loop over backgrounds and add to map
  ProcessList backgrounds = ST.Get(kBkg);
  //backgrounds = backgrounds.Remove("QCD");
  //backgrounds = backgrounds.Filter("ttbar");
  for(int s = 0; s < int(backgrounds.GetN()); s++){
    ProcessList bkg;
    bkg += backgrounds[s];
    vec_samples.push_back(std::make_pair(backgrounds[s].Name(), bkg));
  }

  // loop over signals and add to map
  for(auto p = map_vsignals.begin(); p != map_vsignals.end(); p++){
    ProcessList signals;
    for(auto s : p->second){
      signals += ST.Get(s);
    }
    vec_samples.push_back(std::make_pair(p->first, signals));
  }

  // loop over data and add to map
  ProcessList data = ST.Get(kData);
  for(int s = 0; s < int(data.GetN()); s++){
    ProcessList datum;
    datum += data[s];
//    vec_samples.push_back(std::make_pair(data[s].Name(), datum));
  }

  // vectors to hold 'generic' plotting objects
  vector<vector<TH1*>*> hist_stacks; // vector for holding hist stacks
  vector<TH1*> hist_stack_MET; // example of vector for stacking all MET hists
  hist_stacks.push_back(&hist_stack_MET); // example of pushing back stack onto stacks

  // hists for holding number of events
  const int EC_bins = vec_samples.size() + 1;
  const int Zbi_bins = map_vsignals.size();
  int vec_samples_index = 0;
  int Zbi_samples_index = 0;
  TH2D* hist_EventCount = new TH2D("EventCount", "EventCount", 22, 0, 22, EC_bins, 0, EC_bins);
  hist_EventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_Zbi = new TH2D("Zbi", "Zbi", 22, 0, 22, Zbi_bins, 0, Zbi_bins);

  // Cut flows
  vector<TH1*> vect_hist_cutflow;
  const int CF_bins = 10;

  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    hist_EventCount->GetYaxis()->SetBinLabel(vec_samples_index+1, FP.getTitle(p->first).c_str());

    g_PlotTitle = p->first; 
    int Nsample = p->second.GetN();
    string title = p->first;

    // vectors to hold 'generic' plotting objects
    vector<TH1*> hists1;
    vector<TH2*> hists2;
    vector<TEfficiency*> effs;

    // Cut flows
    TH1D* hist_CutFlow = new TH1D((title+"_CutFlow").c_str(), (title+"_CutFlow").c_str(), CF_bins, 0, CF_bins);
    vect_hist_cutflow.push_back(hist_CutFlow);
    hists1.push_back(hist_CutFlow);
    int CF_bin = 0;

    // Declare hists here
    // push_back hists that you want to plot at the end (hists are filled regardless of whether or not you push_back)
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET [GeV]").c_str(), g_NX/2, 0., 1000.);
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot

    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 1000.);
    hists2.push_back(hist_RISR_PTISR);

    TEfficiency* eff_METtrig = new TEfficiency((title+"_eff_METtrig").c_str(), "Efficiency of MET trigger;Eff;MET [GeV]", g_NX, 0., 500.);
    effs.push_back(eff_METtrig);

    bool is_data = false;
    bool is_bkg = false;
    bool is_signal = false;

    for(int s = 0; s < Nsample; s++){
      Process proc = p->second[s];
      
      is_data   = (proc.Type() == kData);
      is_bkg    = (proc.Type() == kBkg);
      is_signal = (proc.Type() == kSig);
    
      if(is_signal)
        hist_Zbi->GetYaxis()->SetBinLabel(Zbi_samples_index+1, FP.getTitle(p->first).c_str());
      
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
        
        cout << "   Processing file " << file << " w/ tree " << tree << " for sample: " << p->first << endl;
        //cout << "      Sample weight is " << sample_weight << endl;
        if(is_FastSim)
          cout << "      Is FastSim" << endl;
        if(do_FilterDilepton)
          cout << "      Filter Out dilepton events" << endl;
        
        TChain* chain = ST.Tree(proc, f);
        
        ReducedBase_V2* base = new ReducedBase_V2(chain);
        
        int Nentry = base->fChain->GetEntries(); 
        if(SKIP < 1.) SKIP = 1.;
        int BKG_SKIP = SKIP;
        if(is_data || is_signal) BKG_SKIP = 1; // only use skip on BKG
        
        // event loop
        for(int e = 0; e < Nentry; e += BKG_SKIP){
          base->GetEntry(e);
          
          if((e/BKG_SKIP)%(std::max(1, int(Nentry/BKG_SKIP/10))) == 0)
            cout << "      event " << e << " | " << Nentry << endl;
          
          double weight = (base->weight != 0.) ? base->weight : 1.;
          if(!is_data && !is_signal)
            weight *= double(BKG_SKIP);

          // keep total in underflow no matter what cuts applied
          CF_bin = 0;
          hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight);
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "NTUPLES (>= 2L)");

          // Apply PreSelection
          int Njet = base->Njet;
          if(base->Njet == 0) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "Njet > 0");
          
          // get variables from root files using base class
          double MET = base->MET;
          double Mperp = base->Mperp;
          double RISR = base->RISR;
          double PTISR = base->PTISR;
          double PISR = base->PISR;

          double MSperpCM0 = base->MSperpCM0;
          double MQperpCM0 = base->MQperpCM0;
          double gammaPerpCM0 = base->gammaPerpCM0;
          double MJ = base->MJ;
          double ML = base->ML;
          double MLa = base->MLa;
          double MLb = base->MLb;
          double MSCM0 = base->MSCM0;
          double MQCM0 = base->MQCM0;
          double gammaCM0 = base->gammaCM0;

          if(MET < 150.) // PreSelection
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "MET > 150");
          
          // apply trigger to data and FullSim events
          if(!base->METORtrigger && !is_FastSim) // PreSelection
          //if(!base->SingleElectrontrigger && !base->SingleMuontrigger && !is_FastSim)
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "MET trigger");

          if(PTISR < 200.) // PreSelection
	    continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "PTISR > 200");
            
          if(RISR < 0.5 || RISR > 1.0) // PreSelection
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "R_{ISR} > 0.5");

          // Cleaning cuts...
          double dphiCMI = base->dphiCMI;
          double PTCM = base->PTCM;
          double x = fabs(dphiCMI);
          
          if(PTCM > 200.)
            continue;
          if(PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
             -2.777*x*x+1.388*x+0.8264 > 0.)
            continue;
          if(PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
             -1.5625*x*x+7.8125*x-8.766 > 0.)
            continue;
          // End of Cleaning cuts...
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "Cleaning Cuts");
            
          double dphiMET_V = base->dphiMET_V;
          if(fabs(base->dphiMET_V) > acos(-1.)/2.) // PreSelection
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "#Delta #phi(MET, V) < #frac{#pi}{2}");

          // Get Physics Objects
          int Nlep     = base->Nlep;
          int NjetS    = base->Njet_S;
          int NbjetS   = base->Nbjet_S;
          int NjetISR  = base->Njet_ISR;
          int Njet_a   = base->Njet_a;
          int Njet_b   = base->Njet_b;
          int Nlep_a   = base->Nlep_a;
          int Nlep_b   = base->Nlep_b;
          int NbjetISR = base->Nbjet_ISR;

          //if(NbjetISR + NbjetS > 1) continue; // B-Veto
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "B-Veto");

          LepList list_a;
          LepList list_b;
          LepList list_leps;  
          int index;
            
          // ID leps as gold, silver, or bronze
          for(int i = 0; i < Nlep_a; i++){
            index = (*base->index_lep_a)[i];
              
            int PDGID = base->PDGID_lep->at(index);
              
            int qual = base->LepQual_lep->at(index);
            LepID id = kBronze;
            if(qual == 0) id = kGold;
            if(qual == 1) id = kSilver;
            LepFlavor flavor;
            if(abs(PDGID) == 11)
              flavor = kElec;
            else
              flavor = kMuon;
            LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
            LepSource source = LepSource(base->SourceID_lep->at(index));
              
            list_a += Lep(flavor, charge, id, source);
            list_leps += Lep(flavor, charge, id, source);
          }
          for(int i = 0; i < Nlep_b; i++){
            index = (*base->index_lep_b)[i];
            
            int PDGID = base->PDGID_lep->at(index);

            int qual = base->LepQual_lep->at(index);
            LepID id = kBronze;
            if(qual == 0) id = kGold;
            if(qual == 1) id = kSilver;
            LepFlavor flavor;
            if(abs(PDGID) == 11)
              flavor = kElec;
            else
              flavor = kMuon;
            LepCharge charge = (base->Charge_lep->at(index) > 0 ? kPos : kNeg);
            LepSource source = LepSource(base->SourceID_lep->at(index));
            
            list_b += Lep(flavor, charge, id, source);
            list_leps += Lep(flavor, charge, id, source);
          }

          // cut on lepton quality
          bool skip = false;
          int nBL = 0; // number of bronze leps
          int nSL = 0; // number of silver leps
          int nGL = 0; // number of gold leps
          for(int i = 0; i < list_leps.GetN(); i++){
            if(list_leps[i].ID() == kBronze) { nBL++; } // 'Bronze'
            if(list_leps[i].ID() == kSilver) { nSL++; } // 'Silver'
            if(list_leps[i].ID() == kGold) { nGL++; } // 'Gold'
          }
          //if(nBL > 0) continue; // no bronze leps
          //if(nSL > 0) continue; // no silver leps
          //if(nGL > 0) continue; // no gold leps
          //if(nSL > 3) continue; // no more than N silver leps
          // ATLAS
          //if(base->PT_lep->at(0) < 28.) continue;
          //if(base->PT_lep->at(1) < 20.) continue;
          //if(Nlep > 2)
          //  if(base->PT_lep->at(2) < 10.) continue;

          // Fill hists, effs, etc.
          hist_MET->Fill(MET, weight);
          hist_RISR_PTISR->Fill(RISR, PTISR, weight);
          eff_METtrig->Fill(base->METtrigger, MET, weight);

          // Event Counting
          int EC_X = 0; // root hists have underflow in bin 0
          int EC_Y = vec_samples_index+1; // root hists have underflow in bin 0
          int Zbi_EC_Y = Zbi_samples_index+1;
          int cat_Nleps = list_leps.GetN();
          if (cat_Nleps > 4) cat_Nleps = 4; // ge4L is upper limit
          //if (cat_Nleps > 5) cat_Nleps = 5; // ge5L is upper limit
          
          // Extract total lepton charge
          std::vector<LepFlavor> flavors;
          std::vector<LepCharge> charges;
          int total_charge = 0;
          for (int i = 0; i < cat_Nleps; i++){
            flavors.push_back(list_leps[i].Flavor());
            charges.push_back(list_leps[i].Charge());
            total_charge += (charges.back() == LepCharge::kPos) ? 1 : -1;
          }
          std::sort(flavors.begin(), flavors.end());
          int num_e = std::count(flavors.begin(), flavors.end(), kElec);
          int abs_charge = abs(total_charge);
          
          if (cat_Nleps == 2) {
            bool same_flavor = (flavors[0] == flavors[1]);
            bool opposite_sign = (charges[0] != charges[1]);

            if (!same_flavor && !opposite_sign)     EC_X = 1; //flavor_category = "SSOF";
            else if (same_flavor && !opposite_sign) EC_X = 2; //flavor_category = "SSSF";
            else if (!same_flavor && opposite_sign) EC_X = 3; //flavor_category = "OSOF";
            else                                    EC_X = 4; //flavor_category = "OSSF";
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
          }
          else if(cat_Nleps == 3) {
            EC_X = 5;
            EC_X += num_e;
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            EC_X = 9;
            EC_X += abs_charge;
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
          }
          else {//if(cat_Nleps == 4){
            EC_X = 13;
            EC_X += num_e;
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            EC_X = 18;
            EC_X += abs_charge;
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
          }
          //else {
          //  EC_X = 23;
          //  hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
          //  if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
          //  if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
          //  if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
          //}
          hist_EventCount->SetBinContent(0,EC_Y,hist_EventCount->GetBinContent(0,EC_Y)+weight); // normalized to selection
          if(is_bkg){ // total SM bkg
            hist_EventCount->SetBinContent(0,EC_bins,hist_EventCount->GetBinContent(0,EC_bins)+weight); // normalized to selection
          }

        }
        delete base;
        delete chain;
        std::cout << "Integral of file: " << file << " is " << hist_MET->Integral() << std::endl;
      } // for(int f = 0; f < Nfile; f++)

    } // for(int s = 0; s < Nsample; s++)

    // call plotting functions after looping over files in samples
    for(int hist1 = 0; hist1 < int(hists1.size()); hist1++){
      Plot_Hist(hists1[hist1], Norm, lumi, signal_boost, is_signal);
    }
    for(int hist2 = 0; hist2 < int(hists2.size()); hist2++){
      Plot_Hist(hists2[hist2], Norm, lumi, signal_boost, is_signal);
    }
    for(int eff = 0; eff < int(effs.size()); eff++)
      Plot_Eff(effs[eff]);
    hists1.clear();
    hists2.clear();
    effs.clear();
    vec_samples_index++;
    if(is_signal)
      Zbi_samples_index++;

  } // for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){

  int index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    std::cout << "Integral of proc: " << p->first << " is " << (*hist_stacks[0])[index]->Integral() << std::endl;
    index++;
  }

  for(int stack_h = 0; stack_h < int(hist_stacks.size()); stack_h++){
     g_PlotTitle = (*hist_stacks[stack_h])[0]->GetName();
     // loop through process names, if process name in title, remove that part from title and add '_stack'
     for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
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
     Plot_Stack(*hist_stacks[stack_h], vec_samples, signal_boost);
  }

  Plot_EventCount((TH2*)hist_EventCount->Clone("EventCount_Scaled"), true, lumi, false, 0., false, false);
  Plot_EventCount((TH2*)hist_EventCount->Clone("EventCount_SoBAllBKG"), true, lumi, false, 0., true, true);
  Plot_EventCount(hist_EventCount, false, lumi, false, 0., false, false);
  Plot_EventCount((TH2*)hist_Zbi->Clone("EventCount_SoB"), true, lumi, false, 0., true, false);
  Plot_EventCount(hist_Zbi, true, lumi, true, 0.2, false, false);

  Plot_CutFlow(vect_hist_cutflow, true, lumi, signal_boost, vec_samples);

  Long64_t end = gSystem->Now();
  std::cout << "Time to process " << (end-start)/1000.0 << " seconds" << std::endl;
  gApplication->Terminate(0);
} // End of macro
