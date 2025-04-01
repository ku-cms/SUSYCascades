#include "PlottingTools.h"

void Plot_Advanced(){

  std::cout << "Saving plots to: " << output_root_file << std::endl;
  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";
  //string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Advanced";
  //g_Label = "PreSelection";
  //g_Label = "2L gold 0J SR";
  //g_Label = "2 lepton SR";
  //g_Label = "2 lepton ttbar CR";
  output_root_file += "_"+g_Label+".root";
  // Replaces spaces in name of output file with _
  std::replace(output_root_file.begin(), output_root_file.end(), ' ', '_');

  int SKIP = 10000; // note that this only applies to BKG
  //double lumi = 138.; // Run 2
  //double lumi = 138.+109+27+34; // Run 2&3
  double lumi = 9.451; // Summer23BPix
  //double lumi = 400.;
  double signal_boost = 1000.;
  bool Norm = true; // scale hists by lumi

  std::vector<std::pair<std::string, ProcessList>> vec_samples;

  std::map<std::string, VS> map_vsignals;
  map_vsignals.insert(std::make_pair("Cascades_220", VS({"Cascades_*_220_*"})));
  map_vsignals.insert(std::make_pair("Cascades_260", VS({"Cascades_*_260_*"})));
  map_vsignals.insert(std::make_pair("Cascades_270", VS({"Cascades_*_270_*"})));
   
  // loop over signals and add to map
  for(auto p = map_vsignals.begin(); p != map_vsignals.end(); p++){
    ProcessList signals;
    for(auto s : p->second){
      signals += ST.Get(s);
    }
//    vec_samples.push_back(std::make_pair(p->first, signals));
  }
  
  // loop over backgrounds and add to map
  ProcessList backgrounds = ST.Get(kBkg);
  //backgrounds = backgrounds.Remove("QCD");
  for(int s = 0; s < int(backgrounds.GetN()); s++){
    ProcessList bkg;
    bkg += backgrounds[s];
    vec_samples.push_back(std::make_pair(backgrounds[s].Name(), bkg));
  }

  // loop over data and add to map
  ProcessList data = ST.Get(kData);
  for(int s = 0; s < int(data.GetN()); s++){
    ProcessList datum;
    datum += data[s];
    vec_samples.push_back(std::make_pair(data[s].Name(), datum));
  }

  // vectors to hold 'generic' plotting objects
  vector<vector<TH1*>*> hist_stacks; // vector for holding hist stacks
  vector<TH1*> hist_stack_MET; // example of vector for stacking all MET hists
  hist_stacks.push_back(&hist_stack_MET);

  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    // vectors to hold 'generic' plotting objects
    vector<TH1*> hists1;
    vector<TH2*> hists2;
    vector<TEfficiency*> effs;

    g_PlotTitle = p->first; 
    int Nsample = p->second.GetN();
    string title = p->first;

    // Declare hists here
    // push_back hists that you want to plot at the end (hists are filled regardless of whether or not you push_back)
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET [GeV]").c_str(), g_NX/4, 100., 500.);
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot
    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 800.);
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
        int BKG_SKIP = SKIP;
        if(is_data || is_signal) BKG_SKIP = 1; // only use skip on BKG
        
        // event loop
        for(int e = 0; e < Nentry; e += BKG_SKIP){
          base->GetEntry(e);
          
          if((e/BKG_SKIP)%(std::max(1, int(Nentry/BKG_SKIP/10))) == 0)
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
          // End of Cleaning cuts...
            
          if(fabs(base->dphiMET_V) > acos(-1.)/2.)
            continue;
            
          //if(base->RISR < 0.5 || base->RISR > 1.0)
          if(base->RISR < 0.4 || base->RISR > 0.7) // CR
          //if(base->RISR < 0.9)
            continue;

          // Get Physics Objects
          int Nlep     = base->Nlep;
          int NjetS    = base->Njet_S;
          int NbjetS   = base->Nbjet_S;
          int NjetISR  = base->Njet_ISR;
          int NbjetISR = base->Nbjet_ISR;

          if(NbjetISR + NbjetS != 2) continue; // CR
          //if(NbjetISR + NbjetS > 1) continue; // SR

          if(Nlep != 2) continue;
          //if(NjetS != 0) continue; //SR

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

          // cut on lepton quality
          bool skip = false;
          for(int i = 0; i < list_leps.GetN(); i++)
            if(list_leps[i].ID() != kGold) skip = true; // SR
          //if(skip) continue; 

          // get variables from root files using base class
          double Mperp = base->Mperp;
          double RISR = base->RISR;
          double PTISR = base->PTISR;
          
          double weight = (base->weight != 0.) ? base->weight : 1.;
          if(!is_data && !is_signal)
            weight *= double(BKG_SKIP);

          // Fill hists, effs, etc.
          hist_MET->Fill(MET, weight);
          hist_RISR_PTISR->Fill(RISR, PTISR, weight);
          eff_METtrig->Fill(base->METtrigger, MET, weight);
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
  gApplication->Terminate(0);
} // End of macro
