#include "PlottingTools.h"

void Plot_Advanced(){

  Long64_t start = gSystem->Now();
  RestFrames::SetStyle();
  //InitRJRtree();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";
  //string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Advanced_Cands_";

  bool CandSameHemi = true;
  if(CandSameHemi) g_Label += "SameHemi_";

  //g_Label = "TESTING";
  g_Label = "PreSelection";
  //g_Label = "PreSelection & Gold";
  //g_Label = "PreSelection & 4 Gold #mu";
  //g_Label = "No Cuts";
  //g_Label = "ATLAS Cuts";
  //g_Label = "MET > 150";
  
  output_root_file += g_Label;
  SanitizeString(output_root_file);
  folder_name = output_root_file;
  output_root_file += ".root";
  if(SavePDF){
    std::cout << "making dir for plots: " << folder_name << std::endl;
    gSystem->Exec(("mkdir -p "+folder_name).c_str());
    gSystem->Exec(("cp macros/Plot_Advanced.C "+folder_name+"/").c_str());
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
  vector<TH1*> hist_stack_RISR;
  hist_stacks.push_back(&hist_stack_RISR);
  vector<TH1*> hist_stack_PTISR;
  hist_stacks.push_back(&hist_stack_PTISR);
  vector<TH1*> hist_stack_gammaPerp;
  hist_stacks.push_back(&hist_stack_gammaPerp);
  vector<TH1*> hist_stack_MQperp;
  hist_stacks.push_back(&hist_stack_MQperp);
  vector<TH1*> hist_stack_MSperpCM0;
  hist_stacks.push_back(&hist_stack_MSperpCM0);
  vector<TH1*> hist_stack_MQperpCM0;
  hist_stacks.push_back(&hist_stack_MQperpCM0);
  vector<TH1*> hist_stack_gammaPerpCM0;
  hist_stacks.push_back(&hist_stack_gammaPerpCM0);
  vector<TH1*> hist_stack_MSCM0;
  hist_stacks.push_back(&hist_stack_MSCM0);
  vector<TH1*> hist_stack_MQCM0;
  hist_stacks.push_back(&hist_stack_MQCM0);
  vector<TH1*> hist_stack_gammaCM0;
  hist_stacks.push_back(&hist_stack_gammaCM0);
  vector<TH1*> hist_stack_ML;
  hist_stacks.push_back(&hist_stack_ML);
  vector<TH1*> hist_stack_MJ;
  hist_stacks.push_back(&hist_stack_MJ);
  vector<TH1*> hist_stack_CandML;
  hist_stacks.push_back(&hist_stack_CandML);
  vector<TH1*> hist_stack_CandBeta;
  hist_stacks.push_back(&hist_stack_CandBeta);
  vector<TH1*> hist_stack_CandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_CandDeltaPhiMET);
  vector<TH1*> hist_stack_CandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_CandCosDecayAngle);
  vector<TH1*> hist_stack_RZPara;
  hist_stacks.push_back(&hist_stack_RZPara);
  vector<TH1*> hist_stack_RZPerp;
  hist_stacks.push_back(&hist_stack_RZPerp);
  vector<TH1*> hist_stack_RZAng;
  hist_stacks.push_back(&hist_stack_RZAng);
  vector<TH1*> hist_stack_PZPara;
  hist_stacks.push_back(&hist_stack_PZPara);
  vector<TH1*> hist_stack_PZPerp;
  hist_stacks.push_back(&hist_stack_PZPerp);
  vector<TH1*> hist_stack_MatchedCandML;
  hist_stacks.push_back(&hist_stack_MatchedCandML);
  vector<TH1*> hist_stack_MatchedCandBeta;
  hist_stacks.push_back(&hist_stack_MatchedCandBeta);
  vector<TH1*> hist_stack_MatchedCandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_MatchedCandDeltaPhiMET);
  vector<TH1*> hist_stack_MatchedCandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_MatchedCandCosDecayAngle);
  vector<TH1*> hist_stack_UnmatchedCandML;
  hist_stacks.push_back(&hist_stack_UnmatchedCandML);
  vector<TH1*> hist_stack_UnmatchedCandBeta;
  hist_stacks.push_back(&hist_stack_UnmatchedCandBeta);
  vector<TH1*> hist_stack_UnmatchedCandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_UnmatchedCandDeltaPhiMET);
  vector<TH1*> hist_stack_UnmatchedCandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_UnmatchedCandCosDecayAngle);

  // hists for holding number of events
  const int EC_bins = vec_samples.size() + 1;
  const int Zbi_bins = map_vsignals.size();
  int vec_samples_index = 0;
  int Zbi_samples_index = 0;
  TH2D* hist_EventCount = new TH2D("EventCount", "EventCount", 22, 0, 22, EC_bins, 0, EC_bins);
  hist_EventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_Zbi = new TH2D("Zbi", "Zbi", 22, 0, 22, Zbi_bins, 0, Zbi_bins);

  // Event Counting by cands
  TH2D* hist_CandsEventCount = new TH2D("CandsEventCount", "CandsEventCount", 23, 0, 23, EC_bins, 0, EC_bins);
  hist_CandsEventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_CandsZbi = new TH2D("CandsZbi", "CandsZbi", 23, 0, 23, Zbi_bins, 0, Zbi_bins);

  // Cut flows
  vector<TH1*> vect_hist_cutflow;
  const int CF_bins = 10;

  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    hist_EventCount->GetYaxis()->SetBinLabel(vec_samples_index+1, FP.getTitle(p->first).c_str());
    hist_CandsEventCount->GetYaxis()->SetBinLabel(vec_samples_index+1, FP.getTitle(p->first).c_str());

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
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET [GeV]").c_str(), g_NX/2, 100., 1000.);
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot
    TH1D* hist_RISR = new TH1D((title+"_RISR").c_str(), (title+"_RISR;R_{ISR}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_RISR);
    hist_stack_RISR.push_back(hist_RISR);    
    TH1D* hist_PTISR = new TH1D((title+"_PTISR").c_str(), (title+"_PTISR;p_{T}^{ISR} [GeV]").c_str(), g_NX/2, 0., 1000.);
    hists1.push_back(hist_PTISR);
    hist_stack_PTISR.push_back(hist_PTISR);
    TH1D* hist_gammaPerp = new TH1D((title+"_gammaPerp").c_str(), (title+"_gammaPerp;#gamma_{#perp}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_gammaPerp);
    hist_stack_gammaPerp.push_back(hist_gammaPerp);
    TH1D* hist_MQperp = new TH1D((title+"_MQperp").c_str(), (title+"_MQperp;M_{#perp}").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_MQperp);
    hist_stack_MQperp.push_back(hist_MQperp);
    TH1D* hist_MSperpCM0 = new TH1D((title+"_MSperpCM0").c_str(), (title+"_MSperpCM0;MS_{#perp CM0}").c_str(), g_NX/2, 0., 600.);
    hists1.push_back(hist_MSperpCM0);
    hist_stack_MSperpCM0.push_back(hist_MSperpCM0);
    TH1D* hist_MQperpCM0 = new TH1D((title+"_MQperpCM0").c_str(), (title+"_MQperpCM0;MQ_{#perp CM0}").c_str(), g_NX/2, 0., 200.);
    hists1.push_back(hist_MQperpCM0);
    hist_stack_MQperpCM0.push_back(hist_MQperpCM0);
    TH1D* hist_gammaPerpCM0 = new TH1D((title+"_gammaPerpCM0").c_str(), (title+"_gammaPerpCM0;#gamma_{#perp CM0}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_gammaPerpCM0);
    hist_stack_gammaPerpCM0.push_back(hist_gammaPerpCM0);
    TH1D* hist_MSCM0 = new TH1D((title+"_MSCM0").c_str(), (title+"_MSCM0;MS_{CM0}").c_str(), g_NX/2, 0., 600.);
    hists1.push_back(hist_MSCM0);
    hist_stack_MSCM0.push_back(hist_MSCM0);
    TH1D* hist_MQCM0 = new TH1D((title+"_MQCM0").c_str(), (title+"_MQCM0;MQ_{CM0}").c_str(), g_NX/2, 0., 200.);
    hists1.push_back(hist_MQCM0);
    hist_stack_MQCM0.push_back(hist_MQCM0);
    TH1D* hist_gammaCM0 = new TH1D((title+"_gammaCM0").c_str(), (title+"_gammaCM0;#gamma_{CM0}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_gammaCM0);
    hist_stack_gammaCM0.push_back(hist_gammaCM0);
    TH1D* hist_ML = new TH1D((title+"_ML").c_str(), (title+"_ML;ML").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_ML);
    hist_stack_ML.push_back(hist_ML);
    TH1D* hist_MJ = new TH1D((title+"_MJ").c_str(), (title+"_MJ;MJ").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_MJ);
    hist_stack_MJ.push_back(hist_MJ);

    TH1D* hist_CandML = new TH1D((title+"_CandML").c_str(), (title+"_CandML;CandML").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_CandML);
    hist_stack_CandML.push_back(hist_CandML);
    TH1D* hist_CandBeta = new TH1D((title+"_CandBeta").c_str(), (title+"_CandBeta;CandBeta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_CandBeta);
    hist_stack_CandBeta.push_back(hist_CandBeta);
    TH1D* hist_CandDeltaPhiMET = new TH1D((title+"_CandDeltaPhiMET").c_str(), (title+"_CandDeltaPhiMET;CandDeltaPhiMET").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_CandDeltaPhiMET);
    hist_stack_CandDeltaPhiMET.push_back(hist_CandDeltaPhiMET);
    TH1D* hist_CandCosDecayAngle = new TH1D((title+"_CandCosDecayAngle").c_str(), (title+"_CandCosDecayAngle;CandCosDecayAngle").c_str(), g_NX/2, -1., 1.);
    hists1.push_back(hist_CandCosDecayAngle);
    hist_stack_CandCosDecayAngle.push_back(hist_CandCosDecayAngle);
    TH1D* hist_RZPara = new TH1D((title+"_RZPara").c_str(), (title+"_RZPara;RZPara").c_str(), g_NX/2, -1.5, 1.5);
    hists1.push_back(hist_RZPara);
    hist_stack_RZPara.push_back(hist_RZPara);
    TH1D* hist_RZPerp = new TH1D((title+"_RZPerp").c_str(), (title+"_RZPerp;RZPerp").c_str(), g_NX/2, 0., 2.);
    hists1.push_back(hist_RZPerp);
    hist_stack_RZPerp.push_back(hist_RZPerp);
    TH1D* hist_RZAng = new TH1D((title+"_RZAng").c_str(), (title+"_RZAng;RZAng").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_RZAng);
    hist_stack_RZAng.push_back(hist_RZAng);
    TH1D* hist_PZPara = new TH1D((title+"_PZPara").c_str(), (title+"_PZPara;PZPara").c_str(), g_NX/2, 0., 1000.);
    hists1.push_back(hist_PZPara);
    hist_stack_PZPara.push_back(hist_PZPara);
    TH1D* hist_PZPerp = new TH1D((title+"_PZPerp").c_str(), (title+"_PZPerp;PZPerp").c_str(), g_NX/2, 0., 1000.);
    hists1.push_back(hist_PZPerp);
    hist_stack_PZPerp.push_back(hist_PZPerp);

    TH1D* hist_MatchedCandML = new TH1D((title+"_MatchedCandML").c_str(), (title+"_MatchedCandML;MatchedCandML").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_MatchedCandML);
    hist_stack_MatchedCandML.push_back(hist_MatchedCandML);
    TH1D* hist_MatchedCandBeta = new TH1D((title+"_MatchedCandBeta").c_str(), (title+"_MatchedCandBeta;MatchedCandBeta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_MatchedCandBeta);
    hist_stack_MatchedCandBeta.push_back(hist_MatchedCandBeta);
    TH1D* hist_MatchedCandDeltaPhiMET = new TH1D((title+"_MatchedCandDeltaPhiMET").c_str(), (title+"_MatchedCandDeltaPhiMET;MatchedCandDeltaPhiMET").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_MatchedCandDeltaPhiMET);
    hist_stack_MatchedCandDeltaPhiMET.push_back(hist_MatchedCandDeltaPhiMET);
    TH1D* hist_MatchedCandCosDecayAngle = new TH1D((title+"_MatchedCandCosDecayAngle").c_str(), (title+"_MatchedCandCosDecayAngle;MatchedCandCosDecayAngle").c_str(), g_NX/2, -1., 1.);
    hists1.push_back(hist_MatchedCandCosDecayAngle);
    hist_stack_MatchedCandCosDecayAngle.push_back(hist_MatchedCandCosDecayAngle);

    TH1D* hist_UnmatchedCandML = new TH1D((title+"_UnmatchedCandML").c_str(), (title+"_UnmatchedCandML;UnmatchedCandML").c_str(), g_NX/2, 0., 150.);
    hists1.push_back(hist_UnmatchedCandML);
    hist_stack_UnmatchedCandML.push_back(hist_UnmatchedCandML);
    TH1D* hist_UnmatchedCandBeta = new TH1D((title+"_UnmatchedCandBeta").c_str(), (title+"_UnmatchedCandBeta;UnmatchedCandBeta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_UnmatchedCandBeta);
    hist_stack_UnmatchedCandBeta.push_back(hist_UnmatchedCandBeta);
    TH1D* hist_UnmatchedCandDeltaPhiMET = new TH1D((title+"_UnmatchedCandDeltaPhiMET").c_str(), (title+"_UnmatchedCandDeltaPhiMET;UnmatchedCandDeltaPhiMET").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_UnmatchedCandDeltaPhiMET);
    hist_stack_UnmatchedCandDeltaPhiMET.push_back(hist_UnmatchedCandDeltaPhiMET);
    TH1D* hist_UnmatchedCandCosDecayAngle = new TH1D((title+"_UnmatchedCandCosDecayAngle").c_str(), (title+"_UnmatchedCandCosDecayAngle;UnmatchedCandCosDecayAngle").c_str(), g_NX/2, -1., 1.);
    hists1.push_back(hist_UnmatchedCandCosDecayAngle);
    hist_stack_UnmatchedCandCosDecayAngle.push_back(hist_UnmatchedCandCosDecayAngle);

    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 1000.);
    hists2.push_back(hist_RISR_PTISR);
    TH2D* hist_RISR_Mperp = new TH2D((title+"_RISR_Mperp").c_str(), (title+"_RISR_Mperp;R_{ISR};M_{#perp} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 100.);
    hists2.push_back(hist_RISR_Mperp);
    TH2D* hist_dphiCMI_PTCM = new TH2D((title+"_dphiCMI_PTCM").c_str(), (title+"_dphiCMI_PTCM;#Delta #phi_{(CM,I)};p_{T}^{CM}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 0., 500.);
    hists2.push_back(hist_dphiCMI_PTCM);
    TH2D* hist_dphiMETV_PTISR = new TH2D((title+"_dphiMETV_PTISR").c_str(), (title+"_dphiMETV_PTISR;#Delta #phi_{(I,V)};p_{T}^{ISR}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 200., 800.);
    hists2.push_back(hist_dphiMETV_PTISR);
    TH2D* hist_gammaPerp_RISR = new TH2D((title+"_gammaPerp_RISR").c_str(), (title+"_gammaPerp_RISR;#gamma_{#perp};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_gammaPerp_RISR);
    TH2D* hist_RISR_mllLEAD = new TH2D((title+"_RISR_mllLEAD").c_str(), (title+"_RISR_mllLEAD;R_{ISR};m_{ll}LEAD").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.);
    hists2.push_back(hist_RISR_mllLEAD);
    TH2D* hist_RISR_mL = new TH2D((title+"_RISR_mL").c_str(), (title+"_RISR_mL;R_{ISR};mL").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 150.); // mass of leptonic system
    hists2.push_back(hist_RISR_mL);

    TH2D* hist_MSperpCM0_RISR = new TH2D((title+"_MSperpCM0_RISR").c_str(), (title+"_MSperpCM0_RISR;MS_{#perp CM0};RISR").c_str(), g_NX/2., 0., 500., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_MSperpCM0_RISR);
    TH2D* hist_MQperpCM0_RISR = new TH2D((title+"_MQperpCM0_RISR").c_str(), (title+"_MQperpCM0_RISR;MQ_{#perp CM0};RISR").c_str(), g_NX/2., 0., 300., g_NX/2., 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
    hists2.push_back(hist_MQperpCM0_RISR);
    TH2D* hist_gammaPerpCM0_RISR = new TH2D((title+"_gammaPerpCM0_RISR").c_str(), (title+"_gammaPerpCM0_RISR;#gamma_{#perp CM0};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_gammaPerpCM0_RISR);

    TH2D* hist_MSperpCM0_gammaPerpCM0 = new TH2D((title+"_MSperpCM0_gammaPerpCM0").c_str(), (title+"_MSperpCM0_gammaPerpCM0;MS_{#perp CM0};#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 1.);
    hists2.push_back(hist_MSperpCM0_gammaPerpCM0);
    TH2D* hist_MSperpCM0_MQperpCM0 = new TH2D((title+"_MSperpCM0_MQperpCM0").c_str(), (title+"_MSperpCM0_MQperpCM0;MS_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 300.);
    hists2.push_back(hist_MSperpCM0_MQperpCM0);
    TH2D* hist_gammaPerpCM0_MQperpCM0 = new TH2D((title+"_gammaPerpCM0_MQperpCM0").c_str(), (title+"_gammaPerpCM0_MQperpCM0;#gamma_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 300.);
    hists2.push_back(hist_gammaPerpCM0_MQperpCM0);
    
    TH2D* hist_MJ_MQperpCM0 = new TH2D((title+"_MJ_MQperpCM0").c_str(), (title+"_MJ_MQperpCM0;MJ;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MJ_MQperpCM0);
    TH2D* hist_MLa_MQperpCM0 = new TH2D((title+"_MLa_MQperpCM0").c_str(), (title+"_MLa_MQperpCM0;MLa;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLa_MQperpCM0);
    TH2D* hist_MLb_MQperpCM0 = new TH2D((title+"_MLb_MQperpCM0").c_str(), (title+"_MLb_MQperpCM0;MLb;MQ_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLb_MQperpCM0);
    TH2D* hist_MJ_MSperpCM0 = new TH2D((title+"_MJ_MSperpCM0").c_str(), (title+"_MJ_MSperpCM0;MJ;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MJ_MSperpCM0);
    TH2D* hist_MLa_MSperpCM0 = new TH2D((title+"_MLa_MSperpCM0").c_str(), (title+"_MLa_MSperpCM0;MLa;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MLa_MSperpCM0);
    TH2D* hist_MLb_MSperpCM0 = new TH2D((title+"_MLb_MSperpCM0").c_str(), (title+"_MLb_MSperpCM0;MLb;MS_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MLb_MSperpCM0);
    TH2D* hist_MJ_gammaPerpCM0 = new TH2D((title+"_MJ_gammaPerpCM0").c_str(), (title+"_MJ_gammaPerpCM0;MJ;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MJ_gammaPerpCM0);
    TH2D* hist_MLa_gammaPerpCM0 = new TH2D((title+"_MLa_gammaPerpCM0").c_str(), (title+"_MLa_gammaPerpCM0;MLa;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MLa_gammaPerpCM0);
    TH2D* hist_MLb_gammaPerpCM0 = new TH2D((title+"_MLb_gammaPerpCM0").c_str(), (title+"_MLb_gammaPerpCM0;MLb;#gamma_{#perp CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MLb_gammaPerpCM0);

    TH2D* hist_MSCM0_RISR = new TH2D((title+"_MSCM0_RISR").c_str(), (title+"_MSCM0_RISR;MS_{CM0};RISR").c_str(), g_NX/2., 0., 500., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_MSCM0_RISR);
    TH2D* hist_MQCM0_RISR = new TH2D((title+"_MQCM0_RISR").c_str(), (title+"_MQCM0_RISR;MQ_{CM0};RISR").c_str(), g_NX/2., 0., 300., g_NX/2., 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
    hists2.push_back(hist_MQCM0_RISR);
    TH2D* hist_gammaCM0_RISR = new TH2D((title+"_gammaCM0_RISR").c_str(), (title+"_gammaCM0_RISR;#gamma_{CM0};RISR").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_gammaCM0_RISR);

    TH2D* hist_MSCM0_gammaCM0 = new TH2D((title+"_MSCM0_gammaCM0").c_str(), (title+"_MSCM0_gammaCM0;MS_{CM0};#gamma_{CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 1.);
    hists2.push_back(hist_MSCM0_gammaCM0);
    TH2D* hist_MSCM0_MQCM0 = new TH2D((title+"_MSCM0_MQCM0").c_str(), (title+"_MSCM0_MQCM0;MS_{CM0};MQ_{CM0}").c_str(), g_NX/2., 0., 500., g_NX/2., 0., 300.);
    hists2.push_back(hist_MSCM0_MQCM0);
    TH2D* hist_gammaCM0_MQCM0 = new TH2D((title+"_gammaCM0_MQCM0").c_str(), (title+"_gammaCM0_MQCM0;#gamma_{CM0};MQ_{CM0}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 300.);
    hists2.push_back(hist_gammaCM0_MQCM0);
    
    TH2D* hist_MJ_MQCM0 = new TH2D((title+"_MJ_MQCM0").c_str(), (title+"_MJ_MQCM0;MJ;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MJ_MQCM0);
    TH2D* hist_MLa_MQCM0 = new TH2D((title+"_MLa_MQCM0").c_str(), (title+"_MLa_MQCM0;MLa;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLa_MQCM0);
    TH2D* hist_MLb_MQCM0 = new TH2D((title+"_MLb_MQCM0").c_str(), (title+"_MLb_MQCM0;MLb;MQ_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLb_MQCM0);
    TH2D* hist_MJ_MSCM0 = new TH2D((title+"_MJ_MSCM0").c_str(), (title+"_MJ_MSCM0;MJ;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MJ_MSCM0);
    TH2D* hist_MLa_MSCM0 = new TH2D((title+"_MLa_MSCM0").c_str(), (title+"_MLa_MSCM0;MLa;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MLa_MSCM0);
    TH2D* hist_MLb_MSCM0 = new TH2D((title+"_MLb_MSCM0").c_str(), (title+"_MLb_MSCM0;MLb;MS_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 500.);
    hists2.push_back(hist_MLb_MSCM0);
    TH2D* hist_MJ_gammaCM0 = new TH2D((title+"_MJ_gammaCM0").c_str(), (title+"_MJ_gammaCM0;MJ;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MJ_gammaCM0);
    TH2D* hist_MLa_gammaCM0 = new TH2D((title+"_MLa_gammaCM0").c_str(), (title+"_MLa_gammaCM0;MLa;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MLa_gammaCM0);
    TH2D* hist_MLb_gammaCM0 = new TH2D((title+"_MLb_gammaCM0").c_str(), (title+"_MLb_gammaCM0;MLb;#gamma_{CM0}").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 1.);
    hists2.push_back(hist_MLb_gammaCM0);

    TH2D* hist_MLb_MJ = new TH2D((title+"_MLb_MJ").c_str(), (title+"_MLb_MJ;MLb;MJ").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLb_MJ);
    TH2D* hist_MLa_MJ = new TH2D((title+"_MLa_MJ").c_str(), (title+"_MLa_MJ;MLa;MJ").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLa_MJ);
    TH2D* hist_MLa_MLb = new TH2D((title+"_MLa_MLb").c_str(), (title+"_MLa_MLb;MLa;MLb").c_str(), g_NX/2., 0., 300., g_NX/2., 0., 300.);
    hists2.push_back(hist_MLa_MLb);

    TH2D* hist_RISR_MJ = new TH2D((title+"_RISR_MJ").c_str(), (title+"_RISR_MJ;R_{ISR};MJ").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
    hists2.push_back(hist_RISR_MJ);
    TH2D* hist_RISR_MLa = new TH2D((title+"_RISR_MLa").c_str(), (title+"_RISR_MLa;R_{ISR};MLa").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
    hists2.push_back(hist_RISR_MLa);
    TH2D* hist_RISR_MLb = new TH2D((title+"_RISR_MLb").c_str(), (title+"_RISR_MLb;R_{ISR};MLb").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 300.);
    hists2.push_back(hist_RISR_MLb);

    TH2D* hist_RISR_CandML = new TH2D((title+"_RISR_CandML").c_str(), (title+"_RISR_CandML;R_{ISR};CandML").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 150.);
    hists2.push_back(hist_RISR_CandML);
    TH2D* hist_RISR_CandBeta = new TH2D((title+"_RISR_CandBeta").c_str(), (title+"_RISR_CandBeta;R_{ISR};CandBeta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_CandBeta);
    TH2D* hist_RISR_CandDeltaPhiMET = new TH2D((title+"_RISR_CandDeltaPhiMET").c_str(), (title+"_RISR_CandDeltaPhiMET;R_{ISR};CandDeltaPhiMET").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_CandDeltaPhiMET);
    TH2D* hist_CandML_CandBeta = new TH2D((title+"_CandML_CandBeta").c_str(), (title+"_CandML_CandBeta;CandML;CandBeta").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_CandBeta);
    TH2D* hist_CandML_CandDeltaPhiMET = new TH2D((title+"_CandML_CandDeltaPhiMET").c_str(), (title+"_CandML_CandDeltaPhiMET;CandML;CandDeltaPhiMET").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandML_CandDeltaPhiMET);
    TH2D* hist_CandML_CandCosDecayAngle = new TH2D((title+"_CandML_CandCosDecayAngle").c_str(), (title+"_CandML_CandCosDecayAngle;CandML;CandCosDecayAngle").c_str(), g_NX/2., 0., 150., g_NX/2, -1., 1.);
    hists2.push_back(hist_CandML_CandCosDecayAngle);
    TH2D* hist_Beta_CandCosDecayAngle = new TH2D((title+"_Beta_CandCosDecayAngle").c_str(), (title+"_Beta_CandCosDecayAngle;Beta;CandCosDecayAngle").c_str(), g_NX/2., 0., 1., g_NX/2, -1., 1.);
    hists2.push_back(hist_Beta_CandCosDecayAngle);
    TH2D* hist_Beta_CandDeltaPhiMET = new TH2D((title+"_Beta_CandDeltaPhiMET").c_str(), (title+"_Beta_CandDeltaPhiMET;Beta;CandDeltaPhiMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_Beta_CandDeltaPhiMET);
    TH2D* hist_CandCosDecayAngle_CandDeltaPhiMET = new TH2D((title+"_CandCosDecayAngle_CandDeltaPhiMET").c_str(), (title+"_CandCosDecayAngle_CandDeltaPhiMET;CandCosDecayAngle;CandDeltaPhiMET").c_str(), g_NX/2., -1., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandCosDecayAngle_CandDeltaPhiMET);
    TH2D* hist_CandML_RZPara = new TH2D((title+"_CandML_RZPara").c_str(), (title+"_CandML_RZPara;CandML;RZPara").c_str(), g_NX/2., 0., 150., g_NX/2, -1.5, 1.5);
    hists2.push_back(hist_CandML_RZPara);
    TH2D* hist_CandML_RZPerp = new TH2D((title+"_CandML_RZPerp").c_str(), (title+"_CandML_RZPerp;CandML;RZPerp").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandML_RZPerp);
    TH2D* hist_CandML_RZAng = new TH2D((title+"_CandML_RZAng").c_str(), (title+"_CandML_RZAng;CandML;RZAng").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandML_RZAng);
    TH2D* hist_RZPara_RZPerp = new TH2D((title+"_RZPara_RZPerp").c_str(), (title+"_RZPara_RZPerp;RZPara;RZPerp").c_str(), g_NX/2., -1.5, 1.5, g_NX/2, 0., 2.);
    hists2.push_back(hist_RZPara_RZPerp);
    TH2D* hist_CandML_PZPara = new TH2D((title+"_CandML_PZPara").c_str(), (title+"_CandML_PZPara;CandML;PZPara").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 1000.);
    hists2.push_back(hist_CandML_PZPara);
    TH2D* hist_CandML_PZPerp = new TH2D((title+"_CandML_PZPerp").c_str(), (title+"_CandML_PZPerp;CandML;PZPerp").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 1000.);
    hists2.push_back(hist_CandML_PZPerp);

    TH2D* hist_RISR_MatchedCandML = new TH2D((title+"_RISR_MatchedCandML").c_str(), (title+"_RISR_MatchedCandML;R_{ISR};MatchedCandML").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 150.);
    hists2.push_back(hist_RISR_MatchedCandML);
    TH2D* hist_RISR_MatchedCandBeta = new TH2D((title+"_RISR_MatchedCandBeta").c_str(), (title+"_RISR_MatchedCandBeta;R_{ISR};MatchedCandBeta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_MatchedCandBeta);
    TH2D* hist_RISR_MatchedCandDeltaPhiMET = new TH2D((title+"_RISR_MatchedCandDeltaPhiMET").c_str(), (title+"_RISR_MatchedCandDeltaPhiMET;R_{ISR};MatchedCandDeltaPhiMET").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_MatchedCandDeltaPhiMET);
    TH2D* hist_MatchedCandML_MatchedCandBeta = new TH2D((title+"_MatchedCandML_MatchedCandBeta").c_str(), (title+"_MatchedCandML_MatchedCandBeta;MatchedCandML;MatchedCandBeta").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandBeta);
    TH2D* hist_MatchedCandML_MatchedCandDeltaPhiMET = new TH2D((title+"_MatchedCandML_MatchedCandDeltaPhiMET").c_str(), (title+"_MatchedCandML_MatchedCandDeltaPhiMET;MatchedCandML;MatchedCandDeltaPhiMET").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 3.15);
    hists2.push_back(hist_MatchedCandML_MatchedCandDeltaPhiMET);
    TH2D* hist_MatchedCandML_MatchedCandCosDecayAngle = new TH2D((title+"_MatchedCandML_MatchedCandCosDecayAngle").c_str(), (title+"_MatchedCandML_MatchedCandCosDecayAngle;MatchedCandML;MatchedCandCosDecayAngle").c_str(), g_NX/2., 0., 150., g_NX/2, -1., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandCosDecayAngle);

    TH2D* hist_RISR_UnmatchedCandML = new TH2D((title+"_RISR_UnmatchedCandML").c_str(), (title+"_RISR_UnmatchedCandML;R_{ISR};UnmatchedCandML").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 150.);
    hists2.push_back(hist_RISR_UnmatchedCandML);
    TH2D* hist_RISR_UnmatchedCandBeta = new TH2D((title+"_RISR_UnmatchedCandBeta").c_str(), (title+"_RISR_UnmatchedCandBeta;R_{ISR};UnmatchedCandBeta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_UnmatchedCandBeta);
    TH2D* hist_RISR_UnmatchedCandDeltaPhiMET = new TH2D((title+"_RISR_UnmatchedCandDeltaPhiMET").c_str(), (title+"_RISR_UnmatchedCandDeltaPhiMET;R_{ISR};UnmatchedCandDeltaPhiMET").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_UnmatchedCandDeltaPhiMET);
    TH2D* hist_UnmatchedCandML_UnmatchedCandBeta = new TH2D((title+"_UnmatchedCandML_UnmatchedCandBeta").c_str(), (title+"_UnmatchedCandML_UnmatchedCandBeta;UnmatchedCandML;UnmatchedCandBeta").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandBeta);
    TH2D* hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET = new TH2D((title+"_UnmatchedCandML_UnmatchedCandDeltaPhiMET").c_str(), (title+"_UnmatchedCandML_UnmatchedCandDeltaPhiMET;UnmatchedCandML;UnmatchedCandDeltaPhiMET").c_str(), g_NX/2., 0., 150., g_NX/2, 0., 3.15);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET);
    TH2D* hist_UnmatchedCandML_UnmatchedCandCosDecayAngle = new TH2D((title+"_UnmatchedCandML_UnmatchedCandCosDecayAngle").c_str(), (title+"_UnmatchedCandML_UnmatchedCandCosDecayAngle;UnmatchedCandML;UnmatchedCandCosDecayAngle").c_str(), g_NX/2., 0., 150., g_NX/2, -1., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandCosDecayAngle);

    TEfficiency* eff_METtrig = new TEfficiency((title+"_eff_METtrig").c_str(), "Efficiency of MET trigger;Eff;MET [GeV]", g_NX, 0., 700.);
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
          
          //if(do_FilterDilepton)
          //  if(SF.DileptonEvent(base))
          //    continue;
          
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

          //if(MET < 50.) // ATLAS
          if(MET < 150.) // PreSelection
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "MET > 150");
          
          // apply trigger to data and FullSim events
          if(!base->METORtrigger && !is_FastSim) // PreSelection
          //if(!base->SingleElectrontrigger && !base->SingleMuontrigger && !is_FastSim) // ATLAS
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "MET trigger");

          if(PTISR < 200.) // PreSelection
          //if(PTISR < 300.) // SR
	    continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "PTISR > 200");
            
          if(RISR < 0.5 || RISR > 1.0) // PreSelection
          //if(RISR < 0.4 || RISR > 0.7) // CR
          //if(RISR < 0.7 || RISR > 1.0)
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "R_{ISR} > 0.5");

          // Cleaning cuts...
          double dphiCMI = base->dphiCMI;
          double PTCM = base->PTCM;
          double x = fabs(dphiCMI);
          
          // PreSelection
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

          //if(NbjetISR + NbjetS != 2) continue; // CR
          //if(NbjetISR + NbjetS > 1) continue; // SR & ATLAS
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "B-Veto");

          //if(Nlep != 3) continue;
          //if(NjetS != 0) continue; // SR

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
              flavor = LepFlavor::kElectron;
            else
              flavor = LepFlavor::kMuon;
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

          TVector3 TV3_MET;
          TV3_MET.SetPtEtaPhi(MET, 0., base->MET_phi);

          // Leptonic Candidates
          std::vector<L_Cand> V_lep_cands;
          for(int i = 0; i < Nlep-1; i++){
            for(int j = i+1; j < Nlep; j++){
              if(list_leps[i].Flavor() == list_leps[j].Flavor() && list_leps[i].Charge() != list_leps[j].Charge()){ // OSSF
                if(!CandSameHemi && ((NjetS == 0 && Nlep == 2) || (inVec((*base->index_lep_a),i) == inVec((*base->index_lep_a),j)))){ // same hemisphere
                  Particle lep_1;
                  lep_1.SetPtEtaPhiM( base->PT_lep->at(i),
                                      base->Eta_lep->at(i),
                                      base->Phi_lep->at(i),
                                      std::max(0.,base->M_lep->at(i)) );
                  Particle lep_2;
                  lep_2.SetPtEtaPhiM( base->PT_lep->at(j),
                    	              base->Eta_lep->at(j),
                    	              base->Phi_lep->at(j),
                    	              std::max(0.,base->M_lep->at(j)) );
                  if((*base->Index_lep)[i] >= 0)
                    lep_1.SetMomPDGID((*base->genMomPDGID_lep)[(*base->Index_lep)[i]]);
                  else
                    lep_1.SetMomPDGID(0);
                  if((*base->Index_lep)[j] >= 0)
                    lep_2.SetMomPDGID((*base->genMomPDGID_lep)[(*base->Index_lep)[j]]);
                  else
                    lep_2.SetMomPDGID(0);
                  ParticleList cand_list_parts;
                  cand_list_parts.push_back(lep_1);
                  cand_list_parts.push_back(lep_2);
                  ConstRestFrameList cand_list_frames;
                  cand_list_frames.Add(La);
                  cand_list_frames.Add(Lb);
                  L_Cand cand(cand_list_parts);
                  //L_Cand cand(cand_list_parts, cand_list_frames);
                  cand.SetFlavor(list_leps[i].Flavor());
                  V_lep_cands.push_back(cand);
                }
              }
            }
          }

          cand_matching(V_lep_cands);
          int N_V_lep_cands = V_lep_cands.size();
          //if(N_V_lep_cands < 1) continue;

          TLorentzVector TLV_LAB, TLV_CM, TLV_S;
          TLV_LAB.SetPtEtaPhiM(base->LAB_Pt, base->LAB_Eta, base->LAB_Phi, base->LAB_M);
          TLV_CM.SetPtEtaPhiM(base->PTCM, base->EtaCM, base->PhiCM, base->MCM);
          TLV_S.SetPtEtaPhiM(base->PTS, base->EtaS, base->PhiS, base->MS);
          TVector3 beta_CM = TLV_CM.BoostVector();
          TLorentzVector S_CM = TLV_S;
          S_CM.Boost(-beta_CM); // S in CM frame

          for(int i = 0; i < N_V_lep_cands; i++){
            TLorentzVector TLV_Cand = V_lep_cands[i].TLV();
            TLorentzVector TLV_Cand_CM = TLV_Cand;
            TLV_Cand_CM.Boost(-beta_CM);
            TVector3 beta_Cand_CM = TLV_Cand_CM.BoostVector();
            Particle part_Cand_Child = V_lep_cands[i].Cand_PartPlus();
            TLorentzVector TLV_Cand_Child;
            TLV_Cand_Child.SetPtEtaPhiM(part_Cand_Child.Pt(), part_Cand_Child.Eta(), part_Cand_Child.Phi(), part_Cand_Child.M());
            TLV_Cand_Child.Boost(-beta_CM);
            TLV_Cand_Child.Boost(-beta_Cand_CM);
            double CosDecayAngle = TLV_Cand_Child.Vect().Unit().Dot(beta_Cand_CM.Unit());
            double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
            double PZPerp = TLV_Cand_CM.Vect().Cross(S_CM.Vect().Unit()).Mag();
            double RZPara = PZPara/S_CM.Vect().Mag();
            double RZPerp = PZPerp/S_CM.Vect().Mag();
            double RZAng = atan(RZPara/RZPerp);
            hist_CandML->Fill(V_lep_cands[i].Mass(), weight);
            hist_CandBeta->Fill(V_lep_cands[i].Beta(), weight);
            hist_CandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandCosDecayAngle->Fill(CosDecayAngle, weight);
            hist_RZPara->Fill(RZPara, weight);
            hist_RZPerp->Fill(RZPerp, weight);
            hist_RISR_CandML->Fill(RISR, V_lep_cands[i].Mass(), weight);
            hist_RISR_CandBeta->Fill(RISR, V_lep_cands[i].Beta(), weight);
            hist_RISR_CandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_CandBeta->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].Beta(), weight);
            hist_CandML_CandDeltaPhiMET->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_CandCosDecayAngle->Fill(V_lep_cands[i].Mass(), CosDecayAngle, weight);
            hist_CandML_RZPara->Fill(V_lep_cands[i].Mass(), RZPara, weight);
            hist_CandML_RZPerp->Fill(V_lep_cands[i].Mass(), RZPerp, weight);
            hist_RZAng->Fill(RZAng, weight);
            hist_PZPara->Fill(RZPara * S_CM.Vect().Mag(), weight);
            hist_PZPerp->Fill(RZPerp * S_CM.Vect().Mag(), weight);
            hist_Beta_CandCosDecayAngle->Fill(V_lep_cands[i].Beta(), CosDecayAngle, weight);
            hist_Beta_CandDeltaPhiMET->Fill(V_lep_cands[i].Beta(), V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandCosDecayAngle_CandDeltaPhiMET->Fill(CosDecayAngle, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_RZAng->Fill(V_lep_cands[i].Mass(), RZAng, weight);
            hist_RZPara_RZPerp->Fill(RZPara, RZPerp, weight);
            hist_CandML_PZPara->Fill(V_lep_cands[i].Mass(), PZPara, weight);
            hist_CandML_PZPerp->Fill(V_lep_cands[i].Mass(), PZPerp, weight);
            if(V_lep_cands[i].Match() == kMatched){
              hist_MatchedCandML->Fill(V_lep_cands[i].Mass(), weight);
              hist_MatchedCandBeta->Fill(V_lep_cands[i].Beta(), weight);
              hist_MatchedCandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandCosDecayAngle->Fill(CosDecayAngle, weight);
              hist_RISR_MatchedCandML->Fill(RISR, V_lep_cands[i].Mass(), weight);
              hist_RISR_MatchedCandBeta->Fill(RISR, V_lep_cands[i].Beta(), weight);
              hist_RISR_MatchedCandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandML_MatchedCandBeta->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].Beta(), weight);
              hist_MatchedCandML_MatchedCandDeltaPhiMET->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandML_MatchedCandCosDecayAngle->Fill(V_lep_cands[i].Mass(), CosDecayAngle, weight);
            }
            else{
              hist_UnmatchedCandML->Fill(V_lep_cands[i].Mass(), weight);
              hist_UnmatchedCandBeta->Fill(V_lep_cands[i].Beta(), weight);
              hist_UnmatchedCandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandCosDecayAngle->Fill(CosDecayAngle, weight);
              hist_RISR_UnmatchedCandML->Fill(RISR, V_lep_cands[i].Mass(), weight);
              hist_RISR_UnmatchedCandBeta->Fill(RISR, V_lep_cands[i].Beta(), weight);
              hist_RISR_UnmatchedCandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandML_UnmatchedCandBeta->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].Beta(), weight);
              hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET->Fill(V_lep_cands[i].Mass(), V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandML_UnmatchedCandCosDecayAngle->Fill(V_lep_cands[i].Mass(), CosDecayAngle, weight);
            }
          }

          // Fill hists, effs, etc.
          hist_MET->Fill(MET, weight);
          hist_RISR_PTISR->Fill(RISR, PTISR, weight);
          hist_RISR_Mperp->Fill(RISR, Mperp, weight);
          hist_gammaPerp_RISR->Fill(base->gammaT, RISR, weight);
          hist_RISR_mL->Fill(RISR, ML, weight);
          hist_dphiCMI_PTCM->Fill(dphiCMI, PTCM, weight);
          hist_dphiMETV_PTISR->Fill(dphiMET_V, PTISR, weight);

          hist_MSperpCM0_RISR->Fill(MSperpCM0, RISR, weight);
          hist_MQperpCM0_RISR->Fill(MQperpCM0, RISR, weight);
          hist_gammaPerpCM0_RISR->Fill(gammaPerpCM0, RISR, weight);

          hist_MSperpCM0_gammaPerpCM0->Fill(MSperpCM0, gammaPerpCM0, weight);
          hist_MSperpCM0_MQperpCM0->Fill(MSperpCM0, MQperpCM0, weight);
          hist_gammaPerpCM0_MQperpCM0->Fill(gammaPerpCM0, MQperpCM0, weight);

          hist_MJ_MQperpCM0->Fill(MJ, MQperpCM0, weight);
          hist_MLa_MQperpCM0->Fill(MLa, MQperpCM0, weight);
          hist_MLb_MQperpCM0->Fill(MLb, MQperpCM0, weight);
          hist_MJ_MSperpCM0->Fill(MJ, MSperpCM0, weight);
          hist_MLa_MSperpCM0->Fill(MLa, MSperpCM0, weight);
          hist_MLb_MSperpCM0->Fill(MLb, MSperpCM0, weight);
          hist_MJ_gammaPerpCM0->Fill(MJ, gammaPerpCM0, weight);
          hist_MLa_gammaPerpCM0->Fill(MLa, gammaPerpCM0, weight);
          hist_MLb_gammaPerpCM0->Fill(MLb, gammaPerpCM0, weight);

          hist_MSCM0_RISR->Fill(MSCM0, RISR, weight);
          hist_MQCM0_RISR->Fill(MQCM0, RISR, weight);
          hist_gammaCM0_RISR->Fill(gammaCM0, RISR, weight);

          hist_MSCM0_gammaCM0->Fill(MSCM0, gammaCM0, weight);
          hist_MSCM0_MQCM0->Fill(MSCM0, MQCM0, weight);
          hist_gammaCM0_MQCM0->Fill(gammaCM0, MQCM0, weight);

          hist_MJ_MQCM0->Fill(MJ, MQCM0, weight);
          hist_MLa_MQCM0->Fill(MLa, MQCM0, weight);
          hist_MLb_MQCM0->Fill(MLb, MQCM0, weight);
          hist_MJ_MSCM0->Fill(MJ, MSCM0, weight);
          hist_MLa_MSCM0->Fill(MLa, MSCM0, weight);
          hist_MLb_MSCM0->Fill(MLb, MSCM0, weight);
          hist_MJ_gammaCM0->Fill(MJ, gammaCM0, weight);
          hist_MLa_gammaCM0->Fill(MLa, gammaCM0, weight);
          hist_MLb_gammaCM0->Fill(MLb, gammaCM0, weight);

          hist_MLa_MJ->Fill(MLa, MJ, weight);
          hist_MLb_MJ->Fill(MLb, MJ, weight);
          hist_MLa_MLb->Fill(MLa, MLb, weight);
          hist_RISR_MJ->Fill(RISR, MJ, weight);
          hist_RISR_MLa->Fill(RISR, MLa, weight);
          hist_RISR_MLb->Fill(RISR, MLb, weight);

          hist_RISR->Fill(RISR, weight);    
          hist_PTISR->Fill(PTISR, weight);
          hist_gammaPerp->Fill(base->gammaT, weight);
          hist_MQperp->Fill(Mperp, weight);
          hist_MSperpCM0->Fill(MSperpCM0, weight);
          hist_MQperpCM0->Fill(MQperpCM0, weight);
          hist_gammaPerpCM0->Fill(gammaPerpCM0, weight);
          hist_MSCM0->Fill(MSCM0, weight);
          hist_MQCM0->Fill(MQCM0, weight);
          hist_gammaCM0->Fill(gammaCM0, weight);
          hist_ML->Fill(ML, weight);
          hist_MJ->Fill(MJ, weight);

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
          int num_e = std::count(flavors.begin(), flavors.end(), LepFlavor::kElectron);
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
            if(EC_X < 4){
              hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
              if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
              if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
              if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            }
            if (EC_X == 4){ // cand counting
              EC_X += N_V_lep_cands;
              hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
              if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
              if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
              if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            }
          }
          else if(cat_Nleps == 3) {
            EC_X = 6;
            EC_X += num_e;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            EC_X = 10;
            EC_X += N_V_lep_cands;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
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
            EC_X = 13;
            EC_X += num_e;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            EC_X = 19;
            EC_X += N_V_lep_cands;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
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
          hist_CandsEventCount->SetBinContent(0,EC_Y,hist_CandsEventCount->GetBinContent(0,EC_Y)+weight); // normalized to selection
          if(is_bkg){ // total SM bkg
            hist_EventCount->SetBinContent(0,EC_bins,hist_EventCount->GetBinContent(0,EC_bins)+weight); // normalized to selection
            hist_CandsEventCount->SetBinContent(0,EC_bins,hist_CandsEventCount->GetBinContent(0,EC_bins)+weight); // normalized to selection
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

  Plot_EventCount((TH2*)hist_CandsEventCount->Clone("CandsEventCount_Scaled"), true, lumi, false, 0., false, false);
  Plot_EventCount((TH2*)hist_CandsEventCount->Clone("CandsEventCount_SoBAllBKG"), true, lumi, false, 0., true, true);
  Plot_EventCount(hist_CandsEventCount, false, lumi, false, 0., false, false);
  Plot_EventCount((TH2*)hist_CandsZbi->Clone("CandsEventCount_SoB"), true, lumi, false, 0., true, false);
  Plot_EventCount(hist_CandsZbi, true, lumi, true, 0.2, false, false);

  Plot_CutFlow(vect_hist_cutflow, true, lumi, signal_boost, vec_samples);

  Long64_t end = gSystem->Now();
  std::cout << "Time to process " << (end-start)/1000.0 << " seconds" << std::endl;
  gApplication->Terminate(0);
} // End of macro
