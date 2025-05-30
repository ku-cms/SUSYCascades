#include "PlottingTools.h"

void Plot_Advanced(){

  Long64_t start = gSystem->Now();
  RestFrames::SetStyle();
  //InitRJRtree();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/"; // v1
  //string NtuplePath = "/local-scratch/zflowers/NTUPLES/old_HADD_Ele8GeV/"; // v0
  //string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Advanced_NTUPLES_v1_Cands_";

  const double CandMinMass = 2.;
  const double CandMaxMass = 30.;
  output_root_file += "CandMinMass"+std::to_string(int(CandMinMass))+"GeV_";
  output_root_file += "CandMaxMass"+std::to_string(int(CandMaxMass))+"GeV_";

  bool CandSameHemi = false;
  if(CandSameHemi) output_root_file += "SameHemi_";

  //g_Label = "TESTING";
  //g_Label = "PreSelection";
  g_Label = "PreSelection CandCosCM<0.8 LowMass";
  //g_Label = "PreSelection & B-Veto & 0S J";
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
  vector<TH1*> hist_stack_CandBetaCM;
  hist_stacks.push_back(&hist_stack_CandBetaCM);
  vector<TH1*> hist_stack_CandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_CandDeltaPhiMET);
  vector<TH1*> hist_stack_CandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_CandCosDecayAngle);
  vector<TH1*> hist_stack_CandCosDecayAngleCM;
  hist_stacks.push_back(&hist_stack_CandCosDecayAngleCM);
  vector<TH1*> hist_stack_RZPara;
  hist_stacks.push_back(&hist_stack_RZPara);
  vector<TH1*> hist_stack_RZPerp;
  hist_stacks.push_back(&hist_stack_RZPerp);
  vector<TH1*> hist_stack_PZAng;
  hist_stacks.push_back(&hist_stack_PZAng);
  vector<TH1*> hist_stack_PZPara;
  hist_stacks.push_back(&hist_stack_PZPara);
  vector<TH1*> hist_stack_PZPerp;
  hist_stacks.push_back(&hist_stack_PZPerp);
  vector<TH1*> hist_stack_RZParaLABMET;
  hist_stacks.push_back(&hist_stack_RZParaLABMET);
  vector<TH1*> hist_stack_RZPerpLABMET;
  hist_stacks.push_back(&hist_stack_RZPerpLABMET);
  vector<TH1*> hist_stack_PZAngLABMET;
  hist_stacks.push_back(&hist_stack_PZAngLABMET);
  vector<TH1*> hist_stack_MRZPara;
  hist_stacks.push_back(&hist_stack_MRZPara);
  vector<TH1*> hist_stack_MRZPerp;
  hist_stacks.push_back(&hist_stack_MRZPerp);
  vector<TH1*> hist_stack_MRZParaLABMET;
  hist_stacks.push_back(&hist_stack_MRZParaLABMET);
  vector<TH1*> hist_stack_MRZPerpLABMET;
  hist_stacks.push_back(&hist_stack_MRZPerpLABMET);

  vector<TH1*> hist_stack_MatchedCandML;
  hist_stacks.push_back(&hist_stack_MatchedCandML);
  vector<TH1*> hist_stack_MatchedCandBeta;
  hist_stacks.push_back(&hist_stack_MatchedCandBeta);
  vector<TH1*> hist_stack_MatchedCandBetaCM;
  hist_stacks.push_back(&hist_stack_MatchedCandBetaCM);
  vector<TH1*> hist_stack_MatchedCandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_MatchedCandDeltaPhiMET);
  vector<TH1*> hist_stack_MatchedCandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_MatchedCandCosDecayAngle);
  vector<TH1*> hist_stack_MatchedCandCosDecayAngleCM;
  hist_stacks.push_back(&hist_stack_MatchedCandCosDecayAngleCM);
  vector<TH1*> hist_stack_UnmatchedCandML;
  hist_stacks.push_back(&hist_stack_UnmatchedCandML);
  vector<TH1*> hist_stack_UnmatchedCandBeta;
  hist_stacks.push_back(&hist_stack_UnmatchedCandBeta);
  vector<TH1*> hist_stack_UnmatchedCandBetaCM;
  hist_stacks.push_back(&hist_stack_UnmatchedCandBetaCM);
  vector<TH1*> hist_stack_UnmatchedCandDeltaPhiMET;
  hist_stacks.push_back(&hist_stack_UnmatchedCandDeltaPhiMET);
  vector<TH1*> hist_stack_UnmatchedCandCosDecayAngle;
  hist_stacks.push_back(&hist_stack_UnmatchedCandCosDecayAngle);
  vector<TH1*> hist_stack_UnmatchedCandCosDecayAngleCM;
  hist_stacks.push_back(&hist_stack_UnmatchedCandCosDecayAngleCM);

  // hists for holding number of events
  const int EC_bins = vec_samples.size() + 1;
  const int Zbi_bins = map_vsignals.size();
  int vec_samples_index = 0;
  int Zbi_samples_index = 0;
  TH2D* hist_EventCount = new TH2D("EventCount", "EventCount", 22, 0, 22, EC_bins, 0, EC_bins);
  hist_EventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_Zbi = new TH2D("Zbi", "Zbi", 22, 0, 22, Zbi_bins, 0, Zbi_bins);

  // Event Counting by cands
  TH2D* hist_CandsEventCount = new TH2D("CandsEventCount", "CandsEventCount", 14, 0, 14, EC_bins, 0, EC_bins);
  hist_CandsEventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_CandsZbi = new TH2D("CandsZbi", "CandsZbi", 14, 0, 14, Zbi_bins, 0, Zbi_bins);

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
    TH1D* hist_MQperp = new TH1D((title+"_MQperp").c_str(), (title+"_MQperp;M_{#perp}").c_str(), g_NX/2, 0., 175.);
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
    TH1D* hist_ML = new TH1D((title+"_ML").c_str(), (title+"_ML;ML").c_str(), g_NX/2, 0., 175.);
    hists1.push_back(hist_ML);
    hist_stack_ML.push_back(hist_ML);
    TH1D* hist_MJ = new TH1D((title+"_MJ").c_str(), (title+"_MJ;MJ").c_str(), g_NX/2, 0., 175.);
    hists1.push_back(hist_MJ);
    hist_stack_MJ.push_back(hist_MJ);

    TH1D* hist_CandML = new TH1D((title+"_CandML").c_str(), (title+"_CandML;CandML").c_str(), g_NX/2, CandMinMass, CandMaxMass);
    hists1.push_back(hist_CandML);
    hist_stack_CandML.push_back(hist_CandML);
    TH1D* hist_CandBeta = new TH1D((title+"_CandBeta").c_str(), (title+"_CandBeta;#beta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_CandBeta);
    hist_stack_CandBeta.push_back(hist_CandBeta);
    TH1D* hist_CandBetaCM = new TH1D((title+"_CandBetaCM").c_str(), (title+"_CandBetaCM;#beta^{CM}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_CandBetaCM);
    hist_stack_CandBetaCM.push_back(hist_CandBetaCM);
    TH1D* hist_CandDeltaPhiMET = new TH1D((title+"_CandDeltaPhiMET").c_str(), (title+"_CandDeltaPhiMET;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_CandDeltaPhiMET);
    hist_stack_CandDeltaPhiMET.push_back(hist_CandDeltaPhiMET);
    TH1D* hist_CandCosDecayAngle = new TH1D((title+"_CandCosDecayAngle").c_str(), (title+"_CandCosDecayAngle;Z* candidate cos#theta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_CandCosDecayAngle);
    hist_stack_CandCosDecayAngle.push_back(hist_CandCosDecayAngle);
    TH1D* hist_CandCosDecayAngleCM = new TH1D((title+"_CandCosDecayAngleCM").c_str(), (title+"_CandCosDecayAngleCM;Z* candidate cos#theta CM").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_CandCosDecayAngleCM);
    hist_stack_CandCosDecayAngleCM.push_back(hist_CandCosDecayAngleCM);
    TH1D* hist_RZPara = new TH1D((title+"_RZPara").c_str(), (title+"_RZPara;RZPara").c_str(), g_NX/2, -1., 1.);
    hists1.push_back(hist_RZPara);
    hist_stack_RZPara.push_back(hist_RZPara);
    TH1D* hist_RZPerp = new TH1D((title+"_RZPerp").c_str(), (title+"_RZPerp;RZPerp").c_str(), g_NX/2, 0., 1.5);
    hists1.push_back(hist_RZPerp);
    hist_stack_RZPerp.push_back(hist_RZPerp);
    TH1D* hist_PZAng = new TH1D((title+"_PZAng").c_str(), (title+"_PZAng;PZAng").c_str(), g_NX/2, -1.6, 1.6);
    hists1.push_back(hist_PZAng);
    hist_stack_PZAng.push_back(hist_PZAng);
    TH1D* hist_PZPara = new TH1D((title+"_PZPara").c_str(), (title+"_PZPara;PZPara").c_str(), g_NX/2, 0., 600.);
    hists1.push_back(hist_PZPara);
    hist_stack_PZPara.push_back(hist_PZPara);
    TH1D* hist_PZPerp = new TH1D((title+"_PZPerp").c_str(), (title+"_PZPerp;PZPerp").c_str(), g_NX/2, 0., 600.);
    hists1.push_back(hist_PZPerp);
    hist_stack_PZPerp.push_back(hist_PZPerp);
    TH1D* hist_RZParaLABMET = new TH1D((title+"_RZParaLABMET").c_str(), (title+"_RZParaLABMET;RZParaLABMET").c_str(), g_NX/2, -1.5, 1.5);
    hists1.push_back(hist_RZParaLABMET);
    hist_stack_RZParaLABMET.push_back(hist_RZParaLABMET);
    TH1D* hist_RZPerpLABMET = new TH1D((title+"_RZPerpLABMET").c_str(), (title+"_RZPerpLABMET;RZPerpLABMET").c_str(), g_NX/2, 0., 1.5);
    hists1.push_back(hist_RZPerpLABMET);
    hist_stack_RZPerpLABMET.push_back(hist_RZPerpLABMET);
    TH1D* hist_PZAngLABMET = new TH1D((title+"_PZAngLABMET").c_str(), (title+"_PZAngLABMET;PZAngLABMET").c_str(), g_NX/2, -1.6, 1.6);
    hists1.push_back(hist_PZAngLABMET);
    hist_stack_PZAngLABMET.push_back(hist_PZAngLABMET);

    TH1D* hist_MRZPara = new TH1D((title+"_MRZPara").c_str(), (title+"_MRZPara;M/RZPara").c_str(), g_NX/2, 0., 500.);
    hists1.push_back(hist_MRZPara);
    hist_stack_MRZPara.push_back(hist_MRZPara);
    TH1D* hist_MRZPerp = new TH1D((title+"_MRZPerp").c_str(), (title+"_MRZPerp;M/RZPerp").c_str(), g_NX/2, 0., 500.);
    hists1.push_back(hist_MRZPerp);
    hist_stack_MRZPerp.push_back(hist_MRZPerp);
    TH1D* hist_MRZParaLABMET = new TH1D((title+"_MRZParaLABMET").c_str(), (title+"_MRZParaLABMET;M/RZParaLABMET").c_str(), g_NX/2, 0., 500.);
    hists1.push_back(hist_MRZParaLABMET);
    hist_stack_MRZParaLABMET.push_back(hist_MRZParaLABMET);
    TH1D* hist_MRZPerpLABMET = new TH1D((title+"_MRZPerpLABMET").c_str(), (title+"_MRZPerpLABMET;M/RZPerpLABMET").c_str(), g_NX/2, 0., 500.);
    hists1.push_back(hist_MRZPerpLABMET);
    hist_stack_MRZPerpLABMET.push_back(hist_MRZPerpLABMET);

    TH1D* hist_MatchedCandML = new TH1D((title+"_MatchedCandML").c_str(), (title+"_MatchedCandML;Matched Z* candidate M [GeV]").c_str(), g_NX/2, CandMinMass, CandMaxMass);
    hists1.push_back(hist_MatchedCandML);
    hist_stack_MatchedCandML.push_back(hist_MatchedCandML);
    TH1D* hist_MatchedCandBeta = new TH1D((title+"_MatchedCandBeta").c_str(), (title+"_MatchedCandBeta;Matched Z* candidate #beta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_MatchedCandBeta);
    hist_stack_MatchedCandBeta.push_back(hist_MatchedCandBeta);
    TH1D* hist_MatchedCandBetaCM = new TH1D((title+"_MatchedCandBetaCM").c_str(), (title+"_MatchedCandBetaCM;Matched Z* candidate #beta^{CM}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_MatchedCandBetaCM);
    hist_stack_MatchedCandBetaCM.push_back(hist_MatchedCandBetaCM);
    TH1D* hist_MatchedCandDeltaPhiMET = new TH1D((title+"_MatchedCandDeltaPhiMET").c_str(), (title+"_MatchedCandDeltaPhiMET;Matched #Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_MatchedCandDeltaPhiMET);
    hist_stack_MatchedCandDeltaPhiMET.push_back(hist_MatchedCandDeltaPhiMET);
    TH1D* hist_MatchedCandCosDecayAngle = new TH1D((title+"_MatchedCandCosDecayAngle").c_str(), (title+"_MatchedCandCosDecayAngle;Matched cos#theta_{Z}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_MatchedCandCosDecayAngle);
    hist_stack_MatchedCandCosDecayAngle.push_back(hist_MatchedCandCosDecayAngle);
    TH1D* hist_MatchedCandCosDecayAngleCM = new TH1D((title+"_MatchedCandCosDecayAngleCM").c_str(), (title+"_MatchedCandCosDecayAngleCM;Matched cos#theta_{Z} CM").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_MatchedCandCosDecayAngleCM);
    hist_stack_MatchedCandCosDecayAngleCM.push_back(hist_MatchedCandCosDecayAngleCM);

    TH1D* hist_UnmatchedCandML = new TH1D((title+"_UnmatchedCandML").c_str(), (title+"_UnmatchedCandML;Unmatched Z* candidate M [GeV]").c_str(), g_NX/2, CandMinMass, CandMaxMass);
    hists1.push_back(hist_UnmatchedCandML);
    hist_stack_UnmatchedCandML.push_back(hist_UnmatchedCandML);
    TH1D* hist_UnmatchedCandBeta = new TH1D((title+"_UnmatchedCandBeta").c_str(), (title+"_UnmatchedCandBeta;Unmatched Z* candidate #beta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_UnmatchedCandBeta);
    hist_stack_UnmatchedCandBeta.push_back(hist_UnmatchedCandBeta);
    TH1D* hist_UnmatchedCandBetaCM = new TH1D((title+"_UnmatchedCandBetaCM").c_str(), (title+"_UnmatchedCandBetaCM;Unmatched Z* candidate #beta^{CM}").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_UnmatchedCandBetaCM);
    hist_stack_UnmatchedCandBetaCM.push_back(hist_UnmatchedCandBetaCM);
    TH1D* hist_UnmatchedCandDeltaPhiMET = new TH1D((title+"_UnmatchedCandDeltaPhiMET").c_str(), (title+"_UnmatchedCandDeltaPhiMET;Unmatched #Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2, 0., 3.15);
    hists1.push_back(hist_UnmatchedCandDeltaPhiMET);
    hist_stack_UnmatchedCandDeltaPhiMET.push_back(hist_UnmatchedCandDeltaPhiMET);
    TH1D* hist_UnmatchedCandCosDecayAngle = new TH1D((title+"_UnmatchedCandCosDecayAngle").c_str(), (title+"_UnmatchedCandCosDecayAngle;Unmatched Z* candidate cos#theta").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_UnmatchedCandCosDecayAngle);
    hist_stack_UnmatchedCandCosDecayAngle.push_back(hist_UnmatchedCandCosDecayAngle);
    TH1D* hist_UnmatchedCandCosDecayAngleCM = new TH1D((title+"_UnmatchedCandCosDecayAngleCM").c_str(), (title+"_UnmatchedCandCosDecayAngleCM;Unmatched Z* candidate cos#theta CM").c_str(), g_NX/2, 0., 1.);
    hists1.push_back(hist_UnmatchedCandCosDecayAngleCM);
    hist_stack_UnmatchedCandCosDecayAngleCM.push_back(hist_UnmatchedCandCosDecayAngleCM);

    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 1000.);
    hists2.push_back(hist_RISR_PTISR);
    TH2D* hist_RISR_Mperp = new TH2D((title+"_RISR_Mperp").c_str(), (title+"_RISR_Mperp;R_{ISR};M_{#perp} [GeV]").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 100.);
    hists2.push_back(hist_RISR_Mperp);
    TH2D* hist_dphiCMI_PTCM = new TH2D((title+"_dphiCMI_PTCM").c_str(), (title+"_dphiCMI_PTCM;#Delta #phi_{(CM,I)};p_{T}^{CM}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 0., 500.);
    hists2.push_back(hist_dphiCMI_PTCM);
    TH2D* hist_dphiMETV_PTISR = new TH2D((title+"_dphiMETV_PTISR").c_str(), (title+"_dphiMETV_PTISR;#Delta #phi_{(I,V)};p_{T}^{ISR}").c_str(), g_NX/2., 0., 3.15, g_NX/2., 200., 800.);
    hists2.push_back(hist_dphiMETV_PTISR);
    TH2D* hist_gammaPerp_RISR = new TH2D((title+"_gammaPerp_RISR").c_str(), (title+"_gammaPerp_RISR;#gamma_{#perp};R_{ISR}").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_gammaPerp_RISR);
    TH2D* hist_RISR_mllLEAD = new TH2D((title+"_RISR_mllLEAD").c_str(), (title+"_RISR_mllLEAD;R_{ISR};m_{ll}LEAD").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 175.);
    hists2.push_back(hist_RISR_mllLEAD);
    TH2D* hist_RISR_mL = new TH2D((title+"_RISR_mL").c_str(), (title+"_RISR_mL;R_{ISR};mL").c_str(), g_NX/2., 0.5, 1., g_NX/2., 0., 175.); // mass of leptonic system
    hists2.push_back(hist_RISR_mL);

    TH2D* hist_MSperpCM0_RISR = new TH2D((title+"_MSperpCM0_RISR").c_str(), (title+"_MSperpCM0_RISR;MS_{#perp CM0};R_{ISR}").c_str(), g_NX/2., 0., 500., g_NX/2., 0.5, 1.);
    hists2.push_back(hist_MSperpCM0_RISR);
    TH2D* hist_MQperpCM0_RISR = new TH2D((title+"_MQperpCM0_RISR").c_str(), (title+"_MQperpCM0_RISR;MQ_{#perp CM0};R_{ISR}").c_str(), g_NX/2., 0., 300., g_NX/2., 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
    hists2.push_back(hist_MQperpCM0_RISR);
    TH2D* hist_gammaPerpCM0_RISR = new TH2D((title+"_gammaPerpCM0_RISR").c_str(), (title+"_gammaPerpCM0_RISR;#gamma_{#perp CM0};R_{ISR}").c_str(), g_NX/2., 0., 1., g_NX/2., 0.5, 1.);
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

    TH2D* hist_RISR_CandML = new TH2D((title+"_RISR_CandML").c_str(), (title+"_RISR_CandML;R_{ISR};Z* candidate M [GeV]").c_str(), g_NX/2., 0.5, 1., g_NX/2, CandMinMass, CandMaxMass);
    hists2.push_back(hist_RISR_CandML);
    TH2D* hist_RISR_CandBeta = new TH2D((title+"_RISR_CandBeta").c_str(), (title+"_RISR_CandBeta;R_{ISR};Z* candidate #beta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_CandBeta);
    TH2D* hist_RISR_CandBetaCM = new TH2D((title+"_RISR_CandBetaCM").c_str(), (title+"_RISR_CandBetaCM;R_{ISR};Z* candidate #beta^{CM}").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_CandBetaCM);
    TH2D* hist_RISR_CandDeltaPhiMET = new TH2D((title+"_RISR_CandDeltaPhiMET").c_str(), (title+"_RISR_CandDeltaPhiMET;R_{ISR};#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_CandDeltaPhiMET);
    TH2D* hist_CandML_CandBeta = new TH2D((title+"_CandML_CandBeta").c_str(), (title+"_CandML_CandBeta;Z* candidate M [GeV];Z* candidate #beta").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_CandBeta);
    TH2D* hist_CandML_CandBetaCM = new TH2D((title+"_CandML_CandBetaCM").c_str(), (title+"_CandML_CandBetaCM;Z* candidate M [GeV];Z* candidate #beta^{CM}").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_CandBetaCM);
    TH2D* hist_CandML_CandDeltaPhiMET = new TH2D((title+"_CandML_CandDeltaPhiMET").c_str(), (title+"_CandML_CandDeltaPhiMET;Z* candidate M [GeV];#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandML_CandDeltaPhiMET);
    TH2D* hist_CandML_CandCosDecayAngle = new TH2D((title+"_CandML_CandCosDecayAngle").c_str(), (title+"_CandML_CandCosDecayAngle;Z* candidate M [GeV];Z* candidate cos#theta").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_CandCosDecayAngle);
    TH2D* hist_Beta_CandCosDecayAngle = new TH2D((title+"_Beta_CandCosDecayAngle").c_str(), (title+"_Beta_CandCosDecayAngle;Z* candidate #beta;Z* candidate cos#theta").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_Beta_CandCosDecayAngle);
    TH2D* hist_BetaCM_CandCosDecayAngle = new TH2D((title+"_BetaCM_CandCosDecayAngle").c_str(), (title+"_BetaCM_CandCosDecayAngle;Z* candidate #beta^{CM};Z* candidate cos#theta").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_BetaCM_CandCosDecayAngle);
    TH2D* hist_CandCosDecayAngle_CandDeltaPhiMET = new TH2D((title+"_CandCosDecayAngle_CandDeltaPhiMET").c_str(), (title+"_CandCosDecayAngle_CandDeltaPhiMET;Z* candidate cos#theta;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandCosDecayAngle_CandDeltaPhiMET);
    TH2D* hist_CandML_CandCosDecayAngleCM = new TH2D((title+"_CandML_CandCosDecayAngleCM").c_str(), (title+"_CandML_CandCosDecayAngleCM;Z* candidate M [GeV];Z* candidate cos#theta CM").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_CandCosDecayAngleCM);
    TH2D* hist_Beta_CandCosDecayAngleCM = new TH2D((title+"_Beta_CandCosDecayAngleCM").c_str(), (title+"_Beta_CandCosDecayAngleCM;Z* candidate #beta;Z* candidate cos#theta CM").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_Beta_CandCosDecayAngleCM);
    TH2D* hist_BetaCM_CandCosDecayAngleCM = new TH2D((title+"_BetaCM_CandCosDecayAngleCM").c_str(), (title+"_BetaCM_CandCosDecayAngleCM;Z* candidate #beta^{CM};Z* candidate cos#theta CM").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_BetaCM_CandCosDecayAngleCM);
    TH2D* hist_CandCosDecayAngleCM_CandDeltaPhiMET = new TH2D((title+"_CandCosDecayAngleCM_CandDeltaPhiMET").c_str(), (title+"_CandCosDecayAngleCM_CandDeltaPhiMET;Z* candidate cos#theta CM;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_CandCosDecayAngleCM_CandDeltaPhiMET);
    TH2D* hist_Beta_CandDeltaPhiMET = new TH2D((title+"_Beta_CandDeltaPhiMET").c_str(), (title+"_Beta_CandDeltaPhiMET;Z* candidate #beta;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_Beta_CandDeltaPhiMET);
    TH2D* hist_BetaCM_CandDeltaPhiMET = new TH2D((title+"_BetaCM_CandDeltaPhiMET").c_str(), (title+"_BetaCM_CandDeltaPhiMET;Z* candidate #beta^{CM};#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_BetaCM_CandDeltaPhiMET);
    TH2D* hist_CandML_RZPara = new TH2D((title+"_CandML_RZPara").c_str(), (title+"_CandML_RZPara;Z* candidate M [GeV];RZPara").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandML_RZPara);
    TH2D* hist_CandML_RZPerp = new TH2D((title+"_CandML_RZPerp").c_str(), (title+"_CandML_RZPerp;Z* candidate M [GeV];RZPerp").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandML_RZPerp);
    TH2D* hist_RZPara_RZPerp = new TH2D((title+"_RZPara_RZPerp").c_str(), (title+"_RZPara_RZPerp;RZPara;RZPerp").c_str(), g_NX/2., -0.5, 0.5, g_NX/2, 0., 0.5);
    hists2.push_back(hist_RZPara_RZPerp);
    TH2D* hist_CandML_PZAng = new TH2D((title+"_CandML_PZAng").c_str(), (title+"_CandML_PZAng;Z* candidate M [GeV];PZAng").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandML_PZAng);
    TH2D* hist_CandML_PZPara = new TH2D((title+"_CandML_PZPara").c_str(), (title+"_CandML_PZPara;Z* candidate M [GeV];PZPara").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_PZPara);
    TH2D* hist_CandML_PZPerp = new TH2D((title+"_CandML_PZPerp").c_str(), (title+"_CandML_PZPerp;Z* candidate M [GeV];PZPerp").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_PZPerp);

    TH2D* hist_CandML_RZParaLABMET = new TH2D((title+"_CandML_RZParaLABMET").c_str(), (title+"_CandML_RZParaLABMET;Z* candidate M [GeV];RZParaLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandML_RZParaLABMET);
    TH2D* hist_CandML_RZPerpLABMET = new TH2D((title+"_CandML_RZPerpLABMET").c_str(), (title+"_CandML_RZPerpLABMET;Z* candidate M [GeV];RZPerpLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandML_RZPerpLABMET);
    TH2D* hist_CandML_PZAngLABMET = new TH2D((title+"_CandML_PZAngLABMET").c_str(), (title+"_CandML_PZAngLABMET;Z* candidate M [GeV];PZAngLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandML_PZAngLABMET);
    TH2D* hist_RZParaLABMET_RZPerpLABMET = new TH2D((title+"_RZParaLABMET_RZPerpLABMET").c_str(), (title+"_RZParaLABMET_RZPerpLABMET;RZParaLABMET;RZPerpLABMET").c_str(), g_NX/2., -0.5, 0.5, g_NX/2, 0., 0.5);
    hists2.push_back(hist_RZParaLABMET_RZPerpLABMET);
    TH2D* hist_CandML_PZParaLABMET = new TH2D((title+"_CandML_PZParaLABMET").c_str(), (title+"_CandML_PZParaLABMET;Z* candidate M [GeV];PZParaLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 400.);
    hists2.push_back(hist_CandML_PZParaLABMET);
    TH2D* hist_CandML_PZPerpLABMET = new TH2D((title+"_CandML_PZPerpLABMET").c_str(), (title+"_CandML_PZPerpLABMET;Z* candidate M [GeV];PZPerpLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 300.);
    hists2.push_back(hist_CandML_PZPerpLABMET);

    TH2D* hist_RZPara_RZParaLABMET = new TH2D((title+"_RZPara_RZParaLABMET").c_str(), (title+"_RZPara_RZParaLABMET;RZPara;RZParaLABMET").c_str(), g_NX/2., -0.75, 0.75, g_NX/2., -0.75, 0.75);
    hists2.push_back(hist_RZPara_RZParaLABMET);
    TH2D* hist_RZPerp_RZPerpLABMET = new TH2D((title+"_RZPerp_RZPerpLABMET").c_str(), (title+"_RZPerp_RZPerpLABMET;RZPerp;RZPerpLABMET").c_str(), g_NX/2., 0., 0.75, g_NX/2., 0., 0.75);
    hists2.push_back(hist_RZPerp_RZPerpLABMET);
    TH2D* hist_PZPara_PZParaLABMET = new TH2D((title+"_PZPara_PZParaLABMET").c_str(), (title+"_PZPara_PZParaLABMET;PZPara;PZParaLABMET").c_str(), g_NX/2., 0., 400., g_NX/2., 0., 400.);
    hists2.push_back(hist_PZPara_PZParaLABMET);
    TH2D* hist_PZPerp_PZPerpLABMET = new TH2D((title+"_PZPerp_PZPerpLABMET").c_str(), (title+"_PZPerp_PZPerpLABMET;PZPerp;PZPerpLABMET").c_str(), g_NX/2., 0., 400., g_NX/2., 0., 400.);
    hists2.push_back(hist_PZPerp_PZPerpLABMET);

    TH2D* hist_CandML_BetaZPara = new TH2D((title+"_CandML_BetaZPara").c_str(), (title+"_CandML_BetaZPara;Z* candidate M [GeV];PZPara/E^{CM}").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, -1., 1.);
    hists2.push_back(hist_CandML_BetaZPara);
    TH2D* hist_CandML_BetaZPerp = new TH2D((title+"_CandML_BetaZPerp").c_str(), (title+"_CandML_BetaZPerp;Z* candidate M [GeV];PZPerp/E^{CM}").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_CandML_BetaZPerp);
    TH2D* hist_CandCosDecayAngle_BetaZPara = new TH2D((title+"_CandCosDecayAngle_BetaZPara").c_str(), (title+"_CandCosDecayAngle_BetaZPara;Z* cos#theta;PZPara/E^{CM}").c_str(), g_NX/2., 0., 1., g_NX/2, -1., 1.);
    hists2.push_back(hist_CandCosDecayAngle_BetaZPara);
    TH2D* hist_CandCosDecayAngle_BetaZPerp = new TH2D((title+"_CandCosDecayAngle_BetaZPerp").c_str(), (title+"_CandCosDecayAngle_BetaZPerp;Z* cos#theta;PZPerp/E^{CM}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_CandCosDecayAngle_BetaZPerp);
    TH2D* hist_CandCosDecayAngleCM_BetaZPara = new TH2D((title+"_CandCosDecayAngleCM_BetaZPara").c_str(), (title+"_CandCosDecayAngleCM_BetaZPara;Z* cos#theta CM;PZPara/E^{CM}").c_str(), g_NX/2., 0., 1., g_NX/2, -1., 1.);
    hists2.push_back(hist_CandCosDecayAngleCM_BetaZPara);
    TH2D* hist_CandCosDecayAngleCM_BetaZPerp = new TH2D((title+"_CandCosDecayAngleCM_BetaZPerp").c_str(), (title+"_CandCosDecayAngleCM_BetaZPerp;Z* cos#theta CM;PZPerp/E^{CM}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_CandCosDecayAngleCM_BetaZPerp);

    TH2D* hist_CandML_MRZPara = new TH2D((title+"_CandML_MRZPara").c_str(), (title+"_CandML_MRZPara;Z* candidate M [GeV];M/RZPara").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_MRZPara);
    TH2D* hist_CandML_MRZPerp = new TH2D((title+"_CandML_MRZPerp").c_str(), (title+"_CandML_MRZPerp;Z* candidate M [GeV];M/RZPerp").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_MRZPerp);
    TH2D* hist_CandML_MRZParaLABMET = new TH2D((title+"_CandML_MRZParaLABMET").c_str(), (title+"_CandML_MRZParaLABMET;Z* candidate M [GeV];M/RZParaLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_MRZParaLABMET);
    TH2D* hist_CandML_MRZPerpLABMET = new TH2D((title+"_CandML_MRZPerpLABMET").c_str(), (title+"_CandML_MRZPerpLABMET;Z* candidate M [GeV];M/RZPerpLABMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 500.);
    hists2.push_back(hist_CandML_MRZPerpLABMET);

    TH2D* hist_CandCosDecayAngle_MPZPara = new TH2D((title+"_CandCosDecayAngle_MPZPara").c_str(), (title+"_CandCosDecayAngle_MPZPara;Z* cos#theta;PZPara/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngle_MPZPara);
    TH2D* hist_CandCosDecayAngle_MPZPerp = new TH2D((title+"_CandCosDecayAngle_MPZPerp").c_str(), (title+"_CandCosDecayAngle_MPZPerp;Z* cos#theta;PZPerp/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngle_MPZPerp);
    TH2D* hist_CandCosDecayAngle_MPZParaLABMET = new TH2D((title+"_CandCosDecayAngle_MPZParaLABMET").c_str(), (title+"_CandCosDecayAngle_MPZParaLABMET;Z* cos#theta;PZParaLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngle_MPZParaLABMET);
    TH2D* hist_CandCosDecayAngle_MPZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_MPZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_MPZPerpLABMET;Z* cos#theta;PZPerpLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngle_MPZPerpLABMET);

    TH2D* hist_CandCosDecayAngleCM_MPZPara = new TH2D((title+"_CandCosDecayAngleCM_MPZPara").c_str(), (title+"_CandCosDecayAngleCM_MPZPara;Z* cos#theta CM;PZPara/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngleCM_MPZPara);
    TH2D* hist_CandCosDecayAngleCM_MPZPerp = new TH2D((title+"_CandCosDecayAngleCM_MPZPerp").c_str(), (title+"_CandCosDecayAngleCM_MPZPerp;Z* cos#theta CM;PZPerp/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngleCM_MPZPerp);
    TH2D* hist_CandCosDecayAngleCM_MPZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_MPZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_MPZParaLABMET;Z* cos#theta CM;PZParaLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngleCM_MPZParaLABMET);
    TH2D* hist_CandCosDecayAngleCM_MPZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_MPZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_MPZPerpLABMET;Z* cos#theta CM;PZPerpLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_CandCosDecayAngleCM_MPZPerpLABMET);

    TH2D* hist_Beta_MPZPara = new TH2D((title+"_Beta_MPZPara").c_str(), (title+"_Beta_MPZPara;Z* #beta;PZPara/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_Beta_MPZPara);
    TH2D* hist_Beta_MPZPerp = new TH2D((title+"_Beta_MPZPerp").c_str(), (title+"_Beta_MPZPerp;Z* #beta;PZPerp/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_Beta_MPZPerp);
    TH2D* hist_Beta_MPZParaLABMET = new TH2D((title+"_Beta_MPZParaLABMET").c_str(), (title+"_Beta_MPZParaLABMET;Z* #beta;PZParaLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_Beta_MPZParaLABMET);
    TH2D* hist_Beta_MPZPerpLABMET = new TH2D((title+"_Beta_MPZPerpLABMET").c_str(), (title+"_Beta_MPZPerpLABMET;Z* #beta;PZPerpLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_Beta_MPZPerpLABMET);

    TH2D* hist_BetaCM_MPZPara = new TH2D((title+"_BetaCM_MPZPara").c_str(), (title+"_BetaCM_MPZPara;Z* #beta^{CM};PZPara/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_BetaCM_MPZPara);
    TH2D* hist_BetaCM_MPZPerp = new TH2D((title+"_BetaCM_MPZPerp").c_str(), (title+"_BetaCM_MPZPerp;Z* #beta^{CM};PZPerp/M").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_BetaCM_MPZPerp);
    TH2D* hist_BetaCM_MPZParaLABMET = new TH2D((title+"_BetaCM_MPZParaLABMET").c_str(), (title+"_BetaCM_MPZParaLABMET;Z* #beta^{CM};PZParaLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_BetaCM_MPZParaLABMET);
    TH2D* hist_BetaCM_MPZPerpLABMET = new TH2D((title+"_BetaCM_MPZPerpLABMET").c_str(), (title+"_BetaCM_MPZPerpLABMET;Z* #beta^{CM};PZPerpLABMET/M_{T}").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 2.);
    hists2.push_back(hist_BetaCM_MPZPerpLABMET);

    TH2D* hist_PZAng_MPZPara = new TH2D((title+"_PZAng_MPZPara").c_str(), (title+"_PZAng_MPZPara;PZAng;PZPara/M").c_str(), g_NX/2., -1.6, 1.6, g_NX/2, 0., 2.);
    hists2.push_back(hist_PZAng_MPZPara);
    TH2D* hist_PZAng_MPZPerp = new TH2D((title+"_PZAng_MPZPerp").c_str(), (title+"_PZAng_MPZPerp;PZAng;PZPerp/M").c_str(), g_NX/2., -1.6, 1.6, g_NX/2, 0., 2.);
    hists2.push_back(hist_PZAng_MPZPerp);
    TH2D* hist_PZAng_MPZParaLABMET = new TH2D((title+"_PZAng_MPZParaLABMET").c_str(), (title+"_PZAng_MPZParaLABMET;PZAng;PZParaLABMET/M_{T}").c_str(), g_NX/2., -1.6, 1.6, g_NX/2, 0., 2.);
    hists2.push_back(hist_PZAng_MPZParaLABMET);
    TH2D* hist_PZAng_MPZPerpLABMET = new TH2D((title+"_PZAng_MPZPerpLABMET").c_str(), (title+"_PZAng_MPZPerpLABMET;PZAng;PZPerpLABMET/M_{T}").c_str(), g_NX/2., -1.6, 1.6, g_NX/2, 0., 2.);
    hists2.push_back(hist_PZAng_MPZPerpLABMET);

    TH2D* hist_PZPara_CandECM = new TH2D((title+"_PZPara_CandECM").c_str(), (title+"_PZPara_CandECM;PZPara;Z* candidate E^{CM} [GeV]").c_str(), g_NX/2., 0., 400., g_NX/2., 0., 400.);
    hists2.push_back(hist_PZPara_CandECM);
    TH2D* hist_PZPerp_CandECM = new TH2D((title+"_PZPerp_CandECM").c_str(), (title+"_PZPerp_CandECM;PZPerp;Z* candidate E^{CM} [GeV]").c_str(), g_NX/2., 0., 400., g_NX/2., 0., 400.);
    hists2.push_back(hist_PZPerp_CandECM);

    TH2D* hist_CandCosDecayAngle_RZPara = new TH2D((title+"_CandCosDecayAngle_RZPara").c_str(), (title+"_CandCosDecayAngle_RZPara;Z* candidate cos#theta;RZPara").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandCosDecayAngle_RZPara);
    TH2D* hist_CandCosDecayAngle_RZPerp = new TH2D((title+"_CandCosDecayAngle_RZPerp").c_str(), (title+"_CandCosDecayAngle_RZPerp;Z* candidate cos#theta;RZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandCosDecayAngle_RZPerp);
    TH2D* hist_CandCosDecayAngle_PZAng = new TH2D((title+"_CandCosDecayAngle_PZAng").c_str(), (title+"_CandCosDecayAngle_PZAng;Z* candidate cos#theta;PZAng").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandCosDecayAngle_PZAng);
    TH2D* hist_CandCosDecayAngle_PZPara = new TH2D((title+"_CandCosDecayAngle_PZPara").c_str(), (title+"_CandCosDecayAngle_PZPara;Z* candidate cos#theta;PZPara").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngle_PZPara);
    TH2D* hist_CandCosDecayAngle_PZPerp = new TH2D((title+"_CandCosDecayAngle_PZPerp").c_str(), (title+"_CandCosDecayAngle_PZPerp;Z* candidate cos#theta;PZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngle_PZPerp);

    TH2D* hist_CandCosDecayAngle_RZParaLABMET = new TH2D((title+"_CandCosDecayAngle_RZParaLABMET").c_str(), (title+"_CandCosDecayAngle_RZParaLABMET;Z* candidate cos#theta;RZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandCosDecayAngle_RZParaLABMET);
    TH2D* hist_CandCosDecayAngle_RZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_RZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_RZPerpLABMET;Z* candidate cos#theta;RZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandCosDecayAngle_RZPerpLABMET);
    TH2D* hist_CandCosDecayAngle_PZAngLABMET = new TH2D((title+"_CandCosDecayAngle_PZAngLABMET").c_str(), (title+"_CandCosDecayAngle_PZAngLABMET;Z* candidate cos#theta;PZAngLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandCosDecayAngle_PZAngLABMET);
    TH2D* hist_CandCosDecayAngle_PZParaLABMET = new TH2D((title+"_CandCosDecayAngle_PZParaLABMET").c_str(), (title+"_CandCosDecayAngle_PZParaLABMET;Z* candidate cos#theta;PZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngle_PZParaLABMET);
    TH2D* hist_CandCosDecayAngle_PZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_PZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_PZPerpLABMET;Z* candidate cos#theta;PZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngle_PZPerpLABMET);

    TH2D* hist_CandCosDecayAngleCM_RZPara = new TH2D((title+"_CandCosDecayAngleCM_RZPara").c_str(), (title+"_CandCosDecayAngleCM_RZPara;Z* candidate cos#theta CM;RZPara").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandCosDecayAngleCM_RZPara);
    TH2D* hist_CandCosDecayAngleCM_RZPerp = new TH2D((title+"_CandCosDecayAngleCM_RZPerp").c_str(), (title+"_CandCosDecayAngleCM_RZPerp;Z* candidate cos#theta CM;RZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandCosDecayAngleCM_RZPerp);
    TH2D* hist_CandCosDecayAngleCM_PZAng = new TH2D((title+"_CandCosDecayAngleCM_PZAng").c_str(), (title+"_CandCosDecayAngleCM_PZAng;Z* candidate cos#theta CM;PZAng").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandCosDecayAngleCM_PZAng);
    TH2D* hist_CandCosDecayAngleCM_PZPara = new TH2D((title+"_CandCosDecayAngleCM_PZPara").c_str(), (title+"_CandCosDecayAngleCM_PZPara;Z* candidate cos#theta CM;PZPara").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngleCM_PZPara);
    TH2D* hist_CandCosDecayAngleCM_PZPerp = new TH2D((title+"_CandCosDecayAngleCM_PZPerp").c_str(), (title+"_CandCosDecayAngleCM_PZPerp;Z* candidate cos#theta CM;PZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngleCM_PZPerp);

    TH2D* hist_CandCosDecayAngleCM_RZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_RZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_RZParaLABMET;Z* candidate cos#theta CM;RZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_CandCosDecayAngleCM_RZParaLABMET);
    TH2D* hist_CandCosDecayAngleCM_RZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_RZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_RZPerpLABMET;Z* candidate cos#theta CM;RZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_CandCosDecayAngleCM_RZPerpLABMET);
    TH2D* hist_CandCosDecayAngleCM_PZAngLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZAngLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZAngLABMET;Z* candidate cos#theta CM;PZAngLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_CandCosDecayAngleCM_PZAngLABMET);
    TH2D* hist_CandCosDecayAngleCM_PZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZParaLABMET;Z* candidate cos#theta CM;PZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngleCM_PZParaLABMET);
    TH2D* hist_CandCosDecayAngleCM_PZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZPerpLABMET;Z* candidate cos#theta CM;PZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_CandCosDecayAngleCM_PZPerpLABMET);

    TH2D* hist_CandCosDecayAngle_CandCosDecayAngleCM = new TH2D((title+"_CandCosDecayAngle_CandCosDecayAngleCM").c_str(), (title+"_CandCosDecayAngle_CandCosDecayAngleCM;Z* candidate cos#theta;Z* candidate cos#theta CM").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 1.);
    hists2.push_back(hist_CandCosDecayAngle_CandCosDecayAngleCM);
    TH2D* hist_CandBeta_CandBetaCM = new TH2D((title+"_CandBeta_CandBetaCM").c_str(), (title+"_CandBeta_CandBetaCM;Z* candidate #beta;Z* candidate #beta^{CM}").c_str(), g_NX/2., 0., 1., g_NX/2., 0., 1.);
    hists2.push_back(hist_CandBeta_CandBetaCM);

    TH2D* hist_Beta_RZPara = new TH2D((title+"_Beta_RZPara").c_str(), (title+"_Beta_RZPara;Z* candidate #beta;RZPara").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_Beta_RZPara);
    TH2D* hist_Beta_RZPerp = new TH2D((title+"_Beta_RZPerp").c_str(), (title+"_Beta_RZPerp;Z* candidate #beta;RZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_Beta_RZPerp);
    TH2D* hist_Beta_PZAng = new TH2D((title+"_Beta_PZAng").c_str(), (title+"_Beta_PZAng;Z* candidate #beta;PZAng").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_Beta_PZAng);
    TH2D* hist_Beta_PZPara = new TH2D((title+"_Beta_PZPara").c_str(), (title+"_Beta_PZPara;Z* candidate #beta;PZPara").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_Beta_PZPara);
    TH2D* hist_Beta_PZPerp = new TH2D((title+"_Beta_PZPerp").c_str(), (title+"_Beta_PZPerp;Z* candidate #beta;PZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_Beta_PZPerp);

    TH2D* hist_Beta_RZParaLABMET = new TH2D((title+"_Beta_RZParaLABMET").c_str(), (title+"_Beta_RZParaLABMET;Z* candidate #beta;RZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_Beta_RZParaLABMET);
    TH2D* hist_Beta_RZPerpLABMET = new TH2D((title+"_Beta_RZPerpLABMET").c_str(), (title+"_Beta_RZPerpLABMET;Z* candidate #beta;RZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_Beta_RZPerpLABMET);
    TH2D* hist_Beta_PZAngLABMET = new TH2D((title+"_Beta_PZAngLABMET").c_str(), (title+"_Beta_PZAngLABMET;Z* candidate #beta;PZAngLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_Beta_PZAngLABMET);
    TH2D* hist_Beta_PZParaLABMET = new TH2D((title+"_Beta_PZParaLABMET").c_str(), (title+"_Beta_PZParaLABMET;Z* candidate #beta;PZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_Beta_PZParaLABMET);
    TH2D* hist_Beta_PZPerpLABMET = new TH2D((title+"_Beta_PZPerpLABMET").c_str(), (title+"_Beta_PZPerpLABMET;Z* candidate #beta;PZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_Beta_PZPerpLABMET);

    TH2D* hist_BetaCM_RZPara = new TH2D((title+"_BetaCM_RZPara").c_str(), (title+"_BetaCM_RZPara;Z* candidate #beta^{CM};RZPara").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_BetaCM_RZPara);
    TH2D* hist_BetaCM_RZPerp = new TH2D((title+"_BetaCM_RZPerp").c_str(), (title+"_BetaCM_RZPerp;Z* candidate #beta^{CM};RZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_BetaCM_RZPerp);
    TH2D* hist_BetaCM_PZAng = new TH2D((title+"_BetaCM_PZAng").c_str(), (title+"_BetaCM_PZAng;Z* candidate #beta^{CM};PZAng").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_BetaCM_PZAng);
    TH2D* hist_BetaCM_PZPara = new TH2D((title+"_BetaCM_PZPara").c_str(), (title+"_BetaCM_PZPara;Z* candidate #beta^{CM};PZPara").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_BetaCM_PZPara);
    TH2D* hist_BetaCM_PZPerp = new TH2D((title+"_BetaCM_PZPerp").c_str(), (title+"_BetaCM_PZPerp;Z* candidate #beta^{CM};PZPerp").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_BetaCM_PZPerp);

    TH2D* hist_BetaCM_RZParaLABMET = new TH2D((title+"_BetaCM_RZParaLABMET").c_str(), (title+"_BetaCM_RZParaLABMET;Z* candidate #beta^{CM};RZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -0.5, 0.5);
    hists2.push_back(hist_BetaCM_RZParaLABMET);
    TH2D* hist_BetaCM_RZPerpLABMET = new TH2D((title+"_BetaCM_RZPerpLABMET").c_str(), (title+"_BetaCM_RZPerpLABMET;Z* candidate #beta^{CM};RZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 0.5);
    hists2.push_back(hist_BetaCM_RZPerpLABMET);
    TH2D* hist_BetaCM_PZAngLABMET = new TH2D((title+"_BetaCM_PZAngLABMET").c_str(), (title+"_BetaCM_PZAngLABMET;Z* candidate #beta^{CM};PZAngLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, -1.6, 1.6);
    hists2.push_back(hist_BetaCM_PZAngLABMET);
    TH2D* hist_BetaCM_PZParaLABMET = new TH2D((title+"_BetaCM_PZParaLABMET").c_str(), (title+"_BetaCM_PZParaLABMET;Z* candidate #beta^{CM};PZParaLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_BetaCM_PZParaLABMET);
    TH2D* hist_BetaCM_PZPerpLABMET = new TH2D((title+"_BetaCM_PZPerpLABMET").c_str(), (title+"_BetaCM_PZPerpLABMET;Z* candidate #beta^{CM};PZPerpLABMET").c_str(), g_NX/2., 0., 1., g_NX/2, 0., 400.);
    hists2.push_back(hist_BetaCM_PZPerpLABMET);

    TH2D* hist_RISR_MatchedCandML = new TH2D((title+"_RISR_MatchedCandML").c_str(), (title+"_RISR_MatchedCandML;R_{ISR};Matched Z* candidate M [GeV]").c_str(), g_NX/2., 0.5, 1., g_NX/2, CandMinMass, CandMaxMass);
    hists2.push_back(hist_RISR_MatchedCandML);
    TH2D* hist_RISR_MatchedCandBeta = new TH2D((title+"_RISR_MatchedCandBeta").c_str(), (title+"_RISR_MatchedCandBeta;R_{ISR};Matched #beta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_MatchedCandBeta);
    TH2D* hist_RISR_MatchedCandBetaCM = new TH2D((title+"_RISR_MatchedCandBetaCM").c_str(), (title+"_RISR_MatchedCandBetaCM;R_{ISR};Matched #beta^{CM}").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_MatchedCandBetaCM);
    TH2D* hist_RISR_MatchedCandDeltaPhiMET = new TH2D((title+"_RISR_MatchedCandDeltaPhiMET").c_str(), (title+"_RISR_MatchedCandDeltaPhiMET;R_{ISR};MatchedCandDeltaPhiMET").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_MatchedCandDeltaPhiMET);
    TH2D* hist_MatchedCandML_MatchedCandBeta = new TH2D((title+"_MatchedCandML_MatchedCandBeta").c_str(), (title+"_MatchedCandML_MatchedCandBeta;MatchedCandML;Matched #beta").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandBeta);
    TH2D* hist_MatchedCandML_MatchedCandBetaCM = new TH2D((title+"_MatchedCandML_MatchedCandBetaCM").c_str(), (title+"_MatchedCandML_MatchedCandBetaCM;MatchedCandML;Matched #beta^{CM}").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandBetaCM);
    TH2D* hist_MatchedCandML_MatchedCandDeltaPhiMET = new TH2D((title+"_MatchedCandML_MatchedCandDeltaPhiMET").c_str(), (title+"_MatchedCandML_MatchedCandDeltaPhiMET;MatchedCandML;MatchedCandDeltaPhiMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 3.15);
    hists2.push_back(hist_MatchedCandML_MatchedCandDeltaPhiMET);
    TH2D* hist_MatchedCandML_MatchedCandCosDecayAngle = new TH2D((title+"_MatchedCandML_MatchedCandCosDecayAngle").c_str(), (title+"_MatchedCandML_MatchedCandCosDecayAngle;MatchedCandML;MatchedCandCosDecayAngle").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandCosDecayAngle);
    TH2D* hist_MatchedCandML_MatchedCandCosDecayAngleCM = new TH2D((title+"_MatchedCandML_MatchedCandCosDecayAngleCM").c_str(), (title+"_MatchedCandML_MatchedCandCosDecayAngleCM;MatchedCandML;MatchedCandCosDecayAngleCM").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_MatchedCandML_MatchedCandCosDecayAngleCM);

    TH2D* hist_RISR_UnmatchedCandML = new TH2D((title+"_RISR_UnmatchedCandML").c_str(), (title+"_RISR_UnmatchedCandML;R_{ISR};UnmatchedCandML").c_str(), g_NX/2., 0.5, 1., g_NX/2, CandMinMass, CandMaxMass);
    hists2.push_back(hist_RISR_UnmatchedCandML);
    TH2D* hist_RISR_UnmatchedCandBeta = new TH2D((title+"_RISR_UnmatchedCandBeta").c_str(), (title+"_RISR_UnmatchedCandBeta;R_{ISR};Unmatched #beta").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_UnmatchedCandBeta);
    TH2D* hist_UnmatchedCandML_UnmatchedCandBeta = new TH2D((title+"_UnmatchedCandML_UnmatchedCandBeta").c_str(), (title+"_UnmatchedCandML_UnmatchedCandBeta;UnmatchedCandML;UnmatchedCandBeta").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandBeta);
    TH2D* hist_RISR_UnmatchedCandBetaCM = new TH2D((title+"_RISR_UnmatchedCandBetaCM").c_str(), (title+"_RISR_UnmatchedCandBetaCM;R_{ISR};Unmatched #beta^{CM}").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 1.);
    hists2.push_back(hist_RISR_UnmatchedCandBetaCM);
    TH2D* hist_UnmatchedCandML_UnmatchedCandBetaCM = new TH2D((title+"_UnmatchedCandML_UnmatchedCandBetaCM").c_str(), (title+"_UnmatchedCandML_UnmatchedCandBetaCM;UnmatchedCandML;UnmatchedCandBetaCM").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandBetaCM);
    TH2D* hist_RISR_UnmatchedCandDeltaPhiMET = new TH2D((title+"_RISR_UnmatchedCandDeltaPhiMET").c_str(), (title+"_RISR_UnmatchedCandDeltaPhiMET;R_{ISR};UnmatchedCandDeltaPhiMET").c_str(), g_NX/2., 0.5, 1., g_NX/2, 0., 3.15);
    hists2.push_back(hist_RISR_UnmatchedCandDeltaPhiMET);
    TH2D* hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET = new TH2D((title+"_UnmatchedCandML_UnmatchedCandDeltaPhiMET").c_str(), (title+"_UnmatchedCandML_UnmatchedCandDeltaPhiMET;UnmatchedCandML;UnmatchedCandDeltaPhiMET").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 3.15);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET);
    TH2D* hist_UnmatchedCandML_UnmatchedCandCosDecayAngle = new TH2D((title+"_UnmatchedCandML_UnmatchedCandCosDecayAngle").c_str(), (title+"_UnmatchedCandML_UnmatchedCandCosDecayAngle;UnmatchedCandML;UnmatchedCandCosDecayAngle").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandCosDecayAngle);
    TH2D* hist_UnmatchedCandML_UnmatchedCandCosDecayAngleCM = new TH2D((title+"_UnmatchedCandML_UnmatchedCandCosDecayAngleCM").c_str(), (title+"_UnmatchedCandML_UnmatchedCandCosDecayAngleCM;UnmatchedCandML;UnmatchedCandCosDecayAngleCM").c_str(), g_NX/2., CandMinMass, CandMaxMass, g_NX/2, 0., 1.);
    hists2.push_back(hist_UnmatchedCandML_UnmatchedCandCosDecayAngleCM);

    TH2D* hist_DecayAngleDiff_PzCM = new TH2D((title+"_DecayAngleDiff_PzCM").c_str(), (title+"_DecayAngleDiff_PzCM;|cos#theta CM| - |cos#theta|;PzCM").c_str(), g_NX/2., -1., 1., g_NX/2, 0., 500.);
    hists2.push_back(hist_DecayAngleDiff_PzCM);
    TH2D* hist_BetaDiff_PzCM = new TH2D((title+"_BetaDiff_PzCM").c_str(), (title+"_BetaDiff_PzCM;#beta CM - #beta;PzCM").c_str(), g_NX/2., -1., 1., g_NX/2, 0., 500.);
    hists2.push_back(hist_BetaDiff_PzCM);
    TH2D* hist_DecayAngleDiff_BetaDiff = new TH2D((title+"_DecayAngleDiff_BetaDiff").c_str(), (title+"_DecayAngleDiff_BetaDiff;|cos#theta CM| - |cos#theta|;#beta CM - #beta").c_str(), g_NX/2., -1., 1., g_NX/2, -1., 1.);
    hists2.push_back(hist_DecayAngleDiff_BetaDiff);

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
          hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "#Delta#phi(MET, V) < #frac{#pi}{2}");

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
          //if(NbjetISR + NbjetS > 1) continue; // SR & ATLAS 'B-Veto'
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "B-Veto");

          //if(Nlep != 2) continue;
          //if(NjetS != 0) continue; // SR
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //hist_CutFlow->GetXaxis()->SetBinLabel(CF_bin, "0S J");

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
          
          TVector3 TV3_MET;
          TV3_MET.SetPtEtaPhi(MET, 0., base->MET_phi);
          TLorentzVector TLV_LAB, TLV_CM, TLV_S;
          TLV_LAB.SetPtEtaPhiM(base->LAB_Pt, base->LAB_Eta, base->LAB_Phi, base->LAB_M);
          TLV_CM.SetPtEtaPhiM(base->PTCM, base->EtaCM, base->PhiCM, base->MCM);
          TLV_S.SetPtEtaPhiM(base->PTS, base->EtaS, base->PhiS, base->MS);
          double PzCM = base->PzCM;
          TVector3 beta_CM = TLV_CM.BoostVector();
          TLorentzVector S_CM = TLV_S;
          S_CM.Boost(-beta_CM); // S in CM frame

          // Leptonic Candidates
          std::vector<L_Cand> V_lep_cands;
          for(int i = 0; i < Nlep-1; i++){
            for(int j = i+1; j < Nlep; j++){
              if(list_leps[i].Flavor() == list_leps[j].Flavor() && list_leps[i].Charge() != list_leps[j].Charge()){ // OSSF
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
                lep_1.SetCharge((base->Charge_lep->at(i) > 0 ? 1 : -1));
                lep_2.SetCharge((base->Charge_lep->at(j) > 0 ? 1 : -1));
                ParticleList cand_list_parts;
                cand_list_parts.push_back(lep_1);
                cand_list_parts.push_back(lep_2);
                L_Cand cand(cand_list_parts);
                if(cand.M() < CandMinMass || cand.M() > CandMaxMass) continue; 

                TLorentzVector TLV_Cand = cand.TLV();
                TLorentzVector TLV_Cand_CM = TLV_Cand;
                TLV_Cand_CM.Boost(-beta_CM);
                TVector3 beta_Cand_CM = TLV_Cand_CM.BoostVector();
                Particle part_Cand_Child = cand.Cand_PartPlus();
                TLorentzVector TLV_Cand_Child;
                TLV_Cand_Child.SetPtEtaPhiM(part_Cand_Child.Pt(), part_Cand_Child.Eta(), part_Cand_Child.Phi(), part_Cand_Child.M());
                TLorentzVector TLV_Cand_Child_CM = TLV_Cand_Child;
                TLV_Cand_Child_CM.Boost(-beta_CM);
                TLV_Cand_Child_CM.Boost(-beta_Cand_CM);
                double CosDecayAngleCM = fabs(TLV_Cand_Child_CM.Vect().Unit().Dot(beta_Cand_CM.Unit())); // 'original' using CM
                if(CosDecayAngleCM > 0.8) continue;
                double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
                double PZPerp = TLV_Cand_CM.Vect().Cross(S_CM.Vect().Unit()).Mag();
                double PZAng = atan(PZPara/PZPerp);
                double MPZPara = fabs(PZPara/cand.M());

                cand.SetFlavor(list_leps[i].Flavor());
                bool i_a = inVec((*base->index_lep_a),i);
                bool j_a = inVec((*base->index_lep_a),j);
                if(i_a == j_a || (Nlep == 2 && NjetS == 0)){
                  cand.SetSameHemi(true);
                }
                else cand.SetSameHemi(false);
                if(CandSameHemi && !cand.IsSameHemi()) continue;
                V_lep_cands.push_back(cand);
              }
            }
          }

          int N_V_lep_cands = V_lep_cands.size();
          cand_matching(V_lep_cands);
          //if(N_V_lep_cands < 1) continue;

          for(int i = 0; i < N_V_lep_cands; i++){
            TLorentzVector TLV_Cand = V_lep_cands[i].TLV();
            TLorentzVector TLV_Cand_CM = TLV_Cand;
            TLV_Cand_CM.Boost(-beta_CM);
            TVector3 beta_Cand_CM = TLV_Cand_CM.BoostVector();
            Particle part_Cand_Child = V_lep_cands[i].Cand_PartPlus();
            TLorentzVector TLV_Cand_Child;
            TLV_Cand_Child.SetPtEtaPhiM(part_Cand_Child.Pt(), part_Cand_Child.Eta(), part_Cand_Child.Phi(), part_Cand_Child.M());
            TLorentzVector TLV_Cand_Child_CM = TLV_Cand_Child;
            TLV_Cand_Child_CM.Boost(-beta_CM);
            TLV_Cand_Child_CM.Boost(-beta_Cand_CM);
            double CosDecayAngleCM = fabs(TLV_Cand_Child_CM.Vect().Unit().Dot(beta_Cand_CM.Unit())); // 'original' using CM
            double CosDecayAngle = fabs(V_lep_cands[i].CosDecayAngle()); // decay angle without using CM
            double CandML = TLV_Cand.M();
            double Beta = V_lep_cands[i].Beta();
            double BetaCM = TLV_Cand_CM.P()/TLV_Cand_CM.E();
            double CandECM = TLV_Cand_CM.E();
            double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
            double PZPerp = TLV_Cand_CM.Vect().Cross(S_CM.Vect().Unit()).Mag();
            double RZPara = PZPara/S_CM.Vect().Mag();
            double RZPerp = PZPerp/S_CM.Vect().Mag();
            double PZAng = atan(PZPara/PZPerp);
            TLorentzVector TLV_Cand_Transverse = TLV_Cand;
            TLV_Cand_Transverse.SetZ(0.);
            TVector3 TV3_Cand_Transverse = TLV_Cand_Transverse.Vect();
            double PZParaLABMET = TV3_Cand_Transverse.Dot(TV3_MET.Unit());
            double PZPerpLABMET = TV3_Cand_Transverse.Cross(TV3_MET.Unit()).Mag();
            double RZParaLABMET = PZParaLABMET/MET;
            double RZPerpLABMET = PZPerpLABMET/MET;
            double PZAngLABMET = atan(PZParaLABMET/PZPerpLABMET);
            double MRZPara = fabs(CandML/RZPara);
            double MRZPerp = fabs(CandML/RZPerp);
            double MRZParaLABMET = fabs(CandML/RZParaLABMET);
            double MRZPerpLABMET = fabs(CandML/RZPerpLABMET);
            double MPZPara = fabs(PZPara/CandML);
            double MPZPerp = fabs(PZPerp/CandML);
            double MPZParaLABMET = fabs(PZParaLABMET/TLV_Cand_Transverse.Mag());
            double MPZPerpLABMET = fabs(PZPerpLABMET/TLV_Cand_Transverse.Mag());
            double BetaZPara = PZPara/CandECM;
            double BetaZPerp = PZPerp/CandECM;
            hist_CandML->Fill(CandML, weight);
            hist_CandBeta->Fill(Beta, weight);
            hist_CandBetaCM->Fill(BetaCM, weight);
            hist_CandBeta_CandBetaCM->Fill(Beta, BetaCM, weight);
            hist_CandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandCosDecayAngle->Fill(CosDecayAngle, weight);
            hist_CandCosDecayAngleCM->Fill(CosDecayAngleCM, weight);
            hist_RZPara->Fill(RZPara, weight);
            hist_RZPerp->Fill(RZPerp, weight);
            hist_RISR_CandML->Fill(RISR, CandML, weight);
            hist_RISR_CandBeta->Fill(RISR, Beta, weight);
            hist_CandML_CandBeta->Fill(CandML, Beta, weight);
            hist_Beta_CandCosDecayAngleCM->Fill(Beta, CosDecayAngleCM, weight);
            hist_Beta_CandCosDecayAngle->Fill(Beta, CosDecayAngle, weight);
            hist_Beta_CandDeltaPhiMET->Fill(Beta, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_RISR_CandBetaCM->Fill(RISR, BetaCM, weight);
            hist_CandML_CandBetaCM->Fill(CandML, BetaCM, weight);
            hist_BetaCM_CandCosDecayAngleCM->Fill(BetaCM, CosDecayAngleCM, weight);
            hist_BetaCM_CandCosDecayAngle->Fill(BetaCM, CosDecayAngle, weight);
            hist_BetaCM_CandDeltaPhiMET->Fill(BetaCM, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_RISR_CandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_CandDeltaPhiMET->Fill(CandML, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_CandCosDecayAngleCM->Fill(CandML, CosDecayAngleCM, weight);
            hist_CandCosDecayAngleCM_CandDeltaPhiMET->Fill(CosDecayAngleCM, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_CandCosDecayAngle->Fill(CandML, CosDecayAngle, weight);
            hist_CandCosDecayAngle_CandDeltaPhiMET->Fill(CosDecayAngle, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
            hist_CandML_RZPara->Fill(CandML, RZPara, weight);
            hist_CandML_RZPerp->Fill(CandML, RZPerp, weight);
            hist_PZAng->Fill(PZAng, weight);
            hist_PZPara->Fill(RZPara, weight);
            hist_PZPerp->Fill(RZPerp, weight);
            hist_CandML_PZAng->Fill(CandML, PZAng, weight);
            hist_RZPara_RZPerp->Fill(RZPara, RZPerp, weight);
            hist_CandML_PZPara->Fill(CandML, PZPara, weight);
            hist_CandML_PZPerp->Fill(CandML, PZPerp, weight);
            hist_CandML_RZParaLABMET->Fill(CandML, RZParaLABMET, weight);
            hist_CandML_RZPerpLABMET->Fill(CandML, RZPerpLABMET, weight);
            hist_PZAngLABMET->Fill(PZAngLABMET, weight);
            hist_CandML_PZAngLABMET->Fill(CandML, PZAngLABMET, weight);
            hist_RZParaLABMET_RZPerpLABMET->Fill(RZParaLABMET, RZPerpLABMET, weight);
            hist_CandML_PZParaLABMET->Fill(CandML, PZParaLABMET, weight);
            hist_CandML_PZPerpLABMET->Fill(CandML, PZPerpLABMET, weight);
            hist_CandCosDecayAngle_RZPara->Fill(CosDecayAngle, RZPara, weight);
            hist_CandCosDecayAngle_RZPerp->Fill(CosDecayAngle, RZPerp, weight); 
            hist_CandCosDecayAngle_PZAng->Fill(CosDecayAngle, PZAng, weight); 
            hist_CandCosDecayAngle_PZPara->Fill(CosDecayAngle, PZPara, weight); 
            hist_CandCosDecayAngle_PZPerp->Fill(CosDecayAngle, PZPara, weight); 
            hist_CandCosDecayAngle_RZParaLABMET->Fill(CosDecayAngle, RZParaLABMET, weight); 
            hist_CandCosDecayAngle_RZPerpLABMET->Fill(CosDecayAngle, RZPerpLABMET, weight); 
            hist_CandCosDecayAngle_PZAngLABMET->Fill(CosDecayAngle, PZAngLABMET, weight); 
            hist_CandCosDecayAngle_PZParaLABMET->Fill(CosDecayAngle, PZParaLABMET, weight); 
            hist_CandCosDecayAngle_PZPerpLABMET->Fill(CosDecayAngle, PZPerpLABMET, weight); 
            hist_CandCosDecayAngleCM_RZPara->Fill(CosDecayAngleCM, RZPara, weight);
            hist_CandCosDecayAngleCM_RZPerp->Fill(CosDecayAngleCM, RZPerp, weight); 
            hist_CandCosDecayAngleCM_PZAng->Fill(CosDecayAngleCM, PZAng, weight); 
            hist_CandCosDecayAngleCM_PZPara->Fill(CosDecayAngleCM, PZPara, weight); 
            hist_CandCosDecayAngleCM_PZPerp->Fill(CosDecayAngleCM, PZPara, weight); 
            hist_CandCosDecayAngleCM_RZParaLABMET->Fill(CosDecayAngleCM, RZParaLABMET, weight); 
            hist_CandCosDecayAngleCM_RZPerpLABMET->Fill(CosDecayAngleCM, RZPerpLABMET, weight); 
            hist_CandCosDecayAngleCM_PZAngLABMET->Fill(CosDecayAngleCM, PZAngLABMET, weight); 
            hist_CandCosDecayAngleCM_PZParaLABMET->Fill(CosDecayAngleCM, PZParaLABMET, weight); 
            hist_CandCosDecayAngleCM_PZPerpLABMET->Fill(CosDecayAngleCM, PZPerpLABMET, weight); 
            hist_Beta_RZPara->Fill(Beta, RZPara, weight); 
            hist_Beta_RZPerp->Fill(Beta, RZPerp, weight); 
            hist_Beta_PZAng->Fill(Beta, PZAng, weight); 
            hist_Beta_PZPara->Fill(Beta, PZPara, weight); 
            hist_Beta_PZPerp->Fill(Beta, PZPerp, weight); 
            hist_Beta_RZParaLABMET->Fill(Beta, RZParaLABMET, weight); 
            hist_Beta_RZPerpLABMET->Fill(Beta, RZPerpLABMET, weight); 
            hist_Beta_PZAngLABMET->Fill(Beta, PZAngLABMET, weight); 
            hist_Beta_PZParaLABMET->Fill(Beta, PZParaLABMET, weight); 
            hist_Beta_PZPerpLABMET->Fill(Beta, PZPerpLABMET, weight); 
            hist_BetaCM_RZPara->Fill(BetaCM, RZPara, weight); 
            hist_BetaCM_RZPerp->Fill(BetaCM, RZPerp, weight); 
            hist_BetaCM_PZAng->Fill(BetaCM, PZAng, weight); 
            hist_BetaCM_PZPara->Fill(BetaCM, PZPara, weight); 
            hist_BetaCM_PZPerp->Fill(BetaCM, PZPerp, weight); 
            hist_BetaCM_RZParaLABMET->Fill(BetaCM, RZParaLABMET, weight); 
            hist_BetaCM_RZPerpLABMET->Fill(BetaCM, RZPerpLABMET, weight); 
            hist_BetaCM_PZAngLABMET->Fill(BetaCM, PZAngLABMET, weight); 
            hist_BetaCM_PZParaLABMET->Fill(BetaCM, PZParaLABMET, weight); 
            hist_BetaCM_PZPerpLABMET->Fill(BetaCM, PZPerpLABMET, weight); 
            hist_RZPara_RZParaLABMET->Fill(RZPara, RZParaLABMET, weight);
            hist_RZPerp_RZPerpLABMET->Fill(RZPerp, RZPerpLABMET, weight);
            hist_PZPara_PZParaLABMET->Fill(PZPara, PZParaLABMET, weight);
            hist_PZPerp_PZPerpLABMET->Fill(PZPerp, PZPerpLABMET, weight);
            hist_PZPara_CandECM->Fill(PZPara, CandECM, weight);
            hist_PZPerp_CandECM->Fill(PZPerp, CandECM, weight);
            hist_MRZPara->Fill(MRZPara, weight);
            hist_MRZPerp->Fill(MRZPerp, weight);
            hist_MRZParaLABMET->Fill(MRZParaLABMET, weight);
            hist_MRZPerpLABMET->Fill(MRZPerpLABMET, weight);
            hist_CandML_MRZPara->Fill(CandML, MRZPara, weight);
            hist_CandML_MRZPerp->Fill(CandML, MRZPerp, weight);
            hist_CandML_MRZParaLABMET->Fill(CandML, MRZParaLABMET, weight);
            hist_CandML_MRZPerpLABMET->Fill(CandML, MRZPerpLABMET, weight);
            hist_CandCosDecayAngle_MPZPara->Fill(CosDecayAngle, MPZPara, weight);
            hist_CandCosDecayAngle_MPZPerp->Fill(CosDecayAngle, MPZPerp, weight);
            hist_CandCosDecayAngle_MPZParaLABMET->Fill(CosDecayAngle, MPZParaLABMET, weight);
            hist_CandCosDecayAngle_MPZPerpLABMET->Fill(CosDecayAngle, MPZPerpLABMET, weight);
            hist_CandCosDecayAngleCM_MPZPara->Fill(CosDecayAngleCM, MPZPara, weight);
            hist_CandCosDecayAngleCM_MPZPerp->Fill(CosDecayAngleCM, MPZPerp, weight);
            hist_CandCosDecayAngleCM_MPZParaLABMET->Fill(CosDecayAngleCM, MPZParaLABMET, weight);
            hist_CandCosDecayAngleCM_MPZPerpLABMET->Fill(CosDecayAngleCM, MPZPerpLABMET, weight);
            hist_Beta_MPZPara->Fill(Beta, MPZPara, weight);
            hist_Beta_MPZPerp->Fill(Beta, MPZPerp, weight);
            hist_Beta_MPZParaLABMET->Fill(Beta, MPZParaLABMET, weight);
            hist_Beta_MPZPerpLABMET->Fill(Beta, MPZPerpLABMET, weight);
            hist_BetaCM_MPZPara->Fill(BetaCM, MPZPara, weight);
            hist_BetaCM_MPZPerp->Fill(BetaCM, MPZPerp, weight);
            hist_BetaCM_MPZParaLABMET->Fill(BetaCM, MPZParaLABMET, weight);
            hist_BetaCM_MPZPerpLABMET->Fill(BetaCM, MPZPerpLABMET, weight);
            hist_PZAng_MPZPara->Fill(PZAng, MPZPara, weight);
            hist_PZAng_MPZPerp->Fill(PZAng, MPZPerp, weight);
            hist_PZAng_MPZParaLABMET->Fill(PZAng, MPZParaLABMET, weight);
            hist_PZAng_MPZPerpLABMET->Fill(PZAng, MPZPerpLABMET, weight);
            hist_CandCosDecayAngle_CandCosDecayAngleCM->Fill(CosDecayAngle, CosDecayAngleCM, weight);
            hist_CandML_BetaZPara->Fill(CandML, BetaZPara, weight);
            hist_CandML_BetaZPerp->Fill(CandML, BetaZPerp, weight);
            hist_CandCosDecayAngle_BetaZPara->Fill(CosDecayAngle, BetaZPara, weight);
            hist_CandCosDecayAngle_BetaZPerp->Fill(CosDecayAngle, BetaZPerp, weight);
            hist_CandCosDecayAngleCM_BetaZPara->Fill(CosDecayAngleCM, BetaZPara, weight);
            hist_CandCosDecayAngleCM_BetaZPerp->Fill(CosDecayAngleCM, BetaZPerp, weight);
            hist_DecayAngleDiff_PzCM->Fill(CosDecayAngleCM - CosDecayAngle, PzCM, weight);
            hist_BetaDiff_PzCM->Fill(BetaCM - Beta, PzCM, weight);
            hist_DecayAngleDiff_BetaDiff->Fill(CosDecayAngleCM - CosDecayAngle, BetaCM - Beta, weight);
            if(V_lep_cands[i].Match() == kMatched){
              hist_MatchedCandML->Fill(CandML, weight);
              hist_MatchedCandBeta->Fill(Beta, weight);
              hist_MatchedCandBetaCM->Fill(BetaCM, weight);
              hist_MatchedCandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandCosDecayAngle->Fill(CosDecayAngle, weight);
              hist_MatchedCandCosDecayAngleCM->Fill(CosDecayAngleCM, weight);
              hist_RISR_MatchedCandML->Fill(RISR, CandML, weight);
              hist_RISR_MatchedCandBeta->Fill(RISR, Beta, weight);
              hist_MatchedCandML_MatchedCandBeta->Fill(CandML, Beta, weight);
              hist_RISR_MatchedCandBetaCM->Fill(RISR, BetaCM, weight);
              hist_MatchedCandML_MatchedCandBetaCM->Fill(CandML, BetaCM, weight);
              hist_RISR_MatchedCandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandML_MatchedCandDeltaPhiMET->Fill(CandML, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_MatchedCandML_MatchedCandCosDecayAngle->Fill(CandML, CosDecayAngle, weight);
              hist_MatchedCandML_MatchedCandCosDecayAngleCM->Fill(CandML, CosDecayAngleCM, weight);
            }
            else{
              hist_UnmatchedCandML->Fill(CandML, weight);
              hist_UnmatchedCandBeta->Fill(Beta, weight);
              hist_UnmatchedCandBetaCM->Fill(BetaCM, weight);
              hist_UnmatchedCandDeltaPhiMET->Fill(V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandCosDecayAngle->Fill(CosDecayAngle, weight);
              hist_UnmatchedCandCosDecayAngleCM->Fill(CosDecayAngleCM, weight);
              hist_RISR_UnmatchedCandML->Fill(RISR, CandML, weight);
              hist_RISR_UnmatchedCandBeta->Fill(RISR, Beta, weight);
              hist_UnmatchedCandML_UnmatchedCandBeta->Fill(CandML, Beta, weight);
              hist_RISR_UnmatchedCandBetaCM->Fill(RISR, BetaCM, weight);
              hist_UnmatchedCandML_UnmatchedCandBetaCM->Fill(CandML, BetaCM, weight);
              hist_RISR_UnmatchedCandDeltaPhiMET->Fill(RISR, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandML_UnmatchedCandDeltaPhiMET->Fill(CandML, V_lep_cands[i].DeltaPhi(TV3_MET), weight);
              hist_UnmatchedCandML_UnmatchedCandCosDecayAngle->Fill(CandML, CosDecayAngle, weight);
              hist_UnmatchedCandML_UnmatchedCandCosDecayAngleCM->Fill(CandML, CosDecayAngleCM, weight);
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
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            if(N_V_lep_cands > 0) EC_X = 5;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
          }
          else if(cat_Nleps == 3) {
            EC_X = 6;
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
            EC_X = 18;
            EC_X += abs_charge;
            hist_EventCount->SetBinContent(EC_X,EC_Y,hist_EventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_Zbi->SetBinContent(EC_X,Zbi_EC_Y,hist_Zbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_EventCount->SetBinContent(EC_X,EC_bins,hist_EventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_Zbi->SetBinContent(EC_X,0,hist_Zbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
            EC_X = 9;
            EC_X += N_V_lep_cands;
            hist_CandsEventCount->SetBinContent(EC_X,EC_Y,hist_CandsEventCount->GetBinContent(EC_X,EC_Y)+weight);
            if(is_signal) hist_CandsZbi->SetBinContent(EC_X,Zbi_EC_Y,hist_CandsZbi->GetBinContent(EC_X,Zbi_EC_Y)+weight);
            if(is_bkg) hist_CandsEventCount->SetBinContent(EC_X,EC_bins,hist_CandsEventCount->GetBinContent(EC_X,EC_bins)+weight);
            if(is_bkg) hist_CandsZbi->SetBinContent(EC_X,0,hist_CandsZbi->GetBinContent(EC_X,0)+weight); // store tot bkg in underflow
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
