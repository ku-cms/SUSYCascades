#include "PlottingTools.h"

void Plot_Advanced(){

  Long64_t start = gSystem->Now();
  RestFrames::SetStyle();
  InitRJRtree();
  //DrawRJRtree(false); gApplication->Terminate(0); // end macro if making tree plot

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD_NTUPLES_v3/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Advanced_NTUPLES_v3_";

  const double CandMinMass = 1.; // 3.5 to get above J/Psi
  const double CandMaxMass = 250.; // 45, 120, 250
  //output_root_file += "CandMinMass"+std::to_string(int(CandMinMass))+"GeV_";
  //output_root_file += "CandMaxMass"+std::to_string(int(CandMaxMass))+"GeV_";

  bool CandSameHemi = false;
  if(CandSameHemi) output_root_file += "SameHemi_";

  //g_Label = "TESTING";
  //g_Label = "Jets RISR0.5";
  //g_Label = "Plotting";
  //g_Label = "NoMET BL";
  g_Label = "3L !Bronze maxSIP3D<3.5 MS>300 4<MVa<65 LepPtTrig Cos2";
  //g_Label = "=3L !Bronze MET150 PTISR150 RISR0.5";
  //g_Label = "=3L !Bronze !0J S";
  //g_Label = "2L PreSelection 0JS";
  //g_Label = "PreSelection 3ObjS HighMass CandCos<0.8 PZPara+ BetaZ<0.9 Exclusive";
  //g_Label = "No_Cuts";
  //g_Label = "ATLAS Cuts";

  bool OSSF = false;
  bool OSOF = false;
  bool SSSF = false;
  bool SSOF = false;
  if(OSSF) g_Label += " OSSF";
  if(OSOF) g_Label += " OSOF";
  if(SSSF) g_Label += " SSSF";
  if(SSOF) g_Label += " SSOF";

  output_root_file += g_Label;
  SanitizeString(output_root_file);
  output_root_file = "plot_outputs/"+output_root_file;
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
  //map_vsignals.insert(std::make_pair("Cascades_220_SMS", VS({"Cascades_SMS_*_220_*"})));
  //map_vsignals.insert(std::make_pair("Cascades_260_SMS", VS({"Cascades_SMS_*_260_*"})));
  //map_vsignals.insert(std::make_pair("Cascades_270_SMS", VS({"Cascades_SMS_*_270_*"})));
  //map_vsignals.insert(std::make_pair("T1bbbb_900_SMS", VS({"T1bbbb_1200_900_*"})));
  //map_vsignals.insert(std::make_pair("T1bbbb_1000_SMS", VS({"T1bbbb_1200_1000_*"})));
  //map_vsignals.insert(std::make_pair("T1bbbb_1175_SMS", VS({"T1bbbb_1200_1176_*"})));
  // loop over signals and add to map
  for(auto p = map_vsignals.begin(); p != map_vsignals.end(); p++){
    ProcessList signals;
    for(auto s : p->second){
      signals += ST.Get(s);
    }
    vec_samples.push_back(std::make_pair(p->first, signals));
  }

  ProcessList signals = ST.Get(kSig);
  signals = signals.Remove("T1").Remove("SMS");
  for(int s = 0; s < int(signals.GetN()); s++){
    ProcessList sig;
    sig += signals[s];
    vec_samples.push_back(std::make_pair(signals[s].Name(), sig));
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
  vector<TH1*> hist_stack_CandMdivPTISR;
  hist_stacks.push_back(&hist_stack_CandMdivPTISR);
  vector<TH1*> hist_stack_MInv;
  hist_stacks.push_back(&hist_stack_MInv);
  vector<TH1*> hist_stack_OneMinusRZPara;
  hist_stacks.push_back(&hist_stack_OneMinusRZPara);
  vector<TH1*> hist_stack_PzCM;
  hist_stacks.push_back(&hist_stack_PzCM);
  vector<TH1*> hist_stack_BetaZLAB;
  hist_stacks.push_back(&hist_stack_BetaZLAB);
  vector<TH1*> hist_stack_RZParaM;
  hist_stacks.push_back(&hist_stack_RZParaM);
  vector<TH1*> hist_stack_RZParaDivRISR;
  hist_stacks.push_back(&hist_stack_RZParaDivRISR);
  vector<TH1*> hist_stack_MCosDecayAngleM;
  hist_stacks.push_back(&hist_stack_MCosDecayAngleM);
  vector<TH1*> hist_stack_MCosDecayAngleP;
  hist_stacks.push_back(&hist_stack_MCosDecayAngleP);
  vector<TH1*> hist_stack_AbsMCosDecayAngleM;
  hist_stacks.push_back(&hist_stack_AbsMCosDecayAngleM);
  vector<TH1*> hist_stack_AbsMCosDecayAngleP;
  hist_stacks.push_back(&hist_stack_AbsMCosDecayAngleP);
  vector<TH1*> hist_stack_SCosDecayAngleP;
  hist_stacks.push_back(&hist_stack_SCosDecayAngleP);

  // hists for holding number of events
  const int EC_bins = vec_samples.size() + 1;
  const int Zbi_bins = map_vsignals.size() + signals.GetN();
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
  const int CF_bins = 25;
  vector<string> vect_str_cutflow_labels(CF_bins, "");

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
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET [GeV]").c_str(), g_NX, 0., 1000.);
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot
    TH1D* hist_RISR = new TH1D((title+"_RISR").c_str(), (title+"_RISR;R_{ISR}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_RISR);
    hist_stack_RISR.push_back(hist_RISR);    
    TH1D* hist_PTISR = new TH1D((title+"_PTISR").c_str(), (title+"_PTISR;p_{T}^{ISR} [GeV]").c_str(), g_NX, 0., 1000.);
    hists1.push_back(hist_PTISR);
    hist_stack_PTISR.push_back(hist_PTISR);
    TH1D* hist_gammaPerp = new TH1D((title+"_gammaPerp").c_str(), (title+"_gammaPerp;#gamma_{#perp}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_gammaPerp);
    hist_stack_gammaPerp.push_back(hist_gammaPerp);
    TH1D* hist_MQperp = new TH1D((title+"_MQperp").c_str(), (title+"_MQperp;M_{#perp}").c_str(), g_NX, 0., 175.);
    hists1.push_back(hist_MQperp);
    hist_stack_MQperp.push_back(hist_MQperp);
    TH1D* hist_MSperpCM0 = new TH1D((title+"_MSperpCM0").c_str(), (title+"_MSperpCM0;MS_{#perp CM0}").c_str(), g_NX, 0., 600.);
    hists1.push_back(hist_MSperpCM0);
    hist_stack_MSperpCM0.push_back(hist_MSperpCM0);
    TH1D* hist_MQperpCM0 = new TH1D((title+"_MQperpCM0").c_str(), (title+"_MQperpCM0;MQ_{#perp CM0}").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_MQperpCM0);
    hist_stack_MQperpCM0.push_back(hist_MQperpCM0);
    TH1D* hist_gammaPerpCM0 = new TH1D((title+"_gammaPerpCM0").c_str(), (title+"_gammaPerpCM0;#gamma_{#perp CM0}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_gammaPerpCM0);
    hist_stack_gammaPerpCM0.push_back(hist_gammaPerpCM0);
    TH1D* hist_MSCM0 = new TH1D((title+"_MSCM0").c_str(), (title+"_MSCM0;MS_{CM0}").c_str(), g_NX, 0., 600.);
    hists1.push_back(hist_MSCM0);
    hist_stack_MSCM0.push_back(hist_MSCM0);
    TH1D* hist_MQCM0 = new TH1D((title+"_MQCM0").c_str(), (title+"_MQCM0;MQ_{CM0}").c_str(), g_NX, 0., 200.);
    hists1.push_back(hist_MQCM0);
    hist_stack_MQCM0.push_back(hist_MQCM0);
    TH1D* hist_gammaCM0 = new TH1D((title+"_gammaCM0").c_str(), (title+"_gammaCM0;#gamma_{CM0}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_gammaCM0);
    hist_stack_gammaCM0.push_back(hist_gammaCM0);
    TH1D* hist_ML = new TH1D((title+"_ML").c_str(), (title+"_ML;ML").c_str(), g_NX, 0., 175.);
    hists1.push_back(hist_ML);
    hist_stack_ML.push_back(hist_ML);
    TH1D* hist_MJ = new TH1D((title+"_MJ").c_str(), (title+"_MJ;MJ").c_str(), g_NX, 0., 175.);
    hists1.push_back(hist_MJ);
    hist_stack_MJ.push_back(hist_MJ);

    TH1D* hist_CandML = new TH1D((title+"_CandML").c_str(), (title+"_CandML;CandML").c_str(), g_NX, CandMinMass, CandMaxMass);
    hists1.push_back(hist_CandML);
    hist_stack_CandML.push_back(hist_CandML);
    TH1D* hist_CandBeta = new TH1D((title+"_CandBeta").c_str(), (title+"_CandBeta;#beta").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_CandBeta);
    hist_stack_CandBeta.push_back(hist_CandBeta);
    TH1D* hist_CandBetaCM = new TH1D((title+"_CandBetaCM").c_str(), (title+"_CandBetaCM;#beta^{CM}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_CandBetaCM);
    hist_stack_CandBetaCM.push_back(hist_CandBetaCM);
    TH1D* hist_CandDeltaPhiMET = new TH1D((title+"_CandDeltaPhiMET").c_str(), (title+"_CandDeltaPhiMET;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0., 3.15);
    hists1.push_back(hist_CandDeltaPhiMET);
    hist_stack_CandDeltaPhiMET.push_back(hist_CandDeltaPhiMET);
    TH1D* hist_CandCosDecayAngle = new TH1D((title+"_CandCosDecayAngle").c_str(), (title+"_CandCosDecayAngle;Z* candidate cos#theta").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_CandCosDecayAngle);
    hist_stack_CandCosDecayAngle.push_back(hist_CandCosDecayAngle);
    TH1D* hist_CandCosDecayAngleCM = new TH1D((title+"_CandCosDecayAngleCM").c_str(), (title+"_CandCosDecayAngleCM;Z* candidate cos#theta CM").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_CandCosDecayAngleCM);
    hist_stack_CandCosDecayAngleCM.push_back(hist_CandCosDecayAngleCM);
    TH1D* hist_RZPara = new TH1D((title+"_RZPara").c_str(), (title+"_RZPara;RZPara").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_RZPara);
    hist_stack_RZPara.push_back(hist_RZPara);
    TH1D* hist_RZPerp = new TH1D((title+"_RZPerp").c_str(), (title+"_RZPerp;RZPerp").c_str(), g_NX, 0., 1.5);
    hists1.push_back(hist_RZPerp);
    hist_stack_RZPerp.push_back(hist_RZPerp);
    TH1D* hist_PZAng = new TH1D((title+"_PZAng").c_str(), (title+"_PZAng;PZAng").c_str(), g_NX, -1.6, 1.6);
    hists1.push_back(hist_PZAng);
    hist_stack_PZAng.push_back(hist_PZAng);
    TH1D* hist_PZPara = new TH1D((title+"_PZPara").c_str(), (title+"_PZPara;PZPara").c_str(), g_NX, 0., 600.);
    hists1.push_back(hist_PZPara);
    hist_stack_PZPara.push_back(hist_PZPara);
    TH1D* hist_PZPerp = new TH1D((title+"_PZPerp").c_str(), (title+"_PZPerp;PZPerp").c_str(), g_NX, 0., 600.);
    hists1.push_back(hist_PZPerp);
    hist_stack_PZPerp.push_back(hist_PZPerp);
    TH1D* hist_RZParaLABMET = new TH1D((title+"_RZParaLABMET").c_str(), (title+"_RZParaLABMET;RZParaLABMET").c_str(), g_NX, -1.5, 1.5);
    hists1.push_back(hist_RZParaLABMET);
    hist_stack_RZParaLABMET.push_back(hist_RZParaLABMET);
    TH1D* hist_RZPerpLABMET = new TH1D((title+"_RZPerpLABMET").c_str(), (title+"_RZPerpLABMET;RZPerpLABMET").c_str(), g_NX, 0., 1.5);
    hists1.push_back(hist_RZPerpLABMET);
    hist_stack_RZPerpLABMET.push_back(hist_RZPerpLABMET);
    TH1D* hist_PZAngLABMET = new TH1D((title+"_PZAngLABMET").c_str(), (title+"_PZAngLABMET;PZAngLABMET").c_str(), g_NX, -1.6, 1.6);
    hists1.push_back(hist_PZAngLABMET);
    hist_stack_PZAngLABMET.push_back(hist_PZAngLABMET);

    TH1D* hist_MRZPara = new TH1D((title+"_MRZPara").c_str(), (title+"_MRZPara;M/RZPara").c_str(), g_NX, 0., 1000.);
    hists1.push_back(hist_MRZPara);
    hist_stack_MRZPara.push_back(hist_MRZPara);
    TH1D* hist_MRZPerp = new TH1D((title+"_MRZPerp").c_str(), (title+"_MRZPerp;M/RZPerp").c_str(), g_NX, 0., 500.);
    hists1.push_back(hist_MRZPerp);
    hist_stack_MRZPerp.push_back(hist_MRZPerp);
    TH1D* hist_MRZParaLABMET = new TH1D((title+"_MRZParaLABMET").c_str(), (title+"_MRZParaLABMET;M/RZParaLABMET").c_str(), g_NX, 0., 500.);
    hists1.push_back(hist_MRZParaLABMET);
    hist_stack_MRZParaLABMET.push_back(hist_MRZParaLABMET);
    TH1D* hist_MRZPerpLABMET = new TH1D((title+"_MRZPerpLABMET").c_str(), (title+"_MRZPerpLABMET;M/RZPerpLABMET").c_str(), g_NX, 0., 500.);
    hists1.push_back(hist_MRZPerpLABMET);
    hist_stack_MRZPerpLABMET.push_back(hist_MRZPerpLABMET);

    TH1D* hist_CandMdivPTISR = new TH1D((title+"_CandMdivPTISR").c_str(), (title+"_CandMdivPTISR;CandM/PTISR").c_str(), g_NX, 0., 0.8);
    hists1.push_back(hist_CandMdivPTISR);
    hist_stack_CandMdivPTISR.push_back(hist_CandMdivPTISR);
    TH1D* hist_MInv = new TH1D((title+"_MInv").c_str(), (title+"_MInv;MInv").c_str(), g_NX, 0., 1000.);
    hists1.push_back(hist_MInv);
    hist_stack_MInv.push_back(hist_MInv);
    TH1D* hist_OneMinusRZPara = new TH1D((title+"_OneMinusRZPara").c_str(), (title+"_OneMinusRZPara;1-RZPara").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_OneMinusRZPara);
    hist_stack_OneMinusRZPara.push_back(hist_OneMinusRZPara);
    TH1D* hist_PzCM = new TH1D((title+"_PzCM").c_str(), (title+"_PzCM;PzCM").c_str(), g_NX, 0., 700.);
    hists1.push_back(hist_PzCM);
    hist_stack_PzCM.push_back(hist_PzCM);
    TH1D* hist_BetaZLAB = new TH1D((title+"_BetaZLAB").c_str(), (title+"_BetaZLAB;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_BetaZLAB);
    hist_stack_BetaZLAB.push_back(hist_BetaZLAB);
    TH1D* hist_RZParaM = new TH1D((title+"_RZParaM").c_str(), (title+"_RZParaM;RZPara/M").c_str(), g_NX, 0., 0.02);
    hists1.push_back(hist_RZParaM);
    hist_stack_RZParaM.push_back(hist_RZParaM);
    TH1D* hist_RZParaDivRISR = new TH1D((title+"_RZParaDivRISR").c_str(), (title+"_RZParaDivRISR;RZPara/RISR").c_str(), g_NX, 0., 2.);
    hists1.push_back(hist_RZParaDivRISR);
    hist_stack_RZParaDivRISR.push_back(hist_RZParaDivRISR);
    TH1D* hist_DeltaRLeptonMatching = new TH1D((title+"_DeltaRLeptonMatching").c_str(), (title+"_DeltaRLeptonMatching;#Delta R Matched Leptons").c_str(), g_NX, 0., 0.01);
    hists1.push_back(hist_DeltaRLeptonMatching);
    TH1D* hist_MCosDecayAngleM = new TH1D((title+"_MCosDecayAngleM").c_str(), (title+"_MCosDecayAngleM;cos#theta_{Pa}^{LEP}*cos#theta_{Pb}^{LEP}").c_str(), g_NX, -1., 1.);
    hists1.push_back(hist_MCosDecayAngleM);
    hist_stack_MCosDecayAngleM.push_back(hist_MCosDecayAngleM);
    TH1D* hist_MCosDecayAngleP = new TH1D((title+"_MCosDecayAngleP").c_str(), (title+"_MCosDecayAngleP;cos#theta_{Pa}^{LEP}-cos#theta_{Pb}^{LEP}").c_str(), g_NX, -2., 2.);
    hists1.push_back(hist_MCosDecayAngleP);
    hist_stack_MCosDecayAngleP.push_back(hist_MCosDecayAngleP);

    TH1D* hist_AbsMCosDecayAngleM = new TH1D((title+"_AbsMCosDecayAngleM").c_str(), (title+"_AbsMCosDecayAngleM;|cos#theta_{Pa}^{LEP}*cos#theta_{Pb}^{LEP}|").c_str(), g_NX, 0., 1.);
    hists1.push_back(hist_AbsMCosDecayAngleM);
    hist_stack_AbsMCosDecayAngleM.push_back(hist_AbsMCosDecayAngleM);
    TH1D* hist_AbsMCosDecayAngleP = new TH1D((title+"_AbsMCosDecayAngleP").c_str(), (title+"_AbsMCosDecayAngleP;|cos#theta_{Pa}^{LEP}-cos#theta_{Pb}^{LEP}|").c_str(), g_NX, 0., 2.);
    hists1.push_back(hist_AbsMCosDecayAngleP);
    hist_stack_AbsMCosDecayAngleP.push_back(hist_AbsMCosDecayAngleP);
    TH1D* hist_SCosDecayAngleP = new TH1D((title+"_SCosDecayAngleP").c_str(), (title+"_SCosDecayAngleP;cos#theta_{Pa}^{LEP}+cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 2.);
    hists1.push_back(hist_SCosDecayAngleP);
    hist_stack_SCosDecayAngleP.push_back(hist_SCosDecayAngleP);

    TH2D* hist_RISR_PTISR = new TH2D((title+"_RISR_PTISR").c_str(), (title+"_RISR_PTISR;R_{ISR};p_{T}^{ISR} [GeV]").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    hists2.push_back(hist_RISR_PTISR);
    TH2D* hist_RISR_Mperp = new TH2D((title+"_RISR_Mperp").c_str(), (title+"_RISR_Mperp;R_{ISR};M_{#perp} [GeV]").c_str(), g_NX, 0., 1., g_NX, 0., 100.);
    //hists2.push_back(hist_RISR_Mperp);
    TH2D* hist_dphiCMI_PTCM = new TH2D((title+"_dphiCMI_PTCM").c_str(), (title+"_dphiCMI_PTCM;#Delta #phi_{(CM,I)};p_{T}^{CM}").c_str(), g_NX, 0., 3.15, g_NX, 0., 500.);
    //hists2.push_back(hist_dphiCMI_PTCM);
    TH2D* hist_dphiMETV_PTISR = new TH2D((title+"_dphiMETV_PTISR").c_str(), (title+"_dphiMETV_PTISR;#Delta #phi_{(I,V)};p_{T}^{ISR}").c_str(), g_NX, 0., 3.15, g_NX, 200., 800.);
    //hists2.push_back(hist_dphiMETV_PTISR);
    TH2D* hist_gammaPerp_RISR = new TH2D((title+"_gammaPerp_RISR").c_str(), (title+"_gammaPerp_RISR;#gamma_{#perp};R_{ISR}").c_str(), g_NX, 0., 1., g_NX, 0.5, 1.);
    //hists2.push_back(hist_gammaPerp_RISR);
    TH2D* hist_RISR_mllLEAD = new TH2D((title+"_RISR_mllLEAD").c_str(), (title+"_RISR_mllLEAD;R_{ISR};m_{ll}LEAD").c_str(), g_NX, 0.5, 1., g_NX, 0., 175.);
    //hists2.push_back(hist_RISR_mllLEAD);
    TH2D* hist_RISR_mL = new TH2D((title+"_RISR_mL").c_str(), (title+"_RISR_mL;R_{ISR};mL").c_str(), g_NX, 0.5, 1., g_NX, 0., 175.); // mass of leptonic system
    //hists2.push_back(hist_RISR_mL);

    TH2D* hist_MSperpCM0_RISR = new TH2D((title+"_MSperpCM0_RISR").c_str(), (title+"_MSperpCM0_RISR;MS_{#perp CM0};R_{ISR}").c_str(), g_NX, 0., 500., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MSperpCM0_RISR);
    TH2D* hist_MQperpCM0_RISR = new TH2D((title+"_MQperpCM0_RISR").c_str(), (title+"_MQperpCM0_RISR;MQ_{#perp CM0};R_{ISR}").c_str(), g_NX, 0., 300., g_NX, 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
    //hists2.push_back(hist_MQperpCM0_RISR);
    TH2D* hist_gammaPerpCM0_RISR = new TH2D((title+"_gammaPerpCM0_RISR").c_str(), (title+"_gammaPerpCM0_RISR;#gamma_{#perp CM0};R_{ISR}").c_str(), g_NX, 0., 1., g_NX, 0.5, 1.);
    //hists2.push_back(hist_gammaPerpCM0_RISR);

    TH2D* hist_MSperpCM0_gammaPerpCM0 = new TH2D((title+"_MSperpCM0_gammaPerpCM0").c_str(), (title+"_MSperpCM0_gammaPerpCM0;MS_{#perp CM0};#gamma_{#perp CM0}").c_str(), g_NX, 0., 500., g_NX, 0., 1.);
    //hists2.push_back(hist_MSperpCM0_gammaPerpCM0);
    TH2D* hist_MSperpCM0_MQperpCM0 = new TH2D((title+"_MSperpCM0_MQperpCM0").c_str(), (title+"_MSperpCM0_MQperpCM0;MS_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX, 0., 500., g_NX, 0., 300.);
    //hists2.push_back(hist_MSperpCM0_MQperpCM0);
    TH2D* hist_gammaPerpCM0_MQperpCM0 = new TH2D((title+"_gammaPerpCM0_MQperpCM0").c_str(), (title+"_gammaPerpCM0_MQperpCM0;#gamma_{#perp CM0};MQ_{#perp CM0}").c_str(), g_NX, 0., 1., g_NX, 0., 300.);
    //hists2.push_back(hist_gammaPerpCM0_MQperpCM0);
    
    TH2D* hist_MJ_MQperpCM0 = new TH2D((title+"_MJ_MQperpCM0").c_str(), (title+"_MJ_MQperpCM0;MJ;MQ_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MJ_MQperpCM0);
    TH2D* hist_MLa_MQperpCM0 = new TH2D((title+"_MLa_MQperpCM0").c_str(), (title+"_MLa_MQperpCM0;MLa;MQ_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLa_MQperpCM0);
    TH2D* hist_MLb_MQperpCM0 = new TH2D((title+"_MLb_MQperpCM0").c_str(), (title+"_MLb_MQperpCM0;MLb;MQ_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLb_MQperpCM0);
    TH2D* hist_MJ_MSperpCM0 = new TH2D((title+"_MJ_MSperpCM0").c_str(), (title+"_MJ_MSperpCM0;MJ;MS_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MJ_MSperpCM0);
    TH2D* hist_MLa_MSperpCM0 = new TH2D((title+"_MLa_MSperpCM0").c_str(), (title+"_MLa_MSperpCM0;MLa;MS_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MLa_MSperpCM0);
    TH2D* hist_MLb_MSperpCM0 = new TH2D((title+"_MLb_MSperpCM0").c_str(), (title+"_MLb_MSperpCM0;MLb;MS_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MLb_MSperpCM0);
    TH2D* hist_MJ_gammaPerpCM0 = new TH2D((title+"_MJ_gammaPerpCM0").c_str(), (title+"_MJ_gammaPerpCM0;MJ;#gamma_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MJ_gammaPerpCM0);
    TH2D* hist_MLa_gammaPerpCM0 = new TH2D((title+"_MLa_gammaPerpCM0").c_str(), (title+"_MLa_gammaPerpCM0;MLa;#gamma_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MLa_gammaPerpCM0);
    TH2D* hist_MLb_gammaPerpCM0 = new TH2D((title+"_MLb_gammaPerpCM0").c_str(), (title+"_MLb_gammaPerpCM0;MLb;#gamma_{#perp CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MLb_gammaPerpCM0);

    TH2D* hist_MSCM0_RISR = new TH2D((title+"_MSCM0_RISR").c_str(), (title+"_MSCM0_RISR;MS_{CM0};RISR").c_str(), g_NX, 0., 500., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MSCM0_RISR);
    TH2D* hist_MQCM0_RISR = new TH2D((title+"_MQCM0_RISR").c_str(), (title+"_MQCM0_RISR;MQ_{CM0};RISR").c_str(), g_NX, 0., 300., g_NX, 0.5, 1.); // quad sum of Ma and Mb / sqrt(2)
    //hists2.push_back(hist_MQCM0_RISR);
    TH2D* hist_gammaCM0_RISR = new TH2D((title+"_gammaCM0_RISR").c_str(), (title+"_gammaCM0_RISR;#gamma_{CM0};RISR").c_str(), g_NX, 0., 1., g_NX, 0.5, 1.);
    //hists2.push_back(hist_gammaCM0_RISR);

    TH2D* hist_MSCM0_gammaCM0 = new TH2D((title+"_MSCM0_gammaCM0").c_str(), (title+"_MSCM0_gammaCM0;MS_{CM0};#gamma_{CM0}").c_str(), g_NX, 0., 500., g_NX, 0., 1.);
    //hists2.push_back(hist_MSCM0_gammaCM0);
    TH2D* hist_MSCM0_MQCM0 = new TH2D((title+"_MSCM0_MQCM0").c_str(), (title+"_MSCM0_MQCM0;MS_{CM0};MQ_{CM0}").c_str(), g_NX, 0., 500., g_NX, 0., 300.);
    //hists2.push_back(hist_MSCM0_MQCM0);
    TH2D* hist_gammaCM0_MQCM0 = new TH2D((title+"_gammaCM0_MQCM0").c_str(), (title+"_gammaCM0_MQCM0;#gamma_{CM0};MQ_{CM0}").c_str(), g_NX, 0., 1., g_NX, 0., 300.);
    //hists2.push_back(hist_gammaCM0_MQCM0);
    
    TH2D* hist_MJ_MQCM0 = new TH2D((title+"_MJ_MQCM0").c_str(), (title+"_MJ_MQCM0;MJ;MQ_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MJ_MQCM0);
    TH2D* hist_MLa_MQCM0 = new TH2D((title+"_MLa_MQCM0").c_str(), (title+"_MLa_MQCM0;MLa;MQ_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLa_MQCM0);
    TH2D* hist_MLb_MQCM0 = new TH2D((title+"_MLb_MQCM0").c_str(), (title+"_MLb_MQCM0;MLb;MQ_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLb_MQCM0);
    TH2D* hist_MJ_MSCM0 = new TH2D((title+"_MJ_MSCM0").c_str(), (title+"_MJ_MSCM0;MJ;MS_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MJ_MSCM0);
    TH2D* hist_MLa_MSCM0 = new TH2D((title+"_MLa_MSCM0").c_str(), (title+"_MLa_MSCM0;MLa;MS_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MLa_MSCM0);
    TH2D* hist_MLb_MSCM0 = new TH2D((title+"_MLb_MSCM0").c_str(), (title+"_MLb_MSCM0;MLb;MS_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 500.);
    //hists2.push_back(hist_MLb_MSCM0);
    TH2D* hist_MJ_gammaCM0 = new TH2D((title+"_MJ_gammaCM0").c_str(), (title+"_MJ_gammaCM0;MJ;#gamma_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MJ_gammaCM0);
    TH2D* hist_MLa_gammaCM0 = new TH2D((title+"_MLa_gammaCM0").c_str(), (title+"_MLa_gammaCM0;MLa;#gamma_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MLa_gammaCM0);
    TH2D* hist_MLb_gammaCM0 = new TH2D((title+"_MLb_gammaCM0").c_str(), (title+"_MLb_gammaCM0;MLb;#gamma_{CM0}").c_str(), g_NX, 0., 300., g_NX, 0., 1.);
    //hists2.push_back(hist_MLb_gammaCM0);

    TH2D* hist_MLb_MJ = new TH2D((title+"_MLb_MJ").c_str(), (title+"_MLb_MJ;MLb;MJ").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLb_MJ);
    TH2D* hist_MLa_MJ = new TH2D((title+"_MLa_MJ").c_str(), (title+"_MLa_MJ;MLa;MJ").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLa_MJ);
    TH2D* hist_MLa_MLb = new TH2D((title+"_MLa_MLb").c_str(), (title+"_MLa_MLb;MLa;MLb").c_str(), g_NX, 0., 300., g_NX, 0., 300.);
    //hists2.push_back(hist_MLa_MLb);

    TH2D* hist_RISR_MJ = new TH2D((title+"_RISR_MJ").c_str(), (title+"_RISR_MJ;R_{ISR};MJ").c_str(), g_NX, 0.5, 1., g_NX, 0., 300.);
    //hists2.push_back(hist_RISR_MJ);
    TH2D* hist_RISR_MLa = new TH2D((title+"_RISR_MLa").c_str(), (title+"_RISR_MLa;R_{ISR};MLa").c_str(), g_NX, 0.5, 1., g_NX, 0., 300.);
    //hists2.push_back(hist_RISR_MLa);
    TH2D* hist_RISR_MLb = new TH2D((title+"_RISR_MLb").c_str(), (title+"_RISR_MLb;R_{ISR};MLb").c_str(), g_NX, 0.5, 1., g_NX, 0., 300.);
    //hists2.push_back(hist_RISR_MLb);

    TH2D* hist_RISR_CandML = new TH2D((title+"_RISR_CandML").c_str(), (title+"_RISR_CandML;R_{ISR};Z* candidate M [GeV]").c_str(), g_NX, 0.5, 1., g_NX, CandMinMass, CandMaxMass);
    //hists2.push_back(hist_RISR_CandML);
    TH2D* hist_RISR_CandBeta = new TH2D((title+"_RISR_CandBeta").c_str(), (title+"_RISR_CandBeta;R_{ISR};Z* candidate #beta").c_str(), g_NX, 0.5, 1., g_NX, 0., 1.);
    //hists2.push_back(hist_RISR_CandBeta);
    TH2D* hist_RISR_CandBetaCM = new TH2D((title+"_RISR_CandBetaCM").c_str(), (title+"_RISR_CandBetaCM;R_{ISR};Z* candidate #beta^{CM}").c_str(), g_NX, 0.5, 1., g_NX, 0., 1.);
    //hists2.push_back(hist_RISR_CandBetaCM);
    TH2D* hist_RISR_CandDeltaPhiMET = new TH2D((title+"_RISR_CandDeltaPhiMET").c_str(), (title+"_RISR_CandDeltaPhiMET;R_{ISR};#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0.5, 1., g_NX, 0., 3.15);
    //hists2.push_back(hist_RISR_CandDeltaPhiMET);
    TH2D* hist_CandML_CandBeta = new TH2D((title+"_CandML_CandBeta").c_str(), (title+"_CandML_CandBeta;Z* candidate M [GeV];Z* candidate #beta").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_CandBeta);
    TH2D* hist_CandML_CandBetaCM = new TH2D((title+"_CandML_CandBetaCM").c_str(), (title+"_CandML_CandBetaCM;Z* candidate M [GeV];Z* candidate #beta^{CM}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_CandBetaCM);
    TH2D* hist_CandML_CandDeltaPhiMET = new TH2D((title+"_CandML_CandDeltaPhiMET").c_str(), (title+"_CandML_CandDeltaPhiMET;Z* candidate M [GeV];#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 3.15);
    //hists2.push_back(hist_CandML_CandDeltaPhiMET);
    TH2D* hist_CandML_CandCosDecayAngle = new TH2D((title+"_CandML_CandCosDecayAngle").c_str(), (title+"_CandML_CandCosDecayAngle;Z* candidate M [GeV];Z* candidate cos#theta").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_CandCosDecayAngle);
    TH2D* hist_Beta_CandCosDecayAngle = new TH2D((title+"_Beta_CandCosDecayAngle").c_str(), (title+"_Beta_CandCosDecayAngle;Z* candidate #beta;Z* candidate cos#theta").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_Beta_CandCosDecayAngle);
    TH2D* hist_BetaCM_CandCosDecayAngle = new TH2D((title+"_BetaCM_CandCosDecayAngle").c_str(), (title+"_BetaCM_CandCosDecayAngle;Z* candidate #beta^{CM};Z* candidate cos#theta").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_BetaCM_CandCosDecayAngle);
    TH2D* hist_CandCosDecayAngle_CandDeltaPhiMET = new TH2D((title+"_CandCosDecayAngle_CandDeltaPhiMET").c_str(), (title+"_CandCosDecayAngle_CandDeltaPhiMET;Z* candidate cos#theta;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0., 1., g_NX, 0., 3.15);
    //hists2.push_back(hist_CandCosDecayAngle_CandDeltaPhiMET);
    TH2D* hist_CandML_CandCosDecayAngleCM = new TH2D((title+"_CandML_CandCosDecayAngleCM").c_str(), (title+"_CandML_CandCosDecayAngleCM;Z* candidate M [GeV];Z* candidate cos#theta CM").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_CandCosDecayAngleCM);
    TH2D* hist_Beta_CandCosDecayAngleCM = new TH2D((title+"_Beta_CandCosDecayAngleCM").c_str(), (title+"_Beta_CandCosDecayAngleCM;Z* candidate #beta;Z* candidate cos#theta CM").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_Beta_CandCosDecayAngleCM);
    TH2D* hist_BetaCM_CandCosDecayAngleCM = new TH2D((title+"_BetaCM_CandCosDecayAngleCM").c_str(), (title+"_BetaCM_CandCosDecayAngleCM;Z* candidate #beta^{CM};Z* candidate cos#theta CM").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_BetaCM_CandCosDecayAngleCM);
    TH2D* hist_CandCosDecayAngleCM_CandDeltaPhiMET = new TH2D((title+"_CandCosDecayAngleCM_CandDeltaPhiMET").c_str(), (title+"_CandCosDecayAngleCM_CandDeltaPhiMET;Z* candidate cos#theta CM;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0., 1., g_NX, 0., 3.15);
    //hists2.push_back(hist_CandCosDecayAngleCM_CandDeltaPhiMET);
    TH2D* hist_Beta_CandDeltaPhiMET = new TH2D((title+"_Beta_CandDeltaPhiMET").c_str(), (title+"_Beta_CandDeltaPhiMET;Z* candidate #beta;#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0., 1., g_NX, 0., 3.15);
    //hists2.push_back(hist_Beta_CandDeltaPhiMET);
    TH2D* hist_BetaCM_CandDeltaPhiMET = new TH2D((title+"_BetaCM_CandDeltaPhiMET").c_str(), (title+"_BetaCM_CandDeltaPhiMET;Z* candidate #beta^{CM};#Delta #phi(#vec{#slash{E}_{T}},Z* cand)").c_str(), g_NX, 0., 1., g_NX, 0., 3.15);
    //hists2.push_back(hist_BetaCM_CandDeltaPhiMET);
    TH2D* hist_CandML_RZPara = new TH2D((title+"_CandML_RZPara").c_str(), (title+"_CandML_RZPara;Z* candidate M [GeV];RZPara").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 0.5);
    //hists2.push_back(hist_CandML_RZPara);
    TH2D* hist_CandML_RZPerp = new TH2D((title+"_CandML_RZPerp").c_str(), (title+"_CandML_RZPerp;Z* candidate M [GeV];RZPerp").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 0.5);
    //hists2.push_back(hist_CandML_RZPerp);
    TH2D* hist_RZPara_RZPerp = new TH2D((title+"_RZPara_RZPerp").c_str(), (title+"_RZPara_RZPerp;RZPara;RZPerp").c_str(), g_NX, 0, 0.5, g_NX, 0., 0.5);
    //hists2.push_back(hist_RZPara_RZPerp);
    TH2D* hist_CandML_PZAng = new TH2D((title+"_CandML_PZAng").c_str(), (title+"_CandML_PZAng;Z* candidate M [GeV];PZAng").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandML_PZAng);
    TH2D* hist_CandML_PZPara = new TH2D((title+"_CandML_PZPara").c_str(), (title+"_CandML_PZPara;Z* candidate M [GeV];PZPara").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 500.);
    //hists2.push_back(hist_CandML_PZPara);
    TH2D* hist_CandML_PZPerp = new TH2D((title+"_CandML_PZPerp").c_str(), (title+"_CandML_PZPerp;Z* candidate M [GeV];PZPerp").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 500.);
    //hists2.push_back(hist_CandML_PZPerp);

    TH2D* hist_CandML_RZParaLABMET = new TH2D((title+"_CandML_RZParaLABMET").c_str(), (title+"_CandML_RZParaLABMET;Z* candidate M [GeV];RZParaLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, -0.5, 0.5);
    //hists2.push_back(hist_CandML_RZParaLABMET);
    TH2D* hist_CandML_RZPerpLABMET = new TH2D((title+"_CandML_RZPerpLABMET").c_str(), (title+"_CandML_RZPerpLABMET;Z* candidate M [GeV];RZPerpLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 0.5);
    //hists2.push_back(hist_CandML_RZPerpLABMET);
    TH2D* hist_CandML_PZAngLABMET = new TH2D((title+"_CandML_PZAngLABMET").c_str(), (title+"_CandML_PZAngLABMET;Z* candidate M [GeV];PZAngLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandML_PZAngLABMET);
    TH2D* hist_RZParaLABMET_RZPerpLABMET = new TH2D((title+"_RZParaLABMET_RZPerpLABMET").c_str(), (title+"_RZParaLABMET_RZPerpLABMET;RZParaLABMET;RZPerpLABMET").c_str(), g_NX, -0.5, 0.5, g_NX, 0., 0.5);
    //hists2.push_back(hist_RZParaLABMET_RZPerpLABMET);
    TH2D* hist_CandML_PZParaLABMET = new TH2D((title+"_CandML_PZParaLABMET").c_str(), (title+"_CandML_PZParaLABMET;Z* candidate M [GeV];PZParaLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_CandML_PZParaLABMET);
    TH2D* hist_CandML_PZPerpLABMET = new TH2D((title+"_CandML_PZPerpLABMET").c_str(), (title+"_CandML_PZPerpLABMET;Z* candidate M [GeV];PZPerpLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 300.);
    //hists2.push_back(hist_CandML_PZPerpLABMET);

    TH2D* hist_RZPara_RZParaLABMET = new TH2D((title+"_RZPara_RZParaLABMET").c_str(), (title+"_RZPara_RZParaLABMET;RZPara;RZParaLABMET").c_str(), g_NX, 0., 0.75, g_NX, -0.75, 0.75);
    //hists2.push_back(hist_RZPara_RZParaLABMET);
    TH2D* hist_RZPerp_RZPerpLABMET = new TH2D((title+"_RZPerp_RZPerpLABMET").c_str(), (title+"_RZPerp_RZPerpLABMET;RZPerp;RZPerpLABMET").c_str(), g_NX, 0., 0.75, g_NX, 0., 0.75);
    //hists2.push_back(hist_RZPerp_RZPerpLABMET);
    TH2D* hist_PZPara_PZParaLABMET = new TH2D((title+"_PZPara_PZParaLABMET").c_str(), (title+"_PZPara_PZParaLABMET;PZPara;PZParaLABMET").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_PZPara_PZParaLABMET);
    TH2D* hist_PZPerp_PZPerpLABMET = new TH2D((title+"_PZPerp_PZPerpLABMET").c_str(), (title+"_PZPerp_PZPerpLABMET;PZPerp;PZPerpLABMET").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_PZPerp_PZPerpLABMET);

    TH2D* hist_CandML_BetaZPara = new TH2D((title+"_CandML_BetaZPara").c_str(), (title+"_CandML_BetaZPara;Z* candidate M [GeV];PZPara/E^{CM}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_BetaZPara);
    TH2D* hist_CandML_BetaZPerp = new TH2D((title+"_CandML_BetaZPerp").c_str(), (title+"_CandML_BetaZPerp;Z* candidate M [GeV];PZPerp/E^{CM}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_BetaZPerp);
    TH2D* hist_CandCosDecayAngle_BetaZPara = new TH2D((title+"_CandCosDecayAngle_BetaZPara").c_str(), (title+"_CandCosDecayAngle_BetaZPara;Z* cos#theta;PZPara/E^{CM}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_BetaZPara);
    TH2D* hist_CandCosDecayAngle_BetaZPerp = new TH2D((title+"_CandCosDecayAngle_BetaZPerp").c_str(), (title+"_CandCosDecayAngle_BetaZPerp;Z* cos#theta;PZPerp/E^{CM}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_BetaZPerp);
    TH2D* hist_CandCosDecayAngleCM_BetaZPara = new TH2D((title+"_CandCosDecayAngleCM_BetaZPara").c_str(), (title+"_CandCosDecayAngleCM_BetaZPara;Z* cos#theta CM;PZPara/E^{CM}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngleCM_BetaZPara);
    TH2D* hist_CandCosDecayAngleCM_BetaZPerp = new TH2D((title+"_CandCosDecayAngleCM_BetaZPerp").c_str(), (title+"_CandCosDecayAngleCM_BetaZPerp;Z* cos#theta CM;PZPerp/E^{CM}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngleCM_BetaZPerp);

    TH2D* hist_CandML_MRZPara = new TH2D((title+"_CandML_MRZPara").c_str(), (title+"_CandML_MRZPara;Z* candidate M [GeV];M/RZPara").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1000.);
    //hists2.push_back(hist_CandML_MRZPara);
    TH2D* hist_CandML_MRZPerp = new TH2D((title+"_CandML_MRZPerp").c_str(), (title+"_CandML_MRZPerp;Z* candidate M [GeV];M/RZPerp").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 500.);
    //hists2.push_back(hist_CandML_MRZPerp);
    TH2D* hist_CandML_MRZParaLABMET = new TH2D((title+"_CandML_MRZParaLABMET").c_str(), (title+"_CandML_MRZParaLABMET;Z* candidate M [GeV];M/RZParaLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 500.);
    //hists2.push_back(hist_CandML_MRZParaLABMET);
    TH2D* hist_CandML_MRZPerpLABMET = new TH2D((title+"_CandML_MRZPerpLABMET").c_str(), (title+"_CandML_MRZPerpLABMET;Z* candidate M [GeV];M/RZPerpLABMET").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 500.);
    //hists2.push_back(hist_CandML_MRZPerpLABMET);

    TH2D* hist_CandCosDecayAngle_MPZPara = new TH2D((title+"_CandCosDecayAngle_MPZPara").c_str(), (title+"_CandCosDecayAngle_MPZPara;Z* cos#theta;PZPara/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngle_MPZPara);
    TH2D* hist_CandCosDecayAngle_MPZPerp = new TH2D((title+"_CandCosDecayAngle_MPZPerp").c_str(), (title+"_CandCosDecayAngle_MPZPerp;Z* cos#theta;PZPerp/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngle_MPZPerp);
    TH2D* hist_CandCosDecayAngle_MPZParaLABMET = new TH2D((title+"_CandCosDecayAngle_MPZParaLABMET").c_str(), (title+"_CandCosDecayAngle_MPZParaLABMET;Z* cos#theta;PZParaLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngle_MPZParaLABMET);
    TH2D* hist_CandCosDecayAngle_MPZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_MPZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_MPZPerpLABMET;Z* cos#theta;PZPerpLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngle_MPZPerpLABMET);

    TH2D* hist_CandCosDecayAngleCM_MPZPara = new TH2D((title+"_CandCosDecayAngleCM_MPZPara").c_str(), (title+"_CandCosDecayAngleCM_MPZPara;Z* cos#theta CM;PZPara/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngleCM_MPZPara);
    TH2D* hist_CandCosDecayAngleCM_MPZPerp = new TH2D((title+"_CandCosDecayAngleCM_MPZPerp").c_str(), (title+"_CandCosDecayAngleCM_MPZPerp;Z* cos#theta CM;PZPerp/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngleCM_MPZPerp);
    TH2D* hist_CandCosDecayAngleCM_MPZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_MPZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_MPZParaLABMET;Z* cos#theta CM;PZParaLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngleCM_MPZParaLABMET);
    TH2D* hist_CandCosDecayAngleCM_MPZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_MPZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_MPZPerpLABMET;Z* cos#theta CM;PZPerpLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_CandCosDecayAngleCM_MPZPerpLABMET);

    TH2D* hist_Beta_MPZPara = new TH2D((title+"_Beta_MPZPara").c_str(), (title+"_Beta_MPZPara;Z* #beta;PZPara/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_Beta_MPZPara);
    TH2D* hist_Beta_MPZPerp = new TH2D((title+"_Beta_MPZPerp").c_str(), (title+"_Beta_MPZPerp;Z* #beta;PZPerp/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_Beta_MPZPerp);
    TH2D* hist_Beta_MPZParaLABMET = new TH2D((title+"_Beta_MPZParaLABMET").c_str(), (title+"_Beta_MPZParaLABMET;Z* #beta;PZParaLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_Beta_MPZParaLABMET);
    TH2D* hist_Beta_MPZPerpLABMET = new TH2D((title+"_Beta_MPZPerpLABMET").c_str(), (title+"_Beta_MPZPerpLABMET;Z* #beta;PZPerpLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_Beta_MPZPerpLABMET);

    TH2D* hist_BetaCM_MPZPara = new TH2D((title+"_BetaCM_MPZPara").c_str(), (title+"_BetaCM_MPZPara;Z* #beta^{CM};PZPara/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_BetaCM_MPZPara);
    TH2D* hist_BetaCM_MPZPerp = new TH2D((title+"_BetaCM_MPZPerp").c_str(), (title+"_BetaCM_MPZPerp;Z* #beta^{CM};PZPerp/M").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_BetaCM_MPZPerp);
    TH2D* hist_BetaCM_MPZParaLABMET = new TH2D((title+"_BetaCM_MPZParaLABMET").c_str(), (title+"_BetaCM_MPZParaLABMET;Z* #beta^{CM};PZParaLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_BetaCM_MPZParaLABMET);
    TH2D* hist_BetaCM_MPZPerpLABMET = new TH2D((title+"_BetaCM_MPZPerpLABMET").c_str(), (title+"_BetaCM_MPZPerpLABMET;Z* #beta^{CM};PZPerpLABMET/M_{T}").c_str(), g_NX, 0., 1., g_NX, 0., 2.);
    //hists2.push_back(hist_BetaCM_MPZPerpLABMET);

    TH2D* hist_PZAng_MPZPara = new TH2D((title+"_PZAng_MPZPara").c_str(), (title+"_PZAng_MPZPara;PZAng;PZPara/M").c_str(), g_NX, -1.6, 1.6, g_NX, 0., 2.);
    //hists2.push_back(hist_PZAng_MPZPara);
    TH2D* hist_PZAng_MPZPerp = new TH2D((title+"_PZAng_MPZPerp").c_str(), (title+"_PZAng_MPZPerp;PZAng;PZPerp/M").c_str(), g_NX, -1.6, 1.6, g_NX, 0., 2.);
    //hists2.push_back(hist_PZAng_MPZPerp);
    TH2D* hist_PZAng_MPZParaLABMET = new TH2D((title+"_PZAng_MPZParaLABMET").c_str(), (title+"_PZAng_MPZParaLABMET;PZAng;PZParaLABMET/M_{T}").c_str(), g_NX, -1.6, 1.6, g_NX, 0., 2.);
    //hists2.push_back(hist_PZAng_MPZParaLABMET);
    TH2D* hist_PZAng_MPZPerpLABMET = new TH2D((title+"_PZAng_MPZPerpLABMET").c_str(), (title+"_PZAng_MPZPerpLABMET;PZAng;PZPerpLABMET/M_{T}").c_str(), g_NX, -1.6, 1.6, g_NX, 0., 2.);
    //hists2.push_back(hist_PZAng_MPZPerpLABMET);

    //TH2D* hist_PZPara_CandECM = new TH2D((title+"_PZPara_CandECM").c_str(), (title+"_PZPara_CandECM;PZPara;Z* candidate E^{CM} [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_PZPara_CandECM);
    //TH2D* hist_PZPerp_CandECM = new TH2D((title+"_PZPerp_CandECM").c_str(), (title+"_PZPerp_CandECM;PZPerp;Z* candidate E^{CM} [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_PZPerp_CandECM);

    TH2D* hist_CandCosDecayAngle_RZPara = new TH2D((title+"_CandCosDecayAngle_RZPara").c_str(), (title+"_CandCosDecayAngle_RZPara;Z* candidate cos#theta;RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngle_RZPara);
    TH2D* hist_CandCosDecayAngle_RZPerp = new TH2D((title+"_CandCosDecayAngle_RZPerp").c_str(), (title+"_CandCosDecayAngle_RZPerp;Z* candidate cos#theta;RZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngle_RZPerp);
    TH2D* hist_CandCosDecayAngle_PZAng = new TH2D((title+"_CandCosDecayAngle_PZAng").c_str(), (title+"_CandCosDecayAngle_PZAng;Z* candidate cos#theta;PZAng").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandCosDecayAngle_PZAng);
    TH2D* hist_CandCosDecayAngle_PZPara = new TH2D((title+"_CandCosDecayAngle_PZPara").c_str(), (title+"_CandCosDecayAngle_PZPara;Z* candidate cos#theta;PZPara").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngle_PZPara);
    TH2D* hist_CandCosDecayAngle_PZPerp = new TH2D((title+"_CandCosDecayAngle_PZPerp").c_str(), (title+"_CandCosDecayAngle_PZPerp;Z* candidate cos#theta;PZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngle_PZPerp);

    TH2D* hist_CandCosDecayAngle_RZParaLABMET = new TH2D((title+"_CandCosDecayAngle_RZParaLABMET").c_str(), (title+"_CandCosDecayAngle_RZParaLABMET;Z* candidate cos#theta;RZParaLABMET").c_str(), g_NX, 0., 1., g_NX, -0.5, 0.5);
    //hists2.push_back(hist_CandCosDecayAngle_RZParaLABMET);
    TH2D* hist_CandCosDecayAngle_RZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_RZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_RZPerpLABMET;Z* candidate cos#theta;RZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngle_RZPerpLABMET);
    TH2D* hist_CandCosDecayAngle_PZAngLABMET = new TH2D((title+"_CandCosDecayAngle_PZAngLABMET").c_str(), (title+"_CandCosDecayAngle_PZAngLABMET;Z* candidate cos#theta;PZAngLABMET").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandCosDecayAngle_PZAngLABMET);
    TH2D* hist_CandCosDecayAngle_PZParaLABMET = new TH2D((title+"_CandCosDecayAngle_PZParaLABMET").c_str(), (title+"_CandCosDecayAngle_PZParaLABMET;Z* candidate cos#theta;PZParaLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngle_PZParaLABMET);
    TH2D* hist_CandCosDecayAngle_PZPerpLABMET = new TH2D((title+"_CandCosDecayAngle_PZPerpLABMET").c_str(), (title+"_CandCosDecayAngle_PZPerpLABMET;Z* candidate cos#theta;PZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngle_PZPerpLABMET);

    //TH2D* hist_CandCosDecayAngleCM_RZPara = new TH2D((title+"_CandCosDecayAngleCM_RZPara").c_str(), (title+"_CandCosDecayAngleCM_RZPara;Z* candidate cos#theta CM;RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngleCM_RZPara);
    //TH2D* hist_CandCosDecayAngleCM_RZPerp = new TH2D((title+"_CandCosDecayAngleCM_RZPerp").c_str(), (title+"_CandCosDecayAngleCM_RZPerp;Z* candidate cos#theta CM;RZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngleCM_RZPerp);
    //TH2D* hist_CandCosDecayAngleCM_PZAng = new TH2D((title+"_CandCosDecayAngleCM_PZAng").c_str(), (title+"_CandCosDecayAngleCM_PZAng;Z* candidate cos#theta CM;PZAng").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZAng);
    //TH2D* hist_CandCosDecayAngleCM_PZPara = new TH2D((title+"_CandCosDecayAngleCM_PZPara").c_str(), (title+"_CandCosDecayAngleCM_PZPara;Z* candidate cos#theta CM;PZPara").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZPara);
    //TH2D* hist_CandCosDecayAngleCM_PZPerp = new TH2D((title+"_CandCosDecayAngleCM_PZPerp").c_str(), (title+"_CandCosDecayAngleCM_PZPerp;Z* candidate cos#theta CM;PZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZPerp);

    //TH2D* hist_CandCosDecayAngleCM_RZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_RZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_RZParaLABMET;Z* candidate cos#theta CM;RZParaLABMET").c_str(), g_NX, 0., 1., g_NX, -0.5, 0.5);
    //hists2.push_back(hist_CandCosDecayAngleCM_RZParaLABMET);
    //TH2D* hist_CandCosDecayAngleCM_RZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_RZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_RZPerpLABMET;Z* candidate cos#theta CM;RZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_CandCosDecayAngleCM_RZPerpLABMET);
    //TH2D* hist_CandCosDecayAngleCM_PZAngLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZAngLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZAngLABMET;Z* candidate cos#theta CM;PZAngLABMET").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZAngLABMET);
    //TH2D* hist_CandCosDecayAngleCM_PZParaLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZParaLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZParaLABMET;Z* candidate cos#theta CM;PZParaLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZParaLABMET);
    //TH2D* hist_CandCosDecayAngleCM_PZPerpLABMET = new TH2D((title+"_CandCosDecayAngleCM_PZPerpLABMET").c_str(), (title+"_CandCosDecayAngleCM_PZPerpLABMET;Z* candidate cos#theta CM;PZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_CandCosDecayAngleCM_PZPerpLABMET);

    TH2D* hist_CandCosDecayAngle_CandCosDecayAngleCM = new TH2D((title+"_CandCosDecayAngle_CandCosDecayAngleCM").c_str(), (title+"_CandCosDecayAngle_CandCosDecayAngleCM;Z* candidate cos#theta;Z* candidate cos#theta CM").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_CandCosDecayAngleCM);
    TH2D* hist_CandBeta_CandBetaCM = new TH2D((title+"_CandBeta_CandBetaCM").c_str(), (title+"_CandBeta_CandBetaCM;Z* candidate #beta;Z* candidate #beta^{CM}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandBeta_CandBetaCM);

    //TH2D* hist_Beta_RZPara = new TH2D((title+"_Beta_RZPara").c_str(), (title+"_Beta_RZPara;Z* candidate #beta;RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_Beta_RZPara);
    //TH2D* hist_Beta_RZPerp = new TH2D((title+"_Beta_RZPerp").c_str(), (title+"_Beta_RZPerp;Z* candidate #beta;RZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_Beta_RZPerp);
    //TH2D* hist_Beta_PZAng = new TH2D((title+"_Beta_PZAng").c_str(), (title+"_Beta_PZAng;Z* candidate #beta;PZAng").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_Beta_PZAng);
    //TH2D* hist_Beta_PZPara = new TH2D((title+"_Beta_PZPara").c_str(), (title+"_Beta_PZPara;Z* candidate #beta;PZPara").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_Beta_PZPara);
    //TH2D* hist_Beta_PZPerp = new TH2D((title+"_Beta_PZPerp").c_str(), (title+"_Beta_PZPerp;Z* candidate #beta;PZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_Beta_PZPerp);

    //TH2D* hist_Beta_RZParaLABMET = new TH2D((title+"_Beta_RZParaLABMET").c_str(), (title+"_Beta_RZParaLABMET;Z* candidate #beta;RZParaLABMET").c_str(), g_NX, 0., 1., g_NX, -0.5, 0.5);
    //hists2.push_back(hist_Beta_RZParaLABMET);
    //TH2D* hist_Beta_RZPerpLABMET = new TH2D((title+"_Beta_RZPerpLABMET").c_str(), (title+"_Beta_RZPerpLABMET;Z* candidate #beta;RZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_Beta_RZPerpLABMET);
    //TH2D* hist_Beta_PZAngLABMET = new TH2D((title+"_Beta_PZAngLABMET").c_str(), (title+"_Beta_PZAngLABMET;Z* candidate #beta;PZAngLABMET").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_Beta_PZAngLABMET);
    //TH2D* hist_Beta_PZParaLABMET = new TH2D((title+"_Beta_PZParaLABMET").c_str(), (title+"_Beta_PZParaLABMET;Z* candidate #beta;PZParaLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_Beta_PZParaLABMET);
    //TH2D* hist_Beta_PZPerpLABMET = new TH2D((title+"_Beta_PZPerpLABMET").c_str(), (title+"_Beta_PZPerpLABMET;Z* candidate #beta;PZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_Beta_PZPerpLABMET);

    //TH2D* hist_BetaCM_RZPara = new TH2D((title+"_BetaCM_RZPara").c_str(), (title+"_BetaCM_RZPara;Z* candidate #beta^{CM};RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_BetaCM_RZPara);
    //TH2D* hist_BetaCM_RZPerp = new TH2D((title+"_BetaCM_RZPerp").c_str(), (title+"_BetaCM_RZPerp;Z* candidate #beta^{CM};RZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_BetaCM_RZPerp);
    //TH2D* hist_BetaCM_PZAng = new TH2D((title+"_BetaCM_PZAng").c_str(), (title+"_BetaCM_PZAng;Z* candidate #beta^{CM};PZAng").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_BetaCM_PZAng);
    //TH2D* hist_BetaCM_PZPara = new TH2D((title+"_BetaCM_PZPara").c_str(), (title+"_BetaCM_PZPara;Z* candidate #beta^{CM};PZPara").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_BetaCM_PZPara);
    //TH2D* hist_BetaCM_PZPerp = new TH2D((title+"_BetaCM_PZPerp").c_str(), (title+"_BetaCM_PZPerp;Z* candidate #beta^{CM};PZPerp").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_BetaCM_PZPerp);

    //TH2D* hist_BetaCM_RZParaLABMET = new TH2D((title+"_BetaCM_RZParaLABMET").c_str(), (title+"_BetaCM_RZParaLABMET;Z* candidate #beta^{CM};RZParaLABMET").c_str(), g_NX, 0., 1., g_NX, -0.5, 0.5);
    //hists2.push_back(hist_BetaCM_RZParaLABMET);
    //TH2D* hist_BetaCM_RZPerpLABMET = new TH2D((title+"_BetaCM_RZPerpLABMET").c_str(), (title+"_BetaCM_RZPerpLABMET;Z* candidate #beta^{CM};RZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_BetaCM_RZPerpLABMET);
    //TH2D* hist_BetaCM_PZAngLABMET = new TH2D((title+"_BetaCM_PZAngLABMET").c_str(), (title+"_BetaCM_PZAngLABMET;Z* candidate #beta^{CM};PZAngLABMET").c_str(), g_NX, 0., 1., g_NX, -1.6, 1.6);
    //hists2.push_back(hist_BetaCM_PZAngLABMET);
    //TH2D* hist_BetaCM_PZParaLABMET = new TH2D((title+"_BetaCM_PZParaLABMET").c_str(), (title+"_BetaCM_PZParaLABMET;Z* candidate #beta^{CM};PZParaLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_BetaCM_PZParaLABMET);
    //TH2D* hist_BetaCM_PZPerpLABMET = new TH2D((title+"_BetaCM_PZPerpLABMET").c_str(), (title+"_BetaCM_PZPerpLABMET;Z* candidate #beta^{CM};PZPerpLABMET").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    //hists2.push_back(hist_BetaCM_PZPerpLABMET);

    TH2D* hist_MRZPara_RZPara = new TH2D((title+"_MRZPara_RZPara").c_str(), (title+"_MRZPara_RZPara;M/RZPara;RZPara").c_str(), g_NX, 0., 1000., g_NX, 0., 0.5);
    //hists2.push_back(hist_MRZPara_RZPara);
    TH2D* hist_MRZPara_CandCosDecayAngleCM = new TH2D((title+"_MRZPara_CandCosDecayAngleCM").c_str(), (title+"_MRZPara_CandCosDecayAngleCM;M/RZPara;Z* cos#theta CM").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    //hists2.push_back(hist_MRZPara_CandCosDecayAngleCM);
    TH2D* hist_MRZPara_BetaCM = new TH2D((title+"_MRZPara_BetaCM").c_str(), (title+"_MRZPara_BetaCM;M/RZPara;Z* #beta^{CM}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    //hists2.push_back(hist_MRZPara_BetaCM);
    TH2D* hist_MRZPara_CandPZzLAB = new TH2D((title+"_MRZPara_CandPZzLAB").c_str(), (title+"_MRZPara_CandPZzLAB;M/RZPara;PZ_{z}^{LAB}").c_str(), g_NX, 0., 1000., g_NX, 0., 500.);
    //hists2.push_back(hist_MRZPara_CandPZzLAB);
    TH2D* hist_RZPara_CandPZzLAB = new TH2D((title+"_RZPara_CandPZzLAB").c_str(), (title+"_RZPara_CandPZzLAB;RZPara;PZ_{z}^{LAB}").c_str(), g_NX, 0., 0.5, g_NX, 0., 500.);
    //hists2.push_back(hist_RZPara_CandPZzLAB);
    TH2D* hist_RISR_RZPara = new TH2D((title+"_RISR_RZPara").c_str(), (title+"_RISR_RZPara;R_{ISR};RZPara").c_str(), g_NX, 0.5, 1., g_NX, 0., 0.5);
    //hists2.push_back(hist_RISR_RZPara);
    TH2D* hist_RISR_MRZPara = new TH2D((title+"_RISR_MRZPara").c_str(), (title+"_RISR_MRZPara;R_{ISR};M/RZPara").c_str(), g_NX, 0.5, 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_RISR_MRZPara);
    TH2D* hist_RZPara_BetaZLAB = new TH2D((title+"_RZPara_BetaZLAB").c_str(), (title+"_RZPara_BetaZLAB;RZPara;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 0.5, g_NX, 0., 1.);
    //hists2.push_back(hist_RZPara_BetaZLAB);
    TH2D* hist_MRZPara_BetaZLAB = new TH2D((title+"_MRZPara_BetaZLAB").c_str(), (title+"_MRZPara_BetaZLAB;M/RZPara;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    //hists2.push_back(hist_MRZPara_BetaZLAB);
    TH2D* hist_RZPerp_CandPZzLAB = new TH2D((title+"_RZPerp_CandPZzLAB").c_str(), (title+"_RZPerp_CandPZzLAB;M/RZPerp;PZ_{z}^{LAB}").c_str(), g_NX, 0., 0.5, g_NX, 0., 500.);
    //hists2.push_back(hist_RZPerp_CandPZzLAB);
    TH2D* hist_RISR_CandPZzLAB = new TH2D((title+"_RISR_CandPZzLAB").c_str(), (title+"_RISR_CandPZzLAB;R_{ISR};CandPZzLAB").c_str(), g_NX, 0.5, 1., g_NX, 0., 500.);
    //hists2.push_back(hist_RISR_CandPZzLAB);
    TH2D* hist_RISR_BetaZLAB = new TH2D((title+"_RISR_BetaZLAB").c_str(), (title+"_RISR_BetaZLAB;R_{ISR};#beta_{Z}^{LAB}").c_str(), g_NX, 0.5, 1.0, g_NX, 0., 1.0);
    //hists2.push_back(hist_RISR_BetaZLAB);
    TH2D* hist_CandML_BetaZLAB = new TH2D((title+"_CandML_BetaZLAB").c_str(), (title+"_CandML_BetaZLAB;Z* candidate mass [GeV];#beta_{Z}^{LAB}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_BetaZLAB);
    TH2D* hist_RZPerp_BetaZLAB = new TH2D((title+"_RZPerp_BetaZLAB").c_str(), (title+"_RZPerp_BetaZLAB;RZPerp;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 0.5, g_NX, 0., 1.);
    //hists2.push_back(hist_RZPerp_BetaZLAB);

    TH2D* hist_CandCosDecayAngle_BetaZLAB_PreCut = new TH2D((title+"_CandCosDecayAngle_BetaZLAB_PreCut").c_str(), (title+"_CandCosDecayAngle_BetaZLAB_PreCut;Z* cos#theta;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_BetaZLAB_PreCut);
    TH2D* hist_CandCosDecayAngle_BetaZLAB = new TH2D((title+"_CandCosDecayAngle_BetaZLAB").c_str(), (title+"_CandCosDecayAngle_BetaZLAB;Z* cos#theta;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_BetaZLAB);
    TH2D* hist_CandCosDecayAngleCM_BetaZLAB = new TH2D((title+"_CandCosDecayAngleCM_BetaZLAB").c_str(), (title+"_CandCosDecayAngleCM_BetaZLAB;Z* cos#theta CM;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngleCM_BetaZLAB);

    TH2D* hist_RZParaM_RZPara = new TH2D((title+"_RZParaM_RZPara").c_str(), (title+"_RZParaM_RZPara;RZPara/M;RZPara").c_str(), g_NX, 0., 0.02, g_NX, 0., 0.5);
    //hists2.push_back(hist_RZParaM_RZPara);
    TH2D* hist_RZParaM_CandCosDecayAngleCM = new TH2D((title+"_RZParaM_CandCosDecayAngleCM").c_str(), (title+"_RZParaM_CandCosDecayAngleCM;RZPara/M;Z* cos#theta CM").c_str(), g_NX, 0., 0.02, g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaM_CandCosDecayAngleCM);
    TH2D* hist_RZParaM_BetaCM = new TH2D((title+"_RZParaM_BetaCM").c_str(), (title+"_RZParaM_BetaCM;RZPara/M;Z* #beta^{CM}").c_str(), g_NX, 0., 0.02, g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaM_BetaCM);
    TH2D* hist_RZParaM_CandPZzLAB = new TH2D((title+"_RZParaM_CandPZzLAB").c_str(), (title+"_RZParaM_CandPZzLAB;RZPara/M;PZ_{z}^{LAB}").c_str(), g_NX, 0., 0.02, g_NX, 0., 500.);
    //hists2.push_back(hist_RZParaM_CandPZzLAB);
    TH2D* hist_RZParaM_BetaZLAB = new TH2D((title+"_RZParaM_BetaZLAB").c_str(), (title+"_RZParaM_BetaZLAB;RZPara/M;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 0.02, g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaM_BetaZLAB);

    //TH2D* hist_MRZPara_OneMinusRZPara = new TH2D((title+"_MRZPara_OneMinusRZPara").c_str(), (title+"_MRZPara_OneMinusRZPara;M/RZPara;1-RZPara").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    //hists2.push_back(hist_MRZPara_OneMinusRZPara);
    //TH2D* hist_OneMinusRZPara_CandPZzLAB = new TH2D((title+"_OneMinusRZPara_CandPZzLAB").c_str(), (title+"_OneMinusRZPara_CandPZzLAB;1-RZPara;PZ_{z}^{LAB}").c_str(), g_NX, 0., 1., g_NX, 0., 500.);
    //hists2.push_back(hist_OneMinusRZPara_CandPZzLAB);
    //TH2D* hist_RISR_OneMinusRZPara = new TH2D((title+"_RISR_OneMinusRZPara").c_str(), (title+"_RISR_OneMinusRZPara;R_{ISR};1-RZPara").c_str(), g_NX, 0.5, 1., g_NX, 0., 1.);
    //hists2.push_back(hist_RISR_OneMinusRZPara);
    //TH2D* hist_OneMinusRZPara_BetaZLAB = new TH2D((title+"_OneMinusRZPara_BetaZLAB").c_str(), (title+"_OneMinusRZPara_BetaZLAB;1-RZPara;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_OneMinusRZPara_BetaZLAB);
    //TH2D* hist_BetaCM_OneMinusRZPara = new TH2D((title+"_BetaCM_OneMinusRZPara").c_str(), (title+"_BetaCM_OneMinusRZPara;Z* candidate #beta^{CM};1-RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_BetaCM_OneMinusRZPara);
    //TH2D* hist_CandCosDecayAngle_OneMinusRZPara = new TH2D((title+"_CandCosDecayAngle_OneMinusRZPara").c_str(), (title+"_CandCosDecayAngle_OneMinusRZPara;Z* candidate cos#theta;1-RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_OneMinusRZPara);
    //TH2D* hist_Beta_OneMinusRZPara = new TH2D((title+"_Beta_OneMinusRZPara").c_str(), (title+"_Beta_OneMinusRZPara;Z* candidate #beta;1-RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_Beta_OneMinusRZPara);
    //TH2D* hist_CandCosDecayAngleCM_OneMinusRZPara = new TH2D((title+"_CandCosDecayAngleCM_OneMinusRZPara").c_str(), (title+"_CandCosDecayAngleCM_OneMinusRZPara;Z* candidate cos#theta CM;1-RZPara").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngleCM_OneMinusRZPara);

    //TH2D* hist_MRZPara_CandMdivPTISR = new TH2D((title+"_MRZPara_CandMdivPTISR").c_str(), (title+"_MRZPara_CandMdivPTISR;M/RZPara;CandM/PTISR").c_str(), g_NX, 0., 1000., g_NX, 0., 0.8);
    //hists2.push_back(hist_MRZPara_CandMdivPTISR);
    //TH2D* hist_CandMdivPTISR_CandPZzLAB = new TH2D((title+"_CandMdivPTISR_CandPZzLAB").c_str(), (title+"_CandMdivPTISR_CandPZzLAB;CandM/PTISR;PZ_{z}^{LAB}").c_str(), g_NX, 0., 0.8, g_NX, 0., 500.);
    //hists2.push_back(hist_CandMdivPTISR_CandPZzLAB);
    //TH2D* hist_RISR_CandMdivPTISR = new TH2D((title+"_RISR_CandMdivPTISR").c_str(), (title+"_RISR_CandMdivPTISR;R_{ISR};CandM/PTISR").c_str(), g_NX, 0.5, 1., g_NX, 0., 1.);
    //hists2.push_back(hist_RISR_CandMdivPTISR);
    //TH2D* hist_CandMdivPTISR_BetaZLAB = new TH2D((title+"_CandMdivPTISR_BetaZLAB").c_str(), (title+"_CandMdivPTISR_BetaZLAB;CandM/PTISR;#beta_{Z}^{LAB}").c_str(), g_NX, 0., .8, g_NX, 0., 1.);
    //hists2.push_back(hist_CandMdivPTISR_BetaZLAB);
    //TH2D* hist_BetaCM_CandMdivPTISR = new TH2D((title+"_BetaCM_CandMdivPTISR").c_str(), (title+"_BetaCM_CandMdivPTISR;Z* candidate #beta^{CM};CandM/PTISR").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_BetaCM_CandMdivPTISR);
    //TH2D* hist_CandCosDecayAngle_CandMdivPTISR = new TH2D((title+"_CandCosDecayAngle_CandMdivPTISR").c_str(), (title+"_CandCosDecayAngle_CandMdivPTISR;Z* candidate cos#theta;CandM/PTISR").c_str(), g_NX, 0., 0.8, g_NX, 0., 1.);
    //hists2.push_back(hist_CandCosDecayAngle_CandMdivPTISR);
    //TH2D* hist_Beta_CandMdivPTISR = new TH2D((title+"_Beta_CandMdivPTISR").c_str(), (title+"_Beta_CandMdivPTISR;Z* candidate #beta;CandM/PTISR").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_Beta_CandMdivPTISR);
    //TH2D* hist_CandCosDecayAngleCM_CandMdivPTISR = new TH2D((title+"_CandCosDecayAngleCM_CandMdivPTISR").c_str(), (title+"_CandCosDecayAngleCM_CandMdivPTISR;Z* candidate cos#theta CM;CandM/PTISR").c_str(), g_NX, 0., 1., g_NX, 0., 0.8);
    //hists2.push_back(hist_CandCosDecayAngleCM_CandMdivPTISR);

    TH2D* hist_MRZPara_MInv = new TH2D((title+"_MRZPara_MInv").c_str(), (title+"_MRZPara_MInv;M/RZPara;MInv").c_str(), g_NX, 0., 1000., g_NX, 0., 1000.);
    //hists2.push_back(hist_MRZPara_MInv);
    TH2D* hist_MInv_CandPZzLAB = new TH2D((title+"_MInv_CandPZzLAB").c_str(), (title+"_MInv_CandPZzLAB;MInv;PZ_{z}^{LAB}").c_str(), g_NX, 0., 1000., g_NX, 0., 500.);
    //hists2.push_back(hist_MInv_CandPZzLAB);
    TH2D* hist_RISR_MInv = new TH2D((title+"_RISR_MInv").c_str(), (title+"_RISR_MInv;R_{ISR};MInv").c_str(), g_NX, 0.5, 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_RISR_MInv);
    TH2D* hist_MInv_BetaZLAB = new TH2D((title+"_MInv_BetaZLAB").c_str(), (title+"_MInv_BetaZLAB;MInv;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    //hists2.push_back(hist_MInv_BetaZLAB);
    TH2D* hist_BetaCM_MInv = new TH2D((title+"_BetaCM_MInv").c_str(), (title+"_BetaCM_MInv;Z* candidate #beta^{CM};MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_BetaCM_MInv);
    TH2D* hist_CandCosDecayAngle_MInv = new TH2D((title+"_CandCosDecayAngle_MInv").c_str(), (title+"_CandCosDecayAngle_MInv;Z* candidate cos#theta;MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_CandCosDecayAngle_MInv);
    TH2D* hist_Beta_MInv = new TH2D((title+"_Beta_MInv").c_str(), (title+"_Beta_MInv;Z* candidate #beta;MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_Beta_MInv);
    TH2D* hist_CandCosDecayAngleCM_MInv = new TH2D((title+"_CandCosDecayAngleCM_MInv").c_str(), (title+"_CandCosDecayAngleCM_MInv;Z* candidate cos#theta CM;MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_CandCosDecayAngleCM_MInv);

    //TH2D* hist_CandML_RZParaM = new TH2D((title+"_CandML_RZParaM").c_str(), (title+"_CandML_RZParaM;Z* candidate M [GeV];RZPara/M").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 0.02);
    //hists2.push_back(hist_CandML_RZParaM);
    //TH2D* hist_CandML_OneMinusRZPara = new TH2D((title+"_CandML_OneMinusRZPara").c_str(), (title+"_CandML_OneMinusRZPara;Z* candidate M [GeV];1-RZPara").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_OneMinusRZPara);
    //TH2D* hist_CandML_CandMdivPTISR = new TH2D((title+"_CandML_CandMdivPTISR").c_str(), (title+"_CandML_CandMdivPTISR;Z* candidate M [GeV];M/PTISR").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 0.8);
    //hists2.push_back(hist_CandML_CandMdivPTISR);
    TH2D* hist_CandML_MInv = new TH2D((title+"_CandML_MInv").c_str(), (title+"_CandML_MInv;Z* candidate M [GeV];MInv").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1000.);
    //hists2.push_back(hist_CandML_MInv);

    //TH2D* hist_RZParaDivRISR_RZParaM = new TH2D((title+"_RZParaDivRISR_RZParaM").c_str(), (title+"_RZParaDivRISR_RZParaM;RZPara/RISR;RZPara/M").c_str(), g_NX, 0., 2., g_NX, 0., 0.02);
    //hists2.push_back(hist_RZParaDivRISR_RZParaM);
    //TH2D* hist_RZParaDivRISR_CandMdivPTISR = new TH2D((title+"_RZParaDivRISR_CandMdivPTISR").c_str(), (title+"_RZParaDivRISR_CandMdivPTISR;RZPara/RISR;M/PTISR").c_str(), g_NX, 0., 2., g_NX, 0., 0.8);
    //hists2.push_back(hist_RZParaDivRISR_CandMdivPTISR);
    //TH2D* hist_RZParaDivRISR_MInv = new TH2D((title+"_RZParaDivRISR_MInv").c_str(), (title+"_RZParaDivRISR_MInv;RZPara/RISR;MInv").c_str(), g_NX, 0., 2., g_NX, 0., 1000.);
    //hists2.push_back(hist_RZParaDivRISR_MInv);
    //TH2D* hist_RZParaDivRISR_BetaZLAB = new TH2D((title+"_RZParaDivRISR_BetaZLAB").c_str(), (title+"_RZParaDivRISR_BetaZLAB;RZPara/RISR;#beta_{Z}^{LAB}").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_BetaZLAB);
    //TH2D* hist_RZParaDivRISR_MRZPara = new TH2D((title+"_RZParaDivRISR_MRZPara").c_str(), (title+"_RZParaDivRISR_MRZPara;RZPara/RISR;M/RZPara").c_str(), g_NX, 0., 2., g_NX, 0., 1000.);
    //hists2.push_back(hist_RZParaDivRISR_MRZPara);
    //TH2D* hist_RZParaDivRISR_RZPara = new TH2D((title+"_RZParaDivRISR_RZPara").c_str(), (title+"_RZParaDivRISR_RZPara;RZPara/RISR;RZPara").c_str(), g_NX, 0., 2., g_NX, 0., 0.5);
    //hists2.push_back(hist_RZParaDivRISR_RZPara);
    //TH2D* hist_RZParaDivRISR_OneMinusRZPara = new TH2D((title+"_RZParaDivRISR_OneMinusRZPara").c_str(), (title+"_RZParaDivRISR_OneMinusRZPara;RZPara/RISR;1-RZPara").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_OneMinusRZPara);
    //TH2D* hist_RZParaDivRISR_Beta = new TH2D((title+"_RZParaDivRISR_Beta").c_str(), (title+"_RZParaDivRISR_Beta;RZPara/RISR;#beta^{LAB}").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_Beta);
    //TH2D* hist_RZParaDivRISR_BetaCM = new TH2D((title+"_RZParaDivRISR_BetaCM").c_str(), (title+"_RZParaDivRISR_BetaCM;RZPara/RISR;#beta^{CM}").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_BetaCM);
    //TH2D* hist_RZParaDivRISR_CandCosDecayAngle = new TH2D((title+"_RZParaDivRISR_CandCosDecayAngle").c_str(), (title+"_RZParaDivRISR_CandCosDecayAngle;RZPara/RISR;cos#theta").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_CandCosDecayAngle);
    //TH2D* hist_RZParaDivRISR_CandCosDecayAngleCM = new TH2D((title+"_RZParaDivRISR_CandCosDecayAngleCM").c_str(), (title+"_RZParaDivRISR_CandCosDecayAngleCM;RZPara/RISR;cos#theta^{CM}").c_str(), g_NX, 0., 2., g_NX, 0., 1.);
    //hists2.push_back(hist_RZParaDivRISR_CandCosDecayAngleCM);
    //TH2D* hist_CandMLDivRISR_CandML = new TH2D((title+"_CandMLDivRISR_CandML").c_str(), (title+"_CandMLDivRISR_CandML;CandML/RISR;CandML").c_str(), g_NX, 0., 2., g_NX, CandMinMass, CandMaxMass);
    //hists2.push_back(hist_CandMLDivRISR_CandML);

    TH2D* hist_DecayAngleDiff_PzCM = new TH2D((title+"_DecayAngleDiff_PzCM").c_str(), (title+"_DecayAngleDiff_PzCM;|cos#theta CM| - |cos#theta|;PzCM").c_str(), g_NX, -1., 1., g_NX, 0., 500.);
    //hists2.push_back(hist_DecayAngleDiff_PzCM);
    TH2D* hist_BetaDiff_PzCM = new TH2D((title+"_BetaDiff_PzCM").c_str(), (title+"_BetaDiff_PzCM;#beta CM - #beta;PzCM").c_str(), g_NX, -1., 1., g_NX, 0., 500.);
    //hists2.push_back(hist_BetaDiff_PzCM);
    TH2D* hist_DecayAngleDiff_BetaDiff = new TH2D((title+"_DecayAngleDiff_BetaDiff").c_str(), (title+"_DecayAngleDiff_BetaDiff;|cos#theta CM| - |cos#theta|;#beta CM - #beta").c_str(), g_NX, -1., 1., g_NX, -1., 1.);
    //hists2.push_back(hist_DecayAngleDiff_BetaDiff);

    TH2D* hist_CandML_MAXcandML = new TH2D((title+"_CandML_MAXcandML").c_str(), (title+"_CandML_MAXcandML;MIN Cand ML [GeV];MAX Cand ML [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, CandMinMass, CandMaxMass);
    //hists2.push_back(hist_CandML_MAXcandML);
    TH2D* hist_BetaZLAB_MAXcandBetaZLAB = new TH2D((title+"_BetaZLAB_MAXcandBetaZLAB").c_str(), (title+"_BetaZLAB_MAXcandBetaZLAB;MIN Cand #beta^{LAB}_{Z};MAX Cand #beta^{LAB}_{Z}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_BetaZLAB_MAXcandBetaZLAB);
    TH2D* hist_CosDecayAngle_MAXcandCosDecayAngle = new TH2D((title+"_CosDecayAngle_MAXcandCosDecayAngle").c_str(), (title+"_CosDecayAngle_MAXcandCosDecayAngle;MIN Cand cos#theta;MAX Cand cos#theta").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CosDecayAngle_MAXcandCosDecayAngle);
    TH2D* hist_MAXcandBetaZLAB_MAXcandCosDecayAngle = new TH2D((title+"_MAXcandBetaZLAB_MAXcandCosDecayAngle").c_str(), (title+"_MAXcandBetaZLAB_MAXcandCosDecayAngle;MAX Cand #beta^{LAB}_{Z};MAX Cand cos#theta").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_MAXcandBetaZLAB_MAXcandCosDecayAngle);
    TH2D* hist_MAXcandML_MT2 = new TH2D((title+"_MAXcandML_MT2").c_str(), (title+"_MAXcandML_MT2;MAX Cand ML [GeV];MT2 [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 250.);
    //hists2.push_back(hist_MAXcandML_MT2);
    TH2D* hist_MAXcandML_MConCM = new TH2D((title+"_MAXcandML_MConCM").c_str(), (title+"_MAXcandML_MConCM;MAX Cand ML [GeV];MCon^{CM} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_MAXcandML_MConCM);
    TH2D* hist_MAXcandML_MConLEPS = new TH2D((title+"_MAXcandML_MConLEPS").c_str(), (title+"_MAXcandML_MConLEPS;MAX Cand ML [GeV];MCon^{LEPS} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_MAXcandML_MConLEPS);
    TH2D* hist_MT2_MConCM = new TH2D((title+"_MT2_MConCM").c_str(), (title+"_MT2_MConCM;MT2 [GeV];MCon^{CM} [GeV]").c_str(), g_NX, 0., 250., g_NX, 0., 400.);
    //hists2.push_back(hist_MT2_MConCM);
    TH2D* hist_MT2_MConLEPS = new TH2D((title+"_MT2_MConLEPS").c_str(), (title+"_MT2_MConLEPS;MT2 [GeV];MCon^{CM} [GeV]").c_str(), g_NX, 0., 250., g_NX, 0., 400.);
    //hists2.push_back(hist_MT2_MConLEPS);
    TH2D* hist_MConCM_MConLEPS = new TH2D((title+"_MConCM_MConLEPS").c_str(), (title+"_MConCM_MConLEPS;MConCM [GeV];MCon^{LEPS} [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_MConCM_MConLEPS);
    TH2D* hist_CandML_MT2 = new TH2D((title+"_CandML_MT2").c_str(), (title+"_CandML_MT2;MIN Cand ML [GeV];MT2 [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 250.);
    //hists2.push_back(hist_CandML_MT2);
    TH2D* hist_CandML_MConCM = new TH2D((title+"_CandML_MConCM").c_str(), (title+"_CandML_MConCM;MIN Cand ML [GeV];MCon^{CM} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_CandML_MConCM);
    TH2D* hist_CandML_MConLEPS = new TH2D((title+"_CandML_MConLEPS").c_str(), (title+"_CandML_MConLEPS;MIN Cand ML [GeV];MCon^{LEPS} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_CandML_MConLEPS);

    TH2D* hist_CandML_MT2MAX = new TH2D((title+"_CandML_MT2MAX").c_str(), (title+"_CandML_MT2MAX;MIN CandML [GeV];MT2 MAX [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 250.);
    //hists2.push_back(hist_CandML_MT2MAX);
    TH2D* hist_CandML_MConLEPSMET = new TH2D((title+"_CandML_MConLEPSMET").c_str(), (title+"_CandML_MConLEPSMET;MIN CandML [GeV];MCon^{LEPS MET} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_CandML_MConLEPSMET);
    TH2D* hist_CandML_CosThetaRazor = new TH2D((title+"_CandML_CosThetaRazor").c_str(), (title+"_CandML_CosThetaRazor;MIN CandML [GeV];cos#theta^{Razor}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_CosThetaRazor);
    TH2D* hist_CandML_MAXCosThetaRazor = new TH2D((title+"_CandML_MAXCosThetaRazor").c_str(), (title+"_CandML_MAXCosThetaRazor;MIN CandML [GeV];cos#theta^{MAXRazor}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_CandML_MAXCosThetaRazor);
    TH2D* hist_CandML_MAXMInv = new TH2D((title+"_CandML_MAXMInv").c_str(), (title+"_CandML_MAXMInv;MIN CandML [GeV];MAX MInv [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1000.);
    //hists2.push_back(hist_CandML_MAXMInv);
    TH2D* hist_MAXCandML_MT2MAX = new TH2D((title+"_MAXCandML_MT2MAX").c_str(), (title+"_MAXCandML_MT2MAX;MAX CandML [GeV];MAX MT2 [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 250.);
    //hists2.push_back(hist_MAXCandML_MT2MAX);
    TH2D* hist_MAXCandML_MCon_LEPSMET = new TH2D((title+"_MAXCandML_MConLEPSMET").c_str(), (title+"_MAXCandML_MCon_LEPSMET;MAX CandML [GeV];MCon^{LEPS MET} [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 400.);
    //hists2.push_back(hist_MAXCandML_MCon_LEPSMET);
    TH2D* hist_MAXCandML_CosThetaRazor = new TH2D((title+"_MAXCandML_CosThetaRazor").c_str(), (title+"_MAXCandML_CosThetaRazor;MAX CandML [GeV];cos#theta^{Razor}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_MAXCandML_CosThetaRazor);
    TH2D* hist_MAXCandML_MAXCosThetaRazor = new TH2D((title+"_MAXCandML_MAXCosThetaRazor").c_str(), (title+"_MAXCandML_MAXCosThetaRazor;MAX CandML [GeV];cos#theta^{MAXRazor}").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_MAXCandML_MAXCosThetaRazor);
    TH2D* hist_MAXCandML_MAXMInv = new TH2D((title+"_MAXCandML_MAXMInv").c_str(), (title+"_MAXCandML_MAXMInv;MAX CandML [GeV];MAX MInv [GeV]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1000.);
    //hists2.push_back(hist_MAXCandML_MAXMInv);
    TH2D* hist_MAXCandML_RISR = new TH2D((title+"_MAXCandML_RISR").c_str(), (title+"_MAXCandML_RISR;MAX CandML [GeV];RISR").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0.5, 1.);
    //hists2.push_back(hist_MAXCandML_RISR);
    TH2D* hist_MAXCandML_MInv = new TH2D((title+"_MAXCandML_MInv").c_str(), (title+"_MAXCandML_MInv;MAX CandML [GeV];MInv [GeV}]").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1000.);
    //hists2.push_back(hist_MAXCandML_MInv);
    TH2D* hist_MConCM_MT2MAX = new TH2D((title+"_MConCM_MT2MAX").c_str(), (title+"_MConCM_MT2MAX;M_{Con}^{CM} [GeV];MAX MT2 [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 250.);
    //hists2.push_back(hist_MConCM_MT2MAX);
    TH2D* hist_MConCM_MConLEPSMET = new TH2D((title+"_MConCM_MConLEPSMET").c_str(), (title+"_MConCM_MConLEPSMET;M_{Con}^{CM} [GeV];MCon^{LEPS MET} [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_MConCM_MConLEPSMET);
    TH2D* hist_MConCM_CosThetaRazor = new TH2D((title+"_MConCM_CosThetaRazor").c_str(), (title+"_MConCM_CosThetaRazor;M_{Con}^{CM} [GeV];cos#theta^{Razor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConCM_CosThetaRazor);
    TH2D* hist_MConCM_MAXCosThetaRazor = new TH2D((title+"_MConCM_MAXCosThetaRazor").c_str(), (title+"_MConCM_MAXCosThetaRazor;M_{Con}^{CM} [GeV];cos#theta^{MAXRazor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConCM_MAXCosThetaRazor);
    TH2D* hist_MConCM_MAXMInv = new TH2D((title+"_MConCM_MAXMInv").c_str(), (title+"_MConCM_MAXMInv;M_{Con}^{CM} [GeV];MAX MInv [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConCM_MAXMInv);
    TH2D* hist_MConCM_RISR = new TH2D((title+"_MConCM_RISR").c_str(), (title+"_MConCM_RISR;M_{Con}^{CM} [GeV];RISR").c_str(), g_NX, 0., 400., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MConCM_RISR);
    TH2D* hist_MConCM_MInv = new TH2D((title+"_MConCM_MInv").c_str(), (title+"_MConCM_MInv;M_{Con}^{CM} [GeV];MInv [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConCM_MInv);
    TH2D* hist_MConLEPS_MT2MAX = new TH2D((title+"_MConLEPS_MT2MAX").c_str(), (title+"_MConLEPS_MT2MAX;M_{Con}^{LEPS} [GeV];MAX MT2 [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 250.);
    //hists2.push_back(hist_MConLEPS_MT2MAX);
    TH2D* hist_MConLEPS_MConLEPSMET = new TH2D((title+"_MConLEPS_MConLEPSMET").c_str(), (title+"_MConLEPS_MConLEPSMET;M_{Con}^{LEPS} [GeV];M_{Con}^{LEPS MET} [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    //hists2.push_back(hist_MConLEPS_MConLEPSMET);
    TH2D* hist_MConLEPS_CosThetaRazor = new TH2D((title+"_MConLEPS_CosThetaRazor").c_str(), (title+"_MConLEPS_CosThetaRazor;M_{Con}^{LEPS} [GeV];cos#theta^{Razor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConLEPS_CosThetaRazor);
    TH2D* hist_MConLEPS_MAXCosThetaRazor = new TH2D((title+"_MConLEPS_MAXCosThetaRazor").c_str(), (title+"_MConLEPS_MAXCosThetaRazor;M_{Con}^{LEPS} [GeV];cos#theta^{MAXRazor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConLEPS_MAXCosThetaRazor);
    TH2D* hist_MConLEPS_MAXMInv = new TH2D((title+"_MConLEPS_MAXMInv").c_str(), (title+"_MConLEPS_MAXMInv;M_{Con}^{LEPS} [GeV];MAX MInv [GeV]").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConLEPS_MAXMInv);
    TH2D* hist_MConLEPS_RISR = new TH2D((title+"_MConLEPS_RISR").c_str(), (title+"_MConLEPS_RISR;M_{Con}^{LEPS} [GeV];RISR").c_str(), g_NX, 0., 400., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MConLEPS_RISR);
    TH2D* hist_MConLEPS_MInv = new TH2D((title+"_MConLEPS_MInv").c_str(), (title+"_MConLEPS_MInv;M_{Con}^{LEPS} [GeV];MInv").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConLEPS_MInv);
    TH2D* hist_MT2_MT2MAX = new TH2D((title+"_MT2_MT2MAX").c_str(), (title+"_MT2_MT2MAX;MT2;MAX MT2").c_str(), g_NX, 0., 250., g_NX, 0., 250.);
    //hists2.push_back(hist_MT2_MT2MAX);
    TH2D* hist_MT2_MConLEPSMET = new TH2D((title+"_MT2_MConLEPSMET").c_str(), (title+"_MT2_MConLEPSMET;MT2;M_{Con}^{LEPS MET} [GeV]").c_str(), g_NX, 0., 250., g_NX, 0., 400.);
    //hists2.push_back(hist_MT2_MConLEPSMET);
    TH2D* hist_MT2_CosThetaRazor = new TH2D((title+"_MT2_CosThetaRazor").c_str(), (title+"_MT2_CosThetaRazor;MT2;cos#theta^{Razor}").c_str(), g_NX, 0., 250., g_NX, 0., 1.);
    //hists2.push_back(hist_MT2_CosThetaRazor);
    TH2D* hist_MT2_MAXCosThetaRazor = new TH2D((title+"_MT2_MAXCosThetaRazor").c_str(), (title+"_MT2_MAXCosThetaRazor;MT2;cos#theta^{MAXRazor}").c_str(), g_NX, 0., 250., g_NX, 0., 1.);
    //hists2.push_back(hist_MT2_MAXCosThetaRazor);
    TH2D* hist_MT2_MAXMInv = new TH2D((title+"_MT2_MAXMInv").c_str(), (title+"_MT2_MAXMInv;MT2;MAX MInv").c_str(), g_NX, 0., 250., g_NX, 0., 1000.);
    //hists2.push_back(hist_MT2_MAXMInv);
    TH2D* hist_MT2_RISR = new TH2D((title+"_MT2_RISR").c_str(), (title+"_MT2_RISR;MT2;RISR").c_str(), g_NX, 0., 250., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MT2_RISR);
    TH2D* hist_MT2_MInv = new TH2D((title+"_MT2_MInv").c_str(), (title+"_MT2_MInv;MT2;MInv").c_str(), g_NX, 0., 250., g_NX, 0., 1000.);
    //hists2.push_back(hist_MT2_MInv);
    TH2D* hist_MT2MAX_MConLEPSMET = new TH2D((title+"_MT2MAX_MConLEPSMET").c_str(), (title+"_MT2MAX_MConLEPSMET;MAX MT2;M_{Con}^{LEPS MET} [GeV]").c_str(), g_NX, 0., 250., g_NX, 0., 400.);
    //hists2.push_back(hist_MT2MAX_MConLEPSMET);
    TH2D* hist_MT2MAX_CosThetaRazor = new TH2D((title+"_MT2MAX_CosThetaRazor").c_str(), (title+"_MT2MAX_CosThetaRazor;MAX MT2;cos#theta^{Razor}").c_str(), g_NX, 0., 250., g_NX, 0., 1.);
    //hists2.push_back(hist_MT2MAX_CosThetaRazor);
    TH2D* hist_MT2MAX_MAXCosThetaRazor = new TH2D((title+"_MT2MAX_MAXCosThetaRazor").c_str(), (title+"_MT2MAX_MAXCosThetaRazor;MAX MT2;cos#theta^{MAXRazor}").c_str(), g_NX, 0., 250., g_NX, 0., 1.);
    //hists2.push_back(hist_MT2MAX_MAXCosThetaRazor);
    TH2D* hist_MT2MAX_MAXMInv = new TH2D((title+"_MT2MAX_MAXMInv").c_str(), (title+"_MT2MAX_MAXMInv;MAX MT2;MAX MInv").c_str(), g_NX, 0., 250., g_NX, 0., 1000.);
    //hists2.push_back(hist_MT2MAX_MAXMInv);
    TH2D* hist_MT2MAX_RISR = new TH2D((title+"_MT2MAX_RISR").c_str(), (title+"_MT2MAX_RISR;MAX MT2;RISR").c_str(), g_NX, 0., 250., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MT2MAX_RISR);
    TH2D* hist_MT2MAX_MInv = new TH2D((title+"_MT2MAX_MInv").c_str(), (title+"_MT2MAX_MInv;MAX MT2;MInv").c_str(), g_NX, 0., 250., g_NX, 0., 1000.);
    //hists2.push_back(hist_MT2MAX_MInv);
    TH2D* hist_MConLEPSMET_CosThetaRazor = new TH2D((title+"_MConLEPSMET_CosThetaRazor").c_str(), (title+"_MConLEPSMET_CosThetaRazor;M_{Con}^{LEPS MET};cos#theta^{Razor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConLEPSMET_CosThetaRazor);
    TH2D* hist_MConLEPSMET_MAXCosThetaRazor = new TH2D((title+"_MConLEPSMET_MAXCosThetaRazor").c_str(), (title+"_MConLEPSMET_MAXCosThetaRazor;M_{Con}^{LEPS MET};cos#theta^{MAXRazor}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    //hists2.push_back(hist_MConLEPSMET_MAXCosThetaRazor);
    TH2D* hist_MConLEPSMET_MAXMInv = new TH2D((title+"_MConLEPSMET_MAXMInv").c_str(), (title+"_MConLEPSMET_MAXMInv;M_{Con}^{LEPS MET};MAX MInv").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConLEPSMET_MAXMInv);
    TH2D* hist_MConLEPSMET_RISR = new TH2D((title+"_MConLEPSMET_RISR").c_str(), (title+"_MConLEPSMET_RISR;M_{Con}^{LEPS MET};RISR").c_str(), g_NX, 0., 400., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MConLEPSMET_RISR);
    TH2D* hist_MConLEPSMET_MInv = new TH2D((title+"_MConLEPSMET_MInv").c_str(), (title+"_MConLEPSMET_MInv;M_{Con}^{LEPS MET};MInv").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    //hists2.push_back(hist_MConLEPSMET_MInv);
    TH2D* hist_CosThetaRazor_MAXCosThetaRazor = new TH2D((title+"_CosThetaRazor_MAXCosThetaRazor").c_str(), (title+"_CosThetaRazor_MAXCosThetaRazor;cos#theta^{Razor};cos#theta^{MAXRazor}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    //hists2.push_back(hist_CosThetaRazor_MAXCosThetaRazor);
    TH2D* hist_CosThetaRazor_MAXMInv = new TH2D((title+"_CosThetaRazor_MAXMInv").c_str(), (title+"_CosThetaRazor_MAXMInv;cos#theta^{Razor};MAX MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_CosThetaRazor_MAXMInv);
    TH2D* hist_CosThetaRazor_RISR = new TH2D((title+"_CosThetaRazor_RISR").c_str(), (title+"_CosThetaRazor_RISR;cos#theta^{Razor};RISR").c_str(), g_NX, 0., 1., g_NX, 0.5, 1.);
    //hists2.push_back(hist_CosThetaRazor_RISR);
    TH2D* hist_CosThetaRazor_MInv = new TH2D((title+"_CosThetaRazor_MInv").c_str(), (title+"_CosThetaRazor_MInv;cos#theta^{Razor};MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_CosThetaRazor_MInv);
    TH2D* hist_MAXCosThetaRazor_MAXMInv = new TH2D((title+"_MAXCosThetaRazor_MAXMInv").c_str(), (title+"_MAXCosThetaRazor_MAXMInv;cos#theta^{MAXRazor};MAX MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_MAXCosThetaRazor_MAXMInv);
    TH2D* hist_MAXCosThetaRazor_RISR = new TH2D((title+"_MAXCosThetaRazor_RISR").c_str(), (title+"_MAXCosThetaRazor_RISR;cos#theta^{MAXRazor};RISR").c_str(), g_NX, 0., 1., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MAXCosThetaRazor_RISR);
    TH2D* hist_MAXCosThetaRazor_MInv = new TH2D((title+"_MAXCosThetaRazor_MInv").c_str(), (title+"_MAXCosThetaRazor_MInv;cos#theta^{MAXRazor};MInv").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    //hists2.push_back(hist_MAXCosThetaRazor_MInv);
    TH2D* hist_MAXMInv_RISR = new TH2D((title+"_MAXMInv_RISR").c_str(), (title+"_MAXMInv_RISR;MAX MInv;RISR").c_str(), g_NX, 0., 1000., g_NX, 0.5, 1.);
    //hists2.push_back(hist_MAXMInv_RISR);
    TH2D* hist_MAXMInv_MInv = new TH2D((title+"_MAXMInv_MInv").c_str(), (title+"_MAXMInv_MInv;MAX Inv;MInv").c_str(), g_NX, 0., 1000., g_NX, 0., 1000.);
    //hists2.push_back(hist_MAXMInv_MInv);
    TH2D* hist_MAXCandML_MAXCosDecayAngle = new TH2D((title+"_MAXCandML_MAXCosDecayAngle").c_str(), (title+"_MAXCandML_MAXCosDecayAngle;MAX Cand ML [GeV];MAX cos#theta").c_str(), g_NX, CandMinMass, CandMaxMass, g_NX, 0., 1.);
    //hists2.push_back(hist_MAXCandML_MAXCosDecayAngle);

    TH2D* hist_RISRLEP_RISR = new TH2D((title+"_RISRLEP_RISR").c_str(), (title+"_RISRLEP_RISR;R_{ISR}^{LEP};R_{ISR}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    hists2.push_back(hist_RISRLEP_RISR);
    TH2D* hist_MQLEP_MSLEP = new TH2D((title+"_MQLEP_MSLEP").c_str(), (title+"_MQLEP_MSLEP;MQ^{LEP};MS^{LEP}").c_str(), g_NX, 0., 400., g_NX, 0., 1000.);
    hists2.push_back(hist_MQLEP_MSLEP);
    TH2D* hist_gammaLEP_MQLEP = new TH2D((title+"_gammaLEP_MQLEP").c_str(), (title+"_gammaLEP_MQLEP;#gamma^{LEP};MQ^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    hists2.push_back(hist_gammaLEP_MQLEP);
    TH2D* hist_gammaLEP_MSLEP = new TH2D((title+"_gammaLEP_MSLEP").c_str(), (title+"_gammaLEP_MSLEP;#gamma^{LEP};MS^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    hists2.push_back(hist_gammaLEP_MSLEP);
    TH2D* hist_MPa_CosDecayAnglePaLEP = new TH2D((title+"_MPa_CosDecayAnglePaLEP").c_str(), (title+"_MPa_CosDecayAnglePaLEP;M_{Pa}^{LEP};cos#theta_{Pa}^{LEP}").c_str(), g_NX, 0., 300., g_NX, -1., 1.);
    hists2.push_back(hist_MPa_CosDecayAnglePaLEP);
    TH2D* hist_MSLEP_CosDecayAngleSLEP = new TH2D((title+"_MSLEP_CosDecayAngleSLEP").c_str(), (title+"_MSLEP_CosDecayAngleSLEP;MS^{LEP};cos#theta_{S}^{LEP}").c_str(), g_NX, 0., 1000., g_NX, -1., 1.);
    hists2.push_back(hist_MSLEP_CosDecayAngleSLEP);

    TH2D* hist_MET_RISRLEP = new TH2D((title+"_MET_RISRLEP").c_str(), (title+"_MET_RISRLEP;MET;R_{ISR}^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 1.0);
    hists2.push_back(hist_MET_RISRLEP);
    TH2D* hist_MET_PTISRLEP = new TH2D((title+"_MET_PTISRLEP").c_str(), (title+"_MET_PTISRLEP;MET;p_{T}^{ISR LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 800.);
    hists2.push_back(hist_MET_PTISRLEP);
    TH2D* hist_RISRLEP_PTISRLEP = new TH2D((title+"_RISRLEP_PTISRLEP").c_str(), (title+"_RISRLEP_PTISRLEP;R_{ISR}^{LEP};p_{T}^{ISR LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 800.);
    hists2.push_back(hist_RISRLEP_PTISRLEP);
    TH2D* hist_PTISR_PTISRLEP = new TH2D((title+"_PTISR_PTISRLEP").c_str(), (title+"_PTISR_PTISRLEP;p_{T}^{ISR};p_{T}^{ISR LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 800.);
    hists2.push_back(hist_PTISR_PTISRLEP);
    TH2D* hist_MET_MS = new TH2D((title+"_MET_MS").c_str(), (title+"_MET_MS;MET;MS^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 800.);
    //hists2.push_back(hist_MET_MS);
    TH2D* hist_MET_MPa = new TH2D((title+"_MET_MPa").c_str(), (title+"_MET_MPa;MET;MPa^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 500.);
    //hists2.push_back(hist_MET_MPa);
    TH2D* hist_MET_MPb = new TH2D((title+"_MET_MPb").c_str(), (title+"_MET_MPb;MET;MPb^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 400.);
    //hists2.push_back(hist_MET_MPb);
    TH2D* hist_MET_MVa = new TH2D((title+"_MET_MVa").c_str(), (title+"_MET_MVa;MET;MVa^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 250.);
    //hists2.push_back(hist_MET_MVa);
    TH2D* hist_MET_MVb = new TH2D((title+"_MET_MVb").c_str(), (title+"_MET_MVb;MET;MVb^{LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 250.);
    //hists2.push_back(hist_MET_MVb);
    TH2D* hist_MET_PTSCM = new TH2D((title+"_MET_PTSCM").c_str(), (title+"_MET_PTSCM;MET;P_{T,S}^{CM LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 800.);
    //hists2.push_back(hist_MET_PTSCM);
    TH2D* hist_MET_MSS0 = new TH2D((title+"_MET_MSS0").c_str(), (title+"_MET_MSS0;MET;M_{S}^{S0 LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 1000.);
    //hists2.push_back(hist_MET_MSS0);
    TH2D* hist_MET_MVS0 = new TH2D((title+"_MET_MVS0").c_str(), (title+"_MET_MVS0;MET;M_{V}^{S0 LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 800.);
    //hists2.push_back(hist_MET_MVS0);
    TH2D* hist_MET_MQS0 = new TH2D((title+"_MET_MQS0").c_str(), (title+"_MET_MQS0;MET;M_{Q}^{S0 LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 400.);
    //hists2.push_back(hist_MET_MQS0);
    TH2D* hist_MET_gammaS0 = new TH2D((title+"_MET_gammaS0").c_str(), (title+"_MET_gammaS0;MET;#gamma^{S0 LEP}").c_str(), g_NX, 0., 700., g_NX, 0., 1.);
    //hists2.push_back(hist_MET_gammaS0);
    TH2D* hist_MET_MPTilde = new TH2D((title+"_MET_MPTilde").c_str(), (title+"_MET_MPTilde;MET;#tilde{M}_{P}").c_str(), g_NX, 0., 700., g_NX, 0., 400.);
    //hists2.push_back(hist_MET_MPTilde);
    TH2D* hist_MET_MSTilde = new TH2D((title+"_MET_MSTilde").c_str(), (title+"_MET_MSTilde;MET;#tilde{M}_{S}").c_str(), g_NX, 0., 700., g_NX, 0., 800.);
    //hists2.push_back(hist_MET_MSTilde);
    TH2D* hist_MET_gammaTilde = new TH2D((title+"_MET_gammaTilde").c_str(), (title+"_MET_gammaTilde;MET;#tilde{#gamma}").c_str(), g_NX, 0., 700., g_NX, 0., 1.);
    //hists2.push_back(hist_MET_gammaTilde);
    TH2D* hist_MS_MPa = new TH2D((title+"_MS_MPa").c_str(), (title+"_MS_MPa;M_{S}^{LEP};M_{Pa}^{LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 500.);
    hists2.push_back(hist_MS_MPa);
    TH2D* hist_MS_MPb = new TH2D((title+"_MS_MPb").c_str(), (title+"_MS_MPb;M_{S}^{LEP};M_{Pb}^{LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 400.);
    hists2.push_back(hist_MS_MPb);
    TH2D* hist_MS_CosDecayAnglePaLEP = new TH2D((title+"_MS_CosDecayAnglePaLEP").c_str(), (title+"_MS_CosDecayAnglePaLEP;M_{S}^{LEP};cos#theta_{Pa}^{LEP}").c_str(), g_NX, 0., 800., g_NX, -1, 1.);
    hists2.push_back(hist_MS_CosDecayAnglePaLEP);
    TH2D* hist_MS_MVa = new TH2D((title+"_MS_MVa").c_str(), (title+"_MS_MVa;M_{S}^{LEP};M_{Va}^{LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 250.);
    hists2.push_back(hist_MS_MVa);
    TH2D* hist_MS_MVb = new TH2D((title+"_MS_MVb").c_str(), (title+"_MS_MVb;M_{S}^{LEP};M_{Vb}^{LEP}").c_str(), g_NX, 0., 800., g_NX, 0., 250.);
    //hists2.push_back(hist_MS_MVb);
    TH2D* hist_MS_PTSCM = new TH2D((title+"_MS_PTSCM").c_str(), (title+"_MS_PTSCM;M_{S}^{LEP};p_{T,S}^{CM}").c_str(), g_NX, 0., 800., g_NX, 0., 800.);
    hists2.push_back(hist_MS_PTSCM);
    TH2D* hist_MPa_MPb = new TH2D((title+"_MPa_MPb").c_str(), (title+"_MPa_MPb;M_{Pa}^{LEP};M_{Pb}^{LEP}").c_str(), g_NX, 0., 500., g_NX, 0., 400.);
    hists2.push_back(hist_MPa_MPb);
    TH2D* hist_MPa_MVa = new TH2D((title+"_MPa_MVa").c_str(), (title+"_MPa_MVa;M_{Pa}^{LEP};M_{Va}^{LEP}").c_str(), g_NX, 0., 500., g_NX, 0., 250.);
    hists2.push_back(hist_MPa_MVa);
    TH2D* hist_CosDecayAnglePaLEP_CosDecayAnglePbLEP = new TH2D((title+"_CosDecayAnglePaLEP_CosDecayAnglePbLEP").c_str(), (title+"_CosDecayAnglePaLEP_CosDecayAnglePbLEP;cos#theta_{Pa}^{LEP};cos#theta_{Pb}^{LEP}").c_str(), g_NX, -1., 1., g_NX, -1., 1.);
    hists2.push_back(hist_CosDecayAnglePaLEP_CosDecayAnglePbLEP);
    TH2D* hist_CosDecayAnglePaLEP_CosDecayAngleSLEP = new TH2D((title+"_CosDecayAnglePaLEP_CosDecayAngleSLEP").c_str(), (title+"_CosDecayAnglePaLEP_CosDecayAngleSLEP;cos#theta_{Pa}^{LEP};cos#theta_{S}^{LEP}").c_str(), g_NX, -1., 1., g_NX, -1., 1.);
    hists2.push_back(hist_CosDecayAnglePaLEP_CosDecayAngleSLEP);
    TH2D* hist_MSS0_MVS0 = new TH2D((title+"_MSS0_MVS0").c_str(), (title+"_MSS0_MVS0;M_{S}^{S0};M_{V}^{S0}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MSS0_MVS0);
    TH2D* hist_MSS0_MQS0 = new TH2D((title+"_MSS0_MQS0").c_str(), (title+"_MSS0_MQS0;M_{S}^{S0};M_{Q}^{S0}").c_str(), g_NX, 0., 1000., g_NX, 0., 400.);
    hists2.push_back(hist_MSS0_MQS0);
    TH2D* hist_MSS0_gammaS0 = new TH2D((title+"_MSS0_gammaS0").c_str(), (title+"_MSS0_gammaS0;M_{S}^{S0};#gamma^{S0}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    hists2.push_back(hist_MSS0_gammaS0);
    TH2D* hist_MSS0_MSTilde = new TH2D((title+"_MSS0_MSTilde").c_str(), (title+"_MSS0_MSTilde;M_{S}^{S0};#tilde{M}_{S}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MSS0_MSTilde);
    TH2D* hist_MVS0_MQS0 = new TH2D((title+"_MVS0_MQS0").c_str(), (title+"_MVS0_MQS0;M_{V}^{S0};M_{Q}^{S0}").c_str(), g_NX, 0., 800., g_NX, 0., 400.);
    hists2.push_back(hist_MVS0_MQS0);
    TH2D* hist_MVS0_gammaS0 = new TH2D((title+"_MVS0_gammaS0").c_str(), (title+"_MVS0_gammaS0;M_{V}^{S0};#gamma^{S0}").c_str(), g_NX, 0., 800., g_NX, 0., 1.);
    hists2.push_back(hist_MVS0_gammaS0);
    TH2D* hist_MQS0_gammaS0 = new TH2D((title+"_MQS0_gammaS0").c_str(), (title+"_MQS0_gammaS0;M_{Q}^{S0};#gamma^{S0}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    hists2.push_back(hist_MQS0_gammaS0);
    TH2D* hist_MQS0_MPTilde = new TH2D((title+"_MQS0_MPTilde").c_str(), (title+"_MQS0_MPTilde;M_{Q}^{S0};#tilde{M}_{P}").c_str(), g_NX, 0., 400., g_NX, 0., 400.);
    hists2.push_back(hist_MQS0_MPTilde);
    TH2D* hist_gammaS0_gammaTilde = new TH2D((title+"_gammaS0_gammaTilde").c_str(), (title+"_gammaS0_gammaTilde;#gamma^{S0};#tilde{#gamma}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    hists2.push_back(hist_gammaS0_gammaTilde);
    TH2D* hist_MPTilde_MSTilde = new TH2D((title+"_MPTilde_MSTilde").c_str(), (title+"_MPTilde_MSTilde;#tilde{M}_{P};#tilde{M}_{S}").c_str(), g_NX, 0., 400., g_NX, 0., 800.);
    hists2.push_back(hist_MPTilde_MSTilde);
    TH2D* hist_MPTilde_gammaTilde = new TH2D((title+"_MPTilde_gammaTilde").c_str(), (title+"_MPTilde_gammaTilde;#tilde{M}_{P};#tilde{#gamma}").c_str(), g_NX, 0., 400., g_NX, 0., 1.);
    hists2.push_back(hist_MPTilde_gammaTilde);
    TH2D* hist_MSTilde_gammaTilde = new TH2D((title+"_MSTilde_gammaTilde").c_str(), (title+"_MSTilde_gammaTilde;#tilde{M}_{P};#tilde{#gamma}").c_str(), g_NX, 0., 800., g_NX, 0., 1.);
    hists2.push_back(hist_MSTilde_gammaTilde);

    TH2D* hist_RISRLEP_MS = new TH2D((title+"_RISRLEP_MS").c_str(), (title+"_RISRLEP_MS;R_{ISR}^{LEP};MS^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 800.);
    hists2.push_back(hist_RISRLEP_MS);
    TH2D* hist_RISRLEP_MPa = new TH2D((title+"_RISRLEP_MPa").c_str(), (title+"_RISRLEP_MPa;R_{ISR}^{LEP};MPa^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 500.);
    hists2.push_back(hist_RISRLEP_MPa);
    TH2D* hist_RISRLEP_MPb = new TH2D((title+"_RISRLEP_MPb").c_str(), (title+"_RISRLEP_MPb;R_{ISR}^{LEP};MPb^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    hists2.push_back(hist_RISRLEP_MPb);
    TH2D* hist_RISRLEP_MVa = new TH2D((title+"_RISRLEP_MVa").c_str(), (title+"_RISRLEP_MVa;R_{ISR}^{LEP};MVa^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 250.);
    hists2.push_back(hist_RISRLEP_MVa);
    TH2D* hist_RISRLEP_MVb = new TH2D((title+"_RISRLEP_MVb").c_str(), (title+"_RISRLEP_MVb;R_{ISR}^{LEP};MVb^{LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 250.);
    //hists2.push_back(hist_RISRLEP_MVb);
    TH2D* hist_RISRLEP_PTSCM = new TH2D((title+"_RISRLEP_PTSCM").c_str(), (title+"_RISRLEP_PTSCM;R_{ISR}^{LEP};P_{T,S}^{CM LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 800.);
    hists2.push_back(hist_RISRLEP_PTSCM);
    TH2D* hist_RISRLEP_MSS0 = new TH2D((title+"_RISRLEP_MSS0").c_str(), (title+"_RISRLEP_MSS0;R_{ISR}^{LEP};M_{S}^{S0 LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 1000.);
    hists2.push_back(hist_RISRLEP_MSS0);
    TH2D* hist_RISRLEP_MVS0 = new TH2D((title+"_RISRLEP_MVS0").c_str(), (title+"_RISRLEP_MVS0;R_{ISR}^{LEP};M_{V}^{S0 LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 800.);
    hists2.push_back(hist_RISRLEP_MVS0);
    TH2D* hist_RISRLEP_MQS0 = new TH2D((title+"_RISRLEP_MQS0").c_str(), (title+"_RISRLEP_MQS0;R_{ISR}^{LEP};M_{Q}^{S0 LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    hists2.push_back(hist_RISRLEP_MQS0);
    TH2D* hist_RISRLEP_gammaS0 = new TH2D((title+"_RISRLEP_gammaS0").c_str(), (title+"_RISRLEP_gammaS0;R_{ISR}^{LEP};#gamma^{S0 LEP}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    hists2.push_back(hist_RISRLEP_gammaS0);
    TH2D* hist_RISRLEP_MPTilde = new TH2D((title+"_RISRLEP_MPTilde").c_str(), (title+"_RISRLEP_MPTilde;R_{ISR}^{LEP};#tilde{M}_{P}").c_str(), g_NX, 0., 1., g_NX, 0., 400.);
    hists2.push_back(hist_RISRLEP_MPTilde);
    TH2D* hist_RISRLEP_MSTilde = new TH2D((title+"_RISRLEP_MSTilde").c_str(), (title+"_RISRLEP_MSTilde;R_{ISR}^{LEP};#tilde{M}_{S}").c_str(), g_NX, 0., 1., g_NX, 0., 800.);
    hists2.push_back(hist_RISRLEP_MSTilde);
    TH2D* hist_RISRLEP_gammaTilde = new TH2D((title+"_RISRLEP_gammaTilde").c_str(), (title+"_RISRLEP_gammaTilde;R_{ISR}^{LEP};#tilde{#gamma}").c_str(), g_NX, 0., 1., g_NX, 0., 1.);
    hists2.push_back(hist_RISRLEP_gammaTilde);

    TH2D* hist_MINVLEP_RISRLEP = new TH2D((title+"_MINVLEP_RISRLEP").c_str(), (title+"_MINVLEP_PTISRLEP;M_{INV}^{LEP};R_{ISR}^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 0.5);
    hists2.push_back(hist_MINVLEP_RISRLEP);
    TH2D* hist_MINVLEP_MET = new TH2D((title+"_MINVLEP_MET").c_str(), (title+"_MINVLEP_MET;M_{INV}^{LEP};MET").c_str(), g_NX, 0., 1000., g_NX, 0., 700.);
    hists2.push_back(hist_MINVLEP_MET);
    TH2D* hist_MINVLEP_PTISRLEP = new TH2D((title+"_MINVLEP_PTISRLEP").c_str(), (title+"_MINVLEP_PTISRLEP;M_{INV}^{LEP};p_{T}^{ISR LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MINVLEP_PTISRLEP);
    TH2D* hist_MINVLEP_MS = new TH2D((title+"_MINVLEP_MS").c_str(), (title+"_MINVLEP_MS;M_{INV}^{LEP};MS^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MINVLEP_MS);
    TH2D* hist_MINVLEP_MPa = new TH2D((title+"_MINVLEP_MPa").c_str(), (title+"_MINVLEP_MPa;M_{INV}^{LEP};MPa^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 500.);
    hists2.push_back(hist_MINVLEP_MPa);
    TH2D* hist_MINVLEP_MPb = new TH2D((title+"_MINVLEP_MPb").c_str(), (title+"_MINVLEP_MPb;M_{INV}^{LEP};MPb^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 400.);
    hists2.push_back(hist_MINVLEP_MPb);
    TH2D* hist_MINVLEP_MVa = new TH2D((title+"_MINVLEP_MVa").c_str(), (title+"_MINVLEP_MVa;M_{INV}^{LEP};MVa^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 250.);
    hists2.push_back(hist_MINVLEP_MVa);
    TH2D* hist_MINVLEP_MVb = new TH2D((title+"_MINVLEP_MVb").c_str(), (title+"_MINVLEP_MVb;M_{INV}^{LEP};MVb^{LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 250.);
    //hists2.push_back(hist_MINVLEP_MVb);
    TH2D* hist_MINVLEP_PTSCM = new TH2D((title+"_MINVLEP_PTSCM").c_str(), (title+"_MINVLEP_PTSCM;M_{INV}^{LEP};P_{T,S}^{CM LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MINVLEP_PTSCM);
    TH2D* hist_MINVLEP_MSS0 = new TH2D((title+"_MINVLEP_MSS0").c_str(), (title+"_MINVLEP_MSS0;M_{INV}^{LEP};M_{S}^{S0 LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 1000.);
    hists2.push_back(hist_MINVLEP_MSS0);
    TH2D* hist_MINVLEP_MVS0 = new TH2D((title+"_MINVLEP_MVS0").c_str(), (title+"_MINVLEP_MVS0;M_{INV}^{LEP};M_{V}^{S0 LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MINVLEP_MVS0);
    TH2D* hist_MINVLEP_MQS0 = new TH2D((title+"_MINVLEP_MQS0").c_str(), (title+"_MINVLEP_MQS0;M_{INV}^{LEP};M_{Q}^{S0 LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 400.);
    hists2.push_back(hist_MINVLEP_MQS0);
    TH2D* hist_MINVLEP_gammaS0 = new TH2D((title+"_MINVLEP_gammaS0").c_str(), (title+"_MINVLEP_gammaS0;M_{INV}^{LEP};#gamma^{S0 LEP}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    hists2.push_back(hist_MINVLEP_gammaS0);
    TH2D* hist_MINVLEP_MPTilde = new TH2D((title+"_MINVLEP_MPTilde").c_str(), (title+"_MINVLEP_MPTilde;M_{INV}^{LEP};#tilde{M}_{P}").c_str(), g_NX, 0., 1000., g_NX, 0., 400.);
    hists2.push_back(hist_MINVLEP_MPTilde);
    TH2D* hist_MINVLEP_MSTilde = new TH2D((title+"_MINVLEP_MSTilde").c_str(), (title+"_MINVLEP_MSTilde;M_{INV}^{LEP};#tilde{M}_{S}").c_str(), g_NX, 0., 1000., g_NX, 0., 800.);
    hists2.push_back(hist_MINVLEP_MSTilde);
    TH2D* hist_MINVLEP_gammaTilde = new TH2D((title+"_MINVLEP_gammaTilde").c_str(), (title+"_MINVLEP_gammaTilde;M_{INV}^{LEP};#tilde{#gamma}").c_str(), g_NX, 0., 1000., g_NX, 0., 1.);
    hists2.push_back(hist_MINVLEP_gammaTilde);
    TH2D* hist_RPLEP_MSLEP = new TH2D((title+"_RPLEP_MSLEP").c_str(), (title+"_RPLEP_MSLEP;M_{Pa}^{LEP}/M_{Pb}^{LEP};M_{S}^{LEP}").c_str(), g_NX, 0., 2., g_NX, 0., 500.);
    hists2.push_back(hist_RPLEP_MSLEP);
    TH2D* hist_RPLEP_MQLEP = new TH2D((title+"_RPLEP_MQLEP").c_str(), (title+"_RPLEP_MQLEP;M_{Pa}^{LEP}/M_{Pb}^{LEP};M_{Q}^{LEP}").c_str(), g_NX, 0., 2., g_NX, 0., 500.);
    hists2.push_back(hist_RPLEP_MQLEP);

    TH2D* hist_MVa_CosDecayAngleVaLEP = new TH2D((title+"_MVa_CosDecayAngleVaLEP").c_str(), (title+"_MVa_CosDecayAngleVaLEP;M_{Va}^{LEP};cos#theta_{Va}^{LEP}").c_str(), g_NX, 0., 300., g_NX, -1., 1.);
    hists2.push_back(hist_MVa_CosDecayAngleVaLEP);
    TH2D* hist_MS_CosDecayAngleVaLEP = new TH2D((title+"_MS_CosDecayAngleVaLEP").c_str(), (title+"_MS_CosDecayAngleVaLEP;M_{S}^{LEP};cos#theta_{Va}^{LEP}").c_str(), g_NX, 0., 800., g_NX, -1, 1.);
    hists2.push_back(hist_MS_CosDecayAngleVaLEP);
    TH2D* hist_CosDecayAngleVaLEP_CosDecayAnglePbLEP = new TH2D((title+"_CosDecayAngleVaLEP_CosDecayAnglePbLEP").c_str(), (title+"_CosDecayAngleVaLEP_CosDecayAnglePbLEP;cos#theta_{Va}^{LEP};cos#theta_{Pb}^{LEP}").c_str(), g_NX, -1., 1., g_NX, -1., 1.);
    hists2.push_back(hist_CosDecayAngleVaLEP_CosDecayAnglePbLEP);
    TH2D* hist_CosDecayAngleVaLEP_CosDecayAngleSLEP = new TH2D((title+"_CosDecayAngleVaLEP_CosDecayAngleSLEP").c_str(), (title+"_CosDecayAngleVaLEP_CosDecayAngleSLEP;cos#theta_{Va}^{LEP};cos#theta_{S}^{LEP}").c_str(), g_NX, -1., 1., g_NX, -1., 1.);
    hists2.push_back(hist_CosDecayAngleVaLEP_CosDecayAngleSLEP);

    TH2D* hist_MSLEP_MCosDecayAngleP = new TH2D((title+"_MSLEP_MCosDecayAngleP").c_str(), (title+"_MSLEP_MCosDecayAngleP;M_{S}^{LEP};cos#theta_{Pa}^{LEP}-cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 600., g_NX, -2., 2.);
    hists2.push_back(hist_MSLEP_MCosDecayAngleP);
    TH2D* hist_MQLEP_MCosDecayAngleP = new TH2D((title+"_MQLEP_MCosDecayAngleP").c_str(), (title+"_MQLEP_MCosDecayAngleP;M_{Q}^{LEP};cos#theta_{Pa}^{LEP}-cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 600., g_NX, -2., 2.);
    hists2.push_back(hist_MQLEP_MCosDecayAngleP);

    TH2D* hist_CosDecayAnglePaLEP_MPb = new TH2D((title+"_CosDecayAnglePaLEP_MPb").c_str(), (title+"_CosDecayAnglePaLEP_MPb;cos#theta^{Pa};M_{Pb}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 300.);
    hists2.push_back(hist_CosDecayAnglePaLEP_MPb);
    TH2D* hist_CosDecayAnglePaLEP_MVa = new TH2D((title+"_CosDecayAnglePaLEP_MVa").c_str(), (title+"_CosDecayAnglePaLEP_MVa;cos#theta^{Pa};M_{Va}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 200.);
    hists2.push_back(hist_CosDecayAnglePaLEP_MVa);
    TH2D* hist_CosDecayAnglePbLEP_MPa = new TH2D((title+"_CosDecayAnglePbLEP_MPa").c_str(), (title+"_CosDecayAnglePbLEP_MPa;cos#theta^{Pb};M_{Pa}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 300.);
    hists2.push_back(hist_CosDecayAnglePbLEP_MPa);
    TH2D* hist_CosDecayAnglePbLEP_MPb = new TH2D((title+"_CosDecayAnglePbLEP_MPb").c_str(), (title+"_CosDecayAnglePbLEP_MPb;cos#theta^{Pb};M_{Pb}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 300.);
    hists2.push_back(hist_CosDecayAnglePbLEP_MPb);
    TH2D* hist_CosDecayAnglePbLEP_MVa = new TH2D((title+"_CosDecayAnglePbLEP_MVa").c_str(), (title+"_CosDecayAnglePbLEP_MVa;cos#theta^{Pb};M_{Va}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 200.);
    hists2.push_back(hist_CosDecayAnglePbLEP_MVa);
    TH2D* hist_CosDecayAnglePbLEP_MS = new TH2D((title+"_CosDecayAnglePbLEP_MS").c_str(), (title+"_CosDecayAnglePbLEP_MS;cos#theta^{Pb};M_{S}^{LEP}").c_str(), g_NX, -1., 1., g_NX, 0., 1000.);
    hists2.push_back(hist_CosDecayAnglePbLEP_MS);
    TH2D* hist_MPb_MVa = new TH2D((title+"_MPb_MVa").c_str(), (title+"_MPb_MVa;M_{Pb}^{LEP};M_{Va}^{LEP}").c_str(), g_NX, 0., 300., g_NX, 0., 200.);
    hists2.push_back(hist_MPb_MVa); 
    TH2D* hist_MPa_MCosDecayAngleM = new TH2D((title+"_MPa_MCosDecayAngleM").c_str(), (title+"_MPa_MCosDecayAngleM;M_{Pa}^{LEP};cos#theta_{Pa}^{LEP}*cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 300., g_NX, -1., 1.);
    hists2.push_back(hist_MPa_MCosDecayAngleM);
    TH2D* hist_MPb_MCosDecayAngleM = new TH2D((title+"_MPb_MCosDecayAngleM").c_str(), (title+"_MPb_MCosDecayAngleM;M_{Pb}^{LEP};cos#theta_{Pb}^{LEP}*cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 300., g_NX, -1., 1.);
    hists2.push_back(hist_MPb_MCosDecayAngleM);
    TH2D* hist_MPa_SCosDecayAngleP = new TH2D((title+"_MPa_SCosDecayAngleP").c_str(), (title+"_MPa_SCosDecayAngleP;M_{Pa}^{LEP};cos#theta_{Pa}^{LEP}+cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 300., g_NX, 0., 2.);
    hists2.push_back(hist_MPa_SCosDecayAngleP);
    TH2D* hist_MPb_SCosDecayAngleP = new TH2D((title+"_MPb_SCosDecayAngleP").c_str(), (title+"_MPb_SCosDecayAngleP;M_{Pb}^{LEP};cos#theta_{Pb}^{LEP}*cos#theta_{Pb}^{LEP}").c_str(), g_NX, 0., 300., g_NX, 0., 2.);
    hists2.push_back(hist_MPb_SCosDecayAngleP);

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
        
        gErrorIgnoreLevel = kFatal;
        int Nentry = base->fChain->GetEntries(); 
        gErrorIgnoreLevel = 0;
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
          vect_str_cutflow_labels[CF_bin] = "NTUPLES";

          // Get Physics Objects
          int Njet     = base->Njet;
          int Nlep     = base->Nlep;
          int NjetS    = base->Njet_S;
          int NbjetS   = base->Nbjet_S;
          int NjetISR  = base->Njet_ISR;
          int Njet_a   = base->Njet_a;
          int Njet_b   = base->Njet_b;
          int Nlep_a   = base->Nlep_a;
          int Nlep_b   = base->Nlep_b;
          int NbjetISR = base->Nbjet_ISR;

          // Apply PreSelection
          //if(base->Njet == 0) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "Njet > 0";

          if(Nlep != 3) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "= 3L";
          
          //if(do_FilterDilepton)
          //  if(SF.DileptonEvent(base))
          //    continue;
          
          // get variables from root files using base class
          double MET = base->MET;
          double Mperp = base->Mperp;
          double RISR = base->RISR;
          double PTISR = base->PTISR;
          double PISR = base->PISR;
          // JET
          //RISR = base->RISR_JET_ISR;
          //PTISR = base->PTISR_JET_ISR;

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
          //if(MET < 150.) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "MET > 150";
          
          // apply trigger to data and FullSim events
          if(!base->SingleElectrontrigger && !base->SingleMuontrigger && !base->DoubleElectrontrigger && !base->DoubleMuontrigger && !base->EMutrigger && !is_FastSim) // ATLAS
            continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "Lepton triggers";
          //if(!base->SingleElectrontrigger && !base->SingleMuontrigger && !is_FastSim) // ATLAS
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "1L Lepton trigger";
          //if(!base->DoubleElectrontrigger && !base->DoubleMuontrigger && !base->EMutrigger && !is_FastSim) // ATLAS
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "2L Lepton trigger";

          //if(!base->METORtrigger && !is_FastSim) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "MET trigger";

          //if(PTISR < 250.) // PreSelection
	  //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "PTISR > 250";

          // Cleaning cuts...
          double dphiCMI = base->dphiCMI;
          double PTCM = base->PTCM;
          double x = fabs(dphiCMI);
          
          // PreSelection
          //if(PTCM > 200.)
          //  continue;
          //if(PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
          //   -2.777*x*x+1.388*x+0.8264 > 0.)
          //  continue;
          //if(PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
          //   -1.5625*x*x+7.8125*x-8.766 > 0.)
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "Cleaning Cuts";
          // End of Cleaning cuts...
            
          double dphiMET_V = base->dphiMET_V;
          //if(fabs(base->dphiMET_V) > acos(-1.)/2.) continue; // PreSelection
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "#Delta#phi(MET, V) < #frac{#pi}{2}";
            
          //if(RISR < 0.4 || RISR > 0.7) // CR
          //if(RISR < 0.7 || RISR > 1.0)
          //if(RISR < 0.5 || RISR > 1.0) // PreSelection
          //if(RISR < 0.5 || RISR > 1.0) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "R_{ISR} > 0.5";

          //if(base->RISR_JET < 0.7) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "R_{ISR} > 0.7";

          //if(base->RISR_JET < 0.8) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "R_{ISR} > 0.8";

          //if(base->RISR_JET < 0.9) // PreSelection
          //  continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "R_{ISR} > 0.9";

          //if(Nlep + NjetS < 3) continue; // min 3 objects in S
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "3 Obj S";

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
          if(nBL > 0) continue; // no bronze leps
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "!Bronze";
          //if(nSL > 0) continue; // no silver leps
          //if(nGL > 0) continue; // no gold leps
          //if(nSL > 3) continue; // no more than N silver leps
          if((abs(base->PDGID_lep->at(0)) == 13 && base->PT_lep->at(0) < 26.) || (abs(base->PDGID_lep->at(0)) == 11 && base->PT_lep->at(0) < 30.)) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "p_{T}^{1} > trig";
          //if(base->PT_lep->at(1) < 12.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "p_{T}^{2} > 12.5";
          //if(Nlep > 2)
          //  if(base->PT_lep->at(2) < 7.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "p_{T}^{3} > 7.5";

          //if(NbjetISR + NbjetS != 2) continue; // CR
          if(NbjetISR + NbjetS > 1) continue; // SR & ATLAS 'B-Veto'
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "B-Veto";

          double maxSIP3D = std::max({base->SIP3D_lep->at(0), base->SIP3D_lep->at(1), base->SIP3D_lep->at(2)});
          if(maxSIP3D > 3.5) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "SIP3DMax < 3.5";

          //if(nGL < 1) continue; // at least N gold leps
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "Lep Qual: GXX";
          //if(nGL < 2) continue; // at least N gold leps
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "Lep Qual: GGX";
          //if(nGL < 3) continue; // at least N gold leps
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "Lep Qual: GGG";

          //if(NjetS != 0) continue; // SR
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "0J S";

          //if(NjetS == 0) continue; // SR
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "!0J S";

          for(int i = 0; i < Nlep; i++){
            for(int j = 0; j < base->genNlep; j++){
              TLorentzVector TLV_Cand_Child_Reco;
              TLV_Cand_Child_Reco.SetPtEtaPhiM(base->PT_lep->at(i), base->Eta_lep->at(i), base->Phi_lep->at(i), base->M_lep->at(i));
              TLorentzVector TLV_Cand_Child_Gen;
              TLV_Cand_Child_Gen.SetPtEtaPhiM(
                (*base->genPT_lep)[(*base->Index_lep)  [i]],
                (*base->genEta_lep)[(*base->Index_lep) [i]],
                (*base->genPhi_lep)[(*base->Index_lep) [i]],
                (*base->genM_lep)[(*base->Index_lep)   [i]]
              );
              double DeltaR = TLV_Cand_Child_Gen.DeltaR(TLV_Cand_Child_Reco);
              hist_DeltaRLeptonMatching->Fill(DeltaR, weight);
            }
          }
          
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
          std::vector<L_Cand> V_lep_cands_reco;
          for(int i = 0; i < Nlep-1; i++){
            for(int j = i+1; j < Nlep; j++){
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
              if((*base->Index_lep)[i] >= 0){
                lep_1.SetMomPDGID((*base->genMomPDGID_lep)[(*base->Index_lep)[i]]);
                lep_1.SetGenMomIndex((*base->genMomIndex_lep)[(*base->Index_lep)[i]]);
              }
              else{
                lep_1.SetMomPDGID(0);
                lep_1.SetGenMomIndex(-1);
              }
              if((*base->Index_lep)[j] >= 0){
                lep_2.SetMomPDGID((*base->genMomPDGID_lep)[(*base->Index_lep)[j]]);
                lep_2.SetGenMomIndex((*base->genMomIndex_lep)[(*base->Index_lep)[j]]);
              }
              else{
                lep_2.SetMomPDGID(0);
                lep_2.SetGenMomIndex(-1);
              }
              lep_1.SetPDGID(base->PDGID_lep->at(i));
              lep_2.SetPDGID(base->PDGID_lep->at(j));
              lep_1.SetIndex(i);
              lep_2.SetIndex(j);
              lep_1.SetCharge((base->Charge_lep->at(i) > 0 ? 1 : -1));
              lep_2.SetCharge((base->Charge_lep->at(j) > 0 ? 1 : -1));
              ParticleList cand_list_parts;
              cand_list_parts.push_back(lep_1);
              cand_list_parts.push_back(lep_2);
              L_Cand cand(cand_list_parts);
              cand_matching(cand);
              //cand.SetFlavor(list_leps[i].Flavor());
              //TLorentzVector TLV_Cand = cand.TLV();
              //TLorentzVector TLV_Cand_CM = TLV_Cand;
              //TLV_Cand_CM.Boost(-beta_CM);
              //double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
              // CandCuts
              //if(cand.M() < CandMinMass || cand.M() > CandMaxMass) continue; 
              //hist_CandCosDecayAngle_BetaZLAB_PreCut->Fill(cand.CosDecayAngle(), TLV_Cand.Pz()/TLV_Cand.E(), weight);
              //if(cand.CosDecayAngle() > 0.8) continue;
              //if(TLV_Cand.Pz()/TLV_Cand.E() > 0.9) continue;
              //if(PZPara < 0.) continue;
              //if(cand.Match() != kMatched) continue;

              //bool i_a = inVec((*base->index_lep_a),i);
              //bool j_a = inVec((*base->index_lep_a),j);
              //if(i_a == j_a || (Nlep == 2 && NjetS == 0)){
              //  cand.SetSameHemi(true);
              //}
              //else cand.SetSameHemi(false);
              //if(CandSameHemi && !cand.IsSameHemi()) continue;

              V_lep_cands_reco.push_back(cand);
            }
          }

          std::vector<L_Cand> V_lep_cands_gen;
          for(int i = 0; i < base->genNlep; i++){
            for(int j = i+1; j < base->genNlep; j++){
              Particle lep_1;
              lep_1.SetPtEtaPhiM( base->genPT_lep->at(i),
                                  base->genEta_lep->at(i),
                                  base->genPhi_lep->at(i),
                                  std::max(0.,base->genM_lep->at(i)) );
              Particle lep_2;
              lep_2.SetPtEtaPhiM( base->genPT_lep->at(j),
                	            base->genEta_lep->at(j),
                	            base->genPhi_lep->at(j),
                	            std::max(0.,base->genM_lep->at(j)) );
              lep_1.SetMomPDGID((*base->genMomPDGID_lep)[i]);
              lep_1.SetGenMomIndex((*base->genMomIndex_lep)[i]);
              lep_2.SetMomPDGID((*base->genMomPDGID_lep)[j]);
              lep_2.SetGenMomIndex((*base->genMomIndex_lep)[j]);
              //if(lep_1.GenMomIndex() != lep_2.GenMomIndex() && lep_1.GenMomIndex() >= 0) continue;
              //if(lep_1.MomPDGID() != 23) continue;
              lep_1.SetIndex(i);
              lep_2.SetIndex(j);
              lep_1.SetPDGID(base->genPDGID_lep->at(i));
              lep_2.SetPDGID(base->genPDGID_lep->at(j));
              lep_1.SetCharge((base->genPDGID_lep->at(i) > 0 ? 1 : -1));
              lep_2.SetCharge((base->genPDGID_lep->at(j) > 0 ? 1 : -1));
              ParticleList cand_list_parts;
              cand_list_parts.push_back(lep_1);
              cand_list_parts.push_back(lep_2);
              L_Cand cand(cand_list_parts);
              //cand.SetFlavor(abs(base->genPDGID_lep->at(i)) == 11 ? kElec : kMuon);
              //TLorentzVector TLV_Cand = cand.TLV();
              //TLorentzVector TLV_Cand_CM = TLV_Cand;
              //TLV_Cand_CM.Boost(-beta_CM);
              //double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
              // CandCuts
              //if(cand.M() < CandMinMass || cand.M() > CandMaxMass) continue; 
              //hist_CandCosDecayAngle_BetaZLAB_PreCut->Fill(cand.CosDecayAngle(), TLV_Cand.Pz()/TLV_Cand.E(), weight);
              //if(CosDecayAngleCM > 0.8) continue;
              //if(cand.CosDecayAngle() > 0.8) continue;
              //if(TLV_Cand.Pz()/TLV_Cand.E() > 0.9) continue;
              //if(PZPara < 0.) continue;
              V_lep_cands_gen.push_back(cand);
            }
          }

          std::vector<L_Cand> V_lep_cands = V_lep_cands_reco;
          //std::vector<L_Cand> V_lep_cands = V_lep_cands_gen;
          int N_V_lep_cands = V_lep_cands.size();

          //gen matched candidates
          //for(int i = 0; i < N_V_lep_cands; i++){
          //  if(V_lep_cands[i].Match() == kMatched){
          //    for(int j = 0; j < 2; j++){
          //      Particle Cand_Child_Reco = V_lep_cands[i][j];
          //      TLorentzVector TLV_Cand_Child_Reco;
          //      TLV_Cand_Child_Reco.SetPtEtaPhiM(Cand_Child_Reco.Pt(), Cand_Child_Reco.Eta(), Cand_Child_Reco.Phi(), Cand_Child_Reco.M());
          //      TLorentzVector TLV_Cand_Child_Gen;
          //      TLV_Cand_Child_Gen.SetPtEtaPhiM(
          //        (*base->genPT_lep)[(*base->Index_lep)[Cand_Child_Reco.Index()]],
          //        (*base->genEta_lep)[(*base->Index_lep)[Cand_Child_Reco.Index()]],
          //        (*base->genPhi_lep)[(*base->Index_lep)[Cand_Child_Reco.Index()]],
          //        (*base->genM_lep)[(*base->Index_lep)[Cand_Child_Reco.Index()]]
          //      );
          //      //hist_DeltaRLeptonMatching->Fill(TLV_Cand_Child_Gen.DeltaR(TLV_Cand_Child_Reco), weight);
          //    }
          //  }
          //}

          //if(N_V_lep_cands == 0) continue;
          int min_cand_index = 0;
          int max_cand_index = 0;
          for(int i = 1; i < N_V_lep_cands; i++){
            if(V_lep_cands[i].Mass() < V_lep_cands[min_cand_index].Mass()) // cand with min mass
              min_cand_index = i;
            if(V_lep_cands[i].Mass() > V_lep_cands[max_cand_index].Mass()) // cand with max mass
              max_cand_index = i;
          }
          if(N_V_lep_cands > 0){
            L_Cand cand = V_lep_cands[min_cand_index];
            L_Cand MAXcand = V_lep_cands[max_cand_index];

            // CandCuts
            TLorentzVector TLV_Cand = cand.TLV();
            TLorentzVector TLV_Cand_CM = TLV_Cand;
            TLV_Cand_CM.Boost(-beta_CM);
            double PZPara = TLV_Cand_CM.Vect().Dot(S_CM.Vect().Unit());
            //if(cand.M() <= CandMinMass || cand.M() > CandMaxMass) continue; 
            //if(MAXcand.M() <= CandMinMass || MAXcand.M() > CandMaxMass) continue; 
            hist_CandCosDecayAngle_BetaZLAB_PreCut->Fill(cand.CosDecayAngle(), TLV_Cand.Pz()/TLV_Cand.E(), weight);
            //if(fabs(cand.CosDecayAngle()) > 0.8) continue;
            //if(TLV_Cand.Pz()/TLV_Cand.E() > 0.9) continue;
            //if(PZPara < 0.) continue;
            //if(cand.Match() != kMatched) continue;

            bool i_a = inVec((*base->index_lep_a),cand[0].Index());
            bool j_a = inVec((*base->index_lep_a),cand[1].Index());
            if(i_a == j_a || (Nlep == 2 && NjetS == 0)){
              cand.SetSameHemi(true);
            }
            else cand.SetSameHemi(false);
            //if(CandSameHemi && !cand.IsSameHemi()) continue;
            // Flavor and/or Charge Req
            //if(OSSF && !(abs(cand.PL()[0].PDGID()) == abs(cand.PL()[1].PDGID()) && cand.PL()[0].Charge() != cand.PL()[1].Charge())) continue;
            //if(OSOF && !(abs(cand.PL()[0].PDGID()) != abs(cand.PL()[1].PDGID()) && cand.PL()[0].Charge() != cand.PL()[1].Charge())) continue;
            //if(SSOF && !(abs(cand.PL()[0].PDGID()) != abs(cand.PL()[1].PDGID()) && cand.PL()[0].Charge() == cand.PL()[1].Charge())) continue;
            //if(SSSF && !(abs(cand.PL()[0].PDGID()) == abs(cand.PL()[1].PDGID()) && cand.PL()[0].Charge() == cand.PL()[1].Charge())) continue;

            TVector3 beta_Cand_CM = TLV_Cand_CM.BoostVector();
            Particle part_Cand_Child = cand.Cand_PartPlus();
            TLorentzVector TLV_Cand_Child;
            TLV_Cand_Child.SetPtEtaPhiM(part_Cand_Child.Pt(), part_Cand_Child.Eta(), part_Cand_Child.Phi(), part_Cand_Child.M());
            TLorentzVector TLV_Cand_Child_CM = TLV_Cand_Child;
            TLV_Cand_Child_CM.Boost(-beta_CM);
            TLV_Cand_Child_CM.Boost(-beta_Cand_CM);
            double CosDecayAngleCM = fabs(TLV_Cand_Child_CM.Vect().Unit().Dot(beta_Cand_CM.Unit())); // 'original' using CM
            double CosDecayAngle = fabs(cand.CosDecayAngle()); // decay angle without using CM
            double CandML = TLV_Cand.M();
            double Beta = cand.Beta();
            double BetaCM = TLV_Cand_CM.P()/TLV_Cand_CM.E();
            double CandECM = TLV_Cand_CM.E();
            double PZzLAB = TLV_Cand.Pz();
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
            double MRZPara = CandML/RZPara;
            double MRZPerp = CandML/RZPerp;
            double MRZParaLABMET = CandML/RZParaLABMET;
            double MRZPerpLABMET = CandML/RZPerpLABMET;
            double MPZPara = fabs(PZPara/CandML);
            double MPZPerp = fabs(PZPerp/CandML);
            double MPZParaLABMET = fabs(PZParaLABMET/TLV_Cand_Transverse.Mag());
            double MPZPerpLABMET = fabs(PZPerpLABMET/TLV_Cand_Transverse.Mag());
            double BetaZPara = PZPara/CandECM;
            double BetaZPerp = PZPerp/CandECM;
            double BetaZLAB = TLV_Cand.Pz()/TLV_Cand.E();
            double MInv = CandML*RISR/RZPara;
            double CandMdivPTISR = CandML/PTISR;
            double RZParaM = RZPara/CandML;
            TLorentzVector TLV_MAXCand = MAXcand.TLV();
            double MAXcandML = TLV_MAXCand.M();
            double MAXcandBetaZLAB = TLV_MAXCand.Pz()/TLV_MAXCand.E();
            double MAXCosDecayAngle = fabs(MAXcand.CosDecayAngle());
            TLorentzVector TLV_Cand_Child_a = cand.TLV_Part(0);
            TLorentzVector TLV_Cand_Child_b = cand.TLV_Part(1);
            double pa[3] = {TLV_Cand_Child_a.M(), TLV_Cand_Child_a.Px(), TLV_Cand_Child_a.Py()};
            double pb[3] = {TLV_Cand_Child_b.M(), TLV_Cand_Child_b.Px(), TLV_Cand_Child_b.Py()};
            double pmiss[3] = {0.0, TV3_MET.Px(), TV3_MET.Py()}; // Invisible particle mass = 0
            mt2calc.set_momenta(pa, pb, pmiss);
            mt2calc.set_mn(0.0); // assume massless invisible particles like neutrinos
            double MT2 = mt2calc.get_mt2();
            TLorentzVector MAXTLV_Cand_Child_a = MAXcand.TLV_Part(0);
            TLorentzVector MAXTLV_Cand_Child_b = MAXcand.TLV_Part(1);
            double MAXpa[3] = {MAXTLV_Cand_Child_a.M(), MAXTLV_Cand_Child_a.Px(), MAXTLV_Cand_Child_a.Py()};
            double MAXpb[3] = {MAXTLV_Cand_Child_b.M(), MAXTLV_Cand_Child_b.Px(), MAXTLV_Cand_Child_b.Py()};
            mt2calc.set_momenta(MAXpa, MAXpb, pmiss);
            mt2calc.set_mn(0.0); // assume massless invisible particles like neutrinos
            double MT2MAX = mt2calc.get_mt2();
            double MCon_CM = MAXcand.MCon(TLV_CM);
            TLorentzVector TLV_Cand_LEPS;
            //std::vector<TLorentzVector> vect_TLV_Cand_LEPS;
            //vect_TLV_Cand_LEPS.push_back(cand.TLV_Part(0));
            //vect_TLV_Cand_LEPS.push_back(cand.TLV_Part(1));
            //vect_TLV_Cand_LEPS.push_back(MAXcand.TLV_Part(0));
            //vect_TLV_Cand_LEPS.push_back(MAXcand.TLV_Part(1));
            //std::vector<TLorentzVector> vect_TLV_Cand_LEPS_unq; // unique
            //for (const auto& lep : vect_TLV_Cand_LEPS) {
            //    bool is_unique = true;
            //    for (const auto& u : vect_TLV_Cand_LEPS_unq) {
            //      if ((lep - u).Vect().Mag() < 1e-6 && std::abs(lep.E() - u.E()) < 1e-6) {
            //        is_unique = false;
            //        break;
            //    }
            //  }
            //  if (is_unique) vect_TLV_Cand_LEPS_unq.push_back(lep);
            //}
            //for (const auto& lep : vect_TLV_Cand_LEPS_unq) TLV_Cand_LEPS += lep;
            for(int i = 0; i < Nlep; i++){
              TLorentzVector dummy_TLV;
              dummy_TLV.SetPtEtaPhiM( base->PT_lep->at(i),
                                      base->Eta_lep->at(i),
                                      base->Phi_lep->at(i),
                                      std::max(0.,base->M_lep->at(i)) );
               TLV_Cand_LEPS += dummy_TLV;
            }
            double MCon_LEPS = MAXcand.MCon(TLV_Cand_LEPS);
            TLorentzVector TLV_LEPS_MET = S_CM;
            for(int i = 0; i < NjetS; i++){
              index = (*base->index_jet_S)[i];
              TLorentzVector SJet;
              SJet.SetPtEtaPhiM(  base->PT_jet->at(i),
                                  base->Eta_jet->at(i),
                                  base->Phi_jet->at(i),
                                  std::max(0.,base->M_jet->at(i)) );
              SJet.Boost(-beta_CM);
              TLV_LEPS_MET -= SJet;
            } // end S-jet loop
            double MCon_LEPSMET = MAXcand.MCon(TLV_LEPS_MET);
            TLorentzVector TLV_Cand_Child_a_CM = TLV_Cand_Child_a;
            TLorentzVector TLV_Cand_Child_b_CM = TLV_Cand_Child_b;
            TLorentzVector MAXTLV_Cand_Child_a_CM = MAXTLV_Cand_Child_a;
            TLorentzVector MAXTLV_Cand_Child_b_CM = MAXTLV_Cand_Child_b;
            if(cand.PL()[0].Charge() != cand.PL()[1].Charge()){ // Opposite Sign
              TLV_Cand_Child_a_CM = cand.Cand_TLVPlus();
              TLV_Cand_Child_b_CM = cand.Cand_TLVMinus();
            } // Opposite Sign
            if(MAXcand.PL()[0].Charge() != MAXcand.PL()[1].Charge()){ // Opposite Sign
              MAXTLV_Cand_Child_a_CM = MAXcand.Cand_TLVPlus();
              MAXTLV_Cand_Child_b_CM = MAXcand.Cand_TLVMinus();
            } // Opposite Sign
            TLV_Cand_Child_a_CM.Boost(-beta_CM);
            TLV_Cand_Child_b_CM.Boost(-beta_CM);
            MAXTLV_Cand_Child_a_CM.Boost(-beta_CM);
            MAXTLV_Cand_Child_b_CM.Boost(-beta_CM);
            double CosThetaRazor = fabs(((TLV_Cand_Child_a_CM.E() - TLV_Cand_Child_b_CM.E()) * (TLV_Cand_Child_a_CM.E() - TLV_Cand_Child_b_CM.E())) /
                                   (((TLV_Cand_Child_a_CM.E() - TLV_Cand_Child_b_CM.E()) * (TLV_Cand_Child_a_CM.E() - TLV_Cand_Child_b_CM.E()))
                                   + (TLV_Cand_Child_a_CM + TLV_Cand_Child_b_CM).Mag2()));
            double MAXCosThetaRazor = fabs(((MAXTLV_Cand_Child_a_CM.E() - MAXTLV_Cand_Child_b_CM.E()) * (MAXTLV_Cand_Child_a_CM.E() - MAXTLV_Cand_Child_b_CM.E())) /
                                   (((MAXTLV_Cand_Child_a_CM.E() - MAXTLV_Cand_Child_b_CM.E()) * (MAXTLV_Cand_Child_a_CM.E() - MAXTLV_Cand_Child_b_CM.E()))
                                   + (MAXTLV_Cand_Child_a_CM + MAXTLV_Cand_Child_b_CM).Mag2()));
            TLorentzVector TLV_MAXCand_CM = TLV_MAXCand;
            TLV_MAXCand_CM.Boost(-beta_CM);
            double PZParaMAX = TLV_MAXCand_CM.Vect().Dot(S_CM.Vect().Unit());
            double RZParaMAX = PZParaMAX/S_CM.Vect().Mag();
            double MAXMInv = MAXcandML*RISR/RZParaMAX;
            hist_CandML->Fill(CandML, weight);
            hist_CandBeta->Fill(Beta, weight);
            hist_CandBetaCM->Fill(BetaCM, weight);
            hist_CandBeta_CandBetaCM->Fill(Beta, BetaCM, weight);
            hist_CandDeltaPhiMET->Fill(cand.DeltaPhi(TV3_MET), weight);
            hist_CandCosDecayAngle->Fill(CosDecayAngle, weight);
            hist_CandCosDecayAngleCM->Fill(CosDecayAngleCM, weight);
            hist_CandMdivPTISR->Fill(CandMdivPTISR, weight);
            hist_MInv->Fill(MInv, weight);
            hist_OneMinusRZPara->Fill(1-RZPara, weight);
            hist_PzCM->Fill(PzCM, weight);
            hist_BetaZLAB->Fill(BetaZLAB, weight);
            hist_RZParaM->Fill(RZParaM, weight);
            hist_RZParaDivRISR->Fill(RZPara/RISR, weight);
            hist_RZPara->Fill(RZPara, weight);
            hist_RZPerp->Fill(RZPerp, weight);
            hist_RISR_CandML->Fill(RISR, CandML, weight);
            hist_RISR_CandBeta->Fill(RISR, Beta, weight);
            hist_CandML_CandBeta->Fill(CandML, Beta, weight);
            hist_Beta_CandCosDecayAngleCM->Fill(Beta, CosDecayAngleCM, weight);
            hist_Beta_CandCosDecayAngle->Fill(Beta, CosDecayAngle, weight);
            hist_Beta_CandDeltaPhiMET->Fill(Beta, cand.DeltaPhi(TV3_MET), weight);
            hist_RISR_CandBetaCM->Fill(RISR, BetaCM, weight);
            hist_CandML_CandBetaCM->Fill(CandML, BetaCM, weight);
            hist_BetaCM_CandCosDecayAngleCM->Fill(BetaCM, CosDecayAngleCM, weight);
            hist_BetaCM_CandCosDecayAngle->Fill(BetaCM, CosDecayAngle, weight);
            hist_BetaCM_CandDeltaPhiMET->Fill(BetaCM, cand.DeltaPhi(TV3_MET), weight);
            hist_RISR_CandDeltaPhiMET->Fill(RISR, cand.DeltaPhi(TV3_MET), weight);
            hist_CandML_CandDeltaPhiMET->Fill(CandML, cand.DeltaPhi(TV3_MET), weight);
            hist_CandML_CandCosDecayAngleCM->Fill(CandML, CosDecayAngleCM, weight);
            hist_CandCosDecayAngleCM_CandDeltaPhiMET->Fill(CosDecayAngleCM, cand.DeltaPhi(TV3_MET), weight);
            hist_CandML_CandCosDecayAngle->Fill(CandML, CosDecayAngle, weight);
            hist_CandCosDecayAngle_CandDeltaPhiMET->Fill(CosDecayAngle, cand.DeltaPhi(TV3_MET), weight);
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
            //hist_CandCosDecayAngleCM_RZPara->Fill(CosDecayAngleCM, RZPara, weight);
            //hist_CandCosDecayAngleCM_RZPerp->Fill(CosDecayAngleCM, RZPerp, weight); 
            //hist_CandCosDecayAngleCM_PZAng->Fill(CosDecayAngleCM, PZAng, weight); 
            //hist_CandCosDecayAngleCM_PZPara->Fill(CosDecayAngleCM, PZPara, weight); 
            //hist_CandCosDecayAngleCM_PZPerp->Fill(CosDecayAngleCM, PZPara, weight); 
            //hist_CandCosDecayAngleCM_RZParaLABMET->Fill(CosDecayAngleCM, RZParaLABMET, weight); 
            //hist_CandCosDecayAngleCM_RZPerpLABMET->Fill(CosDecayAngleCM, RZPerpLABMET, weight); 
            //hist_CandCosDecayAngleCM_PZAngLABMET->Fill(CosDecayAngleCM, PZAngLABMET, weight); 
            //hist_CandCosDecayAngleCM_PZParaLABMET->Fill(CosDecayAngleCM, PZParaLABMET, weight); 
            //hist_CandCosDecayAngleCM_PZPerpLABMET->Fill(CosDecayAngleCM, PZPerpLABMET, weight); 
            //hist_Beta_RZPara->Fill(Beta, RZPara, weight); 
            //hist_Beta_RZPerp->Fill(Beta, RZPerp, weight); 
            //hist_Beta_PZAng->Fill(Beta, PZAng, weight); 
            //hist_Beta_PZPara->Fill(Beta, PZPara, weight); 
            //hist_Beta_PZPerp->Fill(Beta, PZPerp, weight); 
            //hist_Beta_RZParaLABMET->Fill(Beta, RZParaLABMET, weight); 
            //hist_Beta_RZPerpLABMET->Fill(Beta, RZPerpLABMET, weight); 
            //hist_Beta_PZAngLABMET->Fill(Beta, PZAngLABMET, weight); 
            //hist_Beta_PZParaLABMET->Fill(Beta, PZParaLABMET, weight); 
            //hist_Beta_PZPerpLABMET->Fill(Beta, PZPerpLABMET, weight); 
            //hist_BetaCM_RZPara->Fill(BetaCM, RZPara, weight); 
            //hist_BetaCM_RZPerp->Fill(BetaCM, RZPerp, weight); 
            //hist_BetaCM_PZAng->Fill(BetaCM, PZAng, weight); 
            //hist_BetaCM_PZPara->Fill(BetaCM, PZPara, weight); 
            //hist_BetaCM_PZPerp->Fill(BetaCM, PZPerp, weight); 
            //hist_BetaCM_RZParaLABMET->Fill(BetaCM, RZParaLABMET, weight); 
            //hist_BetaCM_RZPerpLABMET->Fill(BetaCM, RZPerpLABMET, weight); 
            //hist_BetaCM_PZAngLABMET->Fill(BetaCM, PZAngLABMET, weight); 
            //hist_BetaCM_PZParaLABMET->Fill(BetaCM, PZParaLABMET, weight); 
            //hist_BetaCM_PZPerpLABMET->Fill(BetaCM, PZPerpLABMET, weight); 
            //hist_RZPara_RZParaLABMET->Fill(RZPara, RZParaLABMET, weight);
            //hist_RZPerp_RZPerpLABMET->Fill(RZPerp, RZPerpLABMET, weight);
            //hist_PZPara_PZParaLABMET->Fill(PZPara, PZParaLABMET, weight);
            //hist_PZPerp_PZPerpLABMET->Fill(PZPerp, PZPerpLABMET, weight);
            //hist_PZPara_CandECM->Fill(PZPara, CandECM, weight);
            //hist_PZPerp_CandECM->Fill(PZPerp, CandECM, weight);
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
            hist_MRZPara_RZPara->Fill(MRZPara, RZPara, weight);
            hist_MRZPara_CandCosDecayAngleCM->Fill(MRZPara, CosDecayAngleCM, weight);
            hist_MRZPara_BetaCM->Fill(MRZPara, BetaCM, weight);
            hist_MRZPara_CandPZzLAB->Fill(MRZPara, PZzLAB, weight);
            hist_RZPara_CandPZzLAB->Fill(RZPara, PZzLAB, weight);
            hist_RZPerp_CandPZzLAB->Fill(RZPerp, PZzLAB, weight);
            hist_RISR_RZPara->Fill(RISR, RZPara, weight);
            hist_RISR_MRZPara->Fill(RISR, MRZPara, weight);
            hist_RISR_CandPZzLAB->Fill(RISR, PZzLAB, weight);
            hist_RISR_BetaZLAB->Fill(RISR, BetaZLAB, weight);
            hist_CandML_BetaZLAB->Fill(CandML, BetaZLAB, weight);
            hist_RZPara_BetaZLAB->Fill(RZPara, BetaZLAB, weight);
            hist_MRZPara_BetaZLAB->Fill(MRZPara, BetaZLAB, weight);
            hist_RZPerp_BetaZLAB->Fill(RZPerp, BetaZLAB, weight);
            hist_CandCosDecayAngle_BetaZLAB->Fill(CosDecayAngle, BetaZLAB, weight);
            hist_CandCosDecayAngleCM_BetaZLAB->Fill(CosDecayAngleCM, BetaZLAB, weight);
            hist_RZParaM_RZPara->Fill(RZParaM, RZPara, weight);
            hist_RZParaM_CandCosDecayAngleCM->Fill(RZParaM, CosDecayAngle, weight);
            hist_RZParaM_BetaCM->Fill(RZParaM, BetaCM, weight);
            hist_RZParaM_CandPZzLAB->Fill(RZParaM, PZzLAB, weight);
            hist_RZParaM_BetaZLAB->Fill(RZParaM, BetaZLAB, weight);
            //hist_MRZPara_OneMinusRZPara->Fill(MRZPara, 1.-RZPara, weight);
            //hist_OneMinusRZPara_CandPZzLAB->Fill(1.-RZPara, PZzLAB, weight);
            //hist_RISR_OneMinusRZPara->Fill(RISR, 1.-RZPara, weight);
            //hist_OneMinusRZPara_BetaZLAB->Fill(1.-RZPara, BetaZLAB, weight);
            //hist_BetaCM_OneMinusRZPara->Fill(BetaCM, 1.-RZPara, weight);
            //hist_CandCosDecayAngle_OneMinusRZPara->Fill(CosDecayAngle, 1.-RZPara, weight);
            //hist_Beta_OneMinusRZPara->Fill(Beta, 1-RZPara, weight);
            //hist_CandCosDecayAngleCM_OneMinusRZPara->Fill(CosDecayAngleCM, 1.-RZPara, weight);
            //hist_MRZPara_CandMdivPTISR->Fill(MRZPara, CandMdivPTISR, weight);
            //hist_CandMdivPTISR_CandPZzLAB->Fill(CandMdivPTISR, PZzLAB, weight);
            //hist_RISR_CandMdivPTISR->Fill(RISR, CandMdivPTISR, weight);
            //hist_CandMdivPTISR_BetaZLAB->Fill(CandMdivPTISR, BetaZLAB, weight);
            //hist_BetaCM_CandMdivPTISR->Fill(BetaCM, CandMdivPTISR, weight);
            //hist_CandCosDecayAngle_CandMdivPTISR->Fill(CosDecayAngle, CandMdivPTISR, weight);
            //hist_Beta_CandMdivPTISR->Fill(Beta, CandMdivPTISR, weight);
            //hist_CandCosDecayAngleCM_CandMdivPTISR->Fill(CosDecayAngleCM, CandMdivPTISR, weight);
            hist_MRZPara_MInv->Fill(MRZPara, MInv, weight);
            hist_MInv_CandPZzLAB->Fill(MInv, PZzLAB, weight);
            hist_RISR_MInv->Fill(RISR, MInv, weight);
            hist_MInv_BetaZLAB->Fill(MInv, BetaZLAB, weight);
            hist_BetaCM_MInv->Fill(BetaCM, MInv, weight);
            hist_CandCosDecayAngle_MInv->Fill(CosDecayAngle, MInv, weight);
            hist_Beta_MInv->Fill(Beta, MInv, weight);
            hist_CandCosDecayAngleCM_MInv->Fill(CosDecayAngleCM, MInv, weight);
            //hist_CandML_RZParaM->Fill(CandML, RZParaM, weight);
            //hist_CandML_OneMinusRZPara->Fill(CandML, 1.-RZPara, weight);
            //hist_CandML_CandMdivPTISR->Fill(CandML, CandMdivPTISR, weight);
            hist_CandML_MInv->Fill(CandML, MInv, weight);
            //hist_RZParaDivRISR_RZParaM->Fill(RZPara/RISR, RZParaM, weight);
            //hist_RZParaDivRISR_CandMdivPTISR->Fill(RZPara/RISR, CandMdivPTISR, weight);
            //hist_RZParaDivRISR_MInv->Fill(RZPara/RISR, MInv, weight);
            //hist_RZParaDivRISR_BetaZLAB->Fill(RZPara/RISR, BetaZLAB, weight);
            //hist_RZParaDivRISR_MRZPara->Fill(RZPara/RISR, MRZPara, weight);
            //hist_RZParaDivRISR_RZPara->Fill(RZPara/RISR, RZPara, weight);
            //hist_RZParaDivRISR_OneMinusRZPara->Fill(RZPara/RISR, 1-RZPara, weight);
            //hist_RZParaDivRISR_Beta->Fill(RZPara/RISR, Beta, weight);
            //hist_RZParaDivRISR_BetaCM->Fill(RZPara/RISR, BetaCM, weight);
            //hist_RZParaDivRISR_CandCosDecayAngle->Fill(RZPara/RISR, CosDecayAngle, weight);
            //hist_RZParaDivRISR_CandCosDecayAngleCM->Fill(RZPara/RISR, CosDecayAngleCM, weight);
            //hist_CandMLDivRISR_CandML->Fill(RZPara/RISR, CandML, weight);
            hist_CandML_MAXcandML->Fill(CandML, MAXcandML, weight);
            hist_BetaZLAB_MAXcandBetaZLAB->Fill(BetaZLAB, MAXcandBetaZLAB, weight);
            hist_CosDecayAngle_MAXcandCosDecayAngle->Fill(CosDecayAngle, MAXCosDecayAngle, weight);
            hist_MAXcandBetaZLAB_MAXcandCosDecayAngle->Fill(MAXcandBetaZLAB, MAXCosDecayAngle, weight);
            hist_MAXcandML_MT2->Fill(MAXcandML, MT2, weight);
            hist_MAXcandML_MConCM->Fill(MAXcandML, MCon_CM, weight);
            hist_MAXcandML_MConLEPS->Fill(MAXcandML, MCon_LEPS, weight);
            hist_MT2_MConCM->Fill(MT2, MCon_CM, weight);
            hist_MT2_MConLEPS->Fill(MT2, MCon_LEPS, weight);
            hist_MConCM_MConLEPS->Fill(MCon_CM, MCon_LEPS, weight);
            hist_CandML_MT2->Fill(CandML, MT2, weight);
            hist_CandML_MConCM->Fill(CandML, MCon_CM, weight);
            hist_CandML_MConLEPS->Fill(CandML, MCon_LEPS, weight);
            hist_CandML_MT2MAX->Fill(CandML, MT2MAX, weight);
            hist_CandML_MConLEPSMET->Fill(CandML, MCon_LEPSMET, weight);
            hist_CandML_CosThetaRazor->Fill(CandML, CosThetaRazor, weight);
            hist_CandML_MAXCosThetaRazor->Fill(CandML, MAXCosThetaRazor, weight);
            hist_CandML_MAXMInv->Fill(CandML, MAXMInv, weight);
            hist_MAXCandML_MT2MAX->Fill(MAXcandML, MT2MAX, weight);
            hist_MAXCandML_MCon_LEPSMET->Fill(MAXcandML, MCon_LEPSMET, weight);
            hist_MAXCandML_CosThetaRazor->Fill(MAXcandML, CosThetaRazor, weight);
            hist_MAXCandML_MAXCosThetaRazor->Fill(MAXcandML, MAXCosThetaRazor, weight);
            hist_MAXCandML_MAXMInv->Fill(MAXcandML, MAXMInv, weight);
            hist_MAXCandML_RISR->Fill(MAXcandML, RISR, weight);
            hist_MAXCandML_MInv->Fill(MAXcandML, MInv, weight);
            hist_MConCM_MT2MAX->Fill(MCon_CM, MT2MAX, weight);
            hist_MConCM_MConLEPSMET->Fill(MCon_CM, MCon_LEPSMET, weight);
            hist_MConCM_CosThetaRazor->Fill(MCon_CM, CosThetaRazor, weight);
            hist_MConCM_MAXCosThetaRazor->Fill(MCon_CM, MAXCosThetaRazor, weight);
            hist_MConCM_MAXMInv->Fill(MCon_CM, MAXMInv, weight);
            hist_MConCM_RISR->Fill(MCon_CM, RISR, weight);
            hist_MConCM_MInv->Fill(MCon_CM, MInv, weight);
            hist_MConLEPS_MT2MAX->Fill(MCon_LEPS, MT2MAX, weight);
            hist_MConLEPS_MConLEPSMET->Fill(MCon_LEPS, MCon_LEPSMET, weight);
            hist_MConLEPS_CosThetaRazor->Fill(MCon_LEPS, CosThetaRazor, weight);
            hist_MConLEPS_MAXCosThetaRazor->Fill(MCon_LEPS, MAXCosThetaRazor, weight);
            hist_MConLEPS_MAXMInv->Fill(MCon_LEPS, MAXMInv, weight);
            hist_MConLEPS_RISR->Fill(MCon_LEPS, RISR, weight);
            hist_MConLEPS_MInv->Fill(MCon_LEPS, MInv, weight);
            hist_MT2_MT2MAX->Fill(MT2, MT2MAX, weight);
            hist_MT2_MConLEPSMET->Fill(MT2, MCon_LEPSMET, weight);
            hist_MT2_CosThetaRazor->Fill(MT2, CosThetaRazor, weight);
            hist_MT2_MAXCosThetaRazor->Fill(MT2, MAXCosThetaRazor, weight);
            hist_MT2_MAXMInv->Fill(MT2, MAXMInv, weight);
            hist_MT2_RISR->Fill(MT2, RISR, weight);
            hist_MT2_MInv->Fill(MT2, MInv, weight);
            hist_MT2MAX_MConLEPSMET->Fill(MT2MAX, MCon_LEPSMET, weight);
            hist_MT2MAX_CosThetaRazor->Fill(MT2MAX, CosThetaRazor, weight);
            hist_MT2MAX_MAXCosThetaRazor->Fill(MT2MAX, MAXCosThetaRazor, weight);
            hist_MT2MAX_MAXMInv->Fill(MT2MAX, MAXMInv, weight);
            hist_MT2MAX_RISR->Fill(MT2MAX, RISR, weight);
            hist_MT2MAX_MInv->Fill(MT2MAX, MInv, weight);
            hist_MConLEPSMET_CosThetaRazor->Fill(MCon_LEPSMET, CosThetaRazor, weight);
            hist_MConLEPSMET_MAXCosThetaRazor->Fill(MCon_LEPSMET, MAXCosThetaRazor, weight);
            hist_MConLEPSMET_MAXMInv->Fill(MCon_LEPSMET, MAXMInv, weight);
            hist_MConLEPSMET_RISR->Fill(MCon_LEPSMET, RISR, weight);
            hist_MConLEPSMET_MInv->Fill(MCon_LEPSMET, MInv, weight);
            hist_CosThetaRazor_MAXCosThetaRazor->Fill(CosThetaRazor, MAXCosThetaRazor, weight);
            hist_CosThetaRazor_MAXMInv->Fill(CosThetaRazor, MAXMInv, weight);
            hist_CosThetaRazor_RISR->Fill(CosThetaRazor, RISR, weight);
            hist_CosThetaRazor_MInv->Fill(CosThetaRazor, MInv, weight);
            hist_MAXCosThetaRazor_MAXMInv->Fill(MAXCosThetaRazor, MAXMInv, weight);
            hist_MAXCosThetaRazor_RISR->Fill(MAXCosThetaRazor, RISR, weight);
            hist_MAXCosThetaRazor_MInv->Fill(MAXCosThetaRazor, MInv, weight);
            hist_MAXMInv_RISR->Fill(MAXMInv, RISR, weight);
            hist_MAXMInv_MInv->Fill(MAXMInv, MInv, weight);
            hist_MAXCandML_MAXCosDecayAngle->Fill(MAXcandML, MAXCosDecayAngle, weight);
          } // at least one lep cand

          // Implement new RJR tree
          //LAB.ClearEvent();
          //vector<RFKey> lep_RJR_keys;
          //INV.SetLabFrameThreeVector(TV3_MET);
          //TLorentzVector TLV_RJR_Jets;
          //for(int i = 0; i < Njet; i++){
          //  TLorentzVector jet; jet.SetPtEtaPhiM( base->PT_jet->at(i),
          //                          base->Eta_jet->at(i),
          //                          base->Phi_jet->at(i),
          //                          std::max(0.,base->M_jet->at(i)) );
          //  TLV_RJR_Jets += jet;
          //}
          //ISR.SetLabFrameFourVector(TLV_RJR_Jets);
          //for(int i = 0; i < Nlep; i++){
          //  TLorentzVector lep; lep.SetPtEtaPhiM( base->PT_lep->at(i),
          //                          base->Eta_lep->at(i),
          //                          base->Phi_lep->at(i),
          //                          std::max(0.,base->M_lep->at(i)) );
          //  RFKey key = COMB_L.AddLabFrameFourVector(lep);
          //  lep_RJR_keys.push_back(key);
          //}
          //if(!LAB.AnalyzeEvent()) cout << "Problem with RJR Analyze Event \n";
          //bool BSideIsA = false;
          //if(Lb.GetFourVector().M() > La.GetFourVector().M()) BSideIsA = true;

          //TVector3 vPISR = S.GetFourVector(CM).Vect();
          double PTISR_LEP = base->PTISR_LEP;//vPISR.Pt();
          //TVector3 vPINV = (Ia.GetFourVector(CM)+Ib.GetFourVector(CM)).Vect();

          //double MPa = Pa.GetFourVector().M();
          //double MPb = Pb.GetFourVector().M();
          //if(BSideIsA){ MPa = Pb.GetFourVector().M(); MPb = Pa.GetFourVector().M(); }
          double MPa = base->MPa_LEP;
          double MPb = base->MPb_LEP;
          double RISR_LEP = base->RISR_LEP;//fabs(vPINV.Dot(vPISR.Unit())) / vPISR.Mag();
          double MQ_LEP = base->MQ_LEP;//sqrt(MPa*MPa+MPb*MPb)/sqrt(2.);
          double MS_LEP = base->MS_LEP;//S.GetFourVector().M();
          double gamma_LEP = base->gamma_LEP;//2.*MQ_LEP/MS_LEP;
          double CosDecayAngle_S_LEP = base->CosDecayAngle_S_LEP;//S.GetCosDecayAngle();
          //double CosDecayAngle_Pa_LEP = Pa.GetCosDecayAngle();
          //double CosDecayAngle_Pb_LEP = Pb.GetCosDecayAngle();
          double CosDecayAngle_Va_LEP = 0.;//sLa.GetCosDecayAngle();
          double CosDecayAngle_Vb_LEP = 0.;//sLb.GetCosDecayAngle();
          //if(BSideIsA){ CosDecayAngle_Pa_LEP = Pb.GetCosDecayAngle(); CosDecayAngle_Pb_LEP = Pa.GetCosDecayAngle(); }
          double CosDecayAngle_Pa_LEP = base->CosDecayAngle_Pa_LEP;
          double CosDecayAngle_Pb_LEP = base->CosDecayAngle_Pb_LEP;
          //double MVa = La.GetFourVector().M();
          //double MVb = Lb.GetFourVector().M();
          //if(BSideIsA){ MVa = Lb.GetFourVector().M(); MVb = La.GetFourVector().M(); }
          double MVa = base->MVa;
          double MVb = base->MVb;
          //double PTS_CM = S.GetFourVector(CM).Pt();
          double PTS_CM = base->PTS_CM_LEP;
          //TLorentzVector TLV_L_CMLEP = La.GetFourVector(CM) + Lb.GetFourVector(CM);
          //double RZPara_LEP = TLV_L_CMLEP.Vect().Dot(S.GetFourVector(CM).Vect().Unit())/S.GetFourVector(CM).Vect().Mag();
          //double MINV_LEP = TLV_L_CMLEP.M()*RISR_LEP/RZPara_LEP;
          double MINV_LEP = base->MINV_LEP;

          //TLorentzVector Ia_S = Ia.GetFourVector(S);
          //TLorentzVector Ib_S = Ib.GetFourVector(S);
          //TLorentzVector La_S = La.GetFourVector(S);
          //TLorentzVector Lb_S = Lb.GetFourVector(S);
          //TLorentzVector Ia_S0;
          //TLorentzVector Ib_S0;
          //TLorentzVector La_S0;
          //TLorentzVector Lb_S0;
          //Ia_S0.SetPtEtaPhiM(Ia_S.Pt(), Ia_S.Eta(), Ia_S.Phi(), 0.);
          //Ib_S0.SetPtEtaPhiM(Ib_S.Pt(), Ib_S.Eta(), Ib_S.Phi(), 0.);
          //La_S0.SetPtEtaPhiM(La_S.Pt(), La_S.Eta(), La_S.Phi(), 0.);
          //Lb_S0.SetPtEtaPhiM(Lb_S.Pt(), Lb_S.Eta(), Lb_S.Phi(), 0.);
          //double MS_S0 = (Ia_S0+Ib_S0+La_S0+Lb_S0).M();
          //double MV_S0 = (La_S0+Lb_S0).M();
          //double MQ_S0 = sqrt(((Ia_S0+La_S0).M2()+(Ib_S0+Lb_S0).M2())/2.);
          //double gamma_S0 = 2.*MQ_S0/MS_S0;
          double MS_S0 = base->MS_S0_LEP;
          double MV_S0 = base->MV_S0_LEP;
          double MQ_S0 = base->MQ_S0_LEP;
          double gamma_S0 = base->gamma_S0_LEP;

          //TLorentzVector Ia_Pa = Ia.GetFourVector(Pa);
          //TLorentzVector Ib_Pb = Ib.GetFourVector(Pb);
          double MPTilde = base->MPTilde_LEP;//Ia_Pa.Vect().Mag() + Ib_Pb.Vect().Mag();
          //TLorentzVector Pa_S = Pa.GetFourVector(S);
          double MSTilde = base->MSTilde_LEP;//sqrt(MPTilde*MPTilde+Pa_S.Vect().Mag2());
          double gammaTilde = base->gammaTilde_LEP;//MPTilde/MSTilde;

          //if(RISR_LEP < 0.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "R_{ISR}^{LEP} > 0.5";

          //if(PTISR_LEP < 150.) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "PT_{ISR}^{LEP} > 150";

          if(MVa < 4 || MVa > 65.) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "4 < M_{Va}^{LEP} < 65";

          if(MS_LEP < 300.) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "MS^{LEP} > 300";

          if(CosDecayAngle_Pa_LEP+CosDecayAngle_Pb_LEP > 0.19) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "cos#theta^{Pa}+cos#theta^{Pb} < 0.19";

          if(CosDecayAngle_S_LEP > 0.8) continue;
          CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          vect_str_cutflow_labels[CF_bin] = "cos#theta^{S} < 0.8";

          //if(MQ_LEP > 150.) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "MQ^{LEP} < 150";

          //if(CosDecayAngle_Pa_LEP > 0.4 && fabs(CosDecayAngle_Pb_LEP) < 0.2) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "cos#theta^{A} > 0.4 @ |cos#theta^{B}| < 0.2";

          //if(CosDecayAngle_Pa_LEP*CosDecayAngle_Pb_LEP) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "cos#theta^{A} * cos#theta^{B} < 0.2";

          //if(MPa/MPb > 0.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "M_{Pa}^{LEP}/M_{Pb}^{LEP} > 0.5";

          //if((MS_LEP > 180. && MS_LEP < 250.) && MPa/MPb > 0.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "M_{Pa}^{LEP}/M_{Pb}^{LEP} > 0.5 MS < 250";

          //if((MS_LEP > 180. && MS_LEP < 250.) && MPa/MPb > 1.) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "M_{Pa}^{LEP}/M_{Pb}^{LEP} > 1. MS < 250";

          //if(gamma_LEP < 0.5) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "#gamma^{LEP} > 0.5";

          //if(CosDecayAngle_Pa_LEP > 0.9) continue;
          //CF_bin++; hist_CutFlow->SetBinContent(CF_bin, hist_CutFlow->GetBinContent(CF_bin)+weight); 
          //vect_str_cutflow_labels[CF_bin] = "cos#theta_{Pa}^{LEP} < 0.9";

          hist_RISRLEP_RISR->Fill(RISR_LEP, RISR, weight);
          hist_MQLEP_MSLEP->Fill(MQ_LEP, MS_LEP, weight);
          hist_gammaLEP_MQLEP->Fill(gamma_LEP, MQ_LEP, weight);
          hist_gammaLEP_MSLEP->Fill(gamma_LEP, MS_LEP, weight);
          hist_MPa_CosDecayAnglePaLEP->Fill(MPa, CosDecayAngle_Pa_LEP, weight);
          hist_MSLEP_CosDecayAngleSLEP->Fill(MS_LEP, CosDecayAngle_S_LEP, weight);
          hist_RISRLEP_PTISRLEP->Fill(RISR_LEP, PTISR_LEP, weight);
          hist_PTISR_PTISRLEP->Fill(PTISR, PTISR_LEP, weight); 
          hist_MET_MS->Fill(MET, MS_LEP, weight);
          hist_MET_MPa->Fill(MET, MPa, weight);
          hist_MET_MPb->Fill(MET, MPb, weight);
          hist_MET_MVa->Fill(MET, MVa, weight);
          hist_MET_MVb->Fill(MET, MVb, weight);
          hist_MET_PTSCM->Fill(MET, PTS_CM, weight);
          hist_MET_MSS0->Fill(MET, MS_S0, weight);
          hist_MET_MVS0->Fill(MET, MV_S0, weight);
          hist_MET_MQS0->Fill(MET, MQ_S0, weight);
          hist_MET_gammaS0->Fill(MET, gamma_S0, weight);
          hist_MET_MPTilde->Fill(MET, MPTilde, weight);
          hist_MET_MSTilde->Fill(MET, MSTilde, weight);
          hist_MET_gammaTilde->Fill(MET, gammaTilde, weight);
          hist_MS_MPa->Fill(MS_LEP, MPa, weight);
          hist_MS_MPb->Fill(MS_LEP, MPb, weight);
          hist_MS_CosDecayAnglePaLEP->Fill(MS_LEP, CosDecayAngle_Pa_LEP, weight);
          hist_MS_MVa->Fill(MS_LEP, MVa, weight);
          hist_MS_MVb->Fill(MS_LEP, MVb, weight);
          hist_MS_PTSCM->Fill(MS_LEP, PTS_CM, weight);
          hist_MPa_MPb->Fill(MPa, MPb, weight);
          hist_MPa_MVa->Fill(MPa, MVa, weight);
          hist_CosDecayAnglePaLEP_CosDecayAnglePbLEP->Fill(CosDecayAngle_Pa_LEP, CosDecayAngle_Pb_LEP, weight);
          hist_CosDecayAnglePaLEP_CosDecayAngleSLEP->Fill(CosDecayAngle_Pa_LEP, CosDecayAngle_S_LEP, weight);
          hist_MSS0_MVS0->Fill(MS_S0, MV_S0, weight);
          hist_MSS0_MQS0->Fill(MS_S0, MQ_S0, weight);
          hist_MSS0_gammaS0->Fill(MS_S0, gamma_S0, weight);
          hist_MSS0_MSTilde->Fill(MS_S0, MSTilde, weight);
          hist_MVS0_MQS0->Fill(MV_S0, MQ_S0, weight);
          hist_MVS0_gammaS0->Fill(MV_S0, gamma_S0, weight);
          hist_MQS0_gammaS0->Fill(MQ_S0, gamma_S0, weight);
          hist_MQS0_MPTilde->Fill(MQ_S0, MPTilde, weight);
          hist_gammaS0_gammaTilde->Fill(gamma_S0, gammaTilde, weight);
          hist_MPTilde_MSTilde->Fill(MPTilde, MSTilde, weight);
          hist_MPTilde_gammaTilde->Fill(MPTilde, gammaTilde, weight);
          hist_MSTilde_gammaTilde->Fill(MSTilde, gammaTilde, weight);

          hist_CosDecayAnglePaLEP_MPb->Fill(CosDecayAngle_Pa_LEP, MPb, weight);
          hist_CosDecayAnglePaLEP_MVa->Fill(CosDecayAngle_Pa_LEP, MVa, weight);
          hist_CosDecayAnglePbLEP_MPa->Fill(CosDecayAngle_Pb_LEP, MPa, weight);
          hist_CosDecayAnglePbLEP_MPb->Fill(CosDecayAngle_Pb_LEP, MPb, weight);
          hist_CosDecayAnglePbLEP_MVa->Fill(CosDecayAngle_Pb_LEP, MVa, weight);
          hist_CosDecayAnglePbLEP_MS->Fill(CosDecayAngle_Pb_LEP, MS_LEP, weight);
          hist_MPb_MVa->Fill(MPb, MVa, weight);

          hist_CosDecayAngleVaLEP_CosDecayAnglePbLEP->Fill(CosDecayAngle_Va_LEP, CosDecayAngle_Pb_LEP, weight);
          hist_CosDecayAngleVaLEP_CosDecayAngleSLEP->Fill(CosDecayAngle_Va_LEP, CosDecayAngle_S_LEP, weight);
          hist_MVa_CosDecayAngleVaLEP->Fill(MVa, CosDecayAngle_Va_LEP, weight);
          hist_MS_CosDecayAngleVaLEP->Fill(MS_LEP, CosDecayAngle_Va_LEP, weight);
          hist_MSLEP_MCosDecayAngleP->Fill(MS_LEP, CosDecayAngle_Pa_LEP-CosDecayAngle_Pb_LEP, weight);
          hist_MQLEP_MCosDecayAngleP->Fill(MQ_LEP, CosDecayAngle_Pa_LEP-CosDecayAngle_Pb_LEP, weight);
          hist_MPa_MCosDecayAngleM->Fill(MPa, CosDecayAngle_Pa_LEP*CosDecayAngle_Pb_LEP, weight);
          hist_MPb_MCosDecayAngleM->Fill(MPb, CosDecayAngle_Pa_LEP*CosDecayAngle_Pb_LEP, weight);
          hist_MPa_SCosDecayAngleP->Fill(MPa, CosDecayAngle_Pa_LEP+CosDecayAngle_Pb_LEP, weight);
          hist_MPb_SCosDecayAngleP->Fill(MPb, CosDecayAngle_Pa_LEP+CosDecayAngle_Pb_LEP, weight);

          hist_MCosDecayAngleP->Fill(CosDecayAngle_Pa_LEP-CosDecayAngle_Pb_LEP, weight);
          hist_MCosDecayAngleM->Fill(CosDecayAngle_Pa_LEP*CosDecayAngle_Pb_LEP, weight);
          hist_AbsMCosDecayAngleP->Fill(fabs(CosDecayAngle_Pa_LEP-CosDecayAngle_Pb_LEP), weight);
          hist_AbsMCosDecayAngleM->Fill(fabs(CosDecayAngle_Pa_LEP*CosDecayAngle_Pb_LEP), weight);
          hist_SCosDecayAngleP->Fill(CosDecayAngle_Pa_LEP+CosDecayAngle_Pb_LEP, weight);

          hist_MET_RISRLEP->Fill(MET, RISR_LEP, weight);
          hist_MET_PTISRLEP->Fill(MET, PTISR_LEP, weight);
          hist_RISRLEP_MS->Fill(RISR_LEP, MS_LEP, weight);
          hist_RISRLEP_MPa->Fill(RISR_LEP, MPa, weight);
          hist_RISRLEP_MPb->Fill(RISR_LEP, MPb, weight);
          hist_RISRLEP_MVa->Fill(RISR_LEP, MVa, weight);
          hist_RISRLEP_MVb->Fill(RISR_LEP, MVb, weight);
          hist_RISRLEP_PTSCM->Fill(RISR_LEP, PTS_CM, weight);
          hist_RISRLEP_MSS0->Fill(RISR_LEP, MS_S0, weight);
          hist_RISRLEP_MVS0->Fill(RISR_LEP, MV_S0, weight);
          hist_RISRLEP_MQS0->Fill(RISR_LEP, MQ_S0, weight);
          hist_RISRLEP_gammaS0->Fill(RISR_LEP, gamma_S0, weight);
          hist_RISRLEP_MPTilde->Fill(RISR_LEP, MPTilde, weight);
          hist_RISRLEP_MSTilde->Fill(RISR_LEP, MSTilde, weight);
          hist_RISRLEP_gammaTilde->Fill(RISR_LEP, gammaTilde, weight);

          hist_MINVLEP_RISRLEP->Fill(MINV_LEP, RISR_LEP, weight);
          hist_MINVLEP_MET->Fill(MINV_LEP, MET, weight);
          hist_MINVLEP_PTISRLEP->Fill(MINV_LEP, PTISR_LEP, weight);
          hist_MINVLEP_MS->Fill(MINV_LEP, MS_LEP, weight);
          hist_MINVLEP_MPa->Fill(MINV_LEP, MPa, weight);
          hist_MINVLEP_MPb->Fill(MINV_LEP, MPb, weight);
          hist_MINVLEP_MVa->Fill(MINV_LEP, MVa, weight);
          hist_MINVLEP_MVb->Fill(MINV_LEP, MVb, weight);
          hist_MINVLEP_PTSCM->Fill(MINV_LEP, PTS_CM, weight);
          hist_MINVLEP_MSS0->Fill(MINV_LEP, MS_S0, weight);
          hist_MINVLEP_MVS0->Fill(MINV_LEP, MV_S0, weight);
          hist_MINVLEP_MQS0->Fill(MINV_LEP, MQ_S0, weight);
          hist_MINVLEP_gammaS0->Fill(MINV_LEP, gamma_S0, weight);
          hist_MINVLEP_MPTilde->Fill(MINV_LEP, MPTilde, weight);
          hist_MINVLEP_MSTilde->Fill(MINV_LEP, MSTilde, weight);
          hist_MINVLEP_gammaTilde->Fill(MINV_LEP, gammaTilde, weight);

          hist_RPLEP_MSLEP->Fill(MPa/MPb, MS_LEP, weight);
          hist_RPLEP_MQLEP->Fill(MPa/MPb, MQ_LEP, weight);

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

      for(int i = 1; i < CF_bins; i++)
        hist_CutFlow->GetXaxis()->SetBinLabel(i, vect_str_cutflow_labels[i].c_str());

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

  std::cout << "Finished looping over all events. Making plots..." << std::endl;

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
