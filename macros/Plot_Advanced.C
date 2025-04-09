#include "PlottingTools.h"

void Plot_Advanced(){

  Long64_t start = gSystem->Now();
  RestFrames::SetStyle();

  string NtuplePath = "/local-scratch/zflowers/NTUPLES/HADD/";
  //string NtuplePath = "root://cmseos.fnal.gov//store/user/lpcsusylep/NTUPLES_v2/";

  // 10 is Summer23BPix
  SampleTool ST(NtuplePath, 10);
  
  ScaleFactorTool SF;

  output_root_file += "Advanced_";
  //g_Label = "PreSelection";
  //g_Label = "2L GG 0B inclS SR";
  //g_Label = "2L notGG 0B inclS SR";
  //g_Label = "2 lepton SR";
  //g_Label = "2 lepton ttbar CR";
  g_Label = "No Cuts";
  //g_Label = "MET > 150";
  //g_Label = "TESTING";
  
  output_root_file += g_Label;
  SanitizeString(output_root_file);
  folder_name = output_root_file;
  output_root_file += ".root";
  if(SavePDF){
    std::cout << "making dir for plots: " << folder_name << std::endl;
    gSystem->Exec(("mkdir -p "+folder_name).c_str());
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

  // hists for holding number of events
  const int EC_bins = vec_samples.size() + 1;
  const int Zbi_bins = map_vsignals.size();
  TH2D* hist_EventCount = new TH2D("EventCount", "EventCount", 22, 0, 22, EC_bins, 0, EC_bins);
  hist_EventCount->GetYaxis()->SetBinLabel(EC_bins, "TOT BKG");
  TH2D* hist_Zbi = new TH2D("Zbi", "Zbi", 22, 0, 22, Zbi_bins, 0, Zbi_bins);

  int vec_samples_index = 0;
  int Zbi_samples_index = 0;
  for (auto p = vec_samples.begin(); p != vec_samples.end(); p++){
    hist_EventCount->GetYaxis()->SetBinLabel(vec_samples_index+1, FP.getTitle(p->first).c_str());

    // vectors to hold 'generic' plotting objects
    vector<TH1*> hists1;
    vector<TH2*> hists2;
    vector<TEfficiency*> effs;

    g_PlotTitle = p->first; 
    int Nsample = p->second.GetN();
    string title = p->first;

    // Declare hists here
    // push_back hists that you want to plot at the end (hists are filled regardless of whether or not you push_back)
    TH1D* hist_MET = new TH1D((title+"_MET").c_str(), (title+"_MET;MET [GeV]").c_str(), g_NX/4, 100., 1000.);
    hists1.push_back(hist_MET);
    hist_stack_MET.push_back(hist_MET); // example pushing hist into vector for stack plot
    TH1D* hist_RISR = new TH1D((title+"_RISR").c_str(), (title+"_RISR;R_{ISR}").c_str(), g_NX/4, 0., 1.);
    hists1.push_back(hist_RISR);
    hist_stack_RISR.push_back(hist_RISR);    
    TH1D* hist_PTISR = new TH1D((title+"_PTISR").c_str(), (title+"_PTISR;p_{T}^{ISR} [GeV]").c_str(), g_NX/4, 0., 1000.);
    hists1.push_back(hist_PTISR);
    hist_stack_PTISR.push_back(hist_PTISR);
    TH1D* hist_gammaPerp = new TH1D((title+"_gammaPerp").c_str(), (title+"_gammaPerp;#gamma_{#perp}").c_str(), g_NX/4, 0., 1.);
    hists1.push_back(hist_gammaPerp);
    hist_stack_gammaPerp.push_back(hist_gammaPerp);
    TH1D* hist_MQperp = new TH1D((title+"_MQperp").c_str(), (title+"_MQperp;M_{#perp}").c_str(), g_NX/4, 0., 150.);
    hists1.push_back(hist_MQperp);
    hist_stack_MQperp.push_back(hist_MQperp);
    TH1D* hist_MSperpCM0 = new TH1D((title+"_MSperpCM0").c_str(), (title+"_MSperpCM0;MS_{#perp CM0}").c_str(), g_NX/4, 0., 600.);
    hists1.push_back(hist_MSperpCM0);
    hist_stack_MSperpCM0.push_back(hist_MSperpCM0);
    TH1D* hist_MQperpCM0 = new TH1D((title+"_MQperpCM0").c_str(), (title+"_MQperpCM0;MQ_{#perp CM0}").c_str(), g_NX/4, 0., 200.);
    hists1.push_back(hist_MQperpCM0);
    hist_stack_MQperpCM0.push_back(hist_MQperpCM0);
    TH1D* hist_gammaPerpCM0 = new TH1D((title+"_gammaPerpCM0").c_str(), (title+"_gammaPerpCM0;#gamma_{#perp CM0}").c_str(), g_NX/4, 0., 1.);
    hists1.push_back(hist_gammaPerpCM0);
    hist_stack_gammaPerpCM0.push_back(hist_gammaPerpCM0);
    TH1D* hist_MSCM0 = new TH1D((title+"_MSCM0").c_str(), (title+"_MSCM0;MS_{CM0}").c_str(), g_NX/4, 0., 600.);
    hists1.push_back(hist_MSCM0);
    hist_stack_MSCM0.push_back(hist_MSCM0);
    TH1D* hist_MQCM0 = new TH1D((title+"_MQCM0").c_str(), (title+"_MQCM0;MQ_{CM0}").c_str(), g_NX/4, 0., 200.);
    hists1.push_back(hist_MQCM0);
    hist_stack_MQCM0.push_back(hist_MQCM0);
    TH1D* hist_gammaCM0 = new TH1D((title+"_gammaCM0").c_str(), (title+"_gammaCM0;#gamma_{CM0}").c_str(), g_NX/4, 0., 1.);
    hists1.push_back(hist_gammaCM0);
    hist_stack_gammaCM0.push_back(hist_gammaCM0);
    TH1D* hist_ML = new TH1D((title+"_ML").c_str(), (title+"_ML;ML").c_str(), g_NX/4, 0., 150.);
    hists1.push_back(hist_ML);
    hist_stack_ML.push_back(hist_ML);
    TH1D* hist_MJ = new TH1D((title+"_MJ").c_str(), (title+"_MJ;MJ").c_str(), g_NX/4, 0., 150.);
    hists1.push_back(hist_MJ);
    hist_stack_MJ.push_back(hist_MJ);

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

          // Apply PreSelection
          
          //if(do_FilterDilepton)
          //  if(SF.DileptonEvent(base))
          //    continue;
          
          // apply trigger to data and FullSim events
          //if(!base->METORtrigger && !is_FastSim)
          //  continue;
          	
          // get variables from root files using base class
          double MET = base->MET;
          double Mperp = base->Mperp;
          double RISR = base->RISR;
          double PTISR = base->PTISR;

          //if(MET < 150.)
          //  continue;

          //if(PTISR < 200.)
          //if(PTISR < 300.) // SR
	  //  continue;

          // Cleaning cuts...
          double dphiCMI = base->dphiCMI;
          double PTCM = base->PTCM;
          double x = fabs(dphiCMI);
          
          //if(PTCM > 200.)
          //  continue;
          //if(PTCM > -500.*sqrt(std::max(0.,-2.777*x*x+1.388*x+0.8264))+575. &&
          //   -2.777*x*x+1.388*x+0.8264 > 0.)
          //  continue;
          //if(PTCM > -500.*sqrt(std::max(0.,-1.5625*x*x+7.8125*x-8.766))+600. &&
          //   -1.5625*x*x+7.8125*x-8.766 > 0.)
          //  continue;
          // End of Cleaning cuts...
            
          double dphiMET_V = base->dphiMET_V;
          //if(fabs(base->dphiMET_V) > acos(-1.)/2.)
          //  continue;
            
          //if(RISR < 0.5 || RISR > 1.0)
          //if(RISR < 0.4 || RISR > 0.7) // CR
          //if(RISR < 0.9)
          //  continue;

          // Get Physics Objects
          int Nlep     = base->Nlep;
          int NjetS    = base->Njet_S;
          int NbjetS   = base->Nbjet_S;
          int NjetISR  = base->Njet_ISR;
          int NbjetISR = base->Nbjet_ISR;

          //if(NbjetISR + NbjetS != 2) continue; // CR
          //if(NbjetISR + NbjetS > 1) continue; // SR

          //if(Nlep != 2) continue;
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
          for(int i = 0; i < base->Nlep_b; i++){
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
          int nGL = 0; // number of Gold leps
          for(int i = 0; i < list_leps.GetN(); i++){
            if(list_leps[i].ID() == kGold) nGL++;
          }
          //if(nGL < 2) skip = true; // SR GG
          //if(nGL == 2) skip = true; // CR notGG
          //if(skip) continue; 
          
          double weight = (base->weight != 0.) ? base->weight : 1.;
          if(!is_data && !is_signal)
            weight *= double(BKG_SKIP);

          // grab vars from ntuple
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

            if (same_flavor && opposite_sign)       EC_X = 1; //flavor_category = "OSSF";
            else if (same_flavor && !opposite_sign) EC_X = 2; //flavor_category = "SSSF";
            else if (!same_flavor && opposite_sign) EC_X = 3; //flavor_category = "OSOF";
            else                                    EC_X = 4; //flavor_category = "SSOF";
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
          else {
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
          hist_EventCount->SetBinContent(0,EC_Y,hist_EventCount->GetBinContent(0,EC_Y)+weight); // normalized to selection
          if(is_bkg) // total SM bkg
            hist_EventCount->SetBinContent(0,EC_bins,hist_EventCount->GetBinContent(0,EC_bins)+weight); // normalized to selection

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

  Long64_t end = gSystem->Now();
  std::cout << "Time to process " << (end-start)/1000.0 << " seconds" << std::endl;
  gApplication->Terminate(0);
} // End of macro
