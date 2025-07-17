#include "../include/LepSFTool.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <TFile.h>

using std::string;

LepSFTool::LepSFTool(){
  b_E_bank = false;
  id_Jmu_MC = nullptr;
  id_Jmu_Data = nullptr;
  id_Zmu_MC = nullptr;
  iso_med_Zmu_MC = nullptr;
  sip_isomed_Zmu_MC = nullptr;
  id_Zmu_Data = nullptr;
  iso_med_Zmu_Data = nullptr;
  sip_isomed_Zmu_Data = nullptr;
  vl_Zmu_Data = nullptr;
  vl_Zmu_MC = nullptr;
  id_Zel_MC = nullptr;
  iso_med_Zel_MC = nullptr;
  sip_isomed_Zel_MC = nullptr;
  id_Zel_Data = nullptr;
  iso_med_Zel_Data = nullptr;
  sip_isomed_Zel_Data = nullptr; 
  vl_Zel_Data = nullptr;
  vl_Zel_MC = nullptr;
  id_SimF_el = nullptr;
  iso_med_SimF_el = nullptr;
  sip_isomed_SimF_el = nullptr;
  vl_SimF_el = nullptr;
  id_SimF_mu = nullptr;
  iso_med_SimF_mu = nullptr;
  sip_isomed_SimF_mu = nullptr;
  vl_SimF_mu = nullptr;
}
  
LepSFTool::~LepSFTool(){
  if(!b_E_bank)
    return;
    
  delete id_Jmu_MC;
  delete id_Jmu_Data;
  delete id_Zmu_MC;
  delete iso_med_Zmu_MC;
  delete sip_isomed_Zmu_MC;
  delete id_Zmu_Data;
  delete iso_med_Zmu_Data;
  delete sip_isomed_Zmu_Data;
  delete vl_Zmu_Data;
  delete vl_Zmu_MC;
  delete id_Zel_MC;
  delete iso_med_Zel_MC;
  delete sip_isomed_Zel_MC;
  delete id_Zel_Data;
  delete iso_med_Zel_Data;
  delete sip_isomed_Zel_Data; 
  delete vl_Zel_Data;
  delete vl_Zel_MC;
  delete id_SimF_el;
  delete iso_med_SimF_el;
  delete sip_isomed_SimF_el;
  delete vl_SimF_el;
  delete id_SimF_mu;
  delete iso_med_SimF_mu;
  delete sip_isomed_SimF_mu;
  delete vl_SimF_mu;
}

double LepSFTool::getMCIdEfficiency(double pt, double eta, int pdg, int year){
  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return id_Zel_MC->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    //threshold defaults to 20GeV for jpsi and z
    if(pt < 20.){
      return id_Jmu_MC->getEfficiency(pt,eta,year);
    }
    if(pt > 20.){
      return id_Zmu_MC->getEfficiency(pt,eta,year);
    }
  }
  return -1;
}

double LepSFTool::getDataIdEfficiency(double pt, double eta, int pdg, int year){
  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return id_Zel_Data->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    //threshold defaults to 20GeV for jpsi and z
    if(pt < 20.){
      return id_Jmu_Data->getEfficiency(pt,eta,year);
    }
    if(pt > 20.){
      return id_Zmu_Data->getEfficiency(pt,eta,year);
    }
  }
  return -1;
}
  
double LepSFTool::getMCIsoEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;
  
  if(abs(pdg) == 11){//electron
    return iso_med_Zel_MC->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_Zmu_MC->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getDataIsoEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;
  
  if(abs(pdg) == 11){//electron
    return iso_med_Zel_Data->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_Zmu_Data->getEfficiency(pt,eta,year);
  }
  return -1;
}
  
double LepSFTool::getMCSipEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;
  
  if(abs(pdg) == 11){//electron
    return sip_isomed_Zel_MC->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_Zmu_MC->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getDataSipEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;
  
  if(abs(pdg) == 11){//electron
    return sip_isomed_Zel_Data->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_Zmu_Data->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getMCVLIdEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return vl_Zel_MC->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    return vl_Zmu_MC->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getDataVLIdEfficiency(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return vl_Zel_Data->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    return vl_Zmu_Data->getEfficiency(pt,eta,year);
  }
  return -1;
}
  
double LepSFTool::getMCIdError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return id_Zel_MC->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    //threshold defaults to 20GeV for jpsi and z
    if(pt < 20.){
      return id_Jmu_MC->getError(pt,eta,year);
    }
    if(pt > 20.){
      return id_Zmu_MC->getError(pt,eta,year);
    }
  }
  return -1;
}

double LepSFTool::getDataIdError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return id_Zel_Data->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    //threshold defaults to 20GeV for jpsi and z
    if(pt < 20.){
      return id_Jmu_Data->getError(pt,eta,year);
    }
    if(pt > 20.){
      return id_Zmu_Data->getError(pt,eta,year);
    }
  }
  return -1;
}
  
double LepSFTool::getMCIsoError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return iso_med_Zel_MC->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_Zmu_MC->getError(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getDataIsoError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return iso_med_Zel_Data->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_Zmu_Data->getError(pt,eta,year);
  }
  return -1;
}
  
double LepSFTool::getMCSipError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return sip_isomed_Zel_MC->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_Zmu_MC->getError(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getDataSipError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return sip_isomed_Zel_Data->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_Zmu_Data->getError(pt,eta,year);
  }
  return -1;

}

double LepSFTool::getMCVLIdError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return vl_Zel_MC->getError(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    return vl_Zmu_MC->getError(pt,eta,year);
  }
  return -1;

}

double LepSFTool::getDataVLIdError(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11 || abs(pdg) == 0){//electron
    return vl_Zel_Data->getError(pt,eta,year);
  }
  if(abs(pdg) == 13 || abs(pdg) == 1){//muon
    return vl_Zmu_Data->getError(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getIdFastSimSF(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return id_SimF_el->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return id_SimF_mu->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getIdFastSimSFerr(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return id_SimF_el->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return id_SimF_mu->getError(pt,eta,year);
  }
  return -1;
}
  
double LepSFTool::getIsoFastSimSF(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;
  
  if(abs(pdg) == 11){//electron
    return iso_med_SimF_el->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_SimF_mu->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getIsoFastSimSFerr(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return iso_med_SimF_el->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return iso_med_SimF_mu->getError(pt,eta,year);
  }
  return -1;
}
 
double LepSFTool::getSipFastSimSF(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return sip_isomed_SimF_el->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_SimF_mu->getEfficiency(pt,eta,year);
  }
  return -1;

}

double LepSFTool::getSipFastSimSFerr(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return sip_isomed_SimF_el->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return sip_isomed_SimF_mu->getError(pt,eta,year);
  }
  return -1;
}
  
double LepSFTool::getVLIdFastSimSF(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return vl_SimF_el->getEfficiency(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return vl_SimF_mu->getEfficiency(pt,eta,year);
  }
  return -1;
}

double LepSFTool::getVLIdFastSimSFerr(double pt, double eta, int pdg, int year){
  if(pdg == 0) pdg = 11;
  if(abs(pdg) == 1) pdg = 13;

  if(abs(pdg) == 11){//electron
    return vl_SimF_el->getError(pt,eta,year);
  }
  if(abs(pdg) == 13){//muon
    return vl_SimF_mu->getError(pt,eta,year);
  }
  return -1;
}

void LepSFTool::BuildMap(const std::string& lepSFfolder){
  string path     = lepSFfolder+"/TnP_Muon_Output/";
  string pathEl   = lepSFfolder+"/TnP_Electron_Output/";
  string pathSimF = lepSFfolder+"/FullSim-FastSim_SF/";

  b_E_bank = true;
  
  std::cout<<"Beginning Electron init"<<std::endl;
  id_Zel_MC = new E_bank(2016,
			 pathEl+"TnPZ_susyID_MC2016_Tight.root",
			 pathEl+"TnPZ_susyID_MC2017_Tight.root",
			 pathEl+"TnPZ_susyID_MC2018_tight.root",
			 "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002}, 2018);
  id_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001}, 2018);
  id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.001, 0.005}, 2017);
  id_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
  id_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.003, 0.003}, 2016);
  id_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.002}, 2016);

  iso_med_Zel_MC = new E_bank(2016,
			      pathEl+"TnPZ_susyID_MC2016_Iso.root",
			      pathEl+"TnPZ_susyID_MC2017_Iso.root",
			      pathEl+"TnPZ_susyID_MC2018_Iso.root",
			      "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002}, 2018);
  iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001}, 2018);
  iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.003, 0.003}, 2017);
  iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
  iso_med_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.004, 0.003}, 2016);
  iso_med_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2016);


  sip_isomed_Zel_MC = new E_bank(2016,
				 pathEl+"TnPZ_susyID_MC2016_Sip3D.root",
				 pathEl+"TnPZ_susyID_MC2017_Sip3D.root",
				 pathEl+"TnPZ_susyID_MC2018_Sip3D.root",
				 "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005}, 2018);
  sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2018);
  sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.006, 0.004, 0.007}, 2017);
  sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.0006}, 2017);
  sip_isomed_Zel_MC->applySystematic_ptRange(0., 20., {0.001, 0.007, 0.011}, 2016);
  sip_isomed_Zel_MC->applySystematic_ptRange(20., 999., {0.005, 0.012, 0.003}, 2016);

  id_Zel_Data = new E_bank(2016,
			   pathEl+"TnPZ_susyID_data2016_Tight.root",
			   pathEl+"TnPZ_susyID_data2017_Tight.root",
			   pathEl+"TnPZ_susyID_data2018_tight.root",
			   "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.01, 0.002}, 2018);
  id_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.001, 0.001}, 2018);
  id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.001, 0.005}, 2017);
  id_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
  id_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.003, 0.003}, 2016);
  id_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.002}, 2016);;


  iso_med_Zel_Data = new E_bank(2016,
				pathEl+"TnPZ_susyID_data2016_Iso.root",
				pathEl+"TnPZ_susyID_data2017_Iso.root",
				pathEl+"TnPZ_susyID_data2018_Iso.root",
				"tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.003, 0.009, 0.002}, 2018);
  iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.0003, 0.0007, 0.001}, 2018);
  iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.003, 0.003}, 2017);
  iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.002}, 2017);
  iso_med_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.004, 0.003}, 2016);
  iso_med_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2016);

  sip_isomed_Zel_Data = new E_bank(2016,
				   pathEl+"TnPZ_susyID_data2016_Sip3D.root",
				   pathEl+"TnPZ_susyID_data2017_Sip3D.root",
				   pathEl+"TnPZ_susyID_data2018_Sip3D.root",
				   "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.001, 0.004, 0.005}, 2018);
  sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.001, 0.001}, 2018);
  sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.006, 0.004, 0.007}, 2017);
  sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.002, 0.002, 0.0006}, 2017);
  sip_isomed_Zel_Data->applySystematic_ptRange(0., 20., {0.001, 0.007, 0.011}, 2016);
  sip_isomed_Zel_Data->applySystematic_ptRange(20., 999., {0.005, 0.012, 0.003}, 2016);

  vl_Zel_MC = new E_bank(2016,
			 pathEl+"TnPZ_susyID_MC2016_veryLoose.root",
			 pathEl+"TnPZ_susyID_MC2017_veryLoose.root",
			 pathEl+"TnPZ_susyID_MC2018_veryLoose.root",
			 "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  vl_Zel_MC->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01}, 2018);
  vl_Zel_MC->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003}, 2018);
  vl_Zel_MC->applySystematic_ptRange(0., 20., {0.002, 0.007, 0.03}, 2017);
  vl_Zel_MC->applySystematic_ptRange(20., 999., {0.003, 0.0001, 0.0007}, 2017);
  vl_Zel_MC->applySystematic_ptRange(0., 20., {0.02, 0.04, 0.02}, 2016);
  vl_Zel_MC->applySystematic_ptRange(20., 999., {0.003, 0.005, 0.003}, 2016);

  vl_Zel_Data = new E_bank(2016,
			   pathEl+"TnPZ_susyID_data2016_veryLoose.root",
			   pathEl+"TnPZ_susyID_data2017_veryLoose.root",
			   pathEl+"TnPZ_susyID_data2018_veryLoose.root",
			   "tpTree/TightSUSY_pt_eta/fit_eff_plots/");
  vl_Zel_Data->applySystematic_ptRange(0., 20., {0.01, 0.02, 0.01}, 2018);
  vl_Zel_Data->applySystematic_ptRange(20., 999., {0.001, 0.002, 0.003}, 2018);
  vl_Zel_Data->applySystematic_ptRange(0., 20., {0.002, 0.007, 0.03}, 2017);
  vl_Zel_Data->applySystematic_ptRange(20., 999., {0.003, 0.0001, 0.0007}, 2017);
  vl_Zel_Data->applySystematic_ptRange(0., 20., {0.02, 0.04, 0.02}, 2016);
  vl_Zel_Data->applySystematic_ptRange(20., 999., {0.003, 0.005, 0.003}, 2016);

  std::cout<<"Beginning Muon init"<<std::endl;
  /////MUON ID/////////
  id_Jmu_Data = new E_bank(2016,
			   path+"TnPJ_MuonID_data2016_Mu7p5Tk2_veryLoose__Medium_E.root",
			   path+"TnPJ_MuonID_data2017_Mu7p5Tk2_veryLoose__Medium_E.root",
			   path+"TnPJ_MuonID_data2018_Mu7p5Tk2_veryLoose__Medium_E.root",
			   "tpTree/Medium_pt_eta/fit_eff_plots/");
  //id_Jmu_Data->applySystematics(0.0151);	
  id_Jmu_Data->applySystematics(std::vector<double>{0.001,0.001}, 2018);
  id_Jmu_Data->applySystematics(std::vector<double>{0.001,0.001}, 2017);
  id_Jmu_Data->applySystematics(std::vector<double>{0.002,0.001}, 2016);

  id_Zmu_Data = new E_bank(2016,
			   path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__Medium_A_pt1.root",
			   path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
			   path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
			   "tpTree/Medium_pt_eta/fit_eff_plots/");
  //id_Zmu_Data->applySystematics(0.0054);
  id_Zmu_Data->applySystematics(std::vector<double>{0.0007,0.0003}, 2018);
  id_Zmu_Data->applySystematics(std::vector<double>{0.001,0.0003}, 2017);
  id_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002}, 2016);
    
  id_Jmu_MC = new E_bank(2016,
			 path+"TnPJ_MuonID_mc2016_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
			 path+"TnPJ_MuonID_mc2017_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
			 path+"TnPJ_MuonID_mc2018_weight_Mu7p5Tk2_veryLoose__Medium_E.root",
			 "tpTree/Medium_pt_eta/fit_eff_plots/");
  //id_Jmu_MC->applySystematics(0.0151);	
  id_Jmu_MC->applySystematics(std::vector<double>{0.001,0.001}, 2018);
  id_Jmu_MC->applySystematics(std::vector<double>{0.001,0.001}, 2017);
  id_Jmu_MC->applySystematics(std::vector<double>{0.002,0.001}, 2016);

  id_Zmu_MC = new E_bank(2016,
			 path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__Medium_A_pt1.root",
			 path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
			 path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__Medium_A_pt1.root",
			 "tpTree/Medium_pt_eta/fit_eff_plots/");
  //id_Zmu_MC->applySystematics(0.0054);
  id_Zmu_MC->applySystematics(std::vector<double>{0.0007,0.0003}, 2018);
  id_Zmu_MC->applySystematics(std::vector<double>{0.001,0.0003}, 2017);
  id_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002}, 2016);

  /////MUON ISO///////////

  iso_med_Zmu_Data = new E_bank_fit(2016,
				    path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				    path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				    path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				    "tpTree/Medium_pt_eta/fit_eff_plots/");
  //iso_med_Zmu_Data->applySystematics(0.0057);	
  iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.004},2018);
  iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.002},2017);
  iso_med_Zmu_Data->applySystematics(std::vector<double>{0.007,0.004},2016);
  iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.007,0.002},2018);
  iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.007,0.004},2017);
  iso_med_Zmu_Data->setSystematicsLow(std::vector<double>{0.005,0.0009},2016);
  iso_med_Zmu_Data->doLowPtFit(3.,50.,20.,
			       path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
			       path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
			       path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
			       "tpTree/Medium_pt_eta/fit_eff_plots/");
  //iso_med_Zmu_Data->applySystematicsLow(std::vector<double>{0.008,0.004});
		
  iso_med_Zmu_MC = new E_bank_fit(2016,
				  path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt1.root",
				  path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				  path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt1.root",
				  "tpTree/Medium_pt_eta/fit_eff_plots/");
  //iso_med_Zmu_MC->applySystematics(0.0057);
  iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.004},2018);
  iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.002},2017);
  iso_med_Zmu_MC->applySystematics(std::vector<double>{0.007,0.004},2016);
  iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.007,0.002},2018);
  iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.007,0.004},2017);
  iso_med_Zmu_MC->setSystematicsLow(std::vector<double>{0.005,0.0009},2016);
  iso_med_Zmu_MC->doLowPtFit(3.,50.,20.,
			     path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__ISO_MED_A_pt3.root",
			     path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
			     path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__ISO_MED_A_pt3.root",
			     "tpTree/Medium_pt_eta/fit_eff_plots/");
  //iso_med_Zmu_MC->applySystematicsLow(std::vector<double>{0.008,0.004});	


  //////MUON SIP//////////
  sip_isomed_Zmu_Data = new E_bank_fit(2016,
				       path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				       path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				       path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				       "tpTree/Medium_pt_eta/fit_eff_plots/");
  //sip_isomed_Zmu_Data->applySystematics(0.0057);
  sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.0004,0.001},2018);
  sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002},2017);
  sip_isomed_Zmu_Data->applySystematics(std::vector<double>{0.001,0.002},2016);	
  sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.003,0.002},2018);
  sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.005,0.003},2017);
  sip_isomed_Zmu_Data->setSystematicsLow(std::vector<double>{0.002,0.002},2016);
  sip_isomed_Zmu_Data->doLowPtFit(6.,50.,20.,
				  path+"TnP_MuonID_data2016_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				  path+"TnP_MuonID_data2017_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				  path+"TnP_MuonID_data2018_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				  "tpTree/Medium_pt_eta/fit_eff_plots/");
  //sip_isomed_Zmu_Data->applySystematicsLow(std::vector<double>{0.005,0.003});
				
  sip_isomed_Zmu_MC = new E_bank_fit(2016,
				     path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt1.root",
				     path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				     path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt1.root",
				     "tpTree/Medium_pt_eta/fit_eff_plots/");
  //sip_isomed_Zmu_MC->applySystematics(0.0057);
  sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.0004,0.001},2018);
  sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002},2017);
  sip_isomed_Zmu_MC->applySystematics(std::vector<double>{0.001,0.002},2016);	
  sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.003,0.002},2018);
  sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.005,0.003},2017);
  sip_isomed_Zmu_MC->setSystematicsLow(std::vector<double>{0.002,0.002},2016);
  sip_isomed_Zmu_MC->doLowPtFit(6.,50.,20.,
				path+"TnP_MuonID_mc2016_weight_IsoTkMu22_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_veryLoose__SIP_ISOMED_A_pt3.root",
				"tpTree/Medium_pt_eta/fit_eff_plots/");
  //sip_isomed_Zmu_MC->applySystematicsLow(std::vector<double>{0.005,0.003});

  /////////MUON VERY LOOSE//////////
  vl_Zmu_Data = new E_bank_fit(2016,
			       path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt1.root",
			       path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt1.root",
			       path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt1.root",
			       "tpTree/Medium_pt_eta/fit_eff_plots/");
  //vl_Zmu_Data->applySystematics(0.0054);
  vl_Zmu_Data->applySystematics(std::vector<double>{0.0002, 0.0001}, 2018);
  vl_Zmu_Data->applySystematics(std::vector<double>{0.001, 0.001}, 2017);
  vl_Zmu_Data->applySystematics(std::vector<double>{0.0004, 0.002}, 2016);	
  vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.0008,0.003}, 2018);
  vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.001,0.0003}, 2017);
  vl_Zmu_Data->setSystematicsLow(std::vector<double>{0.002,0.003}, 2016);
  vl_Zmu_Data->doLowPtFit(6.,50.,20.,
			  path+"TnP_MuonID_data2016_IsoTkMu22_VeryLoose_A_pt3.root",
			  path+"TnP_MuonID_data2017_isoMu24eta2p1_VeryLoose_A_pt3.root",
			  path+"TnP_MuonID_data2018_isoMu24eta2p1_VeryLoose_A_pt3.root",
			  "tpTree/Medium_pt_eta/fit_eff_plots/");
  //vl_Zmu_Data->applySystematicsLow(std::vector<double>{0.001,0.001});

  vl_Zmu_MC = new E_bank_fit(2016,
			     path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt1.root",
			     path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
			     path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt1.root",
			     "tpTree/Medium_pt_eta/fit_eff_plots/");
  //vl_Zmu_MC->applySystematics(0.0054);
  vl_Zmu_MC->applySystematics(std::vector<double>{0.0002, 0.0001}, 2018);
  vl_Zmu_MC->applySystematics(std::vector<double>{0.001, 0.001}, 2017);
  vl_Zmu_MC->applySystematics(std::vector<double>{0.0004, 0.002}, 2016);	
  vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.0008,0.003}, 2018);
  vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.001,0.0003}, 2017);
  vl_Zmu_MC->setSystematicsLow(std::vector<double>{0.002,0.003}, 2016);
  vl_Zmu_MC->doLowPtFit(6.,50.,20.,
			path+"TnP_MuonID_mc2016_weight_IsoTkMu22_VeryLoose_A_pt3.root",
			path+"TnP_MuonID_mc2017_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
			path+"TnP_MuonID_mc2018_weight_isoMu24eta2p1_VeryLoose_A_pt3.root",
			"tpTree/Medium_pt_eta/fit_eff_plots/");
  //vl_Zmu_MC->applySystematics(std::vector<double>{0.001,0.001});

  /////fullsim fast sim SF

  id_SimF_el =  new E_bank(2016,
			   pathSimF+"TTJets_FastVsFull_Summer16_el_sf_ID.root",
			   pathSimF+"TTJets_FastVsFull_Fall17_el_sf_ID.root",
			   pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_ID.root",
			   "");


  iso_med_SimF_el = new E_bank(2016,
			       pathSimF+"TTJets_FastVsFull_Summer16_el_sf_iso.root",
			       pathSimF+"TTJets_FastVsFull_Fall17_el_sf_iso.root",
			       pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_iso.root",
			       "");

  sip_isomed_SimF_el = new E_bank(2016,
				  pathSimF+"TTJets_FastVsFull_Summer16_el_sf_sip3D.root",
				  pathSimF+"TTJets_FastVsFull_Fall17_el_sf_sip3D.root",
				  pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_sip3D.root",
				  ""); 
  vl_SimF_el = new E_bank(2016,
			  pathSimF+"TTJets_FastVsFull_Summer16_el_sf_vl.root",
			  pathSimF+"TTJets_FastVsFull_Fall17_el_sf_vl.root",
			  pathSimF+"TTJets_FastVsFull_Autumn18_el_sf_vl.root",
			  "");

  id_SimF_mu =  new E_bank(2016,
			   pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_ID.root",
			   pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_ID.root",
			   pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_ID.root",
			   "");

  iso_med_SimF_mu = new E_bank(2016,
			       pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_iso.root",
			       pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_iso.root",
			       pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_iso.root",
			       "");

  sip_isomed_SimF_mu = new E_bank(2016,
				  pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_sip3D.root",
				  pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_sip3D.root",
				  pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_sip3D.root",
				  "");
  vl_SimF_mu = new E_bank(2016,
			  pathSimF+"TTJets_FastVsFull_Summer16_mu_sf_vl.root",
			  pathSimF+"TTJets_FastVsFull_Fall17_mu_sf_vl.root",
			  pathSimF+"TTJets_FastVsFull_Autumn18_mu_sf_vl.root",
			  "");
}




