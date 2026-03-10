#ifndef LepSFTool_h
#define LepSFTool_h

#include "../include/E_bank.hh"
#include "../include/E_bank_fit.hh"

using std::vector;

class LepSFTool {
  
public:
  LepSFTool();
  virtual ~LepSFTool();

  double getMCIdEfficiency(double pt, double eta, int pdg, int year = 0);
  double getDataIdEfficiency(double pt, double eta, int pdg, int year = 0);
  
  double getMCIsoEfficiency(double pt, double eta, int pdg, int year = 0);
  double getDataIsoEfficiency(double pt, double eta, int pdg, int year = 0);
  
  double getMCSipEfficiency(double pt, double eta, int pdg, int year = 0);
  double getDataSipEfficiency(double pt, double eta, int pdg, int year = 0);

  double getMCVLIdEfficiency(double pt, double eta, int pdg, int year = 0);
  double getDataVLIdEfficiency(double pt, double eta, int pdg, int year = 0);
  
  double getMCIdError(double pt, double eta, int pdg, int year = 0);
  double getDataIdError(double pt, double eta, int pdg, int year = 0);
  
  double getMCIsoError(double pt, double eta, int pdg, int year = 0);
  double getDataIsoError(double pt, double eta, int pdg, int year = 0);
  
  double getMCSipError(double pt, double eta, int pdg, int year = 0);
  double getDataSipError(double pt, double eta, int pdg, int year = 0);

  double getMCVLIdError(double pt, double eta, int pdg, int year = 0);
  double getDataVLIdError(double pt, double eta, int pdg, int year = 0);

  double getIdFastSimSF(double pt, double eta, int pdg, int year = 0);
  double getIdFastSimSFerr(double pt, double eta, int pdg, int year = 0);
  
  double getIsoFastSimSF(double pt, double eta, int pdg, int year = 0);
  double getIsoFastSimSFerr(double pt, double eta, int pdg, int year = 0);
 
  double getSipFastSimSF(double pt, double eta, int pdg, int year = 0);
  double getSipFastSimSFerr(double pt, double eta, int pdg, int year = 0);
  
  double getVLIdFastSimSF(double pt, double eta, int pdg, int year = 0);
  double getVLIdFastSimSFerr(double pt, double eta, int pdg, int year = 0);
  
  void BuildMap(const std::string& lepSFfolder); 

private:

  bool b_E_bank;
  E_bank* id_Jmu_MC;
  E_bank* id_Jmu_Data;
  E_bank* id_Zmu_MC;
  E_bank_fit* iso_med_Zmu_MC;
  E_bank_fit* sip_isomed_Zmu_MC;
  E_bank* id_Zmu_Data;
  E_bank_fit* iso_med_Zmu_Data;
  E_bank_fit* sip_isomed_Zmu_Data;
  E_bank_fit* vl_Zmu_Data;
  E_bank_fit* vl_Zmu_MC;
  E_bank* id_Zel_MC;
  E_bank* iso_med_Zel_MC;
  E_bank* sip_isomed_Zel_MC;
  E_bank* id_Zel_Data;
  E_bank* iso_med_Zel_Data;
  E_bank* sip_isomed_Zel_Data;
  E_bank* vl_Zel_Data;
  E_bank* vl_Zel_MC;
  E_bank* id_SimF_el;
  E_bank* iso_med_SimF_el;
  E_bank* sip_isomed_SimF_el;
  E_bank* vl_SimF_el;
  E_bank* id_SimF_mu;
  E_bank* iso_med_SimF_mu;
  E_bank* sip_isomed_SimF_mu;
  E_bank* vl_SimF_mu;

};

#endif






