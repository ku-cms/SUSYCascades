// LepSFTool_FastSim.hh
#ifndef LepSFTool_FastSim_hh
#define LepSFTool_FastSim_hh

#include "correction.h"
#include "Particle.hh"
#include <string>
#include <memory>

class LepSFTool_FastSim {

public:
  LepSFTool_FastSim(const std::string& filename = "");
  virtual ~LepSFTool_FastSim();

  //void SetEra(const std::string& era)   { era_  = era;  }
  void SetYear(const std::string& year) { year_ = year; }
  void BuildMap(const std::string& filename);

  double get_BLP_COL_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_BLP_COL_SF_fs_err(double pt, double eta, int pdg, LepID qual);

  double get_ID_BLP_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_ID_BLP_SF_fs_err(double pt, double eta, int pdg, LepID qual);

  double get_ISO_ID_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_ISO_ID_SF_fs_err(double pt, double eta, int pdg, LepID qual);

  double get_Prompt_ISOID_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_Prompt_ISOID_SF_fs_err(double pt, double eta, int pdg, LepID qual);

  double get_NOT_Prompt_ISOID_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_NOT_Prompt_ISOID_SF_fs_err(double pt, double eta, int pdg, LepID qual);

  double get_NOT_ID_nor_ISO_SF_fs(double pt, double eta, int pdg, LepID qual);
  double get_NOT_ID_nor_ISO_SF_fs_err(double pt, double eta, int pdg, LepID qual);

private:
  std::shared_ptr<const correction::CorrectionSet> cset_;
  //std::string era_  = "";  // must be set via SetEra() before any lookup
  std::string year_ = "";  // must be set via SetYear() before any lookup

  double Lookup_fs(const std::string& sf_key, double pt, double eta) const;
  std::string GetFlav(int pdg) const;

};

#endif