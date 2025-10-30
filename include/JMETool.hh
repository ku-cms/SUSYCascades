#ifndef JMETool_h
#define JMETool_h

#include <memory>
#include "correction.h"
#include "FormulaBin.hh"
#include "JecApplication.hh"
#include "JecConfigReader.hh"
#include "JvmApplication.hh"
#include "JvmConfigReader.hh"

using std::string;

class JMETool {
  
public:
  JMETool(){}
  
  virtual ~JMETool();

  void BuildMap(const std::string& JMEfolder);
  bool IsJVMEnabled() const { return m_JvmEnabled; }

  void ParseJESUncertainty(const string& input, int year, const string& prefix = "");
  void ParseJEC(const string& input, int year, const string& prefix = "");
  
  double GetJESFactor(int year, const string& name, double pT, double Eta, double A = 0., double rho = 0.) const;

  double GetJESFactorCLIB(const std::string& year, const bool& isData, const int& runNumber, const double& pt, const double& eta, const double& phi, const double& area, const double& rho);
  double GetJESFactorSystCLIB(const std::string& jmeYearKey, const double& ptAfterJes, const double& eta, const double& phi, const double& area, const double& rho, const std::string& clibSystKey, const std::string& var /* "Up" or "Down" */) const;
  double GetJERScaleFactor(const std::string& jmeYearKey, const JecApplication::JesInputs& jAfter, const JecApplication::JerInputs& jerIn, const std::string& var, const bool& isData) const;
  double GetJERResolution(const std::string& jmeYearKey, double eta, double pt) const;

  TLorentzVector GetCorrectedMET(const std::string& jmeYearKey, bool isData, const std::optional<double>& runNumber, double metRawPt, double metRawPhi, const std::vector<JecApplication::JetForMet>& jetsForMet, const std::vector<JecApplication::JerInputs>& jersForMet, double rho, const JecApplication::SystematicOptions& systOpts, bool isDebug = false) const;
  void SetupJVM(const JvmConfigReader::Jvm& jvm);
  bool JetInJVM(const double& pt, const double& eta, const double& phi, const int& jetID, const double& ChEmEF, const double& NeEmEF) const;

  // JER based on https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
  // with values taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution 
  // and https://github.com/cms-jet/JRDatabase/tree/master/textFiles
  void BuildJERMap(const std::string& JMEfolder);

  void ParseJERSFUncertainty(const string& input, int year);
  void ParseJERSF(const string& input, int year, const string& prefix = "");
  void ParseJER(const string& input, int year);
  double GetJERFactor(int year, double pT, double Eta, double rho) const;
  double GetJERSFFactor(int year, double Eta, int updown = 0) const;
  int GetPairIndex(int year, double Eta) const;

private:
  mutable std::map<string, FormulaBinsBins*> m_Factors[3]; // [year]
  mutable std::map<int, FormulaBinsBins*> m_JERFactors[3]; // [year]
  mutable std::map<int, std::pair<double,double>> m_JERSFEtaBins[3]; // <etaIndex, <etaMin,etaMax>> [year]
  mutable std::map<int, std::pair<double,double>> m_JERSFFactors[3]; // <etaBin, <SF,unc>> [year]

  double popdouble(std::string& line);
  std::string popstring(std::string& line);

  bool m_JvmEnabled{false};
  std::unique_ptr<JvmApplication::VetoChecker> m_JvmChecker{nullptr};
};

#endif






