#ifndef XsecTool_h
#define XsecTool_h

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <nlohmann/json.hpp> // JSON lib

using json = nlohmann::json;

class ProcessInfo {
  public:
    double cross_section;
    std::string accuracy;
    std::string MCM;
    double total_uncertainty;
};

class XsecTool {

public:
  XsecTool();
  virtual ~XsecTool();

  double GetXsec_BKG(const std::string& dataset) const;
  double GetXsec_SMS_code(const std::string& dataset, double MP, unsigned int code, bool Run3) const;
  double GetXsec_SMS(const std::string& dataset, double MP) const; 
  double GetXsec_Cascades(const std::string& dataset, double MP, bool Run3) const;
  void SetFileTag(const std::string& filetag);

  static std::map<std::string,double> m_Label2Xsec_BKG;
  static std::map<std::string,double> InitMap_Xsec_BKG();
  void UpdateXsecFromJSON(const std::string& json_file);

private:

  static std::map<std::string,int> m_N_SMS;
  static std::map<std::string,int> InitMap_N_SMS();
  static std::map<std::string,std::vector<double> > m_Label2Mass_SMS;
  static std::map<std::string,std::vector<double> > InitMap_Mass_SMS();
  static std::map<std::string,std::vector<double> > m_Label2Xsec_SMS;
  static std::map<std::string,std::vector<double> > InitMap_Xsec_SMS();
  static std::map<std::string,std::vector<double> > m_Label2XsecUnc_SMS;
  static std::map<std::string,std::vector<double> > InitMap_XsecUnc_SMS();
  static std::map<std::string,double> m_Label2Xsec_Cascades;
  static std::map<std::string,double> InitMap_Xsec_Cascades();
  json m_jsonData;
  // For resolving conflicts from inputted json file
  std::string m_FileTag;
  ProcessInfo ResolveXsecConflict(const std::string& process_name);
  bool IsMCMMatchFileTag(const std::string& MCM);
  
};

#endif



