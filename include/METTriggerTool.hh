#ifndef METTriggerTool_h
#define METTriggerTool_h

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <Math/DistFunc.h>
#include <TMath.h>

using std::string;

class METTriggerTool {
  
public:
  METTriggerTool();
    
  virtual ~METTriggerTool();

  void BuildMap(const string& csv);
  void BuildMapJSON(const std::string& json);

  int Get_Curve_Index(double HT, int year, bool el, bool mu, bool data);
  double Get_EFF(double MET = 0., double HT = 0., int year = 2017,
		 bool el = false, bool mu = false,
		 bool data = false, int updown = 0);

  double Get_SF(double MET = 0., double HT = 0., int year = 2017,
		bool el = false, bool mu = false,
		bool data = false, int updown = 0);
  double Get_EFF_JSON(double MET, int year, int Nele, int updown = 0);
  double Get_SF_JSON (double MET, int year, int Nele, int updown = 0);
 
private:
  string m_filename;
   std::map<string,std::vector<double> > m_map;
  
  void SetCSV(const string& csv);
  void ParseCSV();

  string Get_Name(double HT, int year, bool el, bool mu, bool data);
  double Get_EFF(string name = "HT-Le600--SingleElectrontrigger-E1--Nele-E1_SingleElectron_2016_Electron",
		 double MET = 150., int updown = 0);

  double Get_Value(std::string& line);
  
};

#endif
