#ifndef BtagSFTool_h
#define BtagSFTool_h

#include <TEfficiency.h>
#include <TH1D.h>
#include "../include/FormulaBin.hh"

using std::vector;

class BtagSFTool {
  
public:
  BtagSFTool(){
    for(int i = 0; i < 3; i++){
      m_BtagEff[i] = nullptr;
      m_BtagEff_FastSim[i] = nullptr;
      for(int j = 0; j < 3; j++)
	for(int k = 0; k < 3; k++){
	  m_SFs[i][j][k] = nullptr;
	  m_SFs_FastSim[i][j][k] = nullptr;
	}
    }
  }
  
  virtual ~BtagSFTool(){
    for(int i = 0; i < 3; i++){
      if(m_BtagEff[i] != nullptr)
	delete m_BtagEff[i];
      if(m_BtagEff_FastSim[i] != nullptr)
	delete m_BtagEff_FastSim[i];
      for(int j = 0; j < 3; j++)
	for(int k = 0; k < 3; k++){
	  if(m_SFs[i][j][k] != nullptr)
	    delete m_SFs[i][j][k];
	  if(m_SFs_FastSim[i][j][k] != nullptr)
	    delete m_SFs_FastSim[i][j][k];
	}
    }
  }

  double EFF(double pT, int flavor, bool FastSim = false);
  double SF(double pT, int year, int flavor, int updown = 0, bool FastSim = false);
  void BuildMap(const std::string& btagSFfolder);
  void SetEfficiencies(const std::string& rootfile, const std::string& dataset_filetag);
  void SetSFs(const std::string& csvfile, int year, bool FastSim = false);
  
  // flavor_jet
  // 0 - b
  // 1 - c
  // 2 - other

private:
  TEfficiency* m_BtagEff[3];
  TEfficiency* m_BtagEff_FastSim[3];

  FormulaBins* m_SFs[3][3][3]; // [year][flavor][sys]
  FormulaBins* m_SFs_FastSim[3][3][3]; // [year][flavor][sys]

  void ParseCSV(const std::string& csvfile, int iyear, bool FastSim = false);
  std::string popcomma(std::string& line);
};

#endif
