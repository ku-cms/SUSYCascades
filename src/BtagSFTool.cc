#include "../include/BtagSFTool.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <TFile.h>

using std::string;

double BtagSFTool::SF(double pT, int year, int flavor, int updown, bool FastSim){
  double SF = 1.;

  if(abs(updown) > 1)
    return 1.;
  if(flavor < 0 || flavor > 2)
    return 1;

  if(FastSim){
    if(year == 2016)
      if(m_SFs_FastSim[0][flavor][updown+1] != nullptr)
	return m_SFs_FastSim[0][flavor][updown+1]->SF(pT);
    if(year == 2017)
      if(m_SFs_FastSim[1][flavor][updown+1] != nullptr)
	return m_SFs_FastSim[1][flavor][updown+1]->SF(pT);
    if(year == 2018)
      if(m_SFs_FastSim[2][flavor][updown+1] != nullptr)
	return m_SFs_FastSim[2][flavor][updown+1]->SF(pT);
  } else {
    if(year == 2016)
      if(m_SFs[0][flavor][updown+1] != nullptr)
	return m_SFs[0][flavor][updown+1]->SF(pT);
    if(year == 2017)
      if(m_SFs[1][flavor][updown+1] != nullptr)
	return m_SFs[1][flavor][updown+1]->SF(pT);
    if(year == 2018)
      if(m_SFs[2][flavor][updown+1] != nullptr)
	return m_SFs[2][flavor][updown+1]->SF(pT);
  }
    
  return 1;
}

void BtagSFTool::BuildMap(const std::string& btagSFfolder){
  SetSFs(btagSFfolder+"/DeepJet_2016.csv", 2016);
  SetSFs(btagSFfolder+"/DeepJet_2017.csv", 2017);
  SetSFs(btagSFfolder+"/DeepJet_2018.csv", 2018);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_16SL_18_3_2019.csv", 2016, true);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_17SL_18_3_2019.csv", 2017, true);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_18SL_7_5_2019.csv", 2018, true);
}

double BtagSFTool::EFF(double pT, int flavor, bool FastSim){
  if(flavor < 0 || flavor > 2)
    return 0.;
  TEfficiency* eff = FastSim ?
                     m_BtagEff_FastSim[flavor] :
                     m_BtagEff[flavor];
  if(!eff)
    return 0.;
  int bin = eff->FindFixBin(pT);
  return eff->GetEfficiency(bin);
}

void BtagSFTool::SetEfficiencies(const std::string& rootfile,
                                 const std::string& dataset_filetag)
{
  for(int i = 0; i < 3; i++){
    delete m_BtagEff[i]; m_BtagEff[i] = nullptr;
    delete m_BtagEff_FastSim[i]; m_BtagEff_FastSim[i] = nullptr;
  }
  TFile* input = TFile::Open(rootfile.c_str(), "READ");
  if(!input || !input->IsOpen()){
    std::cerr << "Failed to open " << rootfile << "\n";
    return;
  }
  for(int i = 0; i < 3; i++){
    std::string num_name = "hist_btag__" + dataset_filetag +
                           "_flavor" + std::to_string(i) + "_num";
    std::string den_name = "hist_btag__" + dataset_filetag +
                           "_flavor" + std::to_string(i) + "_den";
    TH1D* h_num = (TH1D*)input->Get(("Histograms/" + num_name).c_str());
    TH1D* h_den = (TH1D*)input->Get(("Histograms/" + den_name).c_str());
    TEfficiency* eff = nullptr;
    if(h_num && h_den){
      if(!TEfficiency::CheckConsistency(*h_num, *h_den)){
        std::cerr << "[BtagSFTool] Inconsistent histograms for "
                  << dataset_filetag << " flavor " << i << "\n";
        continue;
      }
      eff = new TEfficiency(*h_num, *h_den);
    } else {
      std::string fallback = Form("BtagEff_flavor%d", i);
      TEfficiency* base =
        (TEfficiency*)input->Get(fallback.c_str());
      if(!base){
        std::cerr << "[BtagSFTool] Could not find: " 
                  << num_name << " or " << den_name <<  "\n"; 
        std::cerr << "[BtagSFTool] Missing fallback: "
                  << fallback << "\n";
        continue;
      }
      std::cout << "[BtagSFTool] Fallback: " << fallback << "\n";
      eff = (TEfficiency*)base->Clone();
    }
    m_BtagEff[i] = eff;
    m_BtagEff_FastSim[i] = (TEfficiency*)eff->Clone(); // set FastSim to nominal for now
  }
  input->Close();
}

void BtagSFTool::SetSFs(const std::string& csvfile, int year, bool FastSim){
  if(year == 2016)
    ParseCSV(csvfile, 0, FastSim);
  if(year == 2017)
    ParseCSV(csvfile, 1, FastSim);
  if(year == 2018)
    ParseCSV(csvfile, 2, FastSim);
}

void BtagSFTool::ParseCSV(const std::string& csvfile, int iyear, bool FastSim){

  if(FastSim){
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++){
	if(m_SFs_FastSim[iyear][i][j] != nullptr)
	  delete m_SFs_FastSim[iyear][i][j];
	m_SFs_FastSim[iyear][i][j] = new FormulaBins();
      }
    }
  } else {
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++){
	if(m_SFs[iyear][i][j] != nullptr)
	  delete m_SFs[iyear][i][j];
	m_SFs[iyear][i][j] = new FormulaBins();
      }
    }
  }

  std::ifstream ifile(csvfile.c_str());
  if(!ifile.is_open()){
    std::cout << "can't open csv file " << csvfile << std::endl;
    return;
  }

  string line, name;
  size_t found, end;

  //discard first line
  getline(ifile,line);
  while(getline(ifile,line)){
    //remove whitespace
    //line.erase(remove(line.begin(), line.end(), ' '), line.end());
    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

    // first digit is WP - keeoing only medium (1) working point for now
    if(stoi(line.substr(0,1)) != 1)
      continue;

    // pop WP
    popcomma(line);

    // measurement type
    name = popcomma(line);
    if(FastSim){
      if(name.find("fastsim") == string::npos)
	continue;
    } else {
      if(name.find("comb") == string::npos &&
	 name.find("incl") == string::npos)
	continue;
    }

    // sys type
    int isys = -2;
    name = popcomma(line);
    if(name == "down")
      isys = -1;
    if(name == "central")
      isys = 0;
    if(name == "up")
      isys = 1;
    if(isys < -1)
      continue;

    // jet flavor
    int iflavor = -1;
    name = popcomma(line);
    if(name == "0")
      iflavor = 0;
    if(name == "1")
      iflavor = 1;
    if(name == "2")
      iflavor = 2;
    if(iflavor < 0)
      continue;

    // eta min/max
    popcomma(line);
    popcomma(line);

    // pT min
    name = popcomma(line);
    int pt_min = stoi(name);

    // pT max
    name = popcomma(line);
    int pt_max = stoi(name);

     // disc min/max
    popcomma(line);
    popcomma(line);

    // formula remains
    found = line.find("\"");
    line.erase(0,found+1);
    found = line.find("\"");
    name = line.substr(0,found);

    //cout << iyear << " " << iflavor << " " << isys << " " << name << endl;

    if(FastSim)
      m_SFs_FastSim[iyear][iflavor][isys+1]->AddBin(pt_min, pt_max, name);
    else
      m_SFs[iyear][iflavor][isys+1]->AddBin(pt_min, pt_max, name);
  }
  
  ifile.close();
}

std::string BtagSFTool::popcomma(std::string& line){
  if(line.find(",") == string::npos)
    return "";
  
  size_t p = line.find(",");
  string ret = line.substr(0,p);
  line.erase(0,p+1);

  return ret;
}
