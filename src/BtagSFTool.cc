#include "../include/BtagSFTool.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <TFile.h>

using std::string;



double BtagSFTool::EFF(double pT, int year, int flavor, bool FastSim){
  if(flavor < 0 || flavor > 2)
    return 0.;

  if(FastSim){
    if(year == 2016)
      if(m_BtagEff2016_FastSim[flavor] != nullptr)
	return m_BtagEff2016_FastSim[flavor]->GetEfficiency(m_BtagEff2016[flavor]->FindFixBin(pT));
    if(year == 2017)
      if(m_BtagEff2017_FastSim[flavor] != nullptr)
	return m_BtagEff2017_FastSim[flavor]->GetEfficiency(m_BtagEff2017[flavor]->FindFixBin(pT));
    if(year == 2018)
      if(m_BtagEff2018_FastSim[flavor] != nullptr)
	return m_BtagEff2018_FastSim[flavor]->GetEfficiency(m_BtagEff2018[flavor]->FindFixBin(pT));
  } else {
    if(year == 2016)
      if(m_BtagEff2016[flavor] != nullptr)
	return m_BtagEff2016[flavor]->GetEfficiency(m_BtagEff2016[flavor]->FindFixBin(pT));
    if(year == 2017)
      if(m_BtagEff2017[flavor] != nullptr)
	return m_BtagEff2017[flavor]->GetEfficiency(m_BtagEff2017[flavor]->FindFixBin(pT));
    if(year == 2018)
      if(m_BtagEff2018[flavor] != nullptr)
	return m_BtagEff2018[flavor]->GetEfficiency(m_BtagEff2018[flavor]->FindFixBin(pT));
  }
  
  return 0.;
}

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

void BtagSFTool::BuildMap(const std::string& btagSFfolder, const std::string& proc_rootfile, int year){
  SetEfficiencies(btagSFfolder+"/BtagEff.root", proc_rootfile, year);
  SetSFs(btagSFfolder+"/DeepJet_2016.csv", 2016);
  SetSFs(btagSFfolder+"/DeepJet_2017.csv", 2017);
  SetSFs(btagSFfolder+"/DeepJet_2018.csv", 2018);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_16SL_18_3_2019.csv", 2016, true);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_17SL_18_3_2019.csv", 2017, true);
  SetSFs(btagSFfolder+"/DeepFlav_13TEV_18SL_7_5_2019.csv", 2018, true);
}


void BtagSFTool::SetEfficiencies(const std::string& rootfile, const std::string& proc_rootfile, int year){
  for(int i = 0; i < 3; i++){
    if(m_BtagEff2016[i] != nullptr)
      delete m_BtagEff2016[i];
    if(m_BtagEff2017[i] != nullptr)
      delete m_BtagEff2017[i];
    if(m_BtagEff2018[i] != nullptr)
      delete m_BtagEff2018[i];
    if(m_BtagEff2016_FastSim[i] != nullptr)
      delete m_BtagEff2016_FastSim[i];
    if(m_BtagEff2017_FastSim[i] != nullptr)
      delete m_BtagEff2017_FastSim[i];
    if(m_BtagEff2018_FastSim[i] != nullptr)
      delete m_BtagEff2018_FastSim[i];
  }
//efficiency hack - load the normal root file for other years
//but the specified year/process load the file from the ntuple
std::cout<<"opening original BTAG rootfile: "<<rootfile<<"\n";
  TFile* input = new TFile(rootfile.c_str(),"READ");
  if(!input->IsOpen())
    return;
  
 // root://cmsxrootd.fnal.gov/
  TFile* input2;
  //s..td::string s1 = "root://cmsxrootd.fnal.gov/";
  //root://cmseos.fnal.gov/
  //std::string s2 = proc_rootfile.substr(23);
  //s1 = s1+s2;
  if( proc_rootfile != ""){  
  input2 = TFile::Open(proc_rootfile.c_str(),"READ");
   std::cout<<"opening Histograms: "<<proc_rootfile<<"\n";
  //if(!input2->IsOpen())
  //  return;
  }

  //create the TEffs from ntuple
  //Histograms/hist_btag_flavor0_den
  //Histograms/hist_btag_flavor0_num
 /* commenting this out, in favor of additional SF 
  TEfficiency *flav0, *flav1, *flav2;
  TH1D *f0n, *f0d, *f1n, *f1d, *f2n, *f2d;
  if (proc_rootfile != ""){
  std::cout<<"loading histograms from: "<< proc_rootfile<<"\n";
  input2->cd("Histograms");
  f0n = (TH1D*) input2->Get("Histograms/hist_btag_flavor0_num");
  f0d = (TH1D*) input2->Get("Histograms/hist_btag_flavor0_den");
  f1n = (TH1D*) input2->Get("Histograms/hist_btag_flavor1_num");
  f1d = (TH1D*) input2->Get("Histograms/hist_btag_flavor1_den");
  f2n = (TH1D*) input2->Get("Histograms/hist_btag_flavor2_num");
  f2d = (TH1D*) input2->Get("Histograms/hist_btag_flavor2_den");
 
  flav0 = new TEfficiency(*f0n,*f0d);
  flav1 = new TEfficiency(*f1n,*f1d);
  flav2 = new TEfficiency(*f2n,*f2d);
  std::cout<<" TEffciency created for year "<<year<<" \n";
  
  }*/

   ///fill all data structures normally
  for(int i = 0; i < 3; i++){
    m_BtagEff2016[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2016_flavor%d",i)))->Clone("BEff_2016_flavor%d");
    m_BtagEff2017[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2017_flavor%d",i)))->Clone("BEff_2017_flavor%d");
    m_BtagEff2018[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2018_flavor%d",i)))->Clone("BEff_2018_flavor%d");
    m_BtagEff2016_FastSim[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2016_flavor%d",i)))->Clone("BEff_2016_FastSim_flavor%d");
    m_BtagEff2017_FastSim[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2017_flavor%d",i)))->Clone("BEff_2017_FastSim_flavor%d");
    m_BtagEff2018_FastSim[i] = (TEfficiency*)((TEfficiency*)input->Get(Form("BtagEff_2018_flavor%d",i)))->Clone("BEff_2018_FastSim_flavor%d");
  }

  //overwrite the indicated eff datastructure
  /*
  if( proc_rootfile != ""){
  if( year == 2016 ){
	std::cout<<"overwriting 2016 TEfficiency\n";
	m_BtagEff2016[0] = flav0;
	m_BtagEff2016_FastSim[0] = flav0;
	m_BtagEff2016[1] = flav1;
        m_BtagEff2016_FastSim[1] = flav1;
	m_BtagEff2016[2] = flav2;
        m_BtagEff2016_FastSim[2] = flav2;	
  }
  if( year == 2017 ){
	std::cout<<"overwriting 2017 TEfficiency\n";
	m_BtagEff2017[0] = flav0;
        m_BtagEff2017_FastSim[0] = flav0;
        m_BtagEff2017[1] = flav1;
        m_BtagEff2017_FastSim[1] = flav1;
        m_BtagEff2017[2] = flav2;
        m_BtagEff2017_FastSim[2] = flav2;
  }
  if( year == 2018 ){
	std::cout<<"overwriting 2018 TEfficiency\n";
	m_BtagEff2018[0] = flav0;
        m_BtagEff2018_FastSim[0] = flav0;
        m_BtagEff2018[1] = flav1;
        m_BtagEff2018_FastSim[1] = flav1;
        m_BtagEff2018[2] = flav2;
        m_BtagEff2018_FastSim[2] = flav2;
  }
  }
  */
  input->Close();
  /*if( proc_rootfile != ""){
  input2->Close();
  }*/
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

// double BtagSFTool::GetWeight(int Njet, double* PT_jet, double* flavor_jet, bool* Btag_jet, int year, int updown){
//   if(abs(updown) > 1) return 1.;
  
//   // if(year == 2016)
//   //   if(m_PU2016[updown+1] != nullptr)
//   //     return m_PU2016[updown+1]->weight(NPV);
//   // if(year == 2017)
//   //   if(m_PU2017[updown+1] != nullptr)
//   //     return m_PU2017[updown+1]->weight(NPV);
//   // if(year == 2018)
//   //   if(m_PU2018[updown+1] != nullptr)
//   //     return m_PU2018[updown+1]->weight(NPV);

//   return 1.;
// }
