#include "../include/ScaleFactorTool.hh"
#include "ReducedBase.hh"

#include "Math/DistFunc.h"

///////////////////////////////////////////
////////// ScaleFactorTool class
///////////////////////////////////////////
void ScaleFactorTool::init_bLambda(){
std::cout<<"init b corrections\n";
bLambdas["DB_2016_flavor1"] = 0.99965;
bLambdas["DB_2016_flavor2"] = 1.00047;
bLambdas["DB_2016_flavor0"] = 0.993772;
bLambdas["ST_2016_flavor1"] = 0.999787;
bLambdas["ST_2016_flavor2"] = 0.999732;
bLambdas["ST_2016_flavor0"] = 0.999793;
bLambdas["TB_2016_flavor1"] = 0.999819;
bLambdas["TB_2016_flavor2"] = 1.00024;
bLambdas["TB_2016_flavor0"] = 0.991269;
bLambdas["ZDY_2016_flavor1"] = 0.999665;
bLambdas["ZDY_2016_flavor2"] = 1.00037;
bLambdas["ZDY_2016_flavor0"] = 0.993308;
bLambdas["QCD_2016_flavor1"] = 0.999819;
bLambdas["QCD_2016_flavor2"] = 1.00038;
bLambdas["QCD_2016_flavor0"] = 0.993585;
bLambdas["Wjets_2016_flavor1"] = 0.999674;
bLambdas["Wjets_2016_flavor2"] = 1.00016;
bLambdas["Wjets_2016_flavor0"] = 0.98821;
bLambdas["DB_2017_flavor1"] = 0.998919;
bLambdas["DB_2017_flavor2"] = 0.999634;
bLambdas["DB_2017_flavor0"] = 0.997803;
bLambdas["ST_2017_flavor1"] = 0.999641;
bLambdas["ST_2017_flavor2"] = 0.999757;
bLambdas["ST_2017_flavor0"] = 0.974068;
bLambdas["TB_2017_flavor1"] = 0.999892;
bLambdas["TB_2017_flavor2"] = 0.999385;
bLambdas["TB_2017_flavor0"] = 0.981191;
bLambdas["ZDY_2017_flavor1"] = 0.999448;
bLambdas["ZDY_2017_flavor2"] = 1.00007;
bLambdas["ZDY_2017_flavor0"] = 0.964503;
bLambdas["QCD_2017_flavor1"] = 1.00112;
bLambdas["QCD_2017_flavor2"] = 1.00085;
bLambdas["QCD_2017_flavor0"] = 0.970989;
bLambdas["Wjets_Pt-100To250_2017_flavor1"] = 0.998264;
bLambdas["Wjets_Pt-100To250_2017_flavor2"] = 0.999942;
bLambdas["Wjets_Pt-100To250_2017_flavor0"] = 0.905242;
bLambdas["Wjets_Pt-250To400_2017_flavor1"] = 0.99799;
bLambdas["Wjets_Pt-250To400_2017_flavor2"] = 0.999162;
bLambdas["Wjets_Pt-250To400_2017_flavor0"] = 0.91331;
bLambdas["Wjets_Pt-400To600_2017_flavor1"] = 0.998393;
bLambdas["Wjets_Pt-400To600_2017_flavor2"] = 0.998939;
bLambdas["Wjets_Pt-400To600_2017_flavor0"] = 0.916633;
bLambdas["Wjets_Pt-50To100_2017_flavor1"] = 0.999072;
bLambdas["Wjets_Pt-50To100_2017_flavor2"] = 1.00043;
bLambdas["Wjets_Pt-50To100_2017_flavor0"] = 0.901885;
bLambdas["Wjets_Pt-600ToInf_2017_flavor1"] = 0.99876;
bLambdas["Wjets_Pt-600ToInf_2017_flavor2"] = 0.999023;
bLambdas["Wjets_Pt-600ToInf_2017_flavor0"] = 0.921131;
bLambdas["DB_2018_flavor1"] = 0.999868;
bLambdas["DB_2018_flavor2"] = 0.999444;
bLambdas["DB_2018_flavor0"] = 1.01425;
bLambdas["ST_2018_flavor1"] = 0.999682;
bLambdas["ST_2018_flavor2"] = 1.00016;
bLambdas["ST_2018_flavor0"] = 0.995576;
bLambdas["TB_2018_flavor1"] = 1.00046;
bLambdas["TB_2018_flavor2"] = 0.998602;
bLambdas["TB_2018_flavor0"] = 0.987426;
bLambdas["ZDY_2018_flavor1"] = 0.999398;
bLambdas["ZDY_2018_flavor2"] = 0.999751;
bLambdas["ZDY_2018_flavor0"] = 0.985088;
bLambdas["QCD_2018_flavor1"] = 1.00006;
bLambdas["QCD_2018_flavor2"] = 1.00175;
bLambdas["QCD_2018_flavor0"] = 0.989323;
bLambdas["Wjets_Pt-100To250_2018_flavor1"] = 0.999363;
bLambdas["Wjets_Pt-100To250_2018_flavor2"] = 0.999747;
bLambdas["Wjets_Pt-100To250_2018_flavor0"] = 0.955501;
bLambdas["Wjets_Pt-250To400_2018_flavor1"] = 0.999234;
bLambdas["Wjets_Pt-250To400_2018_flavor2"] = 0.9997;
bLambdas["Wjets_Pt-250To400_2018_flavor0"] = 0.958674;
bLambdas["Wjets_Pt-400To600_2018_flavor1"] = 0.999037;
bLambdas["Wjets_Pt-400To600_2018_flavor2"] = 0.999727;
bLambdas["Wjets_Pt-400To600_2018_flavor0"] = 0.961257;
bLambdas["Wjets_Pt-50To100_2018_flavor1"] = 0.998915;
bLambdas["Wjets_Pt-50To100_2018_flavor2"] = 0.999821;
bLambdas["Wjets_Pt-50To100_2018_flavor0"] = 0.957529;
bLambdas["Wjets_Pt-600ToInf_2018_flavor1"] = 0.998975;
bLambdas["Wjets_Pt-600ToInf_2018_flavor2"] = 0.998786;
bLambdas["Wjets_Pt-600ToInf_2018_flavor0"] = 0.960181;
}
std::vector<double> ScaleFactorTool::GetbLambdas(std::string proc, std::string filename, int year){
	std::string y_str = std::to_string(year);
	std::string keyf0;
	std::string keyf1; 
	std::string keyf2; 
	if( ((year == 2017) || (year == 2018)) && proc.find("Wjets") != std::string::npos){
		//do by pt pbin
		std::string ptstring{};
		if(filename.find("Pt-100To250") != std::string::npos)
			ptstring="Pt-100To250";
		if(filename.find("Pt-250To400") != std::string::npos)
			ptstring="Pt-250To400";
		if(filename.find("Pt-400To600") != std::string::npos)
			ptstring="Pt-400To600";
		if(filename.find("Pt-50To100") != std::string::npos)
			ptstring="Pt-50To100";
		if(filename.find("Pt-600ToInf") != std::string::npos)
			ptstring="Pt-600ToInf";
		keyf0= proc+"_"+ptstring+"_"+y_str+"_flavor0";
		keyf1= proc+"_"+ptstring+"_"+y_str+"_flavor1";
		keyf2= proc+"_"+ptstring+"_"+y_str+"_flavor2";
	}
	else{
		keyf0 = proc+"_"+y_str+"_flavor0";
		keyf1 = proc+"_"+y_str+"_flavor1";
		keyf2 = proc+"_"+y_str+"_flavor2";
	}
	std::cout<<"getting keys:"<<keyf0<<" "<<keyf1<<" "<<keyf2<<"\n";
	double lambda0, lambda12;
	lambda0 = bLambdas[keyf0];
	lambda12 = (bLambdas[keyf1]+bLambdas[keyf2])/2.;
	std::cout<<"using bLambda factors: "<<lambda0<<" "<<lambda12<<"\n";
	std::vector<double> lams(2);
	lams[0] = lambda0;
	lams[1] = lambda12;
	return lams;
}

double ScaleFactorTool::GetX20BRSF(const string& filename, const string& treename){

  if(filename.find("TChiWZ") == string::npos ||
     filename.find("ZToLL") == string::npos)
    return 1.;
  
   int MP, MC;
  
  sscanf(treename.c_str(), "SMS_%d_%d", &MP, &MC);
  
  if(MP-MC < 5)
    return 1.;

  else
    return 1. + 0.46*ROOT::Math::lognormal_pdf(double((MP-MC)-5)/30., 0., 1.1);
  
}

double ScaleFactorTool::GetMETEff(double MET, int updown){
  if(MET <= 0.)
    return 0.;

  double nom = 0.9921*ROOT::Math::normal_cdf(MET, 77.5, 78.57);
  
  if(updown == 0)
    return nom;
  if(updown > 0)
    return std::min(1., nom*(1.01 + (MET < 300. ? 0.2e-5*(MET-300.)*(MET-300.) : 0.)));
  if(updown < 0)
    return std::min(1., nom*(0.99 + (MET < 300. ? -0.2e-5*(MET-300.)*(MET-300.) : 0.)));

  return 1.;
}

double ScaleFactorTool::GetMETSF(double MET, int updown){
  if(MET <= 0.)
    return 1.;

  double nom = 0.9921*ROOT::Math::normal_cdf(MET, 77.5, 78.57);
  nom /= 0.9974*ROOT::Math::normal_cdf(MET, 79.22, 65.05);

  if(updown == 0)
    return nom;
  if(updown > 0)
    return nom*(1.01 + (MET < 300. ? 0.2e-5*(MET-300.)*(MET-300.) : 0.));
  if(updown < 0)
    return nom*(0.99 + (MET < 300. ? -0.2e-5*(MET-300.)*(MET-300.) : 0.));

  return 1.;
}

bool ScaleFactorTool::DileptonEvent(ReducedBase* base){

  int NWnu = 0;
  for(int i = 0; i < base->genNnu; i++)
    if(fabs(base->genMomPDGID_nu->at(i)) == 24)
      NWnu++;

  if(NWnu >= 2)
    return true;

  int NZlep = 0;
  for(int i = 0; i < base->genNlep; i++)
    if(fabs(base->genMomPDGID_lep->at(i)) == 23)
      NZlep++;

  if(NZlep >= 2)
    return true;

  return false;

}


void ScaleFactorTool::AddBtagFolder(const string& btagfold, const string& proc_rootfile, int year){
  m_BtagSFTool.BuildMap(btagfold, proc_rootfile, year);
}
SleptonFlavor ScaleFactorTool::SleptonEvent(std::vector<int> &genPDGID_susy){

  int Nsmu = 0;
  int Nsel = 0;

  //std::cout << "CRASH?\n" << std::endl;
  //std::cout << "WITHIN genNsusy: " << genPDGID_susy.size() << std::endl;
  for(int i = 0; i < genPDGID_susy.size() ; i++){
    if(abs(genPDGID_susy.at(i)) == 1000013 || abs(genPDGID_susy.at(i)) == 2000013) Nsmu++;
    //std::cout << "HERE?\n" << std::endl;
    //std::cout << "pdgId: \n" << base->genPDGID_susy->at(i) << std::endl;
    if(abs(genPDGID_susy.at(i)) == 1000011 || abs(genPDGID_susy.at(i)) == 2000011) Nsel++;
  }

  if(Nsmu >= 2)
    return kSmu;

  if(Nsel >= 2)
    return kSel;

  return kSmuSel;

}

//
//void ScaleFactorTool::AddBtagFolder(const string& btagfold){
//  m_BtagSFTool.BuildMap(btagfold);
//}  
                        
//double ScaleFactorTool::GetBtagSFWeight(ReducedBase* base, int year, bool fastsim, bool HForLF, int updown, ParticleIDType tag){
double ScaleFactorTool::GetBtagSFWeight( std::vector<double>& PT_jet, int year, bool fastsim, bool HForLF, int updown, ParticleIDType tag, double v2_weight ){ 

  //std::cout<<"IN OTF FUNCTION\n"; 
//  int Njet = base->Njet;
  int iflavor=-1;
  double EFF, SF, PT;

  double probMC   = 1.;
  double probDATA = 1.;

  double w_f0, w_f1, w_f2;

  std::cout<<"SF btag otf called year:"<< year <<" fastsim:"<< fastsim << " HForLF:"<< HForLF <<" updown:"<< updown <<" tag:"<< tag<<"\n";
 //reverse engineer calculations to find flavor
  std::cout<<"PT_jet vec size "<<PT_jet.size()<<"Target v2 weight "<<v2_weight<<"\n"; 
  std::cout<<"The PT of the jets are: \n";
  for( int i=0; i<PT_jet.size(); i++){
	std::cout<<PT_jet.at(i)<<"\n";	
  }


/*  
  for(int i = 0; i < Njet; i++){
    if(abs(base->Flavor_jet->at(i)) == 5)
      iflavor = 0;
    else if(abs(base->Flavor_jet->at(i)) == 4)
      iflavor = 1;
    else
      iflavor = 2;

    std::cout<<"Processing iflavor: "<<iflavor<<"\n";  
 
    if(HForLF && iflavor == 2)
      continue;
    if(!HForLF && iflavor != 2)
      continue;

    PT = base->PT_jet->at(i);
    
    EFF = m_BtagSFTool.EFF(PT, year, iflavor, fastsim);
    SF  = m_BtagSFTool.SF(PT, year, iflavor, updown);

    std::cout << "year: "<< year<< " fastsim: " << fastsim << " HForLF: " << HForLF << " updown: "<< updown <<"\n";
    std::cout << "iflavor:" <<iflavor << " PT:  " << PT << " EFF: " << EFF << " SF: " << SF << " base:  " << base->BtagID_jet->at(i) << std::endl;
    if(fastsim)
      SF *= m_BtagSFTool.SF(PT, year, iflavor, updown, fastsim);

    if(base->BtagID_jet->at(i) >= tag){
      probMC   *= EFF;
      probDATA *= SF*EFF;
    } else {
      probMC   *= (1.-EFF);
      probDATA *= (1.-SF*EFF);
    }
  }

  if(probMC <= 0. || probDATA <= 0.)
    return 1.;
 */


  return probDATA/probMC;

}
