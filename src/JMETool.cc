#include "../include/JMETool.hh"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using std::string;
using std::cout;
using std::endl;

JMETool::~JMETool(){
  for(int i = 0; i < 3; i++){
    std::map<string, FormulaBinsBins*>::iterator it = m_Factors[i].begin();
    while(it != m_Factors[i].end()){
      delete it->second;
      it++;
    }
    std::map<int, FormulaBinsBins*>::iterator it_JER = m_JERFactors[i].begin();
    while(it_JER != m_JERFactors[i].end()){
      delete it_JER->second;
      it_JER++;
    }
  }
}

void JMETool::SetupJVM(const JvmConfigReader::Jvm& jvm) {
    if (jvm.use) {
        m_JvmChecker = std::make_unique<JvmApplication::VetoChecker>(jvm.ref, jvm.key);
        m_JvmEnabled = true;
    } else {
        m_JvmEnabled = false;
    }
}

double JMETool::GetJESFactor(int year, const string& name, double pT, double Eta, double A, double rho) const {
  if(year < 2016 || year > 2018)
    return 0;
  
  if(m_Factors[year-2016].count(name) == 0)
    return 0;
  
  return m_Factors[year-2016][name]->SF(Eta, pT, A, rho);
}

double JMETool::GetJERFactor(int year, double pT, double Eta, double rho) const {
  if(year < 2016 || year > 2018)
    return 1.;
  if(fabs(Eta) > 4.69) return 1.;

  if(year == 2016 && rho > 40.9) rho = 40.89;
  if(year == 2017 && rho > 42.52) rho = 42.51;
  if(year == 2018 && rho > 90.) rho = 89.99;
  
  if(m_JERFactors[year-2016].count(year) == 0)
    return 1.;
  
  return m_JERFactors[year-2016][year]->JERSF(Eta, pT, rho); 
}

double JMETool::GetJESFactorCLIB(const std::string& year, const bool& isData, const std::string& DataEra, const int& runNumber, const double& pt, const double& eta, const double& phi, const double& area, const double& rho) {
  auto tagsAK4 = JecConfigReader::getTagsAK4(year, isData, isData ? std::optional<std::string>(DataEra) : std::nullopt);
  JecConfigReader::CorrectionRefs refs(tagsAK4);

  JecApplication::EvalContext ctx{.year = year, .refs = refs, .isData = isData, .runNumber = isData ? std::optional<double>(runNumber) : std::nullopt, .isDebug = false};

  JecApplication::JesInputs j{pt, eta, phi, area, rho, 0.0};
  return JecApplication::getJesCorrectionNom(ctx, j);
}

double JMETool::GetJESFactorSystCLIB(const std::string& jmeYearKey, const double& ptAfterJes, const double& eta, const double& phi, const double& area, const double& rho, const std::string& clibSystKey, const bool& isData, const std::string& DataEra, const std::string& var /* "Up" or "Down" */) const{
  // Build tags & refs using their reader (no era for MC here)
  auto tags = JecConfigReader::getTagsAK4(jmeYearKey, isData, isData ? std::optional<std::string>(DataEra) : std::nullopt);
  JecConfigReader::CorrectionRefs refs(tags);

  // Forward to helper - var must be "Up" or "Down"
  return JecApplication::getJesCorrectionSyst(refs, clibSystKey, var, eta, ptAfterJes, /*isDebug=*/false);
}

// Get multiplicative JER factor using JecApplication::getJerCorrectionNomOrSyst
double JMETool::GetJERScaleFactor(const std::string& jmeYearKey, const JecApplication::JesInputs& jAfter, const JecApplication::JerInputs& jerIn, const std::string& var, const bool& isData, const std::string& DataEra) const {
  // Build tags & refs for AK4 (MC)
  const auto tags = JecConfigReader::getTagsAK4(jmeYearKey, isData, isData ? std::optional<std::string>(DataEra) : std::nullopt);
  JecConfigReader::CorrectionRefs refs(tags);

  // Build EvalContext
  JecApplication::EvalContext ctx{.year = jmeYearKey, .refs = refs, .isData = isData, .runNumber = std::nullopt, .isDebug = false};

  if (jerIn.run > 0 && JecApplication::requiresRunBasedResidual(jmeYearKey)) {
      ctx.runNumber = static_cast<double>(jerIn.run);
  }

  // Build systematics options
  JecApplication::SystematicOptions systOpts{};
  if (var == "nom" || var == "up" || var == "down")
      systOpts.jerVar = var;
  else
      systOpts.jerVar = "nom";

  // Call the helper which internally handles gen-matching and smearing behavior
  const double sJer = JecApplication::getJerCorrectionNomOrSyst(ctx, jAfter, jerIn, systOpts);
  return sJer;
}

// Query pt-resolution from the CLIB ptResolution correction (tagNamePtResolution)
double JMETool::GetJERResolution(const std::string& jmeYearKey, double eta, double pt, const bool& isData, const std::string& DataEra) const {
  const auto tags = JecConfigReader::getTagsAK4(jmeYearKey, isData, isData ? std::optional<std::string>(DataEra) : std::nullopt);
  JecConfigReader::CorrectionRefs refs(tags);

  // tagNamePtResolution is present in Tags (see JecConfigAK4.json)
  const std::string &ptResName = tags.tagNamePtResolution;

  // Use safeGet helper that the JecConfigReader exposes (same pattern as JES syst)
  auto corr = JecConfigReader::safeGet(refs.cs, ptResName);
  if (!corr) return 0.0;

  // evaluate expects vector arguments depending on the correction; typical signature: {eta, pt}
  const double res = corr->evaluate({ eta, pt });
  return res;
}

TLorentzVector JMETool::GetCorrectedMET(const std::string& jmeYearKey, bool isData, const std::string& DataEra, const std::optional<double>& runNumber, double metRawPt, double metRawPhi, const std::vector<JecApplication::JetForMet>& jetsForMet, const std::vector<JecApplication::JerInputs>& jersForMet, double rho, const JecApplication::SystematicOptions& systOpts, bool isDebug) const {
  // Build tags & refs for AK4 (allow era/data selection if needed)
  const auto tags = JecConfigReader::getTagsAK4(jmeYearKey, isData, isData ? std::optional<std::string>(DataEra) : std::nullopt);
  JecConfigReader::CorrectionRefs refs(tags);
  // Build EvalContext
  JecApplication::EvalContext ctx{.year = jmeYearKey, .refs = refs, .isData = isData, .runNumber = std::nullopt, .isDebug = isDebug};
  if (isData && runNumber.has_value() && JecApplication::requiresRunBasedResidual(jmeYearKey)) ctx.runNumber = runNumber;
  // Build MetInputs and call the JecApplication helper
  JecApplication::MetInputs metIn{ metRawPt, metRawPhi };
  TLorentzVector p4 = JecApplication::getCorrectedMet(ctx, metIn, jetsForMet, jersForMet, rho, systOpts);
  return p4;
}

bool JMETool::JetInJVM(const double& pt, const double& eta, const double& phi, const int& jetID, const double& ChEmEF, const double& NeEmEF) const{
  if (!m_JvmEnabled || !m_JvmChecker) return false;
  return m_JvmChecker->checkJetInVetoRegion(eta, phi, pt, jetID, ChEmEF, NeEmEF);
}

int JMETool::GetPairIndex(int year, double Eta) const {
  for(std::map<int, std::pair<double,double>>::iterator it = m_JERSFEtaBins[year-2016].begin(); it != m_JERSFEtaBins[year-2016].end(); it++)
  {
   if((Eta - m_JERSFEtaBins[year-2016][it->first].second)*(Eta - m_JERSFEtaBins[year-2016][it->first].first) <= 0.)
    return it->first;
  }
  return -1;
}

double JMETool::GetJERSFFactor(int year, double Eta, int updown) const {
  if(year < 2016 || year > 2018)
    return 1.;
  if(fabs(Eta) > 4.69) return 1.;
  
  double binEtaMin = 0.;
  double binEtaMax = 0.;
  double SF = 1.;
  double SFunc = 0.;
  // function to get eta min and max
  int eta_Index = GetPairIndex(year, Eta); 
  if(eta_Index < 0) std::cout << "Could not find eta bin for eta = " << Eta << std::endl;
  // function to get SF and SF unc (SF-SFdown)
  SF = m_JERSFFactors[year-2016][eta_Index].first;
  SFunc = m_JERSFFactors[year-2016][eta_Index].second;

  // up
  if(updown > 0)
   return SF+SFunc;
  // down
  else if(updown < 0)
   return SF-SFunc;
  // nominal
  else
   return SF;
  
  return 1.;
}

void JMETool::BuildJERMap(const std::string& JMEfolder){

  string JME[3]; // [year]

  JME[0] = "Summer16_25nsV1b_MC";
  JME[1] = "Fall17_V3_MC";
  JME[2] = "Autumn18_V6_MC";
  
  string post = "AK4PFchs.txt";

  cout << "Building JER maps for all three years" << endl;
  
  for(int y = 0; y < 3; y++){
    ParseJERSFUncertainty(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_SF_"+post, 2016+y);
    ParseJER(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_PtResolution_"+post, 2016+y);
  }
  
}

void JMETool::ParseJERSFUncertainty(const string& input, int year){
  if(year < 2016 || year > 2018)
    return;

  std::ifstream ifile(input.c_str());
   if(!ifile.is_open()){
    std::cout << "can't open JER SF file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   double etamin, etamax, factor, unc;
   int eta_Index = 0;
   while(getline(ifile,line)){
     // discard comment line
     if(line.find("{") != string::npos)
       continue;
       
     etamin = popdouble(line);
     etamax = popdouble(line);
     std::pair<double,double> eta_pair(etamin,etamax);
     m_JERSFEtaBins[year-2016][eta_Index] = eta_pair;

     popdouble(line); //not used, just need to remove this unused entry
     factor = popdouble(line);
     unc = factor - popdouble(line);
     std::pair<double,double> factor_pair(factor,unc);
     m_JERSFFactors[year-2016][eta_Index] = factor_pair;
     eta_Index++;
     
   }
   ifile.close();
}

void JMETool::BuildMap(const std::string& JMEfolder){
  string JME[3]; // [year]
  string JME_FS[3]; // [year];

  JME[0] = "Summer16_07Aug2017_V11_MC";
  JME[1] = "Fall17_17Nov2017_V32_MC";
  JME[2] = "Autumn18_V19_MC";
  
  JME_FS[0] = "Summer16_FastSimV1_MC";
  JME_FS[1] = "Fall17_FastSimV1_MC";
  JME_FS[2] = "Autumn18_FastSimV1_MC";

  string post = "AK4PFchs.txt";

  cout << "Building JES maps for all three years" << endl;
  
  for(int y = 0; y < 3; y++){
    // ParseJESUncertainty(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_Uncertainty_"+post, 2016+y);
    ParseJESUncertainty(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_UncertaintySources_"+post, 2016+y);
    ParseJESUncertainty(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_Uncertainty_"+post, 2016+y, "FS");

    ParseJEC(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_L1FastJet_"+post, 2016+y);
    ParseJEC(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_L1FastJet_"+post, 2016+y, "FS");
  
    ParseJEC(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_L2Relative_"+post, 2016+y);
    ParseJEC(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_L2Relative_"+post, 2016+y, "FS");
    //cout << "...Done" << endl;
  }
}

void JMETool::ParseJER(const string& input, int year){
  if(year < 2016 || year > 2018)
    return;
  
  std::ifstream ifile(input.c_str());
   if(!ifile.is_open()){
    std::cout << "can't open JER file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   // first line with variable and equation info
   getline(ifile,line);
   
   found = line.find("sqrt");
   line.erase(0,found-1);

   string formula = popstring(line);
   int Nvar = 1;
   int Nparam = 4;
  
   // get the first line to get parameter numbers
   double etamin, etamax;
   double emi, ema;
   double rhomin, rhomax;
   int Nblah;
   double ptmin, ptmax;
   getline(ifile,line);
   emi = popdouble(line);
   ema = popdouble(line);
   rhomin = popdouble(line);
   rhomax = popdouble(line);
   Nblah = popdouble(line);
   ptmin = popdouble(line);
   ptmax = popdouble(line);

   etamin = emi;
   etamax = emi;
   
   std::vector<double> params;
   for(int i = 0; i < Nparam; i++)
     params.push_back(popdouble(line));
   
   FormulaBinsBins* Bins = new FormulaBinsBins(formula, Nvar, Nparam);
   m_JERFactors[year-2016][year] = Bins;

   const FormulaBins* Bin = &Bins->AddBin(etamin, etamax);
   Bin->AddBin(rhomin, rhomax, params);

   while(getline(ifile,line)){
     emi = popdouble(line);
     ema = popdouble(line);
     
     // new eta bin
     if(fabs(emi-etamin) + fabs(ema-etamax) > 1e-8){
       etamin = emi;
       etamax = ema;
       Bin = &Bins->AddBin(etamin, etamax);
     }
     
     rhomin = popdouble(line);
     rhomax = popdouble(line);
     
     Nblah = popdouble(line);
     ptmin = popdouble(line);     
     ptmax = popdouble(line);     

     params.clear();
     for(int i = 0; i < Nparam; i++)
       params.push_back(popdouble(line));
     
     Bin->AddBin(rhomin, rhomax, params);
   }

   ifile.close();
}

void JMETool::ParseJESUncertainty(const string& input, int year, const string& prefix){
  if(year < 2016 || year > 2018)
    return;
  
  std::ifstream ifile(input.c_str());
   if(!ifile.is_open()){
    std::cout << "can't open JES file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   string label = "JESUncer_Total";
   if(prefix != "")
     label = prefix+"_"+label;
   bool start = false;
   FormulaBinsBins* Bins = nullptr;

   double ptmin, ptmax, factor;
   while(getline(ifile,line)){
     // discard comment line
     if(line.find("#") != string::npos)
       continue;

     // new name, new label
     if(line.find("[") != string::npos){
       found = line.find("[");
       line.erase(0,found+1);
       found = line.find("]");
       label = "JESUncer_"+line.substr(0,found);
       if(prefix != "")
	 label = prefix+"_"+label;
       
       continue;
     }
       
     // start new bin map
     if(line.find("{") != string::npos){
       Bins = new FormulaBinsBins();
       m_Factors[year-2016][label] = Bins;
       continue;
     }

     // new eta bin for existing bin map
     double etamin = popdouble(line);
     double etamax = popdouble(line);
     
     const FormulaBins& Bin = Bins->AddBin(etamin, etamax);

     int N = popdouble(line);

     ptmin  = popdouble(line);
     while(line.length() > 0){
       factor = popdouble(line);
       popdouble(line); // factor repeated for some reason
       if(line.length() > 0)
	 ptmax = popdouble(line);
       else
	 ptmax = ptmin + 1000000;

       Bin.AddBin(ptmin, ptmax, factor);
       
       ptmin = ptmax;
     }
   }

   ifile.close();
}

void JMETool::ParseJEC(const string& input, int year, const string& prefix){
  if(year < 2016 || year > 2018)
    return;
  
  std::ifstream ifile(input.c_str());
   if(!ifile.is_open()){
    std::cout << "can't open JME file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   // first line with variable and equation info
   getline(ifile,line);
   
   // remove bracket
   found = line.find("{");
   line.erase(0,found+1);

   // get how many bins specified per line
   int Nbinvar = popdouble(line);

   for(int i = 0; i < Nbinvar; i++)
     popstring(line);

   // get the number of variables in the formula
   int Nvar = popdouble(line);

   std::vector<std::string> svar;
   for(int i = 0; i < Nvar; i++)
     svar.push_back(popstring(line));
   
   string formula = popstring(line);
  
   if(Nvar == 3){
     std::vector<std::string> var;
     std::vector<std::string> sym;
     std::vector<std::string> rep;
     var.push_back("JetPt");
     var.push_back("JetA");
     var.push_back("Rho");
     sym.push_back("x");
     sym.push_back("y");
     sym.push_back("z");
     rep.push_back("q");
     rep.push_back("v");
     rep.push_back("f");

     if(formula.find("max") != string::npos)
       formula.replace(formula.find("max"), 3, "%%%");
     
     for(int i = 0; i < 3; i++){
       while(formula.find(sym[i]) != string::npos){
	 found = formula.find(sym[i]);
	 formula.replace(found, 1, rep[i]);
       }
     }

     for(int v = 0; v < 3; v++){
       for(int s = 0; s < 3; s++){
	 if(svar[s].find(var[v]) != string::npos){
	   while(formula.find(rep[s]) != string::npos){
	     found = formula.find(rep[s]);
	     formula.replace(found, 1, sym[v]);
	   }
	 }
       }
     }
   }

   if(formula.find("%%%") != string::npos)
       formula.replace(formula.find("%%%"), 3, "max");
   
   popstring(line);

   string label = popstring(line);
   found = label.find("}");
   label.erase(found,found+1);
   if(prefix != "")
     label = prefix+"_"+label;

   // get the first line to get parameter numbers
   double etamin, etamax;
   double emi, ema;
   double ptmin, ptmax;
   getline(ifile,line);
   etamin = popdouble(line);
   etamax = popdouble(line);
   if(Nbinvar == 2){
     ptmin = popdouble(line);
     ptmax = popdouble(line);
   } else {
     ptmin = 0.;
     ptmax = 1e10;
   }

   int Nparam = popdouble(line);
   
   for(int i = 0; i < Nvar*2; i++)
     popdouble(line);
   Nparam -= Nvar*2;

   std::vector<double> params;
   for(int i = 0; i < Nparam; i++)
     params.push_back(popdouble(line));
   
   FormulaBinsBins* Bins = new FormulaBinsBins(formula, Nvar, Nparam);
   m_Factors[year-2016][label] = Bins;

   const FormulaBins* Bin = &Bins->AddBin(etamin, etamax);
   Bin->AddBin(ptmin, ptmax, params);

   while(getline(ifile,line)){
     emi = popdouble(line);
     ema = popdouble(line);
     
     // new eta bin
     if(fabs(emi-etamin) + fabs(ema-etamax) > 1e-8){
       etamin = emi;
       etamax = ema;
       Bin = &Bins->AddBin(etamin, etamax);
     }
     
     if(Nbinvar == 2){
       ptmin = popdouble(line);
       ptmax = popdouble(line);
     } else {
       ptmin = 0.;
       ptmax = 1e10;
     }
     
     int N = popdouble(line);
     for(int i = 0; i < Nvar*2; i++)
       popdouble(line);
     N -= Nvar*2;

     params.clear();
     for(int i = 0; i < N; i++)
       params.push_back(popdouble(line));
     
     Bin->AddBin(ptmin, ptmax, params);
   }

   ifile.close();
}

double JMETool::popdouble(std::string& line){
  // remove leading whitespace
  while(line[0] == string(" ")[0])
    line.erase(0,1);
  
  string num;
  if(line.find(" ") != string::npos){
    size_t p = line.find(" ");
    num = line.substr(0,p);
    line.erase(0,p+1);
    while(line[0] == string(" ")[0])
      line.erase(0,1);
  } else
    num = line;

  return stod(num);
}

std::string JMETool::popstring(std::string& line){
  // remove leading whitespace
  while(line[0] == string(" ")[0])
    line.erase(0,1);
  
  string ret;
  if(line.find(" ") != string::npos){
    size_t p = line.find(" ");
    ret = line.substr(0,p);
    line.erase(0,p+1);
    while(line[0] == string(" ")[0])
      line.erase(0,1);
  } else
    ret = line;

  return ret;
}

