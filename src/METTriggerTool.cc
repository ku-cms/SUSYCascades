#include "../include/METTriggerTool.hh"

METTriggerTool::METTriggerTool(){
  m_filename = "Parameters.csv";
}
    
METTriggerTool::~METTriggerTool(){

}

void METTriggerTool::BuildMap(const string& csv){
  SetCSV(csv);
  ParseCSV();
}

double METTriggerTool::Get_EFF(double MET, double HT, int year,
			       bool el, bool mu,
			       bool data, int updown){
  string name = Get_Name(HT, year, el, mu, data);
  return Get_EFF(name, MET, updown);
}

double METTriggerTool::Get_SF(double MET, double HT, int year,
			       bool el, bool mu,
			       bool data, int updown){
  string data_name, bkg_name;

  double SF;
  
  if(!el && !mu){
    data_name = Get_Name(HT, year, true, mu, true); 
    bkg_name = Get_Name(HT, year, true, mu, false);
    SF = Get_EFF(data_name, MET)/Get_EFF(bkg_name, MET);
    data_name = Get_Name(HT, year, el, mu, true);
    bkg_name = Get_Name(HT, year, el, mu, false);
  } else {
    data_name = Get_Name(HT, year, el, mu, true); 
    bkg_name = Get_Name(HT, year, el, mu, false);
    SF = Get_EFF(data_name, MET)/Get_EFF(bkg_name, MET);
  }
  
  if(updown == 0)
    return SF; 
  else if(updown > 0)
    return SF*(m_map[bkg_name][4] + (MET < m_map[bkg_name][2] ? m_map[bkg_name][0]*(MET-m_map[bkg_name][2])*(MET-m_map[bkg_name][2]) : 0.)); 
  else
    return SF*(m_map[bkg_name][5] + (MET < m_map[bkg_name][3] ? m_map[bkg_name][1]*(MET-m_map[bkg_name][3])*(MET-m_map[bkg_name][3]) : 0.));

  return 1.;  
}

void METTriggerTool::SetCSV(const string& csv){
  m_filename = csv;
}

void METTriggerTool::ParseCSV(){
  std::ifstream ifile(m_filename);
  if(!ifile.is_open()){
    std::cout << "can't open csv file " << std::endl;
    return;
  }
  std::string line;
  //discard first line
  getline(ifile,line);
  while(getline(ifile,line)){
    size_t comma = line.find(",");
    std::string name = "";
    std::vector<double> values;
    name = line.substr(0,comma);
    line.erase(0,(name+",").length());
    while(line.length() > 0) { 
      values.push_back(Get_Value(line)); 
    }
    m_map.insert(std::make_pair(name,values));
  }
}

string METTriggerTool::Get_Name(double HT, int year, bool el, bool mu, bool data){
  string name = "";
  
  if(HT <= 600.)
    name+="HT-Le600--"; 
  else if(HT > 600. && HT < 750.)
    name+="HT-G600--HT-L750--"; 
  else if(HT >= 750.)
    name+="HT-Ge750--";
  
  if(el)
    name+="SingleElectrontrigger-E1--Nele-E1_"; 
  else if(mu)
    name+="SingleMuontrigger-E1--Nmu-E1_"; 
  else
    name+="Nlep-E0_";
  
  if(data){
    if(el)
      name+="SingleElectron_"+std::to_string(year)+"_Electron"; 
    else if(mu)
      name+="SingleMuon_"+std::to_string(year)+"_Muon"; 
    else {
      name+="SingleElectron_"+std::to_string(year)+"_ZeroLepton";
    }
  } else {
    name+="Bkg_"+std::to_string(year)+"_";
    if(el)
      name+="Electron"; 
    else if(mu)
      name+="Muon"; 
    else
      name+="ZeroLepton"; 
  }

  return name;
}

int METTriggerTool::Get_Curve_Index(double HT, int year, bool el, bool mu, bool data){
  int index = 0;

  if(HT <= 600.)
    index += 1;
  else if(HT > 600. && HT < 750.)
    index += 2;
  else if(HT >= 750.)
    index += 3;
  
  if(el)
    index += 10;
  else if(mu)
    index += 20;
  else
    index += 30;
  
  if(data){
    if(el)
      index += 100;
    else if(mu)
      index += 200;
    else
      index += 300;
  } else {
    if(el)
      index += 400;
    else if(mu)
      index += 500;
    else
      index += 600;
  }

  return index;
}

double METTriggerTool::Get_EFF(string name, double MET, int updown){

  double EFF = 0.;
  if(m_map[name][9] == 0 && m_map[name][10] == 0){ 
    EFF = m_map[name][6]*ROOT::Math::normal_cdf(MET, m_map[name][8], m_map[name][7]); 
  } else {
    EFF = m_map[name][6]*(TMath::Cos(m_map[name][10])*TMath::Cos(m_map[name][10])*ROOT::Math::normal_cdf(MET, m_map[name][8], m_map[name][7])+TMath::Sin(m_map[name][10])*TMath::Sin(m_map[name][10])*ROOT::Math::normal_cdf(MET, m_map[name][8]*m_map[name][9], m_map[name][7]));
  }
  
  if(updown == 0)
    return EFF;
  else if(updown > 0)
    return std::min(1., EFF*(m_map[name][4] + (MET < m_map[name][2] ? m_map[name][0]*(MET-m_map[name][2])*(MET-m_map[name][2]) : 0.))); 
  else
    return std::min(1., EFF*(m_map[name][5] + (MET < m_map[name][3] ? m_map[name][1]*(MET-m_map[name][3])*(MET-m_map[name][3]) : 0.))); 

  return 1.;
}

double METTriggerTool::Get_Value(std::string& line){
  size_t comma_pos = line.find(",");
  std::string value = line.substr(0,comma_pos);
  line.erase(0,comma_pos+1);
  
  return std::stod(value);
}
static void JSON_SkipWS(const std::string& s, size_t& p) {
    while (p < s.size() && std::isspace((unsigned char)s[p])) ++p;
}
 
static bool JSON_ReadString(const std::string& s, size_t& p, std::string& out) {
    JSON_SkipWS(s, p);
    if (p >= s.size() || s[p] != '"') return false;
    ++p;
    out.clear();
    while (p < s.size() && s[p] != '"') {
        if (s[p] == '\\') ++p;
        out += s[p++];
    }
    if (p < s.size()) ++p;
    return true;
}
 
static std::string JSON_CaptureValue(const std::string& s, size_t& p) {
    JSON_SkipWS(s, p);
    if (p >= s.size()) return "";
    if (s[p] == '{') {
        size_t start = p; int depth = 0;
        while (p < s.size()) {
            if (s[p] == '{') ++depth;
            else if (s[p] == '}') { if (--depth == 0) { ++p; break; } }
            ++p;
        }
        return s.substr(start, p - start);
    }
    if (s[p] == '"') {
        std::string tmp; JSON_ReadString(s, p, tmp); return "\"" + tmp + "\"";
    }
    size_t start = p;
    while (p < s.size() && s[p] != ',' && s[p] != '}') ++p;
    return s.substr(start, p - start);
}
 
// Parses { "k":v, ... } into a string->string map.
static void JSON_ParseObject(const std::string& s, size_t& p,
                              std::map<std::string,std::string>& out) {
    JSON_SkipWS(s, p);
    if (p >= s.size() || s[p] != '{') return;
    ++p;
    while (true) {
        JSON_SkipWS(s, p);
        if (p >= s.size() || s[p] == '}') { ++p; break; }
        if (s[p] == ',') { ++p; continue; }
        std::string key;
        if (!JSON_ReadString(s, p, key)) break;
        JSON_SkipWS(s, p);
        if (p < s.size() && s[p] == ':') ++p;
        out[key] = JSON_CaptureValue(s, p);
    }
}
 
// Parses a flat { "k": number } object into a double map.
static void JSON_ParseDoubles(const std::string& blob,
                               std::map<std::string,double>& out) {
    size_t p = 0;
    std::map<std::string,std::string> raw;
    JSON_ParseObject(blob, p, raw);
    for (auto& kv : raw) {
        try { out[kv.first] = std::stod(kv.second); } catch (...) {}
    }
}
 
static double JSON_Get(const std::map<std::string,double>& m,
                        const std::string& k, double def) {
    auto it = m.find(k); return it != m.end() ? it->second : def;
}
 
// storage key used by the new overloads.
static std::string JSON_StoreKey(int Nele, int year, bool isData) {
    return std::string("Electron") + std::to_string(Nele) +
           (isData ? "__data_" : "__bkg_") + std::to_string(year);
}
 
// Parse the JSON key to extract (Nele, year, isData, isComb).
// Returns false for unrecognised keys.
static bool JSON_DecodeKey(const std::string& key,
                            int& Nele, int& year, bool& isData, bool& isComb) {
    if (key.size() < 9) return false;
    if (key.substr(0, 8) != "Electron") return false;
    if (!std::isdigit((unsigned char)key[8])) return false;
    Nele = key[8] - '0';
 
    size_t dd = key.find("__");
    if (dd == std::string::npos) return false;
    std::string rest = key.substr(dd + 2);
 
    if (rest.substr(0, 4) == "bkg_") {
        isData = false;
        std::string yearTok = rest.substr(4); // e.g. "2017" or "2016Comb"
        isComb = (yearTok.find("Comb") != std::string::npos);
        size_t i = 0;
        while (i < yearTok.size() && std::isdigit((unsigned char)yearTok[i])) ++i;
        try { year = std::stoi(yearTok.substr(0, i)); } catch (...) { return false; }
        return true;
    }
    if (rest.substr(0, 5) == "Data_") {
        isData = true;
        // rest = "Data_<Stream>_<YearTok>"
        size_t last = rest.rfind('_');
        if (last == std::string::npos) return false;
        std::string yearTok = rest.substr(last + 1);
        isComb = (yearTok.find("Comb") != std::string::npos);
        size_t i = 0;
        while (i < yearTok.size() && std::isdigit((unsigned char)yearTok[i])) ++i;
        try { year = std::stoi(yearTok.substr(0, i)); } catch (...) { return false; }
        return true;
    }
    return false;
}
 
void METTriggerTool::BuildMapJSON(const std::string& jsonFile) {
    std::ifstream ifs(jsonFile);
    if (!ifs.is_open()) {
        std::cout << "[METTriggerTool] Cannot open JSON: " << jsonFile << std::endl;
        return;
    }
    std::ostringstream ss; ss << ifs.rdbuf();
    std::string src = ss.str();
 
    // Collect all top-level entries.
    struct RawEntry {
        int  Nele;
        int  year;
        bool isData;
        bool isComb;
        std::map<std::string,std::string> fields; // "fit", "bands", fGauss/fDGauss keys
    };
    std::vector<RawEntry> entries;
 
    size_t pos = 0;
    JSON_SkipWS(src, pos);
    if (pos >= src.size() || src[pos] != '{') return;
    ++pos;
 
    while (true) {
        JSON_SkipWS(src, pos);
        if (pos >= src.size() || src[pos] == '}') break;
        if (src[pos] == ',') { ++pos; continue; }
 
        std::string entryKey;
        if (!JSON_ReadString(src, pos, entryKey)) break;
        JSON_SkipWS(src, pos);
        if (pos < src.size() && src[pos] == ':') ++pos;
 
        // Capture entry value (a nested object)
        std::string blob = JSON_CaptureValue(src, pos);
 
        RawEntry e;
        if (!JSON_DecodeKey(entryKey, e.Nele, e.year, e.isData, e.isComb))
            continue;
 
        size_t bp = 0;
        JSON_ParseObject(blob, bp, e.fields);
        entries.push_back(std::move(e));
    }
 
    // Two-pass insertion: non-Comb first, then Comb (overwrites) -> Comb wins
    auto insert = [&](bool combPass) {
        for (auto& e : entries) {
            if (e.isComb != combPass) continue;
            if (e.fields.find("fit") == e.fields.end()) continue;
 
            std::string storeKey = JSON_StoreKey(e.Nele, e.year, e.isData);
 
            std::map<std::string,double> fitMap;
            JSON_ParseDoubles(e.fields.at("fit"), fitMap);
 
            // Build the vector<double> in the same layout as ParseCSV:
            // [0]=a1 [1]=a2 [2]=b1 [3]=b2 [4]=c1 [5]=c2
            // [6]=norm [7]=sigma [8]=mean [9]=scale [10]=weight
            std::vector<double> v(11, 0.);
            v[4] = 1.; v[5] = 1.; // c1, c2 default to 1
 
            // bands (bkg only)
            auto bandsIt = e.fields.find("bands");
            if (bandsIt != e.fields.end()) {
                std::map<std::string,double> bm;
                JSON_ParseDoubles(bandsIt->second, bm);
                v[0] = JSON_Get(bm, "a1", 0.);
                v[1] = JSON_Get(bm, "a2", 0.);
                v[2] = JSON_Get(bm, "b1", 0.);
                v[3] = JSON_Get(bm, "b2", 0.);
                v[4] = JSON_Get(bm, "c1", 1.);
                v[5] = JSON_Get(bm, "c2", 1.);
            }
 
            v[6]  = JSON_Get(fitMap, "norm",   1.);
            v[7]  = JSON_Get(fitMap, "sigma",  1.);
            v[8]  = JSON_Get(fitMap, "mean",   0.);
            v[9]  = JSON_Get(fitMap, "scale",  0.);
            v[10] = JSON_Get(fitMap, "weight", 0.);
 
            m_map[storeKey] = v;
        }
    };

    insert(false); // non-Comb first
    insert(true);  // Comb second -> overwrites non-Comb for the same (Nele,year)
}
 
// -----------------------------------------------------------------------------
// JSON_EvalEFF: evaluates the Add-formula double-Gaussian CDF used by the
// fitting code (Double_Gaussian_CDF_Func_Add):
//   norm * ( cos²(w)*CDF(x, sigma, mean) + sin²(w)*CDF(x, sigma, mean+scale) )
// Single Gaussian when scale==0 && weight==0.
// v layout: [0]=a1 [1]=a2 [2]=b1 [3]=b2 [4]=c1 [5]=c2
//           [6]=norm [7]=sigma [8]=mean [9]=scale [10]=weight
// -----------------------------------------------------------------------------
static double JSON_EvalEFF(const std::vector<double>& v, double MET, int updown) {
    double norm   = v[6], sigma = v[7], mean  = v[8];
    double scale  = v[9], w     = v[10];
    double eff;
    if (scale == 0. && w == 0.) {
        eff = norm * ROOT::Math::normal_cdf(MET, sigma, mean);
    } else {
        double cos2 = TMath::Cos(w)*TMath::Cos(w);
        double sin2 = TMath::Sin(w)*TMath::Sin(w);
        eff = norm * (
              cos2 * ROOT::Math::normal_cdf(MET, sigma, mean)
            + sin2 * ROOT::Math::normal_cdf(MET, sigma + scale, mean)
        );
    }
    if (updown == 0) return eff;
    if (updown > 0)
        return std::min(1., eff * (v[4] + (MET < v[2] ? v[0]*(MET-v[2])*(MET-v[2]) : 0.)));
    else
        return std::min(1., eff * (v[5] + (MET < v[3] ? v[1]*(MET-v[3])*(MET-v[3]) : 0.)));
    return 1.;
}
 
double METTriggerTool::Get_EFF_JSON(double MET, int year, int Nele, int updown) {
    if(MET < 150.) return 1.;
    std::string name = JSON_StoreKey(Nele, year, /*isData=*/true);
    auto it = m_map.find(name);
    if (it == m_map.end()) {
        std::cout << "[METTriggerTool] Get_EFF: key not found: " << name << std::endl;
        return 1.;
    }
    return JSON_EvalEFF(it->second, MET, updown);
}
 
double METTriggerTool::Get_SF_JSON(double MET, int year, int Nele, int updown) {
    std::string data_name = JSON_StoreKey(Nele, year, /*isData=*/true);
    std::string bkg_name  = JSON_StoreKey(Nele, year, /*isData=*/false);
 
    auto data_it = m_map.find(data_name);
    auto bkg_it  = m_map.find(bkg_name);
 
    if (data_it == m_map.end()) {
        std::cout << "[METTriggerTool] Get_SF: data key not found: " << data_name << std::endl;
        return 1.;
    }
    if (bkg_it == m_map.end()) {
        std::cout << "[METTriggerTool] Get_SF: bkg key not found: " << bkg_name << std::endl;
        return 1.;
    }
    const std::vector<double>& bv = bkg_it->second;
    if (bv.size() < 11) {
        std::cout << "[METTriggerTool] Get_SF: bkg vector too short: " << bkg_name << std::endl;
        return 1.;
    }
 
    double eff_bkg = JSON_EvalEFF(bv, MET, 0);
    if (eff_bkg == 0.) return 1.;
    double eff_data = JSON_EvalEFF(data_it->second, MET, 0);
    double SF = eff_data / eff_bkg;
 
    if (updown == 0)
        return SF;
    else if (updown > 0)
        return SF * (bv[4] + (MET < bv[2] ? bv[0]*(MET-bv[2])*(MET-bv[2]) : 0.));
    else
        return SF * (bv[5] + (MET < bv[3] ? bv[1]*(MET-bv[3])*(MET-bv[3]) : 0.));
    return 1.;
}
