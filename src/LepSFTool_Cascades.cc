#include "../include/LepSFTool_Cascades.hh"
#include <iostream>
#include <cmath>

using std::string;

LepSFTool::LepSFTool(const std::string& filename) {
    if(filename != "")
        cset_ = correction::CorrectionSet::from_file(filename);
}

LepSFTool::~LepSFTool() {}

void LepSFTool::BuildMap(const std::string& filename) {
    try {
        cset_ = correction::CorrectionSet::from_file(filename);
    } catch(const std::exception& e) {
        std::cerr << "[LepSFTool] ERROR: failed to load SF json: " << filename
                  << "\n  " << e.what() << "\n";
        std::exit(1);
    }
}

std::string LepSFTool::GetFlav(int pdg) const {
    return (std::abs(pdg) == 11 || pdg == 0) ? "electron" : "muon";
}

double LepSFTool::Lookup(const std::string& sf_key, double pt, double eta) const {
    auto sf_map = cset_->at(sf_key);
    return sf_map->evaluate({pt, std::abs(eta), year_, era_});
}

// BLP_COL: all qualities
double LepSFTool::get_BLP_COL_SF(double pt, double eta, int pdg, LepID qual) {
    return Lookup(GetFlav(pdg) + "_BLP_over_COL", pt, eta);
}
double LepSFTool::get_BLP_COL_SF_err(double pt, double eta, int pdg, LepID qual) {
    return Lookup(GetFlav(pdg) + "_BLP_over_COL_err", pt, eta);
}

// ID_BLP: gold and silver only
double LepSFTool::get_ID_BLP_SF(double pt, double eta, int pdg, LepID qual) {
    if(qual == kBronze) return 1.0;
    return Lookup(GetFlav(pdg) + "_ID_over_BLP", pt, eta);
}
double LepSFTool::get_ID_BLP_SF_err(double pt, double eta, int pdg, LepID qual) {
    if(qual == kBronze) return 0.0;
    return Lookup(GetFlav(pdg) + "_ID_over_BLP_err", pt, eta);
}

// ISO_ID: gold and silver only
double LepSFTool::get_ISO_ID_SF(double pt, double eta, int pdg, LepID qual) {
    if(qual == kBronze) return 1.0;
    return Lookup(GetFlav(pdg) + "_ISO_over_ID", pt, eta);
}
double LepSFTool::get_ISO_ID_SF_err(double pt, double eta, int pdg, LepID qual) {
    if(qual == kBronze) return 0.0;
    return Lookup(GetFlav(pdg) + "_ISO_over_ID_err", pt, eta);
}

// Prompt_ISOID: gold only
double LepSFTool::get_Prompt_ISOID_SF(double pt, double eta, int pdg, LepID qual) {
    if(qual == kSilver || qual == kBronze) return 1.0;
    return Lookup(GetFlav(pdg) + "_PROMPT_over_ISOID", pt, eta);
}
double LepSFTool::get_Prompt_ISOID_SF_err(double pt, double eta, int pdg, LepID qual) {
    if(qual == kSilver || qual == kBronze) return 0.0;
    return Lookup(GetFlav(pdg) + "_PROMPT_over_ISOID_err", pt, eta);
}

// NOT_Prompt_ISOID: silver only
double LepSFTool::get_NOT_Prompt_ISOID_SF(double pt, double eta, int pdg, LepID qual) {
    if(qual == kGold || qual == kBronze) return 1.0;
    return Lookup(GetFlav(pdg) + "_NOT_PROMPT_over_ISOID", pt, eta);
}
double LepSFTool::get_NOT_Prompt_ISOID_SF_err(double pt, double eta, int pdg, LepID qual) {
    if(qual == kGold || qual == kBronze) return 0.0;
    return Lookup(GetFlav(pdg) + "_NOT_PROMPT_over_ISOID_err", pt, eta);
}

// NOT_ID_nor_ISO: bronze only
double LepSFTool::get_NOT_ID_nor_ISO_SF(double pt, double eta, int pdg, LepID qual) {
    if(qual == kGold || qual == kSilver) return 1.0;
    return Lookup(GetFlav(pdg) + "_NOT_ID_nor_ISO", pt, eta);
}
double LepSFTool::get_NOT_ID_nor_ISO_SF_err(double pt, double eta, int pdg, LepID qual) {
    if(qual == kGold || qual == kSilver) return 0.0;
    return Lookup(GetFlav(pdg) + "_NOT_ID_nor_ISO_err", pt, eta);
}