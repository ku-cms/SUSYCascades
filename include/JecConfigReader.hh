#pragma once

#include <nlohmann/json.hpp>
#include <correction.h>
#include <TRandom3.h>

#include <fstream>
#include <string>
#include <optional>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include <map>
#include <memory>
#include <utility>

namespace JecConfigReader {

// ===== One-time config path registry + cached JSONs =====
struct ConfigPaths {
    std::string ak4 = "data/JME/JecConfigAK4.json";
    std::string ak8 = "data/JME/JecConfigAK8.json";
};

// Internal shared state
inline ConfigPaths& _paths_mut() {
    static ConfigPaths p;           // single shared instance
    return p;
}
inline bool& _paths_locked() {
    static bool locked = false;     // one-time lock
    return locked;
}

// Set the paths once (optional; defaults above are used otherwise)
inline void setConfigPaths(const ConfigPaths& p) {
    if (_paths_locked()) return;    // ignore subsequent changes
    _paths_mut() = p;
    _paths_locked() = true;
}

// Read-only accessor used everywhere else
inline const ConfigPaths& _paths() {
    return _paths_mut();
}


/** Read a JSON configuration file from disk. */
inline nlohmann::json loadJsonConfig(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open JSON config: " + filename);
    }
    nlohmann::json j;
    f >> j;
    return j;
}

// Internal loader with one-time caching per file path
inline const nlohmann::json& _loadOnce(const std::string& path) {
    static std::unordered_map<std::string, nlohmann::json> cache;
    auto it = cache.find(path);
    if (it != cache.end()) return it->second;
    cache.emplace(path, loadJsonConfig(path));
    return cache.at(path);
}

// Public: get cached configs
inline const nlohmann::json& getCfgAK4() { return _loadOnce(_paths().ak4); }
inline const nlohmann::json& getCfgAK8() { return _loadOnce(_paths().ak8); }

/** Get a string field from a JSON object or throw. */
inline std::string getTagName(const nlohmann::json& obj, const std::string& key) {
    if (!obj.contains(key)) {
        throw std::runtime_error("Missing required key in JSON: " + key);
    }
    return obj.at(key).get<std::string>();
}

struct Tags {
    std::string jercJsonPath;
    std::string tagNameL1FastJet;
    std::string tagNameL2Relative;
    std::string tagNameL3Absolute;
    std::string tagNameL2Residual;
    std::string tagNameL2L3Residual;
    std::string tagNamePtResolution;
    std::string tagNameJerScaleFactor;
};

/** Gather tag names for a given (year, data?, era). */
inline Tags getTagNames(const nlohmann::json& baseJson,
                        const std::string& year,
                        bool isData,
                        const std::optional<std::string>& era)
{
    if (!baseJson.contains(year)) {
        throw std::runtime_error("Year key not found in JSON: " + year);
    }
    const auto& yearObj = baseJson.at(year);
    Tags tags;

    tags.jercJsonPath          = getTagName(yearObj, "jercJsonPath");
    tags.tagNameL1FastJet      = getTagName(yearObj, "tagNameL1FastJet");
    tags.tagNameL2Relative     = getTagName(yearObj, "tagNameL2Relative");
    tags.tagNameL3Absolute     = getTagName(yearObj, "tagNameL3Absolute");
    tags.tagNameL2Residual     = getTagName(yearObj, "tagNameL2Residual");
    tags.tagNameL2L3Residual   = getTagName(yearObj, "tagNameL2L3Residual");
    tags.tagNamePtResolution   = getTagName(yearObj, "tagNamePtResolution");
    tags.tagNameJerScaleFactor = getTagName(yearObj, "tagNameJerScaleFactor");

    if (isData) {
        if (!yearObj.contains("data")) {
            throw std::runtime_error("Requested data but no 'data' section for year: " + year);
        }
        if (!era.has_value()) {
            throw std::runtime_error("Data requested but no era provided for year: " + year);
        }
        const std::string& eraKey = era.value();
        const auto& dataObj = yearObj.at("data");
        if (!dataObj.contains(eraKey)) {
            throw std::runtime_error("Era key not found under data for year " + year + ": " + eraKey);
        }
        const auto& eraObj = dataObj.at(eraKey);
        tags.tagNameL1FastJet    = getTagName(eraObj, "tagNameL1FastJet");
        tags.tagNameL2Relative   = getTagName(eraObj, "tagNameL2Relative");
        tags.tagNameL3Absolute   = getTagName(eraObj, "tagNameL3Absolute");
        tags.tagNameL2Residual   = getTagName(eraObj, "tagNameL2Residual");
        tags.tagNameL2L3Residual = getTagName(eraObj, "tagNameL2L3Residual");
    }

    return tags;
}

/** Retrieve a correction from a CorrectionSet with nice error. */
inline correction::Correction::Ref safeGet(
    const std::shared_ptr<correction::CorrectionSet>& cs,
    const std::string& name)
{
    try {
        return cs->at(name);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to retrieve correction \"" + name + "\": " + e.what());
    }
}

// ===== Helpers using cached JSONs =====
inline Tags getTagsAK4(const std::string& year, bool isData,
                       const std::optional<std::string>& era) {
    return getTagNames(getCfgAK4(), year, isData, era);
}

inline Tags getTagsAK8OrFallbackAK4(const std::string& year, bool isData,
                                    const std::optional<std::string>& era) {
    const auto& cfgAK8 = getCfgAK8();
    if (cfgAK8.contains(year)) return getTagNames(cfgAK8, year, isData, era);
    return getTagsAK4(year, isData, era);
}

struct CorrectionRefs {
    std::shared_ptr<correction::CorrectionSet> cs;
    correction::Correction::Ref corrRefJesL1FastJet;
    correction::Correction::Ref corrRefJesL2Relative;
    correction::Correction::Ref corrRefJesL2ResL3Res;
    correction::Correction::Ref corrRefJerReso;
    correction::Correction::Ref corrRefJerSf;
    TRandom3 randomGen;

    CorrectionRefs() = default;
    CorrectionRefs(const Tags& tags)
        : randomGen(0)
    {
        static std::map<std::string, std::shared_ptr<correction::CorrectionSet>> cache;

        if (cache.count(tags.jercJsonPath)) {
            cs = cache.at(tags.jercJsonPath);
        } else {
            cs = correction::CorrectionSet::from_file(tags.jercJsonPath);
            cache[tags.jercJsonPath] = cs;
        }

        corrRefJesL1FastJet   = safeGet(cs, tags.tagNameL1FastJet);
        corrRefJesL2Relative  = safeGet(cs, tags.tagNameL2Relative);
        corrRefJesL2ResL3Res  = safeGet(cs, tags.tagNameL2L3Residual);
        corrRefJerReso        = safeGet(cs, tags.tagNamePtResolution);
        corrRefJerSf          = safeGet(cs, tags.tagNameJerScaleFactor);
    }
};

struct JerBin {
    std::string label;
    double etaMin{}, etaMax{}, ptMin{}, ptMax{};
};

using JerSetMap = std::map<std::string, std::vector<JerBin>>;

inline JerSetMap getJerUncertaintySets(const nlohmann::json& baseJson, const std::string& year) {
    JerSetMap out;

    auto itYear = baseJson.find(year);
    if (itYear == baseJson.end()) return out;
    const auto& y = *itYear;

    auto itUncertainty = y.find("ForUncertaintyJER");
    if (itUncertainty == y.end() || !itUncertainty->is_object()) return out;
    const auto& j = *itUncertainty;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string setName = it.key();
        const auto& obj = it.value();
        if (!obj.is_object()) continue;

        std::vector<JerBin> bins;
        bins.reserve(obj.size());
        for (auto it2 = obj.begin(); it2 != obj.end(); ++it2) {
            const std::string label = it2.key();
            const auto& arr = it2.value(); // [etaMin, etaMax, ptMin, ptMax]
            if (!arr.is_array() || arr.size() != 4) {
                throw std::runtime_error(
                    "ForUncertaintyJER bin \"" + label + "\" must be an array [etaMin, etaMax, ptMin, ptMax]."
                );
            }
            JerBin b;
            b.label = label;
            b.etaMin = arr.at(0).get<double>();
            b.etaMax = arr.at(1).get<double>();
            b.ptMin  = arr.at(2).get<double>();
            b.ptMax  = arr.at(3).get<double>();
            bins.push_back(std::move(b));
        }
        out.emplace(setName, std::move(bins));
    }
    return out;
}

// Each entry: { fullTag, base/custom name }
using SystPairJES   = std::pair<std::string, std::string>;
using SystSetMapJES = std::map<std::string, std::vector<SystPairJES>>;

inline SystSetMapJES getSystTagNames(const nlohmann::json& baseJson, const std::string& year) {
    SystSetMapJES out;
    if (!baseJson.contains(year)) return out;

    const auto& yearObj = baseJson.at(year);
    if (!yearObj.contains("ForUncertaintyJES")) return out;

    const auto& systs = yearObj.at("ForUncertaintyJES");
    for (auto it = systs.begin(); it != systs.end(); ++it) {
        const std::string setName = it.key();
        const auto& arr = it.value();
        if (!arr.is_array()) continue;

        std::vector<SystPairJES> pairs;
        pairs.reserve(arr.size());

        for (const auto& item : arr) {
            if (item.is_array() && item.size() >= 2 && item.at(0).is_string() && item.at(1).is_string()) {
                pairs.emplace_back(item.at(0).get<std::string>(), item.at(1).get<std::string>());
            }
        }
        out.emplace(setName, std::move(pairs));
    }
    return out;
}

enum class SystKind { Nominal, JES, JER };

struct SystTagDetail {
    std::string setName;
    std::string var;   // "Up"/"Down" (empty for nominal)
    SystKind kind{SystKind::Nominal};
    bool isNominal() const { return kind == SystKind::Nominal; }
    std::string systSetName() const { return isNominal() ? "Nominal" : setName; }
    virtual std::string systName() const { return isNominal() ? "Nominal" : setName + "_" + var; }
    virtual ~SystTagDetail() = default;
};

struct SystTagDetailJES : public SystTagDetail {
    std::string tagAK4;
    std::string tagAK8;
    std::string baseTag;
    std::string systName() const override { return isNominal() ? "Nominal" : baseTag + "_" + var; }
};

struct SystTagDetailJER : public SystTagDetail {
    std::string baseTag;
    JerBin jerRegion;
    std::string systName() const override { return isNominal() ? "Nominal" : baseTag + "_" + var; }
};

inline std::vector<SystTagDetailJER> buildJerTagDetails(const JerSetMap& jerSets) {
    std::vector<SystTagDetailJER> out;
    for (const auto& [setName, bins] : jerSets) {
        for (const auto& b : bins) {
            for (const char* var : {"Up","Down"}) {
                SystTagDetailJER d;
                d.setName   = setName;
                d.var       = var;
                d.kind      = SystKind::JER;
                d.baseTag   = b.label;
                d.jerRegion = b;
                out.push_back(d);
            }
        }
    }
    return out;
}

inline std::vector<SystTagDetailJES> buildSystTagDetailJES(const SystSetMapJES& sAK4,
                                                            const SystSetMapJES& sAK8)
{
    std::vector<SystTagDetailJES> systTagDetails;

    for (const auto& [set, pairsAK4] : sAK4) {
        auto itAK8 = sAK8.find(set);
        if (itAK8 == sAK8.end()) continue;
        const auto& pairsAK8 = itAK8->second;

        std::unordered_map<std::string, std::string> mapAK4, mapAK8;
        mapAK4.reserve(pairsAK4.size());
        mapAK8.reserve(pairsAK8.size());

        for (const auto& p : pairsAK4) mapAK4.emplace(p.second, p.first); // key = base/custom
        for (const auto& p : pairsAK8) mapAK8.emplace(p.second, p.first);

        for (const auto& [base, fullAK4] : mapAK4) {
            auto itFullAK8 = mapAK8.find(base);
            if (itFullAK8 == mapAK8.end()) continue;

            for (const char* var : {"Up","Down"}) {
                SystTagDetailJES d;
                d.setName = set;
                d.var = var;
                d.kind = SystKind::JES;
                d.tagAK4 = fullAK4;
                d.tagAK8 = itFullAK8->second;
                d.baseTag = base;
                systTagDetails.push_back(d);
            }
        }
    }
    return systTagDetails;
}

} // namespace JecConfigReader

