#pragma once
#include <nlohmann/json.hpp>
#include <correction.h>

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <utility>

namespace JvmConfigReader {

// ---------- internals (no mutex) ----------
inline std::string& cfgPathRef() {
    static std::string path = "data/JME/JvmConfig.json"; // default for backward compat
    return path;
}
inline nlohmann::json& cfgRef() {
    static nlohmann::json cfg;
    return cfg;
}
inline bool& loadedRef() {
    static bool loaded = false;
    return loaded;
}

// ---------- public path helpers ----------
inline const std::string& getConfigPath() { return cfgPathRef(); }

// Set config path. If reload=true, next config() call re-reads the JSON.
inline void setConfigPath(std::string path, bool reload = false) {
    cfgPathRef() = std::move(path);
    if (reload) {
        loadedRef() = false;
        cfgRef() = nlohmann::json(); // drop old contents
    }
}

// Explicitly invalidate cached JSON so next config() re-reads file.
inline void resetConfig() {
    loadedRef() = false;
    cfgRef() = nlohmann::json();
}

// ---------- main types ----------
struct Jvm {
    correction::Correction::Ref ref{};  // valid if use == true
    std::string key{};
    bool use{false};
};

// ---------- config loader (lazy, cached) ----------
inline const nlohmann::json& config() {
    if (!loadedRef()) {
        std::ifstream f(getConfigPath());
        if (!f)
            throw std::runtime_error(
                std::string("Cannot open JSON config: ") + getConfigPath());
        f >> cfgRef();
        loadedRef() = true;
    }
    return cfgRef();
}

// ---------- CorrectionSet cache (per file) ----------
inline std::map<std::string, std::shared_ptr<correction::CorrectionSet>>& csCacheRef() {
    static std::map<std::string, std::shared_ptr<correction::CorrectionSet>> cache;
    return cache;
}

inline std::shared_ptr<correction::CorrectionSet> loadCs(const std::string& file) {
    auto& cache = csCacheRef();
    if (auto it = cache.find(file); it != cache.end()) return it->second;

    // from_file returns unique_ptr → promote to shared_ptr
    auto up = correction::CorrectionSet::from_file(file);
    std::shared_ptr<correction::CorrectionSet> cs(std::move(up));
    cache[file] = cs;
    return cs;
}

// Optional: clear the CS cache (useful if you swap configs completely)
inline void clearCsCache() { csCacheRef().clear(); }

// ---------- main API ----------
inline Jvm getJvmForYear(const std::string& year) {
    const auto& j = config();
    if (!j.contains(year)) return {}; // use=false

    const auto& y = j.at(year);
    const std::string file = y.at("jvmFilePath").get<std::string>();
    const std::string tag  = y.at("jvmTagName").get<std::string>();
    const std::string key  = y.at("jvmKeyName").get<std::string>();

    auto cs = loadCs(file);
    return { cs->at(tag), key, true };
}

} // namespace 


