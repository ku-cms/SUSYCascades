#pragma once
#include <correction.h>
#include <cmath>
#include <string>

namespace JvmApplication {

struct JetInputs {
    double eta{0.0};
    double phi{0.0};
    double pt{0.0};
    int    jetId{0};
    double chEmEF{0.0};
    double neEmEF{0.0};
};

class VetoChecker {
public:
    // Selection / map limits (adjust centrally here if guidelines change)
    static constexpr double kMaxEtaInMap = 5.191;
    static constexpr double kMaxPhiInMap = 3.1415926;
    static constexpr double kMinPt       = 15.0;
    static constexpr int    kMinJetId    = 6;      // TightLepVeto ID
    static constexpr double kMaxEmFrac   = 0.90;   // (charged + neutral) EM fraction

    VetoChecker() = default;
    VetoChecker(const correction::Correction::Ref& jvmRef, std::string jvmKeyName) {
        setJvm(jvmRef, std::move(jvmKeyName));
    }

    // Configure once per event/file
    void setJvm(const correction::Correction::Ref& jvmRef, std::string jvmKeyName) {
        jvmRef_     = jvmRef;
        jvmKeyName_ = std::move(jvmKeyName);
        jvmValid_   = static_cast<bool>(jvmRef_);
    }

    bool isConfigured() const { return jvmValid_; }
    const std::string& jvmKeyName() const { return jvmKeyName_; }

    // Stateless check (pass all jet inputs)
    bool checkJetInVetoRegion(double eta, double phi, double pt, int jetId,
                              double chEmEF, double neEmEF) const {
        if (!jvmValid_) return false;
        JetInputs j{eta, phi, pt, jetId, chEmEF, neEmEF};
        if (!passMinimalSelection(j)) return false;
        return insideVeto(j);
    }

private:
    bool passMinimalSelection(const JetInputs& j) const {
        if (std::abs(j.eta) > kMaxEtaInMap) return false;
        if (std::abs(j.phi) > kMaxPhiInMap) return false;
        if (j.jetId < kMinJetId)            return false;
        if (j.pt < kMinPt)                  return false;
        if ((j.chEmEF + j.neEmEF) > kMaxEmFrac) return false;
        return true;
    }

    bool insideVeto(const JetInputs& j) const {
        // By convention: 0.0 outside veto, 100.0 (or >0) inside veto
        const double val = jvmRef_->evaluate({ jvmKeyName_, j.eta, j.phi });
        return (val > 0.0);
    }

    correction::Correction::Ref jvmRef_{};
    std::string jvmKeyName_{};
    bool jvmValid_{false};

    JetInputs jet_{};
};

} // namespace JvmApplication

