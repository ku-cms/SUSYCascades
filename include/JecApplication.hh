#pragma once
#include <correction.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include <TLorentzVector.h>
#include <TVector2.h>

#include "JecConfigReader.hh"

namespace JecApplication {

// ======================================================
// Utilities & switches
// ======================================================
inline bool hasPhiDependentL2(const std::string& year) {
    return (year == "2023Post" || year == "2024");
}

inline bool requiresRunBasedResidual(const std::string& year) {
    return (year == "2023Pre" || year == "2023Post" || year == "2024");
}

inline bool usesPuppiMet(const std::string& year) {
    static const std::unordered_set<std::string> years{
        "2022Pre", "2022Post", "2023Pre", "2023Post", "2024"
    };
    return years.count(year) != 0;
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
    const double dEta = eta1 - eta2;
    const double dPhi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return std::sqrt(dEta * dEta + dPhi * dPhi);
}

// Small helper to print when debugging.
inline void debugPrint(bool isDebug, const std::string& msg) {
    if (isDebug) std::cout << "[JERC DEBUG] " << msg << '\n';
}

// ======================================================
// Input structs (compact, intent-revealing)
// ======================================================

// Jet energy scale inputs (state at evaluation point).
struct JesInputs {
    double pt{0.0};         // If using NanoAOD Jet_pt, this is corrected pt
    double eta{0.0};
    double phi{0.0};
    double area{0.0};
    double rho{0.0};
    double rawFactor{0.0};  // NanoAOD Jet_rawFactor (0..1). If using RAW input, set 0.0
};

// JER extras: gen-match + event seeding.
struct JerInputs {
    bool   hasGen{false};
    double genPt{0.0};
    double genEta{0.0};
    double genPhi{0.0};

    // Unique event seed inputs
    unsigned long long event{0};
    unsigned int       run{0};
    unsigned int       lumi{0};

    // Matching thresholds
    double maxDr{0.2};  // ~0.2 (AK4), ~0.4–0.8 (AK8)
};

// Minimal per-jet info used during MET propagation.
struct JetForMet {
    double phi{0.0};
    double eta{0.0};
    double area{0.0};
    double rawPt{0.0};           // Raw Jet pt (or rawMinusMuon below if you prefer)
    double muonSubtrFactor{0.0}; // NanoAOD Jet_muonSubtrFactor
    double chEmEf{0.0};
    double neEmEf{0.0};
};

// Raw MET snapshot.
struct MetInputs {
    double metPt{0.0};
    double metPhi{0.0};
};

// JES / JER options that travel together.
struct SystematicOptions {
    // JER: "nom" | "up" | "down"
    std::string jerVar{"nom"};
    // If set, only jets in this bin use jerVar; others use "nom".
    std::optional<JecConfigReader::JerBin> jerRegion{std::nullopt};

    // JES (MC only): uncertainty source and direction ("Up" | "Down")
    std::string jesSystName{};  // empty => no JES syst
    std::string jesSystVar{};   // "Up" or "Down"
};

// Shared call context for repeated operations.
struct EvalContext {
    std::string year{};
    JecConfigReader::CorrectionRefs& refs;
    bool isData{false};
    std::optional<double> runNumber{std::nullopt};
    bool isDebug{false};
};

// ======================================================
// 1) Nominal JES factor (pt_corr = pt_input * S)
// ======================================================
inline double getJesCorrectionNom(const EvalContext& ctx,
                            const JesInputs& j) {
    // Start from user's input pt
    const double ptStart = j.pt;
    if (ptStart <= 0.0) return 1.0;

    // Undo NanoAOD correction if present → pt_raw_from_input
    const double rawSf = 1.0 - j.rawFactor;  // in [0,1]; 0 => already raw
    double ptRaw = ptStart * rawSf;
    debugPrint(ctx.isDebug, "JES: ptStart=" + std::to_string(ptStart) +
                            ", rawFactor=" + std::to_string(j.rawFactor) +
                            ", ptRaw=" + std::to_string(ptRaw));

    // L1FastJet
    double ptAfter = ptRaw;
    const double c1 = ctx.refs.corrRefJesL1FastJet->evaluate({ j.area, j.eta, ptAfter, j.rho });
    ptAfter *= c1;
    debugPrint(ctx.isDebug, "JES L1FastJet: c1=" + std::to_string(c1) +
                            ", ptAfter=" + std::to_string(ptAfter));

    // L2Relative (eta or eta,phi depending on era)
    double c2 = 1.0;
    if (hasPhiDependentL2(ctx.year)) {
        c2 = ctx.refs.corrRefJesL2Relative->evaluate({ j.eta, j.phi, ptAfter });
    } else {
        c2 = ctx.refs.corrRefJesL2Relative->evaluate({ j.eta, ptAfter });
    }
    ptAfter *= c2;
    debugPrint(ctx.isDebug, "JES L2Relative: c2=" + std::to_string(c2) +
                            ", ptAfter=" + std::to_string(ptAfter));

    // Residuals for data (optionally run-based)
    if (ctx.isData) {
        double cRes = 1.0;
        if (requiresRunBasedResidual(ctx.year)) {
            const double run = ctx.runNumber.value_or(0.0);
            cRes = ctx.refs.corrRefJesL2ResL3Res->evaluate({ run, j.eta, ptAfter });
        } else {
            cRes = ctx.refs.corrRefJesL2ResL3Res->evaluate({ j.eta, ptAfter });
        }
        ptAfter *= cRes;
        debugPrint(ctx.isDebug, "JES Residuals: cRes=" + std::to_string(cRes) +
                                ", ptAfter=" + std::to_string(ptAfter));
    }

    const double scale = ptAfter / ptStart;
    debugPrint(ctx.isDebug, "JES final scale=" + std::to_string(scale));
    return scale;
}

// ======================================================
// 2) JER nominal / up / down (factor only)
// Multiply *after* JES has been applied.
// ======================================================
inline double getJerCorrectionNomOrSyst(const EvalContext& ctx,
                           const JesInputs& jAfterJes,  // pt/eta/phi after JES
                           const JerInputs& jer,
                           const SystematicOptions& opts) {
    // Region gating
    const auto inRegion = [&](const JecConfigReader::JerBin& b) {
        const double aeta = std::fabs(jAfterJes.eta);
        return (aeta >= b.etaMin && aeta < b.etaMax &&
                jAfterJes.pt >= b.ptMin && jAfterJes.pt < b.ptMax);
    };
    const std::string useVar = (opts.jerRegion && !inRegion(*opts.jerRegion)) ? "nom" : opts.jerVar;

    // Resolution and scale factor
    const double reso = ctx.refs.corrRefJerReso->evaluate({ jAfterJes.eta, jAfterJes.pt, jAfterJes.rho });
    double sf = 1.0;
    if (usesPuppiMet(ctx.year)) {
        sf = ctx.refs.corrRefJerSf->evaluate({ jAfterJes.eta, jAfterJes.pt, useVar });
    } else {
        sf = ctx.refs.corrRefJerSf->evaluate({ jAfterJes.eta, useVar });
    }
    debugPrint(ctx.isDebug, "JER inputs: useVar=" + useVar +
                            ", reso=" + std::to_string(reso) +
                            ", sf=" + std::to_string(sf));

    // Seed RNG with event key
    ctx.refs.randomGen.SetSeed(jer.event + jer.run + jer.lumi);

    // Try gen-matched hybrid method
    bool matched = false;
    if (jer.hasGen) {
        const double dR = deltaR(jAfterJes.eta, jAfterJes.phi, jer.genEta, jer.genPhi);
        if (dR < jer.maxDr && std::abs(jAfterJes.pt - jer.genPt) < 3.0 * reso * jAfterJes.pt) {
            matched = true;
        }
        debugPrint(ctx.isDebug, "JER match: hasGen=1, dR=" + std::to_string(dR) +
                                ", matched=" + std::to_string(matched));
    } else {
        debugPrint(ctx.isDebug, "JER match: hasGen=0");
    }

    if (matched) {
        const double corr = std::max(0.0, 1.0 + (sf - 1.0) * (jAfterJes.pt - jer.genPt) / jAfterJes.pt);
        debugPrint(ctx.isDebug, "JER matched corr=" + std::to_string(corr));
        return corr;
    } else {
        const double gaus = ctx.refs.randomGen.Gaus(0.0, reso);
        const double varp = std::max(sf * sf - 1.0, 0.0);
        const double corr = std::max(0.0, 1.0 + gaus * std::sqrt(varp));
        debugPrint(ctx.isDebug, "JER stochastic: gaus=" + std::to_string(gaus) +
                                ", varp=" + std::to_string(varp) +
                                ", corr=" + std::to_string(corr));
        return corr;
    }
}

// ======================================================
// 3) Single JES component systematic shift (factor)
// pt_syst = pt_nominal * factor
// ======================================================
inline double getJesCorrectionSyst(JecConfigReader::CorrectionRefs& refs,
                              const std::string& systName,  // key in CorrectionSet
                              const std::string& var,       // "Up" | "Down"
                              double eta,
                              double ptAfterJes,
                              bool isDebug = false) {
    auto systCorr = JecConfigReader::safeGet(refs.cs, systName);
    const double scale = systCorr->evaluate({ eta, ptAfterJes });
    const double factor = (var == "Up") ? (1.0 + scale) : (1.0 - scale);
    debugPrint(isDebug, "JES syst: name=" + systName +
                        ", var=" + var +
                        ", scale=" + std::to_string(scale) +
                        ", factor=" + std::to_string(factor));
    return factor;
}

// ======================================================
// 4) Type-1 MET correction
//
// Recomputes per-jet L1→L2→Residual(+JES/JER) on muon-subtracted RAW jet pt,
// then shifts MET by (full − L1) for jets passing the standard Type-1 selection.
//
// Returns the corrected MET p4.
// ======================================================
inline TLorentzVector getCorrectedMet(const EvalContext& ctx,
                                       const MetInputs& met,
                                       const std::vector<JetForMet>& jets,
                                       const std::vector<JerInputs>& jers,  // same size as jets
                                       double rho,
                                       const SystematicOptions& opts = {}) {
    double metPx = met.metPt * std::cos(met.metPhi);
    double metPy = met.metPt * std::sin(met.metPhi);
    debugPrint(ctx.isDebug, "MET start: pt=" + std::to_string(met.metPt) +
                            ", phi=" + std::to_string(met.metPhi));

    if (!ctx.isData && jers.size() != jets.size()) {
        debugPrint(true, "[WARNING] jers.size() != jets.size(); JER will index-assume.");
    }

    for (std::size_t i = 0; i < jets.size(); ++i) {
        const auto& v = jets[i];

        // Muon-subtracted raw pt
        const double ptRawMinusMuon = v.rawPt * (1.0 - v.muonSubtrFactor);
        debugPrint(ctx.isDebug, "Jet[" + std::to_string(i) + "] rawMinusMuon=" + std::to_string(ptRawMinusMuon));

        // --- L1 ---
        double ptCorr = ptRawMinusMuon;
        const double c1 = ctx.refs.corrRefJesL1FastJet->evaluate({ v.area, v.eta, ptCorr, rho });
        ptCorr *= c1;
        const double ptCorrL1 = ptCorr;  // Save for Type-1 delta
        debugPrint(ctx.isDebug, "  L1 c1=" + std::to_string(c1) +
                                ", ptCorrL1=" + std::to_string(ptCorrL1));

        // --- L2 ---
        double c2 = 1.0;
        if (hasPhiDependentL2(ctx.year)) {
            c2 = ctx.refs.corrRefJesL2Relative->evaluate({ v.eta, v.phi, ptCorr });
        } else {
            c2 = ctx.refs.corrRefJesL2Relative->evaluate({ v.eta, ptCorr });
        }
        ptCorr *= c2;
        debugPrint(ctx.isDebug, "  L2 c2=" + std::to_string(c2) +
                                ", ptAfterL2=" + std::to_string(ptCorr));

        // --- Residuals for data ---
        if (ctx.isData) {
            double cRes = 1.0;
            if (requiresRunBasedResidual(ctx.year)) {
                const double run = ctx.runNumber.value_or(0.0);
                cRes = ctx.refs.corrRefJesL2ResL3Res->evaluate({ run, v.eta, ptCorr });
            } else {
                cRes = ctx.refs.corrRefJesL2ResL3Res->evaluate({ v.eta, ptCorr });
            }
            ptCorr *= cRes;
            debugPrint(ctx.isDebug, "  Residuals cRes=" + std::to_string(cRes) +
                                    ", ptAfterRes=" + std::to_string(ptCorr));
        }

        // --- JES uncertainty (MC only) ---
        if (!ctx.isData && !opts.jesSystName.empty()) {
            const double jesShift = getJesCorrectionSyst(ctx.refs, opts.jesSystName, opts.jesSystVar, v.eta, ptCorr, ctx.isDebug);
            ptCorr *= jesShift;
            debugPrint(ctx.isDebug, "  JES syst factor=" + std::to_string(jesShift) +
                                    ", ptAfterJESsyst=" + std::to_string(ptCorr));
        }

        // --- JER (MC only) ---
        if (!ctx.isData) {
            JesInputs jAfter{};
            jAfter.pt        = ptCorr;  // already includes JES (+ JES syst)
            jAfter.eta       = v.eta;
            jAfter.phi       = v.phi;
            jAfter.area      = v.area;
            jAfter.rho       = rho;
            jAfter.rawFactor = 0.0;

            const JerInputs& jer = (i < jers.size()) ? jers[i] : JerInputs{};
            const double sJer = getJerCorrectionNomOrSyst(ctx, jAfter, jer, opts);
            ptCorr *= sJer;
            debugPrint(ctx.isDebug, "  JER factor=" + std::to_string(sJer) +
                                    ", ptAfterJER=" + std::to_string(ptCorr));
        }

        // --- Type-1 selection ---
        const bool passSel = (ptCorr > 15.0 &&
                              std::abs(v.eta) < 5.2 &&
                              (v.chEmEf + v.neEmEf) < 0.9);
        debugPrint(ctx.isDebug, "  passType1Sel=" + std::string(passSel ? "true" : "false"));

        if (!passSel) continue;

        // Apply (full − L1) to MET
        const double dpt = (ptCorr - ptCorrL1);
        metPx -= dpt * std::cos(v.phi);
        metPy -= dpt * std::sin(v.phi);
        debugPrint(ctx.isDebug, "  dpt=" + std::to_string(dpt) +
                                ", metPx=" + std::to_string(metPx) +
                                ", metPy=" + std::to_string(metPy));
    }

    const double metPt  = std::hypot(metPx, metPy);
    const double metPhi = std::atan2(metPy, metPx);
    debugPrint(ctx.isDebug, "MET final: pt=" + std::to_string(metPt) +
                            ", phi=" + std::to_string(metPhi));

    TLorentzVector p4;
    p4.SetPtEtaPhiM(metPt, 0.0, metPhi, 0.0);
    return p4;
}

} // namespace JecApplication

