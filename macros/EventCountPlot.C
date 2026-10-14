// EventCountPlot.C
//
// For each file in the list below, read the "EventCount" tree and make two 2D
// plots, with z = sum of Nevent over all entries with that mass pair:
//   1) parent mass (MP) vs LSP mass (MC)
//   2) parent mass (MP) vs Delta M = MP - MC
// A CSV of events per mass point is also written for each file.
//
// Both axes are true numeric axes (no merging of nearby masses). Bin width on
// each axis is the smallest spacing between distinct mass values, so every
// mass point sits alone in its own bin.
//
// Run:  root -l -b -q macros/EventCountPlot.C
// Or:   root -l -b -q 'macros/EventCountPlot.C("root/EventCount")'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

namespace {

// ---------------------------------------------------------------------------
// Configuration -- edit these
// ---------------------------------------------------------------------------

// Name of the tree inside each file
const std::string kTreeName = "EventCount";

// Branch names inside the EventCount tree
const std::string kXMassBranch = "MP";      // parent mass
const std::string kYMassBranch = "MC";      // LSP (child) mass
const std::string kCountBranch = "Nevent";  // events for that mass pair

// Masses are stored as integers in units of 1/100 GeV (e.g. 30000 -> 300 GeV).
// Set to 1 if the branch is already in GeV.
const double kMassScale = 100.0;

// Smallest meaningful mass step in GeV. All masses and delta-M values are
// re-quantized to this so floating-point noise (0.6000000001 vs 0.5999999999)
// can't create duplicate axis values.
const double kMassQuantum = 0.01;

// Snap parent mass to the nearest multiple of this many GeV before computing
// delta-M, to absorb off-by-one rounding in the inputs (526 -> 525). 0 = off.
const double kParentSnap = 25.0;

// Straggler delta-M cleanup: a delta-M value that appears in fewer than this
// fraction of the parent-mass rows is treated as a rounding artifact and moved
// to the nearest "common" delta-M value within kDMStragglerTol GeV (ties go up).
// Set kDMStragglerTol = 0 to turn off. Assumes real delta-M points exist for
// most parent masses; lower the fraction if a real delta-M only exists for a few.
const double kDMStragglerFraction = 0.5;
const double kDMStragglerTol      = 1.0;

// Mass cuts (in GeV, applied after dividing by kMassScale)
const double kMinParentMass = 100.0;
const double kMaxParentMass = 650.0;
const double kMinLSPMass    = 0.0;
const double kMaxLSPMass    = 1e9;
const double kMinDeltaM     = 0.0;
const double kMaxDeltaM     = 50.0;

// Z-axis: number of color steps, fixed range for comparing samples (0 = auto)
const int    kNContours = 100;
const double kZMin = 0.0;
const double kZMax = 0.0;

// Axis titles
const std::string kXTitle  = "Parent mass [GeV]";
const std::string kYTitle  = "LSP mass [GeV]";
const std::string kDMTitle = "Parent mass - LSP mass [GeV]";
const std::string kZTitle  = "N_{events}";

// Safety cap on bins per axis (warns if the spacing would need more)
const int kMaxBinsPerAxis = 3000;

// Output directory for plots / root file / csv
const std::string kOutDir = "plots";

// Files to process, relative to the input directory.
// {label used in output names / plot title, path relative to input dir}
const std::vector<std::pair<std::string, std::string>> kFiles = {
    {"TChiWZ_Autumn18",
     "Autumn18_102X_SMS/TChiWZ_genHT-160_genMET-80_TuneCP2_13TeV-madgraphMLM-pythia8_Autumn18_102X.root"},
    {"TChiWZ_Summer22",
     "Summer22_130X_SMS/SMS-TChiWZ_genHT-50_genMET-80_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer22_130X.root"},
    // add more:
    // {"TChiWW_Autumn18", "Autumn18_102X_SMS/TChiWW_....root"},
};

// ---------------------------------------------------------------------------

using PairSums = std::map<std::pair<double, double>, double>;

void DrawLabels(const std::string& title)
{
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.04);
    l.DrawLatex(0.15, 0.943, "#bf{CMS} Simulation Preliminary");
    l.SetTextSize(0.035);
    l.DrawLatex(0.62, 0.943, title.c_str());
}

// Choose uniform bin edges so each distinct value gets its own bin, centered
// on the value. Bin width = smallest gap between distinct sorted values.
struct Axis {
    int nbins = 1;
    double lo = 0.0;
    double hi = 1.0;
};

Axis MakeAxis(const std::set<double>& values)
{
    Axis ax;
    const std::vector<double> v(values.begin(), values.end());
    double w = 1.0;
    if (v.size() > 1) {
        w = std::numeric_limits<double>::max();
        for (size_t i = 1; i < v.size(); ++i) {
            const double gap = v[i] - v[i - 1];
            if (gap > 1e-6) w = std::min(w, gap);
        }
    }
    ax.nbins = static_cast<int>(std::lround((v.back() - v.front()) / w)) + 1;
    if (ax.nbins > kMaxBinsPerAxis) {
        std::cerr << "Warning: axis needs " << ax.nbins << " bins at width " << w
                  << " GeV; capping at " << kMaxBinsPerAxis
                  << " (some mass points will share bins)" << std::endl;
        ax.nbins = kMaxBinsPerAxis;
        w = (v.back() - v.front()) / (ax.nbins - 1);
    }
    ax.lo = v.front() - 0.5 * w;
    ax.hi = ax.lo + ax.nbins * w;
    return ax;
}

void DrawMap(const std::string& name, const std::string& label, const PairSums& sums,
             const std::string& xTitle, const std::string& yTitle, TFile& outFile)
{
    std::set<double> xs, ys;
    for (const auto& [m, c] : sums) {
        xs.insert(m.first);
        ys.insert(m.second);
    }
    const Axis ax = MakeAxis(xs);
    const Axis ay = MakeAxis(ys);

    TH2D hist(name.c_str(), (";" + xTitle + ";" + yTitle + ";" + kZTitle).c_str(),
              ax.nbins, ax.lo, ax.hi, ay.nbins, ay.lo, ay.hi);

    double total = 0.0;
    double minPos = std::numeric_limits<double>::max();
    for (const auto& [m, c] : sums) {
        hist.Fill(m.first, m.second, c);
        total += c;
        if (c > 0) minPos = std::min(minPos, c);
    }
    std::cout << "  " << name << ": " << sums.size() << " mass points, total Nevent = "
              << total << ", bins " << ax.nbins << " x " << ay.nbins << std::endl;
    if (minPos < std::numeric_limits<double>::max()) hist.SetMinimum(0.5 * minPos);

    hist.GetXaxis()->CenterTitle();
    hist.GetYaxis()->CenterTitle();
    hist.GetZaxis()->CenterTitle();
    hist.GetZaxis()->SetTitleOffset(1.65);
    hist.GetZaxis()->SetMoreLogLabels();
    hist.GetZaxis()->SetNoExponent();
    if (kZMax > 0.0) {
        hist.SetMinimum(kZMin > 0.0 ? kZMin : 0.5 * minPos);
        hist.SetMaximum(kZMax);
    }
    gStyle->SetNumberContours(kNContours);

    TCanvas can(("can_" + name).c_str(), name.c_str(), 1000, 800);
    can.SetLeftMargin(0.15);
    can.SetRightMargin(0.18);
    can.SetBottomMargin(0.15);
    can.SetLogz();
    hist.Draw("COLZ");
    DrawLabels(label);

    outFile.cd();
    hist.Write();
    can.Write();
    can.SaveAs((kOutDir + "/" + name + ".pdf").c_str());
    can.SaveAs((kOutDir + "/" + name + ".png").c_str());
}

double Canon(double v) { return std::round(v / kMassQuantum) * kMassQuantum; }

// Fold rare ("straggler") delta-M values onto the nearest common delta-M value.
void SnapStragglerDM(PairSums& dmSums)
{
    if (kDMStragglerTol <= 0.0) return;

    std::map<double, std::set<double>> rowsPerDM;
    std::set<double> allRows;
    for (const auto& [m, c] : dmSums) {
        rowsPerDM[m.second].insert(m.first);
        allRows.insert(m.first);
    }
    std::vector<double> common;  // ascending (map order)
    for (const auto& [dm, rows] : rowsPerDM)
        if (rows.size() >= kDMStragglerFraction * allRows.size()) common.push_back(dm);
    if (common.empty()) return;

    PairSums fixedSums;
    int moved = 0;
    for (const auto& [m, c] : dmSums) {
        double dm = m.second;
        if (!std::binary_search(common.begin(), common.end(), dm)) {
            double best = dm, bestDist = std::numeric_limits<double>::max();
            for (double cm : common) {
                const double d = std::abs(cm - dm);
                // "<=" with ascending order makes ties go to the larger value
                if (d <= kDMStragglerTol + 1e-9 && d <= bestDist + 1e-9) {
                    best = cm;
                    bestDist = d;
                }
            }
            if (best != dm) ++moved;
            dm = best;
        }
        fixedSums[{m.first, dm}] += c;
    }
    std::cout << "  snapped " << moved << " straggler delta-M entries onto "
              << common.size() << " common delta-M values" << std::endl;
    dmSums = std::move(fixedSums);
}

std::string FormatMass(double m)
{
    std::ostringstream s;
    s << m;  // default formatting drops trailing zeros (100, 12.5, ...)
    return s.str();
}

// Draw a pivot table image: one row per distinct first-key value (parent mass),
// one column per distinct second-key value, cell = summed Nevent (blank if no
// such mass point exists).
void DrawTable(const std::string& name, const std::string& label, const PairSums& sums,
               const std::string& rowTitle, const std::string& colTitle, TFile& outFile)
{
    std::set<double> rowSet, colSet;
    for (const auto& [m, c] : sums) {
        rowSet.insert(m.first);
        colSet.insert(m.second);
    }
    const std::vector<double> rows(rowSet.begin(), rowSet.end());
    const std::vector<double> cols(colSet.begin(), colSet.end());
    const int nR = rows.size();
    const int nC = cols.size();

    // Layout in pixels
    const double cellW = 95, cellH = 30;
    const double left = 150, right = 20, top = 110, bottom = 20;
    const double W = left + nC * cellW + right;
    const double H = top + nR * cellH + bottom;

    TCanvas can(("can_" + name).c_str(), name.c_str(), static_cast<int>(W), static_cast<int>(H));
    can.SetMargin(0, 0, 0, 0);
    can.SetFillColor(kWhite);

    const auto nx = [&](double px) { return px / W; };
    const auto ny = [&](double py) { return 1.0 - py / H; };

    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);

    // Title strip
    t.SetTextAlign(11);
    t.SetTextSize(0.45 * cellH / H);
    t.DrawLatex(nx(10), ny(25), "#bf{CMS} Simulation Preliminary");
    t.SetTextAlign(31);
    t.DrawLatex(nx(W - right), ny(25), label.c_str());
    t.SetTextAlign(21);
    t.SetTextSize(0.5 * cellH / H);
    t.DrawLatex(nx(left + 0.5 * nC * cellW), ny(top - cellH - 12), (colTitle + "  (cells: N_{events})").c_str());

    // Header row (column values) and header column (row values)
    t.SetTextSize(0.55 * cellH / H);
    t.SetTextAlign(22);
    for (int j = 0; j < nC; ++j)
        t.DrawLatex(nx(left + (j + 0.5) * cellW), ny(top - 0.5 * cellH), ("#bf{" + FormatMass(cols[j]) + "}").c_str());
    for (int i = 0; i < nR; ++i)
        t.DrawLatex(nx(0.5 * left), ny(top + (i + 0.5) * cellH), ("#bf{" + FormatMass(rows[i]) + "}").c_str());
    t.SetTextSize(0.4 * cellH / H);
    t.DrawLatex(nx(0.5 * left), ny(top - 0.5 * cellH), rowTitle.c_str());

    // Cell contents
    t.SetTextSize(0.5 * cellH / H);
    for (const auto& [m, c] : sums) {
        const int i = std::lower_bound(rows.begin(), rows.end(), m.first) - rows.begin();
        const int j = std::lower_bound(cols.begin(), cols.end(), m.second) - cols.begin();
        t.DrawLatex(nx(left + (j + 0.5) * cellW), ny(top + (i + 0.5) * cellH),
                    std::to_string(static_cast<long long>(std::llround(c))).c_str());
    }

    // Grid lines (include header row/column)
    TLine line;
    line.SetNDC();
    line.SetLineColor(kGray + 1);
    for (int i = -1; i <= nR; ++i)
        line.DrawLine(nx(0), ny(top + i * cellH), nx(left + nC * cellW), ny(top + i * cellH));
    for (int j = -1; j <= nC; ++j) {
        const double px = (j < 0) ? 0 : left + j * cellW;
        line.DrawLine(nx(px), ny(top - cellH), nx(px), ny(top + nR * cellH));
    }

    outFile.cd();
    can.Write();
    can.SaveAs((kOutDir + "/" + name + ".pdf").c_str());
    can.SaveAs((kOutDir + "/" + name + ".png").c_str());
}

bool MakePlots(const std::string& label, const std::string& path, TFile& outFile)
{
    std::unique_ptr<TFile> in(TFile::Open(path.c_str(), "READ"));
    if (!in || in->IsZombie()) {
        std::cerr << "Failed to open " << path << std::endl;
        return false;
    }
    auto* tree = dynamic_cast<TTree*>(in->Get(kTreeName.c_str()));
    if (!tree) {
        std::cerr << "No tree '" << kTreeName << "' in " << path << std::endl;
        return false;
    }

    TLeaf* lx = tree->GetLeaf(kXMassBranch.c_str());
    TLeaf* ly = tree->GetLeaf(kYMassBranch.c_str());
    TLeaf* ln = tree->GetLeaf(kCountBranch.c_str());
    if (!lx || !ly || !ln) {
        std::cerr << "Missing branch in " << path << " (need " << kXMassBranch << ", "
                  << kYMassBranch << ", " << kCountBranch << "). Available:" << std::endl;
        tree->Print();
        return false;
    }

    // Sum Nevent per (parent, LSP) pair and per (parent, DeltaM) pair
    PairSums lspSums, dmSums;
    const Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        double x = lx->GetValue() / kMassScale;
        double y = ly->GetValue() / kMassScale;
        const double c = ln->GetValue();
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(c)) continue;
        const auto canon = [](double v) { return std::round(v / kMassQuantum) * kMassQuantum; };
        if (kParentSnap > 0.0) x = std::round(x / kParentSnap) * kParentSnap;
        x = canon(x);
        y = canon(y);
        const double dm = canon(x - y);
        if (x < kMinParentMass || x > kMaxParentMass) continue;
        if (y < kMinLSPMass    || y > kMaxLSPMass)    continue;
        if (dm < kMinDeltaM    || dm > kMaxDeltaM)    continue;
        dmSums[{x, dm}] += c;
    }
    if (dmSums.empty()) {
        std::cerr << "No usable entries in " << path << std::endl;
        return false;
    }

    std::cout << label << ":" << std::endl;
    SnapStragglerDM(dmSums);
    // Rebuild the LSP-mass sums from the cleaned delta-M so every output agrees
    for (const auto& [m, c] : dmSums)
        lspSums[{m.first, Canon(m.first - m.second)}] += c;

    // Dump events per mass point to a CSV so samples can be compared directly
    {
        std::ofstream csv(kOutDir + "/" + label + "_EventCount.csv");
        csv << "parent_mass,lsp_mass,delta_m,nevent\n";
        for (const auto& [m, c] : lspSums)
            csv << m.first << "," << m.second << "," << (m.first - m.second) << "," << c << "\n";
    }

    DrawMap(label + "_EventCount",    label, lspSums, kXTitle, kYTitle,  outFile);
    DrawMap(label + "_dm_EventCount", label, dmSums,  kXTitle, kDMTitle, outFile);
    DrawTable(label + "_EventCount_table",    label, lspSums, "Parent (rows)", kYTitle,  outFile);
    DrawTable(label + "_dm_EventCount_table", label, dmSums,  "Parent (rows)", kDMTitle, outFile);
    return true;
}

} // namespace

void EventCountPlot(const std::string& inputDir = "root/EventCount")
{
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gSystem->mkdir(kOutDir.c_str(), true);

    TFile outFile((kOutDir + "/output_EventCountPlot.root").c_str(), "RECREATE");
    int ok = 0;
    for (const auto& [label, relPath] : kFiles) {
        const std::string path = inputDir + "/" + relPath;
        std::cout << "Processing " << path << std::endl;
        if (MakePlots(label, path, outFile)) ++ok;
    }
    outFile.Close();
    std::cout << "Processed " << ok << "/" << kFiles.size() << " files; output in "
              << kOutDir << "/" << std::endl;
}