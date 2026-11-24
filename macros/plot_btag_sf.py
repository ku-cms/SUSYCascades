#!/usr/bin/env python3
"""Plot b-tag numerator/denominator ratios from ROOT histogram files."""

import argparse
from pathlib import Path
import re
import tarfile

import ROOT


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = [
    Path(
        "root/BtagSF/Summer23BPix_130X_Cascades_HISTS/"
        "SlepSnuCascade_MN1-270_MN2-280_MC1-275_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer23BPix_130X.root"
    ),
    Path(
        "root/BtagSF/Summer24_130X_SMS_HISTS/"
        "SMS-TChiWZ_Bin-genHT-50-genMET-80_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24_130X.root"
    ),
    Path(
        "root/BtagSF/Summer24_130X_HISTS/"
        "TTto2L2Nu-3Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8_Summer24_130X.root"
    ),
]

# Change this to another SMS mass-point fragment, or to "" to plot all pairs.
SMS_HISTOGRAM_FILTER = "SMS_30000_25000"
FLAVORS = (0, 1, 2)
HISTOGRAM_DIRECTORY = "Histograms"


def find_histogram_pairs(histogram_directory, name_filter, flavor):
    numerator_suffix = f"_flavor{flavor}_num"
    denominator_suffix = f"_flavor{flavor}_den"
    pairs = {}
    for key in histogram_directory.GetListOfKeys():
        name = key.GetName()
        if key.GetClassName() != "TH1D" or name_filter not in name:
            continue
        if name.endswith(numerator_suffix):
            pairs.setdefault(name[: -len(numerator_suffix)], {})["num"] = name
        elif name.endswith(denominator_suffix):
            pairs.setdefault(name[: -len(denominator_suffix)], {})["den"] = name

    incomplete = [
        stem for stem, pair in pairs.items() if "num" not in pair or "den" not in pair
    ]
    if incomplete:
        raise RuntimeError(
            "Histogram numerator/denominator pair is incomplete: "
            + ", ".join(incomplete)
        )
    if not pairs:
        raise RuntimeError(
            "No TH1D histogram pairs ending in "
            f"{numerator_suffix!r} and {denominator_suffix!r} "
            f"matched filter {name_filter!r}"
        )
    return pairs


def histogram_tag(stem, flavor):
    mass_point = re.search(r"(SMS_[^_]+_[^_]+)$", stem)
    if mass_point:
        return f"{mass_point.group(1)}_flavor{flavor}"
    return f"flavor{flavor}"


def make_ratio_plot(numerator, denominator, title, output_path):
    if numerator.GetNbinsX() != denominator.GetNbinsX():
        raise RuntimeError(f"Histogram bin counts differ for {title}")
    for bin_index in range(1, numerator.GetNbinsX() + 2):
        numerator_edge = numerator.GetXaxis().GetBinLowEdge(bin_index)
        denominator_edge = denominator.GetXaxis().GetBinLowEdge(bin_index)
        if numerator_edge != denominator_edge:
            raise RuntimeError(f"Histogram bin edges differ for {title}")

    ratio = numerator.Clone("btag_efficiency")
    ratio.SetDirectory(0)
    if not ratio.Divide(numerator, denominator, 1.0, 1.0, "B"):
        raise RuntimeError(f"Could not divide numerator by denominator for {title}")
    ratio.SetTitle(title)
    ratio.GetXaxis().SetTitle("Jet p_{T} [GeV]")
    ratio.GetYaxis().SetTitle("Numerator / denominator")
    ratio.SetMinimum(0.0)
    ratio.SetMaximum(1.05)
    ratio.SetLineColor(ROOT.kBlue + 1)
    ratio.SetMarkerColor(ROOT.kBlue + 1)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(0.8)

    canvas = ROOT.TCanvas("btag_ratio_canvas", title, 800, 600)
    canvas.SetGrid()
    ratio.Draw("E1")
    canvas.SaveAs(str(output_path))
    canvas.Close()
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Failed to write plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "root/BtagSF/RatioPlots",
        help="directory for generated PDF plots",
    )
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for relative_path in ROOT_FILES:
        input_path = REPOSITORY_ROOT / relative_path
        if not input_path.is_file():
            raise FileNotFoundError(f"ROOT input file does not exist: {input_path}")

        root_file = ROOT.TFile.Open(str(input_path))
        if not root_file or root_file.IsZombie():
            raise OSError(f"Could not open ROOT input file: {input_path}")
        try:
            histogram_directory = root_file.GetDirectory(HISTOGRAM_DIRECTORY)
            if not histogram_directory:
                raise RuntimeError(
                    f"Missing {HISTOGRAM_DIRECTORY!r} directory in {input_path}"
                )

            name_filter = (
                SMS_HISTOGRAM_FILTER
                if "_SMS_HISTS/" in relative_path.as_posix()
                else ""
            )
            sample_output_dir = args.output_dir / input_path.stem
            sample_output_dir.mkdir(parents=True, exist_ok=True)

            for flavor in FLAVORS:
                pairs = find_histogram_pairs(
                    histogram_directory, name_filter, flavor
                )
                for stem, pair in pairs.items():
                    numerator = histogram_directory.Get(pair["num"])
                    denominator = histogram_directory.Get(pair["den"])
                    if not numerator or not denominator:
                        raise RuntimeError(f"Could not read histogram pair {stem!r}")

                    tag = histogram_tag(stem, flavor)
                    output_path = sample_output_dir / f"{tag}_ratio.pdf"
                    make_ratio_plot(
                        numerator,
                        denominator,
                        f"{input_path.stem}: {tag}",
                        output_path,
                    )
                    print(f"Wrote {output_path}")
        finally:
            root_file.Close()

    archive_path = args.output_dir.parent / f"{args.output_dir.name}.tar.gz"
    with tarfile.open(archive_path, "w:gz") as archive:
        archive.add(args.output_dir, arcname=args.output_dir.name)
    print(f"Wrote {archive_path}")


if __name__ == "__main__":
    main()
