#include "GitHashTool.hh"

#include <TFile.h>
#include <TParameter.h>

#include <cstdio>
#include <stdexcept>
#include <array>
#include <memory>

std::string GitHashTool::GetGitCommitHash() const {
    std::array<char, 128> buffer;
    std::string result;

    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen("git rev-parse HEAD", "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("Failed to run git command");
    }

    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }

    if (!result.empty() && result.back() == '\n') {
        result.pop_back();  // remove trailing newline
    }

    return result;
}

void GitHashTool::WriteToROOT(const std::string& outputFileName) const {
    std::string gitHash = GetGitCommitHash();

    TFile *f = new TFile(outputFileName.c_str(), "UPDATE");
    if (!f || f->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + outputFileName);
    }

    TNamed param("GitCommitHash", gitHash);
    param.Write();

    f->Close();
    delete f;
}

