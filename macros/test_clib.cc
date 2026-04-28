R__ADD_INCLUDE_PATH(/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-4f091fd2adcff55f05f1d262b3254c25/lib/python3.9/site-packages/correctionlib/include)
R__LOAD_LIBRARY(/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-4f091fd2adcff55f05f1d262b3254c25/lib/python3.9/site-packages/correctionlib/lib/libcorrectionlib.so)

#include "correction.h"
#include <iostream>
#include <map>
#include <vector>
#include <string>

// a file to test the output of various pt etas and make sure the correctionlib file is loaded correctly

double lookup_sf(double pt, double eta, std::string year, std::string era, std::string sf_key){
    auto cset = correction::CorrectionSet::from_file("../root/LepSF/LepSFs_fixed.json");
    auto sf_map = cset->at(sf_key);
    double sf = sf_map->evaluate({pt, std::abs(eta), year, era});
    return sf;
}

void test_clib(){

    std::cout << std::fixed << std::setprecision(8);
    std::string sf_key = "electron_BLP_over_COL";

    std::map<std::string, std::vector<std::string>> year_eras = {
        {"2016", {"preVFP", "postVFP"}},
        {"2017", {"none"}},
        {"2018", {"none"}},
        {"2022", {"preEE", "postEE"}},
        {"2023", {"preBPix", "postBpix"}},
    };

    // test points: {pt, eta, label}
    std::vector<std::tuple<double,double,std::string>> test_points = {
        {2.0,   0.0,    "corner pt=2.0,   eta=0.0   "},
        {1001.0, 0.0,    "corner pt=1001.0, eta=0.0   "},
        {1001.0, 2.4,    "corner pt=1001.0, eta=2.4   "},
        {5.0,   0.4,    "edge   pt=5.0,   eta=0.4   "},
        {7.0,   0.4,    "edge   pt=7.0,   eta=0.4   "},
        {10.0,  0.4,    "edge   pt=10.0,  eta=0.4   "},
        {20.0,  0.4,    "edge   pt=20.0,  eta=0.4   "},
        {45.0,  0.4,    "edge   pt=45.0,  eta=0.4   "},
        {75.0,  0.4,    "edge   pt=75.0,  eta=0.4   "},
        {30.0,  0.8,    "edge   pt=30.0,  eta=0.8   "},
        {30.0,  1.4442, "edge   pt=30.0,  eta=1.4442"},
        {30.0,  2.5,    "edge   pt=30.0,  eta=2.5   "},
        {10.0,  1.0,    "mid    pt=10.0,  eta=1.0   "},
        {50.0,  2.0,    "mid    pt=50.0,  eta=2.0   "},
    };

    for(auto& [year, eras] : year_eras){
        for(auto& era : eras){
            std::cout << "\n=== " << year << " " << era << " ===\n";
            for(auto& [pt, eta, label] : test_points){
                std::cout << label << "  SF=" << lookup_sf(pt, eta, year, era, sf_key) << "\n";
            }
        }
    }
}