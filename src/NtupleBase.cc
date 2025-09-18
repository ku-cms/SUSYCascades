#include <TFile.h>
#include <TError.h>
#include <TParameter.h>

#include "NtupleBase.hh"
#include "SUSYNANOBase.hh"
#include "NANOULBase.hh"
#include "NANORun3.hh"
#include "NeventTool.hh"

template <class Base>
NtupleBase<Base>::NtupleBase(TTree* tree)
  : AnalysisBase<Base>(tree)
{
  
}

template <class Base>
NtupleBase<Base>::~NtupleBase(){
  int Ntree = m_Trees.size();
  for(int i = 0; i < Ntree; i++)
    if(m_Trees[i])
      delete m_Trees[i];
}

template <class Base>
void NtupleBase<Base>::GetChunks(const Long64_t& u_NTOT, Long64_t& N0, Long64_t& N1, int ichunk, int nchunk){
    if(u_NTOT == 0) return;
    Long64_t base = u_NTOT / nchunk;
    Long64_t rem  = u_NTOT % nchunk;
    if(ichunk <= rem){
        N0 = (ichunk - 1) * (base + 1);
        N1 = N0 + (base + 1);
    } else {
        N0 = rem * (base + 1) + (ichunk - rem - 1) * base;
        N1 = N0 + base;
    }
    if(N1 > u_NTOT) N1 = u_NTOT;
}

template <class Base>
bool NtupleBase<Base>::WriteNtuple(const string& filename, int ichunk, int nchunk, bool do_slim, int NDAS, const string& DAS_datasetname, const string& DAS_filename){
  TFile* outfile = new TFile(filename.c_str(),"RECREATE");
  outfile->cd();

  string sample;

  std::pair<int,int> masses(0,0);

  if(nchunk < 1 || ichunk < 1 || ichunk > nchunk){
    ichunk = 1;
    nchunk = 1;
  }

  std::cout << "Getting entries in tree" << std::endl;
  gErrorIgnoreLevel = kFatal;
  Long64_t NTOT = Base::fChain->GetEntries();
  gErrorIgnoreLevel = 0;
  cout << "NTOT: " << NTOT << endl;
  Long64_t N1, N0;
  GetChunks(NTOT, N0, N1, ichunk, nchunk);
  // Initialize Histogram Booking
  vector<TH1D*> histos;
  if(!AnalysisBase<Base>::IsData())
    AnalysisBase<Base>::InitializeHistograms(histos);

  int Nsys = AnalysisBase<Base>::m_Systematics.GetN();

  if(Nsys > 1)
    do_slim = true;
  
  cout << "looping between " << N0 << " " << N1 << endl;
  for(Long64_t i = N0; i < N1 && i < NTOT; i++){
    int mymod = (N1-N0)/10;
    if(mymod < 1)
      mymod = 1;
    if(i%mymod == 0)
      cout << " event = " << i << " : [" << N0 << " , " << N1 << "]" << endl;
    
    sample = AnalysisBase<Base>::GetEntry(i);
    
    if(m_Label2Tree.count(sample) == 0){
      m_Label2Tree[sample] = std::vector<TTree*>();
      
      for(int s = 0; s < Nsys; s++){
        if(AnalysisBase<Base>::m_Systematics[s].Label().compare("METUncer_GenMET") == 0){
	  if(!AnalysisBase<Base>::IsFastSim()) continue;
	  TTree* tree = InitOutputTree(AnalysisBase<Base>::m_Systematics[s].Up().TreeName(sample), do_slim);
	  m_Label2Tree[sample].push_back(tree);
	  m_Trees.push_back(tree);
        }
        else{
	  TTree* tree = InitOutputTree(AnalysisBase<Base>::m_Systematics[s].Up().TreeName(sample), do_slim);
	  m_Label2Tree[sample].push_back(tree);
	  m_Trees.push_back(tree);
	  if(!(!AnalysisBase<Base>::m_Systematics[s])){
	    tree = InitOutputTree(AnalysisBase<Base>::m_Systematics[s].Down().TreeName(sample), do_slim);
	    m_Label2Tree[sample].push_back(tree);
	    m_Trees.push_back(tree);
	  }
        }
      }
    }
    
    outfile->cd();
    int isys = 0;
    for(int s = 0; s < Nsys; s++){
      if(AnalysisBase<Base>::m_Systematics[s].Label().compare("METUncer_GenMET") == 0 && !AnalysisBase<Base>::IsFastSim())
	continue;

      FillOutputTree(m_Label2Tree[sample][isys], AnalysisBase<Base>::m_Systematics[s].Up(), do_slim);
      isys++;
      if(!(!AnalysisBase<Base>::m_Systematics[s]) && AnalysisBase<Base>::m_Systematics[s].Label().compare("METUncer_GenMET") != 0){
	FillOutputTree(m_Label2Tree[sample][isys], AnalysisBase<Base>::m_Systematics[s].Down(), do_slim);
	isys++;
      }
    }

    if(!AnalysisBase<Base>::IsData())
      AnalysisBase<Base>::BookHistograms(histos);

    // event count bookkeeping
    if(AnalysisBase<Base>::IsSMS())
      masses = AnalysisBase<Base>::GetSUSYMasses();
    if(m_mapNevent.count(masses) == 0){
      m_masses.push_back(masses);
      m_mapNevent[masses] = 1.;
    } else {
      m_mapNevent[masses] += 1.;
    }
  }

  int Ntree = m_Trees.size();
  cout << "Ntree " << Ntree << endl;
  for(int i = 0; i < Ntree; i++){
    outfile->cd();
    m_Trees[i]->Write("",TObject::kOverwrite);
    delete m_Trees[i];
    m_Trees[i] = nullptr;
  }

  // event count tree
  outfile->cd();
  TTree* tout = (TTree*) new TTree("EventCount", "EventCount");
  
  string dataset = string(AnalysisBase<Base>::GetDataSet());
  string filetag = string(AnalysisBase<Base>::GetFileTag());
  double Nevent;
  double Nweight = 0.;
  int Nevent_tot = 0;
  int MP;
  int MC;
  std::string saved_DAS_datasetname = DAS_datasetname;
  std::string saved_DAS_filename = DAS_filename;
  tout->Branch("NDAS", &NDAS);
  tout->Branch("Nevent", &Nevent);
  tout->Branch("Nweight", &Nweight);
  tout->Branch("filetag", &filetag);
  tout->Branch("dataset", &dataset);
  tout->Branch("DAS_datasetname", &saved_DAS_datasetname);
  tout->Branch("DAS_filename", &saved_DAS_filename);
  tout->Branch("MP", &MP);
  tout->Branch("MC", &MC);
  int Nmass = m_masses.size();
  // Loop for DAS check
  for(int i = 0; i < Nmass; i++){
    Nevent_tot += m_mapNevent[m_masses[i]];
  }
  bool passed_DAS = true;
  if(NDAS > 0 && Nevent_tot != NDAS) passed_DAS = false;
  // Loop to fill tree
  for(int i = 0; i < Nmass; i++){
    Nevent = m_mapNevent[m_masses[i]];
    NDAS = Nevent; // since we already passed DAS check above, set NDAS to Nevent for filling tree
    MP = m_masses[i].first;
    MC = m_masses[i].second;
    tout->Fill();
  }

  tout->Write("",TObject::kOverwrite);
  delete tout;

  if(!AnalysisBase<Base>::IsData()){
    outfile->mkdir("Histograms");
    outfile->cd("Histograms");
    int Nhisto = histos.size();
    for(int i = 0; i < Nhisto; i++)
      histos[i]->Write("", TObject::kOverwrite);
  }

  TParameter<bool> isData("IsData", AnalysisBase<Base>::IsData()); isData.Write();
  TParameter<bool> isSMS("IsSMS", AnalysisBase<Base>::IsSMS()); isSMS.Write();
  TParameter<bool> isCascades("IsCascades", AnalysisBase<Base>::IsCascades()); isCascades.Write();
  TParameter<bool> isPrivateMC("IsPrivateMC", AnalysisBase<Base>::IsPrivateMC()); isPrivateMC.Write();
  // Example read back from output:
  // auto b = (TParameter<bool>*) f.Get("someBool");
  // bool flag = b->GetVal();
  
  outfile->Close();
  delete outfile;

  m_Trees.clear();
  m_Label2Tree.clear();
  return passed_DAS;
  
}

template class NtupleBase<SUSYNANOBase>;
template class NtupleBase<NANOULBase>;
template class NtupleBase<NANORun3>;

