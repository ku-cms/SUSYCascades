#include <TFile.h>
#include <TError.h>

#include "NtupleBase.hh"
#include "SUSYNANOBase.hh"
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
void NtupleBase<Base>::WriteNtuple(const string& filename, int ichunk, int nchunk, bool do_slim){
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
  cout << NTOT << endl;
  Long64_t N1, N0;
  if(nchunk >= NTOT){
    N1 = ichunk;
    N0 = ichunk-1;
  } else {
    N1 = NTOT/nchunk;
    if(NTOT%nchunk > 0)
      N1++;
    N0 = (ichunk-1)*N1;
    N1 = N0 + N1;
  }

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

      //cout << "Filling " << AnalysisBase<Base>::m_Systematics[s].Label() << " " << do_slim << endl;
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
  int NDAS = 0;
  double Nevent;
  double Nweight = 0.;
  int MP;
  int MC;
  tout->Branch("NDAS", &NDAS);
  tout->Branch("Nevent", &Nevent);
  tout->Branch("Nweight", &Nweight);
  tout->Branch("filetag", &filetag);
  tout->Branch("dataset", &dataset);
  tout->Branch("MP", &MP);
  tout->Branch("MC", &MC);
  // add DAS count
  NeventTool eventTool;
  NDAS = eventTool.EventsInDAS(dataset, filetag);
  int Nmass = m_masses.size();
  for(int i = 0; i < Nmass; i++){
    Nevent = m_mapNevent[m_masses[i]];
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
  
  outfile->Close();
  delete outfile;

  m_Trees.clear();
  m_Label2Tree.clear();

  
}

template class NtupleBase<NANORun3>;

