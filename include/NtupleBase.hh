#ifndef NtupleBase_h
#define NtupleBase_h

#include "AnalysisBase.hh"

template <class Base>
class NtupleBase : public AnalysisBase<Base> {

public:
  NtupleBase(TTree* tree = 0);
  virtual ~NtupleBase();

  bool WriteNtuple(const std::string& filename, int ichunk = 1, int nchunch = 1, bool do_slim = false, int NDAS = 0, const string& DAS_datasetname = "", const string& DAS_filename = "");
  void GetChunks(const Long64_t& NTOT, Long64_t& N0, Long64_t& N1, int ichunk, int nchunk);

protected:
  std::vector<TTree*>     m_Trees;
  std::map<string,vector<TTree*>> m_Label2Tree;
 
private:
  virtual TTree* InitOutputTree(const std::string& sample, bool do_slim = false, bool tree_is_sys = false) = 0;
  virtual void FillOutputTree(TTree* tree, const Systematic& sys = Systematic::Default(), bool do_slim = false) = 0;

  // for event count bookkeeping
  std::vector<std::pair<int,int> > m_masses;
  std::map<std::pair<int,int>,double > m_mapNevent;

};

#endif
