#ifndef SampleTool_h
#define SampleTool_h

#include <iostream>
#include <string>
#include <ctype.h>
#include <map>
#include <vector>
#include <TChain.h>

#include "Process.hh"
#include "ScaleFactorTool.hh"

using std::string;
using std::vector;
using std::pair;

class SampleTool {

public:
  SampleTool(const string& ntuple_path, int year);
  virtual ~SampleTool();

  SampleTool& SetNtuplePath(const string& ntuple_path);
  SampleTool& SetYear(int year);

  double Lumi();
  double HEMLumi(bool HEM_Veto);

  ProcessList Get(const string& name) const;
  ProcessList Get(ProcessType type) const;
  ProcessList GetStrictSignalMatch( const string& name) const;
 
  int NTrees(const Process& proc);
  TChain* Tree(const Process& proc, int itree = -1, std::string treeName=""); // remember to delete TChain!
  string TreeName(const Process& proc, int itree);
  string FileName(const Process& proc, int itree);
  
  bool IsFastSim(const Process& proc, int itree);
  bool FilterDilepton(const Process& proc, int itree);
  SleptonFlavor FilterSleptons(const Process& proc, int itree);
  double GetSampleWeight(const Process& proc, int itree);
  
private:
  string m_Path;
  int    m_iYear;

  int YearMap(int year);

  void InitSMS_treeSys(const string& treeSys, const string& prefix, const string& filename, double weight=1., bool FS=false, bool DL=false, SleptonFlavor kFlavor = kSmuSel);

  void InitSMS(const string& prefix, const string& filename, double weight = 1., bool FS = false, bool DL = false, SleptonFlavor kFlavor = kSmuSel);

  void InitProcMap();
  static bool m_ProcInit;
  static std::map<Process, pair<vector<string>,string> > m_Proc[11];
  static double m_Lumi[11];
  static double m_HEMLumi[2];
  // signal only
  void InitSignalProc(const Process& proc);
  // Note: Need to increase array sizes here to match total number of 'eras'
  static std::map<Process, bool> m_SProcInit[11]; // checked combined normalizations already?
  static std::map<Process, std::map<string,bool> >   m_SProcFS[11]; // FastSim?
  static std::map<Process, std::map<string,bool> >   m_SProcDL[11]; // di-lepton filter (ZToLL or dilepton filter);
  static std::map<Process, std::map<string,SleptonFlavor> >   m_SProcSlepFlavor[11]; // Slepton filter 
  static std::map<Process, std::map<string,double> > m_SProcW[11];  // some additional weight to apply
  
  
};

#endif
