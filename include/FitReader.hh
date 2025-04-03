#ifndef FitReader_HH
#define FitReader_HH

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TColor.h>
#include <iostream>
#include <vector>

#include "Category.hh"
#include "Process.hh"

using std::map;
using std::string;
using std::pair;
using std::basic_string;
///////////////////////////////////////////
////////// FitReader class
///////////////////////////////////////////

class FitReader {
public:
  FitReader(const string& inputfile,
	    const string& otherfile = "",
	    const string& otherfold = "");

  virtual ~FitReader();

  void PrintCategories(bool verbose = false);
  void PrintProcesses(bool verbose = false);

  // checks if filled, loads histogram for use
  bool IsFilled(const Category&   cat,
		const Process&    proc,
		const Systematic& sys = Systematic::Default()) const;

  // checks if filled, deletes histogram after check
  bool IsThere(const Category&   cat,
	       const Process&    proc,
	       const Systematic& sys = Systematic::Default()) const;

  const TH1D* GetHistogram(const Category&   cat,
			   const Process&    proc,
			   const Systematic& sys = Systematic::Default()) const;

  // returns a 1-bin histogram added over all RISR/Mperp bins, all cats, all procs
  TH1D* GetIntegralHist(const string&       name,
			const CategoryList& cats,
			const ProcessList&  procs,
			const Systematic&   sys = Systematic::Default()) const;

  // returns a RISR/Mperp histogram added over all cats, all procs (need to have same binning or errors)
  TH1D* GetAddedHist(const string&       name,
		     const CategoryList& cats,
		     const ProcessList&  procs,
		     const Systematic&   sys = Systematic::Default()) const;

  vector<double> GetAddedHistValues(const CategoryList& cats,
				    const ProcessList&  procs,
				    const Systematic&   sys = Systematic::Default()) const;

  vector<double> IntegrateMperp(const FitBin fitBin, const vector<double> histVec) const;
  TH1D* IntegrateMperp(const TString& name, const FitBin fitBin, TH1D* hist) const;
  double CalculateZbi(double Nsig, double Nbkg, double deltaNbkg = 0.2);

  bool IsFilled2D(const Category&   cat,
		  const Process&    proc,
		  const Systematic& sys = Systematic::Default()) const;

  const TH2D* GetHistogram2D(const Category&   cat,
			     const Process&    proc,
			     const Systematic& sys = Systematic::Default()) const;

  double Integral(const Category&   cat,
		  const Process&    proc,
		  const Systematic& sys = Systematic::Default()) const;

  bool HasSystematic(const Process& proc, const Systematic& sys) const;

  string GetSignalTitle(const string&);

  VS GetChannels() const;  
  const ProcessList&  GetProcesses() const;
  const CategoryList& GetCategories(const string& channel = "") const;         
  const Systematics&  GetSystematics() const;

  map<string,VS> m_Strings;
  mutable map<Process,Systematics> m_ProcSys;

  TColor cms1, cms2, cms3, cms4, cms5, cms6, cms7, cms8, cms9, cms10;
  
  /*
  const TColor cms1(8001, 0.24705882352941178, 0.5647058823529412, 0.8549019607843137);
  const TColor cms2(8002, 1., 0.6627450980392157, 0.054901960784313725);
  const TColor cms3(8004, 0.7411764705882353, 0.12156862745098039, 0.00392156862745098);
  const	TColor cms4(8005, 0.5803921568627451, 0.6431372549019608, 0.6352941176470588);
  const	TColor cms5(8006, 0.5137254901960784, 0.17647058823529413, 0.7137254901960784);
  const	TColor cms6(8007, 0.6627450980392157, 0.4196078431372549, 0.34901960784313724);
  const	TColor cms7(8008, 0.9058823529411765, 0.38823529411764707, 0.);
  const	TColor cms8(8009, 0.7254901960784313, 0.6745098039215687, 0.4392156862745098);
  const	TColor cms9(8010, 0.44313725490196076, 0.4588235294117647, 0.5058823529411764);
  const TColor cms10(8011, 0.5725490196078431, 0.8549019607843137, 0.8666666666666667);
  */
protected:
  mutable TFile  m_File;
  mutable TFile* m_FilePtr;
  string         m_FileFold;

  mutable map<Process,map<Category,TH1D*> > m_ProcHist;
  mutable map<Process,map<Systematic,map<Category,pair<TH1D*,TH1D*> > > >m_ProcHistSys;
  mutable map<Process,map<Category,TH2D*> > m_ProcHist_2D;
  mutable map<Process,map<Systematic,map<Category,pair<TH2D*,TH2D*> > > >m_ProcHistSys_2D;

  ProcessList                      m_Proc;
  mutable map<string,CategoryList> m_Chan;
  mutable map<string,bool>         m_CatLabel;
  CategoryList                     m_Cat;
  Systematics                      m_Sys;
  
  ProcessBranch m_ProcBranch;
  void ReadProcesses();
  
  CategoryBranch m_CatBranch;
  void ReadCategories();
  
};

#endif
