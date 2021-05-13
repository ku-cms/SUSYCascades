#ifndef shapeTemplate_HH
#define shapeTemplate_HH


#include <vector>
#include <iostream>

#include "Category.hh"
#include <TH1D.h>
#include "Criteria.hh"
#include "CategoryTree.hh"
#include "Process.hh"



using std::map;
using std::string;
using std::pair;
using std::basic_string;
using std::vector;
class shapeTemplate{

	public:
		shapeTemplate(TH1D* hist,TH1D* hist_cons);
		virtual ~shapeTemplate();
		void sortBins(vector<pair<int,double>>& array);
		TH1D* replaceHistogram();
		double compareShapes();
		void sortHistograms();
			

	private:
		int m_idx = 4;
		int m_nBins;
		double m_norm;
		TH1D* m_hist_OG;
		TH1D* m_hist_cons;
		TH1D* m_hist_rec;
		TH1D* m_hist_scaled;
		double m_minVal = 10.;
		vector<pair<int,double>> m_binsSorted;
		vector<pair<int,double>> m_binsSorted_cons;
		
		void normalizeHistogram(TH1D* hist);
		void unweightHistogram(TH1D* hist);
		
		void replaceBins(int idx=0, double Prem=1., double Prem_cons=1.);
		void setErrors();


};


class shapeTemplateTool{
public:
	shapeTemplateTool(const CategoryTree& CT, ProcessList proc, const string &inputfile);
	virtual ~shapeTemplateTool();

	void createTemplates(bool smooth=false);
	CategoryTree getCategoryTree();
	VS getProcess();


private:
	string m_file;
	VS m_proc;
	VS ST_fakes;
	VS dom_fakes;
	CategoryTree m_CT;
	vector<std::pair<TH1D*,string>> m_histsAndLabels;

	map<string,VS> m_domToRare;

	map<string,double> m_nameToNorm;
	map<int,TH1D*> m_listToHist;
	map<string,const char*> m_nameToTitle;
};

#endif

