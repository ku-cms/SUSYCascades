#include "LHETool.hh"

void LHETool::GetFromLHAPDF(int& nPDFWeights, std::vector<double>& PDFWeights, int id1, int id2, float x1, float x2, float scalePDF, int year){

#ifdef _CMSSW_
  std::string errorTypeCount = pdfset.errorType();
  size_t npar = std::count(errorTypeCount.begin(),errorTypeCount.end(),'+');
  nPDFWeights = npar;
  PDFWeights.resize(pdfset.size());
  std::fill(PDFWeights.begin(),PDFWeights.end(),0.0);
  double x1_cen = pdfs[0]->xfxQ(id1,x1,scalePDF);
  double x2_cen = pdfs[0]->xfxQ(id2,x2,scalePDF);
  double normweight = x1_cen*x2_cen;
  for(int i = 0; i < pdfset.size(); i++){
    double x1_imem = pdfs[i]->xfxQ(id1,x1,scalePDF);
    double x2_imem = pdfs[i]->xfxQ(id2,x2,scalePDF);
    PDFWeights[i] = x1_imem*x2_imem/normweight;
  }
#endif
 
}

void LHETool::AddLHAPDF(int year)
{
#ifdef _CMSSW_
  if(year == 2016)
    pdfset = LHAPDF::getPDFSet("NNPDF30_lo_as_0130"); 
  else // NNPDF31_nnlo_as_0118 for preUL '17 & '18, need to check for UL & Run3
    pdfset = LHAPDF::getPDFSet("NNPDF31_nnlo_as_0118"); 
  pdfs = pdfset.mkPDFs(); 
#endif
}

double LHETool::GetWeight(unsigned int nLHEPdfWeight, float* LHEPdfWeight, int Generator_id1, int Generator_id2, float Generator_x1, float Generator_x2, float Generator_scalePDF, int year, int updown){
  if(updown == 0) return 1.;
  double w84 = 1.;
  double w16 = 1.;
  double mean = 1.; 
  double err = 0.;
  int nPDFWeights = nLHEPdfWeight;
  std::vector<double> PDFWeights;

  if(nPDFWeights == 0){
#ifdef _CMSSW_
    GetFromLHAPDF(nPDFWeights,PDFWeights,Generator_id1,Generator_id2,Generator_x1,Generator_x2,Generator_scalePDF,year);
#else
    return 1;
#endif
  } else {
    for(int i = 0; i < nPDFWeights; i++) { PDFWeights.push_back(LHEPdfWeight[i]); }
  }
  
  if(nPDFWeights == 101)
  {
   std::sort(PDFWeights.begin()+1,PDFWeights.end());
   w84 = PDFWeights[84];
   w16 = PDFWeights[16];
   mean = (w84+w16)/2.;
   err = (w84-w16)/2.;
  }
  else
  {
   mean = std::accumulate(PDFWeights.begin()+1,PDFWeights.end(),0.0)/(PDFWeights.size()-1.);
   double sq_sum = std::inner_product(PDFWeights.begin()+1,PDFWeights.end(),PDFWeights.begin()+1,0.0);
   err = sqrt(sq_sum/(PDFWeights.size()-1.) - mean*mean);
  }
  if(updown > 0)
    return (1.+err/mean);
  else if(updown < 0)
    return (1.-err/mean);
  else
    return 1.;
}


