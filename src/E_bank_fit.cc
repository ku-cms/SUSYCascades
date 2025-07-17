#include "E_bank_fit.hh"

E_bank_fit::E_bank_fit(int year, std::string f16, std::string f17, std::string f18, std::string histPath){
  //set default year to draw from
  _year = year;
  TFile* f1 = TFile::Open(f16.c_str());
  f1->cd(histPath.c_str());
  std::string name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  populateMap((TCanvas*) f1->Get((histPath+name).c_str()),_map16);
  f1->Close();

  TFile* f2 = TFile::Open(f17.c_str());
  f2->cd(histPath.c_str());
  name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  populateMap((TCanvas*) f2->Get((histPath+name).c_str()),_map17);
  f2->Close();

  TFile* f3 = TFile::Open(f18.c_str());
  f3->cd(histPath.c_str());
  name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  populateMap((TCanvas*) f3->Get((histPath+name).c_str()),_map18);
  f3->Close();

}
void E_bank_fit::extractFit( TCanvas* hCanv, double fitLow, double fitUp, std::map<std::pair<double,double>, std::pair<double,double> >& _fmap, int year){
	
  if(year==2016) _syst_low= _syst_low16;
  if(year==2017) _syst_low= _syst_low17;
  if(year==2018) _syst_low= _syst_low18;
  //if not specified default 2017
  if(year==0) _syst_low= _syst_low17;
	
  TList* t = hCanv->GetListOfPrimitives();
  //t->Print();
  TH2D* h = (TH2D*) hCanv->GetPrimitive(t->At(1)->GetName());
  int nbinsX = h->GetNbinsX();	

  //declare variables to be used inside
  double diff,sig;
  double prob1,prob2;
  int hbinX;
  double sum2=0.;
  //loop and get projection
  double pt_edge,sys;
  for(int i=1; i<=nbinsX; i++){
    TH1D* proj = h->ProjectionY(("hprojY"+std::to_string(i)).c_str(),i,i);
    //loop through and adjust errors according to systematics
    for(int j=1; j<=proj->GetNbinsX(); j++){
      pt_edge = proj->GetXaxis()->GetBinLowEdge(j);
      if(pt_edge < 20) sys = _syst_low[i-1];
      if(pt_edge >= 20) sys = _syst[i-1];
      proj->SetBinError(j,  std::sqrt( proj->GetBinError(j)*proj->GetBinError(j) + sys*sys ) );
    }


    proj->Fit("pol2","Q","",fitLow,fitUp);
    //proj->Fit("pol1","+","",fitLow,fitUp);
    TF1* fit2 = proj->GetFunction("pol2");
    //TF1* fit1 = proj->GetFunction("pol1");

		
    double X[1] = {0.};		
    //choose better fit pol1 or pol2
    prob2 = fit2->GetProb();
    //prob1 = fit1->GetProb();
    TF1* fit=fit2;
    //std::cout<<"pol1 P="<<prob1<<" pol2 P="<<prob2<<" Selecting ";
    //std::cout<<"pol2 P="<<prob2<<std::endl;		
    /*if( prob1 < 0.05){
      fit = fit2;
      std::cout<<" pol2"<<std::endl;
      }
      else{
      fit = fit1;
      std::cout<<" pol1"<<std::endl;
      }*/
    //if( prob2 < 0.01){
    //  std::cout<<"Warning >> pol2 with P<1%"<<std::endl;
    //}
    //TH1D *hint = new TH1D(("hint"+std::to_string(i)).c_str(),"Fitted CI", 17, 3, 20);
    TH1D *hint = new TH1D(("hint"+std::to_string(i)).c_str(),"Fitted CI", 17, 3, 20);//simple binning for gsb plots (this is temporary and not intended for efficiency extraction
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint, 0.68);
    hint->SetFillColor(kRed);
    /*hbinX = proj->GetNbinsX();
      sum2=0.;
      for(int j =1; j< hbinX; j++){
      X[0] = proj->GetBinCenter(j);
      diff = proj->GetBinContent(j) - fit->EvalPar(X);
      sum2+= diff*diff;
      }		
      sig = std::sqrt( sum2/fit->GetNDF() );
		
      std::cout<<"Model systematic error w/ chi2/ndf=1 : "<<sig<<std::endl;
    */
    //debugging histograms///////////
    TCanvas* cnew = new TCanvas();
    proj->Draw();		
    hint->Draw("E3 SAME");
    proj->Draw("SAME");
    //cnew->SaveAs(("test1_"+std::to_string(i)+".pdf").c_str());
    TH1D* clone = (TH1D*) hint->Clone();
    //clone contains the CI with added systematics
    for(int j=1; j<= hint->GetNbinsX(); j++){
      pt_edge = clone->GetXaxis()->GetBinLowEdge(j);
      if(pt_edge < 20) sys = _syst_low[i-1];
      if(pt_edge >= 20) sys = _syst[i-1];
      clone->SetBinError(j,  std::sqrt(clone->GetBinError(j)*clone->GetBinError(j)+_syst_low[i]*_syst_low[i])  );
    }
    TCanvas* canclone = new TCanvas();
    proj->Draw();		
    clone->Draw("E3 SAME");
    proj->Draw("SAME");
    //canclone->SaveAs(("test2_"+std::to_string(i)+".pdf").c_str());
    ////////////////////////////////

    //discretize function and errors into a map (addressed fmap)
    //hardcode all the way down to pt1, these will be unused but its generally safe
    /*double ptedge=1.0; //split threshold better be set here
      sig=0.;
      while(ptedge< _split_threshold){
      X[0] = ptedge + 0.5 ;
      _fmap[ std::make_pair(ptedge,h->GetXaxis()->GetBinLowEdge(i)) ] = std::make_pair( fit->EvalPar(X) , sig );
      ptedge = ptedge + 1.0;
      }*/

    for(int j=1; j<= hint->GetNbinsX(); j++){
      X[0] = hint->GetXaxis()->GetBinCenter(j) ;
      _fmap[ std::make_pair( hint->GetXaxis()->GetBinLowEdge(j), h->GetXaxis()->GetBinLowEdge(i)) ] = std::make_pair( fit->EvalPar( X ), clone->GetBinError(j) );
    }
		

  }
	

}
void E_bank_fit::doLowPtFit(double fitLow, double fitUp, double split_threshold, std::string ffit16, std::string ffit17, std::string ffit18, std::string histPath ){
  _split_threshold = split_threshold;
  //open files similarly to e_bank
  TFile* f1 = TFile::Open(ffit16.c_str());
  f1->cd(histPath.c_str());
  std::string name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  //std::cout<< "Performing 2016 fits ... ";
  extractFit( (TCanvas*) f1->Get((histPath+name).c_str()),fitLow ,fitUp , _fitmap16, 2016);
  f1->Close();

  TFile* f2 = TFile::Open(ffit17.c_str());
  f2->cd(histPath.c_str());
  name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  //std::cout<< "Performing 2017 fits ... ";
  extractFit( (TCanvas*) f2->Get((histPath+name).c_str()),fitLow ,fitUp , _fitmap17, 2017);
  f2->Close();

  TFile* f3 = TFile::Open(ffit18.c_str());
  f3->cd(histPath.c_str());
  name = (gDirectory->GetListOfKeys()->At(0)->GetName());
  //std::cout<< "Performing 2018 fits ... ";
  extractFit( (TCanvas*) f3->Get((histPath+name).c_str()),fitLow ,fitUp , _fitmap18, 2017);
  f3->Close();	
}
std::pair<double,double> E_bank_fit::getPair(double pt, double eta, int year){
  if(year>0)_year=year;
  //check if fits are initialized
  bool fitMapIsEmpty= false;
  if(_fitmap16.size()==0 || _fitmap17.size()==0 || _fitmap18.size()==0) fitMapIsEmpty =true;
  if(pt>= _split_threshold || fitMapIsEmpty){
    if(_year == 2016) return _map16[transformValue(pt,eta, _map16)];
    if(_year == 2017) return _map17[transformValue(pt,eta, _map17)];
    if(_year == 2018) return _map18[transformValue(pt,eta, _map18)];
  }
  else if(pt< _split_threshold){
    if(_year == 2016) return _fitmap16[transformValue(pt,eta, _fitmap16)];
    if(_year == 2017) return _fitmap17[transformValue(pt,eta, _fitmap17)];
    if(_year == 2018) return _fitmap18[transformValue(pt,eta, _fitmap18)];
  }
  return std::make_pair(-1,-1);
}
double E_bank_fit::getEfficiency(double pt, double eta, int year){
  if(year>0)_year=year;
  //check if fits are initialized
  bool fitMapIsEmpty= false;
  if(_fitmap16.size()==0 || _fitmap17.size()==0 || _fitmap18.size()==0) fitMapIsEmpty =true;
  if(pt>= _split_threshold || fitMapIsEmpty){
    if(_year == 2016) return _map16[transformValue(pt,eta, _map16)].first;
    if(_year == 2017) return _map17[transformValue(pt,eta, _map17)].first;
    if(_year == 2018) return _map18[transformValue(pt,eta, _map18)].first;
  }
  else if(pt< _split_threshold){
    if(_year == 2016) return _fitmap16[transformValue(pt,eta, _fitmap16)].first;
    if(_year == 2017) return _fitmap17[transformValue(pt,eta, _fitmap17)].first;
    if(_year == 2018) return _fitmap18[transformValue(pt,eta, _fitmap18)].first;
  }
  return -1;
}
double E_bank_fit::getError(double pt, double eta, int year){
  if(year>0)_year=year;
  //check if fits are initialized
  bool fitMapIsEmpty= false;
  if(_fitmap16.size()==0 || _fitmap17.size()==0 || _fitmap18.size()==0) fitMapIsEmpty =true;
  if(pt>= _split_threshold || fitMapIsEmpty){
    if(_year == 2016) return _map16[transformValue(pt,eta, _map16)].second;
    if(_year == 2017) return _map17[transformValue(pt,eta, _map17)].second;
    if(_year == 2018) return _map18[transformValue(pt,eta, _map18)].second;
  }
  else if(pt< _split_threshold){
    if(_year == 2016) return _fitmap16[transformValue(pt,eta, _fitmap16)].second;
    if(_year == 2017) return _fitmap17[transformValue(pt,eta, _fitmap17)].second;
    if(_year == 2018) return _fitmap18[transformValue(pt,eta, _fitmap18)].second;
  }
  return -1;
}
void E_bank_fit::setSystematicsLow(std::vector<double> errors, int year){
  if(year == 0) _syst_low = errors;
  if(year == 2016) _syst_low16 = errors;
  if(year == 2017) _syst_low17 = errors;
  if(year == 2018) _syst_low18 = errors;
	
}
/*void E_bank_fit::applySystematicsLow(std::vector<double> errors, int year=0){
  _syst_low = errors;
  if(year==0){
  addError(errors,_fitmap16);
  addError(errors,_fitmap17);
  addError(errors,_fitmap18);
  }
  if(year==2016) addError(errors,_fitmap16);
  if(year==2017) addError(errors,_fitmap17);
  if(year==2018) addError(errors,_fitmap18);
  }
*/
