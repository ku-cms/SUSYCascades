#include "../include/CategoryTool.hh"

#include "Leptonic.hh"
#include "Hadronic.hh"

using std::cout;
using std::endl;

///////////////////////////////////////////
////////// CategoryTool class
///////////////////////////////////////////

CategoryTool::CategoryTool() {}

CategoryTool::~CategoryTool() {}

// Get all 0L,1L,2L,3L categories
CategoryList CategoryTool::GetCategories(bool maskSR) const {
  CategoryList Categories;
  Categories += GetCategories_0L(maskSR);
  Categories += GetCategories_1L(maskSR);
  Categories += GetCategories_2L(maskSR);
  Categories += GetCategories_3L(maskSR);

  return Categories;
}

FitBin CategoryTool::GetFitBin(int Nlep, int index, bool maskSR) const {
  vector<RBin> RBins;

  if(Nlep == 0){
    if(index == 0){
      RBins.push_back(RBin(0.95,  0.985, VD().a(0.).a(40.)));
      if(!maskSR)
	RBins.push_back(RBin(0.985, 1.00,  VD().a(0.).a(5.).a(10.).a(40.)));
    }
    if(index == 1){ 
      RBins.push_back(RBin(0.8,   0.90,  VD().a(0.).a(100.)));
      RBins.push_back(RBin(0.90,  0.93,  VD().a(0.).a(100.)));//remove Mperp B101
      if(!maskSR){
	RBins.push_back(RBin(0.93,  0.96,  VD().a(0.).a(20.).a(100.)));
	RBins.push_back(RBin(0.96,  1.0,   VD().a(0.).a(15.).a(100.)));
      }
    }
    if(index == 2){
      RBins.push_back(RBin(0.65,  0.75,  VD().a(0.).a(120.)));//remove Mperp B101
      RBins.push_back(RBin(0.75,  0.85,  VD().a(0.).a(120.)));
      if(!maskSR){
	RBins.push_back(RBin(0.85,  0.9,   VD().a(0.).a(120.)));//remove Mperp B101
	RBins.push_back(RBin(0.9,   0.95,  VD().a(0.).a(20.).a(120.)));
	RBins.push_back(RBin(0.95,  1.0,   VD().a(0.).a(120.)));
      }
    }
    if(index == 3){
      RBins.push_back(RBin(0.55,  0.65, VD().a(0.).a(200.)));//remove Mperp B101
      RBins.push_back(RBin(0.65,  0.75, VD().a(0.).a(200.)));
      if(!maskSR){
	RBins.push_back(RBin(0.75,  0.85, VD().a(0.).a(200.)));//remove Mperp B101
	RBins.push_back(RBin(0.85,  0.9,  VD().a(0.).a(50.).a(200.)));
	RBins.push_back(RBin(0.9,   1.0,  VD().a(0.).a(200.)));
      }
    }
    if(index == 4){
      RBins.push_back(RBin(0.55,  0.65, VD().a(0.).a(300.)));//remove Mperp B101
      RBins.push_back(RBin(0.65,  0.75, VD().a(0.).a(300.)));
      if(!maskSR){
	RBins.push_back(RBin(0.75,  0.85, VD().a(0.).a(80.).a(300.)));
	RBins.push_back(RBin(0.85,  1.0,  VD().a(0.).a(300.)));
      }
    }
    if(index == 5){
      RBins.push_back(RBin(0.5,   0.6,  VD().a(0.).a(400.)));//remove Mperp B101
      RBins.push_back(RBin(0.6,   0.7,  VD().a(0.).a(400.)));
      if(!maskSR){
	RBins.push_back(RBin(0.7,   0.8,  VD().a(0.).a(150.).a(400.)));
	RBins.push_back(RBin(0.8,   1.0,  VD().a(0.).a(400.)));
      }
    }
  }
  
  if(Nlep == 1){
    if(index == 0){
     // RBins.push_back(RBin(0.9,   0.96,   VD().a(0.).a(40.))); //remove RISR B101
      RBins.push_back(RBin(0.96,  0.98,   VD().a(0.).a(10.).a(40.)));
      if(!maskSR){
	RBins.push_back(RBin(0.98,  1.00,   VD().a(0.).a(5.).a(10.).a(40.)));
      }
    }
    if(index == 1){
    //  RBins.push_back(RBin(0.65,  0.75,  VD().a(0.).a(50.).a(120.))); //remove RISR B101
    //  RBins.push_back(RBin(0.75,  0.85,  VD().a(0.).a(40.).a(120.))); //keeping masked region the same
      if(!maskSR){
	RBins.push_back(RBin(0.85,  0.9,   VD().a(0.).a(30.).a(120.)));
	RBins.push_back(RBin(0.9,   0.95,  VD().a(0.).a(20.).a(120.)));
	RBins.push_back(RBin(0.95,  1.0,   VD().a(0.).a(120.)));
      }
    }
    if(index == 2){
    //  RBins.push_back(RBin(0.55,  0.65, VD().a(0.).a(110.).a(200.))); //remove RISR B101
   //   RBins.push_back(RBin(0.65,  0.75, VD().a(0.).a(90.).a(200.)));
      if(!maskSR){
	RBins.push_back(RBin(0.8,  0.85, VD().a(0.).a(70.).a(200.)));
	RBins.push_back(RBin(0.85,  0.9,  VD().a(0.).a(50.).a(200.)));
	RBins.push_back(RBin(0.9,   1.0,  VD().a(0.).a(200.)));
      }
    }
    if(index == 3){
  //   RBins.push_back(RBin(0.55,  0.65, VD().a(0.).a(150.).a(300.))); //remove RISR  B101
      RBins.push_back(RBin(0.65,  0.75, VD().a(0.).a(100.).a(300.)));
      if(!maskSR){
	RBins.push_back(RBin(0.75,  0.85, VD().a(0.).a(80.).a(300.)));
	RBins.push_back(RBin(0.85,  1.0,  VD().a(0.).a(300.)));
      }
    }
    if(index == 4){
  //    RBins.push_back(RBin(0.5,   0.6,  VD().a(0.).a(210.).a(400.))); //remove RISR B101
      RBins.push_back(RBin(0.6,   0.7,  VD().a(0.).a(180.).a(400.)));
      if(!maskSR){
	RBins.push_back(RBin(0.7,   0.8,  VD().a(0.).a(150.).a(400.)));
	RBins.push_back(RBin(0.8,   1.0,  VD().a(0.).a(400.)));
      }
    }
  }

  if(Nlep == 2){//index 0,1,2,3 refer to specifically gold configurations for maskSR
    if(index == 0){//2L 0J slvr CR will match gold, bronz will include 2 bins
       RBins.push_back(RBin(0.6,   0.7,   VD().a(0.).a(50.).a(120)));//adding Extra bins for CR uniformity
       RBins.push_back(RBin(0.7,   0.8,   VD().a(0.).a(40.).a(120.)));
       if(!maskSR){
         RBins.push_back(RBin(0.8,   0.9,   VD().a(0.).a(30.).a(120.)));
	 RBins.push_back(RBin(0.9,   0.95,  VD().a(0.).a(20.).a(120.)));
	 RBins.push_back(RBin(0.95,  1.0,   VD().a(0.).a(15.).a(120.)));
       }
    }
    if(index == 1){
      RBins.push_back(RBin(0.5,   0.6,   VD().a(0.).a(100).a(200.)));//added for uniformity
      RBins.push_back(RBin(0.6,   0.7,  VD().a(0.).a(80.).a(200.)));//changed from 0.55->0.6 for uniformity
      if(!maskSR){
        RBins.push_back(RBin(0.7,   0.8,  VD().a(0.).a(60.).a(200.)));
	RBins.push_back(RBin(0.8,   0.9,  VD().a(0.).a(40.).a(200.)));
	RBins.push_back(RBin(0.9,   1.0,  VD().a(0.).a(30.).a(200.)));
      }
    }
    if(index == 2){
      RBins.push_back(RBin(0.5,   0.65,  VD().a(0.).a(100.).a(300.)));
      RBins.push_back(RBin(0.65, 0.75, VD().a(0.).a(80.).a(300.)));//shifted from blind region into CR, it'll be okay (highest Signal presence though)
      if(!maskSR){
       // RBins.push_back(RBin(0.65,  0.75,  VD().a(0.).a(80.).a(300.)));
	RBins.push_back(RBin(0.75,  0.85,  VD().a(0.).a(60.).a(300.)));
	RBins.push_back(RBin(0.85,  1.0,   VD().a(0.).a(300.)));
      }
    }
    if(index == 3){
      RBins.push_back(RBin(0.5,   0.65,  VD().a(0.).a(130.).a(400.)));
      if(!maskSR){
	RBins.push_back(RBin(0.65,  0.8,   VD().a(0.).a(100.).a(400.)));
	RBins.push_back(RBin(0.8,   1.0,   VD().a(0.).a(400.)));
      }
    }
    if(index == 10){//2L 0J bron
	RBins.push_back(RBin(0.7,   0.8,   VD().a(0.).a(40.).a(120.)));
        RBins.push_back(RBin(0.8,   0.9,   VD().a(0.).a(30.).a(120.)));
       if(!maskSR){
         RBins.push_back(RBin(0.9,   0.95,  VD().a(0.).a(20.).a(120.)));
         RBins.push_back(RBin(0.95,  1.0,   VD().a(0.).a(15.).a(120.)));
       }
		
    }
    if(index == 11){//2L 1J slvr/bron
	RBins.push_back(RBin(0.55,  0.7,  VD().a(0.).a(80.).a(200.)));
        RBins.push_back(RBin(0.7,   0.8,  VD().a(0.).a(60.).a(200.)));
      if(!maskSR){
	RBins.push_back(RBin(0.8,   0.9,  VD().a(0.).a(40.).a(200.)));
        RBins.push_back(RBin(0.9,   1.0,  VD().a(0.).a(30.).a(200.)));
      }

    }
    if(index == 12){//2L ge2J slvr/bron
	RBins.push_back(RBin(0.5,   0.65,  VD().a(0.).a(100.).a(300.)));
        RBins.push_back(RBin(0.65,  0.75,  VD().a(0.).a(80.).a(300.)));
      if(!maskSR){
        RBins.push_back(RBin(0.75,  0.85,  VD().a(0.).a(60.).a(300.)));
        RBins.push_back(RBin(0.85,  1.0,   VD().a(0.).a(300.)));
      }

    }

  }

  if(Nlep == 3){
    if(index == 0){
       RBins.push_back(RBin(0.6,  0.7,  VD().a(0.).a(250.)));
       RBins.push_back(RBin(0.7,  0.8,  VD().a(0.).a(250.)));
       if(!maskSR){
	 RBins.push_back(RBin(0.8,  0.9,  VD().a(0.).a(250.)));
	 RBins.push_back(RBin(0.9,  1.0,  VD().a(0.).a(250.)));
       }
    }
    if(index == 1){
      RBins.push_back(RBin(0.55,   0.7,  VD().a(0.).a(250.)));
      if(!maskSR){
	RBins.push_back(RBin(0.7,    0.85, VD().a(0.).a(250.)));
	RBins.push_back(RBin(0.85,   1.0,  VD().a(0.).a(250.)));
      }
    }
  }


  FitBin Bins(RBins);
  return Bins;

}

Hadronic CategoryTool::GetHadronicRegion(int Njet, int index) const {

  /////////
  // 0 jet
  /////////
  if(Njet == 0){
    // fully inclusive
    if(index == 0){
      Hadronic H(0, 0, 0, "incl");
      H += Hadronic(0, 0, -1);
      H += Hadronic(-1, 0, 0);
      H += Hadronic(-1, 0, -1);
      H += Hadronic(-1, -1, 0);
      H += Hadronic(-1, -1, -1);
      return H;
    }

    // C1
    if(index == 1){
      Hadronic H(0, 0, 0, "0j0sv");
      return H;
    }

    // C2
    if(index == 2){
      Hadronic H(0, 0, -1, "0jge1sv");
      return H;
    }

    // exactly 1 sv
    if(index == 3){
      Hadronic H(0, 0, 1, "0j1sv");
      return H;
    }

    // ge 2 sv
    if(index == 4){
      Hadronic H(0, 0, -2, "0jge2sv");
      return H;
    }

    // 0j
    if(index == 5){
      Hadronic H(0, 0, 0, "0j");
      H += Hadronic(0, 0, -1);
      return H;
    }
  }

  /////////
  // 1 jet
  /////////
  if(Njet == 1){
    // fully inclusive ge1 jet
    if(index == 0){
      Hadronic H(-1, 0, 0, "ge1j");
      H += Hadronic(-1, 0, -1);
      H += Hadronic(-1, -1, 0);
      H += Hadronic(-1, -1, -1);
      return H;
    }

    // inclusive exactly 1 jet
    if(index == 1){
      Hadronic H(1, 0, 0, "1j");
      H += Hadronic(1, 0, -1);
      H += Hadronic(1, -1, 0);
      H += Hadronic(1, -1, -1);
      return H;
    }

    // C3 exactly 1 jet 0b 0sv
    if(index == 2){
      Hadronic H(1, 0, 0, "1j0b0sv");
      return H;
    }

    // C4 exactly 1 jet 1b 0sv
    if(index == 3){
      Hadronic H(1, 1, 0, "1j1b0sv");
      return H;
    }

    // C5 exactly 1 jet 0b ge1sv
    if(index == 4){
      Hadronic H(1, 0, -1, "1j0bge1sv");
      return H;
    }

    // C6 exactly 1 jet 1b ge1sv
    if(index == 5){
      Hadronic H(1, 1, -1, "1j1bge1sv");
      return H;
    }

    // B1 ge1 jet 0 b
    if(index == 6){
      Hadronic H(-1, 0, 0, "ge1j0b");
      H += Hadronic(-1, 0, -1);
      return H;
    }

    // B2 ge1 jet ge1 b
    if(index == 7){
      Hadronic H(-1, -1, 0, "ge1jge1b");
      H += Hadronic(-1, -1, -1);
      return H;
    }

    // exactly 1 jet 0b
    if(index == 8){
      Hadronic H(1, 0, 0, "1j0b");
      H += Hadronic(1, 0, -1);
      return H;
    }

    // exactly 1 jet 1b
    if(index == 9){
      Hadronic H(1, 1, 0, "1j1b");
      H += Hadronic(1, 1, -1);
      return H;
    }

    // not C5+C6, exactly 1 jet incl b 0sv
    if(index == 10){
      Hadronic H(1, 0, 0, "1j0sv");
      H += Hadronic(1, 1, 0);
      return H;
    }

    // C5+C6 exactly 1 jet incl b ge1sv
    if(index == 11){
      Hadronic H(1, 0, -1, "1jge1sv");
      H += Hadronic(1, 1, -1);
      return H;
    }
  }

  /////////
  // 2 jet
  /////////
  if(Njet == 2){
    // fully inclusive ge2 jet
    if(index == 0){
      Hadronic H(-2, 0, 0, "ge2j");
      H += Hadronic(-2, 0, -1);
      H += Hadronic(-2, -1, 0);
      H += Hadronic(-2, -1, -1);
      return H;
    }

    // fully inclusive exactly 2 jet
    if(index == 1){
      Hadronic H(2, 0, 0, "2j");
      H += Hadronic(2, 0, -1);
      H += Hadronic(2, -1, 0);
      H += Hadronic(2, -1, -1);
      return H;
    }

    // C7 2 jet 0b
    if(index == 2){
      Hadronic H(2, 0, 0, "2j0b");
      H += Hadronic(2, 0, -1);
      return H;
    }

    // C8 2 jet ge1b
    if(index == 3){
      Hadronic H(2, -1, 0, "2jge1b");
      H += Hadronic(2, -1, -1);
      return H;
    }
    
    // C8 2 jet 1b
    if(index == 4){
      Hadronic H(2, 1, 0, "2j1b");
      H += Hadronic(2, 1, -1);
      return H;
    }

    // C9 2 jet 2b
    if(index == 5){
      Hadronic H(2, 2, 0, "2j2b");
      H += Hadronic(2, 2, -1);
      return H;
    }

    // ge2j 0b
    if(index == 6){
      Hadronic H(-2, 0, 0, "ge2j0b");
      H += Hadronic(-2, 0, -1);
      return H;
    }

    // ge2j ge1b
    if(index == 7){
      Hadronic H(-2, -1, 0, "ge2jge1b");
      H += Hadronic(-2, -1, -1);
      return H;
    }

    //  ge2j 1b
    if(index == 8){
      Hadronic H(-2, 1, 0, "ge2j1b");
      H += Hadronic(-2, 1, -1);
      return H;
    }

    //  ge2j ge2b
    if(index == 9){
      Hadronic H(-2, -2, 0, "ge2jge2b");
      H += Hadronic(-2, -2, -1);
      return H;
    }
  }

  /////////
  // 3 jet
  /////////
  if(Njet == 3){
    // fully inclusive ge3 jet
    if(index == 0){
      Hadronic H(-3, 0, 0, "ge3j");
      H += Hadronic(-3, 0, -1);
      H += Hadronic(-3, -1, 0);
      H += Hadronic(-3, -1, -1);
      return H;
    }
    
    // fully inclusive exactly 3 jet
    if(index == 1){
      Hadronic H(3, 0, 0, "3j");
      H += Hadronic(3, 0, -1);
      H += Hadronic(3, -1, 0);
      H += Hadronic(3, -1, -1);
      return H;
    }

    // C10 3j 0b
    if(index == 2){
      Hadronic H(3, 0, 0, "3j0b");
      H += Hadronic(3, 0, -1);
      return H;
    }

    // C11 3j ge1b
    if(index == 3){
      Hadronic H(3, -1, 0, "3jge1b");
      H += Hadronic(3, -1, -1);
      return H;
    }

    //  3j 1b
    if(index == 4){
      Hadronic H(3, 1, 0, "3j1b");
      H += Hadronic(3, 1, -1);
      return H;
    }

    //  3j ge2b
    if(index == 5){
      Hadronic H(3, -2, 0, "3jge2b");
      H += Hadronic(3, -2, -1);
      return H;
    }

    // ge3j 0b
    if(index == 6){
      Hadronic H(-3, 0, 0, "ge3j0b");
      H += Hadronic(-3, 0, -1);
      return H;
    }

    // ge3j ge1b
    if(index == 7){
      Hadronic H(-3, -1, 0, "ge3jge1b");
      H += Hadronic(-3, -1, -1);
      return H;
    }

    //  ge3j 1b
    if(index == 8){
      Hadronic H(-3, 1, 0, "ge3j1b");
      H += Hadronic(-3, 1, -1);
      return H;
    }

    //  ge3j ge2b
    if(index == 9){
      Hadronic H(-3, -2, 0, "ge3jge2b");
      H += Hadronic(-3, -2, -1);
      return H;
    } 
  }

  /////////
  // 4 jet
  /////////
  if(Njet == 4){
    // C12 fully inclusive ge4 jet
    if(index == 0){
      Hadronic H(-4, 0, 0, "ge4j");
      H += Hadronic(-4, 0, -1);
      H += Hadronic(-4, -1, 0);
      H += Hadronic(-4, -1, -1);
      return H;
    }
    
    // fully inclusive exactly 4 jet
    if(index == 1){
      Hadronic H(4, 0, 0, "4j");
      H += Hadronic(4, 0, -1);
      H += Hadronic(4, -1, 0);
      H += Hadronic(4, -1, -1);
      return H;
    }

    // 4j 0b
    if(index == 2){
      Hadronic H(4, 0, 0, "4j0b");
      H += Hadronic(4, 0, -1);
      return H;
    }

    // 4j ge1b
    if(index == 3){
      Hadronic H(4, -1, 0, "4jge1b");
      H += Hadronic(4, -1, -1);
      return H;
    }

    //  4j 1b
    if(index == 4){
      Hadronic H(4, 1, 0, "4j1b");
      H += Hadronic(4, 1, -1);
      return H;
    }

    //  4j ge2b
    if(index == 5){
      Hadronic H(4, -2, 0, "4jge2b");
      H += Hadronic(4, -2, -1);
      return H;
    }

    // ge4j 0b
    if(index == 6){
      Hadronic H(-4, 0, 0, "ge4j0b");
      H += Hadronic(-4, 0, -1);
      return H;
    }

    // ge4j ge1b
    if(index == 7){
      Hadronic H(-4, -1, 0, "ge4jge1b");
      H += Hadronic(-4, -1, -1);
      return H;
    }

    //  ge4j 1b
    if(index == 8){
      Hadronic H(-4, 1, 0, "ge4j1b");
      H += Hadronic(-4, 1, -1);
      return H;
    }

    //  ge4j ge2b
    if(index == 9){
      Hadronic H(-4, -2, 0, "ge4jge2b");
      H += Hadronic(-4, -2, -1);
      return H;
    }
  }

  /////////
  // 5 jet
  /////////
  if(Njet == 5){
    // fully inclusive ge5 jet
    if(index == 0){
      Hadronic H(-5, 0, 0, "ge5j");
      H += Hadronic(-5, 0, -1);
      H += Hadronic(-5, -1, 0);
      H += Hadronic(-5, -1, -1);
      return H;
    }
    
    // fully inclusive exactly 5 jet
    if(index == 1){
      Hadronic H(5, 0, 0, "5j");
      H += Hadronic(5, 0, -1);
      H += Hadronic(5, -1, 0);
      H += Hadronic(5, -1, -1);
      return H;
    }

    // 5j 0b
    if(index == 2){
      Hadronic H(5, 0, 0, "5j0b");
      H += Hadronic(5, 0, -1);
      return H;
    }

    // 5j ge1b
    if(index == 3){
      Hadronic H(5, -1, 0, "5jge1b");
      H += Hadronic(5, -1, -1);
      return H;
    }

    //  5j 1b
    if(index == 4){
      Hadronic H(5, 1, 0, "5j1b");
      H += Hadronic(5, 1, -1);
      return H;
    }

    //  5j ge2b
    if(index == 5){
      Hadronic H(5, -2, 0, "5jge2b");
      H += Hadronic(5, -2, -1);
      return H;
    }

    // ge5j 0b
    if(index == 6){
      Hadronic H(-5, 0, 0, "ge5j0b");
      H += Hadronic(-5, 0, -1);
      return H;
    }

    // ge5j ge1b
    if(index == 7){
      Hadronic H(-5, -1, 0, "ge5jge1b");
      H += Hadronic(-5, -1, -1);
      return H;
    }

    //  ge5j 1b
    if(index == 8){
      Hadronic H(-5, 1, 0, "ge5j1b");
      H += Hadronic(-5, 1, -1);
      return H;
    }

    //  ge5j ge2b
    if(index == 9){
      Hadronic H(-5, -2, 0, "ge5jge2b");
      H += Hadronic(-5, -2, -1);
      return H;
    }
    
  }

  // default inclusive
  Hadronic H(0, 0, 0, "incl");
  H += Hadronic(0, 0, -1);
  H += Hadronic(-1, 0, 0);
  H += Hadronic(-1, 0, -1);
  H += Hadronic(-1, -1, 0);
  H += Hadronic(-1, -1, -1);
  return H;

}

float gamTval = 0.0;

CategoryList CategoryTool::GetCategories_0L(bool maskSR) const {
  CategoryList Cats_0L;
  cout << "Creating 0L Categories" << endl;

  // define leptonic categories
  // 0 L
  Leptonic L_0(LepList(), LepList(), "0");

  // Get hadronic ISR region lists - inclusive and b-separated
  vector<Hadronic> H_ISR_noB;
  vector<Hadronic> H_ISR_B;
  H_ISR_noB.push_back(GetHadronicRegion(1, 0)); // ge1j
  H_ISR_B.push_back(GetHadronicRegion(1, 6)); // ge1j0b
  H_ISR_B.push_back(GetHadronicRegion(1, 7)); // ge1jge1b
  
  // Get PTISR binning
  vector<double> PTISR_incl;
  PTISR_incl.push_back(300.);

  vector<double> PTISR_bin2;
  PTISR_bin2.push_back(300.);
  PTISR_bin2.push_back(400.);

  vector<double> PTISR_bin3;
  PTISR_bin3.push_back(300.);
  PTISR_bin3.push_back(450.);

  vector<double> PTISR_bin4;
  PTISR_bin4.push_back(350.);
  PTISR_bin4.push_back(500.);
 
  // Get SV binning
  vector<double> SV_eta;
  SV_eta.push_back(0.);
  SV_eta.push_back(1.5);

  // Get gammaT binning
  vector<double> gamT0;
  gamT0.push_back(0.);
  vector<double> gamT;
  gamT.push_back(gamTval);//Change lower cut on gammaT (originally 0)
  gamT.push_back(0.5);

  ///////////////////////////////////////////////

  cout << "Building 0L 0Sjet SV object regions" << endl;
  CategoryList Cats_0L_1SV;
  CategoryList Cats_0L_ge2SV;
  Cats_0L_1SV += Category(L_0, "Ch0L");
  Cats_0L_ge2SV += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_1SV;
  vector<Hadronic> H_0L_ge2SV;
  H_0L_1SV.push_back(GetHadronicRegion(0, 3)); // 0j1sv
  H_0L_ge2SV.push_back(GetHadronicRegion(0, 4)); // 0jge2sv

  Cats_0L_1SV = Cats_0L_1SV.CreateFitBinRegions(GetFitBin(0, 0, maskSR));
  Cats_0L_1SV = Cats_0L_1SV.CreateHadronicSRegions(H_0L_1SV);
  Cats_0L_1SV = Cats_0L_1SV.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_1SV = Cats_0L_1SV.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_0L_1SV = Cats_0L_1SV.CreateGenericRegions("gamT", gamT0);
  Cats_0L_1SV = Cats_0L_1SV.CreateGenericRegions("SVeta", SV_eta);

  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateFitBinRegions(GetFitBin(0, 0, maskSR));
  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateHadronicSRegions(H_0L_ge2SV);
  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateGenericRegions("gamT", gamT0);
  Cats_0L_ge2SV = Cats_0L_ge2SV.CreateGenericRegions("SVeta", SV_eta);

  Cats_0L += Cats_0L_1SV;
  Cats_0L += Cats_0L_ge2SV;

  ///////////////////////////////////////////////

  cout << "Building 0L 1Sjet 0 SV object regions" << endl;
  CategoryList Cats_0L_1j0SV;
  Cats_0L_1j0SV += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_1j0SV;
  
  H_0L_1j0SV.push_back(GetHadronicRegion(1, 2)); // 1j0b0sv
  H_0L_1j0SV.push_back(GetHadronicRegion(1, 3)); // 1j1b0sv
 

  Cats_0L_1j0SV = Cats_0L_1j0SV.CreateFitBinRegions(GetFitBin(0, 1, maskSR));
  Cats_0L_1j0SV = Cats_0L_1j0SV.CreateHadronicSRegions(H_0L_1j0SV);
  Cats_0L_1j0SV = Cats_0L_1j0SV.CreateHadronicISRRegions(H_ISR_B);
  Cats_0L_1j0SV = Cats_0L_1j0SV.CreateGenericRegions("PTISR", VD().a(400.).a(550.));

  Cats_0L += Cats_0L_1j0SV;
  
  ///////////////////////////////////////////////

  cout << "Building 0L 1Sjet 1 SV object regions" << endl;
  CategoryList Cats_0L_1j1SV;
  Cats_0L_1j1SV += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_1j1SV;
  
  H_0L_1j1SV.push_back(GetHadronicRegion(1, 11)); // 1jge1sv

  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateFitBinRegions(GetFitBin(0, 1, maskSR));
  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateHadronicSRegions(H_0L_1j1SV);
  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateGenericRegions("PTISR", VD().a(400.));
  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateGenericRegions("gamT", gamT0);
  Cats_0L_1j1SV = Cats_0L_1j1SV.CreateGenericRegions("SVeta", SV_eta);
  
  Cats_0L += Cats_0L_1j1SV;
  
  ///////////////////////////////////////////////

  cout << "Building 0L 2j 0/1/2b object regions" << endl;
  CategoryList Cats_0L_2jS;
  CategoryList Cats_0L_2jS_2b;
  Cats_0L_2jS    += Category(L_0, "Ch0L");
  Cats_0L_2jS_2b += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_2jS;
  vector<Hadronic> H_0L_2jS_2b;
  H_0L_2jS.push_back(GetHadronicRegion(2, 2)); // 2j0b
  H_0L_2jS.push_back(GetHadronicRegion(2, 4)); // 2j1b
  H_0L_2jS_2b.push_back(GetHadronicRegion(2, 5)); // 2j2b

  Cats_0L_2jS = Cats_0L_2jS.CreateFitBinRegions(GetFitBin(0, 2, maskSR));
  Cats_0L_2jS = Cats_0L_2jS.CreateHadronicSRegions(H_0L_2jS);
  Cats_0L_2jS = Cats_0L_2jS.CreateHadronicISRRegions(H_ISR_B);
  Cats_0L_2jS = Cats_0L_2jS.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_2jS = Cats_0L_2jS.CreateGenericRegions("gamT", gamT);
  
  Cats_0L_2jS_2b = Cats_0L_2jS_2b.CreateFitBinRegions(GetFitBin(0, 2, maskSR));
  Cats_0L_2jS_2b = Cats_0L_2jS_2b.CreateHadronicSRegions(H_0L_2jS_2b);
  Cats_0L_2jS_2b = Cats_0L_2jS_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_2jS_2b = Cats_0L_2jS_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_2jS_2b = Cats_0L_2jS_2b.CreateGenericRegions("gamT", gamT);

  Cats_0L += Cats_0L_2jS;
  Cats_0L += Cats_0L_2jS_2b;

  ///////////////////////////////////////////////

  cout << "Building 0L 3j 0/1/2b object regions" << endl;
  CategoryList Cats_0L_3jS;
  CategoryList Cats_0L_3jS_2b;
  Cats_0L_3jS    += Category(L_0, "Ch0L");
  Cats_0L_3jS_2b += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_3jS;
  vector<Hadronic> H_0L_3jS_2b;
  H_0L_3jS.push_back(GetHadronicRegion(3, 2)); // 3j0b
  H_0L_3jS.push_back(GetHadronicRegion(3, 4)); // 3j1b
  H_0L_3jS_2b.push_back(GetHadronicRegion(3, 5)); // 3j2b

  Cats_0L_3jS = Cats_0L_3jS.CreateFitBinRegions(GetFitBin(0, 3, maskSR));
  Cats_0L_3jS = Cats_0L_3jS.CreateHadronicSRegions(H_0L_3jS);
  Cats_0L_3jS = Cats_0L_3jS.CreateHadronicISRRegions(H_ISR_B);
  Cats_0L_3jS = Cats_0L_3jS.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_3jS = Cats_0L_3jS.CreateGenericRegions("gamT", gamT);

  Cats_0L_3jS_2b = Cats_0L_3jS_2b.CreateFitBinRegions(GetFitBin(0, 3, maskSR));
  Cats_0L_3jS_2b = Cats_0L_3jS_2b.CreateHadronicSRegions(H_0L_3jS_2b);
  Cats_0L_3jS_2b = Cats_0L_3jS_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_3jS_2b = Cats_0L_3jS_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_3jS_2b = Cats_0L_3jS_2b.CreateGenericRegions("gamT", gamT);

  Cats_0L += Cats_0L_3jS;
  Cats_0L += Cats_0L_3jS_2b;

  ///////////////////////////////////////////////

  cout << "Building 0L 4j 0/1/2b object regions" << endl;
  CategoryList Cats_0L_4jS;
  CategoryList Cats_0L_4jS_2b;
  Cats_0L_4jS    += Category(L_0, "Ch0L");
  Cats_0L_4jS_2b += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_4jS;
  vector<Hadronic> H_0L_4jS_2b;
  H_0L_4jS.push_back(GetHadronicRegion(4, 2)); // 4j0b
  H_0L_4jS.push_back(GetHadronicRegion(4, 4)); // 4j1b
  H_0L_4jS_2b.push_back(GetHadronicRegion(4, 5)); // 4j2b

  Cats_0L_4jS = Cats_0L_4jS.CreateFitBinRegions(GetFitBin(0, 4, maskSR));
  Cats_0L_4jS = Cats_0L_4jS.CreateHadronicSRegions(H_0L_4jS);
  Cats_0L_4jS = Cats_0L_4jS.CreateHadronicISRRegions(H_ISR_B);
  Cats_0L_4jS = Cats_0L_4jS.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_4jS = Cats_0L_4jS.CreateGenericRegions("gamT", gamT);

  Cats_0L_4jS_2b = Cats_0L_4jS_2b.CreateFitBinRegions(GetFitBin(0, 4, maskSR));
  Cats_0L_4jS_2b = Cats_0L_4jS_2b.CreateHadronicSRegions(H_0L_4jS_2b);
  Cats_0L_4jS_2b = Cats_0L_4jS_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_4jS_2b = Cats_0L_4jS_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_0L_4jS_2b = Cats_0L_4jS_2b.CreateGenericRegions("gamT", gamT);

  Cats_0L += Cats_0L_4jS;
  Cats_0L += Cats_0L_4jS_2b;

  ///////////////////////////////////////////////

  cout << "Building 0L 5j 0/ge1b object regions" << endl;
  CategoryList Cats_0L_5jS;
  CategoryList Cats_0L_5jS_2b;
  Cats_0L_5jS += Category(L_0, "Ch0L");
  Cats_0L_5jS_2b += Category(L_0, "Ch0L");
  vector<Hadronic> H_0L_5jS;
  vector<Hadronic> H_0L_5jS_2b;
  H_0L_5jS.push_back(GetHadronicRegion(5, 6)); // ge5j0b
  //H_0L_5jS.push_back(GetHadronicRegion(5, 7)); // ge5jge1b
  H_0L_5jS.push_back(GetHadronicRegion(5, 8)); // ge5j1b
  H_0L_5jS_2b.push_back(GetHadronicRegion(5, 9)); // ge5jge2b

  Cats_0L_5jS = Cats_0L_5jS.CreateFitBinRegions(GetFitBin(0, 5, maskSR));
  Cats_0L_5jS = Cats_0L_5jS.CreateHadronicSRegions(H_0L_5jS);
  Cats_0L_5jS = Cats_0L_5jS.CreateHadronicISRRegions(H_ISR_B);
  Cats_0L_5jS = Cats_0L_5jS.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  //Cats_0L_5jS = Cats_0L_5jS.CreateGenericRegions("gamT", gamT);
  Cats_0L_5jS = Cats_0L_5jS.CreateGenericRegions("gamT",gamT0);

  Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateFitBinRegions(GetFitBin(0, 5, maskSR));
  Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateHadronicSRegions(H_0L_5jS_2b);
  Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  //Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateGenericRegions("gamT", gamT);
  Cats_0L_5jS_2b = Cats_0L_5jS_2b.CreateGenericRegions("gamT", gamT0);

  Cats_0L += Cats_0L_5jS;
  Cats_0L += Cats_0L_5jS_2b;

  return Cats_0L;
}

CategoryList CategoryTool::GetCategories_1L(bool maskSR) const {
  CategoryList Cats_1L;
  cout << "Creating 1L Categories" << endl;

  // define leptons
  Lep elp(kElec, kPos);
  Lep elm(kElec, kNeg);
  Lep mup(kMuon, kPos);
  Lep mum(kMuon, kNeg);

  ///////////////////////
  // separated by charge
  ///////////////////////
  Leptonic L_elp_0(LepList().Add(elp), LepList(), "elp");
  Leptonic L_elm_0(LepList().Add(elm), LepList(), "elm");
  Leptonic L_mup_0(LepList().Add(mup), LepList(), "mup");
  Leptonic L_mum_0(LepList().Add(mum), LepList(), "mum");
  CategoryList Cats_1Lcharge;
  Cats_1Lcharge += Category(L_elp_0, "Ch1L"); 
  Cats_1Lcharge += Category(L_elm_0, "Ch1L");
  Cats_1Lcharge += Category(L_mup_0, "Ch1L");
  Cats_1Lcharge += Category(L_mum_0, "Ch1L");

  // flavor integrated
  Leptonic L_lp_0(LepList().Add(elp), LepList(), "lp");
  L_lp_0 += Leptonic(LepList().Add(mup), LepList());
  Leptonic L_lm_0(LepList().Add(elm), LepList(), "lm");
  L_lm_0 += Leptonic(LepList().Add(mum), LepList());
  CategoryList Cats_1Lcharge_noflavor;
  Cats_1Lcharge_noflavor += Category(L_lp_0, "Ch1L"); 
  Cats_1Lcharge_noflavor += Category(L_lm_0, "Ch1L");

  //////////////////////
  // charge inclusive
  //////////////////////
  Leptonic L_elpm_0(LepList().Add(elp), LepList(), "elpm");
  L_elpm_0 += Leptonic(LepList().Add(elm), LepList());
  Leptonic L_mupm_0(LepList().Add(mup), LepList(), "mupm");
  L_mupm_0 += Leptonic(LepList().Add(mum), LepList());
  CategoryList Cats_1Lincl;
  Cats_1Lincl += Category(L_elpm_0, "Ch1L"); 
  Cats_1Lincl += Category(L_mupm_0, "Ch1L");

  // flavor integrated
  Leptonic L_lpm_0(LepList().Add(elp), LepList(), "lpm");
  L_lpm_0 += Leptonic(LepList().Add(elm), LepList());
  L_lpm_0 += Leptonic(LepList().Add(mup), LepList());
  L_lpm_0 += Leptonic(LepList().Add(mum), LepList());
  
  CategoryList Cats_1Lincl_noflavor;
  Cats_1Lincl_noflavor += Category(L_lpm_0, "Ch1L"); 
 
  // Get hadronic ISR region lists - inclusive and b-separated
  vector<Hadronic> H_ISR_noB;
  vector<Hadronic> H_ISR_B;
  H_ISR_noB.push_back(GetHadronicRegion(1, 0)); // ge1j
  H_ISR_B.push_back(GetHadronicRegion(1, 6)); // ge1j0b
  H_ISR_B.push_back(GetHadronicRegion(1, 7)); // ge1jge1b

   // Get PTISR binning
  vector<double> PTISR_incl;
  PTISR_incl.push_back(300.);

  vector<double> PTISR_bin2;
  PTISR_bin2.push_back(300.);
  PTISR_bin2.push_back(400.);

  vector<double> PTISR_bin3;
  PTISR_bin3.push_back(300.);
  PTISR_bin3.push_back(450.);

  vector<double> PTISR_bin4;
  PTISR_bin4.push_back(350.);
  PTISR_bin4.push_back(500.);
 
  // Get SV binning
  vector<double> SV_eta;
  SV_eta.push_back(0.);
  SV_eta.push_back(1.5);

  // Get gammaT binning
  vector<double> gamT0;
  gamT0.push_back(0.);
  vector<double> gamT;
  gamT.push_back(gamTval);
  gamT.push_back(0.5);

  // Define Lepton ID-based categories:
  vector<LepIDsList> IDs_all;
  vector<LepIDsList> IDs_gold;
  vector<LepIDsList> IDs_slvrbron;
  
  LepIDsList ID_gold("gold");
  LepIDsList ID_silver("slvr");
  LepIDsList ID_bronze("bron");

  ID_gold += LepIDs().a(kGold);
  ID_silver += LepIDs().a(kSilver);
  ID_bronze += LepIDs().a(kBronze);

  IDs_all.push_back(ID_gold);
  IDs_all.push_back(ID_silver);
  IDs_all.push_back(ID_bronze);
  IDs_slvrbron.push_back(ID_silver);
  IDs_slvrbron.push_back(ID_bronze);
  IDs_gold.push_back(ID_gold);

  ///////////////////////////////////////////////
  
  cout << "Building 1L 0S object regions" << endl;
  // flavor and charge separated for gold
  CategoryList Cats_1L_0S_gold = Cats_1Lcharge;
  // flavor separated for silver/bronze
  CategoryList Cats_1L_0S_slvrbron = Cats_1Lincl;
  vector<Hadronic> H_1L_0S;
  H_1L_0S.push_back(GetHadronicRegion(0, 1)); // 0j0sv

  Cats_1L_0S_gold = Cats_1L_0S_gold.CreateFitBinRegions(GetFitBin(1, 0, maskSR)); 
  Cats_1L_0S_gold = Cats_1L_0S_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_0S_gold = Cats_1L_0S_gold.CreateHadronicSRegions(H_1L_0S);
  Cats_1L_0S_gold = Cats_1L_0S_gold.CreateHadronicISRRegions(H_ISR_B);
  Cats_1L_0S_gold = Cats_1L_0S_gold.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
 
  Cats_1L_0S_slvrbron = Cats_1L_0S_slvrbron.CreateFitBinRegions(GetFitBin(1, 0, maskSR));
  Cats_1L_0S_slvrbron = Cats_1L_0S_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_0S_slvrbron = Cats_1L_0S_slvrbron.CreateHadronicSRegions(H_1L_0S);
  Cats_1L_0S_slvrbron = Cats_1L_0S_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_0S_slvrbron = Cats_1L_0S_slvrbron.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
 
  Cats_1L += Cats_1L_0S_gold;
  Cats_1L += Cats_1L_0S_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 1L 0Sjet SV object regions" << endl;
  // charge separated, flavor integrated for gold
  CategoryList Cats_1L_SV0S_gold   = Cats_1Lcharge_noflavor;
  // flavor separated, charge integrated
  CategoryList Cats_1L_SV0S_slvrbron = Cats_1Lincl;
  vector<Hadronic> H_1L_SV0S;
  H_1L_SV0S.push_back(GetHadronicRegion(0, 2)); // 0jge1sv

  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateFitBinRegions(GetFitBin(1, 0, maskSR));
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateHadronicSRegions(H_1L_SV0S);
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateGenericRegions("gamT", gamT0);
  Cats_1L_SV0S_gold = Cats_1L_SV0S_gold.CreateGenericRegions("SVeta", SV_eta);

  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateFitBinRegions(GetFitBin(1, 0, maskSR));
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateHadronicSRegions(H_1L_SV0S);
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateGenericRegions("gamT", gamT0);
  Cats_1L_SV0S_slvrbron = Cats_1L_SV0S_slvrbron.CreateGenericRegions("SVeta", SV_eta);

  Cats_1L += Cats_1L_SV0S_gold;
  Cats_1L += Cats_1L_SV0S_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 1L 1j 0/1b 0sv object regions" << endl;
  // flavor and charge separated for gold
  CategoryList Cats_1L_1j0svS_gold   = Cats_1Lincl_noflavor;
  // charge inclusive for silver/bron
  CategoryList Cats_1L_1j0svS_slvrbron = Cats_1Lincl;
  
  vector<Hadronic> H_1L_1j0svS;
  H_1L_1j0svS.push_back(GetHadronicRegion(1, 2)); // 1j0b0sv
  H_1L_1j0svS.push_back(GetHadronicRegion(1, 3)); // 1j1b0sv
  vector<Hadronic> H_1L_1j0svS_bkg;
  H_1L_1j0svS_bkg.push_back(GetHadronicRegion(1, 10)); // 1j0sv

  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateFitBinRegions(GetFitBin(1, 1, maskSR));
  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateHadronicSRegions(H_1L_1j0svS);
  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateHadronicISRRegions(H_ISR_B);
  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_1j0svS_gold = Cats_1L_1j0svS_gold.CreateGenericRegions("gamT", gamT);

  Cats_1L_1j0svS_slvrbron = Cats_1L_1j0svS_slvrbron.CreateFitBinRegions(GetFitBin(1, 1, maskSR));
  Cats_1L_1j0svS_slvrbron = Cats_1L_1j0svS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_1j0svS_slvrbron = Cats_1L_1j0svS_slvrbron.CreateHadronicSRegions(H_1L_1j0svS_bkg);
  Cats_1L_1j0svS_slvrbron = Cats_1L_1j0svS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_1j0svS_slvrbron = Cats_1L_1j0svS_slvrbron.CreateGenericRegions("PTISR", VD().a(350.).a(500.));

  Cats_1L += Cats_1L_1j0svS_gold;
  Cats_1L += Cats_1L_1j0svS_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 1L 1j 0/1b ge1sv object regions" << endl;
  // charge separated, flavor integrated for gold
  CategoryList Cats_1L_1jge1svS_gold   = Cats_1Lincl_noflavor;
  // flavor separated, charge integrted for silver/bron
  CategoryList Cats_1L_1jge1svS_slvrbron = Cats_1Lincl;

  vector<Hadronic> H_1L_1jge1svS;
  // H_1L_1jge1svS.push_back(GetHadronicRegion(1, 4)); // 1j0bge1sv
  // H_1L_1jge1svS.push_back(GetHadronicRegion(1, 5)); // 1j1bge1sv
  H_1L_1jge1svS.push_back(GetHadronicRegion(1, 11)); // 1jge1sv

  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateFitBinRegions(GetFitBin(1, 1, maskSR));
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateHadronicSRegions(H_1L_1jge1svS);
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateGenericRegions("gamT", gamT0);
  Cats_1L_1jge1svS_gold = Cats_1L_1jge1svS_gold.CreateGenericRegions("SVeta", SV_eta);

  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateFitBinRegions(GetFitBin(1, 1, maskSR));
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateHadronicSRegions(H_1L_1jge1svS);
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateGenericRegions("PTISR", VD().a(350.));
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateGenericRegions("gamT", gamT0);
  Cats_1L_1jge1svS_slvrbron = Cats_1L_1jge1svS_slvrbron.CreateGenericRegions("SVeta", SV_eta);

  Cats_1L += Cats_1L_1jge1svS_gold;
  Cats_1L += Cats_1L_1jge1svS_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 1L 2j 0/1/2b object regions" << endl;
  // flavor and charge integrated for gold
  CategoryList Cats_1L_2jS_gold    = Cats_1Lincl_noflavor;
  CategoryList Cats_1L_2jS_gold_2b = Cats_1Lincl_noflavor;
  // charge separated for silver/bron
  CategoryList Cats_1L_2jS_slvrbron = Cats_1Lincl;
 
  vector<Hadronic> H_1L_2jS;
  vector<Hadronic> H_1L_2jS_2b;
  H_1L_2jS.push_back(GetHadronicRegion(2, 2)); // 2j0b
  H_1L_2jS.push_back(GetHadronicRegion(2, 4)); // 2j1b
  H_1L_2jS_2b.push_back(GetHadronicRegion(2, 5)); // 2j2b

  vector<Hadronic> H_1L_2jS_bkg;
  // H_1L_2jS_bkg.push_back(GetHadronicRegion(2, 2)); // 2j0b
  // H_1L_2jS_bkg.push_back(GetHadronicRegion(2, 3)); // 2jge1b
  H_1L_2jS_bkg.push_back(GetHadronicRegion(2, 1)); // 2j

  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateFitBinRegions(GetFitBin(1, 2, maskSR));
  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateHadronicSRegions(H_1L_2jS);
  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateHadronicISRRegions(H_ISR_B);
  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_2jS_gold = Cats_1L_2jS_gold.CreateGenericRegions("gamT", gamT);

  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateFitBinRegions(GetFitBin(1, 2, maskSR));
  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateHadronicSRegions(H_1L_2jS_2b);
  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_2jS_gold_2b = Cats_1L_2jS_gold_2b.CreateGenericRegions("gamT", gamT);

  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateFitBinRegions(GetFitBin(1, 2, maskSR));
  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateHadronicSRegions(H_1L_2jS_bkg);
  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_2jS_slvrbron = Cats_1L_2jS_slvrbron.CreateGenericRegions("gamT", gamT);

  Cats_1L += Cats_1L_2jS_gold;
  Cats_1L += Cats_1L_2jS_gold_2b;
  Cats_1L += Cats_1L_2jS_slvrbron;

  ///////////////////////////////////////////////

  cout << "Building 1L 3j 0/1/2b object regions" << endl;
  // flavor and charge integrated for gold
  CategoryList Cats_1L_3jS_gold    = Cats_1Lincl_noflavor;
  CategoryList Cats_1L_3jS_gold_2b = Cats_1Lincl_noflavor;
  // flavor separated for silver/bron
  CategoryList Cats_1L_3jS_slvrbron = Cats_1Lincl;
 
  vector<Hadronic> H_1L_3jS;
  vector<Hadronic> H_1L_3jS_2b;
  H_1L_3jS.push_back(GetHadronicRegion(3, 2)); // 3j0b
  H_1L_3jS.push_back(GetHadronicRegion(3, 4)); // 3j1b
  H_1L_3jS_2b.push_back(GetHadronicRegion(3, 5)); // 3jge2b

  vector<Hadronic> H_1L_3jS_bkg;
  // H_1L_3jS_bkg.push_back(GetHadronicRegion(3, 2)); // 3j0b
  // H_1L_3jS_bkg.push_back(GetHadronicRegion(3, 3)); // 3jge1b
  H_1L_3jS_bkg.push_back(GetHadronicRegion(3, 1)); // 3j

  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateFitBinRegions(GetFitBin(1, 3, maskSR));
  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateHadronicSRegions(H_1L_3jS);
  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateHadronicISRRegions(H_ISR_B);
  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_3jS_gold = Cats_1L_3jS_gold.CreateGenericRegions("gamT", gamT);

  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateFitBinRegions(GetFitBin(1, 3, maskSR));
  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateHadronicSRegions(H_1L_3jS_2b);
  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_3jS_gold_2b = Cats_1L_3jS_gold_2b.CreateGenericRegions("gamT", gamT);

  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateFitBinRegions(GetFitBin(1, 3, maskSR));
  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateHadronicSRegions(H_1L_3jS_bkg);
  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_3jS_slvrbron = Cats_1L_3jS_slvrbron.CreateGenericRegions("gamT", gamT);

  Cats_1L += Cats_1L_3jS_gold;
  Cats_1L += Cats_1L_3jS_gold_2b;
  Cats_1L += Cats_1L_3jS_slvrbron;

  ///////////////////////////////////////////////

  cout << "Building 1L 4j 0/1b/2b object regions" << endl;
  // flavor and charge integrated for gold
  CategoryList Cats_1L_4jS_gold    = Cats_1Lincl_noflavor;
  CategoryList Cats_1L_4jS_gold_2b = Cats_1Lincl_noflavor;
  // flavor separated for silver/bron
  CategoryList Cats_1L_4jS_slvrbron = Cats_1Lincl;
 
  vector<Hadronic> H_1L_4jS;
  vector<Hadronic> H_1L_4jS_2b;
  H_1L_4jS.push_back(GetHadronicRegion(4, 6)); // ge4j0b
  //H_1L_4jS.push_back(GetHadronicRegion(4, 7)); // ge4jge1b
  H_1L_4jS.push_back(GetHadronicRegion(4, 8)); // ge4j1b
  H_1L_4jS_2b.push_back(GetHadronicRegion(4, 9)); // ge4jge2b
  vector<Hadronic> H_1L_4jS_bkg;
  H_1L_4jS_bkg.push_back(GetHadronicRegion(4, 0)); // ge4j

  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateFitBinRegions(GetFitBin(1, 4, maskSR));
  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateHadronicSRegions(H_1L_4jS);
  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateHadronicISRRegions(H_ISR_B);
  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_4jS_gold = Cats_1L_4jS_gold.CreateGenericRegions("gamT", gamT);

  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateFitBinRegions(GetFitBin(1, 4, maskSR));
  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateLeptonIDRegions(IDs_gold);
  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateHadronicSRegions(H_1L_4jS_2b);
  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_4jS_gold_2b = Cats_1L_4jS_gold_2b.CreateGenericRegions("gamT", gamT);

  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateFitBinRegions(GetFitBin(1, 4, maskSR));
  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateHadronicSRegions(H_1L_4jS_bkg);
  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateGenericRegions("PTISR", VD().a(350.).a(500.));
  Cats_1L_4jS_slvrbron = Cats_1L_4jS_slvrbron.CreateGenericRegions("gamT", gamT);

  Cats_1L += Cats_1L_4jS_gold;
  Cats_1L += Cats_1L_4jS_gold_2b;
  Cats_1L += Cats_1L_4jS_slvrbron;

  ///////////////////////////////////////////////

  return Cats_1L;
}

CategoryList CategoryTool::GetCategories_2L(bool maskSR) const {
  CategoryList Cats_2L;
  cout << "Creating 2L Categories" << endl;

  // define leptons
  Lep elp(kElec, kPos);
  Lep elm(kElec, kNeg);
  Lep mup(kMuon, kPos);
  Lep mum(kMuon, kNeg);

  // 2 L, opposite hemisphere, charge separated
  Leptonic    L_OS_el_el(LepList().Add(elp), LepList().Add(elm), "OSel^el");
  Leptonic    L_SS_el_el(LepList().Add(elp), LepList().Add(elp), "SSel^el");
  L_SS_el_el += Leptonic(LepList().Add(elm), LepList().Add(elm));
  Leptonic    L_OS_mu_mu(LepList().Add(mup), LepList().Add(mum), "OSmu^mu");
  Leptonic    L_SS_mu_mu(LepList().Add(mup), LepList().Add(mup), "SSmu^mu");
  L_SS_mu_mu += Leptonic(LepList().Add(mum), LepList().Add(mum));
  Leptonic    L_OS_el_mu(LepList().Add(elp), LepList().Add(mum), "OSel^mu");
  L_OS_el_mu += Leptonic(LepList().Add(elm), LepList().Add(mup));
  Leptonic    L_SS_el_mu(LepList().Add(elp), LepList().Add(mup), "SSel^mu");
  L_SS_el_mu += Leptonic(LepList().Add(elm), LepList().Add(mum));
  // 2 L, same hemisphere, charge separated
  Leptonic    L_OS_elel_0(LepList().Add(elp).Add(elm), LepList(), "OSelel^0");
  Leptonic    L_SS_elel_0(LepList().Add(elp).Add(elp), LepList(), "SSelel^0");
  L_SS_elel_0 += Leptonic(LepList().Add(elm).Add(elm), LepList());
  Leptonic    L_OS_mumu_0(LepList().Add(mup).Add(mum), LepList(), "OSmumu^0");
  Leptonic    L_SS_mumu_0(LepList().Add(mup).Add(mup), LepList(), "SSmumu^0");
  L_SS_mumu_0 += Leptonic(LepList().Add(mum).Add(mum), LepList());
  Leptonic    L_OS_elmu_0(LepList().Add(elp).Add(mum), LepList(), "OSelmu^0");
  L_OS_elmu_0 += Leptonic(LepList().Add(elm).Add(mup), LepList());
  Leptonic    L_SS_elmu_0(LepList().Add(elp).Add(mup), LepList(), "SSelmu^0");
  L_SS_elmu_0 += Leptonic(LepList().Add(elm).Add(mum), LepList());
  CategoryList Cats_2Lchargehem;
  Cats_2Lchargehem += Category(L_OS_el_el, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_el_el, "Ch2L");
  Cats_2Lchargehem += Category(L_OS_mu_mu, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_mu_mu, "Ch2L");
  Cats_2Lchargehem += Category(L_OS_el_mu, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_el_mu, "Ch2L");
  Cats_2Lchargehem += Category(L_OS_elel_0, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_elel_0, "Ch2L");
  Cats_2Lchargehem += Category(L_OS_mumu_0, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_mumu_0, "Ch2L");
  Cats_2Lchargehem += Category(L_OS_elmu_0, "Ch2L");
  Cats_2Lchargehem += Category(L_SS_elmu_0, "Ch2L");

  // 2 L, hemipshere separated only
  Leptonic    L_el_el(LepList().Add(elp), LepList().Add(elm), "el^el");
  L_el_el += Leptonic(LepList().Add(elp), LepList().Add(elp));
  L_el_el += Leptonic(LepList().Add(elm), LepList().Add(elm));
  Leptonic    L_mu_mu(LepList().Add(mup), LepList().Add(mum), "mu^mu");
  L_mu_mu += Leptonic(LepList().Add(mup), LepList().Add(mup));
  L_mu_mu += Leptonic(LepList().Add(mum), LepList().Add(mum));
  Leptonic    L_el_mu(LepList().Add(elp), LepList().Add(mum), "el^mu");
  L_el_mu += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_el_mu += Leptonic(LepList().Add(elp), LepList().Add(mup));
  L_el_mu += Leptonic(LepList().Add(elm), LepList().Add(mum));
  Leptonic    L_elel_0(LepList().Add(elp).Add(elm), LepList(), "elel^0");
  L_elel_0 += Leptonic(LepList().Add(elp).Add(elp), LepList());
  L_elel_0 += Leptonic(LepList().Add(elm).Add(elm), LepList());
  Leptonic    L_mumu_0(LepList().Add(mup).Add(mum), LepList(), "mumu^0");
  L_mumu_0 += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_mumu_0 += Leptonic(LepList().Add(mum).Add(mum), LepList());
  Leptonic    L_elmu_0(LepList().Add(elp).Add(mum), LepList(), "elmu^0");
  L_elmu_0 += Leptonic(LepList().Add(elm).Add(mup), LepList());
  L_elmu_0 += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_elmu_0 += Leptonic(LepList().Add(elm).Add(mum), LepList());
  CategoryList Cats_2Lhem;
  Cats_2Lhem += Category(L_el_el, "Ch2L");
  Cats_2Lhem += Category(L_mu_mu, "Ch2L");
  Cats_2Lhem += Category(L_el_mu, "Ch2L");
  Cats_2Lhem += Category(L_elel_0, "Ch2L");
  Cats_2Lhem += Category(L_mumu_0, "Ch2L");
  Cats_2Lhem += Category(L_elmu_0, "Ch2L");
  
  // 2 L, charge separated only
  Leptonic    L_OS_elel(LepList().Add(elp), LepList().Add(elm), "OSelel");
  L_OS_elel += Leptonic(LepList().Add(elp).Add(elm), LepList());
  Leptonic    L_SS_elel(LepList().Add(elp), LepList().Add(elp), "SSelel");
  L_SS_elel += Leptonic(LepList().Add(elm), LepList().Add(elm));
  L_SS_elel += Leptonic(LepList().Add(elp).Add(elp), LepList());
  L_SS_elel += Leptonic(LepList().Add(elm).Add(elm), LepList());
  Leptonic    L_OS_mumu(LepList().Add(mup), LepList().Add(mum), "OSmumu");
  L_OS_mumu += Leptonic(LepList().Add(mup).Add(mum), LepList());
  Leptonic    L_SS_mumu(LepList().Add(mup), LepList().Add(mup), "SSmumu");
  L_SS_mumu += Leptonic(LepList().Add(mum), LepList().Add(mum));
  L_SS_mumu += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_SS_mumu += Leptonic(LepList().Add(mum).Add(mum), LepList());
  Leptonic    L_OS_elmu(LepList().Add(elp), LepList().Add(mum), "OSelmu");
  L_OS_elmu += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_OS_elmu += Leptonic(LepList().Add(elp).Add(mum), LepList());
  L_OS_elmu += Leptonic(LepList().Add(elm).Add(mup), LepList());
  Leptonic    L_SS_elmu(LepList().Add(elp), LepList().Add(mup), "SSelmu");
  L_SS_elmu += Leptonic(LepList().Add(elm), LepList().Add(mum));
  L_SS_elmu += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_SS_elmu += Leptonic(LepList().Add(elm).Add(mum), LepList());
  CategoryList Cats_2Lcharge;
  Cats_2Lcharge += Category(L_OS_elel, "Ch2L");
  Cats_2Lcharge += Category(L_SS_elel, "Ch2L");
  Cats_2Lcharge += Category(L_OS_mumu, "Ch2L");
  Cats_2Lcharge += Category(L_SS_mumu, "Ch2L");
  Cats_2Lcharge += Category(L_OS_elmu, "Ch2L");
  Cats_2Lcharge += Category(L_SS_elmu, "Ch2L");
  CategoryList Cats_2Lcharge_OS;
  Cats_2Lcharge_OS += Category(L_OS_elel, "Ch2L");
  Cats_2Lcharge_OS += Category(L_OS_mumu, "Ch2L");
  Cats_2Lcharge_OS += Category(L_OS_elmu, "Ch2L");
  CategoryList Cats_2Lcharge_SS;
  Cats_2Lcharge_SS += Category(L_SS_mumu, "Ch2L");
  Cats_2Lcharge_SS += Category(L_SS_elmu, "Ch2L");
  Cats_2Lcharge_SS += Category(L_SS_elel, "Ch2L");

  Leptonic   L_elelOSSH(LepList().Add(elp).Add(elm), LepList(), "OSSHelel");
  Leptonic   L_mumuOSSH(LepList().Add(mup).Add(mum), LepList(), "OSSHmumu");
  Leptonic   L_elelOSOH(LepList().Add(elp), LepList().Add(elm), "OSOHelel");
  Leptonic   L_mumuOSOH(LepList().Add(mup), LepList().Add(mum), "OSOHmumu");
  Leptonic   L_elelSS(LepList().Add(elp).Add(elp), LepList(), "SSelel");
  L_elelSS += Leptonic(LepList().Add(elm).Add(elm), LepList());
  L_elelSS += Leptonic(LepList().Add(elp), LepList().Add(elp));
  L_elelSS += Leptonic(LepList().Add(elm), LepList().Add(elm));
  Leptonic   L_mumuSS(LepList().Add(mum).Add(mum), LepList(), "SSmumu");
  L_mumuSS += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_mumuSS += Leptonic(LepList().Add(mup), LepList().Add(mup));
  L_mumuSS += Leptonic(LepList().Add(mum), LepList().Add(mum));
  Leptonic   L_elmuOS(LepList().Add(elp), LepList().Add(mum), "OSelmu");
  L_elmuOS += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_elmuOS += Leptonic(LepList().Add(elp).Add(mum), LepList());
  L_elmuOS += Leptonic(LepList().Add(elm).Add(mup), LepList());
  Leptonic   L_elmuSS(LepList().Add(elp), LepList().Add(mup), "SSelmu");
  L_elmuSS += Leptonic(LepList().Add(elm), LepList().Add(mum));
  L_elmuSS += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_elmuSS += Leptonic(LepList().Add(elm).Add(mum), LepList());
   
  CategoryList Cats_2Lmix_OS;
  CategoryList Cats_2Lmix_SS;
  Cats_2Lmix_OS += Category(L_elelOSSH, "Ch2L");
  Cats_2Lmix_OS += Category(L_mumuOSSH, "Ch2L");
  Cats_2Lmix_OS += Category(L_elelOSOH, "Ch2L");
  Cats_2Lmix_OS += Category(L_mumuOSOH, "Ch2L");
  Cats_2Lmix_OS += Category(L_elmuOS, "Ch2L");
  Cats_2Lmix_SS += Category(L_elmuSS, "Ch2L");
  Cats_2Lmix_SS += Category(L_elelSS, "Ch2L");
  Cats_2Lmix_SS += Category(L_mumuSS, "Ch2L");

  Leptonic   L_Zstar(LepList().Add(elp).Add(elm), LepList(), "Zstar");
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum), LepList());
  Leptonic   L_noZOS(LepList().Add(elp), LepList().Add(elm), "noZ");
  L_noZOS += Leptonic(LepList().Add(mup), LepList().Add(mum));
  L_noZOS += Leptonic(LepList().Add(elp), LepList().Add(mum));
  L_noZOS += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_noZOS += Leptonic(LepList().Add(elp).Add(mum), LepList());
  L_noZOS += Leptonic(LepList().Add(elm).Add(mup), LepList());
  Leptonic   L_SS(LepList().Add(elp).Add(elp), LepList(), "SS");
  L_SS += Leptonic(LepList().Add(elm).Add(elm), LepList());
  L_SS += Leptonic(LepList().Add(elp), LepList().Add(elp));
  L_SS += Leptonic(LepList().Add(elm), LepList().Add(elm));
  L_SS += Leptonic(LepList().Add(mum).Add(mum), LepList());
  L_SS += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_SS += Leptonic(LepList().Add(mup), LepList().Add(mup));
  L_SS += Leptonic(LepList().Add(mum), LepList().Add(mum));
  L_SS += Leptonic(LepList().Add(elp), LepList().Add(mup));
  L_SS += Leptonic(LepList().Add(elm), LepList().Add(mum));
  L_SS += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_SS += Leptonic(LepList().Add(elm).Add(mum), LepList());

  CategoryList Cats_2Lnoflavor_OS;
  CategoryList Cats_2Lnoflavor_SS;
  Cats_2Lnoflavor_OS += Category(L_Zstar, "Ch2L");
  Cats_2Lnoflavor_OS += Category(L_noZOS, "Ch2L");
  Cats_2Lnoflavor_SS += Category(L_SS, "Ch2L");
 
  // 2L, flavor separated only
  Leptonic    L_elel(LepList().Add(elp), LepList().Add(elm), "elel");
  L_elel += Leptonic(LepList().Add(elp).Add(elm), LepList());
  L_elel += Leptonic(LepList().Add(elp), LepList().Add(elp));
  L_elel += Leptonic(LepList().Add(elm), LepList().Add(elm));
  L_elel += Leptonic(LepList().Add(elp).Add(elp), LepList());
  L_elel += Leptonic(LepList().Add(elm).Add(elm), LepList());
  Leptonic    L_mumu(LepList().Add(mup), LepList().Add(mum), "mumu");
  L_mumu += Leptonic(LepList().Add(mup).Add(mum), LepList());
  L_mumu += Leptonic(LepList().Add(mup), LepList().Add(mup));
  L_mumu += Leptonic(LepList().Add(mum), LepList().Add(mum));
  L_mumu += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_mumu += Leptonic(LepList().Add(mum).Add(mum), LepList());
  Leptonic    L_elmu(LepList().Add(elp), LepList().Add(mum), "elmu");
  L_elmu += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_elmu += Leptonic(LepList().Add(elp).Add(mum), LepList());
  L_elmu += Leptonic(LepList().Add(elm).Add(mup), LepList());
  L_elmu += Leptonic(LepList().Add(elp), LepList().Add(mup));
  L_elmu += Leptonic(LepList().Add(elm), LepList().Add(mum));
  L_elmu += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_elmu += Leptonic(LepList().Add(elm).Add(mum), LepList());
  CategoryList Cats_2Lflavor;
  Cats_2Lflavor += Category(L_elel, "Ch2L");
  Cats_2Lflavor += Category(L_mumu, "Ch2L");
  Cats_2Lflavor += Category(L_elmu, "Ch2L");

  // 2L, fully inclusive
  Leptonic    L_2Lincl(LepList().Add(elp), LepList().Add(elm), "ll");
  L_2Lincl += Leptonic(LepList().Add(elp).Add(elm), LepList());
  L_2Lincl += Leptonic(LepList().Add(elp), LepList().Add(elp));
  L_2Lincl += Leptonic(LepList().Add(elm), LepList().Add(elm));
  L_2Lincl += Leptonic(LepList().Add(elp).Add(elp), LepList());
  L_2Lincl += Leptonic(LepList().Add(elm).Add(elm), LepList());
  L_2Lincl += Leptonic(LepList().Add(mup), LepList().Add(mum));
  L_2Lincl += Leptonic(LepList().Add(mup).Add(mum), LepList());
  L_2Lincl += Leptonic(LepList().Add(mup), LepList().Add(mup));
  L_2Lincl += Leptonic(LepList().Add(mum), LepList().Add(mum));
  L_2Lincl += Leptonic(LepList().Add(mup).Add(mup), LepList());
  L_2Lincl += Leptonic(LepList().Add(mum).Add(mum), LepList());
  L_2Lincl += Leptonic(LepList().Add(elp), LepList().Add(mum));
  L_2Lincl += Leptonic(LepList().Add(elm), LepList().Add(mup));
  L_2Lincl += Leptonic(LepList().Add(elp).Add(mum), LepList());
  L_2Lincl += Leptonic(LepList().Add(elm).Add(mup), LepList());
  L_2Lincl += Leptonic(LepList().Add(elp), LepList().Add(mup));
  L_2Lincl += Leptonic(LepList().Add(elm), LepList().Add(mum));
  L_2Lincl += Leptonic(LepList().Add(elp).Add(mup), LepList());
  L_2Lincl += Leptonic(LepList().Add(elm).Add(mum), LepList());
  CategoryList Cats_2Lincl;
  Cats_2Lincl += Category(L_2Lincl, "Ch2L");


  // Get hadronic ISR region lists - inclusive and b-separated
  vector<Hadronic> H_ISR_noB;
  vector<Hadronic> H_ISR_B;
  H_ISR_noB.push_back(GetHadronicRegion(1, 0)); // ge1j
  H_ISR_B.push_back(GetHadronicRegion(1, 6)); // ge1j0b
  H_ISR_B.push_back(GetHadronicRegion(1, 7)); // ge1jge1b
  
 // Get PTISR binning
  vector<double> PTISR_incl;
  PTISR_incl.push_back(250.);

  vector<double> PTISR_bin;
  PTISR_bin.push_back(250.);
  PTISR_bin.push_back(350.);

 
  // Get SV binning
  vector<double> SV_eta;
  SV_eta.push_back(0.);
  SV_eta.push_back(1.5);

  // Get gammaT binning
  vector<double> gamT0;
  gamT0.push_back(0.);
  vector<double> gamT;
  gamT.push_back(gamTval);
  gamT.push_back(0.5);

  // Define Lepton ID-based categories:
  vector<LepIDsList> IDs_all;
  vector<LepIDsList> IDs_gold;
  vector<LepIDsList> IDs_slvrbron;
  vector<LepIDsList> IDs_slvr;
  vector<LepIDsList> IDs_bron;  

  LepIDsList ID_gold("gold");
  LepIDsList ID_silver("slvr");
  LepIDsList ID_bronze("bron");

  ID_gold += LepIDs().a(kGold).a(kGold);
  ID_silver += LepIDs().a(kGold).a(kSilver);
  ID_silver += LepIDs().a(kSilver).a(kSilver);
  ID_bronze += LepIDs().a(kGold).a(kBronze);
  ID_bronze += LepIDs().a(kSilver).a(kBronze);
  // ID_bronze += LepIDs().a(kBronze).a(kBronze); // do we want/need this?

  IDs_all.push_back(ID_gold);
  IDs_all.push_back(ID_silver);
  IDs_all.push_back(ID_bronze);
  IDs_slvrbron.push_back(ID_silver);
  IDs_slvrbron.push_back(ID_bronze);
  IDs_gold.push_back(ID_gold);
  IDs_slvr.push_back(ID_silver);
  IDs_bron.push_back(ID_bronze);
  
  ///////////////////////////////////////////////
  
  cout << "Building 2L 0S object regions" << endl;
  CategoryList Cats_2L_0S_gold_OS = Cats_2Lcharge_OS;
  CategoryList Cats_2L_0S_gold_SS = Cats_2Lnoflavor_SS;
  //CategoryList Cats_2L_0S_slvrbron = Cats_2Lflavor;
  CategoryList Cats_2L_0S_slvr = Cats_2Lflavor;
  CategoryList Cats_2L_0S_bron = Cats_2Lflavor;
 
  vector<Hadronic> H_2L_0S;
  H_2L_0S.push_back(GetHadronicRegion(0, 1)); // 0j0sv
  

  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateFitBinRegions(GetFitBin(2, 0, maskSR)); 
  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateHadronicSRegions(H_2L_0S);
  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateHadronicISRRegions(H_ISR_B);
  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateGenericRegions("PTISR", PTISR_bin, maskSR);
  Cats_2L_0S_gold_OS = Cats_2L_0S_gold_OS.CreateGenericRegions("gamT", gamT);

  Cats_2L_0S_gold_SS = Cats_2L_0S_gold_SS.CreateFitBinRegions(GetFitBin(2, 0, maskSR));
  Cats_2L_0S_gold_SS = Cats_2L_0S_gold_SS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_0S_gold_SS = Cats_2L_0S_gold_SS.CreateHadronicSRegions(H_2L_0S);
  Cats_2L_0S_gold_SS = Cats_2L_0S_gold_SS.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0S_gold_SS = Cats_2L_0S_gold_SS.CreateGenericRegions("PTISR", PTISR_incl);
  
//spliting slvrbron to have different binning for CR
//  Cats_2L_0S_slvrbron = Cats_2L_0S_slvrbron.CreateFitBinRegions(GetFitBin(2, 0, maskSR));
//  Cats_2L_0S_slvrbron = Cats_2L_0S_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
//  Cats_2L_0S_slvrbron = Cats_2L_0S_slvrbron.CreateHadronicSRegions(H_2L_0S);
//  Cats_2L_0S_slvrbron = Cats_2L_0S_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
//  Cats_2L_0S_slvrbron = Cats_2L_0S_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L_0S_slvr = Cats_2L_0S_slvr.CreateFitBinRegions(GetFitBin(2, 0, maskSR));
  Cats_2L_0S_slvr = Cats_2L_0S_slvr.CreateLeptonIDRegions(IDs_slvr);
  Cats_2L_0S_slvr = Cats_2L_0S_slvr.CreateHadronicSRegions(H_2L_0S);
  Cats_2L_0S_slvr = Cats_2L_0S_slvr.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0S_slvr = Cats_2L_0S_slvr.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L_0S_bron = Cats_2L_0S_bron.CreateFitBinRegions(GetFitBin(2, 0, maskSR));//bin uniformity change
  Cats_2L_0S_bron = Cats_2L_0S_bron.CreateLeptonIDRegions(IDs_bron);
  Cats_2L_0S_bron = Cats_2L_0S_bron.CreateHadronicSRegions(H_2L_0S);
  Cats_2L_0S_bron = Cats_2L_0S_bron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0S_bron = Cats_2L_0S_bron.CreateGenericRegions("PTISR", PTISR_incl);




  Cats_2L += Cats_2L_0S_gold_OS;
  Cats_2L += Cats_2L_0S_gold_SS;
  //Cats_2L += Cats_2L_0S_slvrbron;
  Cats_2L += Cats_2L_0S_slvr;
  Cats_2L += Cats_2L_0S_bron;

  ///////////////////////////////////////////////
  
  cout << "Building 2L 0 jet 1SV object regions" << endl;
  CategoryList Cats_2L_0j1sv_gold = Cats_2Lincl;
 // CategoryList Cats_2L_0j1sv_slvrbron = Cats_2Lincl;
  CategoryList Cats_2L_0j1sv_slvr = Cats_2Lincl;
  CategoryList Cats_2L_0j1sv_bron = Cats_2Lincl;

  vector<Hadronic> H_2L_0j1sv;
  H_2L_0j1sv.push_back(GetHadronicRegion(0, 2)); // 0jge1sv

  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateFitBinRegions(GetFitBin(2, 0, maskSR)); 
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateHadronicSRegions(H_2L_0j1sv);
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateGenericRegions("PTISR", PTISR_incl);
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateGenericRegions("gamT", gamT0);
  Cats_2L_0j1sv_gold = Cats_2L_0j1sv_gold.CreateGenericRegions("SVeta", SV_eta);

//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateFitBinRegions(GetFitBin(2, 0, maskSR));
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateHadronicSRegions(H_2L_0j1sv);
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateGenericRegions("gamT", gamT0);
//  Cats_2L_0j1sv_slvrbron = Cats_2L_0j1sv_slvrbron.CreateGenericRegions("SVeta", SV_eta);

  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateFitBinRegions(GetFitBin(2, 0, maskSR));
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateLeptonIDRegions(IDs_slvr);
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateHadronicSRegions(H_2L_0j1sv);
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateGenericRegions("PTISR", PTISR_incl);
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateGenericRegions("gamT", gamT0);
  Cats_2L_0j1sv_slvr = Cats_2L_0j1sv_slvr.CreateGenericRegions("SVeta", SV_eta);

  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateFitBinRegions(GetFitBin(2, 0, maskSR));//bin uniformity change
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateLeptonIDRegions(IDs_bron);
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateHadronicSRegions(H_2L_0j1sv);
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateGenericRegions("PTISR", PTISR_incl);
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateGenericRegions("gamT", gamT0);
  Cats_2L_0j1sv_bron = Cats_2L_0j1sv_bron.CreateGenericRegions("SVeta", SV_eta);




  Cats_2L += Cats_2L_0j1sv_gold;
//  Cats_2L += Cats_2L_0j1sv_slvrbron;
  Cats_2L += Cats_2L_0j1sv_slvr;
  Cats_2L += Cats_2L_0j1sv_bron;

  ///////////////////////////////////////////////

  cout << "Building 2L 1j 0/1b object regions" << endl;
  CategoryList Cats_2L_1jS_gold_OS = Cats_2Lnoflavor_OS;
  CategoryList Cats_2L_1jS_gold_SS = Cats_2Lnoflavor_SS;
  CategoryList Cats_2L_1jS_slvrbron = Cats_2Lflavor;
 
  vector<Hadronic> H_2L_1jS;
  H_2L_1jS.push_back(GetHadronicRegion(1, 8)); // 1j0b
  H_2L_1jS.push_back(GetHadronicRegion(1, 9)); // 1j1b
  vector<Hadronic> H_2L_1jS_bkg;
  H_2L_1jS_bkg.push_back(GetHadronicRegion(1, 1)); // 1j

  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateFitBinRegions(GetFitBin(2, 1, maskSR));
  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateHadronicSRegions(H_2L_1jS);
  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateHadronicISRRegions(H_ISR_B);
  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateGenericRegions("PTISR", PTISR_bin, maskSR);
  Cats_2L_1jS_gold_OS = Cats_2L_1jS_gold_OS.CreateGenericRegions("gamT", gamT);

  Cats_2L_1jS_gold_SS = Cats_2L_1jS_gold_SS.CreateFitBinRegions(GetFitBin(2, 1, maskSR));
  Cats_2L_1jS_gold_SS = Cats_2L_1jS_gold_SS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_1jS_gold_SS = Cats_2L_1jS_gold_SS.CreateHadronicSRegions(H_2L_1jS_bkg);
  Cats_2L_1jS_gold_SS = Cats_2L_1jS_gold_SS.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_1jS_gold_SS = Cats_2L_1jS_gold_SS.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L_1jS_slvrbron = Cats_2L_1jS_slvrbron.CreateFitBinRegions(GetFitBin(2, 1, maskSR));//special CR binning//removed for uniformity
  Cats_2L_1jS_slvrbron = Cats_2L_1jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_2L_1jS_slvrbron = Cats_2L_1jS_slvrbron.CreateHadronicSRegions(H_2L_1jS_bkg);
  Cats_2L_1jS_slvrbron = Cats_2L_1jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_1jS_slvrbron = Cats_2L_1jS_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L += Cats_2L_1jS_gold_OS;
  Cats_2L += Cats_2L_1jS_gold_SS;
  Cats_2L += Cats_2L_1jS_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 2L ge2j 0/ge1b object regions" << endl;
  CategoryList Cats_2L_2jS_gold_OS = Cats_2Lnoflavor_OS;
  CategoryList Cats_2L_2jS_gold_SS = Cats_2Lnoflavor_SS;
  CategoryList Cats_2L_2jS_slvrbron = Cats_2Lflavor;
 
  vector<Hadronic> H_2L_2jS;
   H_2L_2jS.push_back(GetHadronicRegion(2, 6)); // ge2j0b
  // H_2L_2jS.push_back(GetHadronicRegion(2, 8)); // ge2j1b
  // H_2L_2jS.push_back(GetHadronicRegion(2, 9)); // ge2jge2b
   H_2L_2jS.push_back(GetHadronicRegion(2, 7)); // ge2jge1b
 
 // H_2L_2jS.push_back(GetHadronicRegion(2, 2)); // 2j0b 
 // H_2L_2jS.push_back(GetHadronicRegion(2, 3)); // 2jge1b
  
  vector<Hadronic> H_2L_2jS_bkg;
  // H_2L_2jS_bkg.push_back(GetHadronicRegion(2, 6)); // ge2j0b
  // H_2L_2jS_bkg.push_back(GetHadronicRegion(2, 7)); // ge2jge1b
  //H_2L_2jS_bkg.push_back(GetHadronicRegion(2, 1)); // 2j
  H_2L_2jS_bkg.push_back(GetHadronicRegion(2, 0));//ge2j inclusive

  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateFitBinRegions(GetFitBin(2, 2, maskSR));
  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateHadronicSRegions(H_2L_2jS);
  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateHadronicISRRegions(H_ISR_B);
  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateGenericRegions("PTISR", PTISR_bin, maskSR);
  Cats_2L_2jS_gold_OS = Cats_2L_2jS_gold_OS.CreateGenericRegions("gamT", gamT);

  Cats_2L_2jS_gold_SS = Cats_2L_2jS_gold_SS.CreateFitBinRegions(GetFitBin(2, 2, maskSR));
  Cats_2L_2jS_gold_SS = Cats_2L_2jS_gold_SS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_2jS_gold_SS = Cats_2L_2jS_gold_SS.CreateHadronicSRegions(H_2L_2jS_bkg);
  Cats_2L_2jS_gold_SS = Cats_2L_2jS_gold_SS.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_2jS_gold_SS = Cats_2L_2jS_gold_SS.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L_2jS_slvrbron = Cats_2L_2jS_slvrbron.CreateFitBinRegions(GetFitBin(2, 2, maskSR));//special CR binning// removed for uniformity
  Cats_2L_2jS_slvrbron = Cats_2L_2jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_2L_2jS_slvrbron = Cats_2L_2jS_slvrbron.CreateHadronicSRegions(H_2L_2jS_bkg);
  Cats_2L_2jS_slvrbron = Cats_2L_2jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_2jS_slvrbron = Cats_2L_2jS_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L += Cats_2L_2jS_gold_OS;
  Cats_2L += Cats_2L_2jS_gold_SS;
  Cats_2L += Cats_2L_2jS_slvrbron;
 
  ///////////////////////////////////////////////
/* Merging 2L ge3J into 2J 
  cout << "Building 2L ge3j 0/ge1b object regions" << endl;
  CategoryList Cats_2L_3jS_gold_OS = Cats_2Lnoflavor_OS;
  CategoryList Cats_2L_3jS_gold_SS = Cats_2Lnoflavor_SS;
  CategoryList Cats_2L_3jS_slvrbron = Cats_2Lflavor;
 
  vector<Hadronic> H_2L_3jS;
  H_2L_3jS.push_back(GetHadronicRegion(3, 6)); // ge3j0b
  H_2L_3jS.push_back(GetHadronicRegion(3, 7)); // ge3jge1b

  vector<Hadronic> H_2L_3jS_bkg;
  H_2L_3jS_bkg.push_back(GetHadronicRegion(3, 0)); // ge3j

  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateFitBinRegions(GetFitBin(2, 3, maskSR));
  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateHadronicSRegions(H_2L_3jS);
  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateHadronicISRRegions(H_ISR_B);
  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateGenericRegions("PTISR", PTISR_bin);
  Cats_2L_3jS_gold_OS = Cats_2L_3jS_gold_OS.CreateGenericRegions("gamT", gamT);

  Cats_2L_3jS_gold_SS = Cats_2L_3jS_gold_SS.CreateFitBinRegions(GetFitBin(2, 3, maskSR));
  Cats_2L_3jS_gold_SS = Cats_2L_3jS_gold_SS.CreateLeptonIDRegions(IDs_gold);
  Cats_2L_3jS_gold_SS = Cats_2L_3jS_gold_SS.CreateHadronicSRegions(H_2L_3jS_bkg);
  Cats_2L_3jS_gold_SS = Cats_2L_3jS_gold_SS.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_3jS_gold_SS = Cats_2L_3jS_gold_SS.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L_3jS_slvrbron = Cats_2L_3jS_slvrbron.CreateFitBinRegions(GetFitBin(2, 3, maskSR));
  Cats_2L_3jS_slvrbron = Cats_2L_3jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_2L_3jS_slvrbron = Cats_2L_3jS_slvrbron.CreateHadronicSRegions(H_2L_3jS_bkg);
  Cats_2L_3jS_slvrbron = Cats_2L_3jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_2L_3jS_slvrbron = Cats_2L_3jS_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_2L += Cats_2L_3jS_gold_OS;
  Cats_2L += Cats_2L_3jS_gold_SS;
  Cats_2L += Cats_2L_3jS_slvrbron;
*/
  ///////////////////////////////////////////////
 
  return Cats_2L;
}

CategoryList CategoryTool::GetCategories_3L(bool maskSR) const {
  CategoryList Cats_3L;
  cout << "Creating 3L Categories" << endl;

  // define leptons
  Lep elp(kElec, kPos);
  Lep elm(kElec, kNeg);
  Lep mup(kMuon, kPos);
  Lep mum(kMuon, kNeg);

  ////////////////////
  // flavor inclusive 
  ////////////////////
  Leptonic    L_Zstar(LepList().Add(elp).Add(elm), LepList().Add(elp), "Zstar");
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(elm));
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm).Add(elp), LepList());
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm).Add(elm), LepList());
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mup));
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mum));
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum).Add(mup), LepList());
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum).Add(mum), LepList());
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mup));
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mum));
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm).Add(mup), LepList());
  L_Zstar += Leptonic(LepList().Add(elp).Add(elm).Add(mum), LepList());
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elp), "OSSFmumu^el");
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elm));
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum).Add(elp), LepList());
  L_Zstar += Leptonic(LepList().Add(mup).Add(mum).Add(elm), LepList());

  Leptonic    L_noZ(LepList().Add(elp).Add(elp), LepList().Add(elm), "noZ");
  L_noZ += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elp));
  L_noZ += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mum));
  L_noZ += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mup));
  L_noZ += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mum));
  L_noZ += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mup));
  L_noZ += Leptonic(LepList().Add(elp).Add(elp).Add(mum), LepList());
  L_noZ += Leptonic(LepList().Add(elm).Add(elm).Add(mup), LepList());
  L_noZ += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elm));
  L_noZ += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elm));
  L_noZ += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elp));
  L_noZ += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elp));
  L_noZ += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elp));
  L_noZ += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elm));
  L_noZ += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elm));
  L_noZ += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elp));
  L_noZ += Leptonic(LepList().Add(mup).Add(mup).Add(elm), LepList());
  L_noZ += Leptonic(LepList().Add(mum).Add(mum).Add(elp), LepList());
  L_noZ += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mum));
  L_noZ += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mum));
  L_noZ += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mup));
  L_noZ += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mup));
  L_noZ += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mup));
  L_noZ += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mum));

  Leptonic    L_SS(LepList().Add(elp).Add(elp), LepList().Add(elp), "SS");
  L_SS += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elm));
  L_SS += Leptonic(LepList().Add(elp).Add(elp).Add(elp), LepList());
  L_SS += Leptonic(LepList().Add(elm).Add(elm).Add(elm), LepList());
  L_SS += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mup));
  L_SS += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mum));
  L_SS += Leptonic(LepList().Add(mup).Add(mup).Add(mup), LepList());
  L_SS += Leptonic(LepList().Add(mum).Add(mum).Add(mum), LepList());
  L_SS += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mup));
  L_SS += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mum));
  L_SS += Leptonic(LepList().Add(elp).Add(elp).Add(mup), LepList());
  L_SS += Leptonic(LepList().Add(elm).Add(elm).Add(mum), LepList());
  L_SS += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elp));
  L_SS += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elm));
  L_SS += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elp));
  L_SS += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elm));
  L_SS += Leptonic(LepList().Add(mup).Add(mup).Add(elp), LepList());
  L_SS += Leptonic(LepList().Add(mum).Add(mum).Add(elm), LepList());
  L_SS += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mup));
  L_SS += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mum));
  
  // 3 L, elelel - OSSF hem separated
  Leptonic    L_OSSF_elel_el(LepList().Add(elp).Add(elm), LepList().Add(elp), "OSSFelel^el");
  L_OSSF_elel_el += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(elm));
  L_OSSF_elel_el += Leptonic(LepList().Add(elp).Add(elm).Add(elp), LepList());
  L_OSSF_elel_el += Leptonic(LepList().Add(elp).Add(elm).Add(elm), LepList());
  Leptonic    L_nOSSF_elel_el(LepList().Add(elp).Add(elp), LepList().Add(elp), "nOSSFelel^el");
  L_nOSSF_elel_el += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elm));
  L_nOSSF_elel_el += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(elm));
  L_nOSSF_elel_el += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elp));
  L_nOSSF_elel_el += Leptonic(LepList().Add(elp).Add(elp).Add(elp), LepList());
  L_nOSSF_elel_el += Leptonic(LepList().Add(elm).Add(elm).Add(elm), LepList());
  // 3 L, elelel - OSSF incl
  Leptonic    L_OSSF_elelel(LepList().Add(elp).Add(elm), LepList().Add(elp), "OSSFelelel");
  L_OSSF_elelel += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(elm));
  L_OSSF_elelel += Leptonic(LepList().Add(elp).Add(elm).Add(elp), LepList());
  L_OSSF_elelel += Leptonic(LepList().Add(elp).Add(elm).Add(elm), LepList());
  L_OSSF_elelel += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(elm));
  L_OSSF_elelel += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elp));
  Leptonic    L_nOSSF_elelel(LepList().Add(elp).Add(elp), LepList().Add(elp), "nOSSFelelel");
  L_nOSSF_elelel += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elm));
  L_nOSSF_elelel += Leptonic(LepList().Add(elp).Add(elp).Add(elp), LepList());
  L_nOSSF_elelel += Leptonic(LepList().Add(elm).Add(elm).Add(elm), LepList());
  // 3 L, elelel - incl
  Leptonic    L_elelel(LepList().Add(elp).Add(elm), LepList().Add(elp), "elelel");
  L_elelel += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(elm));
  L_elelel += Leptonic(LepList().Add(elp).Add(elm).Add(elp), LepList());
  L_elelel += Leptonic(LepList().Add(elp).Add(elm).Add(elm), LepList());
  L_elelel += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(elm));
  L_elelel += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elp));
  L_elelel += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(elp));
  L_elelel += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elm));
  L_elelel += Leptonic(LepList().Add(elp).Add(elp).Add(elp), LepList());
  L_elelel += Leptonic(LepList().Add(elm).Add(elm).Add(elm), LepList());

  // 3 L, mumumu - OSSF hem separated
  Leptonic    L_OSSF_mumu_mu(LepList().Add(mup).Add(mum), LepList().Add(mup), "OSSFmumu^mu");
  L_OSSF_mumu_mu += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mum));
  L_OSSF_mumu_mu += Leptonic(LepList().Add(mup).Add(mum).Add(mup), LepList());
  L_OSSF_mumu_mu += Leptonic(LepList().Add(mup).Add(mum).Add(mum), LepList());
  Leptonic    L_nOSSF_mumu_mu(LepList().Add(mup).Add(mup), LepList().Add(mup), "nOSSFmumu^mu");
  L_nOSSF_mumu_mu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mum));
  L_nOSSF_mumu_mu += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mum));
  L_nOSSF_mumu_mu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mup));
  L_nOSSF_mumu_mu += Leptonic(LepList().Add(mup).Add(mup).Add(mup), LepList());
  L_nOSSF_mumu_mu += Leptonic(LepList().Add(mum).Add(mum).Add(mum), LepList());
  // 3 L, mumumu - OSSF incl
  Leptonic    L_OSSF_mumumu(LepList().Add(mup).Add(mum), LepList().Add(mup), "OSSFmumumu");
  L_OSSF_mumumu += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mum));
  L_OSSF_mumumu += Leptonic(LepList().Add(mup).Add(mum).Add(mup), LepList());
  L_OSSF_mumumu += Leptonic(LepList().Add(mup).Add(mum).Add(mum), LepList());
  L_OSSF_mumumu += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mum));
  L_OSSF_mumumu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mup));
  Leptonic    L_nOSSF_mumumu(LepList().Add(mup).Add(mup), LepList().Add(mup), "nOSSFmumumu");
  L_nOSSF_mumumu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mum));
  L_nOSSF_mumumu += Leptonic(LepList().Add(mup).Add(mup).Add(mup), LepList());
  L_nOSSF_mumumu += Leptonic(LepList().Add(mum).Add(mum).Add(mum), LepList());
  // 3 L, mumumu - incl
  Leptonic    L_mumumu(LepList().Add(mup).Add(mum), LepList().Add(mup), "mumumu");
  L_mumumu += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mum));
  L_mumumu += Leptonic(LepList().Add(mup).Add(mum).Add(mup), LepList());
  L_mumumu += Leptonic(LepList().Add(mup).Add(mum).Add(mum), LepList());
  L_mumumu += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mum));
  L_mumumu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mup));
  L_mumumu += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mup));
  L_mumumu += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mum));
  L_mumumu += Leptonic(LepList().Add(mup).Add(mup).Add(mup), LepList());
  L_mumumu += Leptonic(LepList().Add(mum).Add(mum).Add(mum), LepList());

  // 3 L, elelmu - OSSF hem separated
  Leptonic    L_OSSF_elel_mu(LepList().Add(elp).Add(elm), LepList().Add(mup), "OSSFelel^mu");
  L_OSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mum));
  L_OSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elm).Add(mup), LepList());
  L_OSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elm).Add(mum), LepList());
  Leptonic    L_nOSSF_elel_mu(LepList().Add(elp).Add(elp), LepList().Add(mup), "nOSSFelel^mu");
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mum));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mum));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mup));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elp).Add(mup), LepList());
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(elm).Add(mum), LepList());
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(elp).Add(mum), LepList());
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(elm).Add(mup), LepList());
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elm));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elm));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elp));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elp));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elp));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elp));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elm));
  L_nOSSF_elel_mu += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elm));
  // 3 L, elelmu - OSSF incl
  // Leptonic    L_OSSF_elelmu(LepList().Add(elp).Add(elm), LepList().Add(mup), "OSSFelelmu");
  // L_OSSF_elelmu += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mum));
  // L_OSSF_elelmu += Leptonic(LepList().Add(elp).Add(elm).Add(mup), LepList());
  // L_OSSF_elelmu += Leptonic(LepList().Add(elp).Add(elm).Add(mum), LepList());
  // L_OSSF_elelmu += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elm));
  // L_OSSF_elelmu += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elm));
  // L_OSSF_elelmu += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elp));
  // L_OSSF_elelmu += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elp));
  // Leptonic    L_nOSSF_elelmu(LepList().Add(elp).Add(elp), LepList().Add(mup), "nOSSFelelmu");
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mum));
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mum));
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mup));
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elp).Add(elp).Add(mup), LepList());
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elm).Add(elm).Add(mum), LepList());
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elp).Add(elp).Add(mum), LepList());
  // L_nOSSF_elelmu += Leptonic(LepList().Add(elm).Add(elm).Add(mup), LepList());
  // 3 L, elelmu - incl
  Leptonic    L_elelmu(LepList().Add(elp).Add(elm), LepList().Add(mup), "elelmu");
  L_elelmu += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mum));
  L_elelmu += Leptonic(LepList().Add(elp).Add(elm).Add(mup), LepList());
  L_elelmu += Leptonic(LepList().Add(elp).Add(elm).Add(mum), LepList());
  L_elelmu += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elm));
  L_elelmu += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elm));
  L_elelmu += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elp));
  L_elelmu += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elp));
  L_elelmu += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elp));
  L_elelmu += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elp));
  L_elelmu += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elm));
  L_elelmu += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elm));
  L_elelmu += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mup));
  L_elelmu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mum));
  L_elelmu += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mum));
  L_elelmu += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mup));
  L_elelmu += Leptonic(LepList().Add(elp).Add(elp).Add(mup), LepList());
  L_elelmu += Leptonic(LepList().Add(elm).Add(elm).Add(mum), LepList());
  L_elelmu += Leptonic(LepList().Add(elp).Add(elp).Add(mum), LepList());
  L_elelmu += Leptonic(LepList().Add(elm).Add(elm).Add(mup), LepList());

    // 3 L, mumuel - OSSF hem separated
  Leptonic    L_OSSF_mumu_el(LepList().Add(mup).Add(mum), LepList().Add(elp), "OSSFmumu^el");
  L_OSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elm));
  L_OSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mum).Add(elp), LepList());
  L_OSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mum).Add(elm), LepList());
  Leptonic    L_nOSSF_mumu_el(LepList().Add(mup).Add(mup), LepList().Add(elp), "nOSSFmumu^el");
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elm));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elm));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elp));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mup).Add(elp), LepList());
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(mum).Add(elm), LepList());
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(mup).Add(elm), LepList());
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(mum).Add(elp), LepList());
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mum));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mum));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mup));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mup));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mup));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mup));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mum));
  L_nOSSF_mumu_el += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mum));
  // // 3 L, mumuel - OSSF incl
  // Leptonic    L_OSSF_mumuel(LepList().Add(mup).Add(mum), LepList().Add(elp), "OSSFmumuel");
  // L_OSSF_mumuel += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elm));
  // L_OSSF_mumuel += Leptonic(LepList().Add(mup).Add(mum).Add(elp), LepList());
  // L_OSSF_mumuel += Leptonic(LepList().Add(mup).Add(mum).Add(elm), LepList());
  // L_OSSF_mumuel += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mum));
  // L_OSSF_mumuel += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mum));
  // L_OSSF_mumuel += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mup));
  // L_OSSF_mumuel += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mup));
  // Leptonic    L_nOSSF_mumuel(LepList().Add(mup).Add(mup), LepList().Add(elp), "nOSSFmumuel");
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elm));
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elm));
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elp));
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mup).Add(mup).Add(elp), LepList());
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mum).Add(mum).Add(elm), LepList());
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mup).Add(mup).Add(elm), LepList());
  // L_nOSSF_mumuel += Leptonic(LepList().Add(mum).Add(mum).Add(elp), LepList());
  // 3 L, mumuel - incl
  Leptonic    L_mumuel(LepList().Add(mup).Add(mum), LepList().Add(elp), "mumuel");
  L_mumuel += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elm));
  L_mumuel += Leptonic(LepList().Add(mup).Add(mum).Add(elp), LepList());
  L_mumuel += Leptonic(LepList().Add(mup).Add(mum).Add(elm), LepList());
  L_mumuel += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mum));
  L_mumuel += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mum));
  L_mumuel += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mup));
  L_mumuel += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mup));
  L_mumuel += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mup));
  L_mumuel += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mup));
  L_mumuel += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mum));
  L_mumuel += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mum));
  L_mumuel += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elp));
  L_mumuel += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elm));
  L_mumuel += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elm));
  L_mumuel += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elp));
  L_mumuel += Leptonic(LepList().Add(mup).Add(mup).Add(elp), LepList());
  L_mumuel += Leptonic(LepList().Add(mum).Add(mum).Add(elm), LepList());
  L_mumuel += Leptonic(LepList().Add(mup).Add(mup).Add(elm), LepList());
  L_mumuel += Leptonic(LepList().Add(mum).Add(mum).Add(elp), LepList());

  /////////////////////////////
  // Leptonic    L_OSSF_elel_lep(LepList().Add(elp).Add(elm), LepList().Add(elp), "OSSFelel^lep");
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(elm));
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm).Add(elp), LepList());
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm).Add(elm), LepList());
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mup));
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm), LepList().Add(mum));
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm).Add(mup), LepList());
  // L_OSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elm).Add(mum), LepList());
  // Leptonic    L_nOSSF_elel_lep(LepList().Add(elp).Add(elp), LepList().Add(elp), "nOSSFelel^lep");
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elm));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(elm));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(elp));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp).Add(elp), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm).Add(elm), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mup));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mum));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp), LepList().Add(mum));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm), LepList().Add(mup));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp).Add(mup), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm).Add(mum), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(elp).Add(mum), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(elm).Add(mup), LepList());
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(mup), LepList().Add(elm));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elp).Add(mum), LepList().Add(elm));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(mup), LepList().Add(elp));
  // L_nOSSF_elel_lep += Leptonic(LepList().Add(elm).Add(mum), LepList().Add(elp));

  // Leptonic    L_OSSF_mumu_lep(LepList().Add(mup).Add(mum), LepList().Add(mup), "OSSFmumu^lep");
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(mum));
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum).Add(mup), LepList());
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum).Add(mum), LepList());
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elp));
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum), LepList().Add(elm));
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum).Add(elp), LepList());
  // L_OSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mum).Add(elm), LepList());
  // Leptonic    L_nOSSF_mumu_lep(LepList().Add(mup).Add(mup), LepList().Add(mup), "nOSSFmumu^lep");
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mum));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(mum));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(mup));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup).Add(mup), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum).Add(mum), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elp));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elm));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup), LepList().Add(elm));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum), LepList().Add(elp));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup).Add(elp), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum).Add(elm), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(mup).Add(elm), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(mum).Add(elp), LepList());
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(elp), LepList().Add(mum));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mup).Add(elm), LepList().Add(mum));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(elp), LepList().Add(mup));
  // L_nOSSF_mumu_lep += Leptonic(LepList().Add(mum).Add(elm), LepList().Add(mup));

  /////////////////////////////

  CategoryList Cats_3Lnoflavor;
  Cats_3Lnoflavor += Category(L_Zstar, "Ch3L");
  Cats_3Lnoflavor += Category(L_noZ, "Ch3L");
  Cats_3Lnoflavor += Category(L_SS, "Ch3L"); 

  CategoryList Cats_3Lnoflavor_OS;
  Cats_3Lnoflavor_OS += Category(L_Zstar, "Ch3L");
  Cats_3Lnoflavor_OS += Category(L_noZ, "Ch3L");

  CategoryList Cats_3Lnoflavor_SS;
  Cats_3Lnoflavor_SS += Category(L_SS, "Ch3L");
  
  CategoryList Cats_3Lchargehem;
  Cats_3Lchargehem += Category(L_OSSF_elel_el, "Ch3L");
  Cats_3Lchargehem += Category(L_nOSSF_elel_el, "Ch3L");
  Cats_3Lchargehem += Category(L_OSSF_mumu_mu, "Ch3L");
  Cats_3Lchargehem += Category(L_nOSSF_mumu_mu, "Ch3L");
  Cats_3Lchargehem += Category(L_OSSF_elel_mu, "Ch3L");
  Cats_3Lchargehem += Category(L_nOSSF_elel_mu, "Ch3L");
  Cats_3Lchargehem += Category(L_OSSF_mumu_el, "Ch3L");
  Cats_3Lchargehem += Category(L_nOSSF_mumu_el, "Ch3L");

  CategoryList Cats_3Lchargehem_OSSF;
  Cats_3Lchargehem_OSSF += Category(L_OSSF_elel_el, "Ch3L");
  Cats_3Lchargehem_OSSF += Category(L_OSSF_elel_mu, "Ch3L");
  Cats_3Lchargehem_OSSF += Category(L_OSSF_mumu_mu, "Ch3L");
  Cats_3Lchargehem_OSSF += Category(L_OSSF_mumu_el, "Ch3L");
  
  CategoryList Cats_3Lchargehem_nOSSF;
  Cats_3Lchargehem_nOSSF += Category(L_nOSSF_elel_el, "Ch3L");
  Cats_3Lchargehem_nOSSF += Category(L_nOSSF_mumu_mu, "Ch3L");
  Cats_3Lchargehem_nOSSF += Category(L_nOSSF_elel_mu, "Ch3L");
  Cats_3Lchargehem_nOSSF += Category(L_nOSSF_mumu_el, "Ch3L");

  // CategoryList Cats_3Lcharge;
  // Cats_3Lcharge += Category(L_OSSF_elelel, "Ch3L");
  // Cats_3Lcharge += Category(L_nOSSF_elelel, "Ch3L");
  // Cats_3Lcharge += Category(L_OSSF_mumumu, "Ch3L");
  // Cats_3Lcharge += Category(L_nOSSF_mumumu, "Ch3L");
  // Cats_3Lcharge += Category(L_OSSF_elelmu, "Ch3L");
  // Cats_3Lcharge += Category(L_nOSSF_elelmu, "Ch3L");
  // Cats_3Lcharge += Category(L_OSSF_mumuel, "Ch3L");
  // Cats_3Lcharge += Category(L_nOSSF_mumuel, "Ch3L");

  CategoryList Cats_3Lflavor;
  Cats_3Lflavor += Category(L_elelel, "Ch3L");
  Cats_3Lflavor += Category(L_mumumu, "Ch3L");
  Cats_3Lflavor += Category(L_elelmu, "Ch3L");
  Cats_3Lflavor += Category(L_mumuel, "Ch3L");
 
  // Get hadronic ISR region lists - inclusive and b-separated
  vector<Hadronic> H_ISR_noB;
  vector<Hadronic> H_ISR_B;
  H_ISR_noB.push_back(GetHadronicRegion(1, 0)); // ge1j
  H_ISR_B.push_back(GetHadronicRegion(1, 6)); // ge1j0b
  H_ISR_B.push_back(GetHadronicRegion(1, 7)); // ge1jge1b
  
 // Get PTISR binning
  vector<double> PTISR_incl;
  PTISR_incl.push_back(250.); // minimum PTISR_incl

  vector<double> PTISR_bin;
  PTISR_bin.push_back(200.);
  PTISR_bin.push_back(300.);
 
  // Get SV binning
  vector<double> SV_eta;
  SV_eta.push_back(0.);
  SV_eta.push_back(1.5);

  // Get gammaT binning
  vector<double> gamT0;
  gamT0.push_back(0.);
  vector<double> gamT;
  gamT.push_back(gamTval);
  gamT.push_back(0.5);

  // Define Lepton ID-based categories:
  vector<LepIDsList> IDs_all;
  vector<LepIDsList> IDs_gold;
  vector<LepIDsList> IDs_slvrbron;
  
  LepIDsList ID_gold("gold");
  LepIDsList ID_silver("slvr");
  LepIDsList ID_bronze("bron");

  ID_gold += LepIDs().a(kGold).a(kGold).a(kGold);
  ID_silver += LepIDs().a(kGold).a(kGold).a(kSilver);
  ID_silver += LepIDs().a(kGold).a(kSilver).a(kSilver);
  ID_bronze += LepIDs().a(kGold).a(kGold).a(kBronze);
  ID_bronze += LepIDs().a(kGold).a(kSilver).a(kBronze);
  // ID_bronze += LepIDs().a(kGold).a(kBronze).a(kBronze); // do we want/need this?

  IDs_all.push_back(ID_gold);
  IDs_all.push_back(ID_silver);
  IDs_all.push_back(ID_bronze);
  IDs_slvrbron.push_back(ID_silver);
  IDs_slvrbron.push_back(ID_bronze);
  IDs_gold.push_back(ID_gold);

  ///////////////////////////////////////////////
  
  cout << "Building 3L 0S object regions" << endl;
  //CategoryList Cats_3L_0S_gold = Cats_3Lchargehem; //Cats_3Lcharge
  CategoryList Cats_3L_0S_gold = Cats_3Lnoflavor_OS;
  CategoryList Cats_3L_0S_slvrbron = Cats_3Lnoflavor_OS;
 
  vector<Hadronic> H_3L_0S;
  H_3L_0S.push_back(GetHadronicRegion(0, 5)); // 0j

  Cats_3L_0S_gold = Cats_3L_0S_gold.CreateFitBinRegions(GetFitBin(3, 0, maskSR)); // R1
  Cats_3L_0S_gold = Cats_3L_0S_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_3L_0S_gold = Cats_3L_0S_gold.CreateHadronicSRegions(H_3L_0S);
  Cats_3L_0S_gold = Cats_3L_0S_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_0S_gold = Cats_3L_0S_gold.CreateGenericRegions("PTISR", PTISR_incl);
  //Cats_3L_0S_gold = Cats_3L_0S_gold.CreateGenericRegions("gamT", gamT);

  Cats_3L_0S_slvrbron = Cats_3L_0S_slvrbron.CreateFitBinRegions(GetFitBin(3, 0, maskSR));
  Cats_3L_0S_slvrbron = Cats_3L_0S_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_3L_0S_slvrbron = Cats_3L_0S_slvrbron.CreateHadronicSRegions(H_3L_0S);
  Cats_3L_0S_slvrbron = Cats_3L_0S_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_0S_slvrbron = Cats_3L_0S_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_3L += Cats_3L_0S_gold;
  Cats_3L += Cats_3L_0S_slvrbron;
  
  ///////////////////////////////////////////////

  cout << "Building 3L ge1j 0/1b object regions" << endl;
  CategoryList Cats_3L_1jS_gold = Cats_3Lnoflavor_OS;
  CategoryList Cats_3L_1jS_slvrbron = Cats_3Lnoflavor_OS;
  
  vector<Hadronic> H_3L_1jS;
  H_3L_1jS.push_back(GetHadronicRegion(1, 6)); // ge1j0b
  H_3L_1jS.push_back(GetHadronicRegion(1, 7)); // ge1jge1b

  vector<Hadronic> H_3L_1jS_bkg;
  H_3L_1jS_bkg.push_back(GetHadronicRegion(1, 0)); // ge1j

  Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateFitBinRegions(GetFitBin(3, 1, maskSR));
  Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateHadronicSRegions(H_3L_1jS_bkg);
  Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateGenericRegions("PTISR", PTISR_incl);
  //Cats_3L_1jS_gold = Cats_3L_1jS_gold.CreateGenericRegions("gamT", gamT);
 
  Cats_3L_1jS_slvrbron = Cats_3L_1jS_slvrbron.CreateFitBinRegions(GetFitBin(3, 1, maskSR));
  Cats_3L_1jS_slvrbron = Cats_3L_1jS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_3L_1jS_slvrbron = Cats_3L_1jS_slvrbron.CreateHadronicSRegions(H_3L_1jS_bkg);
  Cats_3L_1jS_slvrbron = Cats_3L_1jS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_1jS_slvrbron = Cats_3L_1jS_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_3L += Cats_3L_1jS_gold;
  Cats_3L += Cats_3L_1jS_slvrbron;
 
  ///////////////////////////////////////////////

  cout << "Building 3L SSS regions" << endl;
  CategoryList Cats_3L_SSS_gold = Cats_3Lnoflavor_SS;
  CategoryList Cats_3L_SSS_slvrbron = Cats_3Lnoflavor_SS;
  
  vector<Hadronic> H_3L_SSS_bkg;
  H_3L_SSS_bkg.push_back(GetHadronicRegion(0, 0)); // ge0j (fully inclusive)

  Cats_3L_SSS_gold = Cats_3L_SSS_gold.CreateFitBinRegions(GetFitBin(3, 0, maskSR));
  Cats_3L_SSS_gold = Cats_3L_SSS_gold.CreateLeptonIDRegions(IDs_gold);
  Cats_3L_SSS_gold = Cats_3L_SSS_gold.CreateHadronicSRegions(H_3L_SSS_bkg);
  Cats_3L_SSS_gold = Cats_3L_SSS_gold.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_SSS_gold = Cats_3L_SSS_gold.CreateGenericRegions("PTISR", PTISR_incl);
 
  Cats_3L_SSS_slvrbron = Cats_3L_SSS_slvrbron.CreateFitBinRegions(GetFitBin(3, 0, maskSR));
  Cats_3L_SSS_slvrbron = Cats_3L_SSS_slvrbron.CreateLeptonIDRegions(IDs_slvrbron);
  Cats_3L_SSS_slvrbron = Cats_3L_SSS_slvrbron.CreateHadronicSRegions(H_3L_SSS_bkg);
  Cats_3L_SSS_slvrbron = Cats_3L_SSS_slvrbron.CreateHadronicISRRegions(H_ISR_noB);
  Cats_3L_SSS_slvrbron = Cats_3L_SSS_slvrbron.CreateGenericRegions("PTISR", PTISR_incl);

  Cats_3L += Cats_3L_SSS_gold;
  Cats_3L += Cats_3L_SSS_slvrbron;
 
  ///////////////////////////////////////////////

  // cout << "Building 3L ge2j 0/ge1b object regions" << endl;
  // CategoryList Cats_3L_2jS_gold = Cats_3Lchargehem;
  // CategoryList Cats_3L_2jS_silver = Cats_3Lflavor;
  // CategoryList Cats_3L_2jS_bronze = Cats_3Lflavor;
  // vector<Hadronic> H_3L_2jS;
  // H_3L_2jS.push_back(GetHadronicRegion(2, 6)); // ge2j0b
  // H_3L_2jS.push_back(GetHadronicRegion(2, 7)); // ge2jge1b

  // Cats_3L_2jS_gold = Cats_3L_2jS_gold.CreateFitBinRegions(GetFitBin(4));
  // Cats_3L_2jS_gold = Cats_3L_2jS_gold.CreateLeptonIDRegions(IDs_gold, 1);
  // Cats_3L_2jS_gold = Cats_3L_2jS_gold.CreateHadronicSRegions(H_3L_2jS);
  // Cats_3L_2jS_gold = Cats_3L_2jS_gold.CreateHadronicISRRegions(H_ISR_B);
  // Cats_3L_2jS_gold = Cats_3L_2jS_gold.CreateGenericRegions("PTISR", PTISR_incl);

  // Cats_3L_2jS_silver = Cats_3L_2jS_silver.CreateFitBinRegions(GetFitBin(4));
  // Cats_3L_2jS_silver = Cats_3L_2jS_silver.CreateLeptonIDRegions(IDs_silver, 1);
  // Cats_3L_2jS_silver = Cats_3L_2jS_silver.CreateHadronicSRegions(H_3L_2jS);
  // Cats_3L_2jS_silver = Cats_3L_2jS_silver.CreateHadronicISRRegions(H_ISR_noB);
  // Cats_3L_2jS_silver = Cats_3L_2jS_silver.CreateGenericRegions("PTISR", PTISR_incl);

  // Cats_3L_2jS_bronze = Cats_3L_2jS_bronze.CreateFitBinRegions(GetFitBin(4));
  // Cats_3L_2jS_bronze = Cats_3L_2jS_bronze.CreateLeptonIDRegions(IDs_bronze, 1);
  // Cats_3L_2jS_bronze = Cats_3L_2jS_bronze.CreateHadronicSRegions(H_3L_2jS);
  // Cats_3L_2jS_bronze = Cats_3L_2jS_bronze.CreateHadronicISRRegions(H_ISR_noB);
  // Cats_3L_2jS_bronze = Cats_3L_2jS_bronze.CreateGenericRegions("PTISR", PTISR_incl);

  // Cats_3L += Cats_3L_2jS_gold;
  // Cats_3L += Cats_3L_2jS_silver;
  // Cats_3L += Cats_3L_2jS_bronze;

  ///////////////////////////////////////////////
  
  return Cats_3L;
}
