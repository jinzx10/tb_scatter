#include "../modelconst.h"
#include "../mathconst.h"

using namespace arma;
using namespace std;

// number of trajectories
const uword cNumTraj = 20;

// if true, time and total potential are saved to txt files with an additional "_test" postfix
const bool cConvergenceTest = true;

// fixed-terminal
const double cTerminalZ = 15.001;
const uword cMaxNumStep = 6000;
const double cDt = 1.0;

// friction switch for EFLD
const bool cFricSwitch = false;

// BCME correction force switch, becomes CME if false
const bool cBroadenSwitch = false;

// freeze the total lattice
const bool cLatFreeze = false;

// restrict the movements of lattice sites to the first 2 dimensions
const bool cLatZFreeze = false;

// print intermediate information
const bool cPrintSwitch = true;

// 'g' (grand canonical) or 'c' (canonical)
const char cEnsemble = 'g'; // 'c' only works in EFLD

// thermodynamic quantity 
const double cThermBeta = 1000.0;
const double cChemPot = 0.0; // unused in canonical ensemble
const double cNumElec = 1.0; // unused in grand canonical ensemble

// broaden the energy level to estimate DOS
const double cLevelBroadenWidth = 0.05;

// initial condition of the center of the molecule
//const vec cMolInitCenterCoor = {-2.625, 7.1447, 20}; // for 45 degree incidence
const vec cMolInitCenterCoor = {13.75, 4.7631, 15}; // on-top (on site 6)
//const vec cMolInitCenterCoor = {13.75, 7.9386, 15}; // empty hollow (among site 6, 9, 10)
//const vec cMolInitCenterCoor = {11.00, 6.3509, 15}; // full hollow (among site 5, 6, 9)
//const vec cMolInitCenterCoor = {12.375, 7.1447, 15}; // bridge (between site 6, 9)
//const vec cMolInitCenterCoor = {0, 0, 18.5};
//const vec cMolInitCenterCoor = {12, 0, 0};
//const vec cMolInitCenterVelo = {0.007071, 0.0, -0.007071}; // Wodtke2015 H incident kinetic energy 2.7 eV
const vec cMolInitCenterVelo = {0.0, 0.0, -0.01}; // Wodtke2015 H incident kinetic energy 2.7 eV

namespace monoatomic
{
    const vector<vec> cMolCoor = { vec{0.0, 0.0, 0.0} };
    const vec cMolMassList = { 1836.0 };

    const double cMolPotCoef = 0;
    const double cMolPotEqDist = 0;
    const double cMolPotExpDecayLen = 0;

    const double cImageC = 0.24;
    const double cImageD = 0.063;
    const double cRefZ = 3.6; // Tully2009 z_image = 1.15 A
    const double cWorkFunc = 0.1874;
    const double cElecAffinity = 0.0277;
}

namespace hexagonal
{
    // UnitCell
    const vector<mat> cLatOnSiteE = { mat{0.05} };
    // Chang2008 gold p-orbital on-site 1.4188 eV ~ 0.05 a.u.
    const vector<vec> cLatInCellCoor = { vec{0,0,0} };
    const vec cLatMassList = {360000}; // mass of gold atom ~ 197 u ~ 197*1840 a.u.

    // BravLat
    // bond length of gold ~ 288 pm ~ 5.5 a.u.
    // [111] plane
    const vector<vec> cLatVec = {   5.5*vec{1.0, 0.0, 0.0},
				    5.5*vec{0.5, 0.5*ROOT3, 0.0},
    				    5.5*vec{0.5, 0.5/ROOT3, -ROOT2/ROOT3} };
    // [110] plane
//    const vector<vec> cLatVec = {   5.5*vec{1.0, 0.0, 0.0},
//				    5.5*vec{0.0, ROOT2, 0.0},
//				    5.5*vec{0.5, 1.0/ROOT2, -0.5}  };
    
    // LatSupCell
    // [111]
    const uvec cNumUnitCell = {4, 4, 3};
    // [110]
    // const uvec cNumUnitCell = {4, 4, 2};
    const uvec cIsPeriodic = {1, 1, 0};
    
    // BrilZone
    const uvec cNumK = {1, 1};

    // LatHop
    const arma::cx_mat cLatHopCoef = cx_mat{ mat{-0.05}, zeros(1,1) };
    // Chang2008 gold hopping amplitude -0.3056 eV ~ -0.01 a.u.
    const arma::mat cLatHopDecayLen = {3.0};
    const arma::mat cLatHopEqDist = {5.5};
    const arma::uword cLatHopMaxDistOrder = 1;

    // CplHop
    // adjust cCplHopCoef to see different Gamma limit!
    const arma::cx_mat cCplHopCoef = cx_mat{ mat{-0.03}, zeros(1,1)};
    const arma::mat cCplHopDecayLen = mat{1};
    const arma::mat cCplHopEqDist = mat{4};
    const double cCplHopCutoff = 45.0;
    const uword cCplHopMaxExtOrder = 3;
 
    // LatPot
    //const arma::mat cLatPotSpringConst = mat{0.01}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotSpringConst = mat{1}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotEqDist = mat{5.5};
    const arma::uword cLatPotMaxDistOrder = 1;

    // CplPot
//    const mat cCplPotCoef = {5.8, 87.0}; // row: lat; col: mol
//    const mat cCplPotCoef = {11.7, 174.1}; // row: lat; col: mol
//    const mat cCplPotExpDecayLen = {0.63, 0.50};
//    const mat cCplPotCoef = {0.33};
    //const mat cCplPotCoef = {0.0159}; // estimated from Olander1971, use D = 10 kcal/mol
    //const mat cCplPotExpDecayLen = {0.75}; // use m = 4 ( l = re / m )
    const mat cCplPotCoef = {0.010}; // estimated from Olander1971, use D = 10 kcal/mol
    const mat cCplPotExpDecayLen = {0.8};
    const mat cCplPotEqDist = {5.67}; // use re = 3 A
    const double cCplPotCutoff = 45.0;
    const arma::uword cCplPotMaxExtOrder = 3;
}

namespace honeycomb
{
    // UnitCell
    //const vector<mat> cLatOnSiteE = { mat{0.00}, mat{0.00} };
    const vector<mat> cLatOnSiteE = { randu(3,3), randu(2,2) };
    const vector<vec> cLatInCellCoor = { vec{0.0, 0.0, 0.0},
					 vec{0.0, -1.0/ROOT3, 0.0}*4.64 };
    const vec cLatMassList = {22000, 21999};
    
    // BravLat
    const vector<vec> cLatVec = {   4.64*vec{1.0, 0.0, 0.0},
				    4.64*vec{0.5, 0.5*ROOT3, 0.0}   };
    
    // LatSupCell
    const uvec cNumUnitCell = {4, 4};
    const uvec cIsPeriodic = {1, 1};
    
    // BrilZone
    const uvec cNumK = {2, 2};

    // LatHop
    const arma::cx_mat cLatHopCoef = cx_mat{ mat{-0.01}, zeros(1,1) };
    // Chang2008 gold hopping amplitude -0.3056 eV ~ -0.01 a.u.
    const arma::mat cLatHopDecayLen = {3.0};
    const arma::mat cLatHopEqDist = {5.5};
    const arma::uword cLatHopMaxDistOrder = 1;

    // CplHop
    // adjust cCplHopCoef to see different Gamma limit!
    const arma::cx_mat cCplHopCoef = cx_mat{ mat{-0.003}, zeros(1,1)};
    const arma::mat cCplHopDecayLen = mat{3.0};
    const arma::mat cCplHopEqDist = mat{6.8}; // NO-Au geometry optimization yields N-Au 2.8 A,
					      // NO bond length 1.17 A, 47 degrees to normal
    const double cCplHopCutoff = 20.0;
    const uword cCplHopMaxExtOrder = 3;
 
    // LatPot
    const arma::mat cLatPotSpringConst = mat{0.01}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotEqDist = mat{5.5};
    const arma::uword cLatPotMaxDistOrder = 1;

    // CplPot
    const mat cCplPotCoef = {11.7, 174.1}; // row: lat; col: mol
    const mat cCplPotExpDecayLen = {0.63, 0.50};
    const double cCplPotCutoff = 20.0;
    const arma::uword cCplPotMaxExtOrder = 3;
   
}

namespace test
{
    const vector<vec> cMolCoor = { vec{0.00, 0.00, 0.00},
        			   vec{0.00, 0.00, 2.17} };
    const vec cMolMassList = {25704, 29376};

//    const vector<mat> cMolOnSiteE = { mat{}, mat{}, mat{} };
//    const vector<vec> cMolCoor = { vec{0.0, 0.0, 0.0},
//				   vec{0.0, 0.0, 0.2},
//				   vec{0.0, 0.0, 0.3} };
//    const vec cMolMassList = {1, 2, 3};
//    const vector<mat> cMolOnSiteVarCoef = { mat{}, mat{}, mat{} };
//    const vector<mat> cMolOnSiteVarLen = { mat{}, mat{}, mat{} };

    // UnitCell
    //const vector<mat> cLatOnSiteE = { mat{{0.05,0.01},{0.01, 0.03}} };
    const vector<mat> cLatOnSiteE = { mat{0.05} };
    // Chang2008 gold p-orbital on-site 1.4188 eV ~ 0.05 a.u.
    const vector<vec> cLatInCellCoor = { vec{0,0,0} };
    const vec cLatMassList = {361692}; // mass of gold atom ~ 197 u ~ 197*1836 a.u.

    // BravLat
    // bond length of gold ~ 288 pm ~ 5.5 a.u.
    const vector<vec> cLatVec = {   5.5*vec{1.0, 0.0, 0.0},
				    5.5*vec{0.5, 0.5*ROOT3, 0.0},
    				    5.5*vec{0.5, 0.5/ROOT3, -ROOT2/ROOT3} };
    
    // LatSupCell
    const uvec cNumUnitCell = {4, 4, 3};
    const uvec cIsPeriodic = {1, 1, 0};
    
    // BrilZone
    const uvec cNumK = {1, 1}; // number of elements shall equals to

    // LatHop
    //const arma::cx_mat cLatHopCoef = cx_mat{ mat{{-0.01, -0.03}, {-0.03, -0.02}}, zeros(2,2) };
    const arma::cx_mat cLatHopCoef = cx_mat{ mat{-0.01}, zeros(1,1) };
    // Chang2008 gold hopping amplitude -0.3056 eV ~ -0.01 a.u.
    const arma::mat cLatHopDecayLen = {3.0};
    //const arma::mat cLatHopDecayLen = {{3.0, 2.0}, {2.0, 1.0}};
    const arma::mat cLatHopEqDist = {5.5};
    //const arma::mat cLatHopEqDist = {{5.5, 4.4}, {4.4, 6.6}};
    const arma::uword cLatHopMaxDistOrder = 1;

    // CplHop
    // adjust cCplHopCoef to see different Gamma limit!
    const arma::cx_mat cCplHopCoef = cx_mat{ mat{-0.01}, zeros(1,1)};
    //const arma::cx_mat cCplHopCoef = cx_mat{ vec{-0.01, -0.02}, zeros(2,1)};
    const arma::mat cCplHopDecayLen = mat{3.0};
    //const arma::mat cCplHopDecayLen = vec{3.0, 4.0};
    const arma::mat cCplHopEqDist = mat{6.8}; // NO-Au geometry optimization N-Au 2.8 A,
    //const arma::mat cCplHopEqDist = vec{6.8, 7.7};
					      // NO bond length 1.17 A, 47 degrees to normal
    const double cCplHopCutoff = 20;
    const uword cCplHopMaxExtOrder = 3;
    
    // MolOnSite
    const double cMorse0Coef = 0.24;
    const double cMorse0EqDist = 2.15;
    const double cMorse0ExpDecayLen = 0.69;
    const double cMorse1Coef = 0.19;
    const double cMorse1EqDist = 2.44;
    const double cMorse1ExpDecayLen = 0.76;

    const double cImageC = 2.34;
    const double cImageD = 0.13;
    const double cRefZ = 0.00;
    const double cWorkFunc = 0.19;
    const double cElecAffinity = -0.00;

    // LatPot
    const arma::mat cLatPotSpringConst = mat{0.01}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotEqDist = mat{5.5};
    const arma::uword cLatPotMaxDistOrder = 1;

    // MolPot
    const double cMolPotCoef = 0.24;
    const double cMolPotEqDist = 2.17;
    const double cMolPotExpDecayLen = 0.69;
    // vibrational period 2*pi/sqrt(2*D/m/l/l) = 732

    // CplPot
    const mat cCplPotCoef = {11.7, 174.1}; // row: lat; col: mol
//    const mat cCplPotCoef = {0, 0}; // row: lat; col: mol
    const mat cCplPotExpDecayLen = {0.63, 0.50};
    const double cCplPotCutoff = 20.0;
    const arma::uword cCplPotMaxExtOrder = 3;
}
