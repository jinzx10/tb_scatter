#include "modelconst.h"
#include "mathconst.h"

using namespace arma;
using namespace std;

// thermodynamic quantity 
const double cThermBeta = 1000.0;
const double cChemPot = 0.0; // unused in canonical ensemble

// broaden the energy level to estimate DOS
const double cLevelBroadenWidth = 0.02;

// freeze the whole lattice
const bool cLatFreeze = true;

// restrict the movements of lattice sites to the first 2 dimensions
const bool cLatZFreeze = false;

namespace hexagonal
{
    // UnitCell
    const vector<vec> cLatInCellCoor = { vec{0,0,0} };
    const vec cLatMassList = {359070}; // mass of gold atom
    // [111] plane, bond length of gold ~ 288 pm ~ 5.5 a.u.
    // Tully2009 bond length ~ 295 pm ~ 5.58 a.u.
    const vector<vec> cLatVec = {   5.5*vec{1.0, 0.0, 0.0},
				    5.5*vec{0.5, 0.5*ROOT3, 0.0},
    				    5.5*vec{0.5, 0.5/ROOT3, -ROOT2/ROOT3} };
    // [110] plane
//    const vector<vec> cLatVec = {   5.5*vec{1.0, 0.0, 0.0},
//				    5.5*vec{0.0, ROOT2, 0.0},
//				    5.5*vec{0.5, 1.0/ROOT2, -0.5}  };

    // [111]
    const uvec cNumUnitCell = {4, 4, 3};
    // [110]
    // const uvec cNumUnitCell = {4, 4, 2};

    const uvec cIsPeriodic = {1, 1, 0};
    const uvec cNumK = {4, 4};

    // see Koskinen2006 (and Chang2008 ?)
    //const vector<cx_mat> cLatOnSiteE = { cx_mat{mat{{-0.215, 0}, {0, -0.01}}, zeros(2,2)} };
    const vector<cx_mat> cLatOnSiteE = { cx_mat{-0.01} };

    // LatHop
    //const arma::cx_mat cLatHopCoef = cx_mat{mat{ {-0.06, 0.09}, {0.09, 0.06} }, zeros(2, 2) };
    //const arma::mat cLatHopDecayLen = { {2.0, 2.0}, {2.0, 2.0} };
    //const arma::mat cLatHopEqDist = { {5.5, 5.5}, {5.5, 5.5} };
    const arma::cx_mat cLatHopCoef = cx_mat{ mat{0.06}, zeros(1,1) };
    const arma::mat cLatHopDecayLen = {2.0};
    const arma::mat cLatHopEqDist = {5.5};
    const arma::uword cLatHopMaxDistOrder = 1;

    // LatPot
    //const arma::mat cLatPotSpringConst = mat{0.01}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotSpringConst = mat{1}; // 1 a.u. = 1557 N/m
    const arma::mat cLatPotEqDist = mat{5.5};
    const arma::uword cLatPotMaxDistOrder = 1;

    // CplHop
    //const std::vector<arma::cx_vec> cCplHopCoef = {cx_vec{vec{0.02, 0.03}, zeros<vec>(2)}};
    //const std::vector<arma::vec> cCplHopDecayLen = {vec{0.4, 0.6}};
    //const std::vector<arma::vec> cCplHopEqDist = {vec{5.0, 5.5}};
    const std::vector<arma::cx_vec> cCplHopCoef = std::vector<cx_vec>{cx_vec{0.07}};
    const std::vector<arma::vec> cCplHopDecayLen = std::vector<vec>{vec{0.5}};
    const std::vector<arma::vec> cCplHopEqDist = std::vector<vec>{vec{5.5}};
    const double cCplHopCutoff = 25.0;
    const uword cCplHopMaxExtOrder = 3;

    namespace Tully2009_NO_Au
    {
	const vector<vec> cMolCoor = { vec{0.0, 0.0, 0.0},
				       vec{0.0, 0.0, 2.1753} };
    	const vec cMolMassList = { 25704.0, 29376.0 };

	const double cNeutralMorseCoef = 0.2432;
	const double cNeutralMorseEqDist = 2.1753;
	const double cNeutralMorseExpDecayLen = 0.6891;
	const double cNeutralAuNCoef = 11.6957;
	const double cNeutralAuOCoef = 174.0792;
	const double cNeutralAuNExpDecayLen = 0.6284;
	const double cNeutralAuOExpDecayLen = 0.5028;
	const double cIonicMorseCoef = 0.1889;
	const double cIonicMorseEqdist = 2.4393;
	const double cIonicMorseExpDecayLen = 0.7595;
	const double cIonicAuNCoef = 0.0092;
	const double cIonicAuOCoef = 174.0792;
	const double cIonicAuNExpDecayLen = 0.9621;
	const double cIonicAuOExpDecayLen = 0.5028;
	const double cIonicAuNEqDist = 4.4406;
	const double cImageCoef = 0.25;
	const double cImageRegLength = 2.3484;
	const double cImageRefZ = 2.1807;
	const double cWorkFunc = 0.1874;
	const double cElecAffinity = 0.0003;
	const double cCplPotCutoff = 25.0;
	const arma::uword cCplPotMaxExtOrder = 3;
    }
}
