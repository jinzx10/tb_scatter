#ifndef __MODEL_CONSTANTS_H__
#define __MODEL_CONSTANTS_H__

#include <armadillo>
#include <vector>

extern const double cThermBeta;
extern const double cChemPot;
extern const double cLevelBroadenWidth;
extern const bool cLatFreeze;
extern const bool cLatZFreeze;

namespace hexagonal
{
    extern const std::vector<arma::vec> cLatInCellCoor;
    extern const arma::vec cLatMassList;
    extern const std::vector<arma::vec> cLatVec;
    extern const arma::uvec cNumUnitCell;

    extern const arma::uvec cIsPeriodic;
    extern const arma::uvec cNumK;

    extern const std::vector<arma::cx_mat> cLatOnSiteE;
    extern const arma::cx_mat cLatHopCoef;
    extern const arma::mat cLatHopDecayLen;
    extern const arma::mat cLatHopEqDist;
    extern const arma::uword cLatHopMaxDistOrder;

    extern const arma::mat cLatPotSpringConst;
    extern const arma::mat cLatPotEqDist;
    extern const arma::uword cLatPotMaxDistOrder;

    extern const std::vector<arma::cx_vec> cCplHopCoef;
    extern const std::vector<arma::vec> cCplHopDecayLen;
    extern const std::vector<arma::vec> cCplHopEqDist;
    extern const double cCplHopCutoff;
    extern const arma::uword cCplHopMaxExtOrder;

    namespace monoatomic
    {
	extern const std::vector<arma::vec> cMolCoor;
    	extern const arma::vec cMolMassList;

	// should contain some other CplPot stuff
	extern const double cCplPotCutoff;
	extern const arma::uword cCplPotMaxExtOrder;

    	extern const double cImageCoef;
    	extern const double cImageRegLength;
    	extern const double cImageRefZ;
    	extern const double cWorkFunc;
    	extern const double cElecAffinity;
    }

    namespace Tully2009_NO_Au
    {
	extern const std::vector<arma::vec> cMolCoor;
	extern const arma::vec cMolMassList;

	extern const double cNeutralMorseCoef;
	extern const double cNeutralMorseEqDist;
	extern const double cNeutralMorseExpDecayLen;
	extern const double cNeutralAuNCoef;
	extern const double cNeutralAuOCoef;
	extern const double cNeutralAuNExpDecayLen;
	extern const double cNeutralAuOExpDecayLen;
	extern const double cIonicMorseCoef;
	extern const double cIonicMorseEqdist;
	extern const double cIonicMorseExpDecayLen;
	extern const double cIonicAuNCoef;
	extern const double cIonicAuOCoef;
	extern const double cIonicAuNExpDecayLen;
	extern const double cIonicAuOExpDecayLen;
	extern const double cIonicAuNEqDist;
	extern const double cImageCoef;
	extern const double cImageRegLength;
	extern const double cImageRefZ;
	extern const double cWorkFunc;
	extern const double cElecAffinity;
	extern const double cCplPotCutoff;
	extern const arma::uword cCplPotMaxExtOrder;
    }
}

#endif
