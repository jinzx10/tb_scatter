#include <armadillo>
#include <chrono>
#include <iomanip>
#include "bcme_all.h"
#include <omp.h>
#define ARMA_DONT_USE_OPENMP

using namespace arma;
using namespace std;
using namespace hexagonal;
using namespace monoatomic;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    myclock::time_point start = myclock::now();
	
#pragma omp parallel num_threads(4)
{
    cout << omp_get_num_threads() << endl;
#pragma omp for schedule(dynamic, 5)
    for (uword traj = 0; traj < cNumTraj; ++traj) {
        arma_rng::set_seed_random();
    
        BCME bcme( cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		   cMolOnSiteE, cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo,
		   cMolOnSiteVarCoef, cMolOnSiteVarLen, cIsPeriodic, cNumK, 
		   cHopCoef, cHopDecayLen, cHopEqDist, cLatMaxNbOrder,
	       	   cCplHopCutoff, cCplHopExtOrder, cCplPotCoef, cCplPotPowerDecay,
		   cCplPotExpDecayLen, cCplPotCharLen, cCplPotCutoff, cCplPotExtOrder,
		   cLatFreeze, cLatZFreeze, cThermBeta, cChemPot, cBroadenSwitch, cDt,
		   cMaxNumStep, cFixedTerminal, cTerminalZ, cPrintSwitch );

        double init_kinE = bcme.nucl.molKinE();
        bcme.propagate();
        cout << "initial kinE: " << setw(10) << init_kinE
	     << "  final kinE: " << setw(10) << bcme.nucl.molKinE()
	     << "  final state: " << bcme.data.state_t[bcme.data.state_t.size()-1]
	     << "  success: " << bcme.data.success << endl;
    }
}
        chrono::duration<double> dur = myclock::now() - start;
        cout << "time spent: " << dur.count() << endl;

//    bcme.setInit();

//    cout.precision(10);
//    bcme.sys.H.molOnSite().raw_print(std::cout); cout << endl;
//    bcme.sys.H.rawMolOnSite().raw_print(std::cout); cout << endl;
//    bcme.sys.var.var(48).raw_print(std::cout); cout << endl;
//    bcme.sys.var.diff(48).raw_print(std::cout); cout << endl;


//    if (cNumTimeStep<= 1000) {
//	cRunTime.save("time.txt", raw_ascii);
//    	conv_to<vec>::from(bcme.data.totE_t).save("energy.txt", raw_ascii);
//    } else {
//	cRunTime.save("time_test.txt", raw_ascii);
//    	conv_to<vec>::from(bcme.data.totE_t).save("energy_test.txt", raw_ascii);
//    }

//    conv_to<vec>::from(bcme.data.state_t).save("state.txt", raw_ascii);
//    extract(bcme.data.coor_t, 0).save("xcoor.txt", raw_ascii);
//    extract(bcme.data.coor_t, 1).save("ycoor.txt", raw_ascii);
//    extract(bcme.data.coor_t, 2).save("zcoor.txt", raw_ascii);

    return 0;
}
