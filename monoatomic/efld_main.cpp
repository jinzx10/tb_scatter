#include <armadillo>
#include <chrono>
#include <iomanip>
#include "efld_all.h"
#include <omp.h>

using namespace arma;
using namespace std;
using namespace hexagonal;
using namespace monoatomic;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    myclock::time_point main_start = myclock::now();

    arma_rng::set_seed_random();

//////////////////////////////////////////////
//
//
//	main
//
//
//////////////////////////////////////////////

    uword num_trajs = 1; // cNumTraj
    uword num_trapped_trajs = 0;

    mat data = zeros(num_trajs, 7);
    // 0: init kinE
    // 1: kinE loss
    // 2: kinE loss percent
    // 3: time elapsed
    // 4: minimal z
    // 5: angle of reflection
    // 6: trapped or not
//    mat vib_energy_quanta = zeros(num_trajs, 2);
//    vec rot_energy = zeros(num_trajs);

//    #pragma omp parallel for num_threads(4) shared(num_trapped_trajs) schedule(dynamic, 1)
    for (uword traj_idx = 0; traj_idx < num_trajs; ++traj_idx) {
	arma_rng::set_seed_random();
	arma_rng::set_seed(1);
        EFLD efld(
		cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo,
    		cIsPeriodic, cNumK,
    		cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, cLatHopMaxDistOrder,
    		cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, cCplHopCutoff, cCplHopMaxExtOrder,
    		cImageC, cImageD, cRefZ, cWorkFunc, cElecAffinity, 
    		cLatPotSpringConst, cLatPotEqDist, cLatPotMaxDistOrder, 
    		cMolPotCoef, cMolPotEqDist, cMolPotExpDecayLen,
    		cCplPotCoef, cCplPotExpDecayLen, cCplPotEqDist, cCplPotCutoff, cCplPotMaxExtOrder,
    		cLatFreeze, cLatZFreeze, cLevelBroadenWidth, cEnsemble, cThermBeta, cChemPot,
    		cFricSwitch, cDt, /*cNumTimeStep*/ cMaxNumStep ,
    		cTerminalZ, cPrintSwitch );

	double init_mol_kinE = efld.sys.sup_cell.molKinE();
	efld.setInitLatBoltzmann();
	//efld.setInitMolVib(15.0);
	//vib_energy_quanta(traj_idx, 1) = efld.sys.numVibQuanta();
	//cout << "dt: " << efld.dt << endl;
	//cout << "init molpot: " << efld.sys.pot.mol_pot.energy() << endl;
	//cout << "numVibQuanta: " << efld.sys.numVibQuanta() << endl;
	//cout << "vib T: " << 2*PI/efld.sys.molVibOmega() << endl;

	vec init_posit = genInUnitCell( vec{12.375, 7.1447, 15}, cLatVec[0], cLatVec[1]);
	init_posit(2) = 15.0;
	efld.sys.sup_cell.coor.tail_cols(1) = init_posit;
	efld.sys.sup_cell.coor.tail_cols(1) = vec{11.00, 6.3509, 15};

	myclock::time_point start = myclock::now();
	efld.computeAndSaveInit();
	efld.propagate();
	chrono::duration<double> dur = myclock::now() - start;

	double minz = 20;
    	for (auto& coor : efld.data.coor_t) {
    	    double nowz = coor(2, efld.sys.sup_cell.numLatSites());
    	    if ( nowz < minz )
    	        minz = nowz;
    	}

	// 0: kinE loss
    	// 1: init kinE 
    	// 2: kinE loss percent
    	// 3: time elapsed
    	// 4: minimal z
    	// 5: angle of reflection
    	// 6: trapped or not
	data(traj_idx, 0) = init_mol_kinE - efld.sys.sup_cell.molKinE();
	data(traj_idx, 1) = init_mol_kinE;
	// column 2 is computed from column 0 and 1

	data(traj_idx, 3) = dur.count();
	data(traj_idx, 4) = minz;

	vec v = efld.data.velo_t[efld.data.velo_t.size()-1].tail_cols(1);
	data(traj_idx, 5) = asin( norm(v.head_rows(2)) / norm(v) ) * 180.0 / PI;
	data(traj_idx, 6) = efld.data.success; // scattered: true; trapped: false

	//vib_energy_quanta(traj_idx, 0) = efld.sys.numVibQuanta();

	if ( !efld.data.success )
	    num_trapped_trajs += 1;

	cout << endl 
	     << "dur: " << dur.count() << endl
	     << "init position: " << init_posit.t()
	     << "num steps: " << efld.data.coor_t.size() << endl
	     << "minz : " << minz << endl
	     << "kin_loss: " << data(traj_idx, 0) / init_mol_kinE * 100 << "%" << endl
	     //<< "trapped: " << num_trapped_trajs << endl
	     << endl;

	if (!cConvergenceTest) {
    	    conv_to<vec>::from(efld.data.run_time).save("time.txt", raw_ascii);
    		conv_to<vec>::from(efld.data.totE_t).save("energy.txt", raw_ascii);
    	} else {
    	    conv_to<vec>::from(efld.data.run_time).save("time_test.txt", raw_ascii);
    		conv_to<vec>::from(efld.data.totE_t).save("energy_test.txt", raw_ascii);
    	}

    	conv_to<vec>::from(efld.data.molCharge_t).save("charge.txt", raw_ascii);
    	extract(efld.data.coor_t, 0).save("xcoor.txt", raw_ascii);
    	extract(efld.data.coor_t, 1).save("ycoor.txt", raw_ascii);
    	extract(efld.data.coor_t, 2).save("zcoor.txt", raw_ascii);
    }
  
    data.col(2) = data.col(0) / data.col(1);
    data.save("data.txt", raw_ascii);
    //vib_energy_quanta.save("vib_energy.txt", raw_ascii);
    //cout << "num_trapped_trajs = " << num_trapped_trajs << endl;


    chrono::duration<double> main_dur = myclock::now() - main_start;
    cout << "main_dur = " << main_dur.count() << endl;

    return 0;
}
