#include <armadillo>
#include <chrono>
#include <iomanip>
#include "efld_all.h"
#include "omp.h"

using namespace arma;
using namespace std;
//using namespace test;
using namespace hexagonal;
using namespace monoatomic;
//using namespace diatomic;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    myclock::time_point main_start = myclock::now();

    arma_rng::set_seed_random();
//    arma_rng::set_seed(1);

//////////////////////////////////////////////
//
//
//	main (parallel)
//
//
//////////////////////////////////////////////

    uword num_trajs = 4; // cNumTraj
    uword num_trapped_trajs = 0;

    mat kin_loss = zeros(num_trajs, 3);
    vec time_elapsed = zeros(num_trajs);
    vec minimal_z = zeros(num_trajs);
    vec out_velo_angle = zeros(num_trajs);
//    mat vib_energy_quanta = zeros(num_trajs, 2);
//    vec rot_energy = zeros(num_trajs);

    #pragma omp parallel for schedule(dynamic)
    for (uword traj_idx = 0; traj_idx < num_trajs; ++traj_idx) {
        EFLD efld(
		cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo,
    		cIsPeriodic, cNumK,
    		cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, cLatHopMaxDistOrder,
    		cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, cCplHopCutoff, cCplHopMaxExtOrder,
    		cMorse0Coef, cMorse0EqDist, cMorse0ExpDecayLen,
    		cMorse1Coef, cMorse1EqDist, cMorse1ExpDecayLen,
    		cImageC, cImageD, cRefZ, cWorkFunc, cElecAffinity, 
    		cLatPotSpringConst, cLatPotEqDist, cLatPotMaxDistOrder, 
    		cMolPotCoef, cMolPotEqDist, cMolPotExpDecayLen,
    		cCplPotCoef, cCplPotExpDecayLen, cCplPotEqDist, cCplPotCutoff, cCplPotMaxExtOrder,
    		cLatFreeze, cLatZFreeze, cLevelBroadenWidth, cEnsemble, cThermBeta, cChemPot,
    		cFricSwitch, cDt, /*cNumTimeStep*/ cMaxNumStep ,
    		cFixedTerminal, cTerminalZ, cPrintSwitch );

	double init_mol_kinE = efld.sys.sup_cell.molKinE();
	efld.setInitLatBoltzmann();
	//efld.setInitMolVib(15.0);
	//vib_energy_quanta(traj_idx, 1) = efld.sys.numVibQuanta();
	//cout << "dt: " << efld.dt << endl;
	//cout << "init molpot: " << efld.sys.pot.mol_pot.energy() << endl;
	//cout << "numVibQuanta: " << efld.sys.numVibQuanta() << endl;
	//cout << "vib T: " << 2*PI/efld.sys.molVibOmega() << endl;

	myclock::time_point start = myclock::now();
	efld.propagate();
	chrono::duration<double> dur = myclock::now() - start;

	double minz = 20;
    	for (auto& coor : efld.data.coor_t) {
    	    double nowz = coor(2, efld.sys.sup_cell.numLatSites());
    	    if ( nowz < minz )
    	        minz = nowz;
    	}
	minimal_z(traj_idx) = minz;

	kin_loss(traj_idx, 0) = init_mol_kinE - efld.sys.sup_cell.molKinE();
	kin_loss(traj_idx, 1) = init_mol_kinE;

	time_elapsed(traj_idx) = dur.count();

	vec v = efld.data.velo_t[efld.data.velo_t.size()-1].tail_cols(1);
	out_velo_angle(traj_idx) = asin( norm(v.head_rows(2)) / norm(v) ) * 180.0 / PI;

	//vib_energy_quanta(traj_idx, 0) = efld.sys.numVibQuanta();

	if ( !efld.data.success )
	    num_trapped_trajs += 1;

	cout << endl << endl;
	cout << "minz : " << minz << endl
	     << "kin_loss: " << kin_loss(traj_idx, 0) / init_mol_kinE * 100 << "%" << endl
	     << "trapped: " << num_trapped_trajs << endl
	     << endl;

	//conv_to<vec>::from(efld.data.run_time).save("time.txt", raw_ascii);
    	//conv_to<vec>::from(efld.data.totE_t).save("energy.txt", raw_ascii);
	//extract(efld.data.coor_t, 0).save("xcoor.txt", raw_ascii);
    	//extract(efld.data.coor_t, 1).save("ycoor.txt", raw_ascii);
    	//extract(efld.data.coor_t, 2).save("zcoor.txt", raw_ascii);
	
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
  
    minimal_z.save("minimal_z.txt", raw_ascii);
    kin_loss.col(2) = kin_loss.col(0) / kin_loss.col(1);
    kin_loss.save("kin_loss.txt", raw_ascii);
    time_elapsed.save("time_elapsed.txt", raw_ascii);
    out_velo_angle.save("out_angle.txt", raw_ascii);
    //vib_energy_quanta.save("vib_energy.txt", raw_ascii);
    cout << "num_trapped_trajs = " << num_trapped_trajs << endl;



    chrono::duration<double> main_dur = myclock::now() - main_start;
    cout << "main_dur = " << main_dur.count() << endl;

    return 0;
}

