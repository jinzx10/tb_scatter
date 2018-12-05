#include <armadillo>
#include <chrono>
#include <mpi.h>
#include <sstream>
#include "common.h"
#include "NO_Au_Diabats.h"
#include "CME.h"

using namespace arma;
using namespace std;
using namespace hexagonal::Tully2009_NO_Au;
using namespace hexagonal;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    const bool broaden_switch = true;
    const double dt = 5.0;
    const uword max_num_steps = 10000;
    const double terminal_com_z = 15.001;
    const double mol_init_com_z = 15;
    const bool print_switch = true;

    const int incident = 0;
    const int num_trajs = 1;
    const double mol_init_trans_kinE = 0.02;
    const double mol_init_vib_quanta = 15.0;
    const double mol_init_rot_quanta = 0;

    const double cpl = abs(cCplHopCoef[0](0));

    string dir_name = "";
    if (num_trajs == 1)
	dir_name = "test_data/onetraj/";
    else
	dir_name = "data/";
    string data_file_name = dir_name;

    int num_procs;
    int id;

    int* is_scattered = nullptr;
    int* final_state = nullptr; // 1 = ionic; 0 = neutral
    int* num_hops = nullptr;
    int* num_steps = nullptr;
    double* time_elapsed = nullptr;
    double* walltime_elapsed = nullptr;
    double* init_trans_kinE = nullptr;
    double* init_rot_quanta = nullptr;
    double* init_vib_quanta = nullptr;
    double* final_trans_kinE = nullptr;
    double* final_rot_quanta = nullptr;
    double* final_vib_quanta = nullptr;
    double* minimal_com_z = nullptr;
    double* reflection_angle = nullptr;
    double* init_com_x = nullptr;
    double* init_com_y = nullptr;
    double* init_com_z = nullptr;
    double* init_com_vx = nullptr;
    double* init_com_vy = nullptr;
    double* init_com_vz = nullptr;

    mat data;
    
    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    const int num_trajs_per_proc = num_trajs / num_procs;
    assert( num_trajs % num_procs == 0 );

    myclock::time_point main_start; // let the master proc do the data collection and timing
    chrono::duration<double> main_dur;

    if (id == 0) {
	main_start = myclock::now();

	is_scattered = new int[num_trajs];
	final_state = new int[num_trajs];
    	num_hops = new int[num_trajs];
	num_steps = new int[num_trajs];
    	time_elapsed = new double[num_trajs];
    	walltime_elapsed = new double[num_trajs];
	init_trans_kinE = new double[num_trajs];
	init_rot_quanta = new double[num_trajs];
	init_vib_quanta = new double[num_trajs];
	final_trans_kinE = new double[num_trajs];
	final_rot_quanta = new double[num_trajs];
	final_vib_quanta = new double[num_trajs];
	minimal_com_z = new double[num_trajs];
	reflection_angle = new double[num_trajs];
	init_com_x = new double[num_trajs];
	init_com_y = new double[num_trajs];
	init_com_z = new double[num_trajs];
	init_com_vx = new double[num_trajs];
	init_com_vy = new double[num_trajs];
	init_com_vz = new double[num_trajs];
    }

    int local_is_scattered[num_trajs_per_proc];
    int local_final_state[num_trajs_per_proc];
    int local_num_hops[num_trajs_per_proc];
    int local_num_steps[num_trajs_per_proc];
    double local_time_elapsed[num_trajs_per_proc];
    double local_walltime_elapsed[num_trajs_per_proc];
    double local_init_trans_kinE[num_trajs_per_proc];
    double local_init_rot_quanta[num_trajs_per_proc];
    double local_init_vib_quanta[num_trajs_per_proc];
    double local_final_trans_kinE[num_trajs_per_proc];
    double local_final_rot_quanta[num_trajs_per_proc];
    double local_final_vib_quanta[num_trajs_per_proc];
    double local_minimal_com_z[num_trajs_per_proc];
    double local_reflection_angle[num_trajs_per_proc];
    double local_init_com_x[num_trajs_per_proc];
    double local_init_com_y[num_trajs_per_proc];
    double local_init_com_z[num_trajs_per_proc];
    double local_init_com_vx[num_trajs_per_proc];
    double local_init_com_vy[num_trajs_per_proc];
    double local_init_com_vz[num_trajs_per_proc];


    arma_rng::set_seed_random();

    for (int local_traj_idx = 0; local_traj_idx < num_trajs_per_proc; ++local_traj_idx) {
	System< NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
		cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		cMolCoor, cMolMassList, vec{0,0,0}, vec{0,0,0}, cIsPeriodic, cNumK,
		LatTB( cLatOnSiteE, cLatHopCoef, cLatHopDecayLen, cLatHopEqDist,
		       nullptr, mat{}, {}, cLatHopMaxDistOrder ),
		CplHop( cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, 
			nullptr, nullptr, nullptr, {},
	    	        cCplHopCutoff, cCplHopMaxExtOrder ),
	    	LatPot( cLatPotSpringConst, cLatPotEqDist,
	    	        nullptr, mat{}, {}, cLatPotMaxDistOrder ),
	    	NO_Au_Diabats(
		    cNeutralMorseCoef, cNeutralMorseEqDist, cNeutralMorseExpDecayLen, 
	    	    cNeutralAuNCoef, cNeutralAuOCoef,
	    	    cNeutralAuNExpDecayLen, cNeutralAuOExpDecayLen,
		    cIonicMorseCoef, cIonicMorseEqdist, cIonicMorseExpDecayLen, 
		    cIonicAuNCoef, cIonicAuOCoef,
	    	    cIonicAuNExpDecayLen, cIonicAuOExpDecayLen, cIonicAuNEqDist, 
	    	    cImageCoef, cImageRegLength, cImageRefZ, cWorkFunc, cElecAffinity, 
	    	    nullptr, nullptr, nullptr, nullptr, {}, cCplPotCutoff, 1 ),
		cLevelBroadenWidth );

	sys.randSetLatVelo( cThermBeta );

	vec target = ( sys.latCoor(6) + sys.latCoor(9) ) / 2.0;
	sys.randSetMolComCoorVelo( mol_init_com_z, target, mol_init_trans_kinE, incident ); 

	sys.diabats.randSetMolOrientation();
	sys.diabats.setMolVibQuanta(mol_init_vib_quanta);
	sys.diabats.setMolRotQuanta(mol_init_rot_quanta);


	CME< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > cme(
		&sys, cThermBeta, cChemPot, broaden_switch, cLatFreeze, cLatZFreeze, 
		dt, max_num_steps, terminal_com_z, print_switch );


	local_init_com_x[local_traj_idx] = sys.molComCoor()(0);
	local_init_com_y[local_traj_idx] = sys.molComCoor()(1);
	local_init_com_z[local_traj_idx] = sys.molComCoor()(2);
	local_init_com_vx[local_traj_idx] = sys.molComVelo()(0); 
	local_init_com_vy[local_traj_idx] = sys.molComVelo()(1); 
	local_init_com_vz[local_traj_idx] = sys.molComVelo()(2); 

    	myclock::time_point start = myclock::now();
    	cme.propagate();
    	chrono::duration<double> dur = myclock::now() - start;

	local_is_scattered[local_traj_idx] = cme.is_scattered;
	local_final_state[local_traj_idx] = cme.state;
	local_num_hops[local_traj_idx] = cme.num_hops;
	local_num_steps[local_traj_idx] = cme.run_time.size();
	local_time_elapsed[local_traj_idx] = cme.curr_time;
	local_walltime_elapsed[local_traj_idx] = dur.count();

	local_init_trans_kinE[local_traj_idx] = mol_init_trans_kinE;
	local_final_trans_kinE[local_traj_idx] = sys.molComKinE();
	local_init_rot_quanta[local_traj_idx] = mol_init_rot_quanta;
	local_final_rot_quanta[local_traj_idx] = sys.diabats.numRotQuanta(cme.state);
	local_init_vib_quanta[local_traj_idx] = mol_init_vib_quanta;
	local_final_vib_quanta[local_traj_idx] = sys.diabats.numVibQuanta(cme.state);
	local_minimal_com_z[local_traj_idx] = cme.minimal_z;

	vec v_out = sys.molComVelo();
	local_reflection_angle[local_traj_idx] = acos( v_out(2) / norm(v_out) ) * 180.0 / PI;

	if (num_trajs == 1) {
    	    extract(cme.lat_coor_t, 0).save( dir_name+"cme_lat_coor_x.txt", raw_ascii );
    	    extract(cme.lat_coor_t, 1).save( dir_name+"cme_lat_coor_y.txt", raw_ascii );
    	    extract(cme.lat_coor_t, 2).save( dir_name+"cme_lat_coor_z.txt", raw_ascii );
    	    extract(cme.mol_coor_t, 0).save( dir_name+"cme_mol_coor_x.txt", raw_ascii );
    	    extract(cme.mol_coor_t, 1).save( dir_name+"cme_mol_coor_y.txt", raw_ascii );
    	    extract(cme.mol_coor_t, 2).save( dir_name+"cme_mol_coor_z.txt", raw_ascii );
	    conv_to<vec>::from(cme.run_time).save( dir_name+"cme_run_time.txt", raw_ascii );
	}
    }

    MPI_Gather( local_is_scattered, num_trajs_per_proc, MPI_INT, 
		is_scattered, num_trajs_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather( local_final_state, num_trajs_per_proc, MPI_INT, 
		final_state, num_trajs_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather( local_num_hops, num_trajs_per_proc, MPI_INT, 
		num_hops, num_trajs_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather( local_num_steps, num_trajs_per_proc, MPI_INT, 
		num_steps, num_trajs_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather( local_time_elapsed, num_trajs_per_proc, MPI_DOUBLE, 
		time_elapsed, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_walltime_elapsed, num_trajs_per_proc, MPI_DOUBLE, 
		walltime_elapsed, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_trans_kinE, num_trajs_per_proc, MPI_DOUBLE, 
		init_trans_kinE, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_rot_quanta, num_trajs_per_proc, MPI_DOUBLE, 
		init_rot_quanta, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_vib_quanta, num_trajs_per_proc, MPI_DOUBLE, 
		init_vib_quanta, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_final_trans_kinE, num_trajs_per_proc, MPI_DOUBLE, 
		final_trans_kinE, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_final_rot_quanta, num_trajs_per_proc, MPI_DOUBLE, 
		final_rot_quanta, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_final_vib_quanta, num_trajs_per_proc, MPI_DOUBLE, 
		final_vib_quanta, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_minimal_com_z, num_trajs_per_proc, MPI_DOUBLE, 
		minimal_com_z, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_reflection_angle, num_trajs_per_proc, MPI_DOUBLE, 
		reflection_angle, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_x, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_x, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_y, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_y, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_z, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_z, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_vx, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_vx, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_vy, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_vy, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( local_init_com_vz, num_trajs_per_proc, MPI_DOUBLE, 
		init_com_vz, num_trajs_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
    if (id == 0) {
	data = zeros(num_trajs, 20);
	data.col(0) = conv_to<vec>::from( Col<int>(is_scattered, num_trajs) );
	data.col(1) = conv_to<vec>::from( Col<int>(final_state, num_trajs) );
	data.col(2) = conv_to<vec>::from( Col<int>(num_hops, num_trajs) );
	data.col(3) = conv_to<vec>::from( Col<int>(num_steps, num_trajs) );
	data.col(4) = vec(time_elapsed, num_trajs);
	data.col(5) = vec(walltime_elapsed, num_trajs);
	data.col(6) = vec(init_trans_kinE, num_trajs);
	data.col(7) = vec(init_rot_quanta, num_trajs);
	data.col(8) = vec(init_vib_quanta, num_trajs);
	data.col(9) = vec(final_trans_kinE, num_trajs);
	data.col(10) = vec(final_rot_quanta, num_trajs);
	data.col(11) = vec(final_vib_quanta, num_trajs);
	data.col(12) = vec(minimal_com_z, num_trajs);
	data.col(13) = vec(reflection_angle, num_trajs);
	data.col(14) = vec(init_com_x, num_trajs);
	data.col(15) = vec(init_com_y, num_trajs);
	data.col(16) = vec(init_com_z, num_trajs);
	data.col(17) = vec(init_com_vx, num_trajs);
	data.col(18) = vec(init_com_vy, num_trajs);
	data.col(19) = vec(init_com_vz, num_trajs);

	delete[] is_scattered;
	delete[] final_state;
    	delete[] num_hops;
    	delete[] num_steps;
    	delete[] time_elapsed;
    	delete[] walltime_elapsed;
    	delete[] init_trans_kinE;
    	delete[] init_rot_quanta;
    	delete[] init_vib_quanta;
    	delete[] final_trans_kinE;
    	delete[] final_rot_quanta;
    	delete[] final_vib_quanta;
	delete[] minimal_com_z;
	delete[] reflection_angle;
	delete[] init_com_x;
	delete[] init_com_y;
	delete[] init_com_z;
	delete[] init_com_vx;
	delete[] init_com_vy;
	delete[] init_com_vz;


	if (broaden_switch)
	    data_file_name += "BCME_";
	else
	    data_file_name += "CME_";

	stringstream sstream;

	sstream.str("");
	sstream << incident;
	data_file_name += "angle" + sstream.str() + "_";
	
	sstream.str("");
	sstream << mol_init_trans_kinE;
	data_file_name += "trans" + sstream.str() + "_";
	
	sstream.str("");
	sstream << mol_init_vib_quanta;
	data_file_name += "vib" + sstream.str() + "_";
	
	sstream.str("");
	sstream << mol_init_rot_quanta;
	data_file_name += "rot" + sstream.str() + "_";
	
	sstream.str("");
	sstream << cpl;
	data_file_name += "cpl" + sstream.str() + "_";
	
	sstream.str("");
	sstream << dt;
	data_file_name += "dt" + sstream.str() + "_";
	
	sstream.str("");
	sstream << num_trajs;
	data_file_name += "traj" + sstream.str() + ".txt";

	data.save(data_file_name, raw_ascii);

	main_dur = myclock::now() - main_start;
	cout << "main_dur = " << main_dur.count() << endl;
    }

    MPI_Finalize();

    return 0;
}
