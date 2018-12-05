#include "common.h"
#include <armadillo>
#include <sstream>
#include "EFLD.h"
#include "CME.h"
#include "NO_Au_Diabats.h"

using namespace arma;
using namespace std;
using namespace hexagonal;
using namespace hexagonal::Tully2009_NO_Au;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    myclock::time_point main_start = myclock::now();
    chrono::duration<double> main_dur;
    vec zerovec = {0,0,0};
    cout << "start test " << endl;

    bool Atoms_Switch = 0;
    bool BravLat_Switch = 0;
    bool LatSupCell_Switch = 0;
    bool SupCell_Switch = 0;
    bool LatTB_Switch = 0;
    bool CplHop_Switch = 0;
    bool LatPot_Switch = 0;
    bool NO_Au_Diabats_Switch = 0;
    bool BrilZone_Switch = 0;
    bool System_Switch = 0;
    bool E_dep_Switch = 0;
    bool Z_dep_Switch = 0;
    bool Two_Diabats_Switch = 0;
    bool Check_Diabats_Switch = 0;
    bool Broadened_Diabats_Switch = 0;
    bool Adiabats_Switch = 0;
    bool CME_Switch = 0;
    bool EFLD_Switch = 1;
    bool CME_Dynamics_Switch = 0;
    bool BO_Dynamics_Switch = 0;
    bool Sample_Traj_Switch = 0;

    std::vector<arma::vec> sup_lat_vec_list{};
    for (arma::uword d = 0; d != cIsPeriodic.size(); ++d)
	if (cIsPeriodic[d])
	    sup_lat_vec_list.push_back( cNumUnitCell[d]*cLatVec[d] );

///////////////////////////////////////////////////////////////
//
//	    Atoms
//
///////////////////////////////////////////////////////////////

if (Atoms_Switch) {
    cout << "Atoms start" << endl;

    Atoms unit_cell;
    unit_cell = Atoms( cLatInCellCoor, cLatMassList );

    cout << unit_cell.numSites() << endl;
    cout << unit_cell.coor() << endl;
    cout << unit_cell.coor(0) << endl;
    cout << unit_cell.mass() << endl;
    cout << unit_cell.mass(0) << endl;
    cout << unit_cell.totMass() << endl;
    //cout << unit_cell.reducedMass() << endl;


    Atoms mol;
    mol = Atoms( cMolCoor, cMolMassList );
    cout << mol.numSites() << endl;
    cout << mol.coor() << endl;
    cout << mol.coor(1) << endl;
    cout << mol.mass() << endl;
    cout << mol.mass(1) << endl;
    cout << mol.totMass() << endl;
    cout << mol.reducedMass() << endl;

    cout << "Atoms end" << endl << endl << endl;;
}


///////////////////////////////////////////////////////////////
//
//	    BravLat
//
///////////////////////////////////////////////////////////////

if (BravLat_Switch) {
    cout << "BravLat start" << endl;

    BravLat bv;
    bv = BravLat( cLatVec );

    cout << bv.latVec(0) << endl;
    cout << bv.rcpVec(2) << endl;
    cout << bv.dim() << endl;
    cout << bv.latVec() << endl;
    cout << bv.rcpVec() << endl;
    cout << dot(bv.latVec(0), bv.rcpVec(1)) << endl;
    cout << dot(bv.latVec(0), bv.rcpVec(2)) << endl;
    cout << dot(bv.latVec(1), bv.rcpVec(0)) << endl;
    cout << dot(bv.latVec(1), bv.rcpVec(2)) << endl;
    cout << dot(bv.latVec(2), bv.rcpVec(0)) << endl;
    cout << dot(bv.latVec(2), bv.rcpVec(1)) << endl;

    cout << "BravLat end" << endl << endl << endl;;
}


///////////////////////////////////////////////////////////////
//
//	    LatSupCell
//
///////////////////////////////////////////////////////////////

if (LatSupCell_Switch) {
    cout << "LatSupCell start" << endl;

    LatSupCell lsc;
    lsc = LatSupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell );

    cout << lsc.numSites() << endl;
    cout << lsc.numUnitCells() << endl;
    cout << lsc.numUnitCells(2) << endl;
    cout << lsc.coor() << endl;
    cout << lsc.coor(2) << endl;
    cout << lsc.mass() << endl;
    cout << lsc.mass(1) << endl;
    cout << lsc.subLatIdx(3) << endl;
    cout << lsc.unit_cell.coor(0) << endl;
    cout << lsc.unit_cell.mass(0) << endl;
    cout << lsc.unit_cell.mass() << endl;
    cout << lsc.unit_cell.totMass() << endl;
    cout << lsc.brav_lat.dim() << endl;
    cout << lsc.brav_lat.latVec(1) << endl;
    cout << lsc.brav_lat.rcpVec(1) << endl;

    cout << "LatSupCell end" << endl << endl << endl;;
}


///////////////////////////////////////////////////////////////
//
//	    SupCell
//
///////////////////////////////////////////////////////////////

if (SupCell_Switch) {
    cout << "SupCell start" << endl;

    SupCell sc;
    sc = SupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, zerovec, zerovec );

    sc.setMolComCoorVelo(15, {0,0,0}, 0.1, 30, 0);
    sc.mol_coor.print(); cout << endl;
    sc.mol_velo.print(); cout << endl;
    sc.lat_coor.print(); cout << endl;
    sc.lat_velo.print(); cout << endl;

    cout << endl;
    cout << sc.numLatSites() << endl;
    cout << sc.numMolSites() << endl;
    cout << sc.latMassMat() << endl;
    cout << sc.molMassMat() << endl;
    cout << sc.latCoor(1) << endl;
    cout << sc.latVelo(2) << endl;
    cout << sc.molCoor(1) << endl;
    cout << sc.molVelo(0) << endl;
    cout << sc.molComVelo() << endl;
    cout << sc.molComCoor() << endl;
    cout << sc.molBondLength() << endl;
    cout << sc.molMomInertia() << endl;
    cout << sc.kinE() << endl;
    cout << sc.latKinE() << endl;
    cout << sc.molKinE() << endl;
    cout << sc.molComKinE() << endl;
    cout << sc.molVibKinE() << endl;
    cout << sc.molRotKinE() << endl;
    cout << endl;

    cout << "SupCell end" << endl << endl << endl;;
}

///////////////////////////////////////////////////////////////
//
//	    LatTB 
//
///////////////////////////////////////////////////////////////

if (LatTB_Switch) {
    cout << "LatTB start" << endl;
    LatTB lat_tb;
    lat_tb = LatTB( cLatOnSiteE, cLatHopCoef, cLatHopDecayLen, cLatHopEqDist,
		    nullptr, mat{}, {}, 1 );
    SupCell sc;
    sc = SupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, zerovec, zerovec );

    cout << lat_tb.isComplete() << endl;
    lat_tb.reset(&sc.lat_coor, sc.lat_sup_cell.coor(), sup_lat_vec_list);
    cout << lat_tb.isComplete() << endl;

    cout << lat_tb.ampl(0, 1, zerovec) << endl;
    cout << lat_tb.diff(0, 1, 0, zerovec) << endl;
    cout << lat_tb.ampl(0, 3, zerovec) << endl;
    cout << lat_tb.diff(0, 3, 0, zerovec) << endl;
    cout << lat_tb.H(zerovec).col(0) << endl;
    cout << lat_tb.subLatIdx(3) << endl;
    cout << lat_tb.numSites() << endl;
    cout << lat_tb.numBaseSites() << endl;
    cout << lat_tb.numUnitCells() << endl;
    cout << lat_tb.numOrbs() << endl;
    cout << lat_tb.numBasisOrbs() << endl;
    cout << lat_tb.numOrbs(1) << endl;
    cout << lat_tb.idxOrbs(2) << endl;
    cout << lat_tb.idxBasisOrbs(2) << endl;
    cout << lat_tb.neighbor(1) << endl;
    cout << lat_tb.neighbor.vecToBeNeighbor(1, 3) << endl;
    cout << lat_tb.neighbor.vecToBeNeighbor(1, 5) << endl;

    cout << "LatTB end" << endl << endl << endl;;
}

///////////////////////////////////////////////////////////////
//
//	    CplHop
//
///////////////////////////////////////////////////////////////

if (CplHop_Switch) {
    cout << "CplHop start" << endl;

    CplHop cpl_hop;
    cpl_hop = CplHop( cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, 
		      nullptr, nullptr, nullptr,
		      {}, cCplHopCutoff, cCplHopMaxExtOrder );

    SupCell sc;
    sc = SupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, zerovec, zerovec );

    sc.setMolComCoor(vec{0,0,5});
    cout << cpl_hop.isComplete() << endl;
    cpl_hop.reset(&sc.lat_coor, &sc.mol_coor, &sc.mol.mass(), sup_lat_vec_list);
    cout << cpl_hop.isComplete() << endl;
    cout << cpl_hop.ampl(0, zerovec) << endl;
    cout << cpl_hop.diff(0, 2, zerovec) << endl;
    cout << cpl_hop.diff(0, 0, 2, zerovec) << endl;
    cout << cpl_hop.diff(1, 0, 2, zerovec) << endl;
    cout << cpl_hop.H(zerovec) << endl;
    cout << cpl_hop.numLatSites() << endl;
    cout << cpl_hop.numMolSites() << endl;
    cout << cpl_hop.numLatBaseSites() << endl;
    cout << cpl_hop.numLatUnitCells() << endl;
    cout << cpl_hop.numLatOrbs() << endl;
    cout << cpl_hop.numLatBasisOrbs() << endl;
    cout << cpl_hop.subLatIdx(5) << endl;
    cout << cpl_hop.numLatOrbs(3) << endl;
    cout << cpl_hop.idxLatOrbs(2) << endl;
    cout << cpl_hop.idxLatBasisOrbs(4) << endl;
    cout << cpl_hop.cpl_orb_neighbor.lat(0) << endl;
    cout << cpl_hop.cpl_orb_neighbor.mol() << endl;

    cout << "CplHop end" << endl << endl << endl;;
}


///////////////////////////////////////////////////////////////
//
//	    LatPot
//
///////////////////////////////////////////////////////////////

if (LatPot_Switch) {
    cout << "LatPot start" << endl;

    SupCell sc;
    sc = SupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, zerovec, zerovec );

    LatPot lat_pot;
    lat_pot = LatPot( cLatPotSpringConst, cLatPotEqDist,
		      nullptr, mat{}, {}, 1 );

    cout << lat_pot.isComplete() << endl;
    lat_pot.reset( &sc.lat_coor, sc.lat_sup_cell.coor(), sup_lat_vec_list );
    cout << lat_pot.isComplete() << endl;

    cout << lat_pot.energy(0, 1) << endl;
    cout << lat_pot.force() << endl;
    cout << lat_pot.energy() << endl;
    cout << lat_pot.neighbor(0) << endl;

    cout << "LatPot end" << endl << endl << endl;;
}

///////////////////////////////////////////////////////////////
//
//	    NO_Au_Diabats
//
///////////////////////////////////////////////////////////////

if (NO_Au_Diabats_Switch) {
    cout << "NO_Au_Diabats start" << endl;

    NO_Au_Diabats diab(
	    cNeutralMorseCoef, cNeutralMorseEqDist, cNeutralMorseExpDecayLen, 
	    cNeutralAuNCoef, cNeutralAuOCoef,
	    cNeutralAuNExpDecayLen, cNeutralAuOExpDecayLen, 
	    cIonicMorseCoef, cIonicMorseEqdist, cIonicMorseExpDecayLen, 
	    cIonicAuNCoef, cIonicAuOCoef,
	    cIonicAuNExpDecayLen, cIonicAuOExpDecayLen, cIonicAuNEqDist, 
	    cImageCoef, cImageRegLength, cImageRefZ, cWorkFunc, cElecAffinity, 
	    nullptr, nullptr, nullptr, nullptr, 
	    {}, cCplPotCutoff, 1
    );

    SupCell sc;
    sc = SupCell( cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, zerovec, zerovec );

    cout << diab.isComplete() << endl;
    diab.reset( &sc.lat_coor, &sc.mol_coor, &sc.mol_velo, &sc.mol.mass(),
		sup_lat_vec_list );
    sc.setMolComCoor(vec{12, 9, 6});
    sc.setMolComVelo(vec{0, 0, -0.1});
    cout << diab.isComplete() << endl;
    cout << diab.neutralVibOmega() << endl;
    cout << diab.ionicVibOmega() << endl;
    cout << diab.neutralMomInertia() << endl;
    cout << diab.ionicMomInertia() << endl;
    cout << diab.neutralMorseEigEnergy(0) << endl;
    cout << diab.ionicMorseEigEnergy(0) << endl;
    cout << diab.numVibQuanta() << endl;
    cout << diab.numRotQuanta() << endl;

    arma_rng::set_seed_random();
    diab.randSetMolOrientation();
    diab.setMolVibQuanta(5);
    diab.setMolRotQuanta(3);
    cout << sc.molComCoor() << endl;
    cout << sc.molComVelo() << endl;
    cout << sc.mol_coor << endl;
    cout << sc.mol_velo << endl;
    cout << diab.numVibQuanta() << endl;
    cout << diab.numRotQuanta() << endl;

    cout << "NO_Au_Diabats end" << endl << endl << endl;;
}



///////////////////////////////////////////////////////////////
//
//	    BrilZone
//
///////////////////////////////////////////////////////////////

if (BrilZone_Switch) {
    cout << "BrilZone start" << endl;

    BrilZone bz;
    BravLat sup_brav(sup_lat_vec_list);
    bz = BrilZone(sup_brav, cNumK);

    cout << bz.dim() << endl;
    cout << bz.numk() << endl;
    for (auto& sv : sup_lat_vec_list)
	sv.print();
    cout << bz.rcpVec() << endl;
    cout << dot(sup_lat_vec_list[0], bz.rcpVec(1)) << endl;
    cout << dot(sup_lat_vec_list[1], bz.rcpVec(0)) << endl;
    cout << bz.k(2) << endl;
    cout << bz.k() << endl;
    bz.k().save("test_data/k.txt", raw_ascii);

    cout << "BrilZone end" << endl << endl << endl;;
}    



///////////////////////////////////////////////////////////////
//
//	    System
//
///////////////////////////////////////////////////////////////

if (System_Switch) {
    cout << "System start" << endl;

    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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

    sys.setMolComCoor(vec{0,0,5});

    sys.bril_zone.k().print(); cout << endl;
    sys.sup_brav.latVec().print(); cout << endl;
    sys.lat_H.H(zerovec).col(0).print(); cout << endl;
    sys.lat_pot.force().print(); cout << endl;
    sys.H(0).print();
    sys.H(10).print();

    cout << sys.diabats.numVibQuanta(0) << endl;

    cout << "System end" << endl << endl << endl;;
}

// {13.75, 4.7631, 15}; // on-top (on site 6)
// {13.75, 7.9386, 15}; // empty hollow (among site 6, 9, 10)
// {11.00, 6.3509, 15}; // full hollow (among site 5, 6, 9)
// {12.375, 7.1447, 15}; // bridge (between site 6, 9)


///////////////////////////////////////////////////////////////
//
//	    E_Gamma and E_dos
//
///////////////////////////////////////////////////////////////

if (E_dep_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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
    mat E_Gamma = zeros(100, 4);
    mat E_dos = zeros(100, 2);
    E_Gamma.col(0) = linspace(sys.latBandMinE(), sys.latBandMaxE());
    E_dos.col(0) = linspace(sys.latBandMinE(), sys.latBandMaxE());
    double z = 4;

    // on-top
    sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z});
    sys.updCplNeighbor();
    E_Gamma.col(1) = sys.Gamma();
    E_dos.col(1) = sys.latBandDos();

    // bridge
    sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z});
    sys.updCplNeighbor();
    E_Gamma.col(2) = sys.Gamma();

    // hollow
    sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3+vec{0, 0, z});
    sys.updCplNeighbor();
    E_Gamma.col(3) = sys.Gamma();

    string E_Gamma_file_name = "test_data/E_Gamma_z";
    string E_dos_file_name = "test_data/E_dos_z";
    stringstream sstream;
    sstream.str("");
    sstream << z;
    E_Gamma_file_name += sstream.str() + "_";
    E_dos_file_name += sstream.str() + "_";
    sstream.str("");
    sstream << sys.bril_zone.numk();
    E_Gamma_file_name += "nk" + sstream.str() + ".txt";
    E_dos_file_name += "nk" + sstream.str() + ".txt";
    E_Gamma.save(E_Gamma_file_name, raw_ascii);
    E_dos.save(E_dos_file_name, raw_ascii);

    cout << "E_Gamma and E_DOS finished" << endl;

}


///////////////////////////////////////////////////////////////
//
//	    z_Gamma
//
///////////////////////////////////////////////////////////////

if (Z_dep_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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
    int n_grids = 50;
    vec z = linspace(3, 10, n_grids);
    mat z_Gamma = zeros(n_grids, 4);
    z_Gamma.col(0) = z;
    for (int i = 0; i < n_grids; ++i) {
	// on-top
	sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z(i)});
	sys.updCplNeighbor();
	z_Gamma(i, 1) = sys.Gamma_h();

	// bridge
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	z_Gamma(i, 2) = sys.Gamma_h();

	// hollow
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3
			    +vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	z_Gamma(i, 3) = sys.Gamma_h();
    }
    string z_Gamma_file_name = "test_data/z_Gamma_nk";
    stringstream sstream;
    sstream.str("");
    sstream << sys.bril_zone.numk();
    z_Gamma_file_name += sstream.str() + ".txt";
    z_Gamma.save(z_Gamma_file_name, raw_ascii);

    cout << "z_Gamma completed" << endl;
}


///////////////////////////////////////////////////////////////
//
//	    Two_Diabats
//
///////////////////////////////////////////////////////////////

if (Two_Diabats_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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
    int n_grids = 50;
    vec z = linspace(3, 10, n_grids);
    mat diabats = zeros(n_grids, 7);
    diabats.col(0) = z;
    double bond_length = sys.diabats.neutralMaxBondLength(10);
    sys.setMolBondLength(bond_length);
    //sys.computeLatBand();
    double latE = 0;

    auto entropy = [] (vec& v) {
	double S = 0;
	v.for_each( [&S] (double& val) {
		if ( (val>EPS) && ((1-val)>EPS) )
		    S -= val*log(val) + (1-val)*log(1-val);
	} );
	return S;
    };


    for (int i = 0; i < n_grids; ++i) {
	// on-top
	sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z(i)});
	sys.updCplNeighbor();
	diabats(i, 1) = sys.diabats.neutral();
	diabats(i, 4) = sys.diabats.ionic();

	// bridge
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	diabats(i, 2) = sys.diabats.neutral();
	diabats(i, 5) = sys.diabats.ionic();

	// hollow
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3
			    +vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	diabats(i, 3) = sys.diabats.neutral();
	diabats(i, 6) = sys.diabats.ionic();
    }

    vec occ;
    for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	occ = fermi(sys.latBand(k_idx), cThermBeta, cChemPot);
	latE += accu(sys.latBand(k_idx) % occ) - entropy(occ) / cThermBeta - cChemPot*accu(occ);
    }
    latE /= sys.bril_zone.numk();

    diabats.tail_cols(6) += latE;
    //cout << bond_length << endl;
    string diabats_file_name = "test_data/diabats.txt";
    diabats.save(diabats_file_name, raw_ascii);

    cout << "Two diabats finished" << endl;
}


///////////////////////////////////////////////////////////////
//
//	    Check_Diabats
//
///////////////////////////////////////////////////////////////

if (Check_Diabats_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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

    bool broaden_switch = 1;
    double dt = 1;
    uword max_num_steps = 10000;
    double terminal_z = 15.001;
    bool print_switch = 1;
    CME< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > cme(
	    &sys, cThermBeta, cChemPot, broaden_switch, cLatFreeze, cLatZFreeze,
	    dt, max_num_steps, terminal_z, print_switch );

    int n_grids = 50;
    vec z = linspace(3, 10, n_grids);
    double dz = z(1) - z(0);
    mat diabats_var = zeros(n_grids, 7); // constructed by integrating forces
    mat f_z_N = zeros(n_grids, 6); //
    mat f_z_O = zeros(n_grids, 6); //
    diabats_var.col(0) = z;
    double bond_length = sys.diabats.neutralMaxBondLength(10);
    sys.setMolBondLength(bond_length);

    double latE = 0;

    auto entropy = [] (vec& v) {
	double S = 0;
	v.for_each( [&S] (double& val) {
		if ( (val>EPS) && ((1-val)>EPS) )
		    S -= val*log(val) + (1-val)*log(1-val);
	} );
	return S;
    };

    rowvec ref = zeros<rowvec>(6);
    mat f = zeros<mat>(3, 2);
    for (int i = 0; i < n_grids; ++i) {
	// on-top
	sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z(i)});
	sys.updCplNeighbor();

	f = sys.diabats.neutralForceMol();
	f_z_N(i, 0) = f.col(0)(2);
	f_z_O(i, 0) = f.col(1)(2);

	f = sys.diabats.ionicForceMol();
	f_z_N(i, 3) = f.col(0)(2);
	f_z_O(i, 3) = f.col(1)(2);

	if ( i == (n_grids-1) ) {
	    ref(0) = sys.diabats.neutral();
	    ref(3) = sys.diabats.ionic();
	}

	// bridge
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z(i)});
    	sys.updCplNeighbor();

	f = sys.diabats.neutralForceMol();
	f_z_N(i, 1) = f.col(0)(2);
	f_z_O(i, 1) = f.col(1)(2);

	f = sys.diabats.ionicForceMol();
	f_z_N(i, 4) = f.col(0)(2);
	f_z_O(i, 4) = f.col(1)(2);

	if ( i == (n_grids-1) ) {
	    ref(1) = sys.diabats.neutral();
	    ref(4) = sys.diabats.ionic();
	}

	// hollow
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3
			    +vec{0, 0, z(i)});
    	sys.updCplNeighbor();

	f = sys.diabats.neutralForceMol();
	f_z_N(i, 2) = f.col(0)(2);
	f_z_O(i, 2) = f.col(1)(2);

	f = sys.diabats.ionicForceMol();
	f_z_N(i, 5) = f.col(0)(2);
	f_z_O(i, 5) = f.col(1)(2);

	if ( i == (n_grids-1) ) {
	    ref(2) = sys.diabats.neutral();
	    ref(5) = sys.diabats.ionic();
	}
    }

    vec occ;
    for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	occ = fermi(sys.latBand(k_idx), cThermBeta, cChemPot);
	latE += accu(sys.latBand(k_idx) % occ) - entropy(occ) / cThermBeta - cChemPot*accu(occ);
    }
    latE /= sys.bril_zone.numk();

    diabats_var.tail_cols(6) += latE;
    diabats_var.tail_cols(6).each_row() += ref;

    for (int i = 0; i < n_grids; ++i) {
	diabats_var(i, 1) += dz * ( accu(f_z_N.col(0).tail(n_grids-i)) +
					  accu(f_z_O.col(0).tail(n_grids-i)) );
	diabats_var(i, 2) += dz * ( accu(f_z_N.col(1).tail(n_grids-i)) +
					  accu(f_z_O.col(1).tail(n_grids-i)) );
	diabats_var(i, 3) += dz * ( accu(f_z_N.col(2).tail(n_grids-i)) +
					  accu(f_z_O.col(2).tail(n_grids-i)) );
	diabats_var(i, 4) += dz * ( accu(f_z_N.col(3).tail(n_grids-i)) +
					  accu(f_z_O.col(3).tail(n_grids-i)) );
	diabats_var(i, 5) += dz * ( accu(f_z_N.col(4).tail(n_grids-i)) +
					  accu(f_z_O.col(4).tail(n_grids-i)) );
	diabats_var(i, 6) += dz * ( accu(f_z_N.col(5).tail(n_grids-i)) +
					  accu(f_z_O.col(5).tail(n_grids-i)) );
    }

    string file_name = "test_data/diabats_var.txt";
    diabats_var.save(file_name, raw_ascii);

}

///////////////////////////////////////////////////////////////
//
//	    Broadened_Diabats
//
///////////////////////////////////////////////////////////////

if (Broadened_Diabats_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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

    bool broaden_switch = 1;
    double dt = 1;
    uword max_num_steps = 10000;
    double terminal_z = 15.001;
    bool print_switch = 1;
    CME< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > cme(
	    &sys, cThermBeta, cChemPot, broaden_switch, cLatFreeze, cLatZFreeze,
	    dt, max_num_steps, terminal_z, print_switch );

    int n_grids = 50;
    vec z = linspace(3, 10, n_grids);
    double dz = z(1) - z(0);
    mat broadened_diabats = zeros(n_grids, 7);
    mat f_corr_z_N = zeros(n_grids, 6); // neutral and ionic BCME correction forces for N
    mat f_corr_z_O = zeros(n_grids, 6); // for O
    broadened_diabats.col(0) = z;
    double bond_length = sys.diabats.neutralMaxBondLength(10);
    sys.setMolBondLength(bond_length);

    double latE = 0;

    auto entropy = [] (vec& v) {
	double S = 0;
	v.for_each( [&S] (double& val) {
		if ( (val>EPS) && ((1-val)>EPS) )
		    S -= val*log(val) + (1-val)*log(1-val);
	} );
	return S;
    };

    mat f_corr = zeros<mat>(3, 2);
    for (int i = 0; i < n_grids; ++i) {
	// on-top
	sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z(i)});
	sys.updCplNeighbor();
	broadened_diabats(i, 1) = sys.diabats.neutral();
	broadened_diabats(i, 4) = sys.diabats.ionic();

	cme.state = 0;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 0) = f_corr.col(0)(2);
	f_corr_z_O(i, 0) = f_corr.col(1)(2);

	cme.state = 1;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 3) = f_corr.col(0)(2);
	f_corr_z_O(i, 3) = f_corr.col(1)(2);

	// bridge
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	broadened_diabats(i, 2) = sys.diabats.neutral();
	broadened_diabats(i, 5) = sys.diabats.ionic();

	cme.state = 0;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 1) = f_corr.col(0)(2);
	f_corr_z_O(i, 1) = f_corr.col(1)(2);

	cme.state = 1;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 4) = f_corr.col(0)(2);
	f_corr_z_O(i, 4) = f_corr.col(1)(2);

	// hollow
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3
			    +vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	broadened_diabats(i, 3) = sys.diabats.neutral();
	broadened_diabats(i, 6) = sys.diabats.ionic();

	cme.state = 0;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 2) = f_corr.col(0)(2);
	f_corr_z_O(i, 2) = f_corr.col(1)(2);

	cme.state = 1;
	f_corr = cme.corrForce();
	f_corr_z_N(i, 5) = f_corr.col(0)(2);
	f_corr_z_O(i, 5) = f_corr.col(1)(2);
    }

    vec occ;
    for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	occ = fermi(sys.latBand(k_idx), cThermBeta, cChemPot);
	latE += accu(sys.latBand(k_idx) % occ) - entropy(occ) / cThermBeta - cChemPot*accu(occ);
    }
    latE /= sys.bril_zone.numk();

    broadened_diabats.tail_cols(6) += latE;

    for (int i = 0; i < n_grids; ++i) {
	broadened_diabats(i, 1) += dz * ( accu(f_corr_z_N.col(0).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(0).tail(n_grids-i)) );
	broadened_diabats(i, 2) += dz * ( accu(f_corr_z_N.col(1).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(1).tail(n_grids-i)) );
	broadened_diabats(i, 3) += dz * ( accu(f_corr_z_N.col(2).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(2).tail(n_grids-i)) );
	broadened_diabats(i, 4) += dz * ( accu(f_corr_z_N.col(3).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(3).tail(n_grids-i)) );
	broadened_diabats(i, 5) += dz * ( accu(f_corr_z_N.col(4).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(4).tail(n_grids-i)) );
	broadened_diabats(i, 6) += dz * ( accu(f_corr_z_N.col(5).tail(n_grids-i)) +
					  accu(f_corr_z_O.col(5).tail(n_grids-i)) );
    }

    string file_name = "test_data/broadened_diabats.txt";
    broadened_diabats.save(file_name, raw_ascii);
    
    cout << "broadened diabats finished" << endl;
}


///////////////////////////////////////////////////////////////
//
//	    adiabats
//
///////////////////////////////////////////////////////////////

if (Adiabats_Switch) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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

    int n_grids = 50;
    vec z = linspace(3, 10, n_grids);
    mat adiabats = zeros(n_grids, 4);
    adiabats.col(0) = z;
    double bond_length = sys.diabats.neutralMaxBondLength(10);
    sys.setMolBondLength(bond_length);

    auto entropy = [] (vec& v) {
	double S = 0;
	v.for_each( [&S] (double& val) {
		if ( (val>EPS) && ((1-val)>EPS) )
		    S -= val*log(val) + (1-val)*log(1-val);
	} );
	return S;
    };
	
    cx_mat eigvec;
    vec eigval;
    vec occ;
    double therm_dyn_pot = 0;
    
    for (int i = 0; i < n_grids; ++i) {
	// on-top
	sys.setMolComCoor(sys.latCoor(6)+vec{0, 0, z(i)});
	sys.updCplNeighbor();
	therm_dyn_pot = 0;
	for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	    eigval = eig_sym( sys.H(k_idx) );
	    occ = fermi(eigval, cThermBeta, cChemPot);
	    therm_dyn_pot += accu(occ%eigval) - entropy(occ) / cThermBeta- cChemPot*accu(occ);
	}
	adiabats(i, 1) = therm_dyn_pot/sys.bril_zone.numk() + sys.diabats.neutral();

	// bridge
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9))/2 + vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	therm_dyn_pot = 0;
	for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	    eigval = eig_sym( sys.H(k_idx) );
	    occ = fermi(eigval, cThermBeta, cChemPot);
	    therm_dyn_pot += accu(occ%eigval) - entropy(occ) / cThermBeta- cChemPot*accu(occ);
	}
	adiabats(i, 2) = therm_dyn_pot/sys.bril_zone.numk() + sys.diabats.neutral();

	// hollow
	sys.setMolComCoor((sys.latCoor(6)+sys.latCoor(9)+sys.latCoor(10))/3
			    +vec{0, 0, z(i)});
    	sys.updCplNeighbor();
	therm_dyn_pot = 0;
	for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
	    eigval = eig_sym( sys.H(k_idx) );
	    occ = fermi(eigval, cThermBeta, cChemPot);
	    therm_dyn_pot += accu(occ%eigval) - entropy(occ) / cThermBeta- cChemPot*accu(occ);
	}
	adiabats(i, 3) = therm_dyn_pot/sys.bril_zone.numk() + sys.diabats.neutral();
    }
    string adiabats_file_name = "test_data/adiabats.txt";
    adiabats.save(adiabats_file_name, raw_ascii);

    cout << "adiabats finished" << endl;
}



///////////////////////////////////////////////////////////////
//
//	    CME
//
///////////////////////////////////////////////////////////////

if (CME_Switch) {
    bool broaden_switch = 1;
    double dt = 1;
    uword max_num_steps = 10000;
    double terminal_z = 15.001;
    bool print_switch = 1;

    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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

    //sys.randSetMolComCoorVelo(15, vec{12, 9, 0}, 0.01, 0);
    sys.setMolComCoorVelo(15, vec{13.75, 7.9386, 0}, 0.01, 0, 0);
    sys.diabats.randSetMolOrientation();
    sys.diabats.setMolVibQuanta(15);
    sys.updCplNeighbor();
    CME< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > cme(
	    &sys, cThermBeta, cChemPot, broaden_switch, cLatFreeze, cLatZFreeze,
	    dt, max_num_steps, terminal_z, print_switch );
    cout << sys.molComCoor() << endl;
    cout << sys.molComVelo() << endl;
    cout << sys.mol_coor << endl;
    cout << sys.mol_velo << endl;
    cout << sys.diabats.numRotQuanta() << endl;
    cout << sys.diabats.numVibQuanta() << endl;
    //cme.propagate();


}



///////////////////////////////////////////////////////////////
//
//	    EFLD
//
///////////////////////////////////////////////////////////////

if (EFLD_Switch) {
    bool fric_switch = 1;
    double dt = 1;
    uword max_num_steps = 10000;
    double terminal_z = 15.001;
    bool print_switch = 1;
    //vec all_width = {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
    vec all_width = regspace(0.004, 0.002, 0.08);
    vec fric;
    cout << "numUnitCells = " << endl << cNumUnitCell << endl;
    cout << "numk = " << endl << cNumK << endl;
for (auto& width : all_width) {
    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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
	    width);


    //sys.setMolComCoorVelo(15, vec{13.75, 7.9386, 0}, 0.01, 0, 0);
    //sys.diabats.randSetMolOrientation();
    //sys.diabats.setMolVibQuanta(15);
    //mat mol_coor_x;
    //mat mol_coor_y;
    //mat mol_coor_z;
    //mol_coor_x.load("test_data/onetraj2/cme_mol_coor_x.txt");
    //mol_coor_y.load("test_data/onetraj2/cme_mol_coor_y.txt");
    //mol_coor_z.load("test_data/onetraj2/cme_mol_coor_z.txt");
    mat mol_coor = { {13.75, 13.75}, {7.9386, 7.9386}, {3.53, 6.53} };
    //mat mol_coor = { {0, 0}, {0, 0}, {3.53, 6.53} };
    //sys.computeLatBand();
    //cout << "bond length = " << norm(mol_coor.col(0) - mol_coor.col(1) ) << endl;
    //cout << "band max E = " << sys.latBandMaxE() << endl;
    //cout << "band min E = " << sys.latBandMinE() << endl;
    //mat mol_coor = zeros(3,2);
    //mol_coor.row(0) = mol_coor_x.row(13823);
    //mol_coor.row(1) = mol_coor_y.row(13823);
    //mol_coor.row(2) = mol_coor_z.row(13823);

    sys.setMolCoor(mol_coor);
    //arma_rng::set_seed_random();
    arma_rng::set_seed(1);
    sys.lat_coor += 1*randu(size(sys.lat_coor));
    sys.updCplNeighbor();

    //cout << sys.molComCoor() << endl; cout << endl;
    //sys.H(0).tail_cols(1).print(); cout << endl;
    //vec eigval = eig_sym(sys.H(0));
    //eigval.print(); cout << endl;

    // check if h and dh match
    //double h0 = sys.diabats.h();
    //double dh0 = sys.diabats.dhMol(0, 2);
    //sys.mol_coor(2, 0) += 0.01;
    //double h1 = sys.diabats.h();
    //double dh1 = sys.diabats.dhMol(0, 2);
    //cout << (h1-h0) / 0.01 << endl;
    //cout << (dh0+dh1)/2 << endl;
    //assert(false);

    EFLD< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > efld(
	    &sys, cThermBeta, cChemPot, fric_switch, true, cLatZFreeze,
	    dt, max_num_steps, terminal_z, print_switch );

    //efld.update();
    cout << efld.mol_rand_force_corr[5](0,0) << endl;
    fric.insert_rows(fric.n_rows, vec{efld.mol_rand_force_corr[5](0,0)});
} // end of level broadening test for-loop
    mat fric_bd = zeros(fric.n_rows, 2);
    fric_bd.col(0) = all_width;
    fric_bd.col(1) = fric;
    string data_file_name = "test_data/fric_bd";
    stringstream ss;
    ss.str("");
    ss << prod(cNumK);
    data_file_name += ss.str() + ".txt";
    fric_bd.save(data_file_name, raw_ascii);
    cout << "EFLD finished" << endl;
}

///////////////////////////////////////////////////////////////
//
//	    CME_Dynamics
//
///////////////////////////////////////////////////////////////

if (CME_Dynamics_Switch) {
    double dt = 20;
    double height = 15;
    double terminal_z = 15.001;

    System < NO_Au_Diabats, LatTB, CplHop, LatPot > sys(
	    cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
	    cMolCoor, cMolMassList, zerovec, zerovec, cIsPeriodic, cNumK,
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


    sys.setMolComCoorVelo(height, vec{13.75, 7.9386, 0}, 0.02, 0, 0);
    sys.diabats.randSetMolOrientation();
    sys.setMolBondLength(2.5);
    sys.diabats.setMolVibQuanta(15);
    cout << "init vib quanta = " << sys.diabats.numVibQuanta() << endl;

    cout << "com coor: " << endl;
    sys.molComCoor().print(); cout << endl;
    cout << "com velo: " << endl;
    sys.molComVelo().print(); cout << endl;
    cout << "mol coor" << endl;
    sys.mol_coor.raw_print(std::cout); cout << endl;
    cout << "mol velo" << endl;
    sys.mol_velo.raw_print(std::cout); cout << endl;
    cout << "mol mass mat" << endl;
    sys.molMassMat().print();
    cout << "bond_length = " << sys.molBondLength() << endl;

    sys.lat_velo.zeros();
    sys.updCplNeighbor();

    mat mol_force(size(sys.mol_coor));
    mol_force = sys.diabats.neutralForceMol();

    mat mol_force_dt(size(sys.mol_coor));
    double curr_time = 0;
    std::vector<double> energy;
    std::vector<double> time;
    energy.push_back(sys.diabats.neutral() + sys.molKinE());
    time.push_back(curr_time);
    
    double min_z = height;
    while ( sys.molComCoor()(2) < terminal_z ) {
	sys.mol_coor += dt * (sys.mol_velo + 0.5 * dt * mol_force / sys.molMassMat() );
	sys.updCplNeighbor();
	mol_force_dt = sys.diabats.neutralForceMol();
	sys.mol_velo += 0.5 * dt * (mol_force + mol_force_dt) / sys.molMassMat();
	mol_force = mol_force_dt;
	curr_time += dt;
	energy.push_back(sys.diabats.neutral() + sys.molKinE());
	time.push_back(curr_time);
	if (sys.molComCoor()(2) < min_z)
	    min_z= sys.molComCoor()(2);
    }
    cout << "min_z = " << min_z << endl;

    stringstream sstream;
    sstream << dt;
    conv_to<vec>::from(energy).save("test_data/energy_dt" + sstream.str() + ".txt", raw_ascii);
    conv_to<vec>::from(time).save("test_data/time_dt" + sstream.str() + ".txt", raw_ascii);
}


///////////////////////////////////////////////////////////////
//
//	    BO_Dynamics
//
///////////////////////////////////////////////////////////////

if (BO_Dynamics_Switch) {
    const bool fric_switch = false;
    const double dt = 2.0;
    const uword max_num_steps = 20000;
    const double terminal_com_z = 15.001;
    const double mol_init_com_z = 15;
    const bool print_switch = true;
    const bool conv_test = true;

    const int incident = 0;
    const double mol_init_trans_kinE = 0.02;
    const double mol_init_vib_quanta = 15.0;
    const double mol_init_rot_quanta = 0;

    arma_rng::set_seed(1);
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
    
    EFLD< System< NO_Au_Diabats, LatTB, CplHop, LatPot> > efld(
	    &sys, cThermBeta, cChemPot, fric_switch, cLatFreeze, cLatZFreeze,
	    dt, max_num_steps, terminal_com_z, print_switch );

    efld.propagate();

    string dir = "test_data/onetraj/";
    extract(efld.lat_coor_t, 0).save( dir+"efld_lat_coor_x.txt", raw_ascii );
    extract(efld.lat_coor_t, 1).save( dir+"efld_lat_coor_y.txt", raw_ascii );
    extract(efld.lat_coor_t, 2).save( dir+"efld_lat_coor_z.txt", raw_ascii );
    extract(efld.mol_coor_t, 0).save( dir+"efld_mol_coor_x.txt", raw_ascii );
    extract(efld.mol_coor_t, 1).save( dir+"efld_mol_coor_y.txt", raw_ascii );
    extract(efld.mol_coor_t, 2).save( dir+"efld_mol_coor_z.txt", raw_ascii );

    if (!conv_test) {
	conv_to<vec>::from(efld.totE_t).save( dir+"efld_totE.txt", raw_ascii );
        conv_to<vec>::from(efld.run_time).save( dir+"efld_time.txt", raw_ascii);
    } else {
        conv_to<vec>::from(efld.totE_t).save( dir+"efld_totE_test.txt", raw_ascii );
        conv_to<vec>::from(efld.run_time).save( dir+"efld_time_test.txt", raw_ascii);
    }
}

if (Sample_Traj_Switch) {
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
    // read data
    mat mol_coor_x;    
    mat mol_coor_y;    
    mat mol_coor_z;    


//    mol_coor_x.load("");
//    mol_coor_y.load("");
//    mol_coor_z.load("");

    uword sz = mol_coor_x.n_rows;

    // diabats along the trajectory
    mat neutral = zeros(sz, 1);
    mat ionic = zeros(sz, 1);
    mat Gamma = zeros(sz, 1);


    mat mol_coor = zeros(3, 2);
    for (uword i = 0; i < sz; ++i) {
	mol_coor.row(0) = mol_coor_x.row(i);
	mol_coor.row(1) = mol_coor_y.row(i);
	mol_coor.row(2) = mol_coor_z.row(i);
	sys.setMolCoor(mol_coor);
	sys.updCplNeighbor();

	neutral(i, 0) = sys.diabats.neutral();
	ionic(i, 0) = sys.diabats.ionic();
	Gamma(i, 0) = sys.Gamma_h();
    }
    neutral.save("", raw_ascii);
    ionic.save("", raw_ascii);
    Gamma.save("", raw_ascii);


}

    cout << "end test" << endl;
    main_dur = myclock::now() - main_start; 
    cout << "time elapsed = " << main_dur.count() << endl;
    return 0;
}
