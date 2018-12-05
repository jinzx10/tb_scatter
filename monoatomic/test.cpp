#include <iostream>
#include <vector>
#include <string>
#include <armadillo>
#include <complex>
#include <cmath>
//#include "bcme_all.h"
#include "efld_all.h"
#include <chrono>

using namespace arma;
using namespace std;
//using namespace test;
using namespace hexagonal;
using namespace monoatomic;
using myclock = std::chrono::high_resolution_clock;

int main()
{
    myclock::time_point main_start = myclock::now();
    chrono::duration<double> main_dur;

    bool timing_switch = true;
    bool E_gamma_dos_switch = true;
    bool z_on_site_E_switch = true;
    bool z_gamma_switch = true;
    bool band_struct_switch = false;
    bool cpl_pot_switch = true;

    double z = 5.0; // for E_gamma_dos

    uword num_z = 100; // for z-dependence
    vec vz = linspace(3.4, 8, num_z);

    double height = 5.0; // for cpl_pot plot

    arma_rng::set_seed_random();

    EFLD efld_test(
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
            cTerminalZ, cPrintSwitch 
    );

//////////////////////////////////////////////
//
//
//	E dependence: Gamma, DOS
//
//
//////////////////////////////////////////////

    if (E_gamma_dos_switch) {
    	efld_test.sys.sup_cell.coor.tail_cols(1)(2) = z;
    	efld_test.computeAndSaveInit();
    	efld_test.sys.lat_band.compute();

    	uword n_egrid = 100;
    	mat E_dep = zeros(n_egrid, 3);

    	E_dep.col(0) = linspace( efld_test.sys.lat_band.minEnergy(), 
    	    		    efld_test.sys.lat_band.maxEnergy(), 100);
    	E_dep.col(1) = efld_test.sys.Gamma();
    	E_dep.col(2) = efld_test.sys.lat_band.dos();
    	//cout << "on site h: " << efld_test.sys.mol_on_site.h() << endl
    	//     << "bandmin: " << efld_test.sys.bandMin() << endl
    	//     << "bandmax: " << efld_test.sys.bandMax() << endl;
    	vec k0 = {0,0,0};
    	//efld_test.sys.Vk(0).print(); cout << endl;
    	//efld_test.sys.lat_band.energy(0).print();
    	//cout << z << "  " << efld_test.sys.Gamma() << "  " << endl;

    	E_dep.save("E_dep.txt", raw_ascii);
    }


//////////////////////////////////////////////////////////////////
//
//
//	on-site energy vs. z
//
//
/////////////////////////////////////////////////////////////////

    if (z_on_site_E_switch) {
    	mat z_dep = zeros(num_z, 9);
    	z_dep.col(0) = vz;
    	for (uword idx = 0; idx < num_z; ++idx) {
    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{13.75, 4.7631, 18}; // on-top
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_dep(idx, 1) = efld_test.sys.pot.cpl_pot.energy(48);
    	    z_dep(idx, 2) = efld_test.sys.pot.cpl_pot.energy(48) + efld_test.sys.mol_on_site.h();

    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{13.75, 7.9386, 20}; // empty hollow
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_dep(idx, 3) = efld_test.sys.pot.cpl_pot.energy(48);
    	    z_dep(idx, 4) = efld_test.sys.pot.cpl_pot.energy(48) + efld_test.sys.mol_on_site.h();
    	    
    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{12.375, 7.1447, 10}; // bridge
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_dep(idx, 5) = efld_test.sys.pot.cpl_pot.energy(48);
    	    z_dep(idx, 6) = efld_test.sys.pot.cpl_pot.energy(48) + efld_test.sys.mol_on_site.h();

    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{11.00, 6.3509, 20}; // full hollow
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_dep(idx, 7) = efld_test.sys.pot.cpl_pot.energy(48);
    	    z_dep(idx, 8) = efld_test.sys.pot.cpl_pot.energy(48) + efld_test.sys.mol_on_site.h();
    	}

    	z_dep.save("z_dep.txt", raw_ascii);
    }


//////////////////////////////////////////////////////////////////
//
//
//	Gamma0 vs. z
//
//
/////////////////////////////////////////////////////////////////

    if (z_gamma_switch) {
	mat z_gamma = zeros(num_z, 5);
    	z_gamma.col(0) = vz;
    	efld_test.sys.lat_band.compute();
    	for (uword idx = 0; idx < num_z; ++idx) {
    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{13.75, 4.7631, 18}; // on-top
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_gamma(idx, 1) = efld_test.sys.Gamma0();

    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{13.75, 7.9386, 20}; // empty hollow
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_gamma(idx, 2) = efld_test.sys.Gamma0();
    	    
    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{12.375, 7.1447, 10}; // bridge
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_gamma(idx, 3) = efld_test.sys.Gamma0();

    	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{11.00, 6.3509, 20}; // full hollow
    	    efld_test.sys.sup_cell.coor.tail_cols(1)(2) = vz(idx);
    	    z_gamma(idx, 4) = efld_test.sys.Gamma0();
    	}

    	z_gamma.save("z_gamma.txt", raw_ascii);
    }


//////////////////////////////////////////////
//
//
//	band structure
//
//
//////////////////////////////////////////////

    if (band_struct_switch) {
	uword numk = efld_test.sys.bril_zone.numk();
    	efld_test.sys.lat_band.compute();
    	efld_test.sys.bril_zone.kMat().save("k.txt", raw_ascii);
    	mat band_struct = zeros(1, numk);
    	for (uword k_idx = 0; k_idx < numk; ++k_idx) {
    	    band_struct.col(k_idx) = efld_test.sys.lat_band.energy(k_idx);
    	}
    	band_struct.save("bandE.txt", raw_ascii);
    }


//////////////////////////////////////////////
//
//
//	coupling potential plot
//
//
//////////////////////////////////////////////

    if (cpl_pot_switch) {
	uvec vn = {50, 50}; // number of points in x and y direction
	double meshsize = 20.0;
	mat meshgrid = mesh(vn) * meshsize; //
	rowvec cplpotval(meshgrid.n_cols);
	for (uword i = 0; i < meshgrid.n_cols; ++i) {
	    efld_test.sys.sup_cell.coor.tail_cols(1) = vec{ meshgrid.col(i)(0), 
							    meshgrid.col(i)(1), height };
	    cplpotval(i) = efld_test.sys.pot.cpl_pot.energy();
	}
	mat joint = join_cols(meshgrid, cplpotval);
	joint.save("cpl_pot.txt", raw_ascii);
	cout << cplpotval.max() - cplpotval.min() << endl;
    }









/*
//    arma_rng::set_seed_random();
    vec k0 = {0,0,0};

//    mat mc = joinCols<mat>({ones(3,2),2*ones(1,2),3*ones(2,2)});
//    mc.print();
//    mat nc = joinRows<mat>({ones(3,2),2*ones(3,1),3*ones(3,3)});
//    nc.print();
//    mat mn = join<mat>({ones(2,1), 2*ones(2,3), 3*ones(2,3), 4*ones(2,1)}, 2, 2, 'r');
//    mn.print();

//    GroupOfAtoms uc(cMolOnSiteE, cMolCoor, cMolMassList);
//    GroupOfAtoms uc = GroupOfAtoms(cLatOnSiteE, cLatInCellCoor, cLatMassList);
//    cout << uc.numOrbs() << endl
//	 << uc.numSites() << endl
//	 << uc.idxOrbs() << endl;
//    cout << "numOrbsOnSite: " << uc.numOrbsOnSite(0) << endl
//	 << "idxOrbsOnSite: " << uc.idxOrbsOnSite(0) << endl
//	 << "onSiteE: " << uc.onSiteE(0) << endl;
	 //<< "numOrbsOnSite: " << uc.numOrbsOnSite(1) << endl
	 //<< "idxOrbsOnSite: " << uc.idxOrbsOnSite(1) << endl << endl
	 //<< uc.mass(1) << endl
	 //<< uc.coor(1) << endl;

//    BravLat bbv;
//    cout << bbv.dim() << endl;
//    bbv = BravLat(cLatVec);
//    cout << bbv.dim() << endl;
//    for (auto& l : bbv.latVecList())
//	l.print();
//    cout << endl;
//    for (auto& l : bbv.rcpVecList())
//	l.print();
//    cout << endl;
//    bbv.latVec(1).print(); cout << endl;
//    bbv.rcpVec(0).print(); cout << endl;

    LatSupCell lsc;
    lsc = LatSupCell(cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell);
//    cout << lsc.numSites() << endl
//	 << lsc.numUnitCells() << endl
//	 << lsc.eqCoor() << endl
//	 << lsc.massMat() << endl
//	 << lsc.numBasisOrbs() <<endl
//	 << lsc.idxOrbsOnSite(22) << endl
//	 << lsc.idxOrbsOnSite(8) << endl << endl
//	 << lsc.subLatIdxOfSite(3) << endl
//	 << lsc.subLatIdxOfSite(22) << endl << endl
//	 << lsc.idxBasisOrbsOnSite(7) << endl
//	 << lsc.idxBasisOrbsOnSite(30) << endl;
	 
    Mol mol;
    mol = Mol( cMolCoor, cMolMassList, lsc.numSites(), lsc.numBasisOrbs() );
//    cout << mol.startSiteIdx() << endl
//	 << mol.endSiteIdx() << endl
//	 << mol.massMat() << endl
////	 << mol.idxOrbsOnSite(49) << endl
////	 << mol.idxOrbsOnSite(50) << endl
////	 << mol.idxBasisOrbsOnSite(48) << endl
////	 << mol.numOrbsOnSite(49) << endl
//	 << mol.numSites() << endl
//	 << mol.numOrbs() << endl;
//    mol.idxOrbs().print(); cout << endl;
//    mol.idxBasisOrbs().print(); cout << endl;

    SupCell sc;
    sc = SupCell( cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		  cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo );
//    sc.coor.print();
//    sc.velo.print();
//    cout << sc.numSites() << endl
//	 << sc.numLatSites() << endl
//	 << sc.numMolSites() << endl
//	 << sc.numBasisOrbs() << endl
//	 << sc.numLatBasisOrbs() << endl
//	 << sc.numOrbsOnSite(0) << endl
//	 << sc.numOrbsOnSite(32) << endl
//	 << sc.mol.startSiteIdx() << endl
//	 << sc.mol.endSiteIdx() << endl
//	 << sc.numOrbsOnSite(32) << endl
//	 << sc.idxOrbsOnSite(5) << endl
//	 << sc.idxOrbsOnSite(20) << endl
//	 << sc.idxOrbsOnSite(32) << endl
//	 << sc.idxBasisOrbsOnSite(9) << endl
//	 << sc.idxBasisOrbsOnSite(27) << endl
//	 << sc.idxBasisOrbsOnSite(32) << endl
//	 << sc.massMat() << endl;
//    sc.idxLatOrbs().print(); cout << endl;
//    sc.idxMolOrbs().print(); cout << endl;
//    sc.idxLatBasisOrbs().print(); cout << endl;
//    sc.idxMolBasisOrbs().print();
//    cout << "molComVelo " << endl << sc.molComVelo() << endl
//	 << "molRotE " << sc.molRotE() << endl
//	 << "molComCoor " << sc.molComCoor() << endl
//	 << "molTotMass " << sc.molTotMass() << endl
////	 << "molReducedMass " << sc.molReducedMass() << endl
//	 << "molKinE " << sc.molKinE() << endl;
//    sc.adjCoor().tail_cols(1).print(); cout << endl;

//    BrilZone bz;
//    bz = BrilZone(bbv.rcpVecList(), cNumK);
//    cout << bz.kMat() << endl
//	 << bz.numk() << endl
//	 << bz.dim() << endl
//	 << bz.k(3) << endl
//	 << bz.numkInDim(1) << endl;


    vector<vec> sup_lat_vec_list{};
    for (uword d = 0; d != cIsPeriodic.size(); ++d)
	if (cIsPeriodic[d])
	    sup_lat_vec_list.push_back( cNumUnitCell[d]*cLatVec[d] );
    BravLat spbv(sup_lat_vec_list);
//    spbv.latVec(0).print();
//    cout << endl;

    LatHop lh;
    lh = LatHop( cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, sc.lat_sup_cell.eqCoor(), 
		 sup_lat_vec_list, cLatHopMaxDistOrder, &sc.coor, &sc.lat_sup_cell );
//    lh.ampl(0, 16, vec{0,0,0}).print(); cout << endl;
//    lh.ampl(0, 1, vec{0,0,0}).print(); cout << endl;
//    cout << lh.diff(0, 1, 0, vec{0,0,0}); cout << endl;
//    cout << lh.diff(1, 0, 0, vec{0,0,0}); cout << endl;
//    cout << lh.diff(0, 1, 1, vec{0,0,0}); cout << endl;
//    lh.lat_neighbor.siteIdxNeighborOfSite(0).print();
//    auto v2b = lh.lat_neighbor.vecToBeNeighborOfSite(0);
//    for (auto& v : v2b) {
//	v.print(); cout << endl;
//    }

//    auto nb0 = lh.lat_neighbor.siteIdxNeighborOfSite(12);
//    for (uword i = 0; i < nb0.size(); ++i) {
//	cout << "site = " << nb0[i] << endl;
//	cout << "ampl: " << endl;
//	lh.ampl(12, nb0[i], k0).print(); cout << endl;
//	cout << "diff: " << endl;
//	lh.diff(12, nb0[i], 0, k0).print(); cout << endl;
//	cout << "v2b: " << endl;
//	lh.lat_neighbor.vecToBeNeighborOfSite(12)[i].print();
//    }


    CplHop ch;
    ch = CplHop( cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, &sc.coor, sup_lat_vec_list,
		 cCplHopCutoff, cCplHopMaxExtOrder, &sc );
    vector<vec> vv = ch.cpl_orb_neighbor.vecToBeNeighborOfSite(48);
    cout << vv.size() << endl;
//    ch.diff(48, 0, 0, k0).print(); cout << endl;
//    ch.diff(0, 0, k0).print(); cout << endl;
//    ch.cpl_orb_neighbor.siteIdxNeighborOfSite(0).print(); cout << endl;
//    ch.cpl_orb_neighbor.siteIdxNeighborOfSite(48).print();
//    ch.ampl(0, vec{0,0,0}).print(); cout << endl;
//    ch.diff(0, 2, vec{0,0,0}).print(); cout << endl;
//    ch.diff(48, 0, 2, vec{0,0,0}).print(); cout << endl;
//    ch.diff(49, 0, 2, vec{0,0,0}).print(); cout << endl;

    Hop hop;
    hop = Hop( cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, sup_lat_vec_list, cLatHopMaxDistOrder,
	       cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, cCplHopCutoff, cCplHopMaxExtOrder,
	       &sc );
//    hop.uniqIdxNeighborOfSite(0).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(48).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(49).print(); cout << endl;
//    hop.diff(0, 48, 2, k0).print(); cout << endl;
//    hop.diff(0, 49, 2, k0).print(); cout << endl;
//    hop.diff(48, 0, 2, k0).print(); cout << endl;
//    hop.diff(49, 0, 2, k0).print(); cout << endl;
//    hop.ampl(0, k0).print(); cout << endl;
//    hop.diff(48, 0, 0, k0).print(); cout << endl;
//    hop.diff(0, 0, k0).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(0).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(48).print(); cout << endl;
//    auto v2b = hop.cpl_hop.cpl_orb_neighbor.vecToBeNeighborOfSite(0);
//    uvec ll = hop.cpl_hop.cpl_orb_neighbor.siteIdxNeighborOfSite(48);
//    cout << "ll" << endl;
//    ll.print(); cout << endl;
//    for (auto& v : v2b) {
//	v.print(); cout << endl;
//    }

//    auto v2b48 = hop.cpl_hop.cpl_orb_neighbor.vecToBeNeighborOfSite(48);
//    for (auto& v : v2b48) {
//	v.print(); cout << endl;
//    }
//    hop.uniqIdxNeighborOfSite(0).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(49).print(); cout << endl;
//    hop.uniqIdxNeighborOfSite(48).print(); cout << endl;
//    hop.ampl(0, 1, vec{0,0,0}).print(); cout << endl;
//    hop.ampl(0, vec{0,0,0}).print(); cout << endl;
//    hop.diff(48, 0, 2, vec{0,0,0}).print(); cout << endl;
//    hop.diff(49, 0, 2, vec{0,0,0}).print(); cout << endl;
//    hop.diff(0, 49, 2, vec{0,0,0}).print(); cout << endl;
//    hop.diff(0, 48, 2, vec{0,0,0}).print(); cout << endl;


    MolOnSite mos( &sc, cMorse0Coef, cMorse0EqDist, cMorse0ExpDecayLen, 
		   cMorse1Coef, cMorse1EqDist, cMorse1ExpDecayLen, cImageC, cImageD, 
		   cRefZ, cWorkFunc, cElecAffinity );
//    cout << mos.h() << endl;
//    cout << mos.dh(48, 2) << endl;

    Hamiltonian H;
    H = Hamiltonian(&sc, &hop, &mos);
//    H(k0).col(0).print(); cout << endl;
//    vec eigval;
//    cx_mat eigvec;
//    eig_sym(eigval, eigvec, H(k0));
//    cx_vec v1 = eigvec.col(0);
//    H(k0).save("H0.txt", raw_ascii);
//    cout << as_scalar(v1.t() * H(k0) * v1) << endl;
//    eigval.print();
//    v1.print(); cout << endl;
//    H(vec{0,0,0}).print(); cout << endl;
//    cout << abs(H(vec{0,0,0})).max() << endl;
//    H.mol().print(); cout << endl; cout << endl;
//    cout << H(k)(0, 48) << endl; cout << endl;
//    cout << H(k)(48, 0) << endl; cout << endl;
//    cout << H(vec{0,0,0})(0,0) << endl
//	 << H.latOnSite()(0,0) << endl;
//    H(vec{0,0,0})(sc.idxBasisOrbsOnSite(0), sc.idxBasisOrbsOnSite(1)).print(); cout << endl;
//    H(vec{0,0,0})(sc.idxBasisOrbsOnSite(0), sc.idxBasisOrbsOnSite(48)).print(); cout << endl;


    LatPot lpot( cLatPotSpringConst, cLatPotEqDist, sc.lat_sup_cell.eqCoor(), sup_lat_vec_list,
		 cLatPotMaxDistOrder, &sc.coor, &sc.lat_sup_cell );
//    lpot.force().print(); cout << endl;
//    cout << lpot.energy(0, 1) << endl;
//    cout << lpot.force(0, 1) << endl;
//    cout << lpot.energy() << endl;
//    lpot.force().print(); cout << endl;

    MolPot mpot( &sc, cMolPotCoef, cMolPotEqDist, cMolPotExpDecayLen );
//    mpot.force().print();
//    cout << mpot.energy() << endl;
//    sc.coor.tail_cols(1) += vec{0,0,0.2};
//    cout << mpot.energy() << endl;
//    mpot.force().print(); cout << endl;

    CplPot cpot;
    cpot = CplPot(&sc, sup_lat_vec_list, cCplPotCoef, cCplPotExpDecayLen, cCplPotCutoff,
		     cCplPotMaxExtOrder );
//    cpot.force().print();

//    cpot.cpl_site_neighbor.siteIdxNeighborOfSite(0).print(); cout << endl;
//    cout << cpot.energy(0, 48) << endl;
//    cpot.force(0, 48).print(); cout << endl;

    System sys( cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo,
		cIsPeriodic, cNumK,
		cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, cLatHopMaxDistOrder,
		cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, cCplHopCutoff, cCplHopMaxExtOrder,
		cMorse0Coef, cMorse0EqDist, cMorse0ExpDecayLen,
		cMorse1Coef, cMorse1EqDist, cMorse1ExpDecayLen,
		cImageC, cImageD, cRefZ, cWorkFunc, cElecAffinity, 
		cLatPotSpringConst, cLatPotEqDist, cLatPotMaxDistOrder, 
		cMolPotCoef, cMolPotEqDist, cMolPotExpDecayLen,
		cCplPotCoef, cCplPotExpDecayLen, cCplPotCutoff, cCplPotMaxExtOrder,
		cLatFreeze, cLatZFreeze );

    EFLD efld(  cLatOnSiteE, cLatInCellCoor, cLatMassList, cLatVec, cNumUnitCell,
		cMolCoor, cMolMassList, cMolInitCenterCoor, cMolInitCenterVelo,
		cIsPeriodic, cNumK,
		cLatHopCoef, cLatHopDecayLen, cLatHopEqDist, cLatHopMaxDistOrder,
		cCplHopCoef, cCplHopDecayLen, cCplHopEqDist, cCplHopCutoff, cCplHopMaxExtOrder,
		cMorse0Coef, cMorse0EqDist, cMorse0ExpDecayLen,
		cMorse1Coef, cMorse1EqDist, cMorse1ExpDecayLen,
		cImageC, cImageD, cRefZ, cWorkFunc, cElecAffinity, 
		cLatPotSpringConst, cLatPotEqDist, cLatPotMaxDistOrder, 
		cMolPotCoef, cMolPotEqDist, cMolPotExpDecayLen,
		cCplPotCoef, cCplPotExpDecayLen, cCplPotCutoff, cCplPotMaxExtOrder,
		cLatFreeze, cLatZFreeze, cEnsemble, cThermBeta, cChemPot,
		cFricSwitch, cDt, cNumTimeStep, cFixedTerminal, cTerminalZ, cPrintSwitch );

//    cout << efld.sys.mol_on_site.h() << endl;
//    efld.elec_struc.val_list[0].print();
//    efld.sys.hop.diff(0, 1, 0, k0).print();
//    efld.elec_struc.DHDR_k[0].slice(0).print();
//    efld.sys.hop.ampl(0, k0).print();

//    sys.H(vec{0,0,0}).print(); cout << endl;
//    sys.sup_cell.coor.print();
//    sys.brav_lat.latVec(0).print();
//
//    cout << sys.brav_lat.dim() << endl
//	 << sys.brav_lat.rcpVec(1) << endl
//	 << sys.bril_zone.kMat() << endl
//	 << sys.bril_zone.numk() << endl
//	 << sys.bril_zone.numkInDim(1) << endl
//	 << sys.sup_cell.idxBasisOrbsOnSite(20) << endl
//	 << sys.sup_cell.idxBasisOrbsOnSite(50) << endl
//	 << sys.sup_cell.idxBasisOrbsOnSite(49) << endl
//	 << sys.sup_cell.massMat() << endl
//	 << sys.sup_cell.numSites() << endl
//	 << sys.sup_cell.numOrbs() << endl
//	 << sys.sup_cell.numLatSites() << endl
//	 << sys.sup_cell.numBasisOrbs() << endl
//	 << sys.sup_cell.numLatBasisOrbs() << endl
//	 << sys.sup_cell.idxLatOrbs() << endl << endl
//	 << sys.sup_cell.idxMolOrbs() << endl << endl
////	 << sys.sup_cell.idxLatBasisOrbs() << endl << endl
//	 << sys.sup_cell.idxMolBasisOrbs() << endl << endl
//	 << sys.sup_cell.idxOrbsOnSite(4) << endl
//	 << sys.sup_cell.idxOrbsOnSite(50) << endl
//	 << sys.sup_cell.idxOrbsOnSite(49) << endl
//	 << sys.sup_cell.numOrbsOnSite(2) << endl
//	 << sys.sup_cell.numOrbsOnSite(50) << endl
//	 << sys.sup_cell.numOrbsOnSite(49) << endl;

//    sys.setInit(randu(size(sys.sup_cell.coor)), randu(size(sys.sup_cell.coor)));
//    sys.sup_cell.coor.print(); cout << endl;
//    sys.hop.lat_hop.ptr_dyn_coor->print(); cout << endl;
//    sys.hop.cpl_hop.ptr_dyn_coor->print(); cout << endl;
//    sys.hop.mol_hop.ptr_dyn_coor->print(); cout << endl;
//
//
//    sys.hop.uniqIdxNeighborOfSite(0).print(); cout << endl;
//    sys.hop.uniqIdxNeighborOfSite(20).print(); cout << endl;
//    sys.hop.uniqIdxNeighborOfSite(32).print(); cout << endl;

//    sys.hop.ampl(49,50).print(); cout << endl;
//    sys.hop.ampl(49,49).print(); cout << endl;
//    sys.hop.ampl(48,50).print(); cout << endl;
//    sys.H.rawMolOnSite().print(); cout << endl;
//    sys.H.molOnSite().raw_print(); cout << endl;
//    sys.H.mol().raw_print(); cout << endl;
//    sys.H.cplLower(vec{0,0,0}).raw_print(); cout << endl;
//    sys.sup_cell.coor.cols(uvec{0, 16}).print(); cout << endl;
//    sys.hop.hop_basics.coef_mat.print(); cout << endl;
//
//    sys.H(vec{0,0,0})( sys.sup_cell.idxBasisOrbsOnSite(0),
//		       sys.sup_cell.idxBasisOrbsOnSite(20) ).print(); cout << endl;
//    sys.H(vec{0,0,0})( sys.sup_cell.idxBasisOrbsOnSite(0),
//		       sys.sup_cell.idxBasisOrbsOnSite(32) ).print(); cout << endl;
//    cout << sys.cpl_pot.energy(0, 32) << endl;
//    vec coor1 = sys.sup_cell.coor.col(0);
//    vec coor2 = sys.sup_cell.coor.col(32);
//    double r = arma::norm(coor2-coor1);
//    cout << cCplPotCoef*pow(cCplPotCharLen/r, cCplPotPowerDecay) * exp( -(r - cCplPotCharLen) / cCplPotExpDecayLen) << endl << endl;

//    mat band = zeros(sys.sup_cell.numLatBasisOrbs(), sys.bril_zone.numk());
//    for (uword k_idx = 0; k_idx < sys.bril_zone.numk(); ++k_idx) {
//	band.col(k_idx) = eig_sym(H.lat(sys.bril_zone.k(k_idx)));
//    }
//
//    band.save("band.txt", raw_ascii);
//    sys.bril_zone.kMat().save("k.txt", raw_ascii);
//
//
//    ElecStruc es(&sys, cEnsemble, cThermBeta, cChemPot, cFricSwitch, cDt);
//    vec k = {0,0,0};
//    auto hh = H(k);
//    cout << abs(hh).max() << endl;
//    vec val = eig_sym(sys.H(k));
//    val.print();
//    es.compute();


//    EFLD efld(&sys, cEnsemble, cThermBeta, cMuOrNe);
//    efld.setInit();
//    efld.ptr_sys->cpl_pot.cpl_neighbor.siteIdxNeighborOfSite(48).print();
//    cout << endl;
//    for (auto& v : efld.ptr_sys->cpl_pot.cpl_neighbor.vecToBeNeighborOfSite(48)) {
//	v.print();
//	cout << endl;
//    }
//    efld.ptr_sys->cpl_pot.cpl_neighbor.isInHomeNeighborOfSite(48).print();
//    cout << efld.ptr_sys->cpl_pot.energy(48, 0)  << endl;
//    cout << efld.ptr_sys->cpl_pot.force(48, 0)  << endl;
//    cout << endl;
//
//    cout << efld.ptr_sys->cpl_pot.energy(48, 3)  << endl;
//    cout << efld.ptr_sys->cpl_pot.force(48, 3)  << endl;

//    myclock::time_point start = myclock::now();
//    efld.propagate();
//    chrono::duration<double> dur = myclock::now() - start;
//    cout << dur.count() << endl;

//    efld.ptr_sys->sup_cell.coor.print();
//    efld.ptr_sys->sup_cell.velo.print();
//    efld.ptr_sys->sup_cell.coor.print();
//    cout << endl;
//    efld.ptr_sys->sup_cell.velo.print();


//    auto a = lat_hop.amplList(0, 23, vec{0,0,0});
//    for (auto& e : a)
//	e.print();
//    lat_hop.ampl(0,23, sys.bril_zone.k(0)).print(); cout << endl;
//    sys.bril_zone.kMat().print();
//    sys.brav_lat.latVecList()[0].print();
//    lat_hop.ampl(0,23, sys.bril_zone.k(3)).print(); cout << endl;
//    auto amp = lat_hop.ampl(0,23, sys.bril_zone.k(3));
//    amp /= exp(-pow(4.64/ROOT3-1, 2));
//    amp.print();
//    lat_hop.diff( 0, 23, sys.bril_zone.k(3), 0).print();
//    lat_hop.diff( 0, 25, sys.bril_zone.k(3), 0).print();
//    lat_hop.diff( 0, 11, sys.bril_zone.k(3), 0).print();

//    MolHop mol_hop( &hb, &sys.sup_cell.coor, &sys.sup_cell.mol);
//    mol_hop.ampl(32, 33).print();
//    mol_hop.ampl(32, 34).print();
//    mol_hop.diff(32, 33, 1).print();
//    mol_hop.diff(32, 34, 0).print();

//    CplHop cpl_hop(&hb, &sys.sup_cell.coor, sys.brav_lat.latVecList(), cCplHopCutoff, 1, &sys.sup_cell);
//    cpl_hop.ampl(3, 32, sys.bril_zone.k(3)).print(); cout << endl;
//    cpl_hop.ampl(3, 32, sys.bril_zone.k(0)).print(); cout << endl;
//    cpl_hop.diff(3, 32, sys.bril_zone.k(3), 0).print(); cout << endl;
//    cpl_hop.diff(3, 32, sys.bril_zone.k(0), 0).print(); cout << endl;


//    cHopCoef.print();

//    vec center_coor = vec{15, 5, 10};
//    mat mol_init_coor = joinRows(mol_coor);                          
//    mol_init_coor.each_col() += center_coor;                         
//    mat init_coor = join_rows(sys.ptrLatSupCell->eqCoor(), mol_init_coor);
//    mat init_velo = zeros(3, sys.nSite());
//    init_velo(2,sys.nSite()-1) = -0.01;
//    model.setInit(init_coor, init_velo);
    
    
//    model.propagate();

//    if (cNumTimeStep<= 1000) {
//	cRunTime.save("time.txt", raw_ascii);
//    	conv_to<vec>::from(efld.data.totE_t).save("energy.txt", raw_ascii);
//    } else {
//	cRunTime.save("time_test.txt", raw_ascii);
//    	conv_to<vec>::from(efld.data.totE_t).save("energy_test.txt", raw_ascii);
//    }
//
//    conv_to<vec>::from(efld.data.molCharge_t).save("charge.txt", raw_ascii);
//    extract(efld.data.coor_t, 0).save("xcoor.txt", raw_ascii);
//    extract(efld.data.coor_t, 1).save("ycoor.txt", raw_ascii);
//    extract(efld.data.coor_t, 2).save("zcoor.txt", raw_ascii);



//    conv_to<vec>::from(model.hist.molCOMCoor_t[0]).save("molCOMxcoor.txt", raw_ascii);
//    conv_to<vec>::from(model.hist.molCOMCoor_t[1]).save("molCOMycoor.txt", raw_ascii);
//    conv_to<vec>::from(model.hist.molCOMCoor_t[2]).save("molCOMzcoor.txt", raw_ascii);
//    conv_to<vec>::from(model.hist.molCOMVelo_t[2]).save("molCOMzvelo.txt", raw_ascii);
*/
    main_dur = myclock::now() - main_start;
    if (timing_switch) {
	cout << "main_dur = " << main_dur.count() << endl;
    }

    return 0;
}
