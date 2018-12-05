#include "CME.h"
#include "auxmath.h"
#include <cmath>
#include <iostream>
#include <iomanip>

template <typename T_sys>
CME<T_sys>::CME(    T_sys*	    const&	ptr_sys_,
		    double	    const&	therm_beta_,
		    double	    const&	chem_pot_,
		    bool	    const&	broaden_switch_,
		    bool	    const&	lat_freeze_,
		    bool	    const&	lat_freeze_z_,
		    double	    const&	dt_,
		    arma::uword	    const&	max_num_steps_,
		    double	    const&	terminal_z_,
		    bool	    const&	print_switch_	    ):
    run_time(), lat_coor_t(), mol_coor_t(), lat_velo_t(), mol_velo_t(), state_t(),
    hop_coor(), hop_h(), hop_state_from(), num_hops(0), minimal_z(0),
    state(false), is_scattered(false), curr_time(0.0),
    ptr_sys(ptr_sys_), therm_beta(therm_beta_), chem_pot(chem_pot_),
    broaden_switch(broaden_switch_), lat_freeze(lat_freeze_), lat_freeze_z(lat_freeze_z_),
    dt(dt_), dt0(dt_), max_num_steps(max_num_steps_), terminal_z(terminal_z_),
    print_switch(print_switch_), 
    lat_force(), lat_force_dt(),
    mol_force(arma::zeros<arma::mat>(3, ptr_sys->numMolSites())),
    mol_force_dt(arma::zeros<arma::mat>(3, ptr_sys->numMolSites())),
    DHDR_diag()
{
    if ( !lat_freeze ) {
	lat_force = arma::zeros<arma::mat>(3, ptr_sys->numLatSites());
	lat_force_dt = arma::zeros<arma::mat>(3, ptr_sys->numLatSites());
    }
    if ( broaden_switch )
	DHDR_diag.resize( 3, arma::zeros<arma::mat>(
		    ptr_sys->numBasisOrbs(), ptr_sys->numMolSites() ) );

    update();
    mol_force = mol_force_dt;
    if ( !lat_freeze )
	lat_force = lat_force_dt;
    save();
    minimal_z = ptr_sys->molComCoor()(2);
}

template <typename T>
void CME<T>::propagate() {
    arma::uword step_count = 1;
    while ( step_count < max_num_steps && ptr_sys->molComCoor()(2) < terminal_z &&
	    ptr_sys->molBondLength() < 3.5 /* about 1.85 A */ ) {
	velocityVerlet();
	hop();
	curr_time += dt;
	save();
	adjDt();
	++step_count;
	if (ptr_sys->molComCoor()(2) < minimal_z)
	    minimal_z = ptr_sys->molComCoor()(2);
	if (print_switch)
	    print();
    }
    if (ptr_sys->molComCoor()(2) >= terminal_z)
	is_scattered = true;
}

template <typename T>
void CME<T>::update() {
    ptr_sys->updCplNeighbor();
    if (!lat_freeze) { // lattice force from DHDR contribution is neglected
	ptr_sys->computeLatBand();
	lat_force_dt = ptr_sys->lat_pot.force();
	lat_force_dt += (state) ? ptr_sys->diabats.ionicForceLat() :
				  ptr_sys->diabats.neutralForceLat();
	if (lat_freeze_z) {
	    lat_force_dt.row(2).zeros();
	    ptr_sys->lat_velo.row(2).zeros();
	}
    }
    mol_force_dt = (state) ? ptr_sys->diabats.ionicForceMol() :
			     ptr_sys->diabats.neutralForceMol();
    if (broaden_switch)
	mol_force_dt += corrForce();
}

template <typename T>
void CME<T>::velocityVerlet() {
    if (!lat_freeze) {
	ptr_sys->lat_coor += dt*(ptr_sys->lat_velo + 0.5*dt*lat_force/ptr_sys->latMassMat());
    }
    ptr_sys->mol_coor += dt*(ptr_sys->mol_velo + 0.5*dt*mol_force/ptr_sys->molMassMat());
    update();
    if (!lat_freeze) {
	ptr_sys->lat_velo += 0.5*dt*(lat_force+lat_force_dt)/ptr_sys->latMassMat();
	lat_force = lat_force_dt;
    }
    ptr_sys->mol_velo += 0.5*dt*(mol_force+mol_force_dt)/ptr_sys->molMassMat();
    mol_force = mol_force_dt;
}

template <typename T>
void CME<T>::hop() {
    if ( (!state) && arma::randu() < dt * ptr_sys->Gamma_h() *
	    fermi(ptr_sys->diabats.h(), therm_beta, chem_pot) ) {
	hop_coor.push_back(ptr_sys->mol_coor);
	hop_h.push_back(ptr_sys->diabats.h());
	hop_state_from.push_back(state);
	++num_hops;
	state = true;
    }
    if ( state && arma::randu() < dt * ptr_sys->Gamma_h() * 
	    ( 1.0 - fermi(ptr_sys->diabats.h(), therm_beta, chem_pot) ) ) {
	hop_coor.push_back(ptr_sys->mol_coor);
	hop_h.push_back(ptr_sys->diabats.h());
	hop_state_from.push_back(state);
	++num_hops;
	state = false;
    }
}

template <typename T>
void CME<T>::save() {
    run_time.push_back(curr_time);
    lat_coor_t.push_back(ptr_sys->lat_coor);
    mol_coor_t.push_back(ptr_sys->mol_coor);
    lat_velo_t.push_back(ptr_sys->lat_velo);
    mol_velo_t.push_back(ptr_sys->mol_velo);
    state_t.push_back(state);
}

template <typename T>
void CME<T>::adjDt() {
    double crit_h = std::abs(ptr_sys->diabats.h());
    double crit_z = ptr_sys->molComCoor()(2);
    if ( crit_h < 0.005 && crit_z < 5.5 ) {
	dt = dt0 / 2.0;
    } else {
	if ( crit_h > 0.005 && crit_z > 5.5 )
	    dt = dt0 * 2.0;
	else
	    dt = dt0;
    }
}

template <typename T>
arma::mat CME<T>::corrForce() {
    return broadenedMeanForce() - unbroadenedMeanForce();
}

template <typename T>
arma::mat CME<T>::unbroadenedMeanForce() {
    arma::mat f_u = arma::zeros(3, ptr_sys->numMolSites());
    for (arma::uword mol_idx = 0; mol_idx < ptr_sys->numMolSites(); ++mol_idx)
	for (arma::uword xyz = 0; xyz < 3; ++xyz)
	    f_u(xyz, mol_idx) -= ptr_sys->diabats.dhMol(mol_idx, xyz);
    return f_u * ( fermi( ptr_sys->diabats.h(), therm_beta, chem_pot )  );
}

template <typename T>
arma::mat CME<T>::broadenedMeanForce() {
    arma::mat f_b = arma::zeros(3, ptr_sys->numMolSites());
    arma::cx_mat U = arma::zeros<arma::cx_mat>( ptr_sys->numBasisOrbs(), ptr_sys->numBasisOrbs() );
    arma::vec band = arma::zeros<arma::vec>( ptr_sys->numBasisOrbs() );
    arma::cx_vec tmp_D = arma::zeros<arma::cx_vec>( ptr_sys->numBasisOrbs() );
    arma::cx_mat tmp_BC = arma::zeros<arma::cx_mat>(1, ptr_sys->numBasisOrbs());

    for (arma::uword k_idx = 0; k_idx < ptr_sys->bril_zone.numk(); ++k_idx) {
	arma::vec k = ptr_sys->bril_zone.k(k_idx);
	arma::eig_sym(band, U, ptr_sys->H(k_idx));
	
	for (auto& each : DHDR_diag)
	    each.zeros();

	for (arma::uword xyz = 0; xyz < 3; ++xyz) {
	    for (arma::uword mol_idx = 0; mol_idx < ptr_sys->numMolSites(); ++mol_idx) {
		// DHDR from cpl hop
		tmp_BC.zeros();
		for (auto& nb_idx : ptr_sys->cpl_H.cpl_orb_neighbor.uniqMol())
		    tmp_BC += ptr_sys->cpl_H.diff(mol_idx, nb_idx, xyz, k) *
			U.rows(ptr_sys->lat_H.idxBasisOrbs(nb_idx));
    	    	tmp_D = ( arma::conj(U.tail_rows(1)) % tmp_BC ).st(); 
    	    	DHDR_diag[xyz].col(mol_idx) += arma::real( tmp_D ) * 2.0;

		// DHDR from on-site h
		DHDR_diag[xyz].col(mol_idx) += arma::real( arma::sum( U.tail_rows(1).t() % 
			(ptr_sys->diabats.dhMol(mol_idx, xyz) * U.tail_rows(1)).st(), 1 ) );
	    }
	}

	for (arma::uword xyz = 0; xyz < 3; ++xyz) {
	    f_b.row(xyz) -= fermi(band, therm_beta, chem_pot).t() * DHDR_diag[xyz];
	}
    }

    return f_b / ptr_sys->bril_zone.numk();
}

template <typename T>
void CME<T>::print() {
    std::cout << "curr z = " << std::setw(8) << ptr_sys->molComCoor()(2) << "    "
	      //<< "target z = " << terminal_z << "    "
	      << "velo z = " << ptr_sys->molComVelo()(2) << "    "
	      //<< "Gamma = " << ptr_sys->Gamma_h() << "    "
	      << "h = " << ptr_sys->diabats.h() << "    "
	      << "state = " << state << "    "
	      << "step count = " << mol_coor_t.size() << "    "
	      << std::endl;
}
