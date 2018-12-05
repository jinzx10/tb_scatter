#include "EFLD.h"
#include "auxmath.h"
#include "mathconst.h"
#include "join.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

template <typename T_sys>
EFLD<T_sys>::EFLD(  T_sys*	    const&	ptr_sys_,
		    double	    const&	therm_beta_,
		    double	    const&	chem_pot_,
		    bool	    const&	fric_switch_,
		    bool	    const&	lat_freeze_,
		    bool	    const&	lat_freeze_z_,
		    double	    const&	dt_,
		    arma::uword	    const&	max_num_steps_,
		    double	    const&	terminal_z_,
		    bool	    const&	print_switch_	):
    run_time(), lat_coor_t(), mol_coor_t(), lat_velo_t(), mol_velo_t(),
    mol_charge_t(), totE_t(),
    minimal_z(0), is_scattered(false), curr_time(0.0), therm_dyn_pot(0), mol_charge(0),
    lat_rand_force_corr(), mol_rand_force_corr(), cpl_rand_force_corr(),
    ptr_sys(ptr_sys_), therm_beta(therm_beta_), chem_pot(chem_pot_), 
    fric_switch(fric_switch_), lat_freeze(lat_freeze_), lat_freeze_z(lat_freeze_z_), 
    dt(dt_), dt0(dt_), max_num_steps(max_num_steps_), terminal_z(terminal_z_),
    print_switch(print_switch_), 
    lat_force(), lat_force_dt(),
    mol_force(arma::zeros<arma::mat>(3, ptr_sys->numMolSites())),
    mol_force_dt(arma::zeros<arma::mat>(3, ptr_sys->numMolSites())),
    eig_val(), eig_vec(), occ(), bff(), lat_DHDR(),
    mol_DHDR(3, arma::zeros<arma::cx_cube>(
		ptr_sys->numBasisOrbs(), ptr_sys->numBasisOrbs(), ptr_sys->numLatSites() ) )
{
    if ( !lat_freeze ) {
	lat_force = arma::zeros<arma::mat>(3, ptr_sys->numLatSites());
	lat_force_dt = arma::zeros<arma::mat>(3, ptr_sys->numLatSites());
	lat_DHDR.resize( lat_freeze_z ? 2 : 3, arma::zeros<arma::cx_cube>(
		    ptr_sys->numBasisOrbs(), ptr_sys->numBasisOrbs(), ptr_sys->numLatSites() ) );
	if (fric_switch) {
	    lat_rand_force_corr.resize( lat_freeze_z ? 3 : 6,
		    arma::zeros( ptr_sys->numLatSites(), ptr_sys->numLatSites() ) );
	    cpl_rand_force_corr.resize( lat_freeze_z ? 6 : 9,
		    arma::zeros( ptr_sys->numLatSites(), ptr_sys->numMolSites() ) );
	}
    }

    if (fric_switch)
	mol_rand_force_corr.resize(6, arma::zeros( ptr_sys->numMolSites(),
						   ptr_sys->numMolSites() ) );

    update();
    mol_force = mol_force_dt;
    if ( !lat_freeze )
	lat_force = lat_force_dt;
    save();
    minimal_z = ptr_sys->molComCoor()(2);
}

template <typename T>
void EFLD<T>::propagate() {
    arma::uword step_count = 1;
    while ( step_count < max_num_steps && ptr_sys->molComCoor()(2) < terminal_z &&
	    ptr_sys->molBondLength() < 3.5 /* 1.85 A */ ) {
	velocityVerlet();
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
void EFLD<T>::velocityVerlet() {
    ptr_sys->mol_coor += dt*(ptr_sys->mol_velo + 0.5*dt*mol_force/ptr_sys->molMassMat());
    if (!lat_freeze) {
	ptr_sys->lat_coor += dt*(ptr_sys->lat_velo + 0.5*dt*lat_force/ptr_sys->latMassMat());
    }
    update();
    ptr_sys->mol_velo += 0.5*dt*(mol_force+mol_force_dt)/ptr_sys->molMassMat();
    mol_force = mol_force_dt;
    if (!lat_freeze) {
	ptr_sys->lat_velo += 0.5*dt*(lat_force+lat_force_dt)/ptr_sys->latMassMat();
	lat_force = lat_force_dt;
    }
}

template <typename T>
void EFLD<T>::update() {
    clear();
    ptr_sys->updCplNeighbor();

    for (arma::uword k_idx = 0; k_idx < ptr_sys->bril_zone.numk(); ++k_idx) {
	arma::eig_sym( eig_val, eig_vec, ptr_sys->H(k_idx) );
	//std::cout << "tot min E = " << eig_val.min() << std::endl;
	//std::cout << "tot max E = " << eig_val.max() << std::endl;
	//arma::uvec idx_near_Ef = find( abs(eig_val) < 10.0 / therm_beta );
	//std::cout << "energy spacing = " << 20.0 / therm_beta / idx_near_Ef.size() << std::endl;
	updOcc();
	addToThermDynPot();
	addToMolCharge();
	updDHDR(ptr_sys->bril_zone.k(k_idx));
	addToMeanForce();
	if (fric_switch) {
	    updBff();
	    addToRandForceCorr();
	}
    }
    
    if (!lat_freeze) {
	lat_force_dt += ptr_sys->lat_pot.force();
	lat_force_dt += ptr_sys->diabats.neutralForceLat();
	if (lat_freeze_z) {
	    lat_force_dt.row(2).zeros();
	    ptr_sys->lat_velo.row(2).zeros();
	}
    }
    mol_force_dt += ptr_sys->diabats.neutralForceMol();

    if (fric_switch)
	addFricAndFluc();
}

template <typename T>
void EFLD<T>::save() {
    run_time.push_back(curr_time);
    lat_coor_t.push_back(ptr_sys->lat_coor);
    mol_coor_t.push_back(ptr_sys->mol_coor);
    lat_velo_t.push_back(ptr_sys->lat_velo);
    mol_velo_t.push_back(ptr_sys->mol_velo);
    mol_charge_t.push_back(mol_charge);
    totE_t.push_back( therm_dyn_pot + ptr_sys->kinE() + ptr_sys->diabats.neutral() +
		      ptr_sys->lat_pot.energy() );
}

template <typename T>
void EFLD<T>::adjDt() {
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
void EFLD<T>::print() {
    std::cout << "curr z = " << std::setw(8) << ptr_sys->molComCoor()(2) << "    "
	      //<< "target z = " << terminal_z << "    "
	      << "velo z = " << ptr_sys->molComVelo()(2) << "    "
	      << "totE = " << totE_t[totE_t.size()-1] << "    "
	      << "step count = " << mol_coor_t.size() << "    "
	      << std::endl;
}

template <typename T>
void EFLD<T>::clear() {
    therm_dyn_pot = 0.0;
    mol_charge = 0.0;
    if (fric_switch) {
	for (auto it = mol_rand_force_corr.begin(); it != mol_rand_force_corr.end(); ++it)
	    it->zeros();
	if (!lat_freeze) {
	    for (auto it = lat_rand_force_corr.begin(); it != lat_rand_force_corr.end(); ++it)
	        it->zeros();
	    for (auto it = cpl_rand_force_corr.begin(); it != cpl_rand_force_corr.end(); ++it)
	        it->zeros();
	}
    }
    lat_force_dt.zeros();
    mol_force_dt.zeros();
}

template <typename T>
double EFLD<T>::entropy() {
    double S = 0.0;
    occ.for_each( [&S] (double& val) { if ( (val>EPS) && ( (1.0-val) > EPS) )
	    S -= val*std::log(val) + (1.0-val)*std::log(1.0-val);} );
    return S;
}

template <typename T>
void EFLD<T>::updOcc() {
    occ = fermi(eig_val, therm_beta, chem_pot);
}


template <typename T>
void EFLD<T>::updDHDR(const arma::vec& k) {
    for (auto it = mol_DHDR.begin(); it != mol_DHDR.end(); ++it)
	it->zeros();

    // see notes for the meaning
    arma::cx_mat tmp_D = arma::zeros<arma::cx_mat>( ptr_sys->numBasisOrbs(),
						    ptr_sys->numBasisOrbs() );
    arma::cx_mat tmp_BC;

    // mol_DHDR
    for (arma::uword xyz = 0; xyz < 3; ++xyz) {
	for (arma::uword mol_idx = 0; mol_idx < ptr_sys->numMolSites(); ++mol_idx) {
	    // mol_DHDR from cpl hop
	    tmp_BC.zeros( 1, ptr_sys->numBasisOrbs() );
	    for (auto& nb_idx : ptr_sys->cpl_H.cpl_orb_neighbor.uniqMol())
		tmp_BC += ptr_sys->cpl_H.diff(mol_idx, nb_idx, xyz, k) * 
		    eig_vec.rows(ptr_sys->idxLatBasisOrbs(nb_idx));
	    tmp_D = eig_vec.tail_rows(1).t() * tmp_BC;
	    mol_DHDR[xyz].slice(mol_idx) += tmp_D + tmp_D.t();

	    // mol_DHDR from on-site energy h
	    mol_DHDR[xyz].slice(mol_idx) += eig_vec.tail_rows(1).t() * 
		ptr_sys->diabats.dhMol(mol_idx, xyz) * eig_vec.tail_rows(1);
	}
    }

    if (lat_freeze)
	return;

    for (auto it = lat_DHDR.begin(); it != lat_DHDR.end(); ++it)
	it->zeros();

    // lat_DHDR
    for (arma::uword xyz = 0; xyz < (lat_freeze_z ? 2 : 3); ++xyz) {
        for (arma::uword lat_idx = 0; lat_idx < ptr_sys->numLatSites(); ++lat_idx) {
    	// lat_DHDR from lat hop
    	tmp_BC.zeros( ptr_sys->lat_H.numOrbs(lat_idx), ptr_sys->numBasisOrbs() );
    	for (auto& nb_idx : ptr_sys->lat_H.neighbor.uniq(lat_idx))
    	    tmp_BC += ptr_sys->lat_H.diff(lat_idx, nb_idx, xyz, k) *
    		eig_vec.rows(ptr_sys->idxLatBasisOrbs(nb_idx));
    	tmp_D = eig_vec.rows(ptr_sys->idxLatBasisOrbs(lat_idx)).t() * tmp_BC;
    	lat_DHDR[xyz].slice(lat_idx) += tmp_D + tmp_D.t();
    
    	// lat_DHDR from cpl hop
    	tmp_BC.zeros( ptr_sys->lat_H.numOrbs(lat_idx), ptr_sys->numBasisOrbs() );
    	tmp_BC += ptr_sys->cpl_H.diff(lat_idx, xyz, k) * eig_vec.tail_rows(1);
    	tmp_D = eig_vec.rows(ptr_sys->idxLatBasisOrbs(lat_idx)).t() * tmp_BC;
    	lat_DHDR[xyz].slice(lat_idx) += tmp_D + tmp_D.t();
        } // end of lat_idx for-loop
    } // end of xyz for-loop
}

template <typename T>
void EFLD<T>::addToThermDynPot() {
    therm_dyn_pot += ( accu(occ%eig_val) - entropy()/therm_beta - chem_pot*accu(occ) )
	/ ptr_sys->bril_zone.numk();
}

template <typename T>
void EFLD<T>::addToMolCharge() {
    if (!ptr_sys->numMolSites())
	return;
    arma::subview< std::complex<double> > v = eig_vec.tail_rows(1);
    mol_charge = arma::as_scalar( ( v%arma::conj(v) ) * occ ).real() / ptr_sys->bril_zone.numk();
    // per unit cell
}

template <typename T>
void EFLD<T>::addToMeanForce() {
    arma::mat mol_DHDR_diag = arma::zeros<arma::mat>( ptr_sys->numBasisOrbs(),
						      ptr_sys->numMolSites() );
    for (arma::uword xyz = 0; xyz != 3; ++xyz) {
	for (arma::uword mol_idx = 0; mol_idx != ptr_sys->numMolSites(); ++mol_idx)
	    mol_DHDR_diag.col(mol_idx) = real( mol_DHDR[xyz].slice(mol_idx).diag() );
	mol_force_dt.row(xyz) -= occ.t() * mol_DHDR_diag / ptr_sys->bril_zone.numk();
    }
    
    if (lat_freeze)
	return;

    arma::mat lat_DHDR_diag = arma::zeros<arma::mat>( ptr_sys->numBasisOrbs(),
						      ptr_sys->numLatSites() );
    for (arma::uword xyz = 0; xyz < (lat_freeze_z ? 2 : 3); ++xyz) {
        for (arma::uword lat_idx = 0; lat_idx != ptr_sys->numLatSites(); ++lat_idx)
	    lat_DHDR_diag.col(lat_idx) = real( lat_DHDR[xyz].slice(lat_idx).diag() );
        lat_force_dt.row(xyz) -= occ.t() * lat_DHDR_diag / ptr_sys->bril_zone.numk();
    }
}

template <typename T>
void EFLD<T>::updBff() {
    double sigma = ptr_sys->level_broaden_width;
//    bff =  arma::conv_to<arma::cx_mat>::from(
//	    ( (1.0-occ) * occ.t() ) / sqrt(2.0*PI) / sigma % 
//	    exp( -arma::pow( eig_val*arma::ones<arma::rowvec>(ptr_sys->numBasisOrbs()) -
//			     arma::ones<arma::vec>(ptr_sys->numBasisOrbs())*eig_val.t(), 2 )
//		/ 2.0 / sigma / sigma ) * 2.0 / 
//	    ( 1.0 - arma::erf( eig_val*arma::ones<arma::rowvec>(ptr_sys->numBasisOrbs()) -
//			       arma::ones<arma::vec>(ptr_sys->numBasisOrbs())*eig_val.t() ) ) );

    bff =  arma::conv_to<arma::cx_mat>::from(
	    ( (1.0-occ) * occ.t() ) / sqrt(2.0*PI) / sigma % 
	    exp( -arma::pow( eig_val*arma::ones<arma::rowvec>(ptr_sys->numBasisOrbs()) -
			     arma::ones<arma::vec>(ptr_sys->numBasisOrbs())*eig_val.t(), 2 )
		/ 2.0 / sigma / sigma ) );
    //bff.for_each([](std::complex<double>& val){if (std::abs(val) < EPS) val = 0;});
}

template <typename T>
void EFLD<T>::addToRandForceCorr() {
    arma::uword dim = ptr_sys->numBasisOrbs() * ptr_sys->numBasisOrbs();
    double coef = PI / ptr_sys->bril_zone.numk();

    // mol_rand_force_corr
//    arma::uword d = ptr_sys->numBasisOrbs();
//    double sigma = ptr_sys->level_broaden_width;
//    arma::uword nsite = ptr_sys->numMolSites();
//    for (arma::uword w = 0; w != 3; ++w) // row xyz
//	for (arma::uword v = w; v != 3; ++v)// col xyz
//	    for (arma::uword r = 0; r < nsite; ++r)
//		for (arma::uword c = 0; c < nsite; ++c)
//		    for (arma::uword i = 0; i < d; ++i)
//			for (arma::uword j = 0; j < d; ++j)
//			    mol_rand_force_corr[w+(v+1)*v/2](r,c) += std::real(
//				    coef * mol_DHDR[w].slice(r)(i, j) *
//				    (1.0 - occ(j)) * mol_DHDR[v].slice(c)(j, i) *
//				    occ(i) / sqrt(2*PI) / sigma * 
//				    exp( -(eig_val(j)-eig_val(i))*(eig_val(j)-eig_val(i))
//					/ 2.0 / sigma / sigma ) );


//    for (arma::uword w = 0; w != 3; ++w) // row
//	for (arma::uword v = w; v != 3; ++v) // col
//	    mol_rand_force_corr[w+(v+1)*v/2] += coef * arma::real(
//		    arma::cx_mat( mol_DHDR[w].memptr(), dim, ptr_sys->numMolSites(), false ).t() *
//		    arma::cx_mat( arma::cx_cube( mol_DHDR[v].each_slice() % bff ).memptr(),
//				  dim, ptr_sys->numMolSites(), false ) );
//    mol_DHDR[2].slice(0).print();
//    std::cout << arma::real(mol_DHDR[2].slice(0)).index_max() << std::endl;
//    assert( false);

/*
    double kT = 1.0 / therm_beta;
    double n_var = 5.0;
    arma::uword n_Egrid = 50;
    arma::vec Egrid = arma::linspace(-n_var*kT, n_var*kT, n_Egrid);
    double dE = Egrid(1) - Egrid(0);
    arma::uword d = ptr_sys->numBasisOrbs();
    double sigma = ptr_sys->level_broaden_width;

    for (auto& E : Egrid) {
        double Tr = 0;
        for (arma::uword i = 0; i < d; ++i)
            for (arma::uword j = 0; j < d; ++j)
        	Tr += std::real( mol_DHDR[2].slice(0)(i,j) * 
        		mol_DHDR[2].slice(0)(j,i) * gaussian(E, eig_val(i), sigma) *
        		gaussian(E, eig_val(j), sigma) );
        mol_rand_force_corr[5](0,0) += coef * fermi(E, therm_beta, chem_pot) *
            (1.0-fermi(E, therm_beta, chem_pot)) * dE * Tr;
    }
*/
    
    arma::vec tmp = gaussian(eig_val, chem_pot, ptr_sys->level_broaden_width);
    arma::cx_mat bb = arma::conv_to<arma::cx_mat>::from(tmp * tmp.t());
    mol_rand_force_corr[5] += coef * arma::real(
		    arma::cx_mat( mol_DHDR[2].memptr(), dim, ptr_sys->numMolSites(), false ).st() * 
		    arma::conj( arma::cx_mat( arma::cx_cube( mol_DHDR[2].each_slice() % bb ).memptr(), dim, ptr_sys->numMolSites(), false ) ) );

    arma::mat DHDR_diag = arma::zeros(ptr_sys->numBasisOrbs(), ptr_sys->numMolSites());
    for (arma::uword col_idx = 0; col_idx < ptr_sys->numMolSites(); ++col_idx)
	DHDR_diag.col(col_idx) = arma::real(mol_DHDR[2].slice(col_idx).diag());
    mol_rand_force_corr[5] -= coef * DHDR_diag.st() * ( DHDR_diag.each_col() % (tmp%tmp));



//    for (arma::uword w = 0; w < 3; ++w)
//	for (arma::uword v = w; v < 3; ++v)
//	    mol_rand_force_corr[w+(v+1)*v/2] += coef * arma::real(
//		    arma::cx_mat( mol_DHDR[w].memptr(), dim, ptr_sys->numMolSites(), false ).st() * 
//		    arma::conj( arma::cx_mat( arma::cx_cube( mol_DHDR[v].each_slice() % bb ).memptr(), dim, ptr_sys->numMolSites(), false ) ) );

//    eig_val.print();
//    assert(false);

    if (lat_freeze)
	return;

    // lat_rand_force_corr
    for (arma::uword w = 0; w < (lat_freeze_z ? 2 : 3); ++w)
	for (arma::uword v = w; v < (lat_freeze_z ? 2 : 3); ++v)
	    lat_rand_force_corr[w+(v+1)*v/2] += coef * arma::real(
		    arma::cx_mat( lat_DHDR[w].memptr(), dim, ptr_sys->numLatSites(), false ).t() *
		    arma::cx_mat( arma::cx_cube( lat_DHDR[v].each_slice() % bff ).memptr(),
				  dim, ptr_sys->numLatSites(), false ) );
    
    // cpl_rand_force_corr
    for (arma::uword w = 0; w < (lat_freeze_z ? 2 : 3); ++w) // row-lat
        for (arma::uword v = 0; v < 3; ++v) // col-mol
	    cpl_rand_force_corr[w+v*(lat_freeze_z ? 2 : 3)] += coef * arma::real(
		    arma::cx_mat( lat_DHDR[w].memptr(), dim, ptr_sys->numLatSites(), false ).t() *
		    arma::cx_mat( arma::cx_cube( mol_DHDR[v].each_slice() % bff ).memptr(),
				  dim, ptr_sys->numMolSites(), false ) );
}

template <typename T>
void EFLD<T>::addFricAndFluc() {
    arma::mat mol_corr = join(std::vector<arma::mat>{
	    mol_rand_force_corr[0],	mol_rand_force_corr[1],	    mol_rand_force_corr[3],
	    mol_rand_force_corr[1].t(), mol_rand_force_corr[2],	    mol_rand_force_corr[4],
	    mol_rand_force_corr[3].t(), mol_rand_force_corr[4].t(), mol_rand_force_corr[5] },
	    3, 3, 'r');
    arma::mat eigvec;
    arma::vec eigval;

    arma::eig_sym(eigval, eigvec, mol_corr);
    eigval.for_each([](double& val) {if (val < 0) val = 0;});
    arma::vec fluc = eigvec * ( sqrt(eigval/dt) % arma::randn(3*ptr_sys->numMolSites()) );
    arma::vec fric = -(therm_beta*mol_corr) * arma::vectorise( ptr_sys->mol_velo.t() );
    mol_force_dt += reshape(fric+fluc, ptr_sys->numMolSites(), 3).t();

    if (lat_freeze)
	return;

    if (lat_freeze_z) {
	arma::mat lat_corr = join(std::vector<arma::mat>{
		lat_rand_force_corr[0],	    lat_rand_force_corr[1],
		lat_rand_force_corr[1].t(), lat_rand_force_corr[2] },
		2, 2, 'r');
	arma::mat cpl_corr = join(cpl_rand_force_corr, 2, 3, 'c');
        arma::mat full_corr = join(std::vector<arma::mat>{
		lat_corr, cpl_corr, cpl_corr.t(), mol_corr }, 2, 2, 'r');

        arma::eig_sym(eigval, eigvec, full_corr);
        eigval.for_each([](double& val) {if (val < 0) val = 0;});
        arma::vec fluc = eigvec * ( sqrt(eigval/dt) % arma::randn( 2*ptr_sys->numLatSites() + 
								   3*ptr_sys->numMolSites() ) );
        arma::vec fric = -(therm_beta*full_corr) * arma::join_cols(
		arma::vectorise( ptr_sys->lat_velo.head_rows(2).t() ),
		arma::vectorise( ptr_sys->mol_velo.t() ) );
        arma::vec tot = fluc + fric;
        lat_force_dt.head_rows(2) += reshape( tot.head_rows(2*ptr_sys->numLatSites()),
					      ptr_sys->numLatSites(), 2 ).t();
        mol_force_dt += reshape( tot.tail_rows(3*ptr_sys->numMolSites()),
				 ptr_sys->numMolSites(), 3 ).t();
    } else {
        arma::mat lat_corr = join(std::vector<arma::mat>{
		lat_rand_force_corr[0],	    lat_rand_force_corr[1],	lat_rand_force_corr[3],
		lat_rand_force_corr[1].t(), lat_rand_force_corr[2],	lat_rand_force_corr[4],
		lat_rand_force_corr[3].t(), lat_rand_force_corr[4].t(), lat_rand_force_corr[5] },
		3, 3, 'r');
	arma::mat cpl_corr = join(cpl_rand_force_corr, 3, 3, 'c');
        arma::mat full_corr = join(std::vector<arma::mat>{
		lat_corr, cpl_corr, cpl_corr.t(), mol_corr }, 2, 2, 'r');
        arma::eig_sym(eigval, eigvec, full_corr);
        eigval.for_each([](double& val) {if (val < 0) val = 0;});
        arma::vec fluc = eigvec * ( sqrt(eigval/dt) % arma::randn( 3*ptr_sys->numLatSites() +
								   3*ptr_sys->numMolSites() ) );
        arma::vec fric = -(therm_beta*full_corr) * arma::join_cols(
		arma::vectorise( ptr_sys->lat_velo.t() ),
		arma::vectorise( ptr_sys->mol_velo.t() ) );
        arma::vec tot = fluc + fric;
        lat_force_dt += reshape( tot.head_rows(3*ptr_sys->numLatSites()),
				 ptr_sys->numLatSites(), 3 ).t();
        mol_force_dt += reshape( tot.tail_rows(3*ptr_sys->numMolSites()),
				 ptr_sys->numMolSites(), 3 ).t();
    }
}
