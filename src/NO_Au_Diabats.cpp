#include <cassert>
#include "NO_Au_Diabats.h"
#include "auxmath.h"
#include "mathconst.h"

using namespace std;
using namespace arma;

NO_Au_Diabats::NO_Au_Diabats():
    cpl_site_neighbor(),
    neutral_morse_coef(0), neutral_morse_eq_dist(0), neutral_morse_exp_decay_len(0),
    neutral_au_n_coef(0), neutral_au_o_coef(0), 
    neutral_au_n_exp_decay_len(0), neutral_au_o_exp_decay_len(0),
    ionic_morse_coef(0), ionic_morse_eq_dist(0), ionic_morse_exp_decay_len(0),
    ionic_au_n_coef(0), ionic_au_o_coef(0), 
    ionic_au_n_exp_decay_len(0), ionic_au_o_exp_decay_len(0),
    ionic_au_n_eq_dist(0),
    image_coef(0), image_reg_length(0), image_ref_z(0),
    work_func(0), elec_affinity(0), cpl_cutoff_len(0),
    ptr_lat_coor(nullptr), ptr_mol_coor(nullptr), ptr_mol_velo(nullptr),
    ptr_mol_mass(nullptr), is_complete(false)
{}

NO_Au_Diabats::NO_Au_Diabats(	double	    const&	neutral_morse_coef_, 
				double	    const&	neutral_morse_eq_dist_,
		    	    	double	    const&      neutral_morse_exp_decay_len_,
		    	    	double	    const&      neutral_au_n_coef_,
		    	    	double      const&      neutral_au_o_coef_,
		    	    	double      const&      neutral_au_n_exp_decay_len_,
		    	    	double      const&      neutral_au_o_exp_decay_len_,
			    	double	    const&	ionic_morse_coef_, 
			    	double	    const&	ionic_morse_eq_dist_,
		    	    	double	    const&      ionic_morse_exp_decay_len_,
		    	    	double      const&      ionic_au_n_coef_,
		    	    	double      const&      ionic_au_o_coef_,
		    	    	double      const&      ionic_au_n_exp_decay_len_,
		    	    	double      const&      ionic_au_o_exp_decay_len_,
				double	    const&	ionic_au_n_eq_dist_,
			    	double	    const&	image_coef_,
			    	double	    const&	image_reg_length_,
			    	double	    const&	image_ref_z_,
			    	double	    const&	work_func_,
			    	double	    const&	elec_affinity_,
				mat*	    const&	ptr_lat_coor_,
				mat*	    const&	ptr_mol_coor_,
				mat*	    const&	ptr_mol_velo_,
				vec const*  const&	ptr_mol_mass_,
				vector<vec> const&	sup_lat_vec_list_,
		    	    	double      const&      cpl_cutoff_len_,
		    	    	uword	    const&	cpl_max_ext_order_	    ):
    cpl_site_neighbor( ptr_lat_coor_, ptr_mol_coor_, sup_lat_vec_list_, 
		       cpl_cutoff_len_, cpl_max_ext_order_ ),
    neutral_morse_coef(neutral_morse_coef_),
    neutral_morse_eq_dist(neutral_morse_eq_dist_),
    neutral_morse_exp_decay_len(neutral_morse_exp_decay_len_),
    neutral_au_n_coef(neutral_au_n_coef_),
    neutral_au_o_coef(neutral_au_o_coef_), 
    neutral_au_n_exp_decay_len(neutral_au_n_exp_decay_len_), 
    neutral_au_o_exp_decay_len(neutral_au_o_exp_decay_len_),
    ionic_morse_coef(ionic_morse_coef_),
    ionic_morse_eq_dist(ionic_morse_eq_dist_),
    ionic_morse_exp_decay_len(ionic_morse_exp_decay_len_),
    ionic_au_n_coef(ionic_au_n_coef_),
    ionic_au_o_coef(ionic_au_o_coef_), 
    ionic_au_n_exp_decay_len(ionic_au_n_exp_decay_len_), 
    ionic_au_o_exp_decay_len(ionic_au_o_exp_decay_len_),
    ionic_au_n_eq_dist(ionic_au_n_eq_dist_),
    image_coef(image_coef_), image_reg_length(image_reg_length_), image_ref_z(image_ref_z_),
    work_func(work_func_), elec_affinity(elec_affinity_), cpl_cutoff_len(cpl_cutoff_len_),
    ptr_lat_coor(ptr_lat_coor_), ptr_mol_coor(ptr_mol_coor_), ptr_mol_velo(ptr_mol_velo_),
    ptr_mol_mass(ptr_mol_mass_), is_complete(false)
{
    initialize();
}

void NO_Au_Diabats::initialize() {
    if ( cpl_site_neighbor.isComplete() && ptr_lat_coor && ptr_mol_coor &&
	 ptr_mol_mass && ptr_mol_coor->n_cols == 2 && ptr_mol_velo->n_cols == 2)
	is_complete = true;
}

void NO_Au_Diabats::reset(mat* const& ptr_lat_coor_, mat* const& ptr_mol_coor_, mat* ptr_mol_velo_, vec const* const& ptr_mol_mass_, const vector<vec>& sup_lat_vec_list_) {
    ptr_lat_coor = ptr_lat_coor_;
    ptr_mol_coor = ptr_mol_coor_;
    ptr_mol_velo = ptr_mol_velo_;
    ptr_mol_mass = ptr_mol_mass_;
    cpl_site_neighbor.reset(ptr_lat_coor_, ptr_mol_coor_, sup_lat_vec_list_);
    is_complete = false;
    initialize();
}

bool NO_Au_Diabats::isComplete() const {
    return is_complete;
}

double NO_Au_Diabats::h() const {
    return ionic() - neutral();
}

double NO_Au_Diabats::neutral() const {
    return neutralMorse() + neutralCpl();
}

double NO_Au_Diabats::ionic() const {
    return ionicMorse() + ionicCpl() + image() + work_func - elec_affinity;
}

double NO_Au_Diabats::neutralMorse() const {
    return morse( norm(rNO()), neutral_morse_coef, neutral_morse_eq_dist,
	    neutral_morse_exp_decay_len );
}

double NO_Au_Diabats::ionicMorse() const {
    return morse( norm(rNO()), ionic_morse_coef, ionic_morse_eq_dist, ionic_morse_exp_decay_len);
}

double NO_Au_Diabats::image() const {
    return -image_coef / sqrt( pow(image_reg_length, 2) + 
			       pow(molComCoor()(2) - image_ref_z, 2) );
}

double NO_Au_Diabats::neutralCpl() const {
    double E = 0.0;
    uvec n_nb = cpl_site_neighbor.mol(0);
    uvec o_nb = cpl_site_neighbor.mol(1);
    for (uword n_nb_idx = 0; n_nb_idx < n_nb.size(); ++n_nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(0, n_nb_idx);
	double dist = norm( ptr_mol_coor->col(0) - (ptr_lat_coor->col(n_nb[n_nb_idx]) + a) );
	E += neutralAuN(dist);
    }
    for (uword o_nb_idx = 0; o_nb_idx < o_nb.size(); ++o_nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(1, o_nb_idx);
	double dist = norm( ptr_mol_coor->col(1) - (ptr_lat_coor->col(o_nb[o_nb_idx]) + a) );
	E += neutralAuO(dist);
    }
    return E;
}

double NO_Au_Diabats::ionicCpl() const {
    double E = 0.0;
    uvec n_nb = cpl_site_neighbor.mol(0);
    uvec o_nb = cpl_site_neighbor.mol(1);
    for (uword n_nb_idx = 0; n_nb_idx < n_nb.size(); ++n_nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(0, n_nb_idx);
	double dist = norm( ptr_mol_coor->col(0) - (ptr_lat_coor->col(n_nb[n_nb_idx]) + a) );
	E += ionicAuN(dist, rNO());
    }
    for (uword o_nb_idx = 0; o_nb_idx < o_nb.size(); ++o_nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(1, o_nb_idx);
	double dist = norm( ptr_mol_coor->col(1) - (ptr_lat_coor->col(o_nb[o_nb_idx]) + a) );
	E += ionicAuO(dist);
    }
    return E;
}

double NO_Au_Diabats::neutralAuN(const double& r) const {
    return neutral_au_n_coef * ( exp( - r / neutral_au_n_exp_decay_len ) - 
				 exp( - cpl_cutoff_len / neutral_au_n_exp_decay_len ) );
}

double NO_Au_Diabats::neutralAuO(const double& r) const {
    return neutral_au_o_coef * ( exp( - r / neutral_au_o_exp_decay_len ) - 
				 exp( - cpl_cutoff_len / neutral_au_o_exp_decay_len ) );
}

double NO_Au_Diabats::ionicAuN(const double& r, const vec& r_no) const {
    // differs by a sign from that of Tully2009, but ok since only cos^2 matters
    double cos_theta = r_no(2) / norm(r_no);
    double x = exp( -(r - ionic_au_n_eq_dist) / ionic_au_n_exp_decay_len );
    double xc = exp( -(cpl_cutoff_len - ionic_au_n_eq_dist) / ionic_au_n_exp_decay_len );
    return ionic_au_n_coef * (x - xc) * (x + xc - 2.0 * cos_theta * cos_theta);
}

double NO_Au_Diabats::ionicAuO(const double& r) const {
    return ionic_au_o_coef * ( exp( - r / ionic_au_o_exp_decay_len ) - 
			       exp( - cpl_cutoff_len / ionic_au_o_exp_decay_len ) );
}

double NO_Au_Diabats::dhLat(const uword& lat_idx, const uword& xyz) const {
    return dIonicLat(lat_idx, xyz) - dNeutralLat(lat_idx, xyz);
}

double NO_Au_Diabats::dhMol(const uword& mol_idx, const uword& xyz) const {
    return dIonicMol(mol_idx, xyz) - dNeutralMol(mol_idx, xyz);
}

double NO_Au_Diabats::dNeutralLat(const uword& lat_idx, const uword& xyz) const {
    return dNeutralCplLat(lat_idx, xyz);
}

double NO_Au_Diabats::dNeutralMol(const uword& mol_idx, const uword& xyz) const {
    return dNeutralMorse(mol_idx, xyz) + dNeutralCplMol(mol_idx, xyz);
}

double NO_Au_Diabats::dIonicLat(const uword& lat_idx, const uword& xyz) const {
    return dIonicCplLat(lat_idx, xyz);
}

double NO_Au_Diabats::dIonicMol(const uword& mol_idx, const uword& xyz) const {
    return dIonicMorse(mol_idx, xyz) + dIonicCplMol(mol_idx, xyz) + dImage(mol_idx, xyz);
}

double NO_Au_Diabats::dNeutralMorse(const uword& mol_idx, const uword& xyz) const {
    assert( isComplete() && mol_idx < numMolSites() && xyz < 3 );
    vec r = ptr_mol_coor->col(mol_idx) - ptr_mol_coor->col(1-mol_idx);
    double dist = norm(r);
    return r(xyz) / dist * 
	dMorse(dist, neutral_morse_coef, neutral_morse_eq_dist, neutral_morse_exp_decay_len);
}

double NO_Au_Diabats::dIonicMorse(const uword& mol_idx, const uword& xyz) const {
    assert( isComplete() && mol_idx < numMolSites() && xyz < 3 );
    vec r = ptr_mol_coor->col(mol_idx) - ptr_mol_coor->col(1-mol_idx);
    double dist = norm(r);
    return r(xyz) / dist * 
	dMorse(dist, ionic_morse_coef, ionic_morse_eq_dist, ionic_morse_exp_decay_len);
}

double NO_Au_Diabats::dImage(const uword& mol_idx, const uword& xyz) const {
    assert( isComplete() && mol_idx < numMolSites() && xyz < 3 );
    if ( xyz != 2)
	return 0.0;
    double mass = (*ptr_mol_mass)(mol_idx);
    double mol_com_z = molComCoor()(2);
    return image_coef * pow( pow(image_reg_length, 2) + pow(mol_com_z - image_ref_z, 2), -1.5 ) * 
	(mol_com_z - image_ref_z) * ( mass / molTotMass() );
}

double NO_Au_Diabats::dNeutralCplLat(const uword& lat_idx, const uword& xyz) const {
    assert( isComplete() && lat_idx < numLatSites() && xyz < 3 );
    double dE = 0.0;
    uvec nb = cpl_site_neighbor.lat(lat_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborLat(lat_idx, nb_idx);
	vec r = ptr_lat_coor->col(lat_idx) - (ptr_mol_coor->col(nb[nb_idx]) + a);
	if ( nb[nb_idx] == 0 ) // nb is N
	    dE += r(xyz) / norm(r) * dNeutralAuN(norm(r));
	else
	    dE += r(xyz) / norm(r) * dNeutralAuO(norm(r));
    }
    return dE;
}

double NO_Au_Diabats::dNeutralCplMol(const uword& mol_idx, const uword& xyz) const {
    assert( isComplete() && mol_idx < numMolSites() && xyz < 3 );
    double dE = 0.0;
    uvec nb = cpl_site_neighbor.mol(mol_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(mol_idx, nb_idx);
	vec r = ptr_mol_coor->col(mol_idx) - (ptr_lat_coor->col(nb[nb_idx]) + a);
	if ( mol_idx == 0 ) // mol_idx is N
	    dE += r(xyz) / norm(r) * dNeutralAuN(norm(r));
	else
	    dE += r(xyz) / norm(r) * dNeutralAuO(norm(r));
    }
    return dE;
}

double NO_Au_Diabats::dIonicCplLat(const uword& lat_idx, const uword& xyz) const {
    assert( isComplete() && lat_idx < numLatSites() && xyz < 3 );
    double dE = 0.0;
    uvec nb = cpl_site_neighbor.lat(lat_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborLat(lat_idx, nb_idx);
	vec r = ptr_lat_coor->col(lat_idx) - (ptr_mol_coor->col(nb[nb_idx]) + a);
	if ( nb[nb_idx] == 0 ) { // nb is N
	    dE += r(xyz) / norm(r) * dIonicAuN_r(norm(r), rNO());
	} else { // nb is O
	    dE += r(xyz) / norm(r) * dIonicAuO(norm(r));
	}
    }
    return dE;
}

double NO_Au_Diabats::dIonicCplMol(const uword& mol_idx, const uword& xyz) const {
    assert( isComplete() && mol_idx < numMolSites() && xyz < 3 );
    double dE = 0.0;
    uvec nb = cpl_site_neighbor.mol(mol_idx);
    double R = norm(rNO()); // NO bond-length
    double Z = rNO()(2); // length of NO projection in z-direction
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(mol_idx, nb_idx);
	vec r = ptr_mol_coor->col(mol_idx) - (ptr_lat_coor->col(nb[nb_idx]) + a);
	if ( mol_idx == 0 ) { // N
	    dE += r(xyz) / norm(r) * dIonicAuN_r(norm(r), rNO());
	    if (xyz == 2) {
	        dE += (R*R - Z*Z) / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    } else {
	        dE -= rNO()(xyz) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    }
	} else { // O
	    dE += r(xyz) / norm(r) * dIonicAuO(norm(r));
	    // the position of O influence the NO axis angle, thus influence ionicAuN
	    if (xyz == 2) {
	        dE -= (R*R - Z*Z) / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    } else {
	        dE += rNO()(xyz) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    }
	}
    }
    return dE;
}

double NO_Au_Diabats::dNeutralAuN(const double& r) const {
    return -neutral_au_n_coef / neutral_au_n_exp_decay_len * exp( -r / neutral_au_n_exp_decay_len );
}

double NO_Au_Diabats::dNeutralAuO(const double& r) const {
    return -neutral_au_o_coef / neutral_au_o_exp_decay_len * exp( -r / neutral_au_o_exp_decay_len );
}

double NO_Au_Diabats::dIonicAuO(const double& r) const {
    return -ionic_au_o_coef / ionic_au_o_exp_decay_len * exp( -r / ionic_au_o_exp_decay_len );
}

double NO_Au_Diabats::dIonicAuN_r(const double& r, const vec& r_no) const {
    double cos_theta = r_no(2) / norm(r_no);
    double x = exp( -(r-ionic_au_n_eq_dist) / ionic_au_n_exp_decay_len );
    return 2.0 * ionic_au_n_coef / ionic_au_n_exp_decay_len * x * (cos_theta*cos_theta - x); 
}

double NO_Au_Diabats::dIonicAuN_cos(const double& r, const vec& r_no) const {
    double cos_theta = r_no(2) / norm(r_no);
    return -4.0 * ionic_au_n_coef * cos_theta * 
	( exp( -(r-ionic_au_n_eq_dist) / ionic_au_n_exp_decay_len ) - 
	  exp( -(cpl_cutoff_len-ionic_au_n_eq_dist) / ionic_au_n_exp_decay_len ) );
}

vec NO_Au_Diabats::neutralForceLat(const uword& lat_idx) const {
    return neutralCplForceLat(lat_idx);
}

vec NO_Au_Diabats::neutralForceMol(const uword& mol_idx) const {
    return neutralMorseForce(mol_idx) + neutralCplForceMol(mol_idx);
}

vec NO_Au_Diabats::ionicForceLat(const uword& lat_idx) const {
    return ionicCplForceLat(lat_idx);
}

vec NO_Au_Diabats::ionicForceMol(const uword& mol_idx) const {
    return ionicMorseForce(mol_idx) + ionicCplForceMol(mol_idx) + imageForce(mol_idx);
}

vec NO_Au_Diabats::neutralMorseForce(const uword& mol_idx) const {
    assert( isComplete() && mol_idx < numMolSites() );
    vec r = ptr_mol_coor->col(mol_idx) - ptr_mol_coor->col(1-mol_idx);
    double dist = norm(r);
    return -r / dist *
	dMorse(dist, neutral_morse_coef, neutral_morse_eq_dist, neutral_morse_exp_decay_len);
}

vec NO_Au_Diabats::ionicMorseForce(const uword& mol_idx) const {
    assert( isComplete() && mol_idx < numMolSites() );
    vec r = ptr_mol_coor->col(mol_idx) - ptr_mol_coor->col(1-mol_idx);
    double dist = norm(r);
    return -r / dist * 
	dMorse(dist, ionic_morse_coef, ionic_morse_eq_dist, ionic_morse_exp_decay_len);
}

vec NO_Au_Diabats::imageForce(const uword& mol_idx) const {
    return vec{0, 0, -dImage(mol_idx, 2)};
}

vec NO_Au_Diabats::neutralCplForceLat(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numLatSites() );
    vec f = zeros(3, 1);
    uvec nb = cpl_site_neighbor.lat(lat_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborLat(lat_idx, nb_idx);
	vec r = ptr_lat_coor->col(lat_idx) - (ptr_mol_coor->col(nb[nb_idx]) + a);
	if ( nb[nb_idx] == 0 ) { // nb is N
	    f -= r / norm(r) * dNeutralAuN(norm(r));
	} else { // nb is O
	    f -= r / norm(r) * dNeutralAuO(norm(r));
	}
    }
    return f;
}

vec NO_Au_Diabats::neutralCplForceMol(const uword& mol_idx) const {
    assert( isComplete() && mol_idx < numMolSites() );
    vec f = zeros(3, 1);
    uvec nb = cpl_site_neighbor.mol(mol_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(mol_idx, nb_idx);
	vec r = ptr_mol_coor->col(mol_idx) - (ptr_lat_coor->col(nb[nb_idx]) + a);
	if ( mol_idx == 0 ) { // mol_idx is N
	    f -= r / norm(r) * dNeutralAuN(norm(r));
	} else {
	    f -= r / norm(r) * dNeutralAuO(norm(r));
	}
    }
    return f;
}


vec NO_Au_Diabats::ionicCplForceLat(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numLatSites() );
    vec f = zeros(3, 1);
    uvec nb = cpl_site_neighbor.lat(lat_idx);
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborLat(lat_idx, nb_idx);
	vec r = ptr_lat_coor->col(lat_idx) - (ptr_mol_coor->col(nb[nb_idx]) + a);
	if ( nb[nb_idx] == 0 ) { // nb is N
	    f -= r / norm(r) * dIonicAuN_r(norm(r), rNO());
	} else { // nb is O
	    f -= r / norm(r) * dIonicAuO(norm(r));
	}
    }
    return f;
}

vec NO_Au_Diabats::ionicCplForceMol(const uword& mol_idx) const {
    assert( isComplete() && mol_idx < numMolSites() );
    vec f = zeros(3, 1);
    uvec nb = cpl_site_neighbor.mol(mol_idx);
    double R = norm(rNO()); // NO bond-length
    double Z = rNO()(2); // length of NO projection in z-direction
    for (uword nb_idx = 0; nb_idx < nb.size(); ++nb_idx) {
	vec a = cpl_site_neighbor.vecToBeNeighborMol(mol_idx, nb_idx);
	vec r = ptr_mol_coor->col(mol_idx) - (ptr_lat_coor->col(nb[nb_idx]) + a);
	if ( mol_idx == 0 ) { // N
	    f -= r / norm(r) * dIonicAuN_r(norm(r), rNO());
	    f(0) += rNO()(0) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    f(1) += rNO()(1) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    f(2) -= (R*R - Z*Z) / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	} else {
	    f -= r / norm(r) * dIonicAuO(norm(r));
	    // the position of O influence the NO axis angle, thus influence ionicAuN
	    f(0) -= rNO()(0) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    f(1) -= rNO()(1) * Z / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	    f(2) += (R*R - Z*Z) / pow(R, 3) * dIonicAuN_cos(norm(r), rNO());
	}
    }
    return f;
}

mat NO_Au_Diabats::neutralForceLat() const {
    mat f = zeros(size(*ptr_lat_coor));
    for (uword lat_idx = 0; lat_idx < f.n_cols; ++lat_idx)
	f.col(lat_idx) = neutralForceLat(lat_idx);
    return f;
}

mat NO_Au_Diabats::neutralForceMol() const {
    mat f = zeros(size(*ptr_mol_coor));
    f.col(0) = neutralForceMol(0);
    f.col(1) = neutralForceMol(1);
    return f;
}

mat NO_Au_Diabats::ionicForceLat() const {
    mat f = zeros(size(*ptr_lat_coor));
    for (uword lat_idx = 0; lat_idx < f.n_cols; ++lat_idx)
	f.col(lat_idx) = ionicForceLat(lat_idx);
    return f;
}

mat NO_Au_Diabats::ionicForceMol() const {
    mat f = zeros(size(*ptr_mol_coor));
    f.col(0) = ionicForceMol(0);
    f.col(1) = ionicForceMol(1);
    return f;
}

vec NO_Au_Diabats::molComCoor() const {
    assert( isComplete() );
    mat mol_mass_mat = ones(3) * (*ptr_mol_mass).t();
    return sum( mol_mass_mat % (*ptr_mol_coor), 1) / molTotMass();
}

uword NO_Au_Diabats::numLatSites() const {
    assert( isComplete() );
    return ptr_lat_coor->n_cols;
}

uword NO_Au_Diabats::numMolSites() const {
    assert( isComplete() );
    return ptr_mol_coor->n_cols;
}

vec NO_Au_Diabats::rNO() const {
    assert( isComplete() );
    return ptr_mol_coor->col(0) - ptr_mol_coor->col(1);
}

double NO_Au_Diabats::molTotMass() const {
    assert( isComplete() );
    return accu( (*ptr_mol_mass) );
}

double NO_Au_Diabats::molReducedMass() const {
    assert( isComplete() );
    return (*ptr_mol_mass)(0) * (*ptr_mol_mass)(1) / molTotMass();
}

void NO_Au_Diabats::updCplNeighbor() {
    cpl_site_neighbor.update();
}

double NO_Au_Diabats::neutralVibOmega() const {
    assert( isComplete() );
    return sqrt( 2.0 * neutral_morse_coef / neutral_morse_exp_decay_len / neutral_morse_exp_decay_len / molReducedMass() );
}

double NO_Au_Diabats::ionicVibOmega() const {
    assert( isComplete() );
    return sqrt( 2.0 * ionic_morse_coef / ionic_morse_exp_decay_len / ionic_morse_exp_decay_len / molReducedMass() );
}

double NO_Au_Diabats::neutralMomInertia() const {
    assert( isComplete() );
    return molReducedMass() * neutral_morse_eq_dist * neutral_morse_eq_dist;
}

double NO_Au_Diabats::ionicMomInertia() const {
    assert( isComplete() );
    return molReducedMass() * ionic_morse_eq_dist * ionic_morse_eq_dist;
}

double NO_Au_Diabats::neutralMorseEigEnergy(const double& vib_quanta) const {
    double E_HO = neutralVibOmega() * (vib_quanta + 0.5);
    return E_HO * (1.0 - E_HO / 4.0 / neutral_morse_coef);
}

double NO_Au_Diabats::ionicMorseEigEnergy(const double& vib_quanta) const {
    double E_HO = ionicVibOmega() * (vib_quanta + 0.5);
    return E_HO * (1.0 - E_HO / 4.0 / ionic_morse_coef);
}

double NO_Au_Diabats::neutralMaxBondLength(const double& vib_quanta) const {
    return neutral_morse_eq_dist - neutral_morse_exp_decay_len * 
	log( 1.0 - sqrt(neutralMorseEigEnergy(vib_quanta) / neutral_morse_coef) );
}

double NO_Au_Diabats::numVibQuanta(const bool& state) const {
    if (state) {
	if ( ionic_morse_coef < molVibKinE() + ionicMorse() ) // dissociated
	    return -1;
	return 2.0 * ( ionic_morse_coef - sqrt( ionic_morse_coef * ( ionic_morse_coef - molVibKinE() - ionicMorse() ) ) ) / ionicVibOmega() - 0.5;
    }
    if ( neutral_morse_coef < molVibKinE() + neutralMorse() ) // dissociated
	return -1;
    return 2.0 * ( neutral_morse_coef - sqrt( neutral_morse_coef * ( neutral_morse_coef - molVibKinE() - neutralMorse() ) ) ) / neutralVibOmega() - 0.5;
}

double NO_Au_Diabats::numRotQuanta(const bool& state) const {
    if (state)
	return 0.5 * ( -1.0 + sqrt( 1.0 + 8.0*ionicMomInertia()*molRotKinE() ) );
    return 0.5 * ( -1.0 + sqrt( 1.0 + 8.0*neutralMomInertia()*molRotKinE() ) );
}

double NO_Au_Diabats::maxVibQuanta(const bool& state) const {
    if (state)
	return 2.0 * ionic_morse_coef / ionicVibOmega() - 1.0;
    return 2.0 * neutral_morse_coef / neutralVibOmega() - 1.0;
}

void NO_Au_Diabats::randSetMolOrientation() {
    double azimuth = randu() * 2 * PI;
    double polar = randu() * PI;
    vec unit = {cos(azimuth)*sin(polar), sin(azimuth)*sin(polar), cos(polar)};
    vec curr_com_coor = molComCoor();
    double R_N = norm( ptr_mol_coor->col(0) - curr_com_coor );
    double R_O = norm( ptr_mol_coor->col(1) - curr_com_coor );
    ptr_mol_coor->col(0) = curr_com_coor + unit*R_N;
    ptr_mol_coor->col(1) = curr_com_coor - unit*R_O;
}

void NO_Au_Diabats::setMolVibQuanta(const double& vib_quanta, const bool& state) {
    double phase = 2.0 * PI * randu();
    double kin_ratio = pow(sin(phase), 2);
    double pot_ratio = pow(cos(phase), 2);
    vec unit = rNO() / norm(rNO());

    double vib_kinE = (state) ? ionicMorseEigEnergy(vib_quanta) * kin_ratio:
				neutralMorseEigEnergy(vib_quanta) * kin_ratio;
    double v_r = sqrt( 2.0 * vib_kinE / molReducedMass() );
    double v_sign = (sin(phase) >= 0) ? -1.0 : 1.0;
    ptr_mol_velo->col(0) += v_sign * unit * v_r * (*ptr_mol_mass)(1) / molTotMass();
    ptr_mol_velo->col(1) -= v_sign * unit * v_r * (*ptr_mol_mass)(0) / molTotMass();

    double vib_potE = (state) ? ionicMorseEigEnergy(vib_quanta) * pot_ratio:
				neutralMorseEigEnergy(vib_quanta) * pot_ratio;
    double pot_sign = (cos(phase) >= 0) ? 1.0 : -1.0;
    double bond_length = (state) ? ionic_morse_eq_dist - ionic_morse_exp_decay_len * log(1.0 - pot_sign * sqrt(vib_potE / ionic_morse_coef)) : neutral_morse_eq_dist - neutral_morse_exp_decay_len * log(1.0 - pot_sign * sqrt(vib_potE / neutral_morse_coef));
    double curr_length = norm(rNO());
    ptr_mol_coor->col(0) += (bond_length - curr_length) * (*ptr_mol_mass)(1) / molTotMass() * unit;
    ptr_mol_coor->col(1) -= (bond_length - curr_length) * (*ptr_mol_mass)(0) / molTotMass() * unit;
}

void NO_Au_Diabats::setMolRotQuanta(const double& rot_quanta, const bool& state) {
    double rot_kinE = rot_quanta * (rot_quanta + 1.0) / 2.0;
    if (state) rot_kinE /= ionicMomInertia();
    else rot_kinE /= neutralMomInertia();
    vec unit_radial = rNO() / norm(rNO());
    vec unit_polar;
    double polar = acos( unit_radial(2) );
    if ( (polar > EPS) && (PI-polar > EPS) ) { // tan(polar) is not zero
	unit_polar = unit_radial % vec{1.0/tan(polar), 1.0/tan(polar), -tan(polar)};
    } else { // polar angle is 0 or pi
	unit_polar = {cos(polar), 0, 0};
    }
    double v_tmp = sqrt( 2.0 * rot_kinE / molTotMass() );
    ptr_mol_velo->col(0) += unit_polar * v_tmp * sqrt( (*ptr_mol_mass)(1) / (*ptr_mol_mass)(0) );
    ptr_mol_velo->col(1) -= unit_polar * v_tmp * sqrt( (*ptr_mol_mass)(0) / (*ptr_mol_mass)(1) );
}

double NO_Au_Diabats::molVibKinE() const {
    vec r01 = ptr_mol_coor->col(0) - ptr_mol_coor->col(1);
    vec v01 = ptr_mol_velo->col(0) - ptr_mol_velo->col(1);
    return 0.5 * molReducedMass() * pow( dot( v01, r01/norm(r01) ), 2 );

}

double NO_Au_Diabats::molRotKinE() const {
    vec r01 = ptr_mol_coor->col(0) - ptr_mol_coor->col(1);
    vec v01 = ptr_mol_velo->col(0) - ptr_mol_velo->col(1);
    return 0.5 * molReducedMass() * ( dot( v01, v01 ) - pow( dot( v01, r01/norm(r01) ), 2 ) );
}
