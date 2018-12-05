#include "CplHop.h"
#include "mathconst.h"
#include <chrono>
#include <cassert>
#include <numeric>

using namespace arma;
using namespace std;
using myclock = std::chrono::high_resolution_clock;

CplHop::CplHop():
    cpl_orb_neighbor(), coef_list(), decay_len_list(), eq_dist_list(),
    ptr_lat_coor(nullptr), ptr_mol_coor(nullptr), ptr_mol_mass(nullptr), is_complete(false)
{}

CplHop::CplHop(	    vector<cx_vec>  const&	coef_list_,
		    vector<vec>	    const&	decay_len_list_,
		    vector<vec>	    const&	eq_dist_list_,
		    mat*	    const&	ptr_lat_coor_,
		    mat*	    const&	ptr_mol_coor_,
		    vec const*	    const&	ptr_mol_mass_,
		    vector<vec>     const&	sup_lat_vec_list_,
		    double	    const&	cutoff_len_,
		    uword	    const&	max_ext_order_		):
    cpl_orb_neighbor( ptr_lat_coor_, ptr_mol_coor_, ptr_mol_mass_, 
		      sup_lat_vec_list_, cutoff_len_, max_ext_order_ ),
    coef_list(coef_list_), decay_len_list(decay_len_list_), eq_dist_list(eq_dist_list_),
    ptr_lat_coor(ptr_lat_coor_), ptr_mol_coor(ptr_mol_coor_), ptr_mol_mass(ptr_mol_mass_),
    is_complete(false)
{
    // should have some check for coef_list, decay_len_list and eq_dist_list
    initialize();
}

void CplHop::initialize() {
    if ( cpl_orb_neighbor.isComplete() )
	is_complete = true;
}

void CplHop::reset(mat* const& ptr_lat_coor_, mat* const& ptr_mol_coor_, vec const* const& ptr_mol_mass_, const vector<vec>& sup_lat_vec_list_) {
    ptr_lat_coor = ptr_lat_coor_;
    ptr_mol_coor = ptr_mol_coor_;
    ptr_mol_mass = ptr_mol_mass_;
    cpl_orb_neighbor.reset(ptr_lat_coor_, ptr_mol_coor_, ptr_mol_mass_, sup_lat_vec_list_);
    is_complete = false;
    initialize();
}

bool CplHop::isComplete() const {
    return is_complete;
}

uword CplHop::numLatSites() const {
    assert( isComplete() );
    return ptr_lat_coor->n_cols;
}

uword CplHop::numMolSites() const {
    assert( isComplete() );
    return ptr_mol_coor->n_cols;
}

uword CplHop::numLatBaseSites() const {
    return coef_list.size();
}

uword CplHop::numLatUnitCells() const {
    assert( isComplete() );
    return numLatSites() / numLatBaseSites();
}

uword CplHop::numLatOrbs() const {
    uword num = 0;
    for (auto& each : coef_list)
	num += each.n_rows;
    return num;
}

uword CplHop::numLatBasisOrbs() const {
    assert( isComplete() );
    return numLatOrbs() * numLatUnitCells();
}

uword CplHop::subLatIdx(const uword& lat_idx) const{
    assert( isComplete() && lat_idx < numLatSites() );
    return lat_idx / numLatUnitCells();
}

uword CplHop::numLatOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numLatSites() );
    return coef_list[subLatIdx(lat_idx)].n_rows;
}

uvec CplHop::idxLatOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numLatSites() );
    if ( numLatOrbs(lat_idx) == 0 )
	return {};
    uword start_idx = 0;
    for (uword i = 0; i < subLatIdx(lat_idx); ++i)
	start_idx += coef_list[i].n_rows;
    return regspace<uvec>(start_idx, start_idx+numLatOrbs(lat_idx)-1 );
}

uvec CplHop::idxLatBasisOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numLatSites() );
    if ( numLatOrbs(lat_idx) == 0 )
	return {};
    uword start_idx = 0;
    for (uword i = 0; i < subLatIdx(lat_idx); ++i)
	start_idx += coef_list[i].n_rows * numLatUnitCells();
    start_idx += lat_idx % numLatUnitCells();
    return regspace<uvec>( start_idx, numLatUnitCells(), 
			   start_idx + (numLatOrbs(lat_idx)-1) * numLatUnitCells() );
}

cx_vec CplHop::H(const vec& k) const {
    assert( isComplete() && k.size() == 3 );
    cx_vec H_cpl_upper = zeros<cx_vec>(numLatBasisOrbs());
    for (uword lat_idx = 0; lat_idx < numLatSites(); ++lat_idx)
	H_cpl_upper(idxLatBasisOrbs(lat_idx)) += ampl(lat_idx, k);
    return H_cpl_upper;
}

cx_vec CplHop::ampl(const double& r, const cx_vec& coef, const vec& decay_len, const vec& eq_dist) const {
    assert(size(coef) == size(decay_len) && size(coef) == size(eq_dist));
    return coef / (1.0 + exp((r-eq_dist) / decay_len));
}

cx_vec CplHop::ampl(const uword& lat_to, const vec& k) const {
    vector<cx_vec> vcv = amplListTo(lat_to, k);
    if (vcv.size() == 0)
	return zeros<cx_vec>( numLatOrbs(lat_to) );
    cx_vec cv = zeros<cx_vec>(size(vcv[0]));
    return std::accumulate(vcv.begin(), vcv.end(), cv);
}

vector<cx_vec> CplHop::amplListTo(const uword& lat_to, const vec& k) const {
    assert( isComplete() && k.size() == 3 );
    assert( lat_to < numLatSites() );
    uvec nb = cpl_orb_neighbor.lat(lat_to);
    vector<cx_vec> vcv(nb.size(), cx_vec{});
    for (uword i = 0; i < nb.size(); ++i) {
	vec a = cpl_orb_neighbor.vecToBeNeighborLat(lat_to, i);
	double dist = norm( ptr_lat_coor->col(lat_to) - ( molComCoor() + a) );
	uword sub_lat_idx = subLatIdx(lat_to);
	vcv[i] = ampl( dist, coef_list[sub_lat_idx], decay_len_list[sub_lat_idx],
		       eq_dist_list[sub_lat_idx] ) * exp(I*dot(k,a));
    }
    return vcv;
}

vector<cx_rowvec> CplHop::amplListFrom(const uword& lat_from, const vec& k) const {
    assert( isComplete() && k.size() == 3 );
    uvec idx = find( cpl_orb_neighbor.mol() == lat_from );
    vector<cx_rowvec> vcv(idx.size(), cx_mat{});
    for (uword i = 0; i < idx.size(); ++i) {
	vec a = cpl_orb_neighbor.vecToBeNeighborMol(idx[i]);
	double dist = norm( molComCoor() - (ptr_lat_coor->col(lat_from) + a) );
	uword sub_lat_idx = subLatIdx(lat_from);
	vcv[i] = ampl( dist, coef_list[sub_lat_idx], decay_len_list[sub_lat_idx], 
		       eq_dist_list[sub_lat_idx] ).t() * exp(I*dot(k,a));
    }
    return vcv;
}

cx_vec CplHop::diff(const uword& lat_to, const uword& xyz, const vec& k) const {
    assert( isComplete() && k.size() == 3 && xyz < 3 );
    vector<cx_vec> vcv = amplListTo(lat_to, k);
    if (vcv.size() == 0)
	return zeros<cx_vec>( numLatOrbs(lat_to) );
    for (uword i = 0; i < vcv.size(); ++i) {
	vec a = cpl_orb_neighbor.vecToBeNeighborLat(lat_to, i);
	vec disp = ptr_lat_coor->col(lat_to) - (molComCoor() + a);
	uword sub_lat_idx = subLatIdx(lat_to);
	vcv[i] %= conv_to<cx_vec>::from( -1.0 / decay_len_list[sub_lat_idx] / 
		( 1.0 + exp( -( norm(disp) - eq_dist_list[sub_lat_idx] ) / 
			    decay_len_list[sub_lat_idx] ) ) * disp[xyz] / norm(disp) );
    }
    cx_vec cv = zeros<cx_vec>(size(vcv[0]));
    return std::accumulate(vcv.begin(), vcv.end(), cv);
}

cx_rowvec CplHop::diff(const uword& mol_site_to, const uword& lat_from, const uword& xyz, const vec& k) const {
    assert( isComplete() && k.size() == 3 && xyz < 3 );
    assert( mol_site_to < numMolSites() );

    vector<cx_rowvec> vcv = amplListFrom(lat_from, k);
    uvec idx = find( cpl_orb_neighbor.mol() == lat_from );
    if (idx.size() == 0)
	return zeros<cx_rowvec>( numLatOrbs(lat_from) );
    for (uword i = 0; i < idx.size(); ++i) {
	vec a = cpl_orb_neighbor.vecToBeNeighborMol(idx[i]);
	vec disp = molComCoor() - (ptr_lat_coor->col(lat_from) + a);
	uword sub_lat_idx = subLatIdx(lat_from);
	vcv[i] %= conv_to<cx_rowvec>::from( -1.0 / decay_len_list[sub_lat_idx].t() / 
		( 1.0 + exp( -( norm(disp) - eq_dist_list[sub_lat_idx].t() ) / 
			    decay_len_list[sub_lat_idx].t() ) ) * disp[xyz] / norm(disp) *
		(*ptr_mol_mass)(mol_site_to) / molTotMass() );
    }
    cx_rowvec cv = zeros<cx_rowvec>(size(vcv[0]));
    return std::accumulate(vcv.begin(), vcv.end(), cv);
}


vec CplHop::molComCoor() const {
    mat mol_mass_mat = ones(3) * (*ptr_mol_mass).t();
    return sum( mol_mass_mat % (*ptr_mol_coor), 1) / molTotMass();
}

double CplHop::molTotMass() const {
    return accu(*ptr_mol_mass);
}

void CplHop::updCplNeighbor() {
    cpl_orb_neighbor.update();
}
