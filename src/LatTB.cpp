#include "LatTB.h"
#include "join.h"
#include "mathconst.h"
#include <cassert>
#include <numeric>

using namespace arma;
using namespace std;

LatTB::LatTB():
    neighbor(), on_site_E_list(), coef_mat(), decay_len_mat(), eq_dist_mat(),
    ptr_coor(nullptr), num_ref_sites(0), on_site_H(), is_complete(false)
{}

LatTB::LatTB(	vector<cx_mat>	    const&	on_site_E_list_,
		cx_mat		    const&	coef_mat_,
		mat		    const&	decay_len_mat_,
		mat		    const&	eq_dist_mat_,
		mat*		    const&	ptr_coor_,
		arma::mat	    const&	ref_coor_,
		vector<vec>	    const&	sup_lat_vec_list_,
		uword		    const&	max_dist_order_		):
    neighbor(ref_coor_, sup_lat_vec_list_, max_dist_order_), on_site_E_list(on_site_E_list_),
    coef_mat(coef_mat_), decay_len_mat(decay_len_mat_), eq_dist_mat(eq_dist_mat_),
    ptr_coor(ptr_coor_), num_ref_sites(ref_coor_.n_cols), on_site_H(), is_complete(false)
{
    assert( [this](){for (auto& e : on_site_E_list) if (!e.is_square()) return 0; return 1;}() );
    assert( coef_mat.is_square() && decay_len_mat.is_square() && eq_dist_mat.is_square() );
    uword num_orbs = 0;
    [this, &num_orbs](){for (auto& e : on_site_E_list) num_orbs += e.n_cols;}();
    assert( coef_mat.n_cols == num_orbs && decay_len_mat.n_cols == num_orbs &&
	    eq_dist_mat.n_cols == num_orbs );
    initialize();
}

void LatTB::initialize() {
    if ( !( neighbor.isComplete() && ptr_coor && 
	    ptr_coor->n_cols == num_ref_sites && ptr_coor->n_rows == 3 ) )
	return;
    is_complete = true;
    vector<cx_mat> each_on_site(numBaseSites(), cx_mat{});
    for (uword idx = 0; idx < each_on_site.size(); ++idx)
        each_on_site[idx] = kron( on_site_E_list[idx], eye( numUnitCells(), numUnitCells() ) );
    on_site_H = joinDiag(each_on_site);
}

void LatTB::reset(mat* const& ptr_coor_, const mat& ref_coor_, const vector<vec>& sup_lat_vec_list_) {
    ptr_coor = ptr_coor_;
    neighbor.reset(ref_coor_, sup_lat_vec_list_);
    num_ref_sites = ref_coor_.n_cols;
    initialize();
}

bool LatTB::isComplete() const {
    return is_complete;
}

uword LatTB::numSites() const {
    assert( isComplete() );
    return ptr_coor->n_cols;
}

uword LatTB::numBaseSites() const {
    assert( isComplete() );
    return on_site_E_list.size();
}

uword LatTB::numUnitCells() const {
    assert( isComplete() );
    return numSites() / numBaseSites();
}

uword LatTB::numOrbs() const {
    assert( isComplete() );
    return coef_mat.n_cols;
}

uword LatTB::numBasisOrbs() const {
    assert( isComplete() );
    return numOrbs() * numUnitCells();
}

uword LatTB::subLatIdx(const uword& lat_idx) const{
    assert( isComplete() && lat_idx < numSites() );
    return lat_idx / numUnitCells();
}

uword LatTB::numOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numSites() );
    return on_site_E_list[subLatIdx(lat_idx)].n_cols;
}

uvec LatTB::idxOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numSites() );
    if ( numOrbs(lat_idx) == 0 )
	return {};
    uword start_idx = 0;
    for (uword i = 0; i < subLatIdx(lat_idx); ++i)
	start_idx += on_site_E_list[i].n_cols;
    return regspace<uvec>(start_idx, start_idx+numOrbs(lat_idx)-1 );
}

uvec LatTB::idxBasisOrbs(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < numSites() );
    if ( numOrbs(lat_idx) == 0 )
	return {};
    uword start_idx = 0;
    for (uword i = 0; i < subLatIdx(lat_idx); ++i)
	start_idx += on_site_E_list[i].n_cols * numUnitCells();
    start_idx += lat_idx % numUnitCells();
    return regspace<uvec>( start_idx, numUnitCells(), 
			   start_idx + (numOrbs(lat_idx)-1) * numUnitCells() );
}

cx_mat LatTB::H(const vec& k) const {
    assert( isComplete() && k.size() == 3 );
    cx_mat H_lat = on_site_H;
    for (uword lat_to = 0; lat_to < numSites(); ++lat_to)
	for (auto& lat_from : neighbor.uniq(lat_to))
	    H_lat( idxBasisOrbs(lat_to), idxBasisOrbs(lat_from) ) += ampl(lat_to, lat_from, k);
    return H_lat;
}

cx_mat LatTB::ampl(const double& r, const cx_mat& coef, const mat& decay_len, const mat& eq_dist) const {
    assert(size(coef) == size(decay_len) && size(coef) == size(eq_dist));
    return coef % exp( -(r-eq_dist) / decay_len );
}

cx_mat LatTB::ampl(const uword& lat_to, const uword& lat_from, const vec& k) const {
    vector<cx_mat> vcm = amplList(lat_to, lat_from, k);
    if (vcm.size() == 0)
	return zeros<cx_mat>( numOrbs(lat_to), numOrbs(lat_from) );
    cx_mat cm = zeros<cx_mat>(size(vcm[0]));
    return std::accumulate(vcm.begin(), vcm.end(), cm);
}

vector<cx_mat> LatTB::amplList(const uword& lat_to, const uword& lat_from, const vec& k) const {
    assert( isComplete() && k.size() == 3 );
    assert( lat_to < numSites() && lat_from < numSites() );
    uvec nb_idx = find( neighbor(lat_to) == lat_from );
    vector<cx_mat> vcm(nb_idx.size(), cx_mat{});
    for (uword i = 0; i < nb_idx.size(); ++i) {
	vec a = neighbor.vecToBeNeighbor(lat_to, nb_idx(i));
	double dist = norm( ptr_coor->col(lat_to) - ( ptr_coor->col(lat_from) + a ) );
	uvec idx_to = idxOrbs(lat_to);
	uvec idx_from = idxOrbs(lat_from);
	vcm[i] = ampl( dist, coef_mat(idx_to, idx_from), decay_len_mat(idx_to, idx_from),
		       eq_dist_mat(idx_to, idx_from) ) * exp(I*dot(k,a));
    }
    return vcm;
}

cx_mat LatTB::diff(const uword& lat_to, const uword& lat_from, const uword& xyz, const vec& k) const {
    assert( xyz < 3 );
    vector<cx_mat> vcm = amplList(lat_to, lat_from, k);
    if (vcm.size() == 0)
	return zeros<cx_mat>( numOrbs(lat_to), numOrbs(lat_from) );
    uvec nb_idx = find( neighbor(lat_to) == lat_from );
    for (uword i = 0; i < nb_idx.size(); ++i) {
	vec a = neighbor.vecToBeNeighbor(lat_to, nb_idx(i));
	vec disp = ptr_coor->col(lat_to) - ( ptr_coor->col(lat_from) + a );
	vcm[i] %= conv_to<cx_mat>::from( -1.0 / decay_len_mat(idxOrbs(lat_to), idxOrbs(lat_from))
					* disp[xyz] / norm(disp) );
    }
    cx_mat cm = zeros<cx_mat>(size(vcm[0]));
    return std::accumulate(vcm.begin(), vcm.end(), cm);
}
