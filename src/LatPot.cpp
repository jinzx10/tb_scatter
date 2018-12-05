#include <cassert>
#include "LatPot.h"

using namespace arma;
using namespace std;

LatPot::LatPot():
    neighbor(), spring_const(), eq_dist(), ptr_coor(nullptr), is_complete(false)
{}

LatPot::LatPot(	    mat		    const&	spring_const_,
		    mat	    	    const&	eq_dist_,
		    mat*	    const&	ptr_coor_,
		    mat	    	    const&	ref_coor_,
		    vector<vec>	    const&	sup_lat_vec_list_,
		    uword	    const&	max_dist_order_		):
    neighbor(ref_coor_, sup_lat_vec_list_, max_dist_order_),
    spring_const(spring_const_), eq_dist(eq_dist_), ptr_coor(ptr_coor_), is_complete(false)
{
    // should have some check for spring_const and eq_dist
    initialize();
}

void LatPot::initialize() {
    if ( !( neighbor.isComplete() && ptr_coor ) )
	return;
    is_complete = true;
}

void LatPot::reset(mat* const& ptr_coor_, const mat& ref_coor_, const vector<vec>& sup_lat_vec_list_) {
    ptr_coor = ptr_coor_;
    neighbor.reset(ref_coor_, sup_lat_vec_list_);
    initialize();
}

bool LatPot::isComplete() const {
    return is_complete;
}

double LatPot::energy() const {
    assert( isComplete() );
    double E = 0.0;
    for (uword lat_idx = 0; lat_idx != numSites(); ++lat_idx)
	for (auto& nb_idx : neighbor.uniq(lat_idx))
	    E += energy(lat_idx, nb_idx);
    return 0.5*E;
}
	
mat LatPot::force() const {
    mat f = zeros(size(*ptr_coor));
    for (uword lat_idx = 0; lat_idx != numSites(); ++lat_idx)
	for (auto& nb_idx : neighbor.uniq(lat_idx)) {
	    vec f_to_lat_idx = force(lat_idx, nb_idx);
	    f.col(lat_idx) += f_to_lat_idx;
	    f.col(nb_idx) -= f_to_lat_idx;
	}
    return 0.5*f;
}

double LatPot::energy(const uword& lat_to, const uword& lat_from) const {
    assert( isComplete() );
    assert( lat_to < numSites() && lat_from < numSites() );
    double E = 0.0;
    uvec idx = find( neighbor(lat_to) == lat_from );
    for (uword i = 0; i < idx.size(); ++i) {
	vec a = neighbor.vecToBeNeighbor(lat_to, idx[i]);
	uword to_idx = subLatIdx(lat_to);
	uword from_idx = subLatIdx(lat_from);
	double var = norm( ptr_coor->col(lat_to) - (ptr_coor->col(lat_from) + a) ) -
		      eq_dist(to_idx, from_idx);
	E += 0.5 * spring_const(to_idx, from_idx) * var * var;
    }
    return E;
}

vec LatPot::force(const uword& lat_to, const uword& lat_from) const {
    assert( isComplete() );
    assert( lat_to < numSites() && lat_from < numSites() );
    vec f = zeros(3,1);
    uvec idx = find( neighbor(lat_to) == lat_from );
    for (uword i = 0; i < idx.size(); ++i) {
	vec a = neighbor.vecToBeNeighbor(lat_to, idx[i]);
	vec r = ptr_coor->col(lat_to) - (ptr_coor->col(lat_from) + a);
	uword to_idx = subLatIdx(lat_to);
	uword from_idx = subLatIdx(lat_from);
	vec disp = r / norm(r) * ( norm(r) - eq_dist(to_idx, from_idx) );
	f -= spring_const(to_idx, from_idx) * disp;
    }
    return f;
}

uword LatPot::numSites() const {
    assert( isComplete() );
    return ptr_coor->n_cols;
}

uword LatPot::numBaseSites() const {
    assert( isComplete() );
    return spring_const.n_cols;
}

uword LatPot::numUnitCells() const {
    assert( isComplete() );
    return numSites() / numBaseSites();
}

uword LatPot::subLatIdx(const uword& lat_idx) const{
    assert( isComplete() && lat_idx < numSites() );
    return lat_idx / numUnitCells();
}
