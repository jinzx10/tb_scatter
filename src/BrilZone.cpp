#include <cassert>
#include "BrilZone.h"
#include "join.h"
#include "auxmath.h"

using namespace arma;
using namespace std;

BrilZone::BrilZone():
    rcp_vec_list(), num_k_list(), k_list()
{}

BrilZone::BrilZone(	BravLat	    const&	brav_lat_,
			uvec	    const&	num_k_list_		    ):
    rcp_vec_list(brav_lat_.rcpVec()), num_k_list(num_k_list_)
{
    assert( num_k_list.size() == rcp_vec_list.n_cols );
    k_list = kGrid(rcp_vec_list, num_k_list);
}

BrilZone::BrilZone(	mat	    const&	rcp_vec_list_,
			uvec	    const&	num_k_list_		    ):
    rcp_vec_list(rcp_vec_list_), num_k_list(num_k_list_)
{
    assert( num_k_list.size() == rcp_vec_list.n_cols );
    k_list = kGrid(rcp_vec_list, num_k_list);
}

mat const& BrilZone::rcpVec() const {
    return rcp_vec_list;
}

subview_col<double> const BrilZone::rcpVec(const uword& idx) const {
    assert( idx < dim() );
    return rcp_vec_list.col(idx);
}

uword BrilZone::numk() const {
    return prod(num_k_list);
}

uword BrilZone::dim() const {
    return num_k_list.size();
}

subview_col<double> const BrilZone::k(const uword& k_idx) const {
    assert( k_idx < numk() );
    return k_list.col(k_idx);
}

mat const& BrilZone::k() const {
    return k_list;
}

mat kGrid(const mat& G, const uvec& num_k) {
    assert( G.n_cols == num_k.size() );

    mat klist = G * mesh(num_k, 'n');
    if (num_k.size() == 1)
	return klist;

    uword nk = prod(num_k);
    mat G_nb = join_rows(G, -G);
    for (uword k_idx = 0; k_idx < nk; ++k_idx) {
	auto k = klist.col(k_idx);
	moveTo1stBZ(k, G_nb);
    }
    return klist;
}

void moveTo1stBZ(subview_col<double>& k, const mat& G_nb) {
    vec k_shifted_norm = zeros(G_nb.n_cols+1);
    mat k_shifted = zeros(3, G_nb.n_cols+1);
    do {
	k_shifted.tail_cols(G_nb.n_cols) = G_nb.each_col()+k;
    	k_shifted.col(0) = k;
    	uword idx = 0;
    	k_shifted_norm.for_each([&k_shifted, &idx](double& val){val = norm(k_shifted.col(idx++));});
	k = k_shifted.col(k_shifted_norm.index_min());
    } while (k_shifted_norm.index_min());
}
