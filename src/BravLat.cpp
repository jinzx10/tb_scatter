#include <cassert>
#include "BravLat.h"
#include "mathconst.h"
#include "join.h"

using namespace arma;
using namespace std;

BravLat::BravLat():
    lat_vec_list(), rcp_vec_list()
{}

BravLat::BravLat(   vector<vec>	    const&	lat_vec_list_	    )
{
    assert( lat_vec_list_.size() <= 3);
    assert( [&lat_vec_list_](){ for (auto& v : lat_vec_list_) 
	    if (v.size() != 3) return 0; return 1;}() );
    lat_vec_list = joinRows(lat_vec_list_);
    rcp_vec_list = rcp(lat_vec_list);
}

mat const& BravLat::latVec() const {
    return lat_vec_list;
}

mat const& BravLat::rcpVec() const {
    return rcp_vec_list;
}

subview_col<double> const BravLat::latVec(const uword& idx) const {
    assert( idx < dim() );
    return lat_vec_list.col(idx);
}

subview_col<double> const BravLat::rcpVec(const uword& idx) const {
    assert( idx < dim() );
    return rcp_vec_list.col(idx);
}

uword BravLat::dim() const {
    return lat_vec_list.n_cols;
}

void BravLat::setLatVec(const uword& idx, const vec& lat_vec) {
    assert( idx < dim() && lat_vec.size() == 3 );
    lat_vec_list.col(idx) = lat_vec;
    rcp_vec_list = rcp(lat_vec_list);
}

void BravLat::addLatVec(const vec& lat_vec) {
    assert( dim() < 3 && lat_vec.size() == 3 );
    lat_vec_list.insert_cols(dim(), lat_vec);
    rcp_vec_list = rcp(lat_vec_list);
}

mat rcp(const mat& a) {
    assert( a.n_cols <= 3 );
    if (a.size() == 0)
	return {};
    a.each_col( [](const vec& v) {assert( v.size() == 3 );} );

    if (a.n_cols == 1) {
	assert( norm(a) > EPS );
	return { 2.0*PI*a / dot(a, a) };
    }

    mat G(size(a));
    const vec a0 = a.col(0);
    const vec a1 = a.col(1);
    const vec a2 = (a.n_cols == 2) ? vec{cross(a0, a1)} : a.col(2);
    double V = dot( a0, cross(a1, a2) );
    assert( abs(V) > EPS );
    G.col(0) = 2.0 * PI / V * cross(a1, a2);
    G.col(1) = 2.0 * PI / V * cross(a2, a0);
    if (G.n_cols == 3)
	G.col(2) = 2.0 * PI / V * cross(a0, a1);
    return G;
}
