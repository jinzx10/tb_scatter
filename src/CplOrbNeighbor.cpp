#include <cassert>
#include <vector>
#include "CplOrbNeighbor.h"
#include "auxmath.h"
#include "join.h"

using namespace arma;
using namespace std;

CplOrbNeighbor::CplOrbNeighbor():
    ptr_lat_coor(nullptr), ptr_mol_coor(nullptr), ptr_mol_mass(nullptr), sup_lat_vec_list(),
    cutoff_len(0), max_ext_order(0), sup_brav_coor(), is_complete(false),
    lat_neighbor(), mol_neighbor()
{}

CplOrbNeighbor::CplOrbNeighbor(	    mat*	const&	    ptr_lat_coor_,
				    mat*	const&	    ptr_mol_coor_,
				    vec const*	const&	    ptr_mol_mass_,
			    	    vector<vec>	const&	    sup_lat_vec_list_,	
			    	    double	const&	    cutoff_len_,
			    	    uword	const&	    max_ext_order_	):
    ptr_lat_coor(ptr_lat_coor_), ptr_mol_coor(ptr_mol_coor_), ptr_mol_mass(ptr_mol_mass_),
    sup_lat_vec_list(sup_lat_vec_list_), cutoff_len(cutoff_len_), max_ext_order(max_ext_order_), 
    sup_brav_coor(), is_complete(false), lat_neighbor(), mol_neighbor()
{
    initialize();
}

void CplOrbNeighbor::initialize() {
    if ( !( ptr_lat_coor && ptr_mol_coor && ptr_mol_mass ) )
	return;

    // need more test with non-periodic case (sup_lat_vec_list.size() == 0)
    assert( sup_lat_vec_list.size() <= 3 );
    assert( [&](){for (auto& v : sup_lat_vec_list) if (v.size() != 3) return 0; return 1;}() );

    is_complete = true;
    sup_brav_coor = bravCoor( sup_lat_vec_list,
        		      ones<uvec>(sup_lat_vec_list.size()) * (2*max_ext_order+1), 's');
    lat_neighbor = vector< vector<uword> >(ptr_lat_coor->n_cols, vector<uword>{});
    update();
}

void CplOrbNeighbor::update() {
    assert( isComplete() );
    for (auto& each : lat_neighbor)
	each.clear();
    mol_neighbor.clear();

    mat ext_lat_coor = extCoor( *ptr_lat_coor, sup_brav_coor );
    for (uword i = 0; i < ext_lat_coor.n_cols; ++i)
	if ( norm( molComCoor() - ext_lat_coor.col(i) ) < cutoff_len )
	    mol_neighbor.push_back(i);

    mat ext_mol_coor = extCoor( molComCoor(), sup_brav_coor );
    for (uword lat_idx = 0; lat_idx < ptr_lat_coor->n_cols; ++lat_idx)
	for (uword i = 0; i < ext_mol_coor.n_cols; ++i)
	    if ( norm( ptr_lat_coor->col(lat_idx) - ext_mol_coor.col(i) ) < cutoff_len )
		lat_neighbor[lat_idx].push_back(i);
}

void CplOrbNeighbor::reset(mat* const& ptr_lat_coor_, mat* const& ptr_mol_coor_, vec const* const& ptr_mol_mass_, const vector<vec>& sup_lat_vec_list_) {
    ptr_lat_coor = ptr_lat_coor_;
    ptr_mol_coor = ptr_mol_coor_;
    ptr_mol_mass = ptr_mol_mass_;
    sup_lat_vec_list = sup_lat_vec_list_;
    is_complete = false;
    initialize();
}

uvec CplOrbNeighbor::lat(const uword& lat_idx) const {
    assert( isComplete() );
    assert( lat_idx < ptr_lat_coor->n_cols );
    return zeros<uvec>( lat_neighbor[lat_idx].size() );
}

uvec CplOrbNeighbor::mol() const {
    assert( isComplete() );
    vector<uword> nb = mol_neighbor;
    for (auto& i : nb)
	i %= ptr_lat_coor->n_cols;
    return conv_to<uvec>::from(nb);
}

uvec CplOrbNeighbor::uniqLat(const uword& lat_idx) const {
    return unique( lat(lat_idx) );
}

uvec CplOrbNeighbor::uniqMol() const {
    return unique( mol() );
}

vec CplOrbNeighbor::vecToBeNeighborLat(const uword& lat_idx, const uword& nb_idx) const {
    const vector<uword>& raw_nb = lat_neighbor[lat_idx];
    assert( nb_idx < raw_nb.size() );
    return sup_brav_coor.col( raw_nb[nb_idx] );
}

vec CplOrbNeighbor::vecToBeNeighborMol(const uword& nb_idx) const {
    assert( nb_idx < mol_neighbor.size() );
    return sup_brav_coor.col( mol_neighbor[nb_idx] / ptr_lat_coor->n_cols );
}

vec CplOrbNeighbor::molComCoor() const {
    mat mol_mass_mat = ones(3) * (*ptr_mol_mass).t();
    return sum( mol_mass_mat % (*ptr_mol_coor), 1) / accu(*ptr_mol_mass);
}

bool CplOrbNeighbor::isComplete() const {
    return is_complete;
}
