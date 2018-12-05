#include <cassert>
#include <vector>
#include "CplSiteNeighbor.h"
#include "auxmath.h"

using namespace arma;
using namespace std;

CplSiteNeighbor::CplSiteNeighbor():
    ptr_lat_coor(nullptr), ptr_mol_coor(nullptr), sup_lat_vec_list(),
    cutoff_len(0), max_ext_order(0), sup_brav_coor(), is_complete(false),
    lat_neighbor(), mol_neighbor()
{}

CplSiteNeighbor::CplSiteNeighbor(   mat*	const&	    ptr_lat_coor_,
				    mat*	const&	    ptr_mol_coor_,
			    	    vector<vec>	const&	    sup_lat_vec_list_,	
			    	    double	const&	    cutoff_len_,
			    	    uword	const&	    max_ext_order_	):
    ptr_lat_coor(ptr_lat_coor_), ptr_mol_coor(ptr_mol_coor_),
    sup_lat_vec_list(sup_lat_vec_list_), cutoff_len(cutoff_len_), max_ext_order(max_ext_order_),
    sup_brav_coor(), is_complete(false), lat_neighbor(), mol_neighbor()
{
    initialize();
}

void CplSiteNeighbor::initialize() {
    if ( !( ptr_lat_coor && ptr_mol_coor ) )
	return;

    // need more test with non-periodic case (sup_lat_vec_list.size() == 0)
    assert( sup_lat_vec_list.size() <= 3 );
    assert( [&](){for (auto& v : sup_lat_vec_list) if (v.size() != 3) return 0; return 1;}() );

    is_complete = true;
    sup_brav_coor = bravCoor( sup_lat_vec_list,
    			  ones<uvec>(sup_lat_vec_list.size()) * (2*max_ext_order+1), 's');
    lat_neighbor = vector< vector<uword> >(ptr_lat_coor->n_cols, vector<uword>{});
    mol_neighbor = vector< vector<uword> >(ptr_mol_coor->n_cols, vector<uword>{});
    update();
}

void CplSiteNeighbor::update() {
    assert( isComplete() );
    for (auto& each: lat_neighbor)
	each.clear();
    for (auto& each: mol_neighbor)
	each.clear();
    
    mat ext_lat_coor = extCoor(*ptr_lat_coor, sup_brav_coor);
    for (uword mol_idx = 0; mol_idx < ptr_mol_coor->n_cols; ++mol_idx)
	for (uword i = 0; i < ext_lat_coor.n_cols; ++i)
	    if ( norm( ptr_mol_coor->col(mol_idx) - ext_lat_coor.col(i) ) < cutoff_len )
		mol_neighbor[mol_idx].push_back(i);

    mat ext_mol_coor = extCoor(*ptr_mol_coor, sup_brav_coor);
    for (uword lat_idx = 0; lat_idx < ptr_mol_coor->n_cols; ++lat_idx)
	for (uword i = 0; i < ext_mol_coor.n_cols; ++i)
	    if ( norm( ptr_lat_coor->col(lat_idx) - ext_mol_coor.col(i) ) < cutoff_len )
		lat_neighbor[lat_idx].push_back(i);
}

void CplSiteNeighbor::reset(mat* const& ptr_lat_coor_, mat* const& ptr_mol_coor_, const vector<vec>& sup_lat_vec_list_) {
    ptr_lat_coor = ptr_lat_coor_;
    ptr_mol_coor = ptr_mol_coor_;
    sup_lat_vec_list = sup_lat_vec_list_;
    is_complete = false;
    initialize();
}

bool CplSiteNeighbor::isComplete() const {
    return is_complete;
}

uvec CplSiteNeighbor::lat(const uword& lat_idx) const {
    assert( isComplete() && lat_idx < ptr_lat_coor->n_cols );
    vector<uword> nb = lat_neighbor[lat_idx];
    for (auto& i : nb)
	i %= ptr_mol_coor->n_cols;
    return conv_to<uvec>::from(nb);
}

uvec CplSiteNeighbor::mol(const uword& mol_idx) const {
    assert( isComplete() && mol_idx < ptr_mol_coor->n_cols );
    vector<uword> nb = mol_neighbor[mol_idx];
    for (auto& i : nb)
	i %= ptr_lat_coor->n_cols;
    return conv_to<uvec>::from(nb);
}

uvec CplSiteNeighbor::uniqLat(const uword& lat_idx) const {
    return unique( lat(lat_idx) );
}

uvec CplSiteNeighbor::uniqMol(const uword& mol_idx) const {
    return unique( mol(mol_idx) );
}

vec CplSiteNeighbor::vecToBeNeighborLat(const uword& lat_idx, const uword& nb_idx) const {
    assert( isComplete() && lat_idx < ptr_lat_coor->n_cols );
    const vector<uword>& raw_nb = lat_neighbor[lat_idx];
    assert( nb_idx < raw_nb.size() );
    return sup_brav_coor.col( raw_nb[nb_idx] / ptr_mol_coor->n_cols );
}

vec CplSiteNeighbor::vecToBeNeighborMol(const uword& mol_idx, const uword& nb_idx) const {
    assert( isComplete() && mol_idx < ptr_mol_coor->n_cols );
    const vector<uword>& raw_nb = mol_neighbor[mol_idx];
    assert( nb_idx < raw_nb.size() );
    return sup_brav_coor.col( raw_nb[nb_idx] / ptr_lat_coor->n_cols );
}
