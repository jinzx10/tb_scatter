#include <cassert>
#include "LatNeighbor.h"
#include "mathconst.h"
#include "auxmath.h"

using namespace arma;
using namespace std;

LatNeighbor::LatNeighbor():
    ref_coor(), sup_lat_vec_list(), max_dist_order(0), sup_brav_coor(),
    is_complete(false), lat_neighbor()
{}

LatNeighbor::LatNeighbor(   arma::mat		    const&	ref_coor_,
			    std::vector<arma::vec>  const&	sup_lat_vec_list_,
			    arma::uword		    const&	max_dist_order_	    ):
    ref_coor(ref_coor_), sup_lat_vec_list(sup_lat_vec_list_), max_dist_order(max_dist_order_),
    sup_brav_coor(), is_complete(false), lat_neighbor()
{
    initialize();
}

void LatNeighbor::initialize() {
    if ( ref_coor.size() == 0  )
	return;

    // need more test with non-periodic case (sup_lat_vec_list.size() == 0)
    assert( ref_coor.n_rows == 3 && sup_lat_vec_list.size() <= 3 );
    assert( [&](){for (auto& v : sup_lat_vec_list) if (v.size() != 3) return 0; return 1;}() );
    
    sup_brav_coor = bravCoor( sup_lat_vec_list,
        		      ones<uvec>(sup_lat_vec_list.size()) * (2*max_dist_order+1), 's' );
    mat ext_coor = extCoor(ref_coor, sup_brav_coor);
    // major_coor is intra-cell coordiate (coor_mat)
    
    lat_neighbor = vector< vector< vector<uword> > >( ref_coor.n_cols,
            vector< vector<uword> >(max_dist_order+1, vector<uword>{}) );
    
    for (uword lat_idx = 0; lat_idx != ref_coor.n_cols; ++lat_idx) {
        vec dist = zeros<vec>(ext_coor.n_cols);
        for (uword idx = 0; idx < ext_coor.n_cols; ++idx)
            dist[idx] = norm(ext_coor.col(idx) - ref_coor.col(lat_idx));
        vec sdist = sort(dist);
        uvec idx_sort = sort_index(dist);
        for (uword nb_order = 0; nb_order <= max_dist_order; ++nb_order) {
            uvec shell = find(abs(sdist-sdist(0)) < EPS);
            lat_neighbor[lat_idx][nb_order] = conv_to< vector<uword> >::from(idx_sort(shell));
            sdist.shed_rows(shell.min(), shell.max());
            idx_sort.shed_rows(shell.min(), shell.max());
        }
    }
    is_complete = true;
}

void LatNeighbor::reset(const mat& ref_coor_, const vector<vec>& sup_lat_vec_list_) {
    ref_coor = ref_coor_;
    sup_lat_vec_list = sup_lat_vec_list_;
    is_complete = false;
    initialize();
}

bool LatNeighbor::isComplete() const {
    return is_complete;
}

uvec LatNeighbor::operator()(const uword& lat_idx) const {
    assert( isComplete() );
    assert( lat_idx < ref_coor.n_cols );
    auto& raw_nb = lat_neighbor[lat_idx];
    uvec nb{};
    for (uword nb_order = 1; nb_order <= max_dist_order; ++nb_order) {
	uvec uv = conv_to<uvec>::from(raw_nb[nb_order]);
	uv.for_each([this](uword& i) {i %= ref_coor.n_cols;});
	nb.insert_rows(nb.n_rows, uv);
    }
    return nb;
}

uvec LatNeighbor::uniq(const uword& lat_idx) const {
    return arma::unique(this->operator()(lat_idx));
}

vec LatNeighbor::vecToBeNeighbor(const uword& lat_idx, const uword& nb_idx) const {
    assert( isComplete() );
    auto& raw_nb = lat_neighbor[lat_idx];
    uword num_nb = 0;
    for (uword od = 1; od <= max_dist_order; ++od)
	num_nb += raw_nb[od].size();
    assert( nb_idx < num_nb );
    uword curr_idx = nb_idx;
    for (uword od = 1; od <= max_dist_order; ++od) {
	if ( curr_idx < raw_nb[od].size() )
	    return sup_brav_coor.col( raw_nb[od][curr_idx] / ref_coor.n_cols);
	else
	    curr_idx -= raw_nb[od].size();
    }
    assert( false );
}
