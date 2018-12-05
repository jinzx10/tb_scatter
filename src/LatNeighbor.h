#ifndef __LATTICE_NEIGHBOR_H__
#define __LATTICE_NEIGHBOR_H__

#include <vector>
#include <armadillo>

class LatNeighbor
{
    public:
	LatNeighbor();
	LatNeighbor(	arma::mat		const&	    ref_coor_,
			std::vector<arma::vec>  const&	    sup_lat_vec_list_,
			arma::uword		const&	    max_dist_order_ = 1	    );

	bool			isComplete() const;
	void			reset(const arma::mat& ref_coor_, const std::vector<arma::vec>& sup_lat_vec_list_);

	arma::uvec		operator()(const arma::uword& lat_idx) const;
	arma::uvec		uniq(const arma::uword& lat_idx) const;
	arma::vec		vecToBeNeighbor(const arma::uword& lat_idx, const arma::uword& nb_idx) const;

    private:
	void			initialize();

	arma::mat		ref_coor;
	std::vector<arma::vec>	sup_lat_vec_list;    
	arma::uword		max_dist_order;
	arma::mat		sup_brav_coor;
	bool			is_complete;

	// first dimension: site index
	// second dimension: neighbor order
	// third dimension: neighbor site index 
	std::vector< std::vector< std::vector<arma::uword> > >	lat_neighbor;
};

#endif
