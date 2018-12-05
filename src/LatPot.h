#ifndef __LATTICE_POTENTIAL_H__
#define __LATTICE_POTENTIAL_H__

#include <armadillo>
#include "LatNeighbor.h"

class LatPot
{
    public:
	LatPot();
	LatPot(	    arma::mat		    const&	spring_const_,
		    arma::mat	    	    const&	eq_dist_,
		    arma::mat*		    const&	ptr_coor_,
		    arma::mat	    	    const&	ref_coor_,
		    std::vector<arma::vec>  const&	sup_lat_vec_list_,
		    arma::uword		    const&	max_dist_order_	= 1	);

	bool		isComplete() const;
	void		reset(arma::mat* const& ptr_coor_, const arma::mat& ref_coor_, const std::vector<arma::vec>& sup_lat_vec_list_);

	double		energy() const;
	double	    	energy(const arma::uword& lat_to, const arma::uword& lat_from) const;
	arma::mat   	force() const;
	arma::vec   	force(const arma::uword& lat_to, const arma::uword& lat_from) const;

	arma::uword	numSites() const;
	arma::uword	numBaseSites() const;
	arma::uword	numUnitCells() const;
	arma::uword	subLatIdx(const arma::uword& lat_idx) const;

	LatNeighbor	neighbor;

    private:
	void		initialize();

	arma::mat	spring_const;
	arma::mat	eq_dist;
	arma::mat*	ptr_coor;
	bool		is_complete;
};

#endif
