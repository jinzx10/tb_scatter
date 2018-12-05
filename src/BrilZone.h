#ifndef __BRILLOUIN_ZONE_H__
#define __BRILLOUIN_ZONE_H__

#include <armadillo>
#include <vector>
#include "BravLat.h"

class BrilZone
{
    public:
	BrilZone();
	BrilZone(   BravLat	    const&	brav_lat_,
		    arma::uvec	    const&	num_k_list_	    );
	BrilZone(   arma::mat	    const&	rcp_vec_list_,
		    arma::uvec	    const&	num_k_list_	    );

	arma::mat		    const&	rcpVec() const;
	arma::subview_col<double>   const	rcpVec(const arma::uword& idx) const;
	arma::uword				dim() const;
	arma::uword		    	    	numk() const;
    	arma::subview_col<double>   const	k(const arma::uword& k_idx) const;
	arma::mat		    const&	k() const;

    private:
	arma::mat		    rcp_vec_list;
	arma::uvec		    num_k_list;
    	arma::mat		    k_list;
};

arma::mat kGrid(const arma::mat& G, const arma::uvec& nk);
void moveTo1stBZ(arma::subview_col<double>& k, const arma::mat& G_nb);

#endif
