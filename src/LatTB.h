#ifndef __LATTICE_TIGHT_BINDING_H__
#define __LATTICE_TIGHT_BINDING_H__

#include "LatNeighbor.h"
#include <vector>
#include <armadillo>

class LatTB // neighbor relation is assumed constant
{
    public:
	LatTB();
	LatTB(	std::vector<arma::cx_mat>   const&	on_site_E_list_,
		arma::cx_mat		    const&	coef_mat_,
		arma::mat		    const&	decay_len_mat_,
		arma::mat		    const&	eq_dist_mat_,
		arma::mat*		    const&	ptr_coor_,
		arma::mat		    const&	ref_coor_,
		std::vector<arma::vec>	    const&	sup_lat_vec_list_,
		arma::uword		    const&	max_dist_order_ = 1	);

	bool			    isComplete() const;
	void			    reset(arma::mat* const& ptr_coor_, const arma::mat& ref_coor_, const std::vector<arma::vec>& sup_lat_vec_list_);

	arma::cx_mat		    H(const arma::vec& k) const;
	arma::cx_mat		    ampl(const arma::uword& lat_to, const arma::uword& lat_from, const arma::vec& k) const;
	arma::cx_mat		    diff(const arma::uword& lat_to, const arma::uword& lat_from, const arma::uword& xyz, const arma::vec& k) const;

	arma::uword		    numSites() const; // total number of lattice sites
	arma::uword		    numBaseSites() const; // number of sites in one unit cell
	arma::uword		    numUnitCells() const;
	arma::uword		    numOrbs() const; // number of orbitals in one unit cell
	arma::uword		    numBasisOrbs() const; // numOrbs()*numUnitCells()

	// coors are viewd as sorted by sub-lattice index first
	arma::uword		    subLatIdx(const arma::uword& lat_idx) const;
	arma::uword		    numOrbs(const arma::uword& lat_idx) const;
	arma::uvec		    idxOrbs(const arma::uword& lat_idx) const;
	// basis orbitals are sorted by sub-lattice first, and orbital type second
	arma::uvec		    idxBasisOrbs(const arma::uword& lat_idx) const;

	LatNeighbor		    neighbor;

    private:
	void			    initialize();
	arma::cx_mat		    ampl(const double& r, const arma::cx_mat& coef, const arma::mat& decay_len, const arma::mat& eq_dist) const;
	std::vector<arma::cx_mat>   amplList(const arma::uword& lat_to, const arma::uword& lat_from, const arma::vec& k) const;

	std::vector<arma::cx_mat>   on_site_E_list;
	arma::cx_mat		    coef_mat;
	arma::mat		    decay_len_mat;
	arma::mat		    eq_dist_mat;
	arma::mat*		    ptr_coor;
	arma::uword		    num_ref_sites;
	arma::cx_mat		    on_site_H;
	bool			    is_complete;
};

#endif
