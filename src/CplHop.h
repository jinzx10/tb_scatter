#ifndef __COUPLING_HOPPING_H__
#define __COUPLING_HOPPING_H__

#include <armadillo>
#include "CplOrbNeighbor.h"

class CplHop
{
    public:
	CplHop();
	CplHop(	    std::vector<arma::cx_vec>	const&	    coef_list_,
		    std::vector<arma::vec>	const&	    decay_len_list_,
		    std::vector<arma::vec>	const&	    eq_dist_list_,
		    arma::mat*	    	    	const&	    ptr_lat_coor_,
		    arma::mat*	    	    	const&	    ptr_mol_coor_,
		    arma::vec const*		const&	    ptr_mol_mass_,
		    std::vector<arma::vec>  	const&	    sup_lat_vec_list_,
		    double		    	const&	    cutoff_len_,
		    arma::uword	    	    	const&	    max_ext_order_ = 1	    );

	bool			    isComplete() const;
	void			    reset(arma::mat* const& ptr_lat_coor_, arma::mat* const& ptr_mol_coor_, arma::vec const* const& ptr_mol_mass, const std::vector<arma::vec>& sup_lat_vec_list_);
	void			    updCplNeighbor();

	arma::cx_vec		    H(const arma::vec& k) const;
	arma::cx_vec		    ampl(const arma::uword& lat_to, const arma::vec& k) const;
	arma::cx_vec		    diff(const arma::uword& lat_to, const arma::uword& xyz, const arma::vec& k) const;
	arma::cx_rowvec		    diff(const arma::uword& mol_to, const arma::uword& lat_from, const arma::uword& xyz, const arma::vec& k) const;

	arma::uword		    numLatSites() const;
	arma::uword		    numMolSites() const;
	arma::uword		    numLatBaseSites() const;
	arma::uword		    numLatUnitCells() const;
	arma::uword		    numLatOrbs() const;
	arma::uword		    numLatBasisOrbs() const;

	// coors are viewd as sorted by sub-lattice index first
	arma::uword		    subLatIdx(const arma::uword& lat_idx) const;
	arma::uword		    numLatOrbs(const arma::uword& lat_idx) const;
	arma::uvec		    idxLatOrbs(const arma::uword& lat_idx) const;
	// basis orbitals are sorted by sub-lattice first, and orbital type second
	arma::uvec		    idxLatBasisOrbs(const arma::uword& lat_idx) const;

	CplOrbNeighbor		    cpl_orb_neighbor;

	arma::vec		    molComCoor() const;
	double			    molTotMass() const;
    private:
	void			    initialize();

	arma::cx_vec		    ampl(const double& r, const arma::cx_vec& coef, const arma::vec& decay_len, const arma::vec& eq_dist) const;
	std::vector<arma::cx_vec>   amplListTo(const arma::uword& lat_to, const arma::vec& k) const;
	std::vector<arma::cx_rowvec> amplListFrom(const arma::uword& lat_from, const arma::vec& k) const;

	std::vector<arma::cx_vec>   coef_list;
	std::vector<arma::vec>	    decay_len_list;
	std::vector<arma::vec>	    eq_dist_list;
	arma::mat*		    ptr_lat_coor;
	arma::mat*		    ptr_mol_coor;
	arma::vec const*	    ptr_mol_mass;
	bool			    is_complete;
};

#endif
