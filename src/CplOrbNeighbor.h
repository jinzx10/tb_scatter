#ifndef __COUPLING_ORBITAL_NEIGHBOR_H__
#define __COUPLING_ORBITAL_NEIGHBOR_H__

#include <armadillo>

class CplOrbNeighbor // assuming only one molecular orbital
{
    public:
	CplOrbNeighbor();
	CplOrbNeighbor(	    arma::mat*		    const&	    ptr_lat_coor_,
			    arma::mat*		    const&	    ptr_mol_coor_,
			    arma::vec const*	    const&	    ptr_mol_mass_,
			    std::vector<arma::vec>  const&	    sup_lat_vec_list_,	
			    double		    const&	    cutoff_len_,
			    arma::uword		    const&	    max_ext_order_ = 1   );
	// the orbital is assumed to be located at the center of mass
	// mol mass is passed for center of mass calculation

	void			update();
	bool			isComplete() const;
	void			reset(arma::mat* const& ptr_lat_coor_, arma::mat* const& ptr_mol_coor_, arma::vec const* const& ptr_mol_mass_, const std::vector<arma::vec>& sup_lat_vec_list_);

	arma::uvec	    	lat(const arma::uword& lat_idx) const;
	arma::uvec		mol() const;
	arma::uvec	    	uniqLat(const arma::uword& lat_idx) const;
	arma::uvec	    	uniqMol() const;
	arma::vec		vecToBeNeighborLat(const arma::uword& lat_idx, const arma::uword& nb_idx) const;
	arma::vec		vecToBeNeighborMol(const arma::uword& nb_idx) const;

    private:
	void			initialize();
	arma::vec		molComCoor() const;

	arma::mat*		ptr_lat_coor;
	arma::mat*		ptr_mol_coor;
	arma::vec const*	ptr_mol_mass;
	std::vector<arma::vec>	sup_lat_vec_list;
	double		    	cutoff_len;
	arma::uword	    	max_ext_order;
	arma::mat	    	sup_brav_coor;
	bool			is_complete;

	// first dimension: lattice site index; second dimension: neighbor molecular site index
	std::vector<std::vector<arma::uword> >	    lat_neighbor;
	std::vector<arma::uword>		    mol_neighbor;
};

#endif
