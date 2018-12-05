#ifndef __COUPLING_SITE_NEIGHBOR_H__
#define __COUPLING_SITE_NEIGHBOR_H__

#include <armadillo>

class CplSiteNeighbor
{
    public:
	CplSiteNeighbor();
	CplSiteNeighbor(    arma::mat*		    const&	ptr_lat_coor_,
			    arma::mat*		    const&	ptr_mol_coor_,
			    std::vector<arma::vec>  const&	sup_lat_vec_list_,	
			    double		    const&	cutoff_len_,
			    arma::uword		    const&	max_ext_order_ = 1   );

	void			update();
	bool		    	isComplete() const;
	void		    	reset(arma::mat* const& ptr_lat_coor_, arma::mat* const& ptr_mol_coor_, const std::vector<arma::vec>& sup_lat_vec_list_);

	arma::uvec		lat(const arma::uword& site_idx) const;
	arma::uvec	    	mol(const arma::uword& site_idx) const;
	arma::uvec	    	uniqLat(const arma::uword& site_idx) const;
	arma::uvec	    	uniqMol(const arma::uword& site_idx) const;
	arma::vec	    	vecToBeNeighborLat(const arma::uword& site_idx, const arma::uword& nb_idx) const;
	arma::vec	    	vecToBeNeighborMol(const arma::uword& site_idx, const arma::uword& nb_idx) const;

    private:
	void			initialize();

	arma::mat*		ptr_lat_coor;
	arma::mat*	    	ptr_mol_coor;
	std::vector<arma::vec>	sup_lat_vec_list;
	double			cutoff_len;
	arma::uword	    	max_ext_order;
	arma::mat	    	sup_brav_coor;
	bool		    	is_complete;

	// first dimension: site index
	// second dimension: neighbor site index
	std::vector< std::vector<arma::uword> >	    lat_neighbor;
	std::vector< std::vector<arma::uword> >	    mol_neighbor;
};

#endif
