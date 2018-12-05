#ifndef __LATTICE_SUPER_CELL_H__
#define __LATTICE_SUPER_CELL_H__

#include <armadillo>
#include "Atoms.h"
#include "BravLat.h"

class LatSupCell
{
    public:
	LatSupCell();
	LatSupCell(	std::vector<arma::vec>	    const&	    coor_list_,
		    	arma::vec		    const&	    mass_list_,
			std::vector<arma::vec>	    const&	    lat_vec_list_,
			arma::uvec		    const&	    num_unit_cell_list_	    );

	Atoms					unit_cell;
	BravLat			    	    	brav_lat;

	arma::uword			    	numSites() const;
	arma::uword		    	    	numUnitCells() const;
	arma::uword		    	    	numUnitCells(const arma::uword& lat_vec_idx) const;
	arma::mat		    const&  	coor() const;
	arma::subview_col<double>   const	coor(const arma::uword& lat_idx) const;
	arma::vec		    const&  	mass() const;
	double			    const&	mass(const arma::uword& lat_idx) const;
	
	arma::uword				subLatIdx(const arma::uword& lat_idx) const;

    private:
	arma::uvec		num_unit_cell_list;
	arma::mat		coor_list; // by construction, sorted by sub-lattice index first
	arma::vec 		mass_list; // in consistent with the order of coor_list
};

#endif
