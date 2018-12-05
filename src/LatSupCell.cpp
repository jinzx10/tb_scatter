#include <cassert>
#include "LatSupCell.h"
#include "join.h"
#include "auxmath.h"

using namespace arma;
using namespace std;

LatSupCell::LatSupCell():
    unit_cell(), brav_lat(), num_unit_cell_list(), coor_list(), mass_list()
{}

LatSupCell::LatSupCell(	    vector<vec>  	const&	    coor_list_,
		    	    vec			const&	    mass_list_,
			    vector<vec>  	const&	    lat_vec_list_,
			    uvec		const&	    num_unit_cell_list_	    ):
    unit_cell(coor_list_, mass_list_), brav_lat(lat_vec_list_),
    // bravCoor as major_coor in extCoor(), i.e. sites are sorted by subLatIdx first
    // within each sub-lattice, the fortran-style array index convention is used
    // bravCoor checks the sizes of lat_vec_list_ and num_unit_cell_list_ should be equal
    num_unit_cell_list(num_unit_cell_list_),
    coor_list( extCoor( bravCoor(lat_vec_list_, num_unit_cell_list_, 'u', 'f'),
			joinRows(coor_list_) ) ),
    mass_list( kron( unit_cell.mass(), ones(numUnitCells()) ) )
{}

uword LatSupCell::numSites() const {
    return unit_cell.numSites() * numUnitCells();
}

uword LatSupCell::numUnitCells() const {
    if (num_unit_cell_list.size() == 0)
	return 0;
    return prod(num_unit_cell_list);
}

uword LatSupCell::numUnitCells(const uword& lat_vec_idx) const {
    assert( lat_vec_idx < num_unit_cell_list.size() );
    return num_unit_cell_list(lat_vec_idx);
}

uword LatSupCell::subLatIdx(const uword& lat_idx) const {
    assert( lat_idx < numSites() );
    return lat_idx / numUnitCells();
}

mat const& LatSupCell::coor() const {
    return coor_list;
}

subview_col<double> const LatSupCell::coor(const uword& lat_idx) const {
    assert( lat_idx < numSites() );
    return coor_list.col(lat_idx);
}

vec const& LatSupCell::mass() const {
    return mass_list;
}

double const& LatSupCell::mass(const uword& lat_idx) const {
    assert( lat_idx < numSites() );
    return mass_list(lat_idx);
}
