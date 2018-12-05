#include <cassert>
#include "Atoms.h"
#include "join.h"

using namespace arma;
using namespace std;

Atoms::Atoms():
    coor_list(), mass_list()
{}

Atoms::Atoms(	vector<vec>	const&	    coor_list_,
		vec		const&	    mass_list_		)
{
    assert( coor_list_.size() == mass_list_.size() );
    assert( [&coor_list_](){for (auto& c : coor_list_) if (c.size() != 3) return 0; return 1;}() );
    assert( [&mass_list_](){for (auto& m : mass_list_) if (m <= 0) return 0; return 1;}() );
    coor_list = joinRows(coor_list_);
    mass_list = mass_list_;
}

mat const& Atoms::coor() const {
    return coor_list;
}

subview_col<double> const Atoms::coor(const uword& site_idx) const {
    assert( site_idx < numSites() );
    return coor_list.col(site_idx);
}

vec const& Atoms::mass() const {
    return mass_list;
}

double Atoms::mass(const uword& site_idx) const {
    assert( site_idx < numSites() );
    return mass_list(site_idx);
}

double Atoms::totMass() const {
    return accu( mass_list );
}

double Atoms::reducedMass() const {
    assert( numSites() == 2 );
    return mass_list(0) * mass_list(1) / ( mass_list(0) + mass_list(1) );
}

uword Atoms::numSites() const {
    return mass_list.size();
}
