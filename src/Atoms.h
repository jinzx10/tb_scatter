#ifndef __ATOMS_H__
#define __ATOMS_H__

#include <armadillo>
#include <vector>

class Atoms
{
    public:
	Atoms();
	Atoms(	std::vector<arma::vec>	const&	    coor_list_,
		arma::vec		const&	    mass_list_	    );

	arma::uword				numSites() const;
	arma::mat		    const&	coor() const;
	arma::subview_col<double>   const	coor(const arma::uword& site_idx) const;
	arma::vec		    const&	mass() const;
	double					mass(const arma::uword& site_idx) const;
	double			    		totMass() const;
	double			    		reducedMass() const; // diatomic only

    private:
	arma::mat				coor_list;
	arma::vec		    		mass_list;
};

#endif
