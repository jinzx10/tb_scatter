#ifndef __SUPER_CELL_H__
#define __SUPER_CELL_H__

#include <armadillo>
#include "LatSupCell.h"

class SupCell
{
    public:
	SupCell();
	SupCell(    std::vector<arma::vec>  const&	lat_coor_list_,
		    arma::vec		    const&	lat_mass_list_,
		    std::vector<arma::vec>  const&	lat_vec_list_,
		    arma::uvec		    const&	num_unit_cell_list_,
		    std::vector<arma::vec>  const&	mol_coor_list_,
		    arma::vec		    const&	mol_mass_list_,
		    arma::vec		    const&	mol_init_ref_coor_,
		    arma::vec		    const&	mol_init_trans_velo_	);
	// initialized mol_coor will be each coor in mol_coor_list_ plus the mol_init_ref_coor_

	LatSupCell			    lat_sup_cell;
	Atoms		    	    	    mol;

	arma::mat		    	    lat_coor;
	arma::mat		    	    mol_coor;
	arma::mat		    	    lat_velo;
	arma::mat		    	    mol_velo;

	arma::uword	    	    	    numSites() const;
	arma::uword	    	    	    numLatSites() const;
	arma::uword	    	    	    numMolSites() const;
	arma::mat		    	    latMassMat() const;
	arma::mat		    	    molMassMat() const;

	arma::vec			    latCoor(const arma::uword& atom_idx) const;
	arma::vec			    latVelo(const arma::uword& atom_idx) const;
	arma::vec		    	    molCoor(const arma::uword& atom_idx) const;
	arma::vec		    	    molVelo(const arma::uword& atom_idx) const;

	arma::vec			    molComCoor() const;
	arma::vec   		    	    molComVelo() const;
	double			    	    molBondLength() const; // instantaneous
	double			    	    molMomInertia() const; // instantaneous

	void				    setMolCoor(const arma::mat& mol_coor);
	void			    	    setMolComCoor(const arma::vec& mol_com_coor);
	void			    	    setMolComVelo(const arma::vec& mol_com_velo);
	void			    	    setMolComCoorVelo(const double& height, const arma::vec& imag_target, const double& kinE, const double& deg_incident, const double& deg_azimuth);
	// imag_target is the point where the molecule is supposed to hit assuming the trajectory is straight, height is the z-coordinate of the molecule
	void			    	    randSetMolComCoorVelo(const double& height, const arma::vec& imag_target_center, const double& kinE, const double& deg_incident);
	// random in azimuth angle and initial x, y position (sampling within a surface unit cell centered at imag_target_center, assuming the surface is z = 0 plane)
	void			    	    randSetLatVelo(const double& therm_beta); // Boltzmann
	void				    setMolBondLength(const double& l); // for diatomic mol only

	double			    	    kinE() const;
	double			    	    latKinE() const;
	double			    	    molKinE() const;
	double			    	    molComKinE() const;
	double			    	    molVibKinE() const;
	double			    	    molRotKinE() const;
};

#endif
