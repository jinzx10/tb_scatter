#ifndef __TULLY2009_NO_AU_DIABATS_H__
#define __TULLY2009_NO_AU_DIABATS_H__

#include <armadillo>
#include "CplSiteNeighbor.h"

class NO_Au_Diabats
{
    public:
	NO_Au_Diabats();
	NO_Au_Diabats(	    double		    const&	neutral_morse_coef_, 
			    double	    	    const&	neutral_morse_eq_dist_,
		    	    double	    	    const&      neutral_morse_exp_decay_len,
		    	    double  	    	    const&      neutral_au_n_coef_,
		    	    double  	    	    const&      neutral_au_o_coef_,
		    	    double  	    	    const&      neutral_au_n_exp_decay_len_,
		    	    double  	    	    const&      neutral_au_o_exp_decay_len_,
			    double	    	    const&	ionic_morse_coef_, 
			    double	    	    const&	ionic_morse_eq_dist_,
		    	    double	    	    const&      ionic_morse_exp_decay_len,
		    	    double  	    	    const&      ionic_au_n_coef_,
		    	    double  	    	    const&      ionic_au_o_coef_,
		    	    double  	    	    const&      ionic_au_n_exp_decay_len_,
		    	    double  	    	    const&      ionic_au_o_exp_decay_len_,
			    double		    const&	ionic_au_n_eq_dist_,
			    double	    	    const&	image_coef_,
			    double	    	    const&	image_reg_length_,
			    double	    	    const&	image_ref_z_,
			    double	    	    const&	work_func_,
			    double	    	    const&	elec_affinity_,
			    arma::mat*		    const&	ptr_lat_coor_,
			    arma::mat*		    const&	ptr_mol_coor_,
			    arma::mat*		    const&	ptr_mol_velo_,
			    arma::vec const*	    const&	ptr_mol_mass_,
			    std::vector<arma::vec>  const&	sup_lat_vec_list_,
		    	    double		    const&      cpl_cutoff_len_,
		    	    arma::uword	    	    const&	cpl_max_ext_order_ = 1	    );

	void		    reset(arma::mat* const& ptr_lat_coor_, arma::mat* const& ptr_mol_coor_, arma::mat* ptr_mol_velo_, arma::vec const* const& ptr_mol_mass_, const std::vector<arma::vec>& sup_lat_vec_list_);
	bool		    isComplete() const;
	void		    updCplNeighbor();

	double		    neutralVibOmega() const;
	double		    ionicVibOmega() const;
	double		    neutralMomInertia() const;
	double		    ionicMomInertia() const;
	double		    neutralMorseEigEnergy(const double& vib_quanta) const;
	double		    ionicMorseEigEnergy(const double& vib_quanta) const;
	double		    neutralMaxBondLength(const double& vib_quanta) const;

	double		    numVibQuanta(const bool& state = false) const;
	double		    numRotQuanta(const bool& state = false) const;

	// Morse potential has a finite number of bound states
	double		    maxVibQuanta(const bool& state = false) const;

	// if called, must before setting vib and rot quanta
	void		    randSetMolOrientation();

	void		    setMolVibQuanta(const double& vib_quanta, const bool& state = false);
	void		    setMolRotQuanta(const double& rot_quanta, const bool& state = false);

	// orbital energy, the difference between ionic and neutral diabat
	double		    h() const; 

	// neutral diabat, contains Morse potential and coupling potential
	double		    neutral() const;
	
	// ionic diabat, contains Morse potential, coupling potential, image potential
	// shifted by work function and electron affinity
	double		    ionic() const;

	// force matrix
	arma::mat	    neutralForceLat() const;
	arma::mat	    neutralForceMol() const;
	arma::mat	    ionicForceLat() const;
	arma::mat	    ionicForceMol() const;

	// individual components of neutral and ionic diabat
	// the coupling potentials are sum of individual pairwise potentials
	double		    neutralMorse() const;
	double		    neutralCpl() const;
	double		    ionicMorse() const;
	double		    ionicCpl() const;
	double		    image() const;

	// derivative of energy/potential 
	// with respect to the (x/y/z)-component of the position of site
	double		    dhLat(const arma::uword& lat_idx, const arma::uword& xyz) const;
	double		    dhMol(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dNeutralLat(const arma::uword& lat_idx, const arma::uword& xyz) const;
	double		    dNeutralMol(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dIonicLat(const arma::uword& lat_idx, const arma::uword& xyz) const;
	double		    dIonicMol(const arma::uword& mol_idx, const arma::uword& xyz) const;

	double		    dNeutralMorse(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dIonicMorse(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dImage(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dNeutralCplLat(const arma::uword& lat_idx, const arma::uword& xyz) const;
	double		    dNeutralCplMol(const arma::uword& mol_idx, const arma::uword& xyz) const;
	double		    dIonicCplLat(const arma::uword& lat_idx, const arma::uword& xyz) const;
	double		    dIonicCplMol(const arma::uword& mol_idx, const arma::uword& xyz) const;

	// force experienced by an individual site
	arma::vec	    neutralForceLat(const arma::uword& lat_idx) const;
	arma::vec	    neutralForceMol(const arma::uword& mol_idx) const;
	arma::vec	    ionicForceLat(const arma::uword& lat_idx) const;
	arma::vec	    ionicForceMol(const arma::uword& mol_idx) const;

	// individual components of forces
	arma::vec	    neutralMorseForce(const arma::uword& mol_idx) const;
	arma::vec	    ionicMorseForce(const arma::uword& mol_idx) const;
	arma::vec	    imageForce(const arma::uword& mol_idx) const;
	arma::vec	    neutralCplForceLat(const arma::uword& lat_idx) const;
	arma::vec	    neutralCplForceMol(const arma::uword& mol_idx) const;
	arma::vec	    ionicCplForceLat(const arma::uword& lat_idx) const;
	arma::vec	    ionicCplForceMol(const arma::uword& mol_idx) const;

	CplSiteNeighbor     cpl_site_neighbor;

    private:
	void		    initialize();
	arma::vec	    molComCoor() const;
	arma::uword	    numLatSites() const;
	arma::uword	    numMolSites() const;
	arma::vec	    rNO() const;
	double		    molTotMass() const;
	double		    molReducedMass() const;
	double		    molVibKinE() const;
	double		    molRotKinE() const;

	// pair-wise coupling potentials as a function of interatomic distance
	// ionicAuN also depends on the angle between NO axis and surface normal
	// r_no = r_n - r_o
	double		    neutralAuN(const double& r) const;
	double		    neutralAuO(const double& r) const;
	double		    ionicAuN(const double& r, const arma::vec& r_no) const;
	double		    ionicAuO(const double& r) const;

	// derivatives of pair-wise potentials taken with respect to the interatomic distance
	double		    dNeutralAuN(const double& r) const;
	double		    dNeutralAuO(const double& r) const;
	double		    dIonicAuO(const double& r) const;

	// dIonicAuN_r is a partial derivative taken with respect to interatomic distance
	// dIonicAuN_cos is a partial derivative taken with respect to cos(theta)
	double		    dIonicAuN_r(const double& r, const arma::vec& r_no) const;
	double		    dIonicAuN_cos(const double& r, const arma::vec& r_no) const;

	double		    neutral_morse_coef; 
	double		    neutral_morse_eq_dist;
	double		    neutral_morse_exp_decay_len;
	double	            neutral_au_n_coef;
	double              neutral_au_o_coef;
	double              neutral_au_n_exp_decay_len;
	double              neutral_au_o_exp_decay_len;
	double		    ionic_morse_coef; 
	double	    	    ionic_morse_eq_dist;
	double		    ionic_morse_exp_decay_len;
	double              ionic_au_n_coef;
	double              ionic_au_o_coef;
	double              ionic_au_n_exp_decay_len;
	double              ionic_au_o_exp_decay_len;
	double		    ionic_au_n_eq_dist;
	double		    image_coef;
	double	    	    image_reg_length;
	double	    	    image_ref_z;
	double	    	    work_func;
	double	    	    elec_affinity;
	double		    cpl_cutoff_len;
	arma::mat*  	    ptr_lat_coor;
	arma::mat*  	    ptr_mol_coor;
	arma::mat*  	    ptr_mol_velo;
	arma::vec const*    ptr_mol_mass;
	bool		    is_complete;
};

#endif
