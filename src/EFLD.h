#ifndef __ELECTRONIC_FRICTION_LANGEVIN_DYNAMICS_H__
#define __ELECTRONIC_FRICTION_LANGEVIN_DYNAMICS_H__

#include <armadillo>
#include "System.h"

template <typename T_sys>
class EFLD
{
    public:
	EFLD(	T_sys*		const&	    ptr_sys_,
		double		const&	    therm_beta_,
		double		const&	    chem_pot_,
		bool		const&	    fric_switch_,
		bool		const&	    lat_freeze_,
		bool		const&	    lat_freeze_z_,
		double		const&	    dt_,
		arma::uword	const&	    max_num_steps_,
		double		const&	    terminal_z_,
		bool		const&	    print_switch_	);

	void			    propagate();
	void			    update();

	std::vector<double>	    run_time;
	std::vector<arma::mat>	    lat_coor_t;
	std::vector<arma::mat>	    mol_coor_t;
	std::vector<arma::mat>	    lat_velo_t;
	std::vector<arma::mat>	    mol_velo_t;
	std::vector<double>	    mol_charge_t;
	std::vector<double>	    totE_t;

	double			    minimal_z;
	bool			    is_scattered;
	double			    curr_time;

	double			    therm_dyn_pot;
	double			    mol_charge;

	// friction add-ons
	std::vector<arma::mat>	    lat_rand_force_corr;
	std::vector<arma::mat>	    mol_rand_force_corr;
	std::vector<arma::mat>	    cpl_rand_force_corr;

    private:
	T_sys*			    ptr_sys;
	double			    therm_beta;
	double			    chem_pot;

	bool			    fric_switch; // if false, becomes Ehrenfest 
	bool			    lat_freeze;
	bool			    lat_freeze_z;

	double	    		    dt; // used in propagation
	double			    dt0; // as a reference
	arma::uword		    max_num_steps;
	double	    		    terminal_z;
	bool	    		    print_switch;

	arma::mat		    lat_force;
	arma::mat		    lat_force_dt;
	arma::mat		    mol_force;
	arma::mat		    mol_force_dt;

	// for a given k
	arma::vec		    eig_val;
	arma::cx_mat		    eig_vec;
	arma::vec		    occ;
	arma::cx_mat		    bff; // (broadened) delta * f * (1-f)
	std::vector<arma::cx_cube>  lat_DHDR;
	std::vector<arma::cx_cube>  mol_DHDR;

	void			    velocityVerlet();
	void			    adjDt();
	void			    save();
	void			    print();

	void			    clear(); // clear everything that computed from "addTo"
	void			    updOcc();
	double			    entropy(); // based on occ
	void			    updDHDR(const arma::vec& k);
	void			    addToThermDynPot();
	void			    addToMolCharge();
	void			    addToMeanForce();

	// friction add-ons
	void			    updBff();
	void			    addToRandForceCorr();
	void			    addFricAndFluc();
};

#include "EFLD.tpp"

#endif
