#ifndef __CLASSICAL_MASTER_EQUATION__
#define __CLASSICAL_MASTER_EQUATION__

#include <armadillo>
#include "System.h"

template <typename T_sys>
class CME
{
    public:
	CME(	T_sys*		const&	    ptr_sys_,
		double		const&	    therm_beta_,
		double		const&	    chem_pot_,
		bool		const&	    broaden_switch_,
		bool		const&	    lat_freeze_,
		bool		const&	    lat_freeze_z_,
		double		const&	    dt_,
		arma::uword	const&	    max_num_steps_,
		double		const&	    terminal_z_,
		bool		const&	    print_switch_	);

	void			    propagate();
	arma::mat		    corrForce();

	std::vector<double>	    run_time;
	std::vector<arma::mat>	    lat_coor_t;
	std::vector<arma::mat>	    mol_coor_t;
	std::vector<arma::mat>	    lat_velo_t;
	std::vector<arma::mat>	    mol_velo_t;
	std::vector<arma::uword>    state_t;

	std::vector<arma::mat>	    hop_coor; // mol_coor at success hops
	std::vector<double>	    hop_h;
	std::vector<arma::uword>    hop_state_from;
	arma::uword		    num_hops;
	double			    minimal_z;
	bool			    state;
	bool			    is_scattered;
	double			    curr_time;

    private:
	T_sys*			    ptr_sys;
	double			    therm_beta;
	double			    chem_pot;

	bool			    broaden_switch; // if true, becomes BCME
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

	// update force_dt, neighbors and lattice band structure (if lat_freeze == false)
	void			    update();

	void			    velocityVerlet();
	void			    hop();
	void			    adjDt();
	void			    save();
	void			    print();

	// BCME add-ons
	std::vector<arma::mat>	    DHDR_diag;
	arma::mat		    broadenedMeanForce();
	arma::mat		    unbroadenedMeanForce();
};

#include "CME.tpp"

#endif
