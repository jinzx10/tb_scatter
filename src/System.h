#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include "SupCell.h"
#include "BravLat.h"
#include "BrilZone.h"
#include <armadillo>
#include <vector>
#include <cassert>

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
class System : public SupCell
{
    public:
	System(	    std::vector<arma::vec>  const&	lat_coor_list_, // SupCell start
		    arma::vec		    const&	lat_mass_list_,
		    std::vector<arma::vec>  const&	lat_vec_list_,
		    arma::uvec		    const&	num_unit_cell_list_,
		    std::vector<arma::vec>  const&	mol_coor_list_,
		    arma::vec		    const&	mol_mass_list_,
		    arma::vec		    const&	mol_init_ref_coor_,
		    arma::vec		    const&	mol_init_trans_velo_, // SupCell end
		    arma::uvec		    const&	is_periodic_list_, // sup_brav
		    arma::uvec		    const&	num_k_list_, // bril_zone
		    T_lat_H		    const&	lat_H_,
		    T_cpl_H		    const&	cpl_H_,
		    T_lat_pot		    const&	lat_pot_,
		    T_diab		    const&	diabats_,
		    double		    const&	level_broaden_width_	);

	BravLat			    sup_brav;
	BrilZone		    bril_zone;
	T_lat_H			    lat_H;
	T_cpl_H			    cpl_H;
	T_lat_pot		    lat_pot;
	T_diab			    diabats;

	arma::vec	const&	    latBand(const arma::uword& k_idx) const;
	arma::cx_mat	const&	    latU(const arma::uword& k_idx) const;

	double			    latBandMinE() const;
	double			    latBandMaxE() const;
	arma::vec		    latBandDos(const arma::vec& E_grid) const;
	arma::vec		    latBandDos(const arma::uword& num_grids = 100) const;
	// use a range specified by latBandMinE() and latBandMaxE()

	arma::cx_mat		    H(const arma::uword& k_idx) const;
	arma::cx_vec		    Vk(const arma::uword& k_idx) const;
	double			    Gamma_h() const;
	arma::vec		    Gamma(const arma::vec& E_grid) const;
	arma::vec		    Gamma(const arma::uword& num_grids = 100) const;

	void			    updCplNeighbor();
	void			    computeLatBand();

	arma::uword		    numBasisOrbs() const;
	arma::uvec		    idxLatBasisOrbs(const arma::uword& lat_idx) const;

	double			    level_broaden_width;

    private:
	std::vector<arma::vec>	    lat_band_energy;
	std::vector<arma::cx_mat>   lat_U_mat; 
};

#include "System.tpp"

#endif
