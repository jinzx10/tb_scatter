#include "System.h"
#include "auxmath.h"
#include <cassert>

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::System(
	std::vector<arma::vec>  const&	    lat_coor_list_, // SupCell start
	arma::vec		const&	    lat_mass_list_,
	std::vector<arma::vec>  const&	    lat_vec_list_,
	arma::uvec		const&	    num_unit_cell_list_,
	std::vector<arma::vec>  const&	    mol_coor_list_,
	arma::vec		const&	    mol_mass_list_,
	arma::vec		const&	    mol_init_ref_coor_,
	arma::vec		const&	    mol_init_trans_velo_, // SupCell end
	arma::uvec		const&	    is_periodic_list_, // sup_brav
	arma::uvec		const&	    num_k_list_, // bril_zone
	T_lat_H			const&	    lat_H_,
	T_cpl_H		    	const&	    cpl_H_,
	T_lat_pot	    	const&	    lat_pot_,
	T_diab		    	const&	    diabats_,
	double		    	const&	    level_broaden_width_	):
    SupCell( lat_coor_list_, lat_mass_list_, lat_vec_list_, num_unit_cell_list_,
	     mol_coor_list_, mol_mass_list_, mol_init_ref_coor_, mol_init_trans_velo_ ),
    sup_brav(), bril_zone(),
    lat_H(lat_H_), cpl_H(cpl_H_), lat_pot(lat_pot_), diabats(diabats_),
    level_broaden_width(level_broaden_width_), lat_band_energy(), lat_U_mat()
{
    assert( is_periodic_list_.size() == lat_vec_list_.size() );

    std::vector<arma::vec> sup_lat_vec_list{};
    for (arma::uword d = 0; d != is_periodic_list_.size(); ++d)
	if (is_periodic_list_[d])
	    sup_lat_vec_list.push_back( num_unit_cell_list_[d]*lat_vec_list_[d] );

    sup_brav = BravLat(sup_lat_vec_list);
    bril_zone = BrilZone(sup_brav, num_k_list_);
    lat_H.reset(&lat_coor, lat_sup_cell.coor(), sup_lat_vec_list);
    cpl_H.reset(&lat_coor, &mol_coor, &mol.mass(), sup_lat_vec_list);
    lat_pot.reset(&lat_coor, lat_sup_cell.coor(), sup_lat_vec_list);
    diabats.reset(&lat_coor, &mol_coor, &mol_velo, &mol.mass(), sup_lat_vec_list);
    lat_band_energy = std::vector<arma::vec>(bril_zone.numk(), arma::vec{});
    lat_U_mat = std::vector<arma::cx_mat>(bril_zone.numk(), arma::cx_mat{});
    computeLatBand();
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
void System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::computeLatBand() {
    for (arma::uword k_idx = 0; k_idx < bril_zone.numk(); ++k_idx)
	eig_sym( lat_band_energy[k_idx], lat_U_mat[k_idx], lat_H.H(bril_zone.k(k_idx)) );
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::vec const& System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latBand(const arma::uword& k_idx) const {
    assert( k_idx < bril_zone.numk() );
    return lat_band_energy[k_idx];
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::cx_mat const& System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latU(const arma::uword& k_idx) const {
    assert( k_idx < lat_U_mat.size() );
    return lat_U_mat[k_idx];
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
double System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latBandMinE() const {
    assert( bril_zone.numk() && lat_band_energy[0].size() );
    double min = lat_band_energy[0].min();
    for (arma::uword k_idx = 1; k_idx < bril_zone.numk(); ++k_idx) {
	double min_k = lat_band_energy[k_idx].min();
	if (min_k < min) min = min_k;
    }
    return min;
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
double System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latBandMaxE() const {
    assert( bril_zone.numk() && lat_band_energy[0].size() );
    double max = lat_band_energy[0].max();
    for (arma::uword k_idx = 1; k_idx < bril_zone.numk(); ++k_idx) {
	double max_k = lat_band_energy[k_idx].max();
	if (max_k > max) max = max_k;
    }
    return max;
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::vec System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latBandDos(const arma::vec& E_grid) const {
    arma::vec dos = arma::zeros<arma::vec>( arma::size(E_grid) );
    for (arma::uword k_idx = 0; k_idx < bril_zone.numk(); ++k_idx)
	for (auto& E0 : latBand(k_idx))
	    dos += gaussian(E_grid, E0, level_broaden_width);
    return dos / bril_zone.numk();
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::vec System < T_diab, T_lat_H, T_cpl_H, T_lat_pot >::latBandDos(const arma::uword& num_grids) const {
    return latBandDos( arma::linspace<arma::vec>(latBandMinE(), latBandMaxE(), num_grids) );
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::cx_mat System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::H(const arma::uword& k_idx) const {
    arma::vec k = bril_zone.k(k_idx);
    arma::uword dim_lat = lat_H.numBasisOrbs();
    arma::cx_mat H = arma::zeros<arma::cx_mat>( dim_lat+1, dim_lat+1 );
    H( arma::span(0, dim_lat-1), arma::span(0, dim_lat-1) ) = lat_H.H(k);
    H( arma::span(0, dim_lat-1), dim_lat ) = cpl_H.H(k);
    H( dim_lat, arma::span(0, dim_lat-1) ) = H( arma::span(0, dim_lat-1), dim_lat ).t();
    H( dim_lat, dim_lat ) = diabats.h();
    return H;
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::cx_vec System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::Vk(const arma::uword& k_idx) const {
    return latU(k_idx).t() * cpl_H.H(bril_zone.k(k_idx));
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
double System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::Gamma_h() const {
    double gamma = 0.0;
    for (arma::uword k_idx = 0; k_idx < bril_zone.numk(); ++k_idx) {
	arma::cx_vec vk_tmp = Vk(k_idx);
	for (arma::uword i = 0; i < lat_H.numBasisOrbs(); ++i)
	    gamma += 2.0 * PI * pow( abs(vk_tmp(i)), 2 ) * 
		     gaussian(diabats.h(), lat_band_energy[k_idx](i), level_broaden_width);
    }
    return gamma / bril_zone.numk();
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::vec System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::Gamma(const arma::vec& E_grid) const {
    arma::vec gamma = arma::zeros<arma::vec>(arma::size(E_grid));
    for (arma::uword k_idx = 0; k_idx < bril_zone.numk(); ++k_idx) {
	arma::cx_vec vk_tmp = Vk(k_idx);
	for (arma::uword i = 0; i < lat_H.numBasisOrbs(); ++i)
	    gamma += 2.0 * PI * pow( abs(vk_tmp(i)), 2 ) * 
		     gaussian(E_grid, lat_band_energy[k_idx](i), level_broaden_width);
    }
    return gamma / bril_zone.numk();
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::vec System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::Gamma(const arma::uword& num_grids) const {
    return Gamma( arma::linspace<arma::vec>(latBandMinE(), latBandMaxE(), num_grids) );
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
void System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::updCplNeighbor() {
    cpl_H.updCplNeighbor();
    diabats.updCplNeighbor();
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::uword System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::numBasisOrbs() const {
    return lat_H.numBasisOrbs() + 1;
}

template < typename T_diab, typename T_lat_H, typename T_cpl_H, typename T_lat_pot >
arma::uvec System< T_diab, T_lat_H, T_cpl_H, T_lat_pot >::idxLatBasisOrbs(const arma::uword& lat_idx) const {
    return lat_H.idxBasisOrbs(lat_idx);
}
