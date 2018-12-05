#include <cassert>
#include "SupCell.h"
#include "join.h"
#include "mathconst.h"

using namespace arma;
using namespace std;

SupCell::SupCell():
    lat_sup_cell(), mol(), lat_coor(), mol_coor(), lat_velo(), mol_velo()
{}

SupCell::SupCell(   vector<vec>		const&	    lat_coor_list_,
		    vec			const&	    lat_mass_list_,
		    vector<vec>		const&	    lat_vec_list_,
		    uvec		const&	    num_unit_cell_list_,
		    vector<vec>		const&	    mol_coor_list_,
		    vec			const&	    mol_mass_list_,
		    vec			const&	    mol_init_ref_coor_,
		    vec			const&	    mol_init_trans_velo_	):
    lat_sup_cell( lat_coor_list_, lat_mass_list_, lat_vec_list_, num_unit_cell_list_ ),
    mol( mol_coor_list_, mol_mass_list_ ),
    lat_coor( lat_sup_cell.coor() ), mol_coor( mol.coor() ),
    lat_velo( zeros( size(lat_coor) ) ), mol_velo( zeros( size(mol_coor) ) )
{
    mol_coor.each_col() += mol_init_ref_coor_;
    mol_velo.each_col() += mol_init_trans_velo_;
}

uword SupCell::numSites() const {
    return numLatSites() + numMolSites();
}

uword SupCell::numLatSites() const {
    return lat_sup_cell.numSites();
}

uword SupCell::numMolSites() const {
    return mol.numSites();
}

mat SupCell::latMassMat() const {
    return ones(3) * lat_sup_cell.mass().t();
}

mat SupCell::molMassMat() const {
    return ones(3) * mol.mass().t();
}

vec SupCell::latCoor(const uword& lat_idx) const {
    assert( lat_idx < numLatSites() );
    return lat_coor.col( lat_idx );
}

vec SupCell::latVelo(const uword& lat_idx) const {
    assert( lat_idx < numLatSites() );
    return lat_velo.col( lat_idx );
}

vec SupCell::molCoor(const uword& mol_idx) const {
    assert( mol_idx < numMolSites() );
    return mol_coor.col( mol_idx );
}

vec SupCell::molVelo(const uword& mol_idx) const {
    assert( mol_idx < numMolSites() );
    return mol_velo.col( mol_idx );
}

vec SupCell::molComCoor() const {
    return sum( molMassMat() % mol_coor, 1 ) / mol.totMass();
}

vec SupCell::molComVelo() const {
    return sum( molMassMat() % mol_velo, 1) / mol.totMass();
}

double SupCell::molBondLength() const {
    assert( numMolSites() == 2 );
    return norm( molCoor(0) - molCoor(1) );
}

double SupCell::molMomInertia() const {
    assert( numMolSites() == 2 );
    return mol.reducedMass() * pow( molBondLength(), 2 );
}

void SupCell::setMolCoor(const mat& mol_coor_) {
    assert( arma::size(mol_coor) == arma::size(mol_coor_) );
    mol_coor = mol_coor_;
}

void SupCell::setMolComCoor(const vec& mol_com_coor) {
    mol_coor.each_col() += mol_com_coor - molComCoor();
}

void SupCell::setMolComVelo(const vec& mol_com_velo) {
    mol_velo.each_col() += mol_com_velo - molComVelo();
}

void SupCell::setMolComCoorVelo(const double& height, const vec& imag_target, const double& kinE, const double& deg_incident, const double& deg_azimuth) {
    assert( kinE > 0 && imag_target.size() == 3 );
    double incident = deg_incident * PI / 180.0;
    double azimuth = deg_azimuth * PI / 180.0;
    vec unit_vec = {cos(azimuth)*sin(incident), sin(azimuth)*sin(incident), -cos(incident)};
    double dist = ( height - imag_target(2) ) / cos(incident);
    double v = sqrt( 2.0 * kinE / mol.totMass() );
    setMolComVelo( v * unit_vec );
    setMolComCoor( imag_target - unit_vec*dist );
}

void SupCell::randSetMolComCoorVelo(const double& height, const vec& imag_target_center, const double& kinE, const double& deg_incident) {
    // make sure z = 0 is a surface
    assert( lat_sup_cell.brav_lat.dim() >= 2 );
    assert( lat_sup_cell.brav_lat.latVec(0)(2) < EPS && lat_sup_cell.brav_lat.latVec(1)(2) < EPS );

    double deg_azimuth = randu() * 360;

    // imag_target is sampled from a surface unit cell
    vec imag_target = imag_target_center + (randu()-0.5)*lat_sup_cell.brav_lat.latVec(0) + 
					   (randu()-0.5)*lat_sup_cell.brav_lat.latVec(1);
    setMolComCoorVelo(height, imag_target, kinE, deg_incident, deg_azimuth);
}

void SupCell::randSetLatVelo(const double& therm_beta) {
    lat_velo = randn(3, numLatSites()) / sqrt( latMassMat() * therm_beta );
    lat_velo.each_row([](rowvec& r){r -= mean(r);}); // eliminate the overall drift
}

void SupCell::setMolBondLength(const double& l) {
    assert( numMolSites() == 2 );
    vec r01 = molCoor(0) - molCoor(1);
    vec unit = r01 / norm(r01);
    double curr_length = molBondLength();
    mol_coor.col(0) += (l-curr_length) * (mol.mass(1)/mol.totMass()) * unit;
    mol_coor.col(1) -= (l-curr_length) * (mol.mass(0)/mol.totMass()) * unit;
}

double SupCell::kinE() const {
    return latKinE() + molKinE();
}

double SupCell::latKinE() const {
    return 0.5 * accu( latMassMat() % lat_velo % lat_velo );
}

double SupCell::molKinE() const {
    return 0.5 * accu( molMassMat() % mol_velo % mol_velo );
}

double SupCell::molComKinE() const {
    return 0.5 * mol.totMass() * accu( molComVelo() % molComVelo() );
}

double SupCell::molVibKinE() const {
    if ( numMolSites() == 1 )
	return 0;
    assert( numMolSites() == 2 );
    vec r01 = molCoor(0) - molCoor(1);
    vec v01 = molVelo(0) - molVelo(1);
    return 0.5 * mol.reducedMass() * pow( dot( v01, r01/norm(r01) ), 2 );
}

double SupCell::molRotKinE() const {
    if ( numMolSites() == 1 )
	return 0;
    assert( numMolSites() == 2 );
    vec r01 = molCoor(0) - molCoor(1);
    vec v01 = molVelo(0) - molVelo(1);
    return 0.5 * mol.reducedMass() * ( dot( v01, v01 ) - pow( dot( v01, r01/norm(r01) ), 2 ) );
}
