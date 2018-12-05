#ifndef __AUXILLARY_MATH_H__
#define __AUXILLARY_MATH_H__

#include <armadillo>
#include <vector>

// kron(vector<Mat>{A,B,C,...}) returns the chained Kronecker product A*B*C*...
// returns ...*C*B*A if reverse == true
template <typename T> arma::Mat<typename T::elem_type> kron(const std::vector<T*>&, const bool& reverse = false);
template <typename T> arma::Mat<typename T::elem_type> kron(const std::vector<T >&, const bool& reverse = false);

// converts std::vector< std::vector<T> > to arma::Mat<T>
// 'r': turns inner std::vector<T> into rows
// 'c': turns inner std::vector<T> into cols
template <typename T> arma::Mat<T> vv2m(const std::vector<std::vector<T> >&, const char& style = 'r');

// linspace(vn, style) returns vn.size() armadillo row vectors (rowvec) 
// the i-th rowvec contains vn[i] evenly-spaced elements, whose range depends on the style:
// 'u': [ 0 , vn[i]-1 ]
// 's': [ -vn[i]/2 , (vn[i]-1)/2 ] (integer division)
// 'n': range of 's' further divided (through floating point division) by vn[i]
std::vector<arma::rowvec> linspace(const arma::uvec& vn, const char& style = 'n');

// mesh(vn, style) returns a matrix whose columns form a meshgrid
// The number of elements in each dimension is provided by vn
// The range depends on the style option
// which is the same as that prescribed in linspace(vn, style)
// layout = 'f' (Fortran-style) or 'c' (C-style)
arma::mat mesh(const arma::uvec& vn, const char& style = 'n', const char& layout = 'f');

// extCoor(major_coor, minor_coor) returns a horizontally concatenated matrix
// each concatenated unit is a shifted major_coor:
// the original major_coor added column-wise by a column of minor_coor 
arma::mat extCoor(const arma::mat& major_coor, const arma::mat& minor_coor);

// generate a set of lattice grid points according to vv and vn
// vv - lattice vectors
// vn - corresponding number of points
// the size of vv and vn must equal
// style = 'u' or 's'
// layout = 'f' or 'c'
arma::mat bravCoor(const std::vector<arma::vec>& vv, const arma::uvec& vn, const char& style = 'u', const char& layout = 'f');

// TBC...
template <typename T> arma::Mat<typename T::elem_type> extract(const std::vector<T>& vm, const arma::uword& dim);
template <typename T> arma::Row<typename T::elem_type> to_row(const std::vector<T>& vm);

double fermi(const double& E, const double& beta, const double& mu); 
double lorentz(const double& E, const double& mean, const double& width);
double gaussian(const double& E, const double& mean, const double& stddev);

arma::vec fermi(const arma::vec& E, const double& beta, const double& mu); 
arma::vec lorentz(const arma::vec& E, const double& mean, const double& width);
arma::vec gaussian(const arma::vec& E, const double& mean, const double& stddev);

double morse(const double& r, const double& coef, const double& eq_dist, const double& exp_decay_len);
double dMorse(const double& r, const double& coef, const double& eq_dist, const double& exp_decay_len);

// temporarily unavailable...
//template <typename T> T fermi(const T& E, const double& beta, const double& mu);
//template <typename T> T lorentz(const T& E, const double& mean, const double& Gamma);

arma::vec genInUnitCell(const arma::vec& center_coor, const arma::vec& r1, const arma::vec& r2);

#include "auxmath.tpp"

#endif
