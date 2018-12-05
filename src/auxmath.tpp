#include <vector>
#include <armadillo>
#include <cassert>
#include "aux.h"
#include "mathconst.h"

template <typename T>
arma::Mat<typename T::elem_type> kron(const std::vector<T*>& vpT, const bool& reverse) {
    using mat_t = arma::Mat<typename T::elem_type>;
    if (vpT.size() == 0)
	return mat_t{};

    assert( [&](){for (auto& pT : vpT) if (pT == nullptr) return 0; return 1;}() );

    if (vpT.size() == 1)
	return *vpT[0];

    mat_t kronprod = *vpT[0];
    for (auto it = vpT.begin()+1; it != vpT.end(); ++it) {
	if (!reverse)
	    kronprod = arma::kron(kronprod, **it);
	else
	    kronprod = arma::kron(**it, kronprod);
    }
    return kronprod;
}

template <typename T>
arma::Mat<typename T::elem_type> kron(const std::vector<T>& rvT, const bool& reverse) {
    return kron(rv2vp(rvT), reverse);
}

template <typename T>
arma::Mat<T> vv2m(const std::vector< std::vector<T> >& vvT, const char& style) {
    using mat_t = arma::Mat<T>;
    if (vvT.size() == 0)
	return mat_t{};

    assert( [&vvT](){for (auto& vT : vvT) if (vT.size() != vvT[0].size()) return 0; return 1;}() );
    assert( style == 'r' || style == 'c' );

    if (style == 'r') {
	using row_t = arma::Row<T>;
    	mat_t m = arma::zeros<mat_t>(vvT.size(), vvT[0].size());
    	arma::uword idx{0};
    	m.each_row([&idx, &vvT](row_t& row) {row = arma::conv_to<row_t>::from(vvT[idx++]);});
    	return m;
    } else {
	using col_t = arma::Col<T>;
    	mat_t m = arma::zeros<mat_t>(vvT[0].size(), vvT.size());
    	arma::uword idx{0};
    	m.each_col([&idx, &vvT](col_t& col) {col = arma::conv_to<col_t>::from(vvT[idx++]);});
    	return m;
    }
}

template <typename T>
arma::Mat<typename T::elem_type> extract(const std::vector<T>& vm, const arma::uword& dim) {
    using mat_t = arma::Mat<typename T::elem_type>;
    if (vm.size() == 0)
	return mat_t{};

    mat_t m = arma::zeros(vm.size(), vm[0].n_cols);
    for (arma::uword i = 0; i < vm.size(); ++i)
	m.row(i) = vm[i].row(dim);
    return m;
}

//template <typename T>
//T fermi(const T& E, const double& beta, const double& mu) {
//    return 1.0 / ( exp(beta*(E-mu)) + 1 );
//}
//
//template <typename T>
//T lorentz(const T& E, const double& mean, const double& Gamma) {
//    return 1.0/PI*(Gamma/2.0) / ( pow(E-mean, 2) + (Gamma/2.0)*(Gamma/2.0) );
//}
