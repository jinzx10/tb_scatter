#ifndef __JOIN_H__
#define __JOIN_H__

#include <armadillo>
#include <vector>

template <typename T> arma::Mat<typename T::elem_type> joinCols(const std::vector<T*>&);
template <typename T> arma::Mat<typename T::elem_type> joinCols(const std::vector<T >&);

template <typename T> arma::Mat<typename T::elem_type> joinRows(const std::vector<T*>&);
template <typename T> arma::Mat<typename T::elem_type> joinRows(const std::vector<T >&);

template <typename T> arma::Mat<typename T::elem_type> joinDiag(const std::vector<T*>&);
template <typename T> arma::Mat<typename T::elem_type> joinDiag(const std::vector<T >&);

template <typename T> arma::Mat<typename T::elem_type> join(const std::vector<T*>&, const arma::uword, const arma::uword, const char& = 'c');
template <typename T> arma::Mat<typename T::elem_type> join(const std::vector<T >&, const arma::uword, const arma::uword, const char& = 'c');

#include "join.tpp"

#endif
