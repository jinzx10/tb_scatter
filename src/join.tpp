#include <vector>
#include <armadillo>
#include <cassert>
#include "aux.h"

template <typename T>
arma::Mat<typename T::elem_type> joinCols(const std::vector<T*>& vpT) {
    using joined_t = arma::Mat<typename T::elem_type>;
    if (vpT.size() == 0)
	return joined_t{};

    assert( [&](){for (auto& pT : vpT) if (pT == nullptr) return 0; return 1;}() &&
	    [&](){for (auto& pT : vpT) if (pT->n_cols != vpT[0]->n_cols) return 0; return 1;}() );

    if (vpT.size() == 1)
	return *vpT[0];

    arma::uword nRow{0};
    for (auto& pT : vpT)
	nRow += pT->n_rows;

    joined_t joined = arma::zeros<joined_t>(nRow, vpT[0]->n_cols);
    arma::uword rowLabel{0};
    for (auto it = vpT.begin(); it != vpT.end(); ++it) {
	joined.rows(rowLabel, rowLabel+(*it)->n_rows-1) = **it;
	rowLabel += (*it)->n_rows;
    }
    return joined;
}

template <typename T>
arma::Mat<typename T::elem_type> joinRows(const std::vector<T*>& vpT) {
    using joined_t = arma::Mat<typename T::elem_type>;
    if (vpT.size() == 0)
	return joined_t{};

    assert( [&](){for (auto& pT : vpT) if (pT == nullptr) return 0; return 1;}() &&
	    [&](){for (auto& pT : vpT) if (pT->n_rows != vpT[0]->n_rows) return 0; return 1;}() );

    if (vpT.size() == 1)
	return *vpT[0];

    arma::uword nCol{0};
    for (auto& pT : vpT)
	nCol += pT->n_cols;

    joined_t joined = arma::zeros<joined_t>(vpT[0]->n_rows, nCol);
    arma::uword colLabel{0};
    for (auto it = vpT.begin(); it != vpT.end(); ++it) {
	joined.cols(colLabel, colLabel+(*it)->n_cols-1) = **it;
	colLabel += (*it)->n_cols;
    }
    return joined;
}

template <typename T>
arma::Mat<typename T::elem_type> joinDiag(const std::vector<T*>& vpT) {
    using joined_t = arma::Mat<typename T::elem_type>;
    if (vpT.size() == 0)
	return joined_t{};
    assert( [&](){for (auto& pT : vpT) if (pT == nullptr) return 0; return 1;}() );

    if (vpT.size() == 1)
	return *vpT[0];

    arma::uword nRow{0}, nCol{0};
    for (auto& pT : vpT) {
	nRow += pT->n_rows;
	nCol += pT->n_cols;
    }

    joined_t joined = arma::zeros<joined_t>(nRow, nCol);
    arma::uword rowLabel{0}, colLabel{0};
    for (auto it = vpT.begin(); it != vpT.end(); ++it) {
	if ( (*it)->n_rows != 0 && (*it)->n_cols != 0 )
	    joined( arma::span(rowLabel, rowLabel+(*it)->n_rows-1),
		    arma::span(colLabel, colLabel+(*it)->n_cols-1) ) = **it;
	rowLabel += (*it)->n_rows;
	colLabel += (*it)->n_cols;
    }
    return joined;
}

template <typename T>
arma::Mat<typename T::elem_type> join(const std::vector<T*>& vpT, const arma::uword nRowBlock, const arma::uword nColBlock, const char& major) {
    using joined_t = arma::Mat<typename T::elem_type>;
    if (vpT.size() == 0)
	return joined_t{};
    assert( [&](){for (auto& pT : vpT) if (pT == nullptr) return 0; return 1;}() &&
	    vpT.size() == nRowBlock*nColBlock &&
	    (major == 'c' || major == 'r') );

    if (vpT.size() == 1)
	return *vpT[0];

    arma::uword nAlongMajor{0}, nAlongMinor{0}, *pnRow(&nAlongMajor), *pnCol(&nAlongMinor);
    arma::uword nBlockAlongMajor{nRowBlock}, dimAlongMajor{0}, dimAlongMinor{1};
    arma::uword rowIdx{0}, colIdx{0}, *pIdxLabelMajor{&rowIdx}, *pIdxLabelMinor{&colIdx};
    if (major == 'r') {
	nBlockAlongMajor = nColBlock;
	dimAlongMajor = 1;
	dimAlongMinor = 0;
	pnRow = &nAlongMinor;
	pnCol = &nAlongMajor;
	pIdxLabelMajor = &colIdx;
	pIdxLabelMinor = &rowIdx;
    }
    
    for (arma::uword blockIdx = 0; blockIdx < vpT.size(); blockIdx += nBlockAlongMajor) {
	std::vector<T*> majorBlock(vpT.begin()+blockIdx, vpT.begin()+blockIdx+nBlockAlongMajor);
	assert( [&]() {for (auto& p : majorBlock) if (arma::size(*p, dimAlongMinor) != arma::size(*majorBlock[0], dimAlongMinor)) return 0; return 1; }() );
	nAlongMinor += arma::size(*majorBlock[0], dimAlongMinor);

	arma::uword nAlongMajor_tmp{0};
	for (auto& pT : majorBlock)
	    nAlongMajor_tmp += arma::size(*pT, dimAlongMajor);
	if (blockIdx == 0) 
	    nAlongMajor = nAlongMajor_tmp;
	assert(nAlongMajor_tmp == nAlongMajor);
    }

    joined_t joined = arma::zeros<joined_t>(*pnRow, *pnCol);
    for (auto it = vpT.begin(); it != vpT.end(); ++it) {
	joined( arma::span(rowIdx, rowIdx+(*it)->n_rows-1),
		arma::span(colIdx, colIdx+(*it)->n_cols-1)) = **it;
	*pIdxLabelMajor += arma::size(**it, dimAlongMajor);
	*pIdxLabelMinor += (*pIdxLabelMajor/nAlongMajor) * arma::size(**it, dimAlongMinor);
	*pIdxLabelMajor %= nAlongMajor;
    }
    return joined;
}

template <typename T>
arma::Mat<typename T::elem_type> joinCols(const std::vector<T>& rvT) {
    return joinCols(rv2vp(rvT));
}

template <typename T>
arma::Mat<typename T::elem_type> joinRows(const std::vector<T>& rvT) {
    return joinRows(rv2vp(rvT));
}

template <typename T>
arma::Mat<typename T::elem_type> joinDiag(const std::vector<T>& rvT) {
    return joinDiag(rv2vp(rvT));
}

template <typename T>
arma::Mat<typename T::elem_type> join(const std::vector<T>& rvT, const arma::uword nRowBlock, const arma::uword nColBlock, const char& major) {
    return join(rv2vp(rvT), nRowBlock, nColBlock, major);
}
