#include <cassert>
#include "auxmath.h"
#include "join.h"
#include "mathconst.h"

using namespace arma;
using namespace std;

vector<rowvec> linspace(const uvec& vn, const char& style) {
    assert( style == 'u' || style == 's' || style == 'n');
    if (vn.size() == 0)
	return {};
    vector<rowvec> ls(vn.size(), rowvec{});
    for (uword i = 0; i < vn.size(); ++i) {
	if (vn(i) == 0) {
	    ls[i] = rowvec{};
	    continue;
	}
	switch(style) {
	    case 'u': {
		ls[i] = linspace<rowvec>(0, vn(i)-1, vn(i));
		break;
	    }
	    case 's': {
		ls[i] = linspace<rowvec>(-static_cast<sword>(vn(i)/2), (vn(i)-1)/2, vn(i));
		break;
	    }
	    default:  {
		ls[i] = linspace<rowvec>(-static_cast<sword>(vn(i)/2), (vn(i)-1)/2, vn(i));
		ls[i] /= vn(i);
	    }
	}
    }
    return ls;
}

mat mesh(const uvec& vn, const char& style, const char& layout) {
    assert( style == 'u' || style == 's' || style == 'n' );
    assert( layout == 'f' || layout == 'c' );
    if (vn.size() == 0)
	return mat{};
    vector<rowvec> ls = linspace(vn, style);
    if (vn.size() == 1)
	return ls[0];

    vector<mat> rows(vn.size(), mat{});
    for (uword i = 0; i < vn.size(); ++i) {
	vector<rowvec> vr(vn.size(), rowvec{});
	for (uword j = 0; j < vn.size(); ++j) {
	    if ( i == j ) {
		vr[j] = ls[j];
		continue;
	    }
	    vr[j] = ones<rowvec>(vn(j));
	}
	if (layout == 'f')
	    rows[i] = kron(vr, true);
	else
	    rows[i] = kron(vr);
    }
    return joinCols(rows);
}

mat extCoor(const mat& major_coor, const mat& minor_coor) {
    assert( major_coor.n_rows == minor_coor.n_rows );
    return kron(ones(1,minor_coor.n_cols), major_coor) + kron(minor_coor, ones(1, major_coor.n_cols));
}

mat bravCoor(const vector<vec>& vv, const uvec& vn, const char& style, const char& layout) {
    assert( vv.size() == vn.size() );
    assert( style == 'u' || style == 's' );
    assert( layout == 'f' || layout == 'c' );
    return joinRows(vv) * mesh(vn, style, layout);
}

vec fermi(const vec& E, const double& beta, const double& mu) {
    if (E.size() == 0)
	return {};
    vec f = 1.0 / ( exp(beta*(E-mu)) + 1.0 );
    f.for_each([](double& val){if (val < EPS) val = 0;});
    return f;
}

vec lorentz(const vec& E, const double& mean, const double& width) {
    if (E.size() == 0)
	return {};
    vec l = 1.0/PI*(width/2.0) / ( (E-mean)%(E-mean) + (width/2.0)*(width/2.0) );
    l.for_each([](double& val){if (val < EPS) val = 0;});
    return l;
}

vec gaussian(const vec& E, const double& mean, const double& stddev) {
    if (E.size() == 0)
	return {};
    vec g = 1.0 / stddev / sqrt(2.0*PI) * exp( -(E-mean) % (E-mean) / 2.0 / stddev / stddev );
    g.for_each([](double& val){if (val < EPS) val = 0;});
    return g;
}

double fermi(const double& E, const double& beta, const double& mu) {
    double f = 1.0 / ( exp(beta*(E-mu)) + 1.0 );
    return (f < EPS) ? 0.0 : f;
}


double lorentz(const double& E, const double& mean, const double& width) {
    double l = 1.0/PI*(width/2.0) / ( (E-mean)*(E-mean) + (width/2.0)*(width/2.0) );
    return (l < EPS) ? 0.0 : l;
}

double gaussian(const double& E, const double& mean, const double& stddev) {
    double g = 1.0 / stddev / sqrt(2.0*PI) * exp( -(E-mean) * (E-mean) / 2.0 / stddev / stddev );
    return (g < EPS) ? 0.0 : g;
}

double morse(const double& r, const double& coef, const double& eq_dist, const double& exp_decay_len) {
    return coef * pow( 1.0 - exp(-(r-eq_dist)/exp_decay_len), 2 );
}

double dMorse(const double& r, const double& coef, const double& eq_dist, const double& exp_decay_len) {
    double x = exp( - (r-eq_dist) / exp_decay_len );
    return 2.0 * coef * (1.0 - x) / exp_decay_len * x;
}

vec genInUnitCell(const vec& center_coor, const vec& r1, const vec& r2) {
    return r1 * (randu()-0.5) + r2 * (randu() - 0.5) + center_coor;
}
