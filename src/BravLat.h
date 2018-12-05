#ifndef __BRAVAIS_LATTICE_H__
#define __BRAVAIS_LATTICE_H__

#include <armadillo>
#include <vector>

class BravLat
{
    public:
	BravLat();
	BravLat		    (   std::vector<arma::vec>  const&   lat_vec_list_	);

	arma::mat		    const&	latVec() const;
	arma::mat		    const&	rcpVec() const;
	arma::subview_col<double>   const	latVec(const arma::uword& idx) const;
	arma::subview_col<double>   const	rcpVec(const arma::uword& idx) const;
	arma::uword				dim() const;
	void					setLatVec(const arma::uword& idx, const arma::vec& lat_vec);
	void					addLatVec(const arma::vec& lat_vec);

    private:
	arma::mat				lat_vec_list;
	arma::mat		    		rcp_vec_list;
};

arma::mat rcp(const arma::mat& lat_vec_list);

#endif
