#ifndef __AUXILLARY_H__
#define __AUXILLARY_H__

#include <vector>

template<typename T> const std::vector<T*> rv2vp(const std::vector<T>&);
template<typename T> std::vector<T> concatenate(const std::vector<T>&, const std::vector<T>&);

#include "aux.tpp"

#endif
