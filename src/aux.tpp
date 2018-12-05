#include <vector>

template <typename T>
const std::vector<T*> rv2vp(const std::vector<T>& rvT) {
    if (!rvT.size())
	return {};
    std::vector<T*> vptr2T{};
    for (decltype(rvT.size()) i = 0; i != rvT.size(); ++i)
	vptr2T.push_back( const_cast<T*>(&rvT[i]) );
    return vptr2T;
}

template <typename T>
std::vector<T> concatenate(const std::vector<T>& vT1, const std::vector<T>& vT2) {
    auto cc = vT1;
    cc.insert(cc.end(), vT2.begin(), vT2.end());
    return cc;
}
