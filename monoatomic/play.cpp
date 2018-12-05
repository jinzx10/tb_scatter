#include "join.h"
#include <armadillo>

using namespace arma;
using namespace std;

int main()
{
    vector<mat> vm(4, ones(2,2));
    for (uword i = 0; i < 4; ++i) {
	vm[i] *= i;
    }
    vm[2] = {{1,2},{3,4}};
    auto m = join(vector<mat>{vm[0], vm[2], vm[2].t(), vm[1]}, 2, 2, 'c');
    auto n = join(vector<mat>{vm[0], vm[2], vm[2].t(), vm[1]}, 2, 2, 'r');
    m.print();
    n.print();
    return 0;
}
