
#ifndef BEM_H_
#define BEM_H_

namespace bem{

void getLameSteklovPoincare(
    double * Sarray,
    int nNodes,
    const double * nodes,
    int nElems,
    const int * elems,
    double nu,
    double E,
    int orderNear,
    int orderFar,
    bool verbose = false
    );

}

#endif
