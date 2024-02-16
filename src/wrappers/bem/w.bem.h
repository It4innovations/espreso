
#ifndef SRC_WRAPPERS_BEM_W_BEM_H_
#define SRC_WRAPPERS_BEM_W_BEM_H_

namespace espreso {

// V: ne * ne;
// K: ne * np;
// D: np * np;
// M: ne * np;

void BEM3DLaplace   (int np, double *points, int ne, int *elemNodes, int order, double *V, double *K, double *D, double *M);
void BEM3DElasticity(int np, double *points, int ne, int *elemNodes, int order, double YoungModulus, double PoissonRatio, double *K);

}



#endif /* SRC_WRAPPERS_BEM_W_BEM_H_ */
