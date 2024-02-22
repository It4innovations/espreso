
#ifndef SRC_WRAPPERS_BEM_W_BEM_H_
#define SRC_WRAPPERS_BEM_W_BEM_H_

namespace espreso {

// V: ne * ne;
// K: ne * np;
// D: np * np;
// M: ne * np;

void BEM3DLaplace   (double *K, int np, double *points, int ne, int *elements, double conductivity);
void BEM3DElasticity(double *K, int np, double *points, int ne, int *elements, double YoungModulus, double PoissonRatio);


void BEM3DLaplaceEval(double *results, int np, double *points, int ne, int *elements, int ni, double *inner, double conductivity, double *dirichlet);
void BEM3DElasticityEval(double *results, int np, double *points, int ne, int *elements, int ni, double *inner, double conductivity, double *dirichlet);



}



#endif /* SRC_WRAPPERS_BEM_W_BEM_H_ */
