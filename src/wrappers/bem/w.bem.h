
#ifndef SRC_WRAPPERS_BEM_W_BEM_H_
#define SRC_WRAPPERS_BEM_W_BEM_H_

namespace espreso {

// V: ne * ne;
// K: ne * np;
// D: np * np;
// M: ne * np;

void BEM3DLaplace   (double *K, esint np, double *points, esint ne, esint *elements, double conductivity);
void BEM3DElasticity(double *K, esint np, double *points, esint ne, esint *elements, double YoungModulus, double PoissonRatio);


void BEM3DLaplaceEval(double *results, esint np, double *points, esint ne, esint *elements, esint ni, double *inner, double conductivity, double *dirichlet);
void BEM3DElasticityEval(double *results, esint np, double *points, esint ne, esint *elements, esint ni, double *inner, double YoungModulus, double PoissonRatio, double *dirichlet);



}



#endif /* SRC_WRAPPERS_BEM_W_BEM_H_ */
