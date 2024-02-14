
#ifndef SRC_WRAPPERS_BEM_W_BEM_H_
#define SRC_WRAPPERS_BEM_W_BEM_H_

namespace espreso {

// V: ne * ne;
// K: ne * np;
// D: np * np;
// M: ne * np;

void BEM3dLaplace(int np, double *points, int ne, int *elemNodes, int order, double *V, double *K, double *D, double *M);

//void getLaplace(
//		BEMData* &bem, double *K,
//		esint nNodes, const double *nodes,
//		esint nElements, const esint *elements,
//		double conductivity);
//
//void evaluateLaplace(
//		BEMData* &bem, double *results,
//		esint nNodes, const double *nodes,
//		esint nElements, const esint *elements,
//		esint nPoints, const double *points,
//		double conductivity, double *dirichlet);
//
//void deleteData(BEMData* &bem);

}



#endif /* SRC_WRAPPERS_BEM_W_BEM_H_ */
