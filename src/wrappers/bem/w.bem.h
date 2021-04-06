
#ifndef SRC_WRAPPERS_BEM_W_BEM_H_
#define SRC_WRAPPERS_BEM_W_BEM_H_

namespace espreso {

struct BEMData;

namespace BEM4I {

bool isLinked();

void getLaplace(
		BEMData* &bem, double *K,
		esint nNodes, const double *nodes,
		esint nElements, const esint *elements,
		double conductivity);

void evaluateLaplace(
		BEMData* &bem, double *results,
		esint nNodes, const double *nodes,
		esint nElements, const esint *elements,
		esint nPoints, const double *points,
		double conductivity, double *dirichlet);

void deleteData(BEMData* &bem);

};

}



#endif /* SRC_WRAPPERS_BEM_W_BEM_H_ */
