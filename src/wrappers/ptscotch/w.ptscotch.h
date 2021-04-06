
#ifndef SRC_WRAPPERS_SCOTCH_W_PTSCOTCH_H_
#define SRC_WRAPPERS_SCOTCH_W_PTSCOTCH_H_

namespace espreso {

struct MPIGroup;

struct PTScotch {
	static bool islinked();

	static esint call(
			MPIGroup &group,
			esint *edistribution,
			esint *eframes, esint *eneighbors,
			esint dimensions, float *coordinates,
			esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
			esint *partition);
};

}



#endif /* SRC_WRAPPERS_SCOTCH_W_PTSCOTCH_H_ */
