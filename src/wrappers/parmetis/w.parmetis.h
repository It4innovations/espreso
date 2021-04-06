
#ifndef SRC_WRAPPERS_METIS_W_PARMETIS_H_
#define SRC_WRAPPERS_METIS_W_PARMETIS_H_

namespace espreso {

struct MPIGroup;

struct ParMETIS {
	static bool islinked();

	enum class METHOD {
		ParMETIS_V3_PartKway,
		ParMETIS_V3_RefineKway,
		ParMETIS_V3_AdaptiveRepart,
		ParMETIS_V3_PartGeom
	};

	static esint call(
			METHOD method,
			MPIGroup &group,
			esint *edistribution,
			esint *eframes, esint *eneighbors,
			esint dimensions, float *coordinates,
			esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
			esint *partition);
};

}



#endif /* SRC_WRAPPERS_METIS_W_PARMETIS_H_ */
