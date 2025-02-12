
#ifndef SRC_WRAPPERS_SCOTCH_W_SCOTCH_H_
#define SRC_WRAPPERS_SCOTCH_W_SCOTCH_H_

namespace espreso {

struct ScotchConfiguration;

struct Scotch {
    static bool islinked();

    static esint call(
            const ScotchConfiguration &options,
            esint verticesCount,
            esint *eframes, esint *eneighbors,
            esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
            esint parts, esint *partition);
};

}



#endif /* SRC_WRAPPERS_SCOTCH_W_SCOTCH_H_ */
