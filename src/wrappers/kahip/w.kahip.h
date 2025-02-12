
#ifndef SRC_WRAPPERS_KAHIP_W_KAHIP_H_
#define SRC_WRAPPERS_KAHIP_W_KAHIP_H_

namespace espreso {

struct KaHIPConfiguration;

struct KaHIP {
    static bool islinked();

    static esint call(
        const KaHIPConfiguration &options,
        esint verticesCount,
        esint *eframes, esint *eneighbors,
        esint verticesWeightCount, esint *verticesWeights, esint *edgeWeights,
        esint parts, esint *partition);
};

}

#endif /* SRC_WRAPPERS_KAHIP_W_KAHIP_H_ */
