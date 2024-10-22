
#ifndef SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_DENSE_H_
#define SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_DENSE_H_

#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"

#include <vector>

namespace espreso {

struct VectorUniformDense {

    template <typename T>
    struct Sync: public Synchronization<Vector_Distributed<Vector_Dense, T> > {
        std::vector<std::vector<T> > sBuffer, rBuffer;
        std::vector<std::vector<esint> > rOffset;
        std::vector<esint> nOffset;
        std::vector<int> neighbors;

        void init(DecompositionDirect &decomposition);
        void gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v);
        void scatterToUpper(Vector_Distributed<Vector_Dense, T> &v);
    };

    VectorUniformDense(HeatTransferLoadStepConfiguration &configuration, int multiplicity);
    VectorUniformDense(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity);

    template <typename T>
    void set(DecompositionDirect &decomposition, Vector_Distributed<Vector_Dense, T> *v)
    {
        v->decomposition = &decomposition;
        v->cluster.resize(decomposition.halo.size() + decomposition.end - decomposition.begin);
    }

    template <typename T>
    void map(Vector_Distributed<Vector_Dense, T> *v)
    {
        v->mapping.elements.resize(elements.offset.size());
        for (size_t i = 0; i < elements.offset.size(); ++i) {
            v->mapping.elements[i].data = v->cluster.vals;
            v->mapping.elements[i].position = elements.permutation.data() + elements.offset[i];
        }

        v->mapping.boundary.resize(boundary.size());
        for (size_t r = 0; r < boundary.size(); ++r) {
            v->mapping.boundary[r].resize(boundary[r].offset.size());
            for (size_t i = 0; i < boundary[r].offset.size(); ++i) {
                v->mapping.boundary[r][i].data = v->cluster.vals;
                v->mapping.boundary[r][i].position = boundary[r].permutation.data() + boundary[r].offset[i];
            }
        }
    }

protected:
    struct RegionInfo {
        std::vector<esint> permutation;
        std::vector<esint> offset;
    };

    int dofs;
    RegionInfo elements;
    std::vector<RegionInfo> boundary;

private:
    void buildPattern(int dofs);
};

}

#endif /* SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_DENSE_H_ */
