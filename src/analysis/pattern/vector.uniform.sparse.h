
#ifndef SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_SPARSE_H_
#define SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_SPARSE_H_

#include "decomposition.direct.h"
#include "synchronization.h"
#include "analysis/math/vector_distributed.h"
#include "esinfo/ecfinfo.h"

#include <vector>

namespace espreso {

// Dirichlet
struct VectorUniformSparse {

    template <typename T>
    struct Sync: public Synchronization<Vector_Distributed<Vector_Sparse, T> > {
        void gatherFromUpper(Vector_Distributed<Vector_Sparse, T> &o) {}
        void scatterToUpper(Vector_Distributed<Vector_Sparse, T> &o) {}
    }; // no need for synchronization

    VectorUniformSparse(HeatTransferLoadStepConfiguration &configuration, int multiplicity);
    VectorUniformSparse(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity);

    template <typename T>
    void set(DecompositionDirect &decomposition, Vector_Distributed<Vector_Sparse, T> *v)
    {
        v->decomposition = &decomposition;
        v->cluster.resize(decomposition.end - decomposition.begin, indices.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            v->cluster.indices[i] = indices[i]; // TODO: shallow copy of indices?
        }
    }

    template <typename T>
    void map(Vector_Distributed<Vector_Sparse, T> *v)
    {
        v->mapping.boundary.resize(dirichlet.size());
        for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
            v->mapping.boundary[r].resize(dirichlet[r].offset.size());
            for (size_t t = 0; t < dirichlet[r].offset.size(); ++t) {
                v->mapping.boundary[r][t].filter = dirichlet[r].filter;
                v->mapping.boundary[r][t].data = v->cluster.vals;
                v->mapping.boundary[r][t].position = dirichlet[r].permutation.data() + dirichlet[r].offset[t];
            }
        }
    }

protected:
    struct RegionInfo {
        int filter;
        std::vector<esint> permutation;
        std::vector<esint> offset;
    };

    std::vector<esint> indices;
    std::vector<RegionInfo> dirichlet;

private:
    void buildPattern(std::vector<int> &isset);
};

}




#endif /* SRC_ANALYSIS_PATTERN_VECTOR_UNIFORM_SPARSE_H_ */
