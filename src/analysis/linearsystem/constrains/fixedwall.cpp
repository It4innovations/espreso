
#include "fixedwall.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"
#include "math/math.h"

#include <climits>

namespace espreso {

template <typename T>
void FixedWall<T>::set(const step::Step &step, FETI<T> &feti)
{
    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D) {
        eslog::error("Fixed wall boundary condition is implemented on for STRUCTURAL_MECHANICS_3D.\n");
    }

    int dim = info::mesh->dimension;
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics_3d.load_steps_settings.at(step.loadstep + 1);
    std::vector<double> &normal = StructuralMechanics::Results::normal->data;

    std::vector<std::vector<esint> > &D2C = feti.D2C;
    std::vector<std::vector<esint> > COLS(feti.K.size());
    std::vector<std::vector<T> > VALS(feti.K.size());

    for (auto wall = loadstep.fixed_wall.begin(); wall != loadstep.fixed_wall.end(); ++wall) {
        const BoundaryRegionStore *region = info::mesh->bregion(wall->first);
        // const FixedWallConfiguration &conf = wall->second;
        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
            auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
            for (int d = 0; d < dim; ++d, ++dmap) {
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        COLS[di->domain - feti.decomposition->dbegin].push_back(di->index);
                        VALS[di->domain - feti.decomposition->dbegin].push_back(normal[*n * dim + d]);
                    }
                }
            }
        }
    }
    std::vector<Matrix_CSR<T> > B1(feti.B1.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        B1[d].resize(D2C[d].size() + COLS[d].size() / 3, feti.B1[d].ncols, feti.B1[d].nnz + COLS[d].size());
        memcpy(B1[d].rows, feti.B1[d].rows, sizeof(int) * (feti.B1[d].nrows + 1));
        memcpy(B1[d].cols, feti.B1[d].cols, sizeof(int) * feti.B1[d].nnz);
        memcpy(B1[d].vals, feti.B1[d].vals, sizeof(T)   * feti.B1[d].nnz);
        for (size_t i = 0; i < COLS[d].size() / 3; ++i) {
            D2C[d].push_back(feti.lambdas.size++);
            B1[d].rows[feti.B1[d].nrows + i + 1] = B1[d].rows[feti.B1[d].nrows + i] + 3;
        }
        memcpy(B1[d].cols + feti.B1[d].nnz, COLS[d].data(), sizeof(int) * COLS[d].size());
        memcpy(B1[d].vals + feti.B1[d].nnz, VALS[d].data(), sizeof(T)   * VALS[d].size());
    }

    swap(B1, feti.B1);

    feti.c.resize(feti.lambdas.size);
    feti.lb.resize(feti.lambdas.size - feti.lambdas.equalities);
    feti.ub.resize(feti.lambdas.size - feti.lambdas.equalities);
}

template <typename T>
void FixedWall<T>::update(const step::Step &step, FETI<T> &feti)
{
    math::set(feti.lb, T{0});
    math::set(feti.ub, std::numeric_limits<T>::max());
}

template struct FixedWall<double>;

}
