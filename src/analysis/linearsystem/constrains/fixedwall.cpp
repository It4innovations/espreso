
#include "fixedwall.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"
#include "math/math.h"

#include <climits>

namespace espreso {

template <typename T>
void FixedWall<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);

    if (loadstep.fixed_wall.empty()) {
        return;
    }
    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS || info::mesh->dimension != 3) {
        eslog::error("Fixed wall boundary condition is implemented on for STRUCTURAL_MECHANICS 3D.\n");
    }

    int dim = info::mesh->dimension;
    struct __lambda__ { int dof, size; double value; };
    std::vector<__lambda__> permutation;

    std::vector<double> &normal = StructuralMechanics::Results::normal->data;

    std::vector<std::vector<esint> > &D2C = feti.D2C;
    std::vector<std::vector<esint> > ROWS(feti.K.size()), COLS(feti.K.size());
    std::vector<std::vector<T> > VALS(feti.K.size());

    for (auto wall = loadstep.fixed_wall.begin(); wall != loadstep.fixed_wall.end(); ++wall) {
        const BoundaryRegionStore *region = info::mesh->bregion(wall->first);
        esint dindex = 0;
        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
            auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
            size_t psize = permutation.size();
            double sum = 0;
            for (int d = 0; d < dim; ++d, ++dmap) {
                while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
                if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != *n * dim + d) {
                    permutation.push_back(__lambda__{ *n * dim + d, (int)dmap->size(), normal[*n * dim + d] });
                    sum += dmap->size() * normal[*n * dim + d] * normal[*n * dim + d];
                }
            }
            double norm = std::sqrt(sum);
            for (size_t pp = psize; pp < permutation.size(); ++pp) {
                permutation[pp].value /= norm;
            }
        }
    }

    std::sort(permutation.begin(), permutation.end(), [&] (const __lambda__ &i, const __lambda__ &j) {
        auto imap = feti.decomposition->dmap->cbegin() + i.dof;
        auto jmap = feti.decomposition->dmap->cbegin() + j.dof;
        int size = std::min(i.size, j.size);
        for (int k = 0; k < size; ++k) {
            if (imap->at(k).domain != jmap->at(k).domain) { return imap->at(k).domain < jmap->at(k).domain; }
        }
        if (i.size != j.size) {
            return i.size < j.size;
        } else {
            return i.dof < j.dof;
        }
    });

    feti.lambdas.nc_halo = 0;
    for (size_t i = 0, cprev = feti.lambdas.cmap.size(), lprev = 0; i < permutation.size(); ++i) {
        auto dmap = feti.decomposition->dmap->cbegin() + permutation[i].dof;
        bool nextL = i == 0 || lprev != (size_t)permutation[i].dof / 3;
        size_t cbegin = feti.lambdas.cmap.size();

        if (nextL) {
            lprev = permutation[i].dof / 3;
            feti.lambdas.cmap.push_back(1);
            feti.lambdas.cmap.push_back(permutation[i].size);
            if (dmap->at(0).domain < feti.decomposition->dbegin) {
                ++feti.lambdas.nc_halo;
            }
            for (int c = 0; c < permutation[i].size; ++c) {
                feti.lambdas.cmap.push_back(dmap->at(c).domain);
                if (feti.decomposition->ismy(dmap->at(c).domain)) {
                    ROWS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(1);
                    COLS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(dmap->at(c).index);
                    VALS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(permutation[i].value);
                    D2C [dmap->at(c).domain - feti.decomposition->dbegin].push_back(feti.lambdas.size);
                }
            }
            if (feti.lambdas.cmap[cprev + 1] == feti.lambdas.cmap[cbegin + 1]) { // domains sizes
                if (std::equal(feti.lambdas.cmap.cbegin() + cbegin + 1, feti.lambdas.cmap.cend(), feti.lambdas.cmap.cbegin() + cprev + 1)) {
                    if (cprev != cbegin) {
                        ++feti.lambdas.cmap[cprev];
                        feti.lambdas.cmap.resize(cbegin);
                    }
                } else {
                    cprev = cbegin;
                }
            } else {
                cprev = cbegin;
            }
            feti.lambdas.size++;
        } else {
            for (int c = 0; c < permutation[i].size; ++c) {
                if (feti.decomposition->ismy(dmap->at(c).domain)) {
                    ROWS[dmap->at(c).domain - feti.decomposition->dbegin].back()++;
                    COLS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(dmap->at(c).index);
                    VALS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(permutation[i].value);
                }
            }
        }
    }
    feti.lambdas.nc_size = feti.lambdas.size - feti.lambdas.equalities - feti.lambdas.nc_halo;

    std::vector<Matrix_CSR<T> > B1(feti.B1.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        B1[d].resize(feti.B1[d].nrows + COLS[d].size() / 3, feti.B1[d].ncols, feti.B1[d].nnz + COLS[d].size());
        memcpy(B1[d].rows, feti.B1[d].rows, sizeof(int) * (feti.B1[d].nrows + 1));
        memcpy(B1[d].cols, feti.B1[d].cols, sizeof(int) * feti.B1[d].nnz);
        memcpy(B1[d].vals, feti.B1[d].vals, sizeof(T)   * feti.B1[d].nnz);
        for (size_t i = 0; i < ROWS[d].size(); ++i) {
            B1[d].rows[feti.B1[d].nrows + i + 1] = B1[d].rows[feti.B1[d].nrows + i] + ROWS[d][i];
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
