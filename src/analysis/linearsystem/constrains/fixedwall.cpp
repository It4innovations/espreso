
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
    std::vector<double> &normal = StructuralMechanics::Results::normal->data;

    if (loadstep.fixed_wall.empty()) {
        return;
    }
    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS || info::mesh->dimension != 3) {
        eslog::error("Fixed wall boundary condition is implemented on for STRUCTURAL_MECHANICS 3D.\n");
    }

    feti.lambdas.intervals.push_back({ 0, 0 });
    int dim = info::mesh->dimension;
    struct __lambda__ { int dof, size, cindex; };
    std::vector<__lambda__> permutation;

    std::vector<std::vector<int> > &D2C = feti.D2C;
    std::vector<std::vector<int> > ROWS(feti.K.size()), COLS(feti.K.size());

    int cc = 0;
    for (auto wall = loadstep.fixed_wall.begin(); wall != loadstep.fixed_wall.end(); ++wall) {
        Point pn, pp;
        pn.x = wall->second.normal.x.evaluator->evaluate();
        pn.y = wall->second.normal.y.evaluator->evaluate();
        pn.z = wall->second.normal.z.evaluator->evaluate();
        pp.x = wall->second.point.x.evaluator->evaluate();
        pp.y = wall->second.point.y.evaluator->evaluate();
        pp.z = wall->second.point.z.evaluator->evaluate();

        const BoundaryRegionStore *region = info::mesh->bregion(wall->first);
        int dindex = 0;
        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n, ++cc) {
            Point nn(normal[*n * dim + 0], normal[*n * dim + 1], normal[*n * dim + 2]);
            Point w = info::mesh->nodes->coordinates->datatarray()[*n] - pp;
            if (pn * w > wall->second.gap || pn * nn > -0.1) { // over the gap or large angle
                continue;
            }
            auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
            for (int d = 0; d < dim; ++d, ++dmap) {
                while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
                if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != *n * dim + d) {
                    permutation.push_back(__lambda__{ (int)(*n * dim + d), (int)dmap->size(), cc });
                }
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

    cindex.clear();
    cindex.resize(feti.K.size());
    std::vector<std::vector<esint> > cperm(feti.K.size());
    for (size_t i = 0, cprev = feti.lambdas.cmap.size(), lprev = 0; i < permutation.size(); ++i) {
        auto dmap = feti.decomposition->dmap->cbegin() + permutation[i].dof;
        bool nextL = i == 0 || lprev != (size_t)permutation[i].dof / 3;
        size_t cbegin = feti.lambdas.cmap.size();

        if (nextL) {
            lprev = permutation[i].dof / 3;
            feti.lambdas.cmap.push_back(1);
            feti.lambdas.cmap.push_back(permutation[i].size);
            if (dmap->at(0).domain < feti.decomposition->dbegin) {
                ++feti.lambdas.intervals.back().halo;
            }
            for (int c = 0; c < permutation[i].size; ++c) {
                feti.lambdas.cmap.push_back(dmap->at(c).domain);
                if (feti.decomposition->ismy(dmap->at(c).domain)) {
                    int dd = dmap->at(c).domain - feti.decomposition->dbegin;
                    cindex[dd].push_back(feti.B1[dd].nrows + ROWS[dd].size());
                    cperm[dd].push_back(permutation[i].cindex);
                    ROWS[dd].push_back(1);
                    COLS[dd].push_back(dmap->at(c).index);
                    D2C [dd].push_back(feti.lambdas.size);
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
                }
            }
        }
    }
    feti.lambdas.intervals.back().size = feti.lambdas.size - feti.lambdas.equalities - feti.lambdas.intervals.back().halo;

    #pragma omp parallel for
    for (size_t d = 0; d < cindex.size(); ++d) {
        std::sort(cindex[d].begin(), cindex[d].end(), [&] (esint i, esint j) {
            return cperm[d][i - feti.B1[d].nrows] < cperm[d][j - feti.B1[d].nrows];
        });
    }

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
    }

    swap(B1, feti.B1);

    feti.c.resize(feti.lambdas.size);
    feti.lb.resize(feti.lambdas.size - feti.lambdas.equalities);
    feti.ub.resize(feti.lambdas.size - feti.lambdas.equalities);
    math::set(feti.lb, T{0});
    math::set(feti.ub, std::numeric_limits<T>::max());
}

template <typename T>
void FixedWall<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);
    std::vector<double> &normal = StructuralMechanics::Results::normal->data;
    std::vector<double> &disp = StructuralMechanics::Results::displacement->data;

    int dim = info::mesh->dimension;
    std::vector<int> cc(cindex.size());
    for (auto wall = loadstep.fixed_wall.begin(); wall != loadstep.fixed_wall.end(); ++wall) {
        Point pn, pp;
        pn.x = wall->second.normal.x.evaluator->evaluate();
        pn.y = wall->second.normal.y.evaluator->evaluate();
        pn.z = wall->second.normal.z.evaluator->evaluate();
        pp.x = wall->second.point.x.evaluator->evaluate();
        pp.y = wall->second.point.y.evaluator->evaluate();
        pp.z = wall->second.point.z.evaluator->evaluate();
        pn.normalize();

        const BoundaryRegionStore *region = info::mesh->bregion(wall->first);
        esint dindex = 0;
        std::vector<double> nv; nv.reserve(3);
        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
            Point nn(normal[*n * dim + 0], normal[*n * dim + 1], normal[*n * dim + 2]);
            Point coo(disp[*n * dim + 0], disp[*n * dim + 1], disp[*n * dim + 2]);
            coo += info::mesh->nodes->coordinates->datatarray()[*n];
            Point w = coo - pp;
            if (pn * w > wall->second.gap || pn * nn > -0.1) { // over the gap or large angle
                continue;
            }
            auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
            double sum = 0;
            nv.clear();
            for (int d = 0; d < dim; ++d, ++dmap) {
                while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
                if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != *n * dim + d) {
                    nv.push_back(nn[d]);
                    sum += dmap->size() * nn[d] * nn[d];
                } else {
                    nn[d] = 0;
                }
            }
            double norm = std::sqrt(sum);
            for (size_t i = 0; i < nv.size(); ++i) {
                nv[i] /= norm;
                nn[i] /= norm;
            }
            if (nv.size()) {
                dmap -= dim;
                int d2c = 0;
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        int dd = di->domain - feti.decomposition->dbegin;
                        int row = cindex[dd][cc[dd]++];
                        for (int c = feti.B1[dd].rows[row], ii = 0; c < feti.B1[dd].rows[row + 1]; ++c, ++ii) {
                            feti.B1[dd].vals[c] = nv[ii];
                        }
                        d2c = feti.D2C[dd][row];
                    }
                }
                feti.c.vals[d2c] = -(pn * w) / (pn * nn);
            }
        }
    }
}

template class FixedWall<double>;

}
