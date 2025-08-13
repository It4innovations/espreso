
#include "fixedtube.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "math/math.h"

#include <climits>
#include <numeric>
#include <fstream>

namespace espreso {


template <typename T>
void FixedTube<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);

    if (loadstep.fixed_tube.empty()) {
        return;
    }
    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS || info::mesh->dimension != 3) {
        eslog::error("Fixed tube boundary condition is implemented only for STRUCTURAL_MECHANICS 3D.\n");
    }

    interval = feti.lambdas.intervals.size();
    cmapsize = feti.lambdas.cmap.size();
    feti.lambdas.intervals.push_back({ 0, 0 });

    dsize.resize(feti.D2C.size());
    for (size_t i = 0; i < dsize.size(); ++i) {
        dsize[i] = feti.D2C[i].size();
    }
}

template <typename T>
void FixedTube<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);
    const std::vector<double> &disp = StructuralMechanics::Results::displacement->data;

    if (loadstep.fixed_tube.empty()) {
        return;
    }
    feti.updated.B = true;

    // we need to reset due to different number of equalities
    feti.lambdas.size = feti.lambdas.equalities;
    for (size_t i = 0; i < dsize.size(); ++i) {
        feti.D2C[i].resize(dsize[i]);
        feti.B1[i].resize(dsize[i], feti.B1[i].ncols, feti.B1[i].rows[dsize[i]]);
    }
    feti.c.resize(feti.lambdas.size);
    feti.lb.clear();
    feti.ub.clear();
    feti.lambdas.intervals[interval] = { 0, 0 };
    feti.lambdas.cmap.resize(cmapsize);

    int dim = info::mesh->dimension;
    struct __lambda__ { esint node, cindex; };
    std::vector<__lambda__> permutation;

    std::vector<std::vector<int> > &D2C = feti.D2C;
    std::vector<std::vector<int> > ROWS(feti.K.size()), COLS(feti.K.size());
    std::vector<std::vector<T> >   VALS(feti.K.size());
    std::vector<T>                 C;

    int index = 0;
    for (auto tube = loadstep.fixed_tube.begin(); tube != loadstep.fixed_tube.end(); ++tube, ++index) {
        Point ul, pp;
        tube->second.direction.x.evaluator->getSubstep() = step.substep;
        tube->second.direction.y.evaluator->getSubstep() = step.substep;
        tube->second.direction.z.evaluator->getSubstep() = step.substep;
        tube->second.center.x.evaluator->getSubstep() = step.substep;
        tube->second.center.y.evaluator->getSubstep() = step.substep;
        tube->second.center.z.evaluator->getSubstep() = step.substep;
        ul.x = tube->second.direction.x.evaluator->evaluate();
        ul.y = tube->second.direction.y.evaluator->evaluate();
        ul.z = tube->second.direction.z.evaluator->evaluate();
        pp.x = tube->second.center.x.evaluator->evaluate();
        pp.y = tube->second.center.y.evaluator->evaluate();
        pp.z = tube->second.center.z.evaluator->evaluate();
        ul.normalize();

        _store(pp, ul, tube->second.radius, index);

        const BoundaryRegionStore *region = info::mesh->bregion(tube->first);
        int dindex = 0;
        _Point<int> free;
        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
            Point nn(region->nodeNormals->data[*n * dim + 0], region->nodeNormals->data[*n * dim + 1], region->nodeNormals->data[*n * dim + 2]);
            Point p = info::mesh->nodes->coordinates->datatarray()[*n] + Point(disp[*n * dim + 0], disp[*n * dim + 1], disp[*n * dim + 2]);
            Point w = p - pp;
            double distance = (w - ul * (w * ul)).length() - tube->second.radius;
            if (distance > tube->second.gap) { // over the gap
                continue;
            }

            double dot = 0;
            { // remove Dirichlet from normal
                auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
                for (int d = 0; d < dim; ++d, ++dmap) {
                    while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
                    if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != *n * dim + d) {
                        free[d] = 1;
                        dot += dmap->size() * nn[d] * nn[d];
                    } else {
                        free[d] = 0;
                        nn[d] = 0;
                    }
                }
            }

            if (free.x || free.y || free.z) { // insert into B1
                // normalize without Dirichlet and for more domains
                nn *= 1 / std::sqrt(dot);
                auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        ROWS[di->domain - feti.decomposition->dbegin].push_back(*n);
                    }
                }
                for (int d = 0; d < dim; ++d, ++dmap) {
                    for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                        if (feti.decomposition->ismy(di->domain)) {
                            COLS[di->domain - feti.decomposition->dbegin].push_back(di->index);
                            VALS[di->domain - feti.decomposition->dbegin].push_back(free[d] ? nn[d] : 0);
                        }
                    }
                }
                permutation.push_back({ *n, (esint)C.size() });
                C.push_back(distance);
            }
        }
    }

    if (permutation.size() == 0) {
        return; // empty
    }

    std::sort(permutation.begin(), permutation.end(), [&] (const __lambda__ &i, const __lambda__ &j) {
        auto imap = feti.decomposition->dmap->cbegin() + dim * i.node;
        auto jmap = feti.decomposition->dmap->cbegin() + dim * j.node;
        int size = std::min(imap->size(), jmap->size());
        for (int k = 0; k < size; ++k) {
            if (imap->at(k).domain != jmap->at(k).domain) { return imap->at(k).domain < jmap->at(k).domain; }
        }
        if (imap->size() != jmap->size()) {
            return imap->size() < jmap->size();
        } else {
            return i.node < j.node;
        }
    });
    for (size_t i = 0, cprev = feti.lambdas.cmap.size(); i < permutation.size(); ++i) {
        auto dmap = feti.decomposition->dmap->cbegin() + dim * permutation[i].node;
        size_t cbegin = feti.lambdas.cmap.size();
        feti.lambdas.cmap.push_back(1);
        feti.lambdas.cmap.push_back(dmap->size());
        for (auto di = dmap->begin(); di != dmap->end(); ++di) {
            feti.lambdas.cmap.push_back(di->domain);
            if (feti.decomposition->ismy(di->domain)) {
                D2C[di->domain - feti.decomposition->dbegin].push_back(feti.lambdas.size);
            }
        }
        if (dmap->front().domain < feti.decomposition->dbegin) {
            ++feti.lambdas.intervals[interval].halo;
        }

        // merge last 2 lambdas together
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
    }

    std::vector<esint> order(info::mesh->nodes->size);
    for (size_t i = 0; i < permutation.size(); ++i) {
        order[permutation[i].node] = i;
    }

    std::vector<Matrix_CSR<T> > B1(feti.B1.size());
    Vector_Dense<T> B1c; B1c.resize(feti.lambdas.size);

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        std::vector<esint> rperm(ROWS[d].size());
        std::iota(rperm.begin(), rperm.end(), 0);
        std::sort(rperm.begin(), rperm.end(), [&] (esint i, esint j) { return order[ROWS[d][i]] < order[ROWS[d][j]]; });

        B1[d].resize(feti.B1[d].nrows + ROWS[d].size(), feti.B1[d].ncols, feti.B1[d].nnz + COLS[d].size());
        memcpy(B1[d].rows, feti.B1[d].rows, sizeof(int) * (feti.B1[d].nrows + 1));
        memcpy(B1[d].cols, feti.B1[d].cols, sizeof(int) * feti.B1[d].nnz);
        memcpy(B1[d].vals, feti.B1[d].vals, sizeof(T)   * feti.B1[d].nnz);
        for (size_t r = 0, ci = feti.B1[d].nnz; r < ROWS[d].size(); ++r) {
            B1[d].rows[feti.B1[d].nrows + r + 1] = B1[d].rows[feti.B1[d].nrows + r] + 3; // always 3 despite some zeros in B1
            for (int c = 0; c < 3; ++c, ++ci) {
                B1[d].cols[ci] = COLS[d][3 * rperm[r] + c];
                B1[d].vals[ci] = VALS[d][3 * rperm[r] + c];
            }
        }
    }
    memcpy(B1c.vals, feti.c.vals, sizeof(T) * feti.c.size);
    for (size_t i = 0, c = feti.c.size; i < permutation.size(); ++i) {
        B1c.vals[c++] = C[permutation[i].cindex];
    }

    swap(B1, feti.B1);
    feti.c.swap(B1c);

    feti.lb.resize(feti.lambdas.size - feti.lambdas.equalities);
    feti.ub.resize(feti.lambdas.size - feti.lambdas.equalities);
    math::set(feti.lb, T{0});
    math::set(feti.ub, std::numeric_limits<T>::max());

    feti.lambdas.intervals[interval].size = feti.lambdas.size - feti.lambdas.equalities - feti.lambdas.intervals.back().halo;
}

template <typename T>
void FixedTube<T>::_store(const Point &center, const Point &direction, double radius, int index)
{
    int layers = 3, steps = 40;

    std::ofstream os(info::ecf->outpath + "/tube." + std::to_string(index) + ".vtk");
    os << "# vtk DataFile Version 2.0\n";
    os << "EXAMPLE\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << layers * steps << " float\n";

    Point n1(1, 0, 0);
    if (Point::cross(direction, n1).length() < 0.1) {
        n1 = Point(0, 1, 0);
    }
    Point n2 = Point::cross(direction, n1);

    for (int i = 0; i < layers; ++i) {
        for (int s = 0; s < steps; ++s) {
            Point p = center + direction * (i - 1) + n1 * std::cos(2 * M_PI * s / steps) * radius + n2 * std::sin(2 * M_PI * s / steps) * radius;
            os << p.x << " " << p.y << " " << p.z << "\n";
        }
    }
    os << "CELLS " << steps * (layers - 1) << " " << 5 * steps * (layers - 1) << " \n";
    for (int i = 0; i < layers - 1; ++i) {
        for (int s = 0; s < steps; ++s) {
            int p0 =  s + 0;
            int p1 = (s + 1) % steps;
            int p2 = steps + p0;
            int p3 = steps + p1;
            os << "4 " << i * steps + p0 << " " << i * steps + p1 << " " << i * steps + p3 << " " << i * steps + p2 << "\n";
        }
    }
    os << "CELL_TYPES " << steps * (layers - 1) << " \n";
    for (int i = 0; i < steps * (layers - 1); ++i) {
        os << "9\n";
    }
}

template class FixedTube<double>;

}



