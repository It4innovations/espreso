
#include "fixedsphere.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "math/math.h"

#include <climits>
#include <numeric>
#include <fstream>

namespace espreso {


template <typename T>
void FixedSphere<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);

    if (loadstep.fixed_sphere.empty()) {
        return;
    }
    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS || info::mesh->dimension != 3) {
        eslog::error("Fixed sphere boundary condition is implemented only for STRUCTURAL_MECHANICS 3D.\n");
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
void FixedSphere<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    StructuralMechanicsLoadStepConfiguration &loadstep = info::ecf->structural_mechanics.load_steps_settings.at(step.loadstep + 1);
    const std::vector<double> &disp = StructuralMechanics::Results::displacement->data;

    if (loadstep.fixed_sphere.empty()) {
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

    for (auto sphere = loadstep.fixed_sphere.begin(); sphere != loadstep.fixed_sphere.end(); ++sphere) {
        Point center(sphere->second.center.x.evaluator->evaluate(), sphere->second.center.y.evaluator->evaluate(), sphere->second.center.z.evaluator->evaluate());
        double radius = sphere->second.radius;

//        std::vector<Point> coos;
//        std::vector<std::pair<Point, double> > gaps;
//        _store(center, radius);

        const BoundaryRegionStore *region = info::mesh->bregion(sphere->first);
        int dindex = 0;

        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
            Point coo = info::mesh->nodes->coordinates->datatarray()[*n] + Point(disp[*n * dim + 0], disp[*n * dim + 1], disp[*n * dim + 2]);
            Point direction = center - coo;
            if (direction.length() > sphere->second.gap) {
                continue;
            }

            bool have_dirichlet = 0;
            double dot = 0;
            {
                auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
                for (int d = 0; d < dim; ++d, ++dmap) {
                    dot += dmap->size() * direction[d] * direction[d];
                    while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
                    if (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] == *n * dim + d) {
                        have_dirichlet = true;
                    }
                }
            }
            if (have_dirichlet) {
                continue;
            }
            direction *= 1 / std::sqrt(dot);
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
                        VALS[di->domain - feti.decomposition->dbegin].push_back(direction[d]);
                    }
                }
            }
            permutation.push_back({ *n, (esint)C.size() });
            C.push_back((center - coo).length() - radius);
//            coos.push_back(coo);
//            gaps.push_back(std::make_pair(direction, C.back()));
        }
//        _store(coos, gaps);
    }

    // CODE FOR COMPUTE GAP FROM NORMALS
//    for (auto sphere = loadstep.fixed_sphere.begin(); sphere != loadstep.fixed_sphere.end(); ++sphere) {
//        Point center(sphere->second.center.x.evaluator->evaluate(), sphere->second.center.y.evaluator->evaluate(), sphere->second.center.z.evaluator->evaluate());
//        double radius = sphere->second.radius;
//
//        _store(center, radius);
//
//        const BoundaryRegionStore *region = info::mesh->bregion(sphere->first);
//        int dindex = 0;
//        _Point<int> free;
//        //  center -- X
//        //     \      |
//        //      \     |
//        //       \    ^
//        //        \   | nn
//        //         \  |
//        //          coo
//
//        for (auto n = region->nodes->datatarray().begin(); n < region->nodes->datatarray().end(); ++n) {
//            Point nn(normal[*n * dim + 0], normal[*n * dim + 1], normal[*n * dim + 2]);
//            Point coo = info::mesh->nodes->coordinates->datatarray()[*n] + Point(disp[*n * dim + 0], disp[*n * dim + 1], disp[*n * dim + 2]);
//            Point X = coo + nn * ((center - coo) * nn);
//            if ((X - center).length() >= radius && (coo - center).length() > sphere->second.gap) { // is crossing or over the gap
//                continue;
//            }
//
//            double dot = 0;
//            { // remove Dirichlet from normal
//                auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
//                for (int d = 0; d < dim; ++d, ++dmap) {
//                    while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < *n * dim + d) { ++dindex; }
//                    if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != *n * dim + d) {
//                        free[d] = 1;
//                        dot += dmap->size() * nn[d] * nn[d];
//                    } else {
//                        eslog::error("implement removing Dirichlet from FixedSphere gap.\n");
//                        // TODO: how to do it?
//                        free[d] = 0;
//                        nn[d] = 0;
//                    }
//                }
//            }
//
//            if (free.x || free.y || free.z) { // insert into B1
//                // normalize without Dirichlet and for more domains
//                nn *= 1 / std::sqrt(dot);
//                auto dmap = feti.decomposition->dmap->cbegin() + *n * dim;
//                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
//                    if (feti.decomposition->ismy(di->domain)) {
//                        ROWS[di->domain - feti.decomposition->dbegin].push_back(*n);
//                    }
//                }
//                for (int d = 0; d < dim; ++d, ++dmap) {
//                    for (auto di = dmap->begin(); di != dmap->end(); ++di) {
//                        if (feti.decomposition->ismy(di->domain)) {
//                            COLS[di->domain - feti.decomposition->dbegin].push_back(di->index);
//                            VALS[di->domain - feti.decomposition->dbegin].push_back(free[d] ? nn[d] : 0);
//                        }
//                    }
//                }
//                permutation.push_back({ *n, (esint)C.size() });
//                C.push_back((X - coo).length());
//                printf("C: %+.4e\n", C.back());
//            }
//        }
//    }

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

    feti.lambdas.intervals.back().size = feti.lambdas.size - feti.lambdas.equalities - feti.lambdas.intervals.back().halo;
}

template <typename T>
void FixedSphere<T>::_store(const Point &center, double radius)
{
    int layers = 20, steps = 30;

    std::ofstream os(info::ecf->outpath + "/sphere.vtk");
    os << "# vtk DataFile Version 2.0\n";
    os << "EXAMPLE\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << layers * steps + 2 << " float\n";

    Point p(center.x, center.y, center.z - radius);
    os << p.x << " " << p.y << " " << p.z << "\n";
    for (int i = 1; i <= layers; ++i) {
        p.z = center.z + radius * std::sin(M_PI * i / (layers + 1) - M_PI_2);
        double scale = std::sqrt(radius * radius - (center.z - p.z) * (center.z - p.z));
        for (int s = 0; s < steps; ++s) {
            p.x = center.x + scale * std::cos(2 * M_PI * s / steps);
            p.y = center.y + scale * std::sin(2 * M_PI * s / steps);
            os << p.x << " " << p.y << " " << p.z << "\n";
        }
    }
    os << center.x << " " << center.y << " " << center.z  + radius << "\n";

    os << "CELLS " << steps * (layers - 1) + 2 * steps << " " << 5 * steps * (layers - 1) + 8 * steps << " \n";
    for (int s = 0; s < steps; ++s) {
        int p0 =  s + 0;
        int p1 = (s + 1) % steps;
        os << "3 " << 0 << " " << 1 + p0 << " " << 1 + p1 << "\n";
    }
    for (int i = 0; i < layers - 1; ++i) {
        for (int s = 0; s < steps; ++s) {
            int p0 =  s + 0;
            int p1 = (s + 1) % steps;
            int p2 = steps + p0;
            int p3 = steps + p1;
            os << "4 " << 1 + i * steps + p0 << " " << 1 + i * steps + p1 << " " << 1 + i * steps + p3 << " " << 1 + i * steps + p2 << "\n";
        }
    }
    for (int s = 0; s < steps; ++s) {
        int p0 =  s + 0;
        int p1 = (s + 1) % steps;
        os << "3 " << layers * steps + 2 - 1 << " " << 1 + (layers - 1) * steps + p0 << " " << 1 + (layers - 1) * steps + p1 << "\n";
    }
    os << "CELL_TYPES " << steps * (layers - 1) + 2 * steps << " \n";
    for (int i = 0; i < steps; ++i) {
        os << "5\n";
    }
    for (int i = 0; i < steps * (layers - 1); ++i) {
        os << "9\n";
    }
    for (int i = 0; i < steps; ++i) {
        os << "5\n";
    }
}

template <typename T>
void FixedSphere<T>::_store(const std::vector<Point> &coords, const std::vector<std::pair<Point, double> > &gaps)
{
    std::ofstream os(info::ecf->outpath + "/gap.vtk");
    os << "# vtk DataFile Version 2.0\n";
    os << "EXAMPLE\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
    os << "POINTS " << 2 * gaps.size() << " float\n";
    for (size_t i = 0; i < gaps.size(); ++i) {
        Point from = coords[i];
        Point to = from + (gaps[i].first * gaps[i].second);
        os << from.x << " " << from.y << " " << from.z << "\n";
        os << to.x << " " << to.y << " " << to.z << "\n";
    }
    os << "CELLS " << gaps.size() << " " << 3 * gaps.size() << "\n";
    for (size_t i = 0; i < gaps.size(); ++i) {
        os << 2 << " " << 2 * i << " " << 2 * i + 1 << "\n";
    }
    os << "CELL_TYPES " << gaps.size() << " \n";
    for (size_t i = 0; i < gaps.size(); ++i) {
        os << "3\n";
    }
}

template class FixedSphere<double>;

}
