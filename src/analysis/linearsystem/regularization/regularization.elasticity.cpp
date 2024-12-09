
#include "regularization.h"

#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/domainsurfacestore.h"
#include "math/wrappers/math.blas.h"
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.lapack.h"
#include "wrappers/metis/w.metis.h"

#include "math/math.h"

#include <algorithm>
#include <numeric>
#include <random>

namespace espreso {

static void getCorners(std::vector<esint> &fixPoints, int domain)
{
    esint begin = info::mesh->domains->elements[domain];
    esint end = info::mesh->domains->elements[domain + 1];
    auto element = info::mesh->elements->nodes->begin() + begin;
    auto *epointer = &info::mesh->elements->epointers->datatarray();

    int ex = info::ecf->generator.grid.elements_x;
    int ey = info::ecf->generator.grid.elements_y;
    int ez = info::ecf->generator.grid.elements_z;

    std::vector<esint> nodes;
    for (esint e = 0; e < end - begin; ++e, ++element) {
        for (int n = 0; n < (*epointer)[begin + e]->coarseNodes; ++n) {
            nodes.push_back(element->at(n));
        }
    }
    utils::sortAndRemoveDuplicates(nodes);
    std::sort(nodes.begin(), nodes.end(), [&] (esint i, esint j) {
        const Point &pi = info::mesh->nodes->coordinates->datatarray()[i];
        const Point &pj = info::mesh->nodes->coordinates->datatarray()[j];
        if (pi.z == pj.z) {
            if (pi.y == pj.y) {
                return pi.x < pj.x;
            }
            return pi.y < pj.y;
        }
        return pi.z < pj.z;
    });

    fixPoints.push_back(nodes[0]);
    fixPoints.push_back(nodes[ex]);
    fixPoints.push_back(nodes[(ex + 1) * ey]);
    fixPoints.push_back(nodes[(ex + 1) * ey + ex]);

    if (info::mesh->dimension == 3) {
        fixPoints.push_back(nodes[ez * (ex + 1) * (ey + 1)]);
        fixPoints.push_back(nodes[ez * (ex + 1) * (ey + 1) + ex]);
        fixPoints.push_back(nodes[ez * (ex + 1) * (ey + 1) + (ex + 1) * ey]);
        fixPoints.push_back(nodes[ez * (ex + 1) * (ey + 1) + (ex + 1) * ey + ex]);
    }

    std::sort(fixPoints.begin(), fixPoints.end());
}

static void getFixPoints(std::vector<esint> &fixPoints, int domain, bool onSurface)
{
    if (
            info::ecf->input_type == ECF::INPUT_TYPE::GENERATOR &&
            info::ecf->generator.uniform_clusters && info::ecf->generator.uniform_domains &&
            info::ecf->generator.shape == INPUT_GENERATOR_SHAPE::GRID) {
        getCorners(fixPoints, domain);
        return;
    }

    esint begin = info::mesh->domains->elements[domain];
    esint end = info::mesh->domains->elements[domain + 1];
    auto element = info::mesh->elements->nodes->begin() + begin;
    auto *epointer = &info::mesh->elements->epointers->datatarray();
    if (onSurface) {
        begin = info::mesh->domainsSurface->edistribution[domain];
        end = info::mesh->domainsSurface->edistribution[domain + 1];
        element = info::mesh->domainsSurface->enodes->begin() + begin;
        epointer = &info::mesh->domainsSurface->epointers->datatarray();
    }

    size_t FIX_POINTS_SIZE = 8;

    auto neighs = [] (std::vector<esint> &neighs, Element::CODE code, int node, const esint* nodes) {
        switch (code) {
        case Element::CODE::HEXA8:
        case Element::CODE::HEXA20:
            if (node < 4) {
                neighs.push_back((nodes[(node + 1) % 4]));
                neighs.push_back((nodes[(node + 3) % 4]));
                neighs.push_back((nodes[node + 4]));
            } else {
                neighs.push_back((nodes[(node + 1) % 4 + 4]));
                neighs.push_back((nodes[(node + 3) % 4 + 4]));
                neighs.push_back((nodes[node - 4]));
            }
            return 3;
        case Element::CODE::TETRA4:
        case Element::CODE::TETRA10:
            neighs.push_back(nodes[(node + 1) % 4]);
            neighs.push_back(nodes[(node + 2) % 4]);
            neighs.push_back(nodes[(node + 3) % 4]);
            return 3;
        case Element::CODE::PRISMA6:
        case Element::CODE::PRISMA15:
            if (node < 3) {
                neighs.push_back(nodes[(node + 1) % 3]);
                neighs.push_back(nodes[(node + 2) % 3]);
                neighs.push_back(nodes[node + 3]);
            } else {
                neighs.push_back(nodes[(node + 1) % 3 + 3]);
                neighs.push_back(nodes[(node + 2) % 3 + 3]);
                neighs.push_back(nodes[node - 3]);
            }
            return 3;

        case Element::CODE::PYRAMID5:
        case Element::CODE::PYRAMID13:
            if (node == 4) {
                neighs.insert(neighs.end(), nodes, nodes + 4);
                return 4;
            } else {
                neighs.push_back(nodes[(node + 1) % 4]);
                neighs.push_back(nodes[(node + 3) % 4]);
                neighs.push_back(nodes[4]);
                return 3;
            }

        case Element::CODE::TRIANGLE3:
        case Element::CODE::TRIANGLE6:
            neighs.push_back(nodes[(node + 1) % 3]);
            neighs.push_back(nodes[(node + 2) % 3]);
            return 2;

        case Element::CODE::SQUARE4:
        case Element::CODE::SQUARE8:
            neighs.push_back(nodes[(node + 1) % 4]);
            neighs.push_back(nodes[(node + 3) % 4]);
            return 2;

        case Element::CODE::LINE2:
        case Element::CODE::LINE3:
            neighs.push_back(nodes[(node + 1) % 2]);
            return 1;
        case Element::CODE::POINT1:
        default:
            return 0;
        }
        return 0;
    };

    std::vector<esint> originnodes, neighsnodes;
    originnodes.reserve((end - begin) * 20);
    for (esint e = 0; e < end - begin; ++e, ++element) {
        for (int n = 0; n < (*epointer)[begin + e]->coarseNodes; ++n) {
            originnodes.insert(
                    originnodes.end(),
                    neighs(neighsnodes, (*epointer)[begin + e]->code, n, element->data()),
                    element->at(n));
        }
    }

    std::vector<esint> permutation(originnodes.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
        return originnodes[i] < originnodes[j];
    });

    std::vector<esint> ids, dist, data;
    dist.push_back(0);
    ids.push_back(originnodes[permutation[0]]);
    data.push_back(neighsnodes[permutation[0]]);
    for (size_t i = 1; i < permutation.size(); i++) {
        if (originnodes[permutation[i]] != originnodes[permutation[i - 1]]) {
            utils::sortAndRemoveDuplicates(data, dist.back());
            dist.push_back(data.size());
            ids.push_back(originnodes[permutation[i]]);
        }
        data.push_back(neighsnodes[permutation[i]]);
    }
    utils::sortAndRemoveDuplicates(data, dist.back());
    dist.push_back(data.size());

    if (ids.size() <= FIX_POINTS_SIZE * FIX_POINTS_SIZE / 2) {
        std::random_device rd;
        std::mt19937 g(rd());

        std::shuffle(ids.begin(), ids.end(), g);
        if (FIX_POINTS_SIZE < ids.size()) {
            ids.resize(FIX_POINTS_SIZE);
        }
        fixPoints = ids;
        utils::sortAndRemoveDuplicates(fixPoints);
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        data[i] = std::lower_bound(ids.begin(), ids.end(), data[i]) - ids.begin();
    }

    std::vector<esint> partition(ids.size());
    METISConfiguration options;
    METIS::call(options, ids.size(), dist.data(), data.data(), 0, NULL, NULL, FIX_POINTS_SIZE, partition.data());

    std::vector<std::vector<int> > pids(FIX_POINTS_SIZE), pdist(FIX_POINTS_SIZE, { Indexing::CSR }), pdata(FIX_POINTS_SIZE);
    for (size_t i = 0; i < partition.size(); i++) {
        pids[partition[i]].push_back(i);
    }
    for (size_t i = 0; i < partition.size(); i++) {
        esint p = partition[i];
        for (esint j = dist[i]; j < dist[i + 1]; j++) {
            if (partition[data[j]] == p) {
                size_t index = std::lower_bound(pids[p].begin(), pids[p].end(), data[j]) - pids[p].begin();
                if (pdist[p].size() <= index) {
                    pdata[p].push_back(index + Indexing::CSR);
                }
            }
        }
        pdist[p].push_back(pdata[p].size() + Indexing::CSR);
    }

    for (size_t p = 0; p < FIX_POINTS_SIZE; p++) {
        if (pids[p].size()) {
            std::vector<float> vals(pdata[p].size(), 1), x(pids[p].size(), 1. / pids[p].size()), y(pids[p].size());
            Matrix_CSR<float> M;
            M.shape = Matrix_Shape::UPPER;
            M.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
            Vector_Dense<float> in, out;
            in.size = out.size = M.nrows = M.ncols = pids[p].size();
            M.nnz = vals.size();
            M.rows = pdist[p].data();
            M.cols = pdata[p].data();
            M.vals = vals.data();
            in.vals = x.data();
            out.vals = y.data();
            SpBLAS<Matrix_CSR, float> spblas(M);

            float last_l = pids[p].size(), l = 1;
            while (fabs((l - last_l) / l) > 1e-6) {
                spblas.apply(out, 1, 0, in);
                last_l = l;
                l = math::blas::norm(out.size, out.vals, 1);
                math::blas::scale(out.size, 1 / l, out.vals, 1);
                out.swap(in);
            }
            fixPoints.push_back(ids[pids[p][std::max_element(in.vals, in.vals + in.size) - in.vals]]);
        }
    }
    utils::sortAndRemoveDuplicates(fixPoints);
}

template <typename T>
static void getNtNNtN(Matrix_Dense<T> &N, Matrix_Dense<T> &NtNNtN)
{
    NtNNtN.resize(N.ncols, N.ncols);
    Matrix_Dense<T> _N(N), NNt; NNt.resize(N.nrows, N.nrows);
    math::blas::AAt(N, NNt);
    math::lapack::solve_sym_upper(NNt, _N);
    math::blas::multiply(T{1}, N, _N, T{0}, NtNNtN, true);
}

template <typename T>
static void setRegMat(Matrix_CSR<T> &K, Matrix_CSR<T> &RegMat, Matrix_Dense<T> &NtNNtN, DecompositionFETI *decomposition, int domain, bool onSurface)
{
    std::vector<esint> fixPoints;
    getFixPoints(fixPoints, domain, onSurface);

    std::vector<int> fixCols, permutation;
    std::vector<double> fixVals;
    if (info::mesh->dimension == 2) {
        fixCols.resize(2 * fixPoints.size());
        fixVals.resize(2 * fixPoints.size() * 3);
        for (size_t i = 0; i < fixPoints.size(); ++i) {
            auto dmap0 = decomposition->dmap->begin() + (2 * fixPoints[i] + 0);
            auto dmap1 = decomposition->dmap->begin() + (2 * fixPoints[i] + 1);
            for (auto di = dmap0->begin(); di != dmap0->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    fixCols[ 2 * i + 0         ] = di->index;
                    fixVals[(2 * i + 0) * 3 + 0] = 1;
                    fixVals[(2 * i + 0) * 3 + 1] = 0;
                    fixVals[(2 * i + 0) * 3 + 2] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
                }
            }
            for (auto di = dmap1->begin(); di != dmap1->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    fixCols[ 2 * i + 1         ] = di->index;
                    fixVals[(2 * i + 1) * 3 + 0] = 0;
                    fixVals[(2 * i + 1) * 3 + 1] = 1;
                    fixVals[(2 * i + 1) * 3 + 2] = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
                }
            }
        }
    }
    if (info::mesh->dimension == 3) {
        fixCols.resize(3 * fixPoints.size());
        fixVals.resize(3 * fixPoints.size() * 6);
        for (size_t i = 0; i < fixPoints.size(); ++i) {
            auto dmap0 = decomposition->dmap->begin() + (3 * fixPoints[i] + 0);
            auto dmap1 = decomposition->dmap->begin() + (3 * fixPoints[i] + 1);
            auto dmap2 = decomposition->dmap->begin() + (3 * fixPoints[i] + 2);
            for (auto di = dmap0->begin(); di != dmap0->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    fixCols[ 3 * i + 0         ] = di->index;
                    fixVals[(3 * i + 0) * 6 + 0] = 1;
                    fixVals[(3 * i + 0) * 6 + 1] = 0;
                    fixVals[(3 * i + 0) * 6 + 2] = 0;
                    fixVals[(3 * i + 0) * 6 + 3] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
                    fixVals[(3 * i + 0) * 6 + 4] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].z;
                    fixVals[(3 * i + 0) * 6 + 5] = 0;
                }
            }
            for (auto di = dmap1->begin(); di != dmap1->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    fixCols[ 3 * i + 1         ] = di->index;
                    fixVals[(3 * i + 1) * 6 + 0] = 0;
                    fixVals[(3 * i + 1) * 6 + 1] = 1;
                    fixVals[(3 * i + 1) * 6 + 2] = 0;
                    fixVals[(3 * i + 1) * 6 + 3] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
                    fixVals[(3 * i + 1) * 6 + 4] = 0;
                    fixVals[(3 * i + 1) * 6 + 5] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].z;
                }
            }
            for (auto di = dmap2->begin(); di != dmap2->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    fixCols[ 3 * i + 2         ] = di->index;
                    fixVals[(3 * i + 2) * 6 + 0] = 0;
                    fixVals[(3 * i + 2) * 6 + 1] = 0;
                    fixVals[(3 * i + 2) * 6 + 2] = 1;
                    fixVals[(3 * i + 2) * 6 + 3] = 0;
                    fixVals[(3 * i + 2) * 6 + 4] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
                    fixVals[(3 * i + 2) * 6 + 5] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
                }
            }
        }
    }

    permutation.resize(fixCols.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(), [&fixCols] (int c1, int c2) { return fixCols[c1] < fixCols[c2]; });

    RegMat.resize(K.nrows, K.nrows, (fixCols.size() - 1) * fixCols.size() / 2 + fixCols.size());
    RegMat.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    RegMat.shape = Matrix_Shape::UPPER;

    RegMat.rows[0] = Indexing::CSR;
    esint r = 0;
    for (size_t i = 0; i < permutation.size(); ++i, ++r) {
        while (r < fixCols[permutation[i]]) {
            RegMat.rows[r + 1] = RegMat.rows[r];
            ++r;
        }
        RegMat.rows[r + 1] = RegMat.rows[r] + fixCols.size() - i;
        for (size_t j = i; j < fixCols.size(); ++j) {
            RegMat.cols[RegMat.rows[r] + j - i - Indexing::CSR] = fixCols[permutation[j]] + Indexing::CSR;
        }
    }
    while (r < RegMat.nrows) {
        RegMat.rows[r + 1] = RegMat.rows[r];
        ++r;
    }
    Matrix_Dense<T> N;
    if (info::mesh->dimension == 2) {
        N.resize(3, fixCols.size());
        for (size_t i = 0; i < permutation.size(); ++i) {
            N.vals[0 * N.ncols + i] = fixVals[3 * permutation[i] + 0];
            N.vals[1 * N.ncols + i] = fixVals[3 * permutation[i] + 1];
            N.vals[2 * N.ncols + i] = fixVals[3 * permutation[i] + 2];
        }
    }
    if (info::mesh->dimension == 3) {
        N.resize(6, fixCols.size());
        for (size_t i = 0; i < permutation.size(); ++i) {
            N.vals[0 * N.ncols + i] = fixVals[6 * permutation[i] + 0];
            N.vals[1 * N.ncols + i] = fixVals[6 * permutation[i] + 1];
            N.vals[2 * N.ncols + i] = fixVals[6 * permutation[i] + 2];
            N.vals[3 * N.ncols + i] = fixVals[6 * permutation[i] + 3];
            N.vals[4 * N.ncols + i] = fixVals[6 * permutation[i] + 4];
            N.vals[5 * N.ncols + i] = fixVals[6 * permutation[i] + 5];
        }
    }
    getNtNNtN(N, NtNNtN);
}

template <typename T>
static void updateRegMat(Matrix_CSR<T> &K, Matrix_CSR<T> &RegMat, Matrix_Dense<T> &NtNNtN)
{
    double max = 0;
    for (esint r = 0; r < K.nrows; ++r) {
        max = std::max(max, K.vals[K.rows[r] - Indexing::CSR]);
    }

    for (esint r = 0, i = 0; r < NtNNtN.nrows; ++r) {
        for (esint c = r; c < NtNNtN.ncols; ++c, ++i) {
            RegMat.vals[i] = max * NtNNtN.vals[r * NtNNtN.ncols + c];
        }
    }
}

template <typename T>
static void setR1(Matrix_CSR<T> &K, Matrix_Dense<T> &R1, DecompositionFETI *decomposition, int domain)
{
    if (info::mesh->dimension == 2) {
        R1.resize(3, K.nrows);
        int i = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++i) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin) {
                    switch (i % 2) {
                    case 0:
                        R1.vals[0 * R1.ncols + di->index] = 1;
                        R1.vals[1 * R1.ncols + di->index] = 0;
                        R1.vals[2 * R1.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 2].y;
                        break;
                    case 1:
                        R1.vals[0 * R1.ncols + di->index] = 0;
                        R1.vals[1 * R1.ncols + di->index] = 1;
                        R1.vals[2 * R1.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 2].x;
                        break;
                    }
                }
            }
        }
    }
    if (info::mesh->dimension == 3) {
        R1.resize(6, K.nrows);
        int i = 0;
        for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++i) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (di->domain == domain + decomposition->dbegin && di->index < K.nrows) {
                    switch (i % 3) {
                    case 0:
                        R1.vals[0 * R1.ncols + di->index] = 1;
                        R1.vals[1 * R1.ncols + di->index] = 0;
                        R1.vals[2 * R1.ncols + di->index] = 0;
                        R1.vals[3 * R1.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].y;
                        R1.vals[4 * R1.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].z;
                        R1.vals[5 * R1.ncols + di->index] = 0;
                        break;
                    case 1:
                        R1.vals[0 * R1.ncols + di->index] = 0;
                        R1.vals[1 * R1.ncols + di->index] = 1;
                        R1.vals[2 * R1.ncols + di->index] = 0;
                        R1.vals[3 * R1.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].x;
                        R1.vals[4 * R1.ncols + di->index] = 0;
                        R1.vals[5 * R1.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].z;
                        break;
                    case 2:
                        R1.vals[0 * R1.ncols + di->index] = 0;
                        R1.vals[1 * R1.ncols + di->index] = 0;
                        R1.vals[2 * R1.ncols + di->index] = 1;
                        R1.vals[3 * R1.ncols + di->index] = 0;
                        R1.vals[4 * R1.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].x;
                        R1.vals[5 * R1.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].y;
                        break;
                    }
                }
            }
        }
    }
}

template <typename T>
void Regularization<T>::set(FETI<T> &feti, StructuralMechanicsLoadStepConfiguration &configuration)
{
    NtNNtN.resize(feti.K.size());
    if (feti.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            if (R1)     setR1    (feti.K[d], feti.R1[d], feti.decomposition, d);
            if (regMat) setRegMat(feti.K[d], feti.RegMat[d], NtNNtN[d], feti.decomposition, d, onSurface);
        }
        if (R1) {
            orthonormalize(feti);
        }
    }
}

template <typename T>
void Regularization<T>::update(FETI<T> &feti, StructuralMechanicsLoadStepConfiguration &configuration)
{
    if (feti.configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            if (regMat && feti.updated.K) updateRegMat(feti.K[d], feti.RegMat[d], NtNNtN[d]);
        }
    } else {
        algebraic(feti, 3 * (info::mesh->dimension - 1), feti.configuration.sc_size);
    }
}

template void Regularization<double>::set    (FETI<double> &feti, StructuralMechanicsLoadStepConfiguration &configuration);
template void Regularization<double>::update (FETI<double> &feti, StructuralMechanicsLoadStepConfiguration &configuration);


}
