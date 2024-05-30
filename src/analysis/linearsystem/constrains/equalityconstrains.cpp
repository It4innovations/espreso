
#include "equalityconstrains.h"
#include "analysis/assembler/structuralmechanics.h"
#include "analysis/math/vector_distributed.h"
#include "basis/containers/allocators.h"
#include "analysis/builder/feti.decomposition.h"

#include <algorithm>
#include <numeric>

namespace espreso {

template <typename T>
void EqualityConstrains<T>::set(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    feti.lambdas.intervals.push_back({ 0, 0 });
    feti.B1.resize(feti.K.size());
    feti.D2C.resize(feti.K.size());
    doffset.resize(feti.K.size());

    size_t maxMultiplicity = 2;
    std::vector<std::vector<int> > &D2C = feti.D2C;
    std::vector<std::vector<int> > COLS(feti.K.size()), FIXED(feti.K.size());
    std::vector<std::vector<T> > VALS(feti.K.size());

    struct __lambda__ { int dof, size; };
    std::vector<__lambda__> permutation;

    int dindex = 0, dof = 0;
    for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap, ++dof) {
        while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < dof) { ++dindex; }
        if (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] == dof) {
            for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                if (feti.decomposition->ismy(di->domain)) {
                    FIXED[di->domain - feti.decomposition->dbegin].push_back(di->index);
                }
            }
        } else {
            maxMultiplicity = std::max(maxMultiplicity, dmap->size());
            size_t i = 0;
            while (i < dmap->size() && !feti.decomposition->ismy(dmap->at(i).domain)) { ++i; }
            for (i = std::max(i, 1UL); i < dmap->size(); ++i) {
                permutation.push_back(__lambda__{ dof, (int)i + 1 });
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

    Matrix_Dense<double> lambdas;
    lambdas.resize(maxMultiplicity - 1, maxMultiplicity);
    math::set(lambdas.nrows * lambdas.ncols, lambdas.vals, 1, T{0});
    for (int r = 0, nc = 1; r < lambdas.nrows; ++r, ++nc) {
        double scale = std::sqrt(1 + (double)nc / (nc * nc));
        for (int c = 0; c < nc; ++c) {
            lambdas.vals[r * maxMultiplicity + c] = scale / (nc + 1);
        }
        lambdas.vals[r * maxMultiplicity + nc] = -scale * nc / (nc + 1);
    }

    feti.lambdas.size = 0;
    for (size_t i = 0, prev = feti.lambdas.cmap.size(); i < permutation.size(); ++i, ++feti.lambdas.size) {
        auto dmap = feti.decomposition->dmap->cbegin() + permutation[i].dof;
        size_t cbegin = feti.lambdas.cmap.size();
        feti.lambdas.cmap.push_back(1);
        feti.lambdas.cmap.push_back(permutation[i].size);
        if (dmap->at(0).domain < feti.decomposition->dbegin) {
            feti.lambdas.intervals.back().halo = feti.lambdas.size + 1;
        }
        for (int c = 0, r = permutation[i].size - 2; c < permutation[i].size; ++c) {
            feti.lambdas.cmap.push_back(dmap->at(c).domain);
            if (feti.decomposition->ismy(dmap->at(c).domain)) {
                D2C [dmap->at(c).domain - feti.decomposition->dbegin].push_back(feti.lambdas.size);
                COLS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(dmap->at(c).index);
                VALS[dmap->at(c).domain - feti.decomposition->dbegin].push_back(lambdas.vals[r * maxMultiplicity + c]);
            }
        }
        if (feti.lambdas.cmap[prev + 1] == feti.lambdas.cmap[cbegin + 1]) { // domains sizes
            if (std::equal(feti.lambdas.cmap.cbegin() + cbegin + 1, feti.lambdas.cmap.cend(), feti.lambdas.cmap.cbegin() + prev + 1)) {
                if (prev != cbegin) {
                    ++feti.lambdas.cmap[prev];
                    feti.lambdas.cmap.resize(cbegin);
                }
            } else {
                prev = cbegin;
            }
        } else {
            prev = cbegin;
        }
    }

    for (size_t d = 0; d < feti.K.size(); ++d) {
        doffset[d] = feti.lambdas.size;
        for (size_t i = 0; i < FIXED[d].size(); ++i) {
            D2C [d].push_back(feti.lambdas.size++);
            COLS[d].push_back(FIXED[d][i]);
            VALS[d].push_back(1);
        }
        if (FIXED[d].size()) {
            feti.lambdas.cmap.push_back(FIXED[d].size());
            feti.lambdas.cmap.push_back(1);
            feti.lambdas.cmap.push_back(feti.decomposition->dbegin + d);
        }
    }
    feti.lambdas.equalities = feti.lambdas.size;
    feti.lambdas.intervals.back().size = feti.lambdas.size - feti.lambdas.intervals.back().halo;

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        feti.B1[d].resize(COLS[d].size(), feti.K[d].nrows, COLS[d].size());
        std::iota(feti.B1[d].rows, feti.B1[d].rows + COLS[d].size() + 1, 0); // B1 is indexed from 0
        std::copy(COLS[d].begin(), COLS[d].end(), feti.B1[d].cols);
        std::copy(VALS[d].begin(), VALS[d].end(), feti.B1[d].vals);
    }

    feti.c.resize(feti.lambdas.size);
}

template <typename T>
void EqualityConstrains<T>::update(const step::Step &step, FETI<T> &feti, const Vector_Distributed<Vector_Sparse, T> &dirichlet)
{
    math::set(feti.c, T{0});
    std::vector<size_t> dindex = doffset;

    if (dirichlet.cluster.nnz) {
        auto map = feti.decomposition->dmap->cbegin();
        for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
            map += dirichlet.cluster.indices[i] - prev;
            for (auto di = map->begin(); di != map->end(); ++di) {
                if (feti.decomposition->ismy(di->domain)) {
                    feti.c.vals[dindex[di->domain - feti.decomposition->dbegin]++] = dirichlet.cluster.vals[i];
                }
            }
        }
    }
}

template struct EqualityConstrains<double>;

}


