
#include "cluster.faces.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"

namespace espreso {


template <typename T>
void ClusterFacesGluing<T>::set(const step::Step &step, FETI<T> &feti)
{
    if (feti.K.size() == 1) {
        return; // TOTAL FETI
    }
    feti.B0.resize(feti.K.size());
    feti.D2C0.resize(feti.K.size());

    struct __ijv__ {
        int i,j; double v;

        bool operator<(const __ijv__ &other) const {
            if (i == other.i) {
                return j < other.j;
            }
            return i < other.i;
        }
    };

    std::vector<esint> rindex;
    std::vector<std::vector<__ijv__> > B0(feti.B0.size());

    auto dual = info::mesh->domains->localDual->begin();
    std::vector<esint> rows(info::mesh->clusters->size);
    for (esint d1 = 0; d1 < info::mesh->domains->size; ++d1, ++dual) {
        esint cluster = info::mesh->domains->cluster[d1];
        for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
            if (d1 < *dit) {
                rindex.push_back(rows[cluster]);
                if (feti.R1[d1].nrows < feti.R1[*dit].nrows) {
                    rows[cluster] += feti.R1[*dit].nrows;
                } else {
                    rows[cluster] += feti.R1[d1].nrows ? feti.R1[d1].nrows : 1;
                }
            } else {
                auto dualbegin = info::mesh->domains->localDual->begin();
                auto it = std::lower_bound((dualbegin + *dit)->begin(), (dualbegin + *dit)->end(), d1);
                rindex.push_back(rindex[it - dualbegin->begin()]);
            }
        }
    }
    feti.cluster.gl_size = rows.front(); // TODO: more clusters

    dual = info::mesh->domains->localDual->begin();
    for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap) {
        for (auto di1 = dmap->begin(); di1 != dmap->end(); ++di1) {
            for (auto di2 = di1 + 1; di2 != dmap->end(); ++di2) {
                if (feti.decomposition->ismy(di1->domain) && feti.decomposition->ismy(di2->domain)) {
                    auto it = std::lower_bound(
                            (dual + (di1->domain - feti.decomposition->dbegin))->begin(),
                            (dual + (di1->domain - feti.decomposition->dbegin))->end(),
                            di2->domain - feti.decomposition->dbegin);

                    if (it != (dual + (di1->domain - feti.decomposition->dbegin))->end() && *it == di2->domain - feti.decomposition->dbegin) {
                        esint d1, d2, d1index, d2index;
                        if (di1->domain < di2->domain) {
                            d1 = di1->domain - feti.decomposition->dbegin;
                            d2 = di2->domain - feti.decomposition->dbegin;
                            d1index = di1->index;
                            d2index = di2->index;
                        } else {
                            d1 = di2->domain - feti.decomposition->dbegin;
                            d2 = di1->domain - feti.decomposition->dbegin;
                            d1index = di2->index;
                            d2index = di1->index;
                        }
                        if (feti.R1[d1].nrows) {
                            for (esint r = 0; r < feti.R1[d1].nrows; ++r) {
                                B0[d1].push_back({rindex[it - dual->begin()] + r, d1index,  feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                                B0[d2].push_back({rindex[it - dual->begin()] + r, d2index, -feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                            }
                        } else {
                            B0[d1].push_back({rindex[it - dual->begin()], d1index,  (double)(d1 + 1) / info::mesh->domains->size});
                            B0[d1].push_back({rindex[it - dual->begin()], d2index, -(double)(d1 + 1) / info::mesh->domains->size});
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (esint d = 0; d < info::mesh->domains->size; ++d) {
        std::sort(B0[d].begin(), B0[d].end());
        if (B0[d].size()) {
            feti.D2C0[d].push_back(B0[d][0].i);
        }
        for (size_t i = 1; i < B0[d].size(); ++i) {
            if (feti.D2C0[d].back() != B0[d][i].i) {
                feti.D2C0[d].push_back(B0[d][i].i);
            }
        }
        feti.B0[d].resize(feti.D2C0[d].size(), feti.K[d].ncols, B0[d].size());
        feti.B0[d].rows[0] = 0;
        feti.B0[d].cols[0] = B0[d][0].j;
        feti.B0[d].vals[0] = B0[d][0].v;
        for (size_t i = 1, r = 1; i < B0[d].size(); ++i) {
            feti.B0[d].cols[i] = B0[d][i].j;
            feti.B0[d].vals[i] = B0[d][i].v;
            if (B0[d][i - 1].i != B0[d][i].i) {
                feti.B0[d].rows[r++] = i;
            }
        }
        feti.B0[d].rows[feti.D2C0[d].size()] = B0[d].size();
    }
}

template <typename T>
void ClusterFacesGluing<T>::update(const step::Step &step, FETI<T> &feti)
{
    feti.B0.clear();
    feti.D2C0.clear();
    feti.B0.resize(feti.K.size());
    feti.D2C0.resize(feti.K.size());

    struct __ijv__ {
        int i,j; double v;

        bool operator<(const __ijv__ &other) const {
            if (i == other.i) {
                return j < other.j;
            }
            return i < other.i;
        }
    };

    std::vector<esint> rindex;
    std::vector<std::vector<__ijv__> > B0(feti.B0.size());

    auto dual = info::mesh->domains->localDual->begin();
    std::vector<esint> rows(info::mesh->clusters->size);
    for (esint d1 = 0; d1 < info::mesh->domains->size; ++d1, ++dual) {
        esint cluster = info::mesh->domains->cluster[d1];
        for (auto dit = dual->begin(); dit != dual->end(); ++dit) {
            if (d1 < *dit) {
                rindex.push_back(rows[cluster]);
                if (feti.R1[d1].nrows < feti.R1[*dit].nrows) {
                    rows[cluster] += feti.R1[*dit].nrows;
                } else {
                    rows[cluster] += feti.R1[d1].nrows ? feti.R1[d1].nrows : 1;
                }
            } else {
                auto dualbegin = info::mesh->domains->localDual->begin();
                auto it = std::lower_bound((dualbegin + *dit)->begin(), (dualbegin + *dit)->end(), d1);
                rindex.push_back(rindex[it - dualbegin->begin()]);
            }
        }
    }
    feti.cluster.gl_size = rows.front(); // TODO: more clusters

    dual = info::mesh->domains->localDual->begin();
    for (auto dmap = feti.decomposition->dmap->cbegin(); dmap != feti.decomposition->dmap->cend(); ++dmap) {
        for (auto di1 = dmap->begin(); di1 != dmap->end(); ++di1) {
            for (auto di2 = di1 + 1; di2 != dmap->end(); ++di2) {
                if (feti.decomposition->ismy(di1->domain) && feti.decomposition->ismy(di2->domain)) {
                    auto it = std::lower_bound(
                            (dual + (di1->domain - feti.decomposition->dbegin))->begin(),
                            (dual + (di1->domain - feti.decomposition->dbegin))->end(),
                            di2->domain - feti.decomposition->dbegin);

                    if (it != (dual + (di1->domain - feti.decomposition->dbegin))->end() && *it == di2->domain - feti.decomposition->dbegin) {
                        esint d1, d2, d1index, d2index;
                        if (di1->domain < di2->domain) {
                            d1 = di1->domain - feti.decomposition->dbegin;
                            d2 = di2->domain - feti.decomposition->dbegin;
                            d1index = di1->index;
                            d2index = di2->index;
                        } else {
                            d1 = di2->domain - feti.decomposition->dbegin;
                            d2 = di1->domain - feti.decomposition->dbegin;
                            d1index = di2->index;
                            d2index = di1->index;
                        }
                        if (feti.R1[d1].nrows) {
                            for (esint r = 0; r < feti.R1[d1].nrows; ++r) {
                                B0[d1].push_back({rindex[it - dual->begin()] + r, d1index,  feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                                B0[d2].push_back({rindex[it - dual->begin()] + r, d2index, -feti.R1[d1].vals[feti.R1[d1].ncols * r + d1index]});
                            }
                        } else {
                            B0[d1].push_back({rindex[it - dual->begin()], d1index,  (double)(d1 + 1) / info::mesh->domains->size});
                            B0[d1].push_back({rindex[it - dual->begin()], d2index, -(double)(d1 + 1) / info::mesh->domains->size});
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (esint d = 0; d < info::mesh->domains->size; ++d) {
        std::sort(B0[d].begin(), B0[d].end());
        if (B0[d].size()) {
            feti.D2C0[d].push_back(B0[d][0].i);
        }
        for (size_t i = 1; i < B0[d].size(); ++i) {
            if (feti.D2C0[d].back() != B0[d][i].i) {
                feti.D2C0[d].push_back(B0[d][i].i);
            }
        }
        feti.B0[d].resize(feti.D2C0[d].size(), feti.K[d].ncols, B0[d].size());
        feti.B0[d].rows[0] = 0;
        feti.B0[d].cols[0] = B0[d][0].j;
        feti.B0[d].vals[0] = B0[d][0].v;
        for (size_t i = 1, r = 1; i < B0[d].size(); ++i) {
            feti.B0[d].cols[i] = B0[d][i].j;
            feti.B0[d].vals[i] = B0[d][i].v;
            if (B0[d][i - 1].i != B0[d][i].i) {
                feti.B0[d].rows[r++] = i;
            }
        }
        feti.B0[d].rows[feti.D2C0[d].size()] = B0[d].size();
    }
}

template struct ClusterFacesGluing<double>;


}
