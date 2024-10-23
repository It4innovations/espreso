
#include "matrix.uniform.feti.h"

#include "basis/containers/serializededata.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactinterfacestore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainsurfacestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

MatrixUniformFETI::MatrixUniformFETI(HeatTransferLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = 1;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    for (auto mat = info::ecf->heat_transfer.material_set.begin(); mat != info::ecf->heat_transfer.material_set.end(); ++mat) {
        if (info::ecf->heat_transfer.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
            shape = Matrix_Shape::FULL;
            type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
        }
    }
    if (configuration.translation_motions.size()) {
        shape = Matrix_Shape::FULL;
        type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
    }

    for (auto disc = info::ecf->heat_transfer.discretization.cbegin(); disc != info::ecf->heat_transfer.discretization.cend(); ++disc) {
        if (disc->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
            for (auto n = info::mesh->domainsSurface->nodes->datatarray().cbegin(); n != info::mesh->domainsSurface->nodes->datatarray().cend(); ++n) {
                surface.push_back(*n);
            }
            break;
        }
    }
    BEM = surface.size();

    fillDecomposition(configuration.feti, dofs);
    #pragma omp parallel for
    for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
        buildPattern(dofs, shape, domain);
    }
}

MatrixUniformFETI::MatrixUniformFETI(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = info::mesh->dimension * multiplicity;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;

    if (configuration.mode == StructuralMechanicsLoadStepConfiguration::MODE::NONLINEAR) {
        type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
    }

    if (configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
        type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    }

    for (auto disc = info::ecf->structural_mechanics.discretization.cbegin(); disc != info::ecf->structural_mechanics.discretization.cend(); ++disc) {
        if (disc->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
            for (auto n = info::mesh->domainsSurface->nodes->datatarray().cbegin(); n != info::mesh->domainsSurface->nodes->datatarray().cend(); ++n) {
                for (int m = 0; m < multiplicity; ++m) {
                    for (int d = 0; d < dofs; ++d) {
                        surface.push_back(*n * dofs * multiplicity + m * dofs + d);
                    }
                }
            }
            break;
        }
    }
    BEM = surface.size();

    fillDecomposition(configuration.feti, dofs * multiplicity);
    #pragma omp parallel for
    for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
        buildPattern(dofs * multiplicity, shape, domain);
    }
}

void MatrixUniformFETI::fillDecomposition(FETIConfiguration &feti, int dofs)
{
    elements.resize(info::mesh->domains->size);
    boundary.resize(info::mesh->domains->size);

    decomposition.dbegin = info::mesh->domains->offset;
    decomposition.dend = info::mesh->domains->offset + info::mesh->domains->size;
    decomposition.dtotal = info::mesh->domains->totalSize;
    decomposition.neighDomain.resize(info::mesh->neighbors.size() + 1);
    auto ddist = info::mesh->domains->gatherProcDistribution(); // remove this
    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        decomposition.neighDomain[n] = ddist[info::mesh->neighbors[n]];
    }
    decomposition.neighDomain.back() = decomposition.dbegin;

    decomposition.begin = dofs * info::mesh->nodes->uniqInfo.offset;
    decomposition.end = dofs * (info::mesh->nodes->uniqInfo.offset + info::mesh->nodes->uniqInfo.size);
    decomposition.totalSize = dofs * info::mesh->nodes->uniqInfo.totalSize;

    decomposition.neighbors = info::mesh->neighbors;
    decomposition.neighDOF.resize(decomposition.neighbors.size() + 1, decomposition.begin); // the last is my offset
    decomposition.halo.clear();
    decomposition.halo.reserve(dofs * info::mesh->nodes->uniqInfo.nhalo);
    for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n) {
        for (int dof = 0; dof < dofs; ++dof) {
            decomposition.halo.push_back(dofs * info::mesh->nodes->uniqInfo.position[n] + dof);
        }
    }

    std::vector<esint> dBuffer = { decomposition.begin };
    if (!Communication::gatherUniformNeighbors(dBuffer, decomposition.neighDOF, decomposition.neighbors)) {
        eslog::internalFailure("cannot exchange matrix distribution info.\n");
    }

    std::vector<esint> distribution(dofs * info::mesh->nodes->size + 1);
    auto domainmap = info::mesh->nodes->domains->begin();
    esint sum = 0;
    for (esint i = 0; i < info::mesh->nodes->size; ++i, ++domainmap) {
        for (int dof = 0; dof < dofs; ++dof, sum += domainmap->size()) {
            distribution[i * dofs + dof] = sum;
        }
    }
    distribution.back() = sum;

    // fill DMAP
    decomposition.dmap = new serializededata<esint, DecompositionFETI::DIndex>(
            tarray<esint>(info::env::threads, distribution),
            tarray<DecompositionFETI::DIndex>(info::mesh->nodes->domains->datatarray().distribution(), dofs, DecompositionFETI::DIndex{0, -1}));

    pattern.resize(info::mesh->domains->size);

    if (BEM) { // with BEM, TODO: mixed FEM & BEM
        // go through surface
        for (auto i = surface.begin(); i != surface.end(); ++i) {
            auto domains = info::mesh->nodes->domains->begin() + *i / dofs;
            auto dmap = decomposition.dmap->begin() + *i;
            auto di = dmap->begin();
            for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                di->domain = *d;
                if (decomposition.ismy(*d)) {
                    di->index = pattern[*d - decomposition.dbegin].size++;
                }
            }
        }

        for (size_t i = 0; i < pattern.size(); i++) {
            pattern[i].surface = pattern[i].size;
        }

        { // go through the rest dofs
            esint index = 0;
            auto surf = surface.begin();
            auto dmap = decomposition.dmap->begin();
            for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++index) {
                for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                    if (surf == surface.end() || *surf != dofs * index + dof) {
                        dmap->begin()->domain = *domains->begin();
                        dmap->begin()->index = pattern[dmap->begin()->domain - decomposition.dbegin].size++;
                    }
                    if (surf != surface.end() && *surf == dofs * index + dof) ++surf;
                }
            }
        }
    } else {
        auto dmap = decomposition.dmap->begin();
        for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains) {
            for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                auto di = dmap->begin();
                for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                    di->domain = *d;
                    if (decomposition.ismy(*d)) {
                        di->index = pattern[*d - decomposition.dbegin].size++;
                    }
                }
            }
        }
    }

    decomposition.dsize.resize(pattern.size());
    for (size_t i = 0; i < pattern.size(); ++i) {
        decomposition.dsize[i] = pattern[i].size;
    }
}

void MatrixUniformFETI::buildPattern(int dofs, Matrix_Shape shape, int domain)
{
    auto ebegin = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain];
    auto eend = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain + 1];

    size_t ni = 0;
    std::vector<esint> imap(dofs * info::mesh->nodes->size, -1);
    for (auto dmap = decomposition.dmap->cbegin(); dmap != decomposition.dmap->cend(); ++dmap, ++ni) {
        for (auto di = dmap->begin(); di != dmap->end(); ++di) {
            if (di->domain == domain + decomposition.dbegin) {
                imap[ni] = di->index;
            }
        }
    }

    std::vector<esint> begin(pattern[domain].size + 1);
    for (auto enodes = ebegin; enodes != eend; ++enodes) {
        for (auto n = enodes->begin(); n != enodes->end(); ++n) {
            for (int rd = 0; rd < dofs; ++rd) {
                begin[imap[*n * dofs + rd]] += dofs * enodes->size();
            }
        }
    }
    utils::sizesToOffsets(begin);

    std::vector<esint> end = begin;
    std::vector<esint, initless_allocator<esint> > indices(begin.back());
    for (auto enodes = ebegin; enodes != eend; ++enodes) {
        for (auto from = enodes->begin(); from != enodes->end(); ++from) {
            for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                for (int rd = 0; rd < dofs; ++rd) {
                    for (int cd = 0; cd < dofs; ++cd) {
                        if (shape == Matrix_Shape::FULL || imap[*from * dofs + rd] <= imap[*to * dofs + cd]) {
                            indices[end[imap[*from * dofs + rd]]++] = imap[*to * dofs + cd];
                        }
                    }
                }
            }
        }
    }

    size_t count = 0;
    for (esint n = 0; n < pattern[domain].size; ++n) {
        std::sort(indices.begin() + begin[n], indices.begin() + end[n]);
        esint unique = begin[n];
        for (auto i = begin[n] + 1; i < end[n]; ++i) {
            if (indices[unique] != indices[i]) {
                indices[++unique] = indices[i];
            }
        }
        end[n] = unique + 1;
        count += end[n] - begin[n];
    }

    pattern[domain].row.reserve(count);
    pattern[domain].column.reserve(count);

    std::vector<esint, initless_allocator<esint> > offset;
    offset.reserve(pattern[domain].size);

    for (esint n = 0, size = 0; n < pattern[domain].size; ++n) {
        offset.push_back(size);
        for (esint i = begin[n]; i < end[n]; ++i) {
            pattern[domain].row.push_back(n);
            pattern[domain].column.push_back(indices[i]);
        }
        size += end[n] - begin[n];
    }

    for (esint i = info::mesh->elements->eintervalsDistribution[domain]; i < info::mesh->elements->eintervalsDistribution[domain + 1]; ++i) {
        elements[domain].offset.push_back(elements[domain].permutation.size());
        const auto &einterval = info::mesh->elements->eintervals[i];
        for (auto enodes = info::mesh->elements->nodes->cbegin() + einterval.begin; enodes != info::mesh->elements->nodes->cbegin() + einterval.end; ++enodes) {
            for (int rd = 0, r = 0; rd < dofs; ++rd) {
                for (auto from = enodes->begin(); from != enodes->end(); ++from, ++r) {
                    if (shape == Matrix_Shape::FULL) {
                        esint ri = imap[*from * dofs + rd];
                        auto ibegin = indices.begin() + begin[ri];
                        for (int cd = 0; cd < dofs; ++cd) {
                            for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                                esint coffset = 0;
                                while (ibegin[coffset] < imap[*to * dofs + cd]) ++coffset;
                                elements[domain].permutation.push_back(offset[ri] + coffset);
                            }
                        }
                    } else {
                        for (int cd = 0, c = 0; cd < dofs; ++cd) {
                            for (auto to = enodes->begin(); to != enodes->end(); ++to, ++c) {
                                if (r <= c) {
                                    esint min = std::min(imap[*from * dofs + rd], imap[*to * dofs + cd]);
                                    esint max = std::max(imap[*from * dofs + rd], imap[*to * dofs + cd]);
                                    auto ibegin = indices.begin() + begin[min];
                                    esint coffset = 0;
                                    while (ibegin[coffset] < max) ++coffset;
                                    elements[domain].permutation.push_back(offset[min] + coffset);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    boundary[domain].resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (region->dimension) {
            for (esint i = region->eintervalsDistribution[domain]; i < region->eintervalsDistribution[domain + 1]; ++i) {
                boundary[domain][r].offset.push_back(boundary[domain][r].permutation.size());
                auto element = region->elements->cbegin() + region->eintervals[i].begin;
                for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e, ++element) {
                    for (int rd = 0; rd < dofs; ++rd) {
                        for (auto from = element->begin(); from != element->end(); ++from) {
                            if (shape == Matrix_Shape::FULL) {
                                esint ri = imap[*from * dofs + rd];
                                auto ibegin = indices.begin() + begin[ri];
                                for (int cd = 0; cd < dofs; ++cd) {
                                    for (auto to = element->begin(); to != element->end(); ++to) {
                                        esint coffset = 0;
                                        while (ibegin[coffset] < imap[*to * dofs + cd]) ++coffset;
                                        boundary[domain][r].permutation.push_back(offset[ri] + coffset);
                                    }
                                }
                            } else {
                                for (int cd = 0; cd < dofs; ++cd) {
                                    for (auto to = from; to != element->end(); ++to) {
                                        esint min = std::min(imap[*from * dofs + rd], imap[*to * dofs + cd]);
                                        esint max = std::max(imap[*from * dofs + rd], imap[*to * dofs + cd]);
                                        auto ibegin = indices.begin() + begin[min];
                                        esint coffset = 0;
                                        while (ibegin[coffset] < max) ++coffset;
                                        boundary[domain][r].permutation.push_back(offset[min] + coffset);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
void MatrixUniformFETI::Apply<T>::init(MatrixUniformFETI &m)
{
    in.decomposition = &m.decomposition;
    out.feti.decomposition = &m.decomposition;
    out.direct.resize(m.decomposition.halo.size() + m.decomposition.end - m.decomposition.begin);

    spblas.resize(m.decomposition.dsize.size());
    in.domains.resize(m.decomposition.dsize.size());
    out.feti.domains.resize(m.decomposition.dsize.size());
    for (size_t d = 0; d < m.decomposition.dsize.size(); ++d) {
        in.domains[d].resize(m.decomposition.dsize[d]);
        out.feti.domains[d].resize(m.decomposition.dsize[d]);
    }
}

template <typename T>
void MatrixUniformFETI::Apply<T>::apply(Matrix_FETI<T> &m, Vector_Distributed<Vector_Dense, T> &y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> &x)
{
    x.copyTo(&in);
    math::copy(out.direct, y.cluster);

    #pragma omp parallel for
    for (size_t i = 0; i < spblas.size(); ++i) {
        spblas[i].insert(m.domains[i]);
        spblas[i].apply(out.feti.domains[i], alpha, 0., in.domains[i]);
    }
    out.feti.sumTo(&y);
    math::add(y.cluster, beta, out.direct);
}

template struct MatrixUniformFETI::Sync<double>;
template struct MatrixUniformFETI::Apply<double>;


