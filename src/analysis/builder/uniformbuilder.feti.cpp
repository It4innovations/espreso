
#include "uniformbuilder.feti.h"

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
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

UniformBuilderFETIPattern::UniformBuilderFETIPattern(HeatTransferLoadStepConfiguration &configuration)
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

    dirichletInfo.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (configuration.temperature.find(region->name) != configuration.temperature.end()) {
            dirichletInfo[r].dirichlet = 1;
        }
    }

    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (configuration.temperature.find(region->name) != configuration.temperature.end()) {
            if (dirichletInfo[r].dirichlet) {
                for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
                    for (int d = 0; d < dofs; ++d) {
                        dirichlet.push_back(*n * dofs + d);
                    }
                }
            }
        }
    }

    dirichletInfo[0].size = dofs * (info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size);
    dirichletInfo[0].f = dirichlet; // use the first region to store indices permutation;
    utils::sortAndRemoveDuplicates(dirichlet);
    dirichletInfo[0].indices = dirichlet;
    for (size_t i = 0; i < dirichletInfo[0].f.size(); ++i) {
        dirichletInfo[0].f[i] = std::lower_bound(dirichlet.begin(), dirichlet.end(), dirichletInfo[0].f[i]) - dirichlet.begin();
    }

    fillDecomposition(configuration.feti, dofs);
    #pragma omp parallel for
    for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
        buildPattern(dofs, shape, domain);
    }
}
UniformBuilderFETIPattern::UniformBuilderFETIPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = info::mesh->dimension * multiplicity;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;

    dirichletInfo.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        auto expr = configuration.displacement.find(region->name);
        if (expr != configuration.displacement.end()) {
            for (int m = 0; m < multiplicity; ++m) {
                for (int d = 0; d < dofs; ++d) {
                    if (expr->second.data[d].isset) {
                        dirichletInfo[r].dirichlet += 1 << (dofs * m + d);
                    }
                }
            }
        }
    }

    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (configuration.displacement.find(region->name) != configuration.displacement.end()) {
            for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
                for (int m = 0; m < multiplicity; ++m) {
                    for (int d = 0; d < dofs; ++d) {
                        if (dirichletInfo[r].dirichlet & (1 << d)) {
                            dirichlet.push_back(*n * dofs * multiplicity + m * dofs + d);
                        }
                    }
                }
            }
        }
    }

    for (auto disp = configuration.fixed_wall.cbegin(); disp != configuration.fixed_wall.cend(); ++disp) {
        const BoundaryRegionStore *region = info::mesh->bregion(disp->first);
        for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
            for (int m = 0; m < multiplicity; ++m) {
                for (int d = 0; d < dofs; ++d) {
                    inequality.push_back(*n * dofs * multiplicity + m * dofs + d);
                }
            }
        }
    }
    utils::sortAndRemoveDuplicates(inequality);

    dirichletInfo[0].size = dofs * multiplicity * (info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size);
    dirichletInfo[0].f = dirichlet; // use the first region to store indices permutation;
    utils::sortAndRemoveDuplicates(dirichlet);
    dirichletInfo[0].indices = dirichlet;
    for (size_t i = 0; i < dirichletInfo[0].f.size(); ++i) {
        dirichletInfo[0].f[i] = std::lower_bound(dirichlet.begin(), dirichlet.end(), dirichletInfo[0].f[i]) - dirichlet.begin();
    }

    fillDecomposition(configuration.feti, dofs * multiplicity);
    #pragma omp parallel for
    for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
        buildPattern(dofs * multiplicity, shape, domain);
    }
}

void UniformBuilderFETIPattern::fillDecomposition(FETIConfiguration &feti, int dofs)
{
    elements.resize(info::mesh->domains->size);
    bregion.resize(info::mesh->domains->size);

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
    decomposition.dmap = new serializededata<esint, DIndex>(
            tarray<esint>(info::env::threads, distribution),
            tarray<DIndex>(info::mesh->nodes->domains->datatarray().distribution(), dofs));

    switch (feti.ordering) {
    case FETIConfiguration::ORDERING::ORDERED:
        // order: inner, lambdas, dirichlet, inequalities
        { // go through the rest dofs
            esint index = 0;
            auto fix = dirichlet.begin();
            auto ieq = inequality.begin();
            auto dmap = decomposition.dmap->begin();
            for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++index) {
                for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                    if (
                            domains->size() == 1 &&
                            (fix == dirichlet.end()  || *fix != dofs * index + dof) &&
                            (ieq == inequality.end() || *ieq != dofs * index + dof)
                            ) {

                        dmap->begin()->domain = *domains->begin();
                        dmap->begin()->index = elements[dmap->begin()->domain - decomposition.dbegin].size++;
                    }
                    if (fix != dirichlet.end()  && *fix == dofs * index + dof) ++fix;
                    if (ieq != inequality.end() && *ieq == dofs * index + dof) ++ieq;
                }
            }
        }
        { // go through shared nodes
            esint index = 0;
            auto fix = dirichlet.begin();
            auto dmap = decomposition.dmap->begin();
            for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++index) {
                for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                    if (fix == dirichlet.end() || *fix != dofs * index + dof) {
                        if (domains->size() > 1) {
                            auto di = dmap->begin();
                            for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                                di->domain = *d;
                                if (decomposition.ismy(*d)) {
                                    di->index = elements[*d - decomposition.dbegin].size++;
                                }
                            }
                        }
                    }
                    if (fix != dirichlet.end() && *fix == dofs * index + dof) ++fix;
                }
            }
        }
        // go through dirichlet
        for (auto i = dirichlet.begin(); i != dirichlet.end(); ++i) {
            auto domains = info::mesh->nodes->domains->begin() + *i / dofs;
            auto dmap = decomposition.dmap->begin() + *i;
            auto di = dmap->begin();
            for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                di->domain = *d;
                if (decomposition.ismy(*d)) {
                    di->index = elements[*d - decomposition.dbegin].size++;
                }
            }
        }
        // go through inequalities
        for (auto i = inequality.begin(); i != inequality.end(); ++i) {
            auto domains = info::mesh->nodes->domains->begin() + *i / dofs;
            if (domains->size() == 1) { // already indexed when domains->size > 1
                auto dmap = decomposition.dmap->begin() + *i;
                auto di = dmap->begin();
                for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                    di->domain = *d;
                    if (decomposition.ismy(*d)) {
                        di->index = elements[*d - decomposition.dbegin].size++;
                    }
                }
            }
        }
        break;
    case FETIConfiguration::ORDERING::NATURAL:
        auto dmap = decomposition.dmap->begin();
        for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains) {
            for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                auto di = dmap->begin();
                for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
                    di->domain = *d;
                    if (decomposition.ismy(*d)) {
                        di->index = elements[*d - decomposition.dbegin].size++;
                    }
                }
            }
        }
        break;
    }
}

void UniformBuilderFETIPattern::buildPattern(int dofs, Matrix_Shape shape, int domain)
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

    std::vector<esint> begin(elements[domain].size + 1);
    for (auto enodes = ebegin; enodes != eend; ++enodes) {
        for (auto n = enodes->begin(); n != enodes->end(); ++n) {
            for (int rd = 0; rd < dofs; ++rd) {
                begin[imap[*n * dofs + rd]] += dofs * enodes->size();
            }
        }
    }
    utils::sizesToOffsets(begin);

    size_t dataSize = 0;
    std::vector<esint> end = begin;
    std::vector<esint, initless_allocator<esint> > indices(begin.back());
    for (auto enodes = ebegin; enodes != eend; ++enodes) {
        for (auto from = enodes->begin(); from != enodes->end(); ++from) {
            for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                for (int rd = 0; rd < dofs; ++rd) {
                    for (int cd = 0; cd < dofs; ++cd) {
                        if (shape == Matrix_Shape::FULL || imap[*from * dofs + rd] <= imap[*to * dofs + cd]) {
                            ++dataSize; indices[end[imap[*from * dofs + rd]]++] = imap[*to * dofs + cd];
                        }
                    }
                }
            }
        }
    }

    size_t count = 0;
    for (esint n = 0; n < elements[domain].size; ++n) {
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

    elements[domain].row.reserve(count);
    elements[domain].column.reserve(count);
    elements[domain].K.reserve(dataSize);
    elements[domain].f.reserve(dataSize);

    std::vector<esint, initless_allocator<esint> > offset;
    offset.reserve(elements[domain].size);

    for (esint n = 0, size = 0; n < elements[domain].size; ++n) {
        offset.push_back(size);
        for (esint i = begin[n]; i < end[n]; ++i) {
            elements[domain].row.push_back(n);
            elements[domain].column.push_back(indices[i]);
        }
        size += end[n] - begin[n];
    }

    for (auto enodes = ebegin; enodes != eend; ++enodes) {
        for (int rd = 0, r = 0; rd < dofs; ++rd) {
            for (auto from = enodes->begin(); from != enodes->end(); ++from, ++r) {
                if (shape == Matrix_Shape::FULL) {
                    esint ri = imap[*from * dofs + rd];
                    auto ibegin = indices.begin() + begin[ri];
                    for (int cd = 0; cd < dofs; ++cd) {
                        for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                            esint coffset = 0;
                            while (ibegin[coffset] < imap[*to * dofs + cd]) ++coffset;
                            elements[domain].K.push_back(offset[ri] + coffset);
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
                                elements[domain].K.push_back(offset[min] + coffset);
                            }
                        }
                    }
                }
                elements[domain].f.push_back(imap[*from * dofs + rd]);
            }
        }
    }

    bregion[domain].resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
                auto element = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
                for (esint e = info::mesh->boundaryRegions[r]->eintervals[i].begin; e < info::mesh->boundaryRegions[r]->eintervals[i].end; ++e, ++element) {
                    for (int rd = 0; rd < dofs; ++rd) {
                        for (auto from = element->begin(); from != element->end(); ++from) {
                            bregion[domain][r].f.push_back(imap[*from * dofs + rd]);
                            if (shape == Matrix_Shape::FULL) {
                                esint ri = imap[*from * dofs + rd];
                                auto ibegin = indices.begin() + begin[ri];
                                for (int cd = 0; cd < dofs; ++cd) {
                                    for (auto to = element->begin(); to != element->end(); ++to) {
                                        esint coffset = 0;
                                        while (ibegin[coffset] < imap[*to * dofs + cd]) ++coffset;
                                        bregion[domain][r].K.push_back(offset[ri] + coffset);
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
                                        bregion[domain][r].K.push_back(offset[min] + coffset);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            // TODO: how to set node regions?
//            pattern->bregion[r].f.reserve(info::mesh->boundaryRegions[r]->nodes->datatarray().size());
//            for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); ++n) {
//                belement.clear();
//                for (int dof = 0; dof < dofs; ++dof) {
//                    belement.push_back(info::mesh->nodes->uniqInfo.position[*n] * dofs + dof);
//                }
//                for (size_t i = 0; i < belement.size(); ++i) {
//                    pattern->bregion[r].f.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), belement[i]) - RHSPattern.begin());
//                }
//            }
        }
    }
}

