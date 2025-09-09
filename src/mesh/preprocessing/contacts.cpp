
#include "meshpreprocessing.h"

#include "analysis/assembler/general/basefunctions.h"
#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/bodystore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/contactstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/contactinterfacestore.h"

#include "basis/containers/point.h"
#include "basis/containers/serializededata.h"
#include "basis/structures/intervaltree.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "basis/logging/timelogger.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "output/visualization/debug.h"
#include "wrappers/mpi/communication.h"

#include "math/primitives/matrix_dense.h"
#include "math/wrappers/math.lapack.h"

#include <algorithm>
#include <map>
#include <numeric>
#include <unordered_set>
#include <unordered_map>

namespace espreso {
namespace mesh {

void computeBodiesSurface(NodeStore *nodes, ElementStore *elements, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<BoundaryRegionStore*> &boundaryRegions, SurfaceStore *surface, std::vector<int> &neighbors)
{
    profiler::syncstart("compute_bodies_surface");
    auto eregion = [elementsRegions] (const std::string &name) -> ElementsRegionStore* {
        for (size_t r = 0; r < elementsRegions.size(); r++) {
            if (StringCompare::caseSensitiveEq(elementsRegions[r]->name, name)) {
                return elementsRegions[r];
            }
        }
        eslog::error("Unknown region of elements with name '%s'.\n", name.c_str());
        return NULL;
    };

    for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
        switch (it->second.detection) {
        case ContactInterfaceConfiguration::DETECTION::ALL_BODIES:
            for (size_t r = 0; r < elementsRegions.size(); ++r) {
                elementsRegions[r]->contact.gap = std::max(it->second.gap, elementsRegions[r]->contact.gap);
                elementsRegions[r]->contact.angle = std::max(it->second.angle, elementsRegions[r]->contact.angle);
                elementsRegions[r]->contact.self_contact |= it->second.self_contact;
            }
            break;
        case ContactInterfaceConfiguration::DETECTION::BODY_LIST:
            for (size_t b = 0; b < it->second.body_list.size(); ++b) {
                eregion(it->second.body_list[b])->contact.gap = std::max(it->second.gap, eregion(it->second.body_list[b])->contact.gap);
                eregion(it->second.body_list[b])->contact.angle = std::max(it->second.angle, eregion(it->second.body_list[b])->contact.angle);
                eregion(it->second.body_list[b])->contact.self_contact |= it->second.self_contact;
            }
            break;
        case ContactInterfaceConfiguration::DETECTION::CONTACT_PAIR:
            eslog::internalFailure("implement CONTACT_PAIR detection.\n");
            break;
        }
    }

    elements->contact = new serializededata<esint, ContactInfo>(1, tarray<ContactInfo>(elements->epointers->datatarray().distribution(), 1, ContactInfo()));
    for (size_t r = 0; r < elementsRegions.size(); ++r) {
        for (auto e = elementsRegions[r]->elements->datatarray().begin(); e != elementsRegions[r]->elements->datatarray().end(); ++e) {
            elements->contact->datatarray()[*e].gap = std::max(elementsRegions[r]->contact.gap, elements->contact->datatarray()[*e].gap);
            elements->contact->datatarray()[*e].angle = std::max(elementsRegions[r]->contact.angle, elements->contact->datatarray()[*e].angle);
            elements->contact->datatarray()[*e].self_contact |= elementsRegions[r]->contact.self_contact;
        }
    }

    size_t threads = info::env::OMP_NUM_THREADS;

    std::vector<std::vector<esint> > faces(threads), facesDistribution(threads), parents(threads), body(threads), ecounters(threads, std::vector<esint>((int)Element::CODE::SIZE));
    std::vector<std::vector<Element*> > fpointers(threads);

    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        auto nodes = elements->nodes->cbegin(t);
        auto neighs = elements->faceNeighbors->cbegin(t);
        const auto &epointers = elements->epointers->datatarray();

        std::vector<esint> fdist, fdata, fparents, fbody, ecounter((int)Element::CODE::SIZE);
        std::vector<Element*> fpointer;
        if (t == 0) {
            fdist.push_back(0);
        }

        for (size_t e = elements->distribution.threads[t]; e < elements->distribution.threads[t + 1]; ++e, ++neighs, ++nodes) {
            for (size_t n = 0; n < neighs->size(); ++n) {
                if (neighs->at(n) == -1) {
                    auto face = epointers[e]->faceList->begin() + n;
                    for (auto f = face->begin(); f != face->end(); ++f) {
                        fdata.push_back(nodes->at(*f));
                    }
                    fdist.push_back(fdata.size());
                    fpointer.push_back(epointers[e]->facepointers->datatarray()[n]);
                    fparents.push_back(e);
                    fbody.push_back(elements->body->datatarray()[e]);
                    ++ecounter[(int)fpointer.back()->code];
                }
            }
        }

        facesDistribution[t].swap(fdist);

        faces[t].swap(fdata);
        fpointers[t].swap(fpointer);
        parents[t].swap(fparents);
        body[t].swap(fbody);
        ecounters[t].swap(ecounter);
    }

    for (size_t t = 1; t < threads; t++) {
        for (size_t e = 0; e < ecounters[0].size(); e++) {
            ecounters[0][e] += ecounters[t][e];
        }
    }

    serializededata<esint, Element*>::balance(1, fpointers);
    surface->epointers = new serializededata<esint, Element*>(1, fpointers);
    surface->ecounters = ecounters[0];

    surface->edistribution = surface->epointers->datatarray().distribution();

    if (surface->ecounters[(int)Element::CODE::TRIANGLE3] == (esint)surface->edistribution.back()) {
        serializededata<esint, esint>::balance(3, faces, &surface->edistribution);
        surface->enodes = new serializededata<esint, esint>(3, faces);
        surface->triangles = surface->enodes;
        surface->tdistribution = surface->edistribution;
    } else {
        utils::threadDistributionToFullDistribution(facesDistribution);
        serializededata<esint, esint>::balance(facesDistribution, faces, &surface->edistribution);
        surface->enodes = new serializededata<esint, esint>(facesDistribution, faces);
    }
    serializededata<esint, esint>::balance(1, parents, &surface->edistribution);
    surface->parents = new serializededata<esint, esint>(1, parents);

    serializededata<esint, esint>::balance(1, body, &surface->edistribution);
    surface->body = new serializededata<esint, esint>(1, body);

    surface->nIDs = new serializededata<esint, esint>(*nodes->IDs);
    surface->coordinates = new serializededata<esint, Point>(*nodes->coordinates);

    std::vector<esint> snodes(surface->enodes->datatarray().begin(), surface->enodes->datatarray().end());
    utils::sortAndRemoveDuplicates(snodes);

    boundaryRegions.push_back(new BoundaryRegionStore("SURFACE"));
    boundaryRegions.back()->originalDimension = boundaryRegions.back()->dimension = info::mesh->dimension - 1;
    boundaryRegions.back()->distribution.threads = surface->edistribution;
    boundaryRegions.back()->epointers = new serializededata<esint, Element*>(*surface->epointers);
    boundaryRegions.back()->elements = new serializededata<esint, esint>(*surface->enodes);

    std::vector<std::vector<ContactInfo> > contact(threads);

    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        std::vector<ContactInfo> tcontact;
        tcontact.reserve(surface->parents->datatarray().size(t));
        for (auto e = surface->parents->datatarray().begin(t); e != surface->parents->datatarray().end(t); ++e) {
            tcontact.push_back(elements->contact->datatarray()[*e]);
        }
        contact[t].swap(tcontact);
    }
    surface->contact = new serializededata<esint, ContactInfo>(1, contact);

    surface->size = surface->edistribution.back();
    surface->offset = surface->edistribution.back();
    surface->totalSize = Communication::exscan(surface->offset);

    DebugOutput::surface("surface.bodies", 1, 1);
    profiler::syncend("compute_bodies_surface");
    eslog::checkpointln("MESH: BODY SURFACE COMPUTED");
}

template <size_t nodes, size_t edim, size_t gps=1> struct SurfaceElementBasis {
    double w[gps], N[gps][nodes], dN[gps][nodes][edim], NN[nodes][nodes], dNN[nodes][nodes][edim];
};

template<Element::CODE code, size_t nodes>
static void surfaceNormal2D(SurfaceStore * surface, size_t e, const double* displacement)
{
    SurfaceElementBasis<nodes, 1> basis;
    BaseFunctions<code, 1>::simd(basis);
}

template<Element::CODE code, size_t nodes, size_t gps>
static void surfaceNormal3D(SurfaceStore * surface, size_t e, const double* displacement)
{
    SurfaceElementBasis<nodes, 2, gps> basis;
    BaseFunctions<code, gps>::simd(basis);

    auto enodes = surface->enodes->begin() + e;
    for (size_t n = 0; n < enodes->size(); ++n) {
        double dND[6] = { 0, 0, 0, 0, 0, 0 };
        for (size_t m = 0; m < enodes->size(); ++m) {
            Point disp;
            if (displacement) {
                int dim = info::mesh->dimension;
                disp.x = displacement[enodes->at(m) * dim + 0]; disp.y = displacement[enodes->at(m) * dim + 1]; if (dim == 3) disp.z = displacement[enodes->at(m) * dim + 2];
            }
            Point c = surface->coordinates->datatarray()[enodes->at(m)] + disp;
            dND[0] += basis.dNN[n][m][0] * c.x; dND[1] += basis.dNN[n][m][0] * c.y; dND[2] += basis.dNN[n][m][0] * c.z;
            dND[3] += basis.dNN[n][m][1] * c.x; dND[4] += basis.dNN[n][m][1] * c.y; dND[5] += basis.dNN[n][m][1] * c.z;
        }
        Point normal(
                dND[1] * dND[5] - dND[2] * dND[4],
                dND[2] * dND[3] - dND[0] * dND[5],
                dND[0] * dND[4] - dND[1] * dND[3]);
        normal.normalize();
        normal /= surface->nodeMultiplicity->data[enodes->at(n)];
        surface->nodeNormals->data[3 * enodes->at(n) + 0] += normal.x;
        surface->nodeNormals->data[3 * enodes->at(n) + 1] += normal.y;
        surface->nodeNormals->data[3 * enodes->at(n) + 2] += normal.z;
    }
}

void computeSurfaceNodeNormals(NodeStore *nodes, SurfaceStore * surface, const std::vector<int> &neighbors, const double* displacement)
{
    if (surface->nodeNormals == NULL) {
        surface->nodeNormals = nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "SURFACE_NORMAL", step::TYPE::TIME, false);
        surface->nodeMultiplicity = nodes->appendData(1, NamedData::DataType::SCALAR);

        for (auto enodes = surface->enodes->begin(); enodes != surface->enodes->end(); ++enodes) {
            for (auto n = enodes->begin(); n != enodes->end(); ++n) {
                surface->nodeMultiplicity->data[*n] += 1;
            }
        }
        surface->nodeMultiplicity->synchronize();
    }
    std::fill(surface->nodeNormals->data.begin(), surface->nodeNormals->data.end(), 0);

    for (size_t e = 0; e < surface->enodes->structures(); ++e) {
        switch (surface->epointers->datatarray()[e]->code) {
        case Element::CODE::LINE2    : surfaceNormal2D<Element::CODE::LINE2    , 2>(surface, e, displacement); break;
        case Element::CODE::LINE3    : surfaceNormal2D<Element::CODE::LINE3    , 3>(surface, e, displacement); break;
        case Element::CODE::TRIANGLE3: surfaceNormal3D<Element::CODE::TRIANGLE3, 3, 6>(surface, e, displacement); break;
        case Element::CODE::TRIANGLE6: surfaceNormal3D<Element::CODE::TRIANGLE6, 6, 6>(surface, e, displacement); break;
        case Element::CODE::SQUARE4  : surfaceNormal3D<Element::CODE::SQUARE4  , 4, 1>(surface, e, displacement); break;
        case Element::CODE::SQUARE8  : surfaceNormal3D<Element::CODE::SQUARE8  , 8, 1>(surface, e, displacement); break;
        default:
            eslog::internalFailure("unknown or not implemented surface element.\n");
        }
    }

    surface->nodeNormals->synchronize();
}

void computeWarpedNormals(SurfaceStore * surface, const double* displacement)
{
    if (surface->normal != NULL) {
        delete surface->normal;
    }
    if (surface->parameters != NULL) {
        delete surface->parameters;
    }
    if (surface->base != NULL) {
        delete surface->base;
    }
    profiler::syncstart("compute_warped_surface_normals");

    std::vector<std::vector<Point> > normal(info::env::OMP_NUM_THREADS), base(info::env::OMP_NUM_THREADS), parameters(info::env::OMP_NUM_THREADS);

    #pragma omp parallel for
    for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
        std::vector<Point> tnormal(surface->edistribution[t + 1] - surface->edistribution[t]);
        std::vector<Point> tparamters(2 * (surface->edistribution[t + 1] - surface->edistribution[t]));
        std::vector<Point> tbase(surface->edistribution[t + 1] - surface->edistribution[t]);

        auto epointers = surface->epointers->datatarray().begin(t);
        auto enodes = surface->enodes->begin(t);
        for (size_t e = surface->edistribution[t], i = 0; e < surface->edistribution[t + 1]; ++e, ++i, ++epointers, ++enodes) {
            switch ((*epointers)->code) {
//            case Element::CODE::LINE2:
//            case Element::CODE::LINE3:
//            {
//                const Point &a = surface->coordinates->datatarray()[enodes->at(0)];
//                const Point &b = surface->coordinates->datatarray()[enodes->at(1)];
//                tbase[i] = (a + b) / 2;
//                tnormal[i] = Point(b - a).normalize();
//                std::swap(tnormal[i].x, tnormal[i].y);
//                tnormal[i].y = -tnormal[i].y;
//            } break;
            case Element::CODE::TRIANGLE3:
//            case Element::CODE::TRIANGLE6:
            {
                Point da, db, dc;
                if (displacement) {
                    int dim = info::mesh->dimension;
                    da.x = displacement[dim * enodes->at(0) + 0]; da.y = displacement[dim * enodes->at(0) + 1]; if (dim == 3) da.z = displacement[dim * enodes->at(0) + 2];
                    db.x = displacement[dim * enodes->at(1) + 0]; db.y = displacement[dim * enodes->at(1) + 1]; if (dim == 3) db.z = displacement[dim * enodes->at(1) + 2];
                    dc.x = displacement[dim * enodes->at(2) + 0]; dc.y = displacement[dim * enodes->at(2) + 1]; if (dim == 3) dc.z = displacement[dim * enodes->at(2) + 2];
                }
                Point a = surface->coordinates->datatarray()[enodes->at(0)] + da;
                Point b = surface->coordinates->datatarray()[enodes->at(1)] + db;
                Point c = surface->coordinates->datatarray()[enodes->at(2)] + dc;
                tbase[i] = a;
                tparamters[2 * i] = b - a;
                tparamters[2 * i + 1] = c - a;
                tnormal[i] = Point::cross(tparamters[2 * i], tparamters[2 * i + 1]).normalize();
            } break;
            case Element::CODE::SQUARE4:
//            case Element::CODE::SQUARE8:
            {
                Point da, db, dc, dd;
                if (displacement) {
                    int dim = info::mesh->dimension;
                    da.x = displacement[dim * enodes->at(0) + 0]; da.y = displacement[dim * enodes->at(0) + 1]; if (dim == 3) da.z = displacement[dim * enodes->at(0) + 2];
                    db.x = displacement[dim * enodes->at(1) + 0]; db.y = displacement[dim * enodes->at(1) + 1]; if (dim == 3) db.z = displacement[dim * enodes->at(1) + 2];
                    dc.x = displacement[dim * enodes->at(2) + 0]; dc.y = displacement[dim * enodes->at(2) + 1]; if (dim == 3) dc.z = displacement[dim * enodes->at(2) + 2];
                    dd.x = displacement[dim * enodes->at(3) + 0]; dd.y = displacement[dim * enodes->at(3) + 1]; if (dim == 3) dd.z = displacement[dim * enodes->at(3) + 2];
                }
                Point a = surface->coordinates->datatarray()[enodes->at(0)] + da;
                Point b = surface->coordinates->datatarray()[enodes->at(1)] + db;
                Point c = surface->coordinates->datatarray()[enodes->at(2)] + dc;
                Point d = surface->coordinates->datatarray()[enodes->at(3)] + dd;
                Point center = (a + b + c + d) / 4;
                tnormal[i] = Point::cross(c - a, d - b).normalize();
                Point plane[4] = {
                        a - tnormal[i] * (tnormal[i] * (a - center)),
                        b - tnormal[i] * (tnormal[i] * (b - center)),
                        c - tnormal[i] * (tnormal[i] * (c - center)),
                        d - tnormal[i] * (tnormal[i] * (d - center))
                };
                Point pp[4] = { plane[1] - plane[0], plane[2] - plane[1], plane[3] - plane[2], plane[0] - plane[3] };
                int minv = 0; double mindiff = 2;
                for (int v = 0; v < 4; ++v) {
                    double _s, _t;
                    plane[(v + 2) % 4].getBarycentric(plane[v], pp[v], -pp[(v + 3) % 4], _s, _t);
                    double diff = 0;
                    if (_s < 0) diff += -_s;
                    if (_t < 0) diff += -_t;
                    if (_s > 0) diff += _s - 1;
                    if (_t > 0) diff += _t - 1;
                    if (mindiff > diff) {
                        mindiff = diff;
                        minv = v;
                    }
                }
                tbase[i] = plane[minv];
                tparamters[2 * i] = pp[minv];
                tparamters[2 * i + 1] = -pp[(minv + 3) % 4];
            } break;
            default:
                eslog::internalFailure("unknown or not implemented surface element.\n");
            }
        }

        normal[t].swap(tnormal);
        parameters[t].swap(tparamters);
        base[t].swap(tbase);
    }

    surface->normal = new serializededata<esint, Point>(1, normal);
    surface->parameters = new serializededata<esint, Point>(2, parameters);
    surface->base = new serializededata<esint, Point>(1, base);

    DebugOutput::warpedNormals("surface.planes", displacement, 1, 1);
    profiler::syncend("compute_warped_surface_normals");
    eslog::checkpointln("MESH: WARPED SURFACE NORMALS COMMPUTED");
}

void exchangeContactHalo(SurfaceStore * surface, ContactStore *contact, const double* displacement)
{
    profiler::syncstart("exchange_contact_halo");
    _Point<float> box[2] = {
            _Point<float>( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()),
            _Point<float>(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max())
    };
    auto enodes = surface->enodes->begin();
    for (esint e = 0; e < surface->size; ++e, ++enodes) {
        if (surface->contact->datatarray()[e].gap > 0) {
            for (auto n = enodes->begin(); n != enodes->end(); ++n) {
                Point disp;
                if (displacement) {
                    int dim = info::mesh->dimension;
                    disp.x = displacement[*n * dim + 0]; disp.y = displacement[*n * dim + 1]; if (dim == 3) disp.z = displacement[*n * dim + 2];
                }
                Point p = surface->coordinates->datatarray()[*n] + disp;
                box[0].x = std::min((float)p.x - surface->contact->datatarray()[e].gap, box[0].x);
                box[0].y = std::min((float)p.y - surface->contact->datatarray()[e].gap, box[0].y);
                box[0].z = std::min((float)p.z - surface->contact->datatarray()[e].gap, box[0].z);
                box[1].x = std::max((float)p.x + surface->contact->datatarray()[e].gap, box[1].x);
                box[1].y = std::max((float)p.y + surface->contact->datatarray()[e].gap, box[1].y);
                box[1].z = std::max((float)p.z + surface->contact->datatarray()[e].gap, box[1].z);
            }
        }
    }

    std::vector<_Point<float> > boxes(2 * info::mpi::size);
    Communication::allGather(box, boxes.data(), 6, MPI_FLOAT);

    auto areIntersected = [&] (_Point<float> *block1, _Point<float> *block2) {
        return    !
                ((block1[1].x < block2[0].x || block2[1].x < block1[0].x) ||
                 (block1[1].y < block2[0].y || block2[1].y < block1[0].y) ||
                 (block1[1].z < block2[0].z || block2[1].z < block1[0].z));
    };

    contact->neighbors.clear();
    for (int r = 0; r < info::mpi::size; r++) {
        if (r != info::mpi::rank) {
            if (areIntersected(box, boxes.data() + 2 * r)) {
                contact->neighbors.push_back(r);
            }
        }
    }

    const auto &coordinates = surface->coordinates->datatarray();
    std::vector<std::vector<esint> > sBuffer(contact->neighbors.size()), rBuffer(contact->neighbors.size());
    for (size_t n = 0; n < contact->neighbors.size(); ++n) {
        std::vector<esint> esend, nsend;
        const auto &min = boxes[2 * contact->neighbors[n]];
        const auto &max = boxes[2 * contact->neighbors[n] + 1];
        auto enodes = surface->enodes->begin();
        for (size_t e = 0; e < surface->enodes->structures(); ++e, ++enodes) {
            for (auto ec = enodes->begin(); ec != enodes->end(); ++ec) {
                Point disp;
                if (displacement) {
                    int dim = info::mesh->dimension;
                    disp.x = displacement[*ec * dim + 0]; disp.y = displacement[*ec * dim + 1]; if (dim == 3) disp.z = displacement[*ec * dim + 2];
                }
                if (
                        (min.x <= coordinates[*ec].x + disp.x && coordinates[*ec].x + disp.x <= max.x) &&
                        (min.y <= coordinates[*ec].y + disp.y && coordinates[*ec].y + disp.y <= max.y) &&
                        (min.z <= coordinates[*ec].z + disp.z && coordinates[*ec].z + disp.z <= max.z)) {

                    esend.push_back(e);
                    for (auto c = enodes->begin(); c != enodes->end(); ++c) {
                        nsend.push_back(*c);
                    }
                    break;
                }
            }
        }
        size_t ssize = 0, nsize = nsend.size();
        ssize += 1 + 3 * esend.size(); // fID, body, epointer size
        ssize += 1 + esend.size() + nsend.size(); // enodes size
        ssize += 1 + esend.size() * (4 * sizeof(Point) / sizeof(esint)); // plane
        utils::sortAndRemoveDuplicates(nsend);
        ssize += 1 + nsend.size() * (1 + sizeof(Point) / sizeof(esint)); // ids + coordinates size

        sBuffer[n].reserve(ssize);

        // send offset and bodies
        sBuffer[n].push_back(esend.size());
        for (size_t e = 0; e < esend.size(); ++e) {
            sBuffer[n].push_back(esend[e]);
        }
        for (size_t e = 0; e < esend.size(); ++e) {
            sBuffer[n].push_back(surface->body->datatarray()[esend[e]]);
        }
        for (size_t e = 0; e < esend.size(); ++e) {
            sBuffer[n].push_back(static_cast<int>(surface->epointers->datatarray()[esend[e]]->code));
        }

        // send enodes in target offsets
        sBuffer[n].push_back(esend.size() + nsize);
        enodes = surface->enodes->begin();
        for (size_t e = 0, prev = 0; e < esend.size(); prev = esend[e++]) {
            enodes += esend[e] - prev;
            sBuffer[n].push_back(enodes->size());
            for (auto c = enodes->begin(); c != enodes->end(); ++c) {
                sBuffer[n].push_back(std::lower_bound(nsend.begin(), nsend.end(), *c) - nsend.begin());
            }
        }

        // send planes
        sBuffer[n].push_back(esend.size());
        for (size_t e = 0; e < esend.size(); ++e) {
            const auto &nn = surface->normal->datatarray()[esend[e]];
            sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nn), reinterpret_cast<const esint*>(&nn) + sizeof(nn) / sizeof(esint));
        }
        for (size_t e = 0; e < esend.size(); ++e) {
            const auto &nn = surface->parameters->datatarray()[2 * esend[e]];
            sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nn), reinterpret_cast<const esint*>(&nn) + 2 * sizeof(nn) / sizeof(esint));
        }
        for (size_t e = 0; e < esend.size(); ++e) {
            const auto &nc = surface->base->datatarray()[esend[e]];
            sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&nc), reinterpret_cast<const esint*>(&nc) + sizeof(nc) / sizeof(esint));
        }

        // send global ids
        sBuffer[n].push_back(nsend.size());
        for (size_t c = 0; c < nsend.size(); ++c) {
            sBuffer[n].push_back(surface->nIDs->datatarray()[nsend[c]]);
        }

        // send coordinates
        for (size_t c = 0; c < nsend.size(); ++c) {
            Point disp;
            if (displacement) {
                int dim = info::mesh->dimension;
                disp.x = displacement[nsend[c] * dim + 0]; disp.y = displacement[nsend[c] * dim + 1]; if (dim == 3) disp.z = displacement[nsend[c] * dim + 2];
            }
            Point p = surface->coordinates->datatarray()[nsend[c]] + disp;
            sBuffer[n].insert(sBuffer[n].end(), reinterpret_cast<const esint*>(&p), reinterpret_cast<const esint*>(&p) + sizeof(p) / sizeof(esint));
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, contact->neighbors)) {
        eslog::internalFailure("cannot exchange contact halo.\n");
    }

    contact->surfaces.resize(contact->neighbors.size());
    for (size_t n = 0, i = 0; n < contact->neighbors.size(); ++n, i = 0) {
        contact->surfaces[n] = new SurfaceStore();
        // receive parents and bodies
        contact->surfaces[n]->fID = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
        memcpy(contact->surfaces[n]->fID->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
        contact->surfaces[n]->body = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
        memcpy(contact->surfaces[n]->body->datatarray().data(), rBuffer[n].data() + i + 1 + rBuffer[n][i], sizeof(esint) * rBuffer[n][i]);
        contact->surfaces[n]->epointers = new serializededata<esint, Element*>(1, tarray<Element*>(info::env::OMP_NUM_THREADS, rBuffer[n][i]));
        for (esint e = 0; e < rBuffer[n][i]; ++e) {
            contact->surfaces[n]->epointers->datatarray()[e] = &Mesh::edata[rBuffer[n][2 * rBuffer[n][i] + e + i + 1]];
        }
        i += 1 + 3 * rBuffer[n][i];

        // receive enodes
        std::vector<std::vector<esint> > edist(info::env::OMP_NUM_THREADS), enodes(info::env::OMP_NUM_THREADS);
        size_t size = i + rBuffer[n][i];
        ++i;
        edist[0].push_back(0);
        while(i < size) {
            size_t nsize = rBuffer[n][i++];
            for (size_t nn = 0; nn < nsize; ++nn) {
                enodes[0].push_back(rBuffer[n][i++]);
            }
            edist[0].push_back(enodes[0].size());
        }
        serializededata<esint, esint>::balance(edist, enodes, &contact->surfaces[n]->fID->datatarray().distribution());
        contact->surfaces[n]->enodes = new serializededata<esint, esint>(edist, enodes);

        // receive warped normals
        size = rBuffer[n][i];
        contact->surfaces[n]->normal = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
        memcpy(contact->surfaces[n]->normal->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i + 1), sizeof(Point) * size);
        i += 1 + size * (sizeof(Point) / sizeof(esint));

        contact->surfaces[n]->parameters = new serializededata<esint, Point>(2, tarray<Point>(info::env::OMP_NUM_THREADS, 2 * size));
        memcpy(contact->surfaces[n]->parameters->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), 2 * sizeof(Point) * size);
        i += size * (2 * sizeof(Point) / sizeof(esint));

        contact->surfaces[n]->base = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
        memcpy(contact->surfaces[n]->base->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), sizeof(Point) * size);
        i += size * (sizeof(Point) / sizeof(esint));

        // receive global ids
        size = rBuffer[n][i];
        contact->surfaces[n]->nIDs = new serializededata<esint, esint>(1, tarray<esint>(info::env::OMP_NUM_THREADS, size));
        memcpy(contact->surfaces[n]->nIDs->datatarray().data(), rBuffer[n].data() + i + 1, sizeof(esint) * rBuffer[n][i]);
        i += 1 + rBuffer[n][i];

        // receive coordinates
        contact->surfaces[n]->coordinates = new serializededata<esint, Point>(1, tarray<Point>(info::env::OMP_NUM_THREADS, size));
        memcpy(contact->surfaces[n]->coordinates->datatarray().data(), reinterpret_cast<const Point*>(rBuffer[n].data() + i), sizeof(Point) * size);
    }
    contact->neighborsWithMe = contact->neighbors;
    contact->neighborsWithMe.push_back(info::mpi::rank);
    contact->surfaces.push_back(surface);

    profiler::syncend("exchange_contact_halo");
    eslog::checkpointln("MESH: CLOSE BOUNDARY EXCHANGED");
}

void findCloseElements(ContactStore *contact, const double* displacement)
{
    profiler::syncstart("find_close_elements");
    // checking
    // 1. distance to plane defined by normal and center
    // 2. test if any point is in coarse element defined by: base + s * parameter.u + t * parameter.v

    auto getCoords = [&] (size_t neigh, size_t n) {
        Point p = contact->surfaces[neigh]->coordinates->datatarray()[n];
        if (displacement && neigh + 1 == contact->surfaces.size()) {
            int dim = info::mesh->dimension;
            p.x += displacement[n * dim + 0];
            p.y += displacement[n * dim + 1];
            if (dim == 3) p.z += displacement[n * dim + 2];
        }
        return p;
    };

    // neighbors, local
    std::vector<Point> nstart, nend, lstart, lend;
    std::vector<esint> neigh, offset = { 0 };

    for (size_t n = 0; n < contact->neighborsWithMe.size() - 1; ++n) {
        offset.push_back(offset.back() + contact->surfaces[n]->enodes->structures());
    }
    // push neighbors
    nstart.reserve(offset.back());
    nend.reserve(offset.back());
    neigh.reserve(offset.back());
    for (size_t r = 0, offset = 0; r < contact->neighborsWithMe.size() - 1; ++r) {
        for (auto e = contact->surfaces[r]->enodes->begin(); e != contact->surfaces[r]->enodes->end(); ++e, ++offset) {
            nstart.push_back(getCoords(r, e->front()));
            nend.push_back(nstart.back());
            for (auto n = e->begin() + 1; n != e->end(); ++n) {
                getCoords(r, *n).minmax(nstart.back(), nend.back());
            }
        }
        neigh.insert(neigh.end(), contact->surfaces[r]->enodes->structures(), r);
    }
    // push local
    size_t local_n = contact->surfaces.size() - 1;
    lstart.reserve(contact->surfaces.back()->enodes->structures());
    lend.reserve(contact->surfaces.back()->enodes->structures());
    for (auto e = contact->surfaces.back()->enodes->begin(); e != contact->surfaces.back()->enodes->end(); ++e) {
        lstart.push_back(getCoords(local_n, e->front()));
        lend.push_back(lstart.back());
        for (auto n = e->begin() + 1; n != e->end(); ++n) {
            getCoords(local_n, *n).minmax(lstart.back(), lend.back());
        }
    }

    IntervalTree ntree(nstart, nend);
    IntervalTree ltree(lstart, lend);

    std::unordered_map<esint, std::unordered_set<esint> > bodyPairs;
    std::vector<esint> dist = { 0 }, data;
    dist.reserve(contact->surfaces.back()->size + 1);
    data.reserve(offset.back() + contact->surfaces.back()->size);
    std::vector<esint> nintervals, lintervals;

    for (esint face = 0; face < contact->surfaces.back()->size; ++face) {
        bool self_contact = contact->surfaces.back()->contact->datatarray()[face].self_contact;
        double gap = contact->surfaces.back()->contact->datatarray()[face].gap;
        double angle = -std::cos(M_PI * contact->surfaces.back()->contact->datatarray()[face].angle / 180);

        int body = contact->surfaces.back()->body->datatarray()[face];
        Point normal = contact->surfaces.back()->normal->datatarray()[face];
        Point base = contact->surfaces.back()->base->datatarray()[face];
        Point u = contact->surfaces.back()->parameters->datatarray()[2 * face];
        Point v = contact->surfaces.back()->parameters->datatarray()[2 * face + 1];
        bool istriangle =
                (contact->surfaces.back()->epointers->datatarray()[face]->code == Element::CODE::TRIANGLE3) ||
                (contact->surfaces.back()->epointers->datatarray()[face]->code == Element::CODE::TRIANGLE6);

        auto cohen_sutherland = [&] (int n, esint index) {
            auto getCode = [&] (const double &d, const double &s, const double &t) {
                int code = 0;
                code |= (s < 0) ? 1 : 0;
                code |= (t < 0) ? 2 : 0;
                if (istriangle) {
                    code |= (1 < t + s) ? 4 : 0;
                } else {
                    code |= (1 < s) ? 4 : 0;
                    code |= (1 < t) ? 8 : 0;
                }
                code |= d < -gap ? 16 : 0;
                code |= gap < d  ? 32 : 0;
                return code;
            };

            int code = 63;
            auto nodes = contact->surfaces[n]->enodes->begin() + index;
            for (auto nn = nodes->begin(); nn != nodes->end(); ++nn) {
                Point p = getCoords(n, *nn);
                double s, t, d = normal * (p - base);
                p -= normal * d;
                p.getBarycentric(base, u, v, s, t);
                code &= getCode(d, s, t);
            }

            return code == 0;
        };

        Point gapdirection = normal;
        gapdirection.abs();
        nintervals.clear();
        lintervals.clear();
        ntree.traverse(1, lstart[face] - gapdirection * gap, lend[face] + gapdirection * gap, nintervals);
        ltree.traverse(1, lstart[face] - gapdirection * gap, lend[face] + gapdirection * gap, lintervals);

        for (size_t i = 0; i < nintervals.size(); ++i) {
            for (auto pp = ntree.permutation.cbegin() + ntree.begin(nintervals[i]); pp != ntree.permutation.cbegin() + ntree.end(nintervals[i]); ++pp) {
                int n = neigh[*pp];
                esint index = *pp - offset[n];
                int nbody = contact->surfaces[n]->body->datatarray()[index];
                const Point &nnormal = contact->surfaces[n]->normal->datatarray()[index];
                if ((self_contact || nbody != body) && nnormal * normal <= angle && cohen_sutherland(n, index)) {
                    bodyPairs[std::min(body, nbody)].insert(std::max(body, nbody));
                    data.push_back(n);
                    data.push_back(index);
                }
            }
        }

        for (size_t i = 0; i < lintervals.size(); ++i) {
            for (auto pp = ltree.permutation.cbegin() + ltree.begin(lintervals[i]); pp != ltree.permutation.cbegin() + ltree.end(lintervals[i]); ++pp) {
                int nbody = contact->surfaces.back()->body->datatarray()[*pp];
                const Point &nnormal = contact->surfaces.back()->normal->datatarray()[*pp];
                if ((self_contact || nbody != body) && nnormal * normal <= angle && cohen_sutherland(local_n, *pp)) {
                    bodyPairs[std::min(body, nbody)].insert(std::max(body, nbody));
                    data.push_back(contact->surfaces.size() - 1);
                    data.push_back(*pp);
                }
            }
        }

        dist.push_back(data.size());
    }

    contact->pairs = new serializededata<esint, esint>(dist, data);

    std::vector<std::pair<esint, esint> > interfaces;
    for (auto i = bodyPairs.begin(); i != bodyPairs.end(); ++i) {
        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            interfaces.push_back(std::make_pair(i->first, *j));
        }
    }
    std::sort(interfaces.begin(), interfaces.end());
    Communication::uniqueAllGatherUnknownSize(interfaces);
    for (auto i = interfaces.begin(); i != interfaces.end(); ++i) {
        contact->interfaces.push_back(Interface(i->first, i->second));
    }

    profiler::syncend("find_close_elements");
    eslog::checkpointln("MESH: CLOSE ELEMENTS FOUND");
}

static void triangulate(std::vector<Point> &face, std::vector<Triangle> &triangles)
{
    triangles.clear();
    auto prev= [&] (const int &i) -> int {
        return (i + face.size() - 1) % face.size();
    };
    auto next = [&] (const int &i) -> int {
        return (i + 1) % face.size();
    };

    if (face.size() == 3) {
        triangles.push_back(Triangle{ face, 0, 1, 2 });
        return;
    }
    if (face.size() == 4) {
        for (int v = 0, size = face.size(); v < size; ++v) {
            Point v1 = face[v] - face[prev(v)], v2 = face[next(v)] - face[prev(v)];
            if (Point::cross2d(v1.normalize(), v2.normalize()) < 0) {
                triangles.push_back(Triangle{ face, v, next(v), next(next(v)) });
                triangles.push_back(Triangle{ face, next(next(v)), prev(v), v });
                return;
            }
        }
        if ((face[0] - face[2]).length() < (face[1] - face[3]).length()) {
            triangles.push_back(Triangle{ face, 0, 1, 2 });
            triangles.push_back(Triangle{ face, 0, 2, 3 });
        } else {
            triangles.push_back(Triangle{ face, 0, 1, 3 });
            triangles.push_back(Triangle{ face, 1, 2, 3 });
        }
        return;
    }
    eslog::internalFailure("cannot compute triangles from a polygon.\n");
}

static void clip(const Point &base, const std::vector<Triangle> &triangles, const std::vector<Point> &polygon, const double &gap, std::vector<Triangle> &output)
{
    double eps = 1e-10;
    enum status { in, out, on, processed };
    std::vector<Triangle> res;

    auto crossPoint = [] (const Point &p0, const Point &p1, const Point &q0, const Point &q1) {
        Point u = p1 - p0, v = q1 - q0, w = p0 - q0;
        double t = Point::cross2d(u, w) / Point::cross2d(u, v);
        return q0 + v * t;
    };

    for (size_t i = 0; i < triangles.size(); ++i) {
        std::vector<Point> in(polygon.rbegin(), polygon.rend()), nextin;
        auto curr = [&] (size_t i) { return (i                ) % in.size(); };
        auto next = [&] (size_t i) { return (i + 1            ) % in.size(); };
        auto prev = [&] (size_t i) { return (i - 1 + in.size()) % in.size(); };

        Point p[3] = { triangles[i].p[0], triangles[i].p[1], triangles[i].p[2] };
        Point line[3] = { (p[1] - p[0]).normalize(), (p[2] - p[1]).normalize(), (p[0] - p[2]).normalize() };

        for (int c = 0; c < 3; ++c) {
            std::vector<status> status(in.size());
            for (size_t j = 0; j < in.size(); ++j) {
                double distance = Point::cross2d(line[c], in[j] - p[c]);
                status[j] = eps < distance ? status::in : distance < -eps ? status::out : status::on;
            }
            for (size_t j = 0; j < in.size(); ++j) {
                if (status[j] == status::in && status[next(j)] == status::on) {
                    while (status[j = next(j)] == status::on) {
                        status[j] = status::in;
                    }
                }
            }
            for (size_t j = 0; j < in.size(); ++j) {
                if (status[j] == status::in && status[prev(j)] == status::on) {
                    while (status[j = prev(j)] == status::on) {
                        status[j] = status::in;
                    }
                }
            }
            for (size_t j = 0; j < in.size(); ++j) {
                if (status[j] == status::in && status[prev(j)] == status::out) {
                    nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j)], in[j]));
                    status[prev(j)] = status::processed;
                    for (size_t k = 0; k < in.size(); ++k) {
                        if (status[curr(j + k)] == status::in || status[curr(j + k)] == status::on) {
                            nextin.push_back(in[curr(j + k)]);
                            status[curr(j + k)] = status::processed;
                        } else {
                            nextin.push_back(crossPoint(p[c], p[(c + 1) % 3], in[prev(j + k)], in[curr(j + k)]));
                            status[curr(j + k)] = status::processed;
                            break;
                        }
                    }
                }
            }
            if (nextin.size()) {
                in = nextin;
                nextin.clear();
            } else if (in.size() && status.front() != status::in) {
                in.clear();
            }
        }
        if (in.size()) {
            std::vector<Point> _res;
            for (size_t j = 0; j < in.size(); ++j) {
                if (std::fabs(in[prev(j)].x - in[j].x) > 1e-6 || std::fabs(in[prev(j)].y - in[j].y) > 1e-6) {
                    _res.push_back(in[j]);
                }
            }
            if (3 <= _res.size()) {
                for (esint j = 2, size = _res.size(); j < size; j++) { // Delaunay triangulation?
                    res.push_back(Triangle{ _res, 0, j - 1, j });
                }
            }
        }
    }

    // remove triangles that are above the limit
    auto max = [&] (const double &z, const Point &p0, const Point &p1) { return p0 + (p1 - p0) * ((z - p0.z) / (p1.z - p0.z)); };

    output.clear();
    for (size_t t = 0; t < res.size(); ++t) {
        bool below[3] = { res[t].p[0].z + gap < base.z, res[t].p[1].z + gap < base.z, res[t].p[2].z + gap < base.z };

        if (!below[0] && !below[1] && !below[2]) {
            output.push_back(res[t]);
        } else if (!below[0] || !below[1] || !below[2]) {
            int i = 0;
            while (below[0] || !below[2]) { std::rotate(below, below + 1, below + 3); ++i; }
            if (below[1]) {
                output.push_back(Triangle(res[t].p[i], max(base.z - gap, res[t].p[i], res[t].p[(i + 1) % 3]), max(base.z - gap, res[t].p[i], res[t].p[(i + 2) % 3])));
            } else {
                Point c0 = max(base.z - gap, res[t].p[i], res[t].p[(i + 2) % 3]);
                Point c1 = max(base.z - gap, res[t].p[(i + 1) % 3], res[t].p[(i + 2) % 3]);
                output.push_back(Triangle(res[t].p[i], res[t].p[(i + 1) % 3], c0));
                output.push_back(Triangle(res[t].p[(i + 1) % 3], c1, c0));
            }
        }
    }
    output.swap(res);
    output.clear();
    for (size_t t = 0; t < res.size(); ++t) {
        bool above[3] = { base.z < res[t].p[0].z - gap, base.z < res[t].p[1].z - gap, base.z < res[t].p[2].z - gap };

        if (!above[0] && !above[1] && !above[2]) {
            output.push_back(res[t]);
        } else if (!above[0] || !above[1] || !above[2]) {
            int i = 0;
            while (above[0] || !above[2]) { std::rotate(above, above + 1, above + 3); ++i; }
            if (above[1]) {
                output.push_back(Triangle(res[t].p[i], max(base.z + gap, res[t].p[i], res[t].p[(i + 1) % 3]), max(base.z + gap, res[t].p[i], res[t].p[(i + 2) % 3])));
            } else {
                Point c0 = max(base.z + gap, res[t].p[i], res[t].p[(i + 2) % 3]);
                Point c1 = max(base.z + gap, res[t].p[(i + 1) % 3], res[t].p[(i + 2) % 3]);
                output.push_back(Triangle(res[t].p[i], res[t].p[(i + 1) % 3], c0));
                output.push_back(Triangle(res[t].p[(i + 1) % 3], c1, c0));
            }
        }
    }
}

void computeContactInterface(SurfaceStore* surface, ContactStore* contact, const double *displacement)
{
    profiler::syncstart("compute_contact_interface");

    const std::vector<SurfaceStore*> &surfaces = contact->surfaces;
    auto getCoords = [&] (size_t neigh, size_t n) {
        Point p = surfaces[neigh]->coordinates->datatarray()[n];
        if (displacement && neigh + 1 == contact->surfaces.size()) {
            int dim = info::mesh->dimension;
            p.x += displacement[n * dim + 0];
            p.y += displacement[n * dim + 1];
            if (dim == 3) p.z += displacement[n * dim + 2];
        }
        return p;
    };

    Point axis;
    double cos, sin;
    auto setRotation = [&] (esint e) {
        const Point &normal = surface->normal->datatarray()[e];
        if (normal.z < -0.999) {
            axis.x = 1;
            sin = axis.y = axis.z = 0;
            cos = -1;
        } else if (normal.z < 0.999) {
            axis = Point::cross(normal, Point(0, 0, 1)).normalize();
            double angle = std::acos(normal.z);
            sin = std::sin(angle);
            cos = std::cos(angle);
        } else {
            sin = axis.x = axis.y = 0;
            cos = axis.z = 1;
        }
    };

    std::vector<Point> plane(20), projected(20);
    std::vector<Triangle> planeTriangles, intersection;
    auto setPolygon = [&] (std::vector<Point> &polygon, esint neigh, esint offset) {
        polygon.clear();
        auto nodes = surfaces[neigh]->enodes->cbegin() + offset;
        const auto &epointer = surfaces[neigh]->epointers->datatarray();
        for (auto n = epointer[offset]->polygon->begin(); n != epointer[offset]->polygon->end(); ++n) {
            polygon.push_back(getCoords(neigh, nodes->at(*n)));
            polygon.back().rodrigues(axis, cos, sin);
        }
    };

    std::vector<Triangle> triangles;

    std::vector<esint> sdist = { 0 }, ddist = { 0 }, pdist = { 0 };
    std::vector<SparseSegment> sparse;
    std::vector<DenseSegment> dense;
    std::vector<Point2D> planeCoordinates;

    auto pairs = contact->pairs->cbegin();
    for (esint e = 0; e < surface->size; ++e, ++pairs) {
        if (pairs->size() == 0) {
            continue;
        }

        setRotation(e);

        double earea = 0;
        setPolygon(plane, contact->neighbors.size(), e);
        Point base = surface->base->datatarray()[e];
        base.rodrigues(axis, cos, sin);
        triangulate(plane, planeTriangles);
        for (size_t t = 0; t < planeTriangles.size(); ++t) {
            earea += planeTriangles[t].area();
        }

        std::unordered_map<esint, double> insertedBodies;
        for (auto other = pairs->begin(); other != pairs->end(); ++other) {
            esint neigh = *other++;
            esint offset = *other;

            setPolygon(projected, neigh, *other);
            clip(base, planeTriangles, projected, surface->contact->datatarray()[e].gap, intersection);
            if (intersection.size()) {
                if (insertedBodies.empty()) {
                    sparse.push_back(SparseSegment(surfaces.back()->body->datatarray()[e], e, planeCoordinates.size(), triangles.size(), dense.size()));
                    planeCoordinates.insert(planeCoordinates.end(), plane.begin(), plane.end());
                }
                if (insertedBodies.count(surfaces[neigh]->body->datatarray()[offset]) == 0) {
                    insertedBodies[surfaces[neigh]->body->datatarray()[offset]] = 0;
                }

                ++sparse.back().denseSegmentEnd;
                dense.push_back(DenseSegment(neigh, surfaces[neigh]->body->datatarray()[offset], offset, planeCoordinates.size(), planeCoordinates.size() + projected.size()));
                planeCoordinates.insert(planeCoordinates.end(), projected.begin(), projected.end());
                for (size_t i = 0; i < intersection.size(); i++) {
                    planeCoordinates.insert(planeCoordinates.end(), intersection[i].p, intersection[i].p + 3);
                }
                dense.back().triangles += intersection.size();
            }
            for (size_t i = 0; i < intersection.size(); i++) {
                double area = intersection[i].area();
                insertedBodies[surfaces[neigh]->body->datatarray()[offset]] += area / earea;
                intersection[i].rotate(axis, cos, -sin);
                triangles.push_back(intersection[i]);
            }
        }
        for (auto ib = insertedBodies.begin(); ib != insertedBodies.end(); ++ib) {
            if (ib->second < MIN_SLAVE_COVER_RATIO) {
                for (esint d = sparse.back().denseSegmentBegin; d < sparse.back().denseSegmentEnd; ++d) {
                    if (surfaces[dense[d].neigh]->body->datatarray()[dense[d].element] == ib->first) {
                        dense[d].skip = true;
                    }
                }
            }
        }
        if (insertedBodies.size()) {
            sdist.push_back(sparse.size());
            ddist.push_back(dense.size());
            pdist.push_back(planeCoordinates.size());
        }
    }

    contact->sparseSide = new serializededata<esint, SparseSegment>(sdist, sparse);
    contact->denseSide = new serializededata<esint, DenseSegment>(ddist, dense);
    contact->planeCoordinates = new serializededata<esint, Point2D>(pdist, planeCoordinates);
    contact->intersections = new serializededata<esint, Triangle>(1, triangles);

    profiler::syncend("compute_contact_interface");
    eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
}

template<Element::CODE code, size_t nodes>
static void normal2D(NodeData *nodeMultiplicity, NodeData *nodeNormals, serializededata<esint, Point> *coordinates, BoundaryRegionStore *region, size_t e, const double* displacement)
{
    SurfaceElementBasis<nodes, 1> basis;
    BaseFunctions<code, 1>::simd(basis);
}

template<Element::CODE code, size_t nodes, size_t gps>
static void normal3D(NodeData *nodeMultiplicity, NodeData *nodeNormals, serializededata<esint, Point> *coordinates, BoundaryRegionStore *region, size_t e, const double* displacement)
{
    SurfaceElementBasis<nodes, 2, gps> basis;
    BaseFunctions<code, gps>::simd(basis);

    auto enodes = region->elements->begin() + e;
    for (size_t n = 0; n < enodes->size(); ++n) {
        double dND[6] = { 0, 0, 0, 0, 0, 0 };
        for (size_t m = 0; m < enodes->size(); ++m) {
            Point disp;
            if (displacement) {
                int dim = info::mesh->dimension;
                disp.x = displacement[enodes->at(m) * dim + 0]; disp.y = displacement[enodes->at(m) * dim + 1]; if (dim == 3) disp.z = displacement[enodes->at(m) * dim + 2];
            }
            Point c = coordinates->datatarray()[enodes->at(m)] + disp;
            dND[0] += basis.dNN[n][m][0] * c.x; dND[1] += basis.dNN[n][m][0] * c.y; dND[2] += basis.dNN[n][m][0] * c.z;
            dND[3] += basis.dNN[n][m][1] * c.x; dND[4] += basis.dNN[n][m][1] * c.y; dND[5] += basis.dNN[n][m][1] * c.z;
        }
        Point normal(
                dND[1] * dND[5] - dND[2] * dND[4],
                dND[2] * dND[3] - dND[0] * dND[5],
                dND[0] * dND[4] - dND[1] * dND[3]);
        normal.normalize();
        normal /= nodeMultiplicity->data[enodes->at(n)];
        nodeNormals->data[3 * enodes->at(n) + 0] += normal.x;
        nodeNormals->data[3 * enodes->at(n) + 1] += normal.y;
        nodeNormals->data[3 * enodes->at(n) + 2] += normal.z;
    }
}

void computeBoundaryRegionNormals(NodeStore *nodes, std::vector<BoundaryRegionStore*> &boundaryRegions, const std::vector<int> &neighbors, const double* displacement)
{
    for (size_t r = 1; r < boundaryRegions.size(); ++r) {
        if (boundaryRegions[r]->dimension) {
            if (boundaryRegions[r]->nodeNormals == nullptr) {
                boundaryRegions[r]->nodeNormals = nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR);
                boundaryRegions[r]->nodeMultiplicity = nodes->appendData(1, NamedData::DataType::SCALAR);
                for (auto enodes = boundaryRegions[r]->elements->begin(); enodes != boundaryRegions[r]->elements->end(); ++enodes) {
                    for (auto n = enodes->begin(); n != enodes->end(); ++n) {
                        boundaryRegions[r]->nodeMultiplicity->data[*n] += 1;
                    }
                }
                boundaryRegions[r]->nodeMultiplicity->synchronize();
            }
            std::fill(boundaryRegions[r]->nodeNormals->data.begin(), boundaryRegions[r]->nodeNormals->data.end(), 0);

            for (size_t e = 0; e < boundaryRegions[r]->elements->structures(); ++e) {
                switch (boundaryRegions[r]->epointers->datatarray()[e]->code) {
                case Element::CODE::LINE2    : normal2D<Element::CODE::LINE2    , 2>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                case Element::CODE::LINE3    : normal2D<Element::CODE::LINE3    , 3>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                case Element::CODE::TRIANGLE3: normal3D<Element::CODE::TRIANGLE3, 3, 6>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                case Element::CODE::TRIANGLE6: normal3D<Element::CODE::TRIANGLE6, 6, 6>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                case Element::CODE::SQUARE4  : normal3D<Element::CODE::SQUARE4  , 4, 1>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                case Element::CODE::SQUARE8  : normal3D<Element::CODE::SQUARE8  , 8, 1>(boundaryRegions[r]->nodeMultiplicity, boundaryRegions[r]->nodeNormals, nodes->coordinates, boundaryRegions[r], e, displacement); break;
                default:
                    eslog::internalFailure("unknown or not implemented surface element.\n");
                }
            }
            DebugOutput::regionNormals(boundaryRegions[r], displacement);
        }
    }
}

void arrangeContactInterfaces(NodeStore *nodes, ContactStore* contact, BodyStore *bodies, std::vector<ElementsRegionStore*> &elementsRegions, std::vector<ContactInterfaceStore*> &contactInterfaces, const double* displacement)
{
    profiler::syncstart("arrange_contact_interface");

    if (contact->nodeMultiplicity == NULL) {
        contact->nodeMultiplicity = nodes->appendData(1, NamedData::DataType::SCALAR);
    }
    if (contact->nodeNormals == NULL) {
        contact->nodeNormals = nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "CONTACT_NORMAL");
    }
    std::fill(contact->nodeMultiplicity->data.begin(), contact->nodeMultiplicity->data.end(), 0);
    std::fill(contact->nodeNormals->data.begin(), contact->nodeNormals->data.end(), 0);

    auto eregion = [elementsRegions] (const std::string &name) -> ElementsRegionStore* {
        for (size_t r = 0; r < elementsRegions.size(); r++) {
            if (StringCompare::caseSensitiveEq(elementsRegions[r]->name, name)) {
                return elementsRegions[r];
            }
        }
        eslog::error("Unknown region of elements with name '%s'.\n", name.c_str());
        return NULL;
    };

    const std::vector<SurfaceStore*> &surfaces = contact->surfaces;

    struct istats { double area; esint faces, triangles; bool skip; istats(): area(0), faces(0), triangles(0), skip(true) {} };
    std::unordered_map<esint, std::unordered_map<esint, istats> > istats;

    auto *sside = contact->sparseSide;
    auto *dside = contact->denseSide;
    auto &coors = contact->planeCoordinates->datatarray();
    for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
        bool add = true;
        for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
            if (d->skip) {
                continue;
            }

            if (add) {
                add = false;
                istats[s->body][d->body].faces += 1;
            }
            istats[s->body][d->body].triangles += d->triangles;

            for (esint t = 0; t < d->triangles; ++t) {
                istats[s->body][d->body].area += Triangle::area(coors.data() + d->triangleOffset + 3 * t);
            }
        }
    }

    std::vector<Interface> interfaces = contact->interfaces;
    auto comp = [] (const Interface &i, const Interface &j) {
        if (i.from.body == j.from.body) {
            return i.to.body < j.to.body;
        }
        return i.from.body < j.from.body;
    };
    auto size = interfaces.size();
    std::vector<double> area(2 * size);
    std::vector<esint> faces(2 * size);
    std::vector<esint> triasum(2 * size), triaoffset(2 * size);
    for (auto i = istats.begin(); i != istats.end(); ++i) {
        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            auto it = std::lower_bound(interfaces.begin(), interfaces.end(), Interface(i->first, j->first), comp);
            if (it == interfaces.end() || it->from.body != i->first || it->to.body != j->first) {
                it = std::lower_bound(interfaces.begin(), interfaces.end(), Interface(j->first, i->first), comp);
                auto offset = it - interfaces.begin();
                area[offset + size] += j->second.area;
                faces[offset + size] += j->second.faces;
                triasum[offset + size] += j->second.triangles;
            } else {
                auto offset = it - interfaces.begin();
                area[offset] += j->second.area;
                faces[offset] += j->second.faces;
                triasum[offset] += j->second.triangles;
            }
        }
    }

    Communication::allReduce(area, Communication::OP::SUM);
    Communication::allReduce(faces, Communication::OP::SUM);
    std::vector<esint> triasize = triaoffset = triasum;
    Communication::exscan(triasum, triaoffset);
    contact->interfaces.clear();
    for (size_t i = 0; i < interfaces.size(); ++i) {
        if (faces[i] || faces[i + size]) {
            contact->interfaces.push_back(Interface(interfaces[i].from.body, interfaces[i].to.body));
            contact->interfaces.back().from.area = area[i];
            contact->interfaces.back().from.faces = faces[i];
            contact->interfaces.back().from.triangleOffset = triaoffset[i];
            contact->interfaces.back().from.triangleSize = triasize[i];
            contact->interfaces.back().from.triangleTotalSize = triasum[i];
            contact->interfaces.back().to.area = area[i + size];
            contact->interfaces.back().to.faces = faces[i + size];
            contact->interfaces.back().to.triangleOffset = triaoffset[i + size];
            contact->interfaces.back().to.triangleSize = triasize[i + size];
            contact->interfaces.back().to.triangleTotalSize = triasum[i + size];
            contact->interfaces.back().setOrientation();
        }
    }
    std::sort(contact->interfaces.begin(), contact->interfaces.end(), comp);

    std::vector<std::string> bnames(bodies->totalSize);
    for (size_t r = 1; r < elementsRegions.size(); ++r) {
        for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
            if (bnames[elementsRegions[r]->bodies[b]].size()) {
                bnames[elementsRegions[r]->bodies[b]] += "_";
            }
            bnames[elementsRegions[r]->bodies[b]] += elementsRegions[r]->name;
        }
    }

    esint myrank = surfaces.size() - 1;
    std::vector<int> assigned(contact->interfaces.size()), current;
    std::vector<std::vector<esint> > sBuffer(contact->neighbors.size(), std::vector<esint>(contact->interfaces.size() + 1)), rBuffer(contact->neighbors.size());

    auto preprocess = [&] (size_t index, esint from, esint to) {
        istats[from][to].skip = false;
        for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
            if (s->body != from) {
                continue;
            }
            for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
                if (d->body == to && d->neigh < myrank) {
                    ++sBuffer[d->neigh][index];
                    sBuffer[d->neigh].push_back(surfaces[d->neigh]->fID->datatarray()[d->element]);
                }
            }
        }
    };

    auto create = [&] (const std::string &name, size_t index, esint from, esint to) {
        contactInterfaces.push_back(new ContactInterfaceStore("CONTACT-" + name + "-" + bnames[from] + "-" + bnames[to], index));

        std::vector<esint> dist = { 0 }, data;
        std::vector<Element*> epointer;

        auto push = [&] (esint e) {
            auto face = surfaces.back()->enodes->begin() + e;
            epointer.push_back(surfaces.back()->epointers->datatarray()[e]);
            for (auto n = face->begin(); n != face->end(); ++n) {
                data.push_back(*n);
                contact->nodeMultiplicity->data[*n] += 1;
            }
            dist.push_back(data.size());
        };

        std::vector<esint> dense;
        for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
            if (s->body != from) {
                continue;
            }
            bool addsparse = true;
            for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
                if (d->body == to) {
                    if (addsparse) {
                        addsparse = false;
                        push(s->element);
                    }
                    if (d->neigh == myrank) {
                        dense.push_back(d->element);
                    }
                }
            }
        }
        for (size_t n = 0; n < rBuffer.size(); ++n) {
            for (esint i = rBuffer[n][index]; i < rBuffer[n][index + 1]; ++i) {
                dense.push_back(rBuffer[n][i]);
            }
        }
        utils::sortAndRemoveDuplicates(dense);
        for (size_t i = 0; i < dense.size(); ++i) {
            push(dense[i]);
        }

        contactInterfaces.back()->originalDimension = 2;
        contactInterfaces.back()->dimension = 2;
        contactInterfaces.back()->elements = new serializededata<esint, esint>(dist, data);
        contactInterfaces.back()->epointers = new serializededata<esint, Element*>(1, epointer);
    };

    for (int insert = 0; insert <= 1; ++insert) {
        for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
            switch (it->second.detection) {
            case ContactInterfaceConfiguration::DETECTION::ALL_BODIES:
                for (size_t c = 0; c < contact->interfaces.size(); ++c) {
                    if (assigned[c] == insert) {
                        if (insert) {
                            it->second.found_interfaces.push_back(contactInterfaces.size());
                            create(it->first, c, contact->interfaces[c].from.body, contact->interfaces[c].to.body);
                        } else {
                            preprocess(c, contact->interfaces[c].from.body, contact->interfaces[c].to.body);
                        }
                        assigned[c] = insert + 1;
                    }
                }
                break;
            case ContactInterfaceConfiguration::DETECTION::BODY_LIST:
                current.clear();
                for (size_t b = 0; b < it->second.body_list.size(); ++b) {
                    const ElementsRegionStore *region = eregion(it->second.body_list[b]);
                    current.insert(current.end(), region->bodies.begin(), region->bodies.end());
                }
                utils::sortAndRemoveDuplicates(current);
                for (size_t c = 0; c < contact->interfaces.size(); ++c) {
                    if (
                            assigned[c] == insert &&
                            std::binary_search(current.begin(), current.end(), contact->interfaces[c].from.body) &&
                            std::binary_search(current.begin(), current.end(), contact->interfaces[c].to.body)) {

                        if (insert) {
                            create(it->first, c, contact->interfaces[c].from.body, contact->interfaces[c].to.body);
                        } else {
                            preprocess(c, contact->interfaces[c].from.body, contact->interfaces[c].to.body);
                        }
                        assigned[c] = insert + 1;
                    }
                }
                break;
            case ContactInterfaceConfiguration::DETECTION::CONTACT_PAIR:
                break;
            }
        }

        if (insert == 0) {
            for (size_t n = 0; n < sBuffer.size(); ++n) {
                esint sum = contact->interfaces.size() + 1;
                for (size_t i = 0; i <= contact->interfaces.size(); ++i) {
                    esint tmp = sBuffer[n][i];
                    sBuffer[n][i] = sum;
                    sum += tmp;
                }
            }
            if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, contact->neighbors)) {
                eslog::internalFailure("cannot exchange contact faces to neighbors.\n");
            }
        }
    }

    for (auto s = sside->datatarray().begin(); s != sside->datatarray().end(); ++s) {
        for (auto d = dside->datatarray().begin() + s->denseSegmentBegin; d != dside->datatarray().begin() + s->denseSegmentEnd; ++d) {
            d->skip = d->skip | istats[surfaces.back()->body->datatarray()[s->element]][surfaces[d->neigh]->body->datatarray()[d->element]].skip;
        }
    }

    contact->nodeMultiplicity->synchronize();
    for (size_t i = 0; i < contactInterfaces.size(); ++i) {
        for (size_t e = 0; e < contactInterfaces[i]->elements->structures(); ++e) {
            switch (contactInterfaces[i]->epointers->datatarray()[e]->code) {
            case Element::CODE::LINE2    : normal2D<Element::CODE::LINE2    , 2>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            case Element::CODE::LINE3    : normal2D<Element::CODE::LINE3    , 3>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            case Element::CODE::TRIANGLE3: normal3D<Element::CODE::TRIANGLE3, 3, 6>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            case Element::CODE::TRIANGLE6: normal3D<Element::CODE::TRIANGLE6, 6, 6>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            case Element::CODE::SQUARE4  : normal3D<Element::CODE::SQUARE4  , 4, 1>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            case Element::CODE::SQUARE8  : normal3D<Element::CODE::SQUARE8  , 8, 1>(contact->nodeMultiplicity, contact->nodeNormals, nodes->coordinates, contactInterfaces[i], e, displacement); break;
            default:
                eslog::internalFailure("unknown or not implemented surface element.\n");
            }
        }
    }

    DebugOutput::contact(1, 1);
    profiler::syncend("arrange_contact_interface");
    eslog::checkpointln("MESH: CONTACT INTERFACE ARRANGED");
}

}
}
