
#include "mortar.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/contactstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "math/math.h"
#include "wrappers/mpi/communication.h"

#include <set>
#include <map>
#include <numeric>

namespace espreso {

static void tria3ReferenceCoords(const Matrix_Dense<double> &vertices, const Matrix_Dense<double> &points, Matrix_Dense<double> &result)
{
    // assume non-zero denominator
    result.resize(points.nrows, 2);

    double ux = vertices.vals[1 * vertices.ncols + 0] - vertices.vals[0 * vertices.ncols + 0];
    double uy = vertices.vals[1 * vertices.ncols + 1] - vertices.vals[0 * vertices.ncols + 1];
    double vx = vertices.vals[2 * vertices.ncols + 0] - vertices.vals[0 * vertices.ncols + 0];
    double vy = vertices.vals[2 * vertices.ncols + 1] - vertices.vals[0 * vertices.ncols + 1];
    double uv = ux * vx + uy * vy, uu = ux * ux + uy * uy, vv = vx * vx + vy * vy;
    double denominator = uv * uv - uu * vv;
    for (esint r = 0; r < points.nrows; ++r) {
        double wx = points.vals[r * points.ncols + 0] - vertices.vals[0 * vertices.ncols + 0];
        double wy = points.vals[r * points.ncols + 1] - vertices.vals[0 * vertices.ncols + 1];
        double wu = wx * ux + wy * uy, wv = wx * vx + wy * vy;
        result.vals[r * result.ncols + 0] = (uv * wv - vv * wu) / denominator;
        result.vals[r * result.ncols + 1] = (uv * wu - uu * wv) / denominator;
    }
}

static void square4ReferenceCoords(const Matrix_Dense<double> &vertices, const Matrix_Dense<double> &points, Matrix_Dense<double> &result)
{
    Matrix_Dense<double> tmp; tmp.resize(4, 2);
    tmp.vals[0 * tmp.ncols + 0] = -vertices.vals[0 * vertices.ncols + 0] + vertices.vals[1 * vertices.ncols + 0] + vertices.vals[2 * vertices.ncols + 0] - vertices.vals[3 * vertices.ncols + 0];
    tmp.vals[1 * tmp.ncols + 0] = -vertices.vals[0 * vertices.ncols + 0] - vertices.vals[1 * vertices.ncols + 0] + vertices.vals[2 * vertices.ncols + 0] + vertices.vals[3 * vertices.ncols + 0];
    tmp.vals[2 * tmp.ncols + 0] =  vertices.vals[0 * vertices.ncols + 0] - vertices.vals[1 * vertices.ncols + 0] + vertices.vals[2 * vertices.ncols + 0] - vertices.vals[3 * vertices.ncols + 0];
    tmp.vals[3 * tmp.ncols + 0] =  vertices.vals[0 * vertices.ncols + 0] + vertices.vals[1 * vertices.ncols + 0] + vertices.vals[2 * vertices.ncols + 0] + vertices.vals[3 * vertices.ncols + 0];
    tmp.vals[0 * tmp.ncols + 1] = -vertices.vals[0 * vertices.ncols + 1] + vertices.vals[1 * vertices.ncols + 1] + vertices.vals[2 * vertices.ncols + 1] - vertices.vals[3 * vertices.ncols + 1]; // s
    tmp.vals[1 * tmp.ncols + 1] = -vertices.vals[0 * vertices.ncols + 1] - vertices.vals[1 * vertices.ncols + 1] + vertices.vals[2 * vertices.ncols + 1] + vertices.vals[3 * vertices.ncols + 1]; // t
    tmp.vals[2 * tmp.ncols + 1] =  vertices.vals[0 * vertices.ncols + 1] - vertices.vals[1 * vertices.ncols + 1] + vertices.vals[2 * vertices.ncols + 1] - vertices.vals[3 * vertices.ncols + 1]; // st
    tmp.vals[3 * tmp.ncols + 1] =  vertices.vals[0 * vertices.ncols + 1] + vertices.vals[1 * vertices.ncols + 1] + vertices.vals[2 * vertices.ncols + 1] + vertices.vals[3 * vertices.ncols + 1]; // 1

    Matrix_Dense<double>  F;  F.resize(1, 2);
    Matrix_Dense<double> dF; dF.resize(2, 2);
    double bx, by, rx, ry, detdF, normF;
    double detdF0 = tmp.vals[0 * tmp.ncols + 0] * tmp.vals[1 * tmp.ncols + 1] - tmp.vals[0 * tmp.ncols + 1] * tmp.vals[1 * tmp.ncols + 0];
    int cnt;
    result.resize(points.nrows, 2);
    for (esint i = 0; i < points.nrows; ++i) {
        bx = 4.0 * points.vals[i * points.ncols + 0] - tmp.vals[3 * tmp.ncols + 0];
        by = 4.0 * points.vals[i * points.ncols + 1] - tmp.vals[3 * tmp.ncols + 1];
        rx = ( tmp.vals[1 * tmp.ncols + 1] * bx - tmp.vals[1 * tmp.ncols + 0] * by) / detdF0; // s
        ry = (-tmp.vals[0 * tmp.ncols + 1] * bx + tmp.vals[0 * tmp.ncols + 0] * by) / detdF0; // t
        cnt = 0;
        F .vals[0 *  F.ncols + 0] = tmp.vals[0 * tmp.ncols + 0] * rx + tmp.vals[1 * tmp.ncols + 0] * ry + tmp.vals[2 * tmp.ncols + 0] * rx * ry - bx;
        F .vals[0 *  F.ncols + 1] = tmp.vals[0 * tmp.ncols + 1] * rx + tmp.vals[1 * tmp.ncols + 1] * ry + tmp.vals[2 * tmp.ncols + 1] * rx * ry - by;
        dF.vals[0 * dF.ncols + 0] = tmp.vals[0 * tmp.ncols + 0] + tmp.vals[2 * tmp.ncols + 0] * ry;
        dF.vals[0 * dF.ncols + 1] = tmp.vals[0 * tmp.ncols + 1] + tmp.vals[2 * tmp.ncols + 1] * ry;
        dF.vals[1 * dF.ncols + 0] = tmp.vals[1 * tmp.ncols + 0] + tmp.vals[2 * tmp.ncols + 0] * rx;
        dF.vals[1 * dF.ncols + 1] = tmp.vals[1 * tmp.ncols + 1] + tmp.vals[2 * tmp.ncols + 1] * rx;
        detdF = tmp.vals[0 * tmp.ncols + 0] * tmp.vals[1 * tmp.ncols + 1] - tmp.vals[0 * tmp.ncols + 1] * tmp.vals[1 * tmp.ncols + 0];
        normF = std::sqrt(F.vals[0 * F.ncols + 0] * F.vals[0 * F.ncols + 0] + F.vals[0 * F.ncols + 1] * F.vals[0 * F.ncols + 1]);
        while (cnt > 0 && normF > 1e-5 && std::abs(detdF / detdF0) > 1e-5) {
            rx -= ( dF.vals[1 * dF.ncols + 1] * F.vals[0 * F.ncols + 0] - dF.vals[1 * dF.ncols + 0] * F.vals[0 * F.ncols + 1]) / detdF;
            ry -= (-dF.vals[1 * dF.ncols + 0] * F.vals[0 * F.ncols + 0] + dF.vals[0 * dF.ncols + 0] * F.vals[0 * F.ncols + 1]) / detdF;
            F .vals[0 *  F.ncols + 0] = tmp.vals[0 * tmp.ncols + 0] * rx + tmp.vals[1 * tmp.ncols + 0] * ry + tmp.vals[2 * tmp.ncols + 0] * rx * ry - bx;
            F .vals[0 *  F.ncols + 1] = tmp.vals[0 * tmp.ncols + 1] * rx + tmp.vals[1 * tmp.ncols + 1] * ry + tmp.vals[2 * tmp.ncols + 1] * rx * ry - by;
            dF.vals[0 * dF.ncols + 0] = tmp.vals[0 * tmp.ncols + 0] + tmp.vals[2 * tmp.ncols + 0] * ry;
            dF.vals[0 * dF.ncols + 1] = tmp.vals[0 * tmp.ncols + 1] + tmp.vals[2 * tmp.ncols + 1] * ry;
            dF.vals[1 * dF.ncols + 0] = tmp.vals[1 * tmp.ncols + 0] + tmp.vals[2 * tmp.ncols + 0] * rx;
            dF.vals[1 * dF.ncols + 1] = tmp.vals[1 * tmp.ncols + 1] + tmp.vals[2 * tmp.ncols + 1] * rx;
            detdF = tmp.vals[0 * tmp.ncols + 0] * tmp.vals[1 * tmp.ncols + 1] - tmp.vals[0 * tmp.ncols + 1] * tmp.vals[1 * tmp.ncols + 0];
            normF = std::sqrt(F.vals[0 * F.ncols + 0] * F.vals[0 * F.ncols + 0] + F.vals[0 * F.ncols + 1] * F.vals[0 * F.ncols + 1]);
            cnt--;
        }
        result.vals[i * result.ncols + 0] = rx;
        result.vals[i * result.ncols + 1] = ry;
    }
}

static void recomputeDetJ(Matrix_Dense<double> &coords, Matrix_Dense<double> &resdetJ)
{
    resdetJ.resize(1, 6);
    for (int gp = 0; gp < 6; ++gp) {
        double j00 = -coords.vals[0 * coords.ncols + 0] + coords.vals[1 * coords.ncols + 0];
        double j01 = -coords.vals[0 * coords.ncols + 1] + coords.vals[1 * coords.ncols + 1];
        double j10 = -coords.vals[0 * coords.ncols + 0] + coords.vals[2 * coords.ncols + 0];
        double j11 = -coords.vals[0 * coords.ncols + 1] + coords.vals[2 * coords.ncols + 1];
        resdetJ.vals[0 * resdetJ.ncols + gp] = j00 * j11 - j01 * j10;
    }
}

static void recomputeDetJN(Element *e, Matrix_Dense<double> &coords, Matrix_Dense<double> &resdetJ, Matrix_Dense<double> &resN, Matrix_Dense<double> &refPoints)
{
    resdetJ.resize(1, refPoints.nrows);
    resN.resize(e->nodes, refPoints.nrows);
    for (esint p = 0; p < refPoints.nrows; p++) {
        switch (e->code) {
        case Element::CODE::TRIANGLE3: {
            resN.vals[0 * resN.ncols + p] =  1 - refPoints.vals[p * refPoints.ncols + 0] - refPoints.vals[p * refPoints.ncols + 1];
            resN.vals[1 * resN.ncols + p] =      refPoints.vals[p * refPoints.ncols + 0];
            resN.vals[2 * resN.ncols + p] =                                                refPoints.vals[p * refPoints.ncols + 1];

            double j00 = -coords.vals[0 * coords.ncols + 0] + coords.vals[1 * coords.ncols + 0];
            double j01 = -coords.vals[0 * coords.ncols + 1] + coords.vals[1 * coords.ncols + 1];
            double j10 = -coords.vals[0 * coords.ncols + 0] + coords.vals[2 * coords.ncols + 0];
            double j11 = -coords.vals[0 * coords.ncols + 1] + coords.vals[2 * coords.ncols + 1];
            resdetJ.vals[0 * resdetJ.ncols + p] = j00 * j11 - j01 * j10;
        } break;
        case Element::CODE::SQUARE4: {
            resN.vals[0 * resN.ncols + p] = 0.25 * (1 - refPoints.vals[p * refPoints.ncols + 0]) * (1 - refPoints.vals[p * refPoints.ncols + 1]);
            resN.vals[1 * resN.ncols + p] = 0.25 * (1 + refPoints.vals[p * refPoints.ncols + 0]) * (1 - refPoints.vals[p * refPoints.ncols + 1]);
            resN.vals[2 * resN.ncols + p] = 0.25 * (1 + refPoints.vals[p * refPoints.ncols + 0]) * (1 + refPoints.vals[p * refPoints.ncols + 1]);
            resN.vals[3 * resN.ncols + p] = 0.25 * (1 - refPoints.vals[p * refPoints.ncols + 0]) * (1 + refPoints.vals[p * refPoints.ncols + 1]);
            double dN00 = 0.25 * ( refPoints.vals[p * refPoints.ncols + 1] - 1); double dN10 = 0.25 * ( refPoints.vals[p * refPoints.ncols + 0] - 1);
            double dN01 = 0.25 * (-refPoints.vals[p * refPoints.ncols + 1] + 1); double dN11 = 0.25 * (-refPoints.vals[p * refPoints.ncols + 0] - 1);
            double dN02 = 0.25 * ( refPoints.vals[p * refPoints.ncols + 1] + 1); double dN12 = 0.25 * ( refPoints.vals[p * refPoints.ncols + 0] + 1);
            double dN03 = 0.25 * (-refPoints.vals[p * refPoints.ncols + 1] - 1); double dN13 = 0.25 * (-refPoints.vals[p * refPoints.ncols + 0] + 1);

            double j00 = dN00 * coords.vals[0 * coords.ncols + 0] + dN01 * coords.vals[1 * coords.ncols + 0] + dN02 * coords.vals[2 * coords.ncols + 0] + dN03 * coords.vals[3 * coords.ncols + 0];
            double j01 = dN00 * coords.vals[0 * coords.ncols + 1] + dN01 * coords.vals[1 * coords.ncols + 1] + dN02 * coords.vals[2 * coords.ncols + 1] + dN03 * coords.vals[3 * coords.ncols + 1];
            double j10 = dN10 * coords.vals[0 * coords.ncols + 0] + dN11 * coords.vals[1 * coords.ncols + 0] + dN12 * coords.vals[2 * coords.ncols + 0] + dN13 * coords.vals[3 * coords.ncols + 0];
            double j11 = dN10 * coords.vals[0 * coords.ncols + 1] + dN11 * coords.vals[1 * coords.ncols + 1] + dN12 * coords.vals[2 * coords.ncols + 1] + dN13 * coords.vals[3 * coords.ncols + 1];
            resdetJ.vals[0 * resdetJ.ncols + p] = j00 * j11 - j01 * j10;
        } break;
        default:
            break;
        }
    }
}

#define TOLERATED_SLAVE_COVER_RATIO_FOR_DUAL_SHAPE_COEFFICIENTS_ON_WHOLE_ELEMENT 10 // never at this versions
#define BE_VALUE_TRESHOLD 1e-15

template <typename T>
void MortarContact<T>::assembleMortarInterface(std::vector<Mortar> &B)
{
    int printrank = -1;

    int tGPs = 6;
    std::vector<double> w = { 0.111690794839005, 0.111690794839005, 0.111690794839005, 0.054975871827661, 0.054975871827661, 0.054975871827661};
    std::vector<double> s = { 0.445948490915965, 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.091576213509771, 0.816847572980459};
    std::vector<double> t = { 0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771};

    Matrix_Dense<double> sCoords, dCoords, tCoords, sRefCoords, dRefCoords;
    sRefCoords.resize(3, 2); dRefCoords.resize(3, 2);
    Matrix_Dense<double> sGpCoords, dGpCoords;
    sGpCoords.resize(tGPs, 2), dGpCoords.resize(tGPs, 2);
    Matrix_Dense<double> tmp, sDetJ, dDetJ, tDetJ, sN, dN;

    const std::vector<SurfaceStore*> &surfaces = info::mesh->contact->surfaces;

    auto getReferenceCoords = [&] (Element *e, Matrix_Dense<double> &eCoords, Matrix_Dense<double> &tCoords, Matrix_Dense<double> &ref) {
        switch (e->code) { // TODO: improve for other types
        case Element::CODE::TRIANGLE3: tria3ReferenceCoords(eCoords, tCoords, ref); break;
        case Element::CODE::SQUARE4: square4ReferenceCoords(eCoords, tCoords, ref); break;
        default: eslog::internalFailure("not implemented mortar element.\n");
        }
    };

    auto getGpCoords = [&] (Matrix_Dense<double> &refCoords, Matrix_Dense<double> &gpCoords) {
        for (int i = 0; i < tGPs; ++i) {
            gpCoords.vals[i * gpCoords.ncols + 0] = (1 - s[i] - t[i]) * refCoords.vals[0 * refCoords.ncols + 0] + s[i] * refCoords.vals[1 * refCoords.ncols + 0] + t[i] * refCoords.vals[2 * refCoords.ncols + 0];
            gpCoords.vals[i * gpCoords.ncols + 1] = (1 - s[i] - t[i]) * refCoords.vals[0 * refCoords.ncols + 1] + s[i] * refCoords.vals[1 * refCoords.ncols + 1] + t[i] * refCoords.vals[2 * refCoords.ncols + 1];
        }
    };

    std::unordered_map<esint, std::unordered_map<esint, esint> > interface;
    for (size_t i = 0; i < info::mesh->contact->interfaces.size(); ++i) {
        interface[info::mesh->contact->interfaces[i].from.body][info::mesh->contact->interfaces[i].to.body] = i;
    }

    auto *sside = info::mesh->contact->sparseSide;
    auto *dside = info::mesh->contact->denseSide;
    double *coors = info::mesh->contact->planeCoordinates->datatarray().data()->data();
    for (auto ss = sside->datatarray().begin(); ss != sside->datatarray().end(); ++ss) {
        for (auto iface = interface[ss->body].begin(); iface != interface[ss->body].end(); ++iface) {
            Element* sElement = surfaces.back()->epointers->datatarray()[ss->element];
            auto sNodes = surfaces.back()->enodes->begin() + ss->element;
            const auto &sIDs = surfaces.back()->nIDs->datatarray();
            sCoords.nrows = sElement->nodes;
            sCoords.ncols = 2;
            sCoords.vals = coors + 2 * ss->coordinateOffset;
            if (printrank == info::mpi::rank) {
                printf("s-%d-coords: ", ss->element);
                for (int n = 0; n < sElement->nodes; ++n) {
                    printf("<%+.4f, %+.4f>", sCoords.vals[n * sCoords.ncols + 0], sCoords.vals[n * sCoords.ncols + 1]);
                }
                printf("\n");
            }

            std::vector<double> sparse(sElement->nodes * sElement->nodes);
            Matrix_Dense<double> D, M;
            D.resize(sElement->nodes, sElement->nodes); M.resize(sElement->nodes, sElement->nodes);
            math::set(D, T{0}); math::set(M, T{0});
            for (int dmReady = 0; dmReady <= 1; ++dmReady) {
                for (auto ds = dside->datatarray().begin() + ss->denseSegmentBegin; ds != dside->datatarray().begin() + ss->denseSegmentEnd; ++ds) {
                    if (ds->skip) {
                        continue;
                    }
                    if (iface->first != ds->body) {
                        continue;
                    }

                    Element* dElement = surfaces[ds->neigh]->epointers->datatarray()[ds->element];
                    auto dNodes = surfaces[ds->neigh]->enodes->begin() + ds->element;
                    const auto &dIDs = surfaces[ds->neigh]->nIDs->datatarray();
                    if (dmReady) {
                        dCoords.nrows = sElement->nodes;
                        dCoords.ncols = 2;
                        dCoords.vals = coors + 2 * ds->coordinateOffset;
                        if (printrank == info::mpi::rank) {
                            printf("d-%d-coords: ", ds->element);
                            for (int n = 0; n < dElement->nodes; ++n) {
                                printf("<%+.4f, %+.4f>", dCoords.vals[n * dCoords.ncols + 0], dCoords.vals[n * dCoords.ncols + 1]);
                            }
                            printf("\n");
                        }
                    }

                    if (printrank == info::mpi::rank) printf("assemble %d-%d\n", ss->element, ds->element);

                    std::vector<double> dense(sElement->nodes * dElement->nodes);
                    for (esint dt = 0; dt < ds->triangles; ++dt) {
                        tCoords.nrows = 3;
                        tCoords.ncols = 2;
                        tCoords.vals = coors + 2 * (ds->triangleOffset + 3 * dt);
                        getReferenceCoords(sElement, sCoords, tCoords, sRefCoords);
                        getGpCoords(sRefCoords, sGpCoords);

                        recomputeDetJ(sRefCoords, tDetJ);
                        recomputeDetJN(sElement, sCoords, sDetJ, sN, sGpCoords);

                        if (!dmReady) {
                            // compute De, Me for coefficients Ae of dual basis functions from formula (4.60)
                            for (int gp = 0; gp < tGPs; ++gp) {
                                double weight = w[gp] * tDetJ.vals[0 * tDetJ.ncols + gp] * sDetJ.vals[0 * tDetJ.ncols + gp];
                                for (int i = 0; i < sElement->nodes; i++) {
                                    D.vals[i * sElement->nodes + i] += weight * sN.vals[i * sN.ncols + gp];
                                    M.vals[i * sElement->nodes + i] += weight * sN.vals[i * sN.ncols + gp] * sN.vals[i * sN.ncols + gp];
                                    for (int j = i + 1; j < sElement->nodes; j++) {
                                        M.vals[j * sElement->nodes + i] += weight * sN.vals[i * sN.ncols + gp] * sN.vals[j * sN.ncols + gp];
                                        M.vals[i * sElement->nodes + j] += weight * sN.vals[i * sN.ncols + gp] * sN.vals[j * sN.ncols + gp];
                                    }
                                }
                            }
                        } else {
                            getReferenceCoords(dElement, dCoords, tCoords, dRefCoords);
                            getGpCoords(dRefCoords, dGpCoords);
                            recomputeDetJN(dElement, dCoords, dDetJ, dN, dGpCoords);

                            if (printrank == info::mpi::rank) {
                                printf("\tt-%d-%d-coords: <%+.4f, %+.4f><%+.4f, %+.4f><%+.4f, %+.4f>\n", ss->element, ds->element, tCoords.vals[0 * tCoords.ncols + 0], tCoords.vals[0 * tCoords.ncols + 1], tCoords.vals[1 * tCoords.ncols + 0], tCoords.vals[1 * tCoords.ncols + 1], tCoords.vals[2 * tCoords.ncols + 0], tCoords.vals[2 * tCoords.ncols + 1]);
                                printf("\tt-%d-%d-s-ref-coords: <%+.4f, %+.4f><%+.4f, %+.4f><%+.4f, %+.4f>\n", ss->element, ds->element, sRefCoords.vals[0 * sRefCoords.ncols + 0], sRefCoords.vals[0 * sRefCoords.ncols + 1], sRefCoords.vals[1 * sRefCoords.ncols + 0], sRefCoords.vals[1 * sRefCoords.ncols + 1], sRefCoords.vals[2 * sRefCoords.ncols + 0], sRefCoords.vals[2 * sRefCoords.ncols + 1]);
                                printf("\tt-%d-%d-d-ref-coords: <%+.4f, %+.4f><%+.4f, %+.4f><%+.4f, %+.4f>\n", ss->element, ds->element, dRefCoords.vals[0 * dRefCoords.ncols + 0], dRefCoords.vals[0 * dRefCoords.ncols + 1], dRefCoords.vals[1 * dRefCoords.ncols + 0], dRefCoords.vals[1 * dRefCoords.ncols + 1], dRefCoords.vals[2 * dRefCoords.ncols + 0], dRefCoords.vals[2 * dRefCoords.ncols + 1]);
                            }

                            Matrix_Dense<double> psi;
                            psi.resize(D.nrows, sN.ncols);
                            math::blas::multiply(T{1}, D, sN, T{0}, psi, true, false);
                            if (printrank == info::mpi::rank) {
                                for (int r = 0; r < psi.nrows; ++r) {
                                    for (int c = 0; c < psi.ncols; ++c) {
                                        printf(" %e", psi.vals[r * psi.ncols + c]);
                                    }
                                    printf("\n");
                                }
                                printf("\n");
                            }

                            for (int gp = 0; gp < tGPs; ++gp) {
                                double weight = w[gp] * tDetJ.vals[0 * tDetJ.ncols + gp] * sDetJ.vals[0 * tDetJ.ncols + gp];
                                for (int i = 0; i < sElement->nodes; i++) {
                                    for (int j = 0; j < sElement->nodes; j++) {
                                        sparse[i * sElement->nodes + j] += weight * sN.vals[j * sN.ncols + gp] * psi.vals[i * psi.ncols + gp];
                                    }
                                    for (int j = 0; j < dElement->nodes; j++) {
                                        dense[i * dElement->nodes + j] -= weight * dN.vals[j * dN.ncols + gp] * psi.vals[i * psi.ncols + gp];
                                    }
                                }
                            }
                        }
                    }
                    if (dmReady) {
                        for (int i = 0; i < sElement->nodes; i++) {
                            for (int j = 0; j < dElement->nodes; j++) {
                                if (std::fabs(dense[i * dElement->nodes + j]) > BE_VALUE_TRESHOLD) {
                                    if (printrank == info::mpi::rank) printf("  <%d,%d>=%f\n", sIDs[sNodes->at(i)], dIDs[dNodes->at(j)], dense[i * dElement->nodes + j]);
                                    B.push_back(Mortar(iface->second, sIDs[sNodes->at(i)], dIDs[dNodes->at(j)], dense[i * dElement->nodes + j]));
                                }
                            }
                        }
                    }
                }
                if (!dmReady) {
                    if (printrank == info::mpi::rank) {
                        printf("M");
                        for (int i = 0; i < M.nnz; ++i) {
                            printf(" %.7f", M.vals[i]);
                        }
                        printf("\n");
                        printf("D");
                        for (int i = 0; i < D.nnz; ++i) {
                            printf(" %.2f", D.vals[i]);
                        }
                        printf("\n");
                    }
                    // inv(M) * D  ... need to be transposed to evaluate psi coefs !!!!
                    math::lapack::solve_general(M, D);
                    if (printrank == info::mpi::rank) {
                        printf("D");
                        for (int i = 0; i < D.nnz; ++i) {
                            printf(" %.1f", D.vals[i]);
                        }
                        printf("\n");
                    }
                }
            }
            for (int i = 0; i < sElement->nodes; i++) {
                for (int j = 0; j < sElement->nodes; j++) {
                    if (std::fabs(sparse[i * sElement->nodes + j]) > BE_VALUE_TRESHOLD) {
                        if (printrank == info::mpi::rank) printf("<%d,%d>=%f\n", sIDs[sNodes->at(i)], sIDs[sNodes->at(j)], sparse[i * sElement->nodes + j]);
                        B.push_back(Mortar(iface->second, sIDs[sNodes->at(i)], sIDs[sNodes->at(j)], sparse[i * sElement->nodes + j]));
                    }
                }
            }
        }
    }
}

template <typename T>
void MortarContact<T>::synchronize(FETI<T> &feti, std::vector<Mortar> &B)
{
    std::vector<esint> idPerm(info::mesh->nodes->size);
    std::iota(idPerm.begin(), idPerm.end(), 0);
    std::sort(idPerm.begin(), idPerm.end(), [&] (esint i, esint j) { return info::mesh->nodes->IDs->datatarray()[i] < info::mesh->nodes->IDs->datatarray()[j]; });

//    Communication::serialize([&] () {
//        printf("%2d :: ", info::mpi::rank);
//        for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
//            printf(" %d", info::mesh->neighbors[n]);
//        }
//        printf("\n");
//    });

//    Communication::serialize([&] () {
//        printf("%2d ::\n", info::mpi::rank);
//        for (size_t i = 0; i < B.size(); ++i) {
//            printf("    :: %3d : %3d == %+.2e\n", B[i].from, B[i].to, B[i].value);
//        }
//    });

    { // get node-to-domains map and duplicate gluing
        std::map<esint, esint*> n2d;
        for (size_t i = 0; i < B.size(); i++) {
            n2d.insert(std::make_pair(B[i].to, nullptr));
        }
        std::vector<esint> nodes;
        for (auto it = n2d.cbegin(); it != n2d.cend(); ++it) {
            nodes.push_back(it->first);
        }

        std::vector<std::vector<esint> > sBuffer(info::mesh->neighborsWithMe.size(), nodes), rBuffer(info::mesh->neighborsWithMe.size());
        if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighborsWithMe)) {
            eslog::internalFailure("cannot exchange mortar nodes.\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            sBuffer[n].clear();
            auto &ids = info::mesh->nodes->IDs->datatarray();
            for (size_t i = 0, j = 0; i < rBuffer[n].size(); ++i) {
                while (j < idPerm.size() && ids[idPerm[j]] < rBuffer[n][i]) { ++j; }
                if (j < idPerm.size() && ids[idPerm[j]] == rBuffer[n][i]) {
                    auto dmap = info::mesh->nodes->domains->begin() + idPerm[j];
                    sBuffer[n].push_back(rBuffer[n][i]);
                    sBuffer[n].push_back(dmap->size());
                    for (auto d = dmap->begin(); d != dmap->end(); ++d) {
                        sBuffer[n].push_back(*d);
                    }
                }
            }
        }

        if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighborsWithMe)) {
            eslog::internalFailure("cannot exchange mortar nodes-to-domain map.\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            for (size_t i = 0; i < rBuffer[n].size(); ) {
                esint id = rBuffer[n][i++];
                esint domains = rBuffer[n][i];
                n2d[id] = rBuffer[n].data() + i;
                i += domains + 1;
            }
        }

        // duplicate data to domains (still per element data)
        std::vector<Mortar> pure; B.swap(pure);
        for (size_t i = 0; i < pure.size(); ++i) {
            esint *domains = n2d[pure[i].to];
            for (int d = 0; d < domains[0]; ++d) {
                B.push_back(pure[i]);
                B.back().value /= domains[0];
                B.back().domain = domains[1 + d];
            }
        }
    }

    { // exchange with common 'sparse' nodes
        std::set<esint> local;
        for (size_t i = 0; i < B.size(); ++i) {
            local.insert(B[i].from);
        }

        std::vector<std::vector<Mortar> > rBuffer(info::mesh->neighbors.size());
        if (!Communication::exchangeUnknownSize(B, rBuffer, info::mesh->neighbors)) {
            eslog::internalFailure("cannot synchronize mortars.\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            for (size_t i = 0; i < rBuffer[n].size(); ++i) {
                if (local.count(rBuffer[n][i].from)) {
                    B.push_back(rBuffer[n][i]);
                }
            }
        }
    } // now all processes that generated mortar have complete mortars


    // merge values
    std::sort(B.begin(), B.end());
    if (B.size()) {
        size_t unique = 0;
        for (size_t i = 1; i < B.size(); i++) {
            if (B[unique] != B[i]) {
                B[++unique] = B[i];
            } else {
                B[unique].value += B[i].value;
            }
        }
        B.resize(unique + 1);
    }

//    Communication::serialize([&] () {
//        printf("%2d ::\n", info::mpi::rank);
//        for (size_t i = 0, j, n = 0; i < B.size(); ++n) {
//            printf("    :: %3d ", B[i].from);
//            for (j = i; j < B.size() && B[i].from == B[j].from; ++j) {
//                if (j && B[j].to == B[j - 1].to) {
//                    printf("[%2d,%+.2e]", B[j].domain, B[j].value);
//                } else {
//                    printf(" %3d[%2d,%+.2e]", B[j].to, B[j].domain, B[j].value);
//                }
//            }
//            printf("\n");
//            i = j;
//        }
//    });

    // normalize values
    for (size_t i = 0, j = 0; i < B.size(); i = j) {
        while (j < B.size() && B[i].from == B[j].from) { ++j; }
        double sum = 0;
        for (size_t k = i; k < j; ++k) {
            sum += B[k].value * B[k].value;
        }
        double norm = 1 / std::sqrt(sum);
        for (size_t k = i; k < j; ++k) {
            B[k].value *= norm;
        }
    }

//    std::map<int, double> check, norm;
//    for (size_t i = 0; i < B.size(); ++i) {
//        check[B[i].from] += B[i].value;
//        norm[B[i].from] += B[i].value * B[i].value;
//    }
//    Communication::serialize([&] () {
//        for (auto x = check.begin(); x != check.end(); ++x) {
//            printf("%3d = %+.5e |%+.5e|\n", x->first, x->second, norm[x->first]);
//        }
//    });

    { // exchange with neighbors-of-neighbors and pick only local part of the mortar
        std::vector<int> neighbors = info::mesh->neighbors;
        std::vector<std::vector<int> > rn(neighbors.size());
        if (!Communication::exchangeUnknownSize(neighbors, rn, info::mesh->neighbors)) {
            eslog::internalFailure("cannot exchange neighbors.\n");
        }
        for (size_t n = 0; n < rn.size(); ++n) {
            neighbors.insert(neighbors.end(), rn[n].begin(), rn[n].end());
        }
        utils::sortAndRemoveDuplicates(neighbors);

        std::set<esint> inserted;

        std::vector<std::vector<Mortar> > rBuffer(neighbors.size());
        if (!Communication::exchangeUnknownSize(B, rBuffer, neighbors)) {
            eslog::internalFailure("cannot synchronize mortars.\n");
        }

        for (size_t n = 0; n < rBuffer.size(); ++n) {
            auto &ids = info::mesh->nodes->IDs->datatarray();
            for (size_t i = 0, j = 0; i < rBuffer[n].size(); i = j) {
                while (j < rBuffer[n].size() && rBuffer[n][i].from == rBuffer[n][j].from) { ++j; }
                if (inserted.count(rBuffer[n][i].from) == 0) { // new mortar (it is possible to receive the mortar with more processes)
                    bool local = false;
                    for (size_t k = i; !local && k < j; ++k) {
                        auto it = std::lower_bound(idPerm.begin(), idPerm.end(), rBuffer[n][k].to, [&] (esint id, esint val) { return ids[id] < val; });
                        if (it != idPerm.end() && ids[*it] == rBuffer[n][k].to) {
                            local = true;
                        }
                    }
                    if (local) {
                        mInfo[rBuffer[n][i].from].begin = mortar.size();
                        for (size_t k = i; k < j; ++k) {
                            mInfo[rBuffer[n][i].from].domains.push_back(rBuffer[n][k].domain);
                            auto it = std::lower_bound(idPerm.begin(), idPerm.end(), rBuffer[n][k].to, [&] (esint id, esint val) { return ids[id] < val; });
                            if (it != idPerm.end() && ids[*it] == rBuffer[n][k].to && feti.decomposition->ismy(rBuffer[n][k].domain)) {
                                mortar.push_back(rBuffer[n][k]);
                                mortar.back().to = *it; // global ID to local offset
                            }
                        }
                        utils::sortAndRemoveDuplicates(mInfo[rBuffer[n][i].from].domains);
                        mInfo[rBuffer[n][i].from].end = mortar.size();
                    }
                    inserted.insert(rBuffer[n][i].from);
                }
            }
        }
    }

//    Communication::serialize([&] () {
//        printf("%2d ::\n", info::mpi::rank);
//        for (size_t i = 0, j, n = 0; i < mortar.size(); ++n) {
//            printf("    :: %3d ", mortar[i].from);
//            for (j = i; j < mortar.size() && mortar[i].from == mortar[j].from; ++j) {
//                if (j && mortar[j].to == mortar[j - 1].to) {
//                    printf("[%2d,%+.2e]", mortar[j].domain, mortar[j].value);
//                } else {
//                    printf(" %3d[%2d,%+.2e]", mortar[j].to, mortar[j].domain, mortar[j].value);
//                }
//            }
//            printf("\n");
//            i = j;
//        }
//    });

    { // update neighbors since there can be new neighbors connected via the sparse side
        std::set<esint> domains;
        for (auto info = mInfo.cbegin(); info != mInfo.cend(); ++info) {
            domains.insert(info->second.domains.begin(), info->second.domains.end());
        }

        auto ddist = info::mesh->domains->gatherProcDistribution(); // remove this
        int ni = 0;
        for (auto n = domains.begin(); n != domains.end(); ++n) {
            while (ddist[ni] <= *n) ++ni;
            if (info::mesh->neighborsWithMe.back() != ni - 1) {
                info::mesh->neighborsWithMe.push_back(ni - 1);
            }
        }
        utils::sortAndRemoveDuplicates(info::mesh->neighborsWithMe);
        info::mesh->neighbors.clear();
        for (size_t i = 0; i < info::mesh->neighborsWithMe.size(); ++i) {
            if (info::mesh->neighborsWithMe[i] != info::mpi::rank) {
                info::mesh->neighbors.push_back(info::mesh->neighborsWithMe[i]);
            }
        }

//        Communication::serialize([&] () {
//            printf("%2d :: ", info::mpi::rank);
//            for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
//                printf(" %d", info::mesh->neighbors[n]);
//            }
//            printf("\n");
//        });

        feti.decomposition->update(info::mesh->neighbors);
    }

    { // compute mortar id
        int id = 0;
        std::vector<esint> idmap;
        for (auto m = mInfo.cbegin(); m != mInfo.cend(); ++m) {
            if (feti.decomposition->ismy(m->second.domains.front())) {
                idmap.push_back(mortar[m->second.begin].from);
                idmap.push_back(id++);
            }
        }

        esint moffset = id;
        Communication::exscan(moffset);;
        for (size_t i = 0; i < idmap.size(); i += 2) {
            idmap[i + 1] = mInfo[idmap[i]].id = moffset++;
        }

        std::vector<std::vector<esint> > rids(info::mesh->neighbors.size() + 1);

        if (!Communication::exchangeUnknownSize(idmap, rids, info::mesh->neighbors)) {
            eslog::internalFailure("cannot synchronize mortar ids.\n");
        }
        for (size_t n = 0; n < rids.size(); ++n) {
            for (size_t i = 0; i < rids[n].size(); i += 2) {
                if (mInfo.count(rids[n][i]) == 1) {
                    mInfo[rids[n][i]].id = rids[n][i + 1];
                }
            }
        }
    }


    { // exchange normals
        // compute info
        struct __mInfo__ {
            int from; Point normal;
            __mInfo__(): from(-1) {}
            __mInfo__(int from): from(from) {}
        };

        std::vector<__mInfo__> info;
        auto &ids = info::mesh->nodes->IDs->datatarray();
        for (auto m = mInfo.cbegin(); m != mInfo.cend(); ++m) {
            auto it = std::lower_bound(ids.begin() + info::mesh->nodes->uniqInfo.nhalo, ids.end(), m->first);
            if (it != ids.end() && *it == m->first) {
                info.push_back(__mInfo__(m->first));
                if (StructuralMechanics::Results::normal) {
                    info.back().normal.x = StructuralMechanics::Results::normal->data[info::mesh->dimension * (it - ids.begin()) + 0];
                    info.back().normal.y = StructuralMechanics::Results::normal->data[info::mesh->dimension * (it - ids.begin()) + 1];
                    if (info::mesh->dimension == 3) {
                        info.back().normal.z = StructuralMechanics::Results::normal->data[info::mesh->dimension * (it - ids.begin()) + 2];
                    }
                }
            }
        }

        // exchange info
        std::vector<std::vector<__mInfo__> > rInfo(info::mesh->neighbors.size() + 1);

        if (!Communication::exchangeUnknownSize(info, rInfo, info::mesh->neighbors)) {
            eslog::internalFailure("cannot synchronize mortar info.\n");
        }
        info.swap(rInfo.back());

        for (size_t n = 0; n < rInfo.size(); ++n) {
            for (size_t i = 0; i < rInfo[n].size(); ++i) {
                if (mInfo.count(rInfo[n][i].from) == 1) {
                    mInfo[rInfo[n][i].from].normal = rInfo[n][i].normal;
                }
            }
        }
    }
}

template <typename T>
void MortarContact<T>::set(const step::Step &step, FETI<T> &feti)
{
    bool run = false;
    for (auto ci = info::ecf->input.contact_interfaces.cbegin(); ci != info::ecf->input.contact_interfaces.cend(); ++ci) {
        if (ci->second.criterion != ContactInterfaceConfiguration::CRITERION::SKIP) {
            run = true;
        }
    }
    if (!run) {
        return;
    }

    std::vector<Mortar> B;
    if (info::mesh->contact->sparseSide != NULL) {
        assembleMortarInterface(B);
    }
    synchronize(feti, B);

//    Communication::serialize([&] () {
//        printf("%2d ::\n", info::mpi::rank);
//        for (size_t i = 0, j, n = 0; i < mortar.size(); ++n) {
//            printf("    :: %3d [%d]<%+.2f> --", mortar[i].from, mInfo[mortar[i].from].id, mInfo[mortar[i].from].normal.z);
//            for (size_t d = 0; d < mInfo[mortar[i].from].domains.size(); ++d) {
//                printf(" %2d", mInfo[mortar[i].from].domains[d]);
//            }
//            printf(" -- ");
//            for (j = i; j < mortar.size() && mortar[i].from == mortar[j].from; ++j) {
//                if (j && mortar[j].to == mortar[j - 1].to) {
//                    printf("[%2d]", mortar[j].domain);
//                } else {
//                    printf(" %3d[%2d]", info::mesh->nodes->IDs->datatarray()[mortar[j].to], mortar[j].domain);
//                }
//            }
//            printf("\n");
//            i = j;
//        }
//    });

    if (info::ecf->input.contact_interfaces.size() > 1) {
        eslog::globalerror("implemented support for more contact pairs\n");
    }

    std::vector<int> bound, gap;
    for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
        switch (it->second.criterion) {
        case ContactInterfaceConfiguration::CRITERION::BOUND:
            bound.insert(bound.end(), it->second.found_interfaces.begin(), it->second.found_interfaces.end()); break;
        case ContactInterfaceConfiguration::CRITERION::GAP:
            gap.insert(gap.end(), it->second.found_interfaces.begin(), it->second.found_interfaces.end()); break;
        default:
            break;
        }
    }
    std::sort(bound.begin(), bound.end());
    std::sort(gap.begin(), gap.end());

    std::vector<std::vector<int> > &D2C = feti.D2C;
    std::vector<std::vector<int> > ROWS(feti.K.size()), COLS(feti.K.size());
    std::vector<std::vector<double> > VALS(feti.K.size()), C(feti.K.size());

    int dofs = info::mesh->dimension;
    if (info::ecf->physics == PhysicsConfiguration::TYPE::HEAT_TRANSFER) {
        dofs = 1;
    }

    std::map<int, int> ids;
    for (auto info = mInfo.cbegin(); info != mInfo.cend(); ++info) {
        ids[info->second.id] = info->first;
    }

    // process mortars according to ID
    feti.lambdas.intervals.push_back({ 0, 0 });
    size_t prevmap = feti.lambdas.cmap.size();
    ineq_begin = feti.lambdas.size;
    for (auto id = ids.cbegin(); id != ids.cend(); ++id) {
        MortarInfo &info = mInfo[id->second];

        auto initrow = [&] () {
            for (size_t d = 0; d < info.domains.size(); ++d) {
                if (feti.decomposition->ismy(info.domains[d])) {
                    D2C [info.domains[d] - feti.decomposition->dbegin].push_back(feti.lambdas.size);
                    ROWS[info.domains[d] - feti.decomposition->dbegin].push_back(0);
                }
            }
        };

        auto push = [&] (const Mortar &m, int dof, double value) {
            auto di = (feti.decomposition->dmap->begin() + dofs * m.to + dof)->begin();
            while (di->domain < m.domain) { ++di; }
            ROWS[m.domain - feti.decomposition->dbegin].back()++;
            COLS[m.domain - feti.decomposition->dbegin].push_back(di->index);
            VALS[m.domain - feti.decomposition->dbegin].push_back(value);
        };

        auto pushcmap = [&] (int multiplier) {
            bool pushcmap = true;
            if (prevmap < feti.lambdas.cmap.size() && feti.lambdas.cmap[prevmap + 1] == (int)info.domains.size()) {
                pushcmap = false;
                for (size_t di = 0; di < info.domains.size(); ++di) {
                    if (info.domains[di] != feti.lambdas.cmap[prevmap + 2 + di]) {
                        pushcmap = true;
                    }
                }
            }
            if (pushcmap) {
                prevmap = feti.lambdas.cmap.size();
                feti.lambdas.cmap.push_back(multiplier);
                feti.lambdas.cmap.push_back(info.domains.size());
                for (size_t di = 0; di < info.domains.size(); ++di) {
                    feti.lambdas.cmap.push_back(info.domains[di]);
                }
            } else {
                feti.lambdas.cmap[prevmap] += multiplier;
            }
            if (info.domains.front() < feti.decomposition->dbegin) {
                feti.lambdas.intervals.back().halo += multiplier;
            }
        };

        if (bound.size()) {
            for (int dof = 0; dof < dofs; ++dof) {
                initrow();
                for (esint i = info.begin; i < info.end; ++i) {
                    push(mortar[i], dof, mortar[i].value);
                }
                feti.lambdas.size++;
                feti.lambdas.intervals.back().size++;
            }
            pushcmap(dofs);
            feti.lambdas.equalities = feti.lambdas.size;
        }
        if (gap.size()) {
            initrow();
            for (esint i = info.begin; i < info.end; ++i) {
                const Mortar &m = mortar[i];
                for (int dof = 0; dof < dofs; ++dof) {
                    if (std::fabs(mInfo[m.from].normal[dof]) > BE_VALUE_TRESHOLD) {
                        push(m, dof, mInfo[m.from].normal[dof] * m.value);
                    }
                }
            }
            feti.lambdas.size++;
            feti.lambdas.intervals.back().size++;
            pushcmap(1);
        }
    }
    feti.lb.resize(feti.lambdas.size - feti.lambdas.equalities);
    feti.ub.resize(feti.lambdas.size - feti.lambdas.equalities);
    math::set(feti.lb, T{0});
    math::set(feti.ub, std::numeric_limits<T>::max());
    ineq_end = feti.lambdas.size;

    std::vector<Matrix_CSR<T> > B1(feti.B1.size());
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        B1[d].resize(feti.B1[d].nrows + ROWS[d].size(), feti.B1[d].ncols, feti.B1[d].nnz + COLS[d].size());
        memcpy(B1[d].rows, feti.B1[d].rows, sizeof(int) * (feti.B1[d].nrows + 1));
        memcpy(B1[d].cols, feti.B1[d].cols, sizeof(int) * feti.B1[d].nnz);
        memcpy(B1[d].vals, feti.B1[d].vals, sizeof(T)   * feti.B1[d].nnz);
        for (size_t i = 0; i < ROWS[d].size(); ++i) {
            B1[d].rows[feti.B1[d].nrows + i + 1] = B1[d].rows[feti.B1[d].nrows + i] + ROWS[d][i];
        }
        memcpy(B1[d].cols + feti.B1[d].nnz, COLS[d].data(), sizeof(int) * COLS[d].size());
        memcpy(B1[d].vals + feti.B1[d].nnz, VALS[d].data(), sizeof(double) * VALS[d].size());
    }

    feti.lambdas.intervals.back().size -= feti.lambdas.intervals.back().halo;

    swap(B1, feti.B1);
    feti.c.resize(feti.lambdas.size);
}

template <typename T>
void MortarContact<T>::update(const step::Step &step, FETI<T> &feti)
{
    bool run = false;
    for (auto ci = info::ecf->input.contact_interfaces.cbegin(); ci != info::ecf->input.contact_interfaces.cend(); ++ci) {
        if (ci->second.criterion != ContactInterfaceConfiguration::CRITERION::SKIP) {
            run = true;
        }
    }
    if (!run) {
        return;
    }

    if (info::ecf->physics != PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS) {
        return;
    }

    std::vector<int> gap;
    for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
        switch (it->second.criterion) {
        case ContactInterfaceConfiguration::CRITERION::GAP:
            gap.insert(gap.end(), it->second.found_interfaces.begin(), it->second.found_interfaces.end()); break;
        default: break;
        }
    }

    if (gap.size() == 0) {
        return;
    }

    int dofs = info::mesh->dimension;
    std::vector<double> disp = StructuralMechanics::Results::displacement->data;
    std::map<esint, double> c;

    std::map<int, int> ids;
    for (auto info = mInfo.cbegin(); info != mInfo.cend(); ++info) {
        ids[info->second.id] = info->first;
    }

    for (auto id = ids.cbegin(); id != ids.cend(); ++id) {
        MortarInfo &info = mInfo[id->second];

        for (esint i = info.begin; i < info.end; ++i) {
            const Mortar &m = mortar[i];
            for (int dof = 0; dof < dofs; ++dof) {
                if (std::fabs(info.normal[dof]) > BE_VALUE_TRESHOLD) {
                    auto di = (feti.decomposition->dmap->begin() + dofs * m.to + dof)->begin();
                    while (di->domain < m.domain) { ++di; }
                    c[id->second] -= info.normal[dof] * m.value * (info::mesh->nodes->coordinates->datatarray()[m.to][dof] + disp[dofs * m.to + dof]);
                }
            }
        }
    }

    std::vector<std::pair<esint, double> > sc;
    std::vector<std::vector<std::pair<esint, double> > > rc(info::mesh->neighbors.size());
    for (auto it = c.cbegin(); it != c.cend(); ++it) {
        sc.push_back(std::make_pair(it->first, it->second));
    }
    if (!Communication::exchangeUnknownSize(sc, rc, info::mesh->neighbors)) {
        eslog::internalFailure("cannot exchange gap.\n");
    }

    for (size_t n = 0; n < rc.size(); ++n) {
        for (size_t i = 0; i < rc[n].size(); ++i) {
            c[rc[n][i].first] += rc[n][i].second;
        }
    }

    double *cc = feti.c.vals + ineq_begin;
    for (auto id = ids.cbegin(); id != ids.cend(); ++id, ++cc) {
        *cc = c[id->second];
    }
}

template struct MortarContact<double>;

}
