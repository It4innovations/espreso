
#include "mortar.h"

#include "analysis/assembler/structuralmechanics.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/contactstore.h"
#include "mesh/store/nodestore.h"
#include "math/math.h"
#include "wrappers/mpi/communication.h"

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
void MortarContact<T>::synchronize(std::vector<Mortar> &B, std::vector<LambdaInfo> &lambdas, std::vector<int> &domains)
{
    // synchronize across neighbors
    std::vector<std::vector<Mortar> > Bs(info::mesh->neighbors.size(), B), Br(info::mesh->neighbors.size());

    if (!Communication::exchangeUnknownSize(Bs, Br, info::mesh->neighbors)) {
        eslog::internalFailure("cannot synchronize mortars.\n");
    }

    for (size_t n = 0; n < Br.size(); ++n) {
        B.insert(B.end(), Br[n].begin(), Br[n].end());
    }

    // unique
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

    // remove non-local lambdas
    for (auto it = B.begin(); it != B.end(); it++) {
        auto &ids = info::mesh->nodes->IDs->datatarray();
        auto nit = std::find(ids.begin(), ids.begin() + info::mesh->nodes->uniqInfo.nhalo, it->to);
        if (nit == ids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != it->to) {
            nit = std::lower_bound(ids.begin() + info::mesh->nodes->uniqInfo.nhalo, ids.end(), it->to);
        }
        if (nit != ids.end() && *nit == it->to) {
            it->to = nit - ids.begin();
            mortar.push_back(*it); // only local part is needed
        }
    }


    // compute lambda info

    // collect my lambdas
    std::vector<Lambda> lmap;
    if (mortar.size()) {
        for (auto it = mortar.begin(), prev = it; it != mortar.end(); ++it) {
            if (it == mortar.begin() || prev->pair != it->pair || prev->from != it->from) {
                auto &myids = info::mesh->nodes->IDs->datatarray();
                auto nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), it->from);
                if (nit != myids.end() && *nit == it->from) {
                    lmap.push_back(Lambda(lmap.size(), it->pair, it->from));
                    if (StructuralMechanics::Results::normal) { // make more robust
                        lmap.back().normal.x = StructuralMechanics::Results::normal->data[info::mesh->dimension * (nit - myids.begin()) + 0];
                        lmap.back().normal.y = StructuralMechanics::Results::normal->data[info::mesh->dimension * (nit - myids.begin()) + 1];
                        lmap.back().normal.z = StructuralMechanics::Results::normal->data[info::mesh->dimension * (nit - myids.begin()) + 2];
                    }
                }
                prev = it;
            }
        }
    }

    esint loffset = lmap.size();
    Communication::exscan(loffset);;
    for (auto lit = lmap.begin(); lit != lmap.end(); ++lit) {
        lit->lambda += loffset;
    }

    // synchronize lambda info
    std::vector<std::vector<Lambda> > rLambdas(info::mesh->neighborsWithMe.size());
    if (!Communication::exchangeUnknownSize(lmap, rLambdas, info::mesh->neighborsWithMe)) {
        eslog::internalFailure("cannot exchange mortar lambdas.\n");
    }

    std::vector<esint> ldomains, sBuffer;
    std::vector<double> sScale;
    auto mit = mortar.begin();
    for (size_t n = 0; n < rLambdas.size(); ++n) {
        for (size_t i = 0; i < rLambdas[n].size(); ++i) {
            ldomains.clear();
            double scale = 0;
            while (mit != mortar.end() && mit->from < rLambdas[n][i].id) { ++mit; }
            if (mit != mortar.end() && mit->from == rLambdas[n][i].id) {
                while (mit != mortar.end() && mit->from == rLambdas[n][i].id) {
                    auto dmap = info::mesh->nodes->domains->begin() + mit->to;
                    for (auto d = dmap->begin(); d != dmap->end(); ++d) {
                        ldomains.push_back(*d);
                    }
                    utils::sortAndRemoveDuplicates(ldomains);
                    scale += dmap->size() * (mit->value / dmap->size()) * (mit->value / dmap->size());
                    ++mit;
                }
            }
            if (ldomains.size()) {
                lambdas.push_back(rLambdas[n][i]);
                sScale.push_back(scale);
                sBuffer.push_back(rLambdas[n][i].id);
                sBuffer.push_back(ldomains.size());
                for (size_t ld = 0; ld < ldomains.size(); ++ld) {
                    sBuffer.push_back(ldomains[ld]);
                }
            }
        }
    }

    std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());
    std::vector<std::vector<double> > rScale(info::mesh->neighborsWithMe.size());
    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighborsWithMe)) {
        eslog::internalFailure("cannot exchange mortar lambda info.\n");
    }
    if (!Communication::exchangeUnknownSize(sScale, rScale, info::mesh->neighborsWithMe)) {
        eslog::internalFailure("cannot exchange mortar lambda scale.\n");
    }

    // try to avoid maps if needed
    std::map<esint, std::vector<esint> > idmap;
    std::map<esint, double> idscale;
    for (size_t n = 0; n < rScale.size(); ++n) {
        for (size_t i = 0, j = 0; i < rScale[n].size(); ++i, j += rBuffer[n][j + 1] + 2) {
            esint id = rBuffer[n][j];
            esint ndomains = rBuffer[n][j + 1];
            idmap[id].insert(idmap[id].end(), rBuffer[n].begin() + j + 2, rBuffer[n].begin() + j + 2 + ndomains);
            idscale[id] += rScale[n][i];
        }
    }

    for (auto ids = idscale.begin(); ids != idscale.end(); ++ids) {
        ids->second = 1 / std::sqrt(ids->second);
    }
    for (auto idm = idmap.begin(); idm != idmap.end(); ++idm) {
        utils::sortAndRemoveDuplicates(idm->second);
    }
    for (size_t i = 0; i < lambdas.size(); ++i) {
        lambdas[i].doffset = domains.size();
        lambdas[i].ndomains = idmap[lambdas[i].id].size();
        for (size_t j = 0; j < idmap[lambdas[i].id].size(); ++j) {
            domains.push_back(idmap[lambdas[i].id][j]);
        }
    }

    for (size_t i = 0; i < mortar.size(); ++i) {
        mortar[i].value *= idscale[mortar[i].from];
    }

//    Communication::serialize([&] () {
//        for (size_t i = 0; i < lambdas.size(); ++i) {
//            printf("%3d :: ", lambdas[i].id);
//            printf("%.3f %.3f %.3f :: ", lambdas[i].normal.x, lambdas[i].normal.y, lambdas[i].normal.z);
//            for (esint j = 0; j < lambdas[i].ndomains; ++j) {
//                printf(" %2d", domains[lambdas[i].doffset + j]);
//            }
//            for (size_t j = 0; j < mortar.size(); ++j) {
//                if (mortar[j].from == lambdas[i].id) {
//                    printf(" <%2d:%2d>=%+.3f", mortar[j].from, mortar[j].to, mortar[j].value);
//                }
//            }
//            printf("\n"); // it does not respect decomposition; hence, it is not 1
//        }
//    });
}

template <typename T>
void MortarContact<T>::set(const step::Step &step, FETI<T> &feti)
{
    if (info::mesh->contactInterfaces.size() == 0) {
        return;
    }

    std::vector<Mortar> B;
    if (info::mesh->contact->sparseSide != NULL) {
        assembleMortarInterface(B);
    }
    std::sort(B.begin(), B.end());
    synchronize(B, lambdas, domains);

    std::vector<int> bound, gap;
    for (auto it = info::ecf->input.contact_interfaces.begin(); it != info::ecf->input.contact_interfaces.end(); ++it) {
        switch (it->second.criterion) {
        case ContactInterfaceConfiguration::CRITERION::BOUND:
            bound.insert(bound.end(), it->second.found_interfaces.begin(), it->second.found_interfaces.end()); break;
        case ContactInterfaceConfiguration::CRITERION::GAP:
            gap.insert(gap.end(), it->second.found_interfaces.begin(), it->second.found_interfaces.end()); break;
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

    auto begin = mortar.begin(), end = begin;
    auto getarray = [&] (const LambdaInfo &lambda, std::vector<int> &filter) {
        if (!std::binary_search(filter.begin(), filter.end(), lambda.pair)) { return false; }
        while (begin != mortar.end() && begin->from  < lambda.id) ++begin;
        end = begin;
        while (  end != mortar.end() &&   end->from == lambda.id) ++end;
        return true;
    };

    auto initrow = [&] (const LambdaInfo &lambda) {
        for (esint di = 0; di < lambda.ndomains; ++di) {
            if (feti.decomposition->ismy(domains[lambda.doffset + di])) {
                esint ldi = domains[lambda.doffset + di] - feti.decomposition->dbegin;
                D2C [ldi].push_back(feti.lambdas.size);
                ROWS[ldi].push_back(0);
            }
        }
    };

    size_t prevmap = feti.lambdas.cmap.size();
    auto pushcmap = [&] (const LambdaInfo &lambda, int multiplier) {
        bool pushcmap = true;
        if (prevmap < feti.lambdas.cmap.size() && feti.lambdas.cmap[prevmap + 1] == lambda.ndomains) {
            pushcmap = false;
            for (esint di = 0; di < lambda.ndomains; ++di) {
                if (domains[lambda.doffset + di] != feti.lambdas.cmap[prevmap + 2 + di]) {
                    pushcmap = true;
                }
            }
        }
        if (pushcmap) {
            prevmap = feti.lambdas.cmap.size();
            feti.lambdas.cmap.push_back(multiplier); // TODO: sum up the lambdas with the same map
            feti.lambdas.cmap.push_back(lambda.ndomains);
            for (esint di = 0; di < lambda.ndomains; ++di) {
                feti.lambdas.cmap.push_back(domains[lambda.doffset + di]);
            }
        } else {
            feti.lambdas.cmap[prevmap] += multiplier;
        }
        if (domains[lambda.doffset] < feti.decomposition->dbegin) {
            feti.lambdas.intervals.back().halo += multiplier;
        }
    };

    if (bound.size()) {
        prevmap = feti.lambdas.cmap.size();
        begin = mortar.begin(), end = begin;
        feti.lambdas.intervals.push_back({ 0, 0 });
        for (auto lambda = lambdas.begin(); lambda != lambdas.end(); ++lambda) { // equality constrains
            if (!getarray(*lambda, bound)) { continue; }

            for (int dof = 0; dof < dofs; ++dof) {
                initrow(*lambda);
                for (auto it = begin; it != end; ++it) {
                    auto dmap = feti.decomposition->dmap->begin() + dofs * it->to + dof;
                    for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                        if (feti.decomposition->ismy(di->domain)) {
                            esint ldi = di->domain - feti.decomposition->dbegin;
                            ROWS[ldi].back()++;
                            COLS[ldi].push_back(di->index);
                            VALS[ldi].push_back(it->value / dmap->size());
                        }
                    }
                }
                feti.lambdas.size++;
                feti.lambdas.intervals.back().size++;
            }
            pushcmap(*lambda, dofs);
        }
        feti.lambdas.equalities = feti.lambdas.size;
    }

    if (gap.size()) {
        ineq_begin = feti.lambdas.size;
        prevmap = feti.lambdas.cmap.size();
        begin = mortar.begin(), end = begin;
        feti.lambdas.intervals.push_back({ 0, 0 });
        for (auto lambda = lambdas.begin(); lambda != lambdas.end(); ++lambda) { // equality constrains
            if (!getarray(*lambda, gap)) { continue; }

            initrow(*lambda);
            for (auto it = begin; it != end; ++it) {
                for (int dof = 0; dof < dofs; ++dof) {
                    if (std::fabs(lambda->normal[dof]) > BE_VALUE_TRESHOLD) {
                        auto dmap = feti.decomposition->dmap->begin() + dofs * it->to + dof;
                        for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                            if (feti.decomposition->ismy(di->domain)) {
                                esint ldi = di->domain - feti.decomposition->dbegin;
                                ROWS[ldi].back()++;
                                COLS[ldi].push_back(di->index);
                                VALS[ldi].push_back(lambda->normal[dof] * it->value / dmap->size());
                            }
                        }
                    }
                }
            }
            feti.lambdas.size++;
            feti.lambdas.intervals.back().size++;
            pushcmap(*lambda, 1);
        }
        feti.lb.resize(feti.lambdas.size - feti.lambdas.equalities);
        feti.ub.resize(feti.lambdas.size - feti.lambdas.equalities);
        math::set(feti.lb, T{0});
        math::set(feti.ub, std::numeric_limits<T>::max());
        ineq_end = feti.lambdas.size;
    }

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
    if (info::mesh->contactInterfaces.size() == 0) {
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
    std::vector<std::pair<esint, double> > c;

    auto begin = mortar.begin(), end = begin;
    for (auto lambda = lambdas.begin(); lambda != lambdas.end(); ++lambda) { // equality constrains
        if (!std::binary_search(gap.begin(), gap.end(), lambda->pair)) { continue; }
        while (begin != mortar.end() && begin->from  < lambda->id) ++begin;
        end = begin;
        while (  end != mortar.end() &&   end->from == lambda->id) ++end;

        c.push_back(std::make_pair(begin->from, 0));
        for (auto it = begin; it != end; ++it) {
            for (int dof = 0; dof < dofs; ++dof) {
                c.back().second -= lambda->normal[dof] * it->value * (info::mesh->nodes->coordinates->datatarray()[it->to][dof] + disp[dofs * it->to + dof]);
            }
        }
    }

    Vector_Dense<T> cc; cc.vals = feti.c.vals + ineq_begin; cc.size = ineq_end - ineq_begin;
    math::set(cc, T{0});

    std::vector<std::vector<std::pair<esint, double> > > rc(info::mesh->neighbors.size());
    if (!Communication::exchangeUnknownSize(c, rc, info::mesh->neighbors)) {
        eslog::internalFailure("cannot exchange gap.\n");
    }
    for (size_t n = 0; n < rc.size(); ++n) {
        for (size_t i = 0, j = 0; i < rc[n].size(); ++i) {
            while (j < c.size() && c[j].first < rc[n][i].first) { ++j; }
            if (j < c.size() && c[j].first == rc[n][i].first) {
                c[j].second += rc[n][i].second;
            }
        }
    }

    for (size_t i = 0; i < c.size(); ++i) {
        cc.vals[i] += c[i].second;
    }
}

template class MortarContact<double>;

}
