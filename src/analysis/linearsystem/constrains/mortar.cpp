
#include "mortar.h"

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
void MortarContact<T>::synchronize(std::vector<Mortar> &B)
{
    // synchronize across neighbors
    std::vector<std::vector<Mortar> > Bs(info::mesh->neighbors.size(), B), Br(info::mesh->neighbors.size());

    if (!Communication::exchangeUnknownSize(Bs, Br, info::mesh->neighbors)) {
        eslog::internalFailure("cannot synchronize mortars.\n");
    }

    for (size_t n = 0; n < Br.size(); ++n) {
        B.insert(B.end(), Br[n].begin(), Br[n].end());
    }

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

    for (auto begin = B.begin(), end = begin; begin != B.end(); begin = end) {
        double scale = 1;
        while (end != B.end() && begin->pair == end->pair && begin->from == end->from) {
            if (end->from == end->to) {
                scale = end->value;
            }
            ++end;
        }
        scale /= M_SQRT1_2;
        auto &ids = info::mesh->nodes->IDs->datatarray();
        bool push = false;
        for (auto it = begin; it != end; ++it) {
            it->value /= scale;
            auto nit = std::find(ids.begin(), ids.begin() + info::mesh->nodes->uniqInfo.nhalo, it->to);
            if (nit == ids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != it->to) {
                nit = std::lower_bound(ids.begin() + info::mesh->nodes->uniqInfo.nhalo, ids.end(), it->to);
            }
            if (nit != ids.end() && *nit == it->to) {
                push = true;
            }
        }
        if (push) {
            mortar.insert(mortar.end(), begin, end);
        }
    }
}

template <typename T>
void MortarContact<T>::set(const step::Step &step, FETI<T> &feti, int dofs)
{
    if (info::mesh->contactInterfaces.size() == 0) {
        return;
    }

    std::vector<Mortar> B;
    if (info::mesh->contact->sparseSide != NULL) {
        assembleMortarInterface(B);
    }
    std::sort(B.begin(), B.end());
    synchronize(B);

    struct __lambda__ { esint lambda, pair, id; };

    esint lambdaCount = 0;
    std::vector<__lambda__> lmap;
    auto &myids = info::mesh->nodes->IDs->datatarray();

    // collect my nodes in B in order to compute mortar lambdas
    if (mortar.size()) {
        for (auto it = mortar.begin(), prev = it; it != mortar.end(); ++it) {
            if (it == mortar.begin() || prev->pair != it->pair || prev->from != it->from) {
                auto nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), it->from);
                if (nit != myids.end() && *nit == it->from) {
                    lmap.push_back({ lambdaCount++, it->pair, it->from });
                }
                prev = it;
            }
        }
    }

    Communication::exscan(lambdaCount);
    for (auto lit = lmap.begin(); lit != lmap.end(); ++lit) {
        lit->lambda += lambdaCount;
    }

    std::vector<std::vector<__lambda__> > sLambdas(info::mesh->neighborsWithMe.size(), lmap), rLambdas(info::mesh->neighborsWithMe.size());
    if (!Communication::exchangeUnknownSize(sLambdas, rLambdas, info::mesh->neighborsWithMe)) {
        eslog::internalFailure("cannot exchange mortar lambdas.\n");
    }

    // lambdas <node, lambda> are sorted according to lambda in order to get correct result from FETI solver
    std::vector<__lambda__> lambdas;
    for (size_t n = 0; n < rLambdas.size(); ++n) {
        lambdas.insert(lambdas.end(), rLambdas[n].begin(), rLambdas[n].end());
    }
    std::sort(lambdas.begin(), lambdas.end(), [] (const __lambda__ &l1, const __lambda__ &l2) { return l1.lambda < l2.lambda; });

//  Communication::serialize([&] () {
//      std::cout << info::mpi::rank << ": \n";
//      for (size_t i = 0; i < lambdas.size(); ++i) {
//          printf("%d %d %d\n", lambdas[i].lambda, lambdas[i].pair, lambdas[i].id);
//      }
//  });

    // found nodes that are on mortar interface -> we have exchange domain map in order to get correct decomposition
    std::vector<esint> mids, mylambdas;
    for (auto it = mortar.begin(); it != mortar.end(); ++it) {
        auto nit = std::find(myids.begin(), myids.begin() + info::mesh->nodes->uniqInfo.nhalo, it->to);
        if (nit == myids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != it->to) {
            nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), it->to);
        }
        if (nit != myids.end() && *nit == it->to) {
            mids.push_back(it->to);
            if (mylambdas.size() == 0 || mylambdas.back() != it->from) {
                mylambdas.push_back(it->from);
            }
        }
    }
    utils::sortAndRemoveDuplicates(mylambdas);
    utils::sortAndRemoveDuplicates(mids);

    std::vector<esint> sBuffer;
    std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());
    sBuffer.push_back(mids.size());
    for (auto it = mids.begin(); it != mids.end(); ++it) {
        auto nit = std::find(myids.begin(), myids.begin() + info::mesh->nodes->uniqInfo.nhalo, *it);
        if (nit == myids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != *it) {
            nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), *it);
        }
        if (nit != myids.end() && *nit == *it) {
            sBuffer.push_back(*it);
            auto dmap = feti.decomposition->dmap->begin() + dofs * (nit - myids.begin());
            sBuffer.push_back(dmap->size());
            for (int dof = 0; dof < dofs; ++dof, ++dmap) {
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        sBuffer.push_back(di->domain);
                        sBuffer.push_back(di->index);
                    }
                }
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighborsWithMe)) {
        eslog::internalFailure("cannot exchange mortar d-map.\n");
    }

    struct npair { esint id, n, offset; };
    std::vector<npair> ninfo; // where we have stored info about a particular node
    for (size_t n = 0; n < rBuffer.size(); ++n) {
        for (size_t i = 1; i < rBuffer[n].size(); i += 2 * dofs * rBuffer[n][i + 1] + 2) {
            ninfo.push_back(npair{ rBuffer[n][i], (esint)n, (esint)i });
        }
    }
    std::sort(ninfo.begin(), ninfo.end(), [&] (const npair &i, const npair &j) {
        if (i.id == j.id) {
            return i.n < j.n;
        }
        return i.id < j.id;
    });

    feti.lambdas.intervals.push_back({ 0, 0 });
    std::vector<std::vector<int> > &D2C = feti.D2C;
    std::vector<std::vector<int> > ROWS(feti.K.size()), COLS(feti.K.size());
    std::vector<std::vector<double> > VALS(feti.K.size());

    // build B1 in correct order
    size_t prevmap = feti.lambdas.cmap.size();
    std::vector<int> domains;
    for (auto lambda = lambdas.begin(); lambda != lambdas.end(); ++lambda) {
        if (!std::binary_search(mylambdas.begin(), mylambdas.end(), lambda->id)) {
            continue;
        }
        auto begin = std::lower_bound(mortar.begin(), mortar.end(), *lambda, [] (const Mortar &b, const __lambda__ &l) { return b.pair == l.pair ? b.from < l.id : b.pair < l.pair; });
        auto end = begin;
        while (end != mortar.end() && lambda->pair == end->pair && lambda->id == end->from) {
            ++end;
        }


        domains.clear();
        for (int dof = 0; dof < dofs; ++dof) {
            for (auto it = begin; it != end; ++it) {
                auto nit = std::lower_bound(ninfo.begin(), ninfo.end(), it->to, [&] (const npair &info, const esint &lambda) { return info.id < lambda; });
                int ndomains = 0;
                for (auto nnit = nit; nnit != ninfo.end() && nnit->id == it->to; ++nnit) {
                    ndomains += rBuffer[nnit->n][nnit->offset + 1];
                    if (dof == 0) {
                        for (esint i = 0; i < rBuffer[nnit->n][nnit->offset + 1]; ++i) {
                            domains.push_back(rBuffer[nit->n][nit->offset + 2 + 2 * i]);
                        }
                    }
                }
                while (nit != ninfo.end() && nit->id == it->to) {
                    if (info::mesh->neighborsWithMe[nit->n] == info::mpi::rank) {
                        for (esint i = 0; i < rBuffer[nit->n][nit->offset + 1]; ++i) {
                            if (D2C [rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].size() == 0 || D2C[rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].back() != feti.lambdas.size) {
                                D2C [rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].push_back(feti.lambdas.size);
                                ROWS[rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].push_back(0);
                            }
                            ROWS[rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].back()++;
                            COLS[rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].push_back(rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i + 1]);
                            VALS[rBuffer[nit->n][nit->offset + 2 * dof * ndomains + 2 + 2 * i] - feti.decomposition->dbegin].push_back(it->value / ndomains);
                        }
                    }
                    ++nit;
                }
            }
            feti.lambdas.size++;
            feti.lambdas.intervals.back().size++;
        }
        utils::sortAndRemoveDuplicates(domains);
        bool pushcmap = true;
        if (prevmap < feti.lambdas.cmap.size() && feti.lambdas.cmap[prevmap + 1] == (int)domains.size()) {
            pushcmap = false;
            for (size_t i = 0; i < domains.size(); ++i) {
                if (domains[i] != feti.lambdas.cmap[prevmap + 2 + i]) {
                    pushcmap = true;
                }
            }
        }
        if (pushcmap) {
            prevmap = feti.lambdas.cmap.size();
            feti.lambdas.cmap.push_back(dofs); // TODO: sum up the lambdas with the same map
            feti.lambdas.cmap.push_back(domains.size());
            for (size_t i = 0; i < domains.size(); ++i) {
                feti.lambdas.cmap.push_back(domains[i]);
            }
        } else {
            feti.lambdas.cmap[prevmap] += dofs;
        }
        if (domains.front() < feti.decomposition->dbegin) {
            feti.lambdas.intervals.back().halo += dofs;
        }
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
}

template class MortarContact<double>;

}
