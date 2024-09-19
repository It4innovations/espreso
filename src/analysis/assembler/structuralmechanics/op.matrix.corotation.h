
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixCorotation: SubKernel {
    const char* name() const { return "MatrixCorotation"; }

    MatrixCorotation()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate()
    {
        this->isactive = 1;
    }
};


template <Element::CODE code, size_t nodes, size_t gps, size_t ndim> struct MatrixCorotationKernel;

template <Element::CODE code, size_t nodes, size_t gps>
struct MatrixCorotationKernel<code, nodes, gps, 2>: MatrixCorotation {
    MatrixCorotationKernel(const MatrixCorotation &base): MatrixCorotation(base)
    {
        BasisKernel<code, nodes, gps, 2>::setCenter(cw, cN, cdN);
        n_div = 1. / nodes;
    }

    double cw, cN[nodes], cdN[nodes][2];
    double n_div;

    template <typename Element>
    void simd(Element &element)
    {
        SIMD JC[2];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    JC[i * 2 + j] = JC[i * 2 + j] + element.coords.node[n][i] * load1(cdN[n][j]);
                }
            }
        }

        SIMD det, invJC[4]; inv22(JC, det, invJC);

        SIMD cdND[nodes][2];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(cdN[n][0]);
            SIMD dNY = load1(cdN[n][1]);
            cdND[n][0] = invJC[0] * dNX + invJC[2] * dNY;
            cdND[n][1] = invJC[1] * dNX + invJC[3] * dNY;
        }

        SIMD F[4];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    F[i * 2 + j] = F[i * 2 + j] + element.displacement[n][i] * cdND[n][j];
                }
            }
        }

        F[0] = F[0] + load1(1.0);
        F[3] = F[3] + load1(1.0);

        // R = F * (F' * F)^(-1/2)
        SIMD FF[4]; multAtB<2, 2, 2>(FF, F, F);
        SIMD eigVal[2], eigVec[2];
        eigSym22(FF, eigVal, eigVec);
        SIMD FD[4] = { sqrt(eigVal[0]), zeros(), zeros(), sqrt(eigVal[1]) };
        set<2, 2>(FF, zeros());
        multAtBA<2, 2>(FF, eigVec, FD);
        SIMD invFF[4]; inv22(FF, invFF);
        SIMD R[4]; multAB<2, 2, 2>(R, F, invFF);

        SIMD cdisp[2], corigin[2];
        for (int d = 0; d < 2; ++d) {
            for (size_t n = 0; n < nodes; ++n) {
                cdisp[d] = cdisp[d] + element.coords.node[n][d] + element.displacement[n][d];
                corigin[d] = corigin[d] + element.coords.node[n][d];
            }
            cdisp[d] = cdisp[d] * load1(n_div);
            corigin[d] = corigin[d] * load1(n_div);
        }

        SIMD u_lie[2 * nodes], S_bar[2 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD ldisp[2], lorigin[2];
            for (int d = 0; d < 2; ++d) {
                ldisp[d] = element.coords.node[n][d] + element.displacement[n][d] - cdisp[d];
                lorigin[d] = element.coords.node[n][d] - corigin[d];
            }
            SIMD Mdisp[2] = { -ldisp[1], ldisp[0] };

            // S_bar = -R' * vektor od teziste (deformovane)
            SIMD S_barn[2]; multAtBt<2, 2, 1>(S_barn, R, Mdisp, load1(-1.));
            for (int i = 0; i < 2; ++i) {
                S_bar[0 * nodes + n] = S_barn[0];
                S_bar[1 * nodes + n] = S_barn[1];
            }

            // u_lie = Re' * vektor od teziste (deformovane) - vektor od teziste (nedef.)
            SIMD u_lien[2]; multAtB<2, 2, 1>(u_lien, R, ldisp);
            u_lie[n + 0 * nodes] = u_lien[0] - lorigin[0];
            u_lie[n + 1 * nodes] = u_lien[1] - lorigin[1];
        }

        // A_transp
        // [  0 -3  2 ]
        // [  3  0 -1 ]
        // [ -2  1  0 ]

        // internal force computation
        SIMD AS[4];
        for (size_t n = 0; n < nodes; ++n) {
            for (size_t i = 0; i < 2; ++i) {
                AS[0 * 2 + i] = AS[0 * 2 + i] + cdND[n][1] * S_bar[0 * nodes + n];
                AS[1 * 2 + i] = AS[1 * 2 + i] - cdND[n][0] * S_bar[1 * nodes + n];
            }
        }
        SIMD invAS[4]; inv22(AS, invAS);

        // G_bar = inv(A_trans * S_bar) * A_transp
        SIMD G_bar[1 * 2 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            for (size_t i = 0; i < 2; ++i) {
                G_bar[i * nodes + n] = invAS[i * 2 + 0] * cdND[n][1] - invAS[i * 2 + 1] * cdND[n][0];
            }
        }

        // P_bar = eye(24) - S_bar * G_bar
        SIMD P_bar[2 * nodes * 2 * nodes];
        for (size_t i = 0; i < 2 * nodes; ++i) {
            P_bar[i * 2 * nodes + i] = load1(1.);
        }
        multAB<2 * nodes, 2, 2 * nodes>(P_bar, S_bar, G_bar, load1(-1.));

        // f_ile = Ke * u_lie
        SIMD f_ile[2 * nodes]; multAB<2 * nodes, 2 * nodes, 1>(f_ile, element.K, u_lie);
        SIMD fiep[2 * nodes]; multAtB<2 * nodes, 2 * nodes, 1>(fiep, P_bar, f_ile);
        for (size_t n = 0; n < nodes; ++n) {
            element.nf[0 * nodes + n] = element.nf[0 * nodes + n] + R[0] * fiep[0 * nodes + n] + R[1] * fiep[1 * nodes + n];
            element.nf[1 * nodes + n] = element.nf[1 * nodes + n] + R[2] * fiep[0 * nodes + n] + R[3] * fiep[1 * nodes + n];
        }

        SIMD Fn_bar[2 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            Fn_bar[0 * nodes + n] = -fiep[n + 1 * nodes];
            Fn_bar[1 * nodes + n] =  fiep[n + 0 * nodes];
        }

        // KM =  Re * P_bar' * Ke * P_bar      * Re'
        // KGR = Re * F_bar * G_bar            * Re'
        // KGP = Re * G_bar' * Fn_bar' * P_bar * Re'
        // Kte = KM - KGR - KGP

        // K = Re * ((P_bar' * Ke - G_bar' * Fn_bar') * P_bar - Fn_bar * G_bar) * Re'
        SIMD Ke[2 * nodes * 2 * nodes];
        multAtB <2 * nodes, 2 * nodes, 2 * nodes>(Ke, P_bar, element.K);
        set<2 * nodes, 2 * nodes>(element.K, load1(0));
        multAB  <2 * nodes, 2        , 2 * nodes>(element.K, Fn_bar, G_bar, load1(-1));
        for (size_t n = 0; n < 2 * nodes; ++n) { // (P_bar' * Ke) - (G_bar' * Fn_bar')
            for (size_t m = 0; m < 2 * nodes; ++m) {
                Ke[n * 2 * nodes + m] = Ke[n * 2 * nodes + m] + element.K[m * 2 * nodes + n];
            }
        }
        multAB  <2 * nodes, 2 * nodes, 2 * nodes>(element.K, Ke, P_bar);

        for (size_t n = 0; n < nodes; ++n) {
            for (size_t m = 0; m < nodes; ++m) {
                SIMD RK11 = R[0] * element.K[(0 * nodes + n) * nodes * 2 + 0 * nodes + m] + R[1] * element.K[(1 * nodes + n) * nodes * 2 + 0 * nodes + m];
                SIMD RK12 = R[0] * element.K[(0 * nodes + n) * nodes * 2 + 1 * nodes + m] + R[1] * element.K[(1 * nodes + n) * nodes * 2 + 1 * nodes + m];

                SIMD RK21 = R[2] * element.K[(0 * nodes + n) * nodes * 2 + 0 * nodes + m] + R[3] * element.K[(1 * nodes + n) * nodes * 2 + 0 * nodes + m];
                SIMD RK22 = R[2] * element.K[(0 * nodes + n) * nodes * 2 + 1 * nodes + m] + R[3] * element.K[(1 * nodes + n) * nodes * 2 + 1 * nodes + m];

                element.K[(0 * nodes + n) * nodes * 2 + 0 * nodes + m] = RK11 * R[0] + RK12 * R[1];
                element.K[(0 * nodes + n) * nodes * 2 + 1 * nodes + m] = RK11 * R[2] + RK12 * R[3];

                element.K[(1 * nodes + n) * nodes * 2 + 0 * nodes + m] = RK21 * R[0] + RK22 * R[1];
                element.K[(1 * nodes + n) * nodes * 2 + 1 * nodes + m] = RK21 * R[2] + RK22 * R[3];
            }
        }
    }
};

template <Element::CODE code, size_t nodes, size_t gps>
struct MatrixCorotationKernel<code, nodes, gps, 3>: MatrixCorotation {
    MatrixCorotationKernel(const MatrixCorotation &base): MatrixCorotation(base)
    {
        BasisKernel<code, nodes, gps, 3>::setCenter(cw, cN, cdN);
        n_div = 1. / nodes;
    }

    double cw, cN[nodes], cdN[nodes][3];
    double n_div;

    template <typename Element>
    void simd(Element &element)
    {
        SIMD JC[9];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    JC[i * 3 + j] = JC[i * 3 + j] + element.coords.node[n][i] * load1(cdN[n][j]);
                }
            }
        }

        SIMD det, invJC[9]; inv33(JC, det, invJC);

        SIMD cdND[nodes][3];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(cdN[n][0]);
            SIMD dNY = load1(cdN[n][1]);
            SIMD dNZ = load1(cdN[n][2]);
            cdND[n][0] = invJC[0] * dNX + invJC[3] * dNY + invJC[6] * dNZ;
            cdND[n][1] = invJC[1] * dNX + invJC[4] * dNY + invJC[7] * dNZ;
            cdND[n][2] = invJC[2] * dNX + invJC[5] * dNY + invJC[8] * dNZ;
        }

        SIMD F[9];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    F[i * 3 + j] = F[i * 3 + j] + element.displacement[n][i] * cdND[n][j];
                }
            }
        }

        F[0] = F[0] + load1(1.0);
        F[4] = F[4] + load1(1.0);
        F[8] = F[8] + load1(1.0);

        // R = F * (F' * F)^(-1/2)
        SIMD FF[9]; multAtB<3, 3, 3>(FF, F, F);
        SIMD eigVal[3], eigVec[9];
        eigSym33(FF, eigVal, eigVec);
        SIMD FD[9] = { sqrt(eigVal[0]), zeros(), zeros(), zeros(), sqrt(eigVal[1]), zeros(), zeros(), zeros(), sqrt(eigVal[2]) };
        set<3, 3>(FF, zeros());
        multAtBA<3, 3>(FF, eigVec, FD);
        SIMD invFF[9]; inv33(FF, invFF);
        SIMD R[9]; multAB<3, 3, 3>(R, F, invFF);

        SIMD cdisp[3], corigin[3];
        for (int d = 0; d < 3; ++d) {
            for (size_t n = 0; n < nodes; ++n) {
                cdisp[d] = cdisp[d] + element.coords.node[n][d] + element.displacement[n][d];
                corigin[d] = corigin[d] + element.coords.node[n][d];
            }
            cdisp[d] = cdisp[d] * load1(n_div);
            corigin[d] = corigin[d] * load1(n_div);
        }

        SIMD u_lie[3 * nodes], S_bar[3 * nodes * 3];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD ldisp[3], lorigin[3];
            for (int d = 0; d < 3; ++d) {
                ldisp[d] = element.coords.node[n][d] + element.displacement[n][d] - cdisp[d];
                lorigin[d] = element.coords.node[n][d] - corigin[d];
            }
            SIMD Mdisp[9] = { zeros(), -ldisp[2], ldisp[1], ldisp[2], zeros(), -ldisp[0], -ldisp[1], ldisp[0], zeros() };

            // S_bar = -R' * vektor od teziste (deformovane) * R
            SIMD S_barn[9]; multAtBA<3, 3>(S_barn, R, Mdisp, load1(-1.));
            for (int i = 0; i < 3; ++i) {
                S_bar[(0 * nodes + n) * 3 + i] = S_barn[3 * 0 + i];
                S_bar[(1 * nodes + n) * 3 + i] = S_barn[3 * 1 + i];
                S_bar[(2 * nodes + n) * 3 + i] = S_barn[3 * 2 + i];
            }

            // u_lie = Re' * vektor od teziste (deformovane) - vektor od teziste (nedef.)
            SIMD u_lien[3]; multAtB<3, 3, 1>(u_lien, R, ldisp);
            u_lie[n + 0 * nodes] = u_lien[0] - lorigin[0];
            u_lie[n + 1 * nodes] = u_lien[1] - lorigin[1];
            u_lie[n + 2 * nodes] = u_lien[2] - lorigin[2];
        }

        // A_transp
        // [  0 -3  2 ]
        // [  3  0 -1 ]
        // [ -2  1  0 ]

        // internal force computation
        SIMD AS[9];
        for (size_t n = 0; n < nodes; ++n) {
            for (size_t i = 0; i < 3; ++i) {
                AS[0 * 3 + i] = AS[0 * 3 + i] - cdND[n][2] * S_bar[(1 * nodes + n) * 3 + i] + cdND[n][1] * S_bar[(2 * nodes + n) * 3 + i];
                AS[1 * 3 + i] = AS[1 * 3 + i] + cdND[n][2] * S_bar[(0 * nodes + n) * 3 + i] - cdND[n][0] * S_bar[(2 * nodes + n) * 3 + i];
                AS[2 * 3 + i] = AS[2 * 3 + i] - cdND[n][1] * S_bar[(0 * nodes + n) * 3 + i] + cdND[n][0] * S_bar[(1 * nodes + n) * 3 + i];
            }
        }
        SIMD invAS[9]; inv33(AS, invAS);

        // G_bar = inv(A_trans * S_bar) * A_transp
        SIMD G_bar[3 * 3 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            for (size_t i = 0; i < 3; ++i) {
                G_bar[i * 3 * nodes + 0 * nodes + n] =  invAS[i * 3 + 1] * cdND[n][2] - invAS[i * 3 + 2] * cdND[n][1];
                G_bar[i * 3 * nodes + 1 * nodes + n] = -invAS[i * 3 + 0] * cdND[n][2] + invAS[i * 3 + 2] * cdND[n][0];
                G_bar[i * 3 * nodes + 2 * nodes + n] =  invAS[i * 3 + 0] * cdND[n][1] - invAS[i * 3 + 1] * cdND[n][0];
            }
        }

        // P_bar = eye(24) - S_bar * G_bar
        SIMD P_bar[3 * nodes * 3 * nodes];
        for (size_t i = 0; i < 3 * nodes; ++i) {
            P_bar[i * 3 * nodes + i] = load1(1.);
        }
        multAB<3 * nodes, 3, 3 * nodes>(P_bar, S_bar, G_bar, load1(-1.));

        // f_ile = Ke * u_lie
        SIMD f_ile[3 * nodes]; multAB<3 * nodes, 3 * nodes, 1>(f_ile, element.K, u_lie);
        SIMD fiep[3 * nodes]; multAtB<3 * nodes, 3 * nodes, 1>(fiep, P_bar, f_ile);
        for (size_t n = 0; n < nodes; ++n) {
            element.nf[0 * nodes + n] = element.nf[0 * nodes + n] + R[0] * fiep[0 * nodes + n] + R[1] * fiep[1 * nodes + n] + R[2] * fiep[2 * nodes + n];
            element.nf[1 * nodes + n] = element.nf[1 * nodes + n] + R[3] * fiep[0 * nodes + n] + R[4] * fiep[1 * nodes + n] + R[5] * fiep[2 * nodes + n];
            element.nf[2 * nodes + n] = element.nf[2 * nodes + n] + R[6] * fiep[0 * nodes + n] + R[7] * fiep[1 * nodes + n] + R[8] * fiep[2 * nodes + n];
        }

        SIMD Fn_bar[3 * nodes * 3];
        for (size_t n = 0; n < nodes; ++n) {
            Fn_bar[(n + 0 * nodes) * 3 + 1] = -fiep[n + 2 * nodes];
            Fn_bar[(n + 0 * nodes) * 3 + 2] =  fiep[n + 1 * nodes];
            Fn_bar[(n + 1 * nodes) * 3 + 0] =  fiep[n + 2 * nodes];
            Fn_bar[(n + 1 * nodes) * 3 + 2] = -fiep[n + 0 * nodes];
            Fn_bar[(n + 2 * nodes) * 3 + 0] = -fiep[n + 1 * nodes];
            Fn_bar[(n + 2 * nodes) * 3 + 1] =  fiep[n + 0 * nodes];
        }

        // KM =  Re * P_bar' * Ke * P_bar      * Re'
        // KGR = Re * F_bar * G_bar            * Re'
        // KGP = Re * G_bar' * Fn_bar' * P_bar * Re'
        // Kte = KM - KGR - KGP

        // K = Re * ((P_bar' * Ke - G_bar' * Fn_bar') * P_bar - Fn_bar * G_bar) * Re'
        SIMD Ke[3 * nodes * 3 * nodes];
        multAtB <3 * nodes, 3 * nodes, 3 * nodes>(Ke, P_bar, element.K);
        set<3 * nodes, 3 * nodes>(element.K, load1(0));
        multAB  <3 * nodes, 3        , 3 * nodes>(element.K, Fn_bar, G_bar, load1(-1));
        for (size_t n = 0; n < 3 * nodes; ++n) { // (P_bar' * Ke) - (G_bar' * Fn_bar')
            for (size_t m = 0; m < 3 * nodes; ++m) {
                Ke[n * 3 * nodes + m] = Ke[n * 3 * nodes + m] + element.K[m * 3 * nodes + n];
            }
        }
        multAB  <3 * nodes, 3 * nodes, 3 * nodes>(element.K, Ke, P_bar);

        for (size_t n = 0; n < nodes; ++n) {
            for (size_t m = 0; m < nodes; ++m) {
                SIMD RK11 = R[0] * element.K[(0 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[1] * element.K[(1 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[2] * element.K[(2 * nodes + n) * nodes * 3 + 0 * nodes + m];
                SIMD RK12 = R[0] * element.K[(0 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[1] * element.K[(1 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[2] * element.K[(2 * nodes + n) * nodes * 3 + 1 * nodes + m];
                SIMD RK13 = R[0] * element.K[(0 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[1] * element.K[(1 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[2] * element.K[(2 * nodes + n) * nodes * 3 + 2 * nodes + m];

                SIMD RK21 = R[3] * element.K[(0 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[4] * element.K[(1 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[5] * element.K[(2 * nodes + n) * nodes * 3 + 0 * nodes + m];
                SIMD RK22 = R[3] * element.K[(0 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[4] * element.K[(1 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[5] * element.K[(2 * nodes + n) * nodes * 3 + 1 * nodes + m];
                SIMD RK23 = R[3] * element.K[(0 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[4] * element.K[(1 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[5] * element.K[(2 * nodes + n) * nodes * 3 + 2 * nodes + m];

                SIMD RK31 = R[6] * element.K[(0 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[7] * element.K[(1 * nodes + n) * nodes * 3 + 0 * nodes + m] + R[8] * element.K[(2 * nodes + n) * nodes * 3 + 0 * nodes + m];
                SIMD RK32 = R[6] * element.K[(0 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[7] * element.K[(1 * nodes + n) * nodes * 3 + 1 * nodes + m] + R[8] * element.K[(2 * nodes + n) * nodes * 3 + 1 * nodes + m];
                SIMD RK33 = R[6] * element.K[(0 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[7] * element.K[(1 * nodes + n) * nodes * 3 + 2 * nodes + m] + R[8] * element.K[(2 * nodes + n) * nodes * 3 + 2 * nodes + m];

                element.K[(0 * nodes + n) * nodes * 3 + 0 * nodes + m] = RK11 * R[0] + RK12 * R[1] + RK13 * R[2];
                element.K[(0 * nodes + n) * nodes * 3 + 1 * nodes + m] = RK11 * R[3] + RK12 * R[4] + RK13 * R[5];
                element.K[(0 * nodes + n) * nodes * 3 + 2 * nodes + m] = RK11 * R[6] + RK12 * R[7] + RK13 * R[8];

                element.K[(1 * nodes + n) * nodes * 3 + 0 * nodes + m] = RK21 * R[0] + RK22 * R[1] + RK23 * R[2];
                element.K[(1 * nodes + n) * nodes * 3 + 1 * nodes + m] = RK21 * R[3] + RK22 * R[4] + RK23 * R[5];
                element.K[(1 * nodes + n) * nodes * 3 + 2 * nodes + m] = RK21 * R[6] + RK22 * R[7] + RK23 * R[8];

                element.K[(2 * nodes + n) * nodes * 3 + 0 * nodes + m] = RK31 * R[0] + RK32 * R[1] + RK33 * R[2];
                element.K[(2 * nodes + n) * nodes * 3 + 1 * nodes + m] = RK31 * R[3] + RK32 * R[4] + RK33 * R[5];
                element.K[(2 * nodes + n) * nodes * 3 + 2 * nodes + m] = RK31 * R[6] + RK32 * R[7] + RK33 * R[8];
            }
        }
    }
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_ */
