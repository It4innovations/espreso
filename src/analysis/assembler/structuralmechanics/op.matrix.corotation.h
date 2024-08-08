
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
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
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
    }

    double cw, cN[nodes], cdN[nodes][2];

    template <typename Element>
    void simd(Element &element)
    {

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
        // S = dND[i]
        // [  0 -3  2 ]
        // [  3  0 -1 ]
        // [ -2  1  0 ]

        std::vector<int> perm(3 * nodes);
        for (size_t i = 0; i < nodes; ++i) {
            for (size_t d = 0; d < 3UL; ++d) {
                perm[i * 3 + d] = d * nodes + i;
            }
        }

        {
            SIMD _KK[3 * nodes * 3 * nodes];
            for (size_t r = 0; r < 3 * nodes; ++r) {
                for (size_t c = 0; c < 3 * nodes; ++c) {
                    _KK[r * 3 * nodes + c] = element.K[perm[r] * 3 * nodes + perm[c]];
                }
            }
            for (size_t r = 0; r < 3 * nodes; ++r) {
                for (size_t c = 0; c < 3 * nodes; ++c) {
                    element.K[r * 3 * nodes + c] = _KK[r * 3 * nodes + c];
                }
            }
        }

//        printf("K\n"); print(3 * nodes, 3 * nodes, element.K);

        SIMD JC[9];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    JC[i * 3 + j] = JC[i * 3 + j] + load1(cdN[n][i]) * element.coords.node[n][j];
                }
            }
        }

        SIMD det, invJC[9]; inv(JC, det, invJC);

        SIMD cdND[nodes][3];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(cdN[n][0]);
            SIMD dNY = load1(cdN[n][1]);
            SIMD dNZ = load1(cdN[n][2]);
            cdND[n][0] = invJC[0] * dNX + invJC[1] * dNY + invJC[2] * dNZ;
            cdND[n][1] = invJC[3] * dNX + invJC[4] * dNY + invJC[5] * dNZ;
            cdND[n][2] = invJC[6] * dNX + invJC[7] * dNY + invJC[8] * dNZ;
        }

        SIMD F[9];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (size_t n = 0; n < nodes; ++n) {
                    F[j * 3 + i] = F[j * 3 + i] + element.displacement[n][j] * cdND[n][i]; // transpose
                }
            }
        }
        F[0] = F[0] + load1(1.0);
        F[4] = F[4] + load1(1.0);
        F[8] = F[8] + load1(1.0);

//        printf("F\n"); print(3, 3, F);

        // R = F * (F' * F)^(-1/2)
        SIMD FF[9]; multAtB<3, 3, 3>(FF, F, F);
        SIMD eigVal[3], eigVec[9];
        eigSym(FF, eigVal, eigVec);
        printf("EIG\n"); print(1, 3, eigVal);
        SIMD FD[9] = { sqrt(eigVal[0]), zeros(), zeros(), zeros(), sqrt(eigVal[1]), zeros(), zeros(), zeros(), sqrt(eigVal[2]) };
        set<3, 3>(FF, zeros());
        multAtBA<3, 3>(FF, eigVec, FD);
        SIMD invFF[9]; inv(FF, invFF);
        SIMD R[9]; multAB<3, 3, 3>(R, F, invFF);

//        printf("R\n"); print(3, 3, R);

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
            for (int i = 0; i < 9; ++i) {
                S_bar[9 * n + i] = S_barn[i];
            }

            // u_lie = Re' * vektor od teziste (deformovane) - vektor od teziste (nedef.)
            SIMD u_lien[3]; multAtB<3, 3, 1>(u_lien, R, ldisp);
            u_lie[3 * n + 0] = u_lien[0] - lorigin[0];
            u_lie[3 * n + 1] = u_lien[1] - lorigin[1];
            u_lie[3 * n + 2] = u_lien[2] - lorigin[2];
        }

        SIMD A_transp[3 * 3 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            A_transp[0 * 3 * nodes + 3 * n + 1] = -cdND[n][2];
            A_transp[0 * 3 * nodes + 3 * n + 2] =  cdND[n][1];
            A_transp[1 * 3 * nodes + 3 * n + 0] =  cdND[n][2];
            A_transp[1 * 3 * nodes + 3 * n + 2] = -cdND[n][0];
            A_transp[2 * 3 * nodes + 3 * n + 0] = -cdND[n][1];
            A_transp[2 * 3 * nodes + 3 * n + 1] =  cdND[n][0];
        }

        // internal force computation
        SIMD AS[9]; multAB<3, 3 * nodes, 3>(AS, A_transp, S_bar);
        SIMD invAS[9]; inv(AS, invAS);
        // G_bar = inv(A_trans * S_bar) * A_transp
        SIMD G_bar[3 * 3 * nodes]; multAB<3, 3, 3 * nodes>(G_bar, invAS, A_transp);
        // P_bar = eye(24) - S_bar * G_bar
        SIMD P_bar[3 * nodes * 3 * nodes];
        for (size_t i = 0; i < 3 * nodes; ++i) {
            P_bar[i * 3 * nodes + i] = load1(1.);
        }
        multAB<3 * nodes, 3, 3 * nodes>(P_bar, S_bar, G_bar, load1(-1.));

        SIMD RR[3 * nodes * 3 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            for (int nn = 0; nn < 3; ++nn) {
                for (int mm = 0; mm < 3; ++mm) {
                    RR[(n * 3 + nn) * 3 * nodes + 3 * n + mm] = R[nn * 3 + mm];
                }
            }
        }

        // f_ile = Ke * u_lie
        SIMD f_ile[3 * nodes]; multAB<3 * nodes, 3 * nodes, 1>(f_ile, element.K, u_lie);
        SIMD fiep[3 * nodes]; multAtB<3 * nodes, 3 * nodes, 1>(fiep, P_bar, f_ile);
        multAB<3 * nodes, 3 * nodes, 1>(element.nf, RR, fiep, load1(-1));
//        printf("nf\n"); print(3 * nodes, 1, element.nf);
//        for (size_t n = 0; n < nodes; ++n) {
//            SIMD fiepn[3]; fiepn[0] = fiep[3 * n + 0]; fiepn[1] = fiep[3 * n + 1]; fiepn[2] = fiep[3 * n + 2];
//            SIMD fie[3]; multAB<3, 3, 1>(fie, R, fiepn);
//            element.nf[0 * nodes + n] = element.nf[3 * n + 0] - fie[0];
//            element.nf[1 * nodes + n] = element.nf[3 * n + 1] - fie[1];
//            element.nf[2 * nodes + n] = element.nf[3 * n + 2] - fie[2];
//        }

        SIMD Fn_bar[3 * nodes * 3];
        for (size_t n = 0; n < nodes; ++n) {
            Fn_bar[(3 * n + 0) * 3 + 1] = -fiep[n * 3 + 2];
            Fn_bar[(3 * n + 0) * 3 + 2] =  fiep[n * 3 + 1];
            Fn_bar[(3 * n + 1) * 3 + 0] =  fiep[n * 3 + 2];
            Fn_bar[(3 * n + 1) * 3 + 2] = -fiep[n * 3 + 0];
            Fn_bar[(3 * n + 2) * 3 + 0] = -fiep[n * 3 + 1];
            Fn_bar[(3 * n + 2) * 3 + 1] =  fiep[n * 3 + 0];
        }

//        printf("Fn_bar\n"); print(3 * nodes, 3, Fn_bar);

        // KM = Re * P_bar' * Ke * P_bar * Re'
        // KGR = -Re * F_bar * G_bar * Re'
        // KGP = -Re * G_bar' * Fn_bar' P_bar * Re'
        // Kte = KM + KGR + KGP

//        printf("RR\n"); print(3 * nodes, 3 * nodes, RR);

        // K = Re * ((P_bar' * Ke - G_bar' * Fn_bar') * P_bar - Fn_bar * G_bar) * Re'
        SIMD Ke[3 * nodes * 3 * nodes];
        multAtB <3 * nodes, 3 * nodes, 3 * nodes>(Ke, P_bar, element.K);
//        printf("1\n"); print(3 * nodes, 3 * nodes, Ke);
        multAtBt<3 * nodes, 3        , 3 * nodes>(Ke, G_bar, Fn_bar, load1(-1));
//        printf("2\n"); print(3 * nodes, 3 * nodes, Ke);
        set<3 * nodes, 3 * nodes>(element.K, load1(0));
        multAB  <3 * nodes, 3 * nodes, 3 * nodes>(element.K, Ke, P_bar);
//        printf("3\n"); print(3 * nodes, 3 * nodes, element.K);
        multAB  <3 * nodes, 3         , 3 * nodes>(element.K, Fn_bar, G_bar, load1(-1));
//        printf("4\n"); print(3 * nodes, 3 * nodes, element.K);

        set<3 * nodes, 3 * nodes>(Ke, load1(0));
        multABt <3 * nodes, 3 * nodes, 3 * nodes>(Ke, element.K, RR);
//        printf("5\n"); print(3 * nodes, 3 * nodes, Ke);
        set<3 * nodes, 3 * nodes>(element.K, load1(0));
        multAB  <3 * nodes, 3 * nodes, 3 * nodes>(element.K, RR, Ke);
//        printf("6\n"); print(3 * nodes, 3 * nodes, element.K);

//        printf("K\n"); print(3 * nodes, 3 * nodes, element.K);

        {
            SIMD _KK[3 * nodes * 3 * nodes];
            for (size_t r = 0; r < 3 * nodes; ++r) {
                for (size_t c = 0; c < 3 * nodes; ++c) {
                    _KK[perm[r] * 3 * nodes + perm[c]] = element.K[r * 3 * nodes + c];
                }
            }
            for (size_t r = 0; r < 3 * nodes; ++r) {
                for (size_t c = 0; c < 3 * nodes; ++c) {
                    element.K[r * 3 * nodes + c] = _KK[r * 3 * nodes + c];
                }
            }
        }

        printf("FINAL\n"); print(3 * nodes, 3 * nodes, element.K);
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_COROTATION_H_ */
