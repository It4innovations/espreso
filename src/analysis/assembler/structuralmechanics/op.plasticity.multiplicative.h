
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_

#include "basis/utilities/utils.h"
#include "analysis/assembler/general/element.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "mesh/store/nameddata.h"
#include "esinfo/meshinfo.h"

#include <memory>
#include <cmath>

namespace espreso {

struct PlasticityMultiplicative: SubKernel {
    const char* name() const { return "PlasticityMultiplicative"; }

    PlasticityMultiplicative()
    : behaviour(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN), configuration(nullptr),
      elements(0)
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }

    void activate(size_t interval, StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour, const PlasticityPropertiesConfiguration *configuration)
    {
        this->behaviour = behaviour;
        this->configuration = configuration;
        this->isconst = false;
        this->isactive = true;
    }

    void init(size_t chunks, size_t gps)
    {
        invCp.resize(chunks * gps * SIMD::size * 6 + 3 * SIMD::size);
        alpha.resize(chunks * gps * SIMD::size * 1 + 3 * SIMD::size);
        double *_invCp = utils::getAligned(SIMD::size, invCp);
        for (size_t c = 0; c < chunks; ++c) {
            for (size_t i = 0; i < gps; ++i) {
                for (size_t s = 0; s < SIMD::size; ++s) { // eye in Voigh notation
                    _invCp[c * gps * SIMD::size * 6 + i * SIMD::size * 6 + 0 * SIMD::size + s] = 1;
                    _invCp[c * gps * SIMD::size * 6 + i * SIMD::size * 6 + 1 * SIMD::size + s] = 1;
                    _invCp[c * gps * SIMD::size * 6 + i * SIMD::size * 6 + 2 * SIMD::size + s] = 1;
                }
            }
        }
    }

    StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour;
    const PlasticityPropertiesConfiguration *configuration;

    size_t elements;
    std::vector<double> invCp, alpha;
};

struct PlasticityMultiplicativeStorage: PlasticityMultiplicative {

    PlasticityMultiplicativeStorage(PlasticityMultiplicative &base, SubKernel::Action action)
    : PlasticityMultiplicative(base),
      invCp(base.invCp.data()),
      alpha(base.alpha.data()),
      save(action == SubKernel::SOLUTION)
    {
        if (isactive) {
            invCp = utils::getAligned(SIMD::size, base.invCp);
            alpha = utils::getAligned(SIMD::size, base.alpha);
        }
    }

    double *invCp, *alpha;
    bool save;
};

template <size_t nodes, size_t gps, size_t ndim> struct PlasticityMultiplicativeKernel: PlasticityMultiplicativeStorage {
    PlasticityMultiplicativeKernel(PlasticityMultiplicative &base, SubKernel::Action action): PlasticityMultiplicativeStorage(base, action) {}

    template <typename Element>
    void simd(Element &element)
    {

    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }
};

template <size_t nodes, size_t gps> struct PlasticityMultiplicativeKernel<nodes, gps, 3>: PlasticityMultiplicativeStorage {
    PlasticityMultiplicativeKernel(PlasticityMultiplicative &base, SubKernel::Action action): PlasticityMultiplicativeStorage(base, action) {}

    SIMD logJbarDivJbar;

    template <typename Element>
    void simd(Element &element)
    {
        SIMD ve, Ve, divJbar, mean_dNdx[3 * nodes];
        for (size_t gp = 0; gp < gps; ++gp) {
            SIMD JX[9];
            for (size_t n = 0; n < nodes; ++n) {
                SIMD coordsX = element.coords.node[n][0];
                SIMD coordsY = element.coords.node[n][1];
                SIMD coordsZ = element.coords.node[n][2];
                SIMD dNX = load1(element.dN[gp][n][0]);
                SIMD dNY = load1(element.dN[gp][n][1]);
                SIMD dNZ = load1(element.dN[gp][n][2]);
                JX[0] = JX[0] + dNX * coordsX;
                JX[3] = JX[3] + dNX * coordsY;
                JX[6] = JX[6] + dNX * coordsZ;
                JX[1] = JX[1] + dNY * coordsX;
                JX[4] = JX[4] + dNY * coordsY;
                JX[7] = JX[7] + dNY * coordsZ;
                JX[2] = JX[2] + dNZ * coordsX;
                JX[5] = JX[5] + dNZ * coordsY;
                JX[8] = JX[8] + dNZ * coordsZ;
            }
            Ve = Ve + determinant33(JX);

            SIMD Jx[9];
            for (size_t n = 0; n < nodes; ++n) {
                SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
                SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
                SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
//                printf("dx %+.10e %+.10e\n", element.displacement[n][0][0], element.displacement[n][0][1]);
//                printf("dy %+.10e %+.10e\n", element.displacement[n][1][0], element.displacement[n][1][1]);
//                printf("dz %+.10e %+.10e\n", element.displacement[n][2][0], element.displacement[n][2][1]);
                SIMD dNX = load1(element.dN[gp][n][0]);
                SIMD dNY = load1(element.dN[gp][n][1]);
                SIMD dNZ = load1(element.dN[gp][n][2]);
                Jx[0] = Jx[0] + dNX * coordsX;
                Jx[3] = Jx[3] + dNX * coordsY;
                Jx[6] = Jx[6] + dNX * coordsZ;
                Jx[1] = Jx[1] + dNY * coordsX;
                Jx[4] = Jx[4] + dNY * coordsY;
                Jx[7] = Jx[7] + dNY * coordsZ;
                Jx[2] = Jx[2] + dNZ * coordsX;
                Jx[5] = Jx[5] + dNZ * coordsY;
                Jx[8] = Jx[8] + dNZ * coordsZ;
            }
            SIMD detJx, invJx[9];
            inv33(Jx, detJx, invJx);
//            printf("detJx %+.10e %+.10e\n", detJx[0], detJx[1]);
            ve = ve + detJx;

            for (size_t n = 0; n < nodes; ++n) {
                SIMD dNX = load1(element.dN[gp][n][0]);
                SIMD dNY = load1(element.dN[gp][n][1]);
                SIMD dNZ = load1(element.dN[gp][n][2]);
                mean_dNdx[0 * nodes + n] = mean_dNdx[0 * nodes + n] + detJx * (invJx[0] * dNX + invJx[3] * dNY + invJx[6] * dNZ);
                mean_dNdx[1 * nodes + n] = mean_dNdx[1 * nodes + n] + detJx * (invJx[1] * dNX + invJx[4] * dNY + invJx[7] * dNZ);
                mean_dNdx[2 * nodes + n] = mean_dNdx[2 * nodes + n] + detJx * (invJx[2] * dNX + invJx[5] * dNY + invJx[8] * dNZ);
            }
        }

//        printf("MEAN\n"); print(3, nodes, mean_dNdx);

//        printf("ve %+.10e %+.10e\n", ve[0], ve[1]);
//        printf("Ve %+.10e %+.10e\n", Ve[0], Ve[1]);
        SIMD Jbar = ve / Ve;
        divJbar = load1(1) / Jbar;
        logJbarDivJbar = log(Jbar) * divJbar;

        SIMD bulk_modulus = element.ecf.youngModulus[0] / (load1(3) - load1(6) * element.ecf.poissonRatio[0]);
        SIMD kappabar = bulk_modulus * divJbar - bulk_modulus * logJbarDivJbar;

        multABt<3 * nodes, 1, 3 * nodes>(element.K, mean_dNdx, mean_dNdx, kappabar / ve);
    }

    void n4_otimes_symallcomb(SIMD cc1[3], SIMD cc3[3], SIMD n[9], SIMD vc4m_hat[6*6], size_t gp)
    {
        for (int a = 0; a < 3; ++a) {
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[0 * 3 + a] * n[0 * 3 + a];
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc1[a] * n[0 * 3 + a] * n[0 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[a] * n[0 * 3 + a] * n[1 * 3 + a] * n[0 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[a] * n[0 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[a] * n[0 * 3 + a] * n[1 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[a] * n[1 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[a] * n[1 * 3 + a] * n[2 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[a] * n[0 * 3 + a] * n[2 * 3 + a] * n[0 * 3 + a] * n[2 * 3 + a];
        }
        for (int k = 3, a, b; k < 6; ++k) {
            switch (k) {
            case 3: a = 0; b = 1; break;
            case 4: a = 1; b = 2; break;
            case 5: a = 0; b = 2; break;
            }
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + a] *  n[0 * 3 + b] * n[0 * 3 + b] * load1(2));
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + a] *  n[1 * 3 + b] * n[1 * 3 + b] + n[0 * 3 + b] * n[0 * 3 + b] * n[1 * 3 + a] * n[1 * 3 + a]);
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + a] *  n[2 * 3 + b] * n[2 * 3 + b] + n[0 * 3 + b] * n[0 * 3 + b] * n[2 * 3 + a] * n[2 * 3 + a]);
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[0 * 3 + a] * n[1 * 3 + b] + n[0 * 3 + b] * n[1 * 3 + a]));
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + a] *  n[1 * 3 + b] * n[2 * 3 + b] + n[0 * 3 + b] * n[0 * 3 + b] * n[1 * 3 + a] * n[2 * 3 + a]);
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]));
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + a] *  n[1 * 3 + b] * n[1 * 3 + b] * load1(2));
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + a] *  n[2 * 3 + b] * n[2 * 3 + b] + n[1 * 3 + b] * n[1 * 3 + b] * n[2 * 3 + a] * n[2 * 3 + a]);
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + b] * n[1 * 3 + a] + n[1 * 3 + b] * n[0 * 3 + a]));
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]));
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + a] *  n[0 * 3 + b] * n[2 * 3 + b] + n[1 * 3 + b] * n[1 * 3 + b] * n[0 * 3 + a] * n[2 * 3 + a]);
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + a] *  n[2 * 3 + b] * n[2 * 3 + b] * load1(2));
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + a] *  n[0 * 3 + b] * n[1 * 3 + b] + n[2 * 3 + b] * n[2 * 3 + b] * n[0 * 3 + a] * n[1 * 3 + a]);
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + b] * n[2 * 3 + a] + n[2 * 3 + b] * n[1 * 3 + a]));
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[0 * 3 + b] * n[2 * 3 + a] + n[2 * 3 + b] * n[0 * 3 + a]));
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + b] *  n[1 * 3 + a] * n[1 * 3 + b] * load1(2));
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]));
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[k] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]));
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[k] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[2 * 3 + a] * n[2 * 3 + b] + n[2 * 3 + b] * n[2 * 3 + a]));
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[0 * 3 + b] + n[1 * 3 + b] * n[0 * 3 + a]));
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[k] * (n[2 * 3 + a] * n[2 * 3 + b] *  n[0 * 3 + a] * n[0 * 3 + b] * load1(2));
        }
        for (int k = 6, a, b; k < 9; ++k) {
            switch (k) {
            case 6: a = 0; b = 1; break;
            case 7: a = 1; b = 2; break;
            case 8: a = 0; b = 2; break;
            }
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] *  n[0 * 3 + a] * n[0 * 3 + b] * load1(4));
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] *  n[1 * 3 + a] * n[1 * 3 + b] * load1(4));
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] *  n[2 * 3 + a] * n[2 * 3 + b] * load1(4));
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[0 * 3 + a] * n[1 * 3 + b] + n[0 * 3 + b] * n[1 * 3 + a]) * load1(2));
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc3[k - 6] * (n[0 * 3 + a] * n[0 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc3[k - 6] * (n[1 * 3 + a] * n[1 * 3 + b] *  n[1 * 3 + a] * n[1 * 3 + b] * load1(4));
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc3[k - 6] * (n[1 * 3 + a] * n[1 * 3 + b] *  n[2 * 3 + a] * n[2 * 3 + b] * load1(4));
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc3[k - 6] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + a] * n[1 * 3 + b] + n[0 * 3 + b] * n[1 * 3 + a]) * load1(2));
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc3[k - 6] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc3[k - 6] * (n[1 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc3[k - 6] * (n[2 * 3 + a] * n[2 * 3 + b] *  n[2 * 3 + a] * n[2 * 3 + b] * load1(4));
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc3[k - 6] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[0 * 3 + a] * n[1 * 3 + b] + n[0 * 3 + b] * n[1 * 3 + a]) * load1(2));
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc3[k - 6] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc3[k - 6] * (n[2 * 3 + a] * n[2 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) * load1(2));
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc3[k - 6] * (n[0 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + a] * n[1 * 3 + b] + n[0 * 3 + b] * n[1 * 3 + a]) + n[0 * 3 + b] * n[1 * 3 + a] * (n[0 * 3 + b] * n[1 * 3 + a] + n[0 * 3 + a] * n[1 * 3 + b]));
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc3[k - 6] * (n[0 * 3 + a] * n[1 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) + n[0 * 3 + b] * n[1 * 3 + a] * (n[1 * 3 + b] * n[2 * 3 + a] + n[1 * 3 + a] * n[2 * 3 + b]));
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc3[k - 6] * (n[0 * 3 + a] * n[1 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) + n[0 * 3 + b] * n[1 * 3 + a] * (n[0 * 3 + b] * n[2 * 3 + a] + n[0 * 3 + a] * n[2 * 3 + b]));
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc3[k - 6] * (n[1 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) + n[1 * 3 + b] * n[2 * 3 + a] * (n[1 * 3 + b] * n[2 * 3 + a] + n[1 * 3 + a] * n[2 * 3 + b]));
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc3[k - 6] * (n[1 * 3 + a] * n[2 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) + n[1 * 3 + b] * n[2 * 3 + a] * (n[0 * 3 + b] * n[2 * 3 + a] + n[0 * 3 + a] * n[2 * 3 + b]));
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc3[k - 6] * (n[0 * 3 + a] * n[2 * 3 + b] * (n[0 * 3 + a] * n[2 * 3 + b] + n[0 * 3 + b] * n[2 * 3 + a]) + n[0 * 3 + b] * n[2 * 3 + a] * (n[0 * 3 + b] * n[2 * 3 + a] + n[0 * 3 + a] * n[2 * 3 + b]));
        }
        vc4m_hat[1 * 6 + 0] = vc4m_hat[0 * 6 + 1];
        vc4m_hat[2 * 6 + 0] = vc4m_hat[0 * 6 + 2];
        vc4m_hat[2 * 6 + 1] = vc4m_hat[1 * 6 + 2];
        vc4m_hat[3 * 6 + 0] = vc4m_hat[0 * 6 + 3];
        vc4m_hat[3 * 6 + 1] = vc4m_hat[1 * 6 + 3];
        vc4m_hat[3 * 6 + 2] = vc4m_hat[2 * 6 + 3];
        vc4m_hat[4 * 6 + 0] = vc4m_hat[0 * 6 + 4];
        vc4m_hat[4 * 6 + 1] = vc4m_hat[1 * 6 + 4];
        vc4m_hat[4 * 6 + 2] = vc4m_hat[2 * 6 + 4];
        vc4m_hat[4 * 6 + 3] = vc4m_hat[3 * 6 + 4];
        vc4m_hat[5 * 6 + 0] = vc4m_hat[0 * 6 + 5];
        vc4m_hat[5 * 6 + 1] = vc4m_hat[1 * 6 + 5];
        vc4m_hat[5 * 6 + 2] = vc4m_hat[2 * 6 + 5];
        vc4m_hat[5 * 6 + 3] = vc4m_hat[3 * 6 + 5];
        vc4m_hat[5 * 6 + 4] = vc4m_hat[4 * 6 + 5];
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD sigmaY0 = load1(1e8);
        SIMD Hisotropic = load1(2.4e9);
        SIMD mu = element.ecf.youngModulus[0] / (load1(2.) + load1(2.) * element.ecf.poissonRatio[0]);
        SIMD bulk_modulus = element.ecf.youngModulus[0] / (load1(3) - load1(6) * element.ecf.poissonRatio[0]);

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                element.F[i * 3 + j] = zeros();
                for (size_t n = 0; n < nodes; ++n) {
                    element.F[i * 3 + j] = element.F[i * 3 + j] + element.displacement[n][i] * element.dND[n * 3 + j];
                }
            }
        }
        element.F[0] = element.F[0] + load1(1.);
        element.F[4] = element.F[4] + load1(1.);
        element.F[8] = element.F[8] + load1(1.);

//        if (gp == 0) { printf("F\n"); print(3, 3, element.F); }

        SIMD detF, invF[9];
        inv33(element.F, detF, invF);
//        if (gp == 0) printf("det F %.10e\n", detF[0]);

        SIMD alpha = load(this->alpha);
        SIMD invCp[9] = {
                load(this->invCp + 0 * SIMD::size), load(this->invCp + 3 * SIMD::size), load(this->invCp + 5 * SIMD::size),
                load(this->invCp + 3 * SIMD::size), load(this->invCp + 1 * SIMD::size), load(this->invCp + 4 * SIMD::size),
                load(this->invCp + 5 * SIMD::size), load(this->invCp + 4 * SIMD::size), load(this->invCp + 2 * SIMD::size)
        };

//        if (gp != 100UL) { printf("alpha\n"); print(1, 1, &alpha); }
//        if (gp != 100UL) { printf("invCp\n"); print(3, 3, invCp); }

        SIMD C2[9]; multAtB<3, 3, 3>(C2, element.F, element.F);
        SIMD detC2, invC2[9];
        inv33(C2, detC2, invC2);

        // 0 1 2   0 3 5
        // 3 4 5     1 4
        // 6 7 8       2

        //  0  1  2  3  4  5
        //  6  7  8  9 10 11
        // 12 13 14 15 16 17
        // 18 19 20 21 22 23
        // 24 25 26 27 28 29
        // 30 31 32 33 34 35

        SIMD pbar = bulk_modulus * logJbarDivJbar;
        SIMD vC4[6*6];
        vC4[ 0] =           pbar * detF * (invC2[0] * invC2[0] -            load1(2) * invC2[0] * invC2[0]);
        vC4[ 1] = vC4[ 6] = pbar * detF * (invC2[0] * invC2[4] -            load1(2) * invC2[1] * invC2[1]);
        vC4[ 2] = vC4[12] = pbar * detF * (invC2[0] * invC2[8] -            load1(2) * invC2[2] * invC2[2]);
        vC4[ 3] = vC4[18] = pbar * detF * (invC2[0] * invC2[1] -            load1(2) * invC2[0] * invC2[1]);
        vC4[ 4] = vC4[24] = pbar * detF * (invC2[0] * invC2[5] -            load1(2) * invC2[1] * invC2[2]);
        vC4[ 5] = vC4[30] = pbar * detF * (invC2[0] * invC2[2] -            load1(2) * invC2[0] * invC2[2]);
        vC4[ 7] =           pbar * detF * (invC2[4] * invC2[4] -            load1(2) * invC2[4] * invC2[4]);
        vC4[ 8] = vC4[13] = pbar * detF * (invC2[4] * invC2[8] -            load1(2) * invC2[5] * invC2[5]);
        vC4[ 9] = vC4[19] = pbar * detF * (invC2[4] * invC2[1] -            load1(2) * invC2[4] * invC2[1]);
        vC4[10] = vC4[25] = pbar * detF * (invC2[4] * invC2[5] -            load1(2) * invC2[4] * invC2[5]);
        vC4[11] = vC4[31] = pbar * detF * (invC2[4] * invC2[2] -            load1(2) * invC2[1] * invC2[5]);
        vC4[14] =           pbar * detF * (invC2[8] * invC2[8] -            load1(2) * invC2[8] * invC2[8]);
        vC4[15] = vC4[20] = pbar * detF * (invC2[8] * invC2[1] -            load1(2) * invC2[2] * invC2[5]);
        vC4[16] = vC4[26] = pbar * detF * (invC2[8] * invC2[5] -            load1(2) * invC2[8] * invC2[5]);
        vC4[17] = vC4[32] = pbar * detF * (invC2[8] * invC2[2] -            load1(2) * invC2[8] * invC2[2]);
        vC4[21] =           pbar * detF * (invC2[1] * invC2[1] - invC2[1] * invC2[1] - invC2[0] * invC2[4]);
        vC4[22] = vC4[27] = pbar * detF * (invC2[1] * invC2[5] - invC2[1] * invC2[5] - invC2[2] * invC2[4]);
        vC4[23] = vC4[33] = pbar * detF * (invC2[1] * invC2[2] - invC2[1] * invC2[2] - invC2[0] * invC2[5]);
        vC4[28] =           pbar * detF * (invC2[5] * invC2[5] - invC2[5] * invC2[5] - invC2[4] * invC2[8]);
        vC4[29] = vC4[34] = pbar * detF * (invC2[5] * invC2[2] - invC2[2] * invC2[5] - invC2[1] * invC2[8]);
        vC4[35] =           pbar * detF * (invC2[2] * invC2[2] - invC2[2] * invC2[2] - invC2[0] * invC2[8]);

//        if (gp == 0) { printf("vC4\n"); print(6, 6, vC4); }

        SIMD be_trial[9]; multABAt<3, 3>(be_trial, element.F, invCp);
        SIMD eigVal[3], n_trial[9], N_trial[9];
        eigSym33B(be_trial, eigVal, n_trial);
        SIMD lambda_trial[3] = { sqrt(eigVal[0]), sqrt(eigVal[1]), sqrt(eigVal[2]) };
//        printf("%+.6e %+.6e %+.6e %+.6e %+.6e %+.6e\n", lambda_trial[0][0], lambda_trial[1][0], lambda_trial[2][0], lambda_trial[0][1], lambda_trial[1][1], lambda_trial[2][1]);
        SIMD N_tmp[9] = {
                lambda_trial[0] * n_trial[0], lambda_trial[1] * n_trial[3], lambda_trial[2] * n_trial[6],
                lambda_trial[0] * n_trial[1], lambda_trial[1] * n_trial[4], lambda_trial[2] * n_trial[7],
                lambda_trial[0] * n_trial[2], lambda_trial[1] * n_trial[5], lambda_trial[2] * n_trial[8]
        };
//        if (gp == 0) { printf("N_tmp\n"); print(3, 3, N_tmp); }
        multAB<3, 3, 3>(N_trial, invF, N_tmp);

        SIMD tau_trial_devdiag[3] = {
                load1(2) * mu * log(lambda_trial[0]) - load1(2/3.) * mu * log(detF),
                load1(2) * mu * log(lambda_trial[1]) - load1(2/3.) * mu * log(detF),
                load1(2) * mu * log(lambda_trial[2]) - load1(2/3.) * mu * log(detF)
        };
        SIMD norm_tau_trial_dev = sqrt(tau_trial_devdiag[0] * tau_trial_devdiag[0] + tau_trial_devdiag[1] * tau_trial_devdiag[1] + tau_trial_devdiag[2] * tau_trial_devdiag[2]);

        SIMD f_trial = load1(std::sqrt(3/2.)) * norm_tau_trial_dev - (sigmaY0 * detF + Hisotropic * alpha);

//        if (gp == 0) { printf("f_trial\n"); print(1, 1, &f_trial); }

        double eye3_13[9] = { 2/3., -1/3., -1/3., -1/3., 2/3., -1/3., -1/3., -1/3., 2/3. };

        double sqrt23 = sqrt(2/3.);
        SIMD nu[3], sigma_devdiag[3], cc1[6], cc3[3];
        SIMD be[9];
        for (size_t s = 0; s < SIMD::size; ++s) {
            auto k2ij = [] (int k, int &i, int &j) {
                switch (k) {
                case 0: i = 0; j = 0; break;
                case 1: i = 1; j = 1; break;
                case 2: i = 2; j = 2; break;
                case 3: i = 0; j = 1; break;
                case 4: i = 1; j = 2; break;
                case 5: i = 0; j = 2; break;
                }
            };

            if (f_trial[s] > 0) {
                double Dgamma = f_trial[s] / (3 * mu[s] + Hisotropic[s]);
                for (int k = 0, i, j; k < 6; ++k) {
                    k2ij(k, i, j);
                    nu[i][s] = tau_trial_devdiag[i][s] / (sqrt23 * norm_tau_trial_dev[s]);
                    sigma_devdiag[i][s] = ((1 - (2 * mu[s] * Dgamma) / (sqrt23 * norm_tau_trial_dev[s])) * tau_trial_devdiag[i][s]) / detF[s];
                    cc1[k][s] = (1 / detF[s]) * (2 * mu[s] * (1 - 2 * mu[s] * Dgamma / (sqrt23 * norm_tau_trial_dev[s])) * eye3_13[i * 3 + j] - 4 * mu[s] * mu[s] * nu[i][s] * nu[j][s] * (1. / (3 * mu[s] + Hisotropic[s]) - sqrt23 * Dgamma / norm_tau_trial_dev[s]));
                    if (k < 3) {
                        cc1[k][s] -= 2 * sigma_devdiag[i][s];
                    }
                }
                if (save) {
                    this->alpha[s] += Dgamma;
                    double lambda_e[3] = { 0, 0, 0 };
                    for (int k = 0; k < 3; ++k) {
                        lambda_e[k] = std::exp(std::log(lambda_trial[k][s]) - Dgamma * nu[k][s]);
                    }
                    for (int i = 0; i < 3; ++i) {
                        for (int j = i; j < 3; ++j) {
                            be[i * 3 + j][s] += lambda_e[0] * lambda_e[0] * n_trial[0 * 3 + i][s] * n_trial[0 * 3 + j][s];
                            be[i * 3 + j][s] += lambda_e[1] * lambda_e[1] * n_trial[1 * 3 + i][s] * n_trial[1 * 3 + j][s];
                            be[i * 3 + j][s] += lambda_e[2] * lambda_e[2] * n_trial[2 * 3 + i][s] * n_trial[2 * 3 + j][s];
                        }
                    }
                }
            } else {
                for (int k = 0, i, j; k < 6; ++k) {
                    k2ij(k, i, j);
                    cc1[k][s] = (1 / detF[s]) * 2 * mu[s] * eye3_13[i * 3 + j];
                    if (k < 3) {
                        sigma_devdiag[i][s] = tau_trial_devdiag[i][s] / detF[s];
                        cc1[k][s] -= 2 * sigma_devdiag[i][s];
                    }
                }
            }

            for (int k = 3, i, j; k < 6; ++k) {
                k2ij(k, i, j);
                if (std::fabs(lambda_trial[i][s] - lambda_trial[j][s]) <= 1e-5 * std::fabs(lambda_trial[j][s])) {
//                    printf("up ");
                    cc3[k - 3][s] = mu[s] * (1 / detF[s]) - sigma_devdiag[i][s];
                } else {
//                    printf("down ");
                    cc3[k - 3][s] = (sigma_devdiag[i][s] * lambda_trial[j][s] * lambda_trial[j][s] - sigma_devdiag[j][s] * lambda_trial[i][s] * lambda_trial[i][s]) / (lambda_trial[i][s] * lambda_trial[i][s] - lambda_trial[j][s] * lambda_trial[j][s]);
                }
            }
        }
        if (save) {
            be[3] = be[1]; be[6] = be[2]; be[7] = be[5]; // make be symmetric
            set<3, 3>(invCp, load1(0));
            multABAt<3, 3>(invCp, invF, be);
            for (size_t s = 0; s < SIMD::size; ++s) {
                if (f_trial[s] > 0) {
                    this->invCp[0 * SIMD::size + s] = invCp[0][s];
                    this->invCp[1 * SIMD::size + s] = invCp[4][s];
                    this->invCp[2 * SIMD::size + s] = invCp[8][s];
                    this->invCp[3 * SIMD::size + s] = invCp[1][s];
                    this->invCp[4 * SIMD::size + s] = invCp[5][s];
                    this->invCp[5 * SIMD::size + s] = invCp[2][s];
                }
            }
        }
//        printf("\n");

        if (std::isnan(cc1[0][1])) {
            printf("cc1\n"); print(1, 6, cc1);
        }
//        if (gp == 0) { printf("cc1\n"); print(1, 6, cc1); }
//        if (gp == 0) { printf("cc3\n"); print(1, 3, cc3); }

        SIMD sigma_diag[3] = { sigma_devdiag[0] + pbar, sigma_devdiag[1] + pbar, sigma_devdiag[2] + pbar };
        SIMD sigma[9];
        multABt<3, 1, 3>(sigma, n_trial + 0, n_trial + 0, sigma_diag[0]);
        multABt<3, 1, 3>(sigma, n_trial + 3, n_trial + 3, sigma_diag[1]);
        multABt<3, 1, 3>(sigma, n_trial + 6, n_trial + 6, sigma_diag[2]);

        for (int i = 0; i < 6; ++i) {
            element.vS[i] = zeros();
        }
        for (int a = 0; a < 3; ++a) {
            element.vS[0] = element.vS[0] + sigma_diag[a] * N_trial[0 * 3 + a] * N_trial[0 * 3 + a];
            element.vS[1] = element.vS[1] + sigma_diag[a] * N_trial[1 * 3 + a] * N_trial[1 * 3 + a];
            element.vS[2] = element.vS[2] + sigma_diag[a] * N_trial[2 * 3 + a] * N_trial[2 * 3 + a];
            element.vS[3] = element.vS[3] + sigma_diag[a] * N_trial[0 * 3 + a] * N_trial[1 * 3 + a];
            element.vS[4] = element.vS[4] + sigma_diag[a] * N_trial[1 * 3 + a] * N_trial[2 * 3 + a];
            element.vS[5] = element.vS[5] + sigma_diag[a] * N_trial[0 * 3 + a] * N_trial[2 * 3 + a];
        }

//        if (gp == 0) { printf("S2\n"); print(1, 6, element.vS); }
//        if (gp == 0) { printf("N_trial\n"); print(3, 3, N_trial); }

        n4_otimes_symallcomb(cc1, cc3, N_trial, vC4, gp);

//        if (gp == 0) { printf("C4\n"); print(6, 6, vC4); }

        SIMD BL[6 * 3 * nodes];
        for (size_t n = 0; n < nodes; n++) {
            for (int j = 0; j < 3; j++) {
                BL[0 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 0];
                BL[1 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n * 3 + 1];
                BL[2 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 2] * element.dND[n * 3 + 2];
                BL[3 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 1] + element.F[j * 3 + 1] * element.dND[n * 3 + 0];
                BL[4 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n * 3 + 2] + element.F[j * 3 + 2] * element.dND[n * 3 + 1];
                BL[5 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 2] + element.F[j * 3 + 2] * element.dND[n * 3 + 0];
            }
        }

        SIMD scale = element.det * load1(element.w[gp]);
        multAtBA<6, 3 * nodes>(element.K, BL, vC4, scale);
        multAtB<3 * nodes, 6, 1>(element.nf, BL, element.vS, scale);

        SIMD S[9]; voigt6ToMatrix33(element.vS, S);
        SIMD KS[nodes * nodes];
        multABAt<nodes, 3>(KS, element.dND, S, scale);
//        printf("KS(0,0) = %+.10e\n", KS[0][0]);
        for (size_t n = 0; n < nodes; ++n) {
            for (size_t m = 0; m < nodes; ++m) {
                element.K[(n + 0 * nodes) * 3 * nodes + m + 0 * nodes] = element.K[(n + 0 * nodes) * 3 * nodes + m + 0 * nodes] + KS[n * nodes + m];
                element.K[(n + 1 * nodes) * 3 * nodes + m + 1 * nodes] = element.K[(n + 1 * nodes) * 3 * nodes + m + 1 * nodes] + KS[n * nodes + m];
                element.K[(n + 2 * nodes) * 3 * nodes + m + 2 * nodes] = element.K[(n + 2 * nodes) * 3 * nodes + m + 2 * nodes] + KS[n * nodes + m];
            }
        }

        // move to the next GP
        this->alpha += SIMD::size;
        this->invCp += 6 * SIMD::size;
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_ */
