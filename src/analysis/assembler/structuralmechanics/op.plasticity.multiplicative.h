
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_

#include "basis/utilities/utils.h"
#include "analysis/assembler/general/element.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "mesh/store/nameddata.h"
#include "esinfo/meshinfo.h"

#include <memory>

namespace espreso {

struct PlasticityMultiplicative: SubKernel {
    PlasticityMultiplicative()
    : behaviour(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN), configuration(nullptr),
      isPlastized(nullptr), isPlastizedEnd(nullptr), elements(0)
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }

    void activate(size_t interval, StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour, const PlasticityPropertiesConfiguration *configuration, NamedData *isPlastized)
    {
        this->behaviour = behaviour;
        this->configuration = configuration;
        this->isconst = false;
        this->isactive = true;
        this->isPlastized = isPlastized->data.data() + info::mesh->elements->eintervals[interval].begin;
        this->isPlastizedEnd = isPlastized->data.data() + isPlastized->data.size();
    }

    void init(size_t chunks, size_t gps)
    {
        invCp.resize(chunks * gps * SIMD::size * 6 + SIMD::size);
        alpha.resize(chunks * gps * SIMD::size * 1 + SIMD::size);
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
    double *isPlastized, *isPlastizedEnd;

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

    SIMD Jbar;

    template <typename Element>
    void simd(Element &element)
    {
        SIMD Ve, ve;
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
            ve = ve + determinant33(Jx);
        }
        Jbar = ve / Ve;
    }

    void otimes(SIMD cc1[3], SIMD cc3[3], SIMD n[9], SIMD vs2[6*3], SIMD vc4m_hat[6*6])
    {
        for (int a = 0; a < 3; ++a) {
            for (int i = 0; i < 3; ++i) {
                vs2[1 * 3 + i] = vs2[1 * 3 + i] + cc1[a] * n[1 * 9 + a * 3 + i] * n[1 * 9 + a * 3 + i];
                vs2[2 * 3 + i] = vs2[2 * 3 + i] + cc1[a] * n[2 * 9 + a * 3 + i] * n[2 * 9 + a * 3 + i];
                vs2[3 * 3 + i] = vs2[3 * 3 + i] + cc1[a] * n[3 * 9 + a * 3 + i] * n[3 * 9 + a * 3 + i];
                vs2[4 * 3 + i] = vs2[4 * 3 + i] + cc1[a] * n[1 * 9 + a * 3 + i] * n[2 * 9 + a * 3 + i];
                vs2[5 * 3 + i] = vs2[5 * 3 + i] + cc1[a] * n[2 * 9 + a * 3 + i] * n[3 * 9 + a * 3 + i];
                vs2[6 * 3 + i] = vs2[6 * 3 + i] + cc1[a] * n[1 * 9 + a * 3 + i] * n[3 * 9 + a * 3 + i];
            }
        }
        for (int a = 0; a < 3; ++a) {
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a];
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[1 * 6 + 6] = vc4m_hat[1 * 6 + 6] + cc1[a] * n[1 * 3 + a] * n[1 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[2 * 6 + 6] = vc4m_hat[2 * 6 + 6] + cc1[a] * n[2 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[a] * n[3 * 3 + a] * n[3 * 3 + a] * n[3 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[a] * n[3 * 3 + a] * n[3 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[a] * n[3 * 3 + a] * n[3 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[3 * 6 + 6] = vc4m_hat[3 * 6 + 6] + cc1[a] * n[3 * 3 + a] * n[3 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[a] * n[1 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[2 * 3 + a];
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[a] * n[1 * 3 + a] * n[2 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[4 * 6 + 6] = vc4m_hat[4 * 6 + 6] + cc1[a] * n[1 * 3 + a] * n[2 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[a] * n[2 * 3 + a] * n[3 * 3 + a] * n[2 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[5 * 6 + 6] = vc4m_hat[5 * 6 + 6] + cc1[a] * n[2 * 3 + a] * n[3 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
            vc4m_hat[6 * 6 + 6] = vc4m_hat[6 * 6 + 6] + cc1[a] * n[1 * 3 + a] * n[3 * 3 + a] * n[1 * 3 + a] * n[3 * 3 + a];
        }
        for (int k = 4, a, b; k < 6; ++k) {
            switch (k) {
            case 4: a = 0; b = 1; break;
            case 5: a = 1; b = 2; break;
            case 6: a = 0; b = 2; break;
            }
            for (int i = 0; i < 3; ++i) {
                vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + a] *  n[1 * 9 + b] * n[1 * 3 + b] * load1(2);
                vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + a] *  n[2 * 9 + b] * n[2 * 3 + b] + n[1 * 3 + b] * n[1 * 3 + b] * n[2 * 3 + a] * n[2 * 3 + a];
                vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + a] *  n[3 * 9 + b] * n[3 * 3 + b] + n[1 * 3 + b] * n[1 * 3 + b] * n[3 * 3 + a] * n[3 * 3 + a];
                vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 9 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]);
                vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + a] *  n[2 * 9 + b] * n[3 * 3 + b] + n[1 * 3 + b] * n[1 * 3 + b] * n[2 * 3 + a] * n[3 * 3 + a];
                vc4m_hat[1 * 6 + 6] = vc4m_hat[1 * 6 + 6] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 9 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]);
                vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + a] *  n[2 * 9 + b] * n[2 * 3 + b] * load1(2);
                vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + a] *  n[3 * 9 + b] * n[3 * 3 + b] + n[2 * 3 + b] * n[2 * 3 + b] * n[3 * 3 + a] * n[3 * 3 + a];
                vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 9 + b] * n[2 * 3 + a] + n[2 * 3 + b] * n[1 * 3 + a]);
                vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + b] * (n[2 * 9 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]);
                vc4m_hat[2 * 6 + 6] = vc4m_hat[2 * 6 + 6] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + a] *  n[1 * 9 + b] * n[3 * 3 + b] + n[2 * 3 + b] * n[2 * 3 + b] * n[1 * 3 + a] * n[3 * 3 + a];
                vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + a] *  n[3 * 9 + b] * n[3 * 3 + b] * load1(2);
                vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + a] *  n[1 * 9 + b] * n[2 * 3 + b] + n[3 * 3 + b] * n[3 * 3 + b] * n[1 * 3 + a] * n[2 * 3 + a];
                vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + b] * (n[2 * 9 + b] * n[3 * 3 + a] + n[3 * 3 + b] * n[2 * 3 + a]);
                vc4m_hat[3 * 6 + 6] = vc4m_hat[3 * 6 + 6] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + b] * (n[1 * 9 + b] * n[3 * 3 + a] + n[3 * 3 + b] * n[1 * 3 + a]);
                vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + b] *  n[2 * 9 + a] * n[2 * 3 + b] * load1(2);
                vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 9 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]);
                vc4m_hat[4 * 6 + 6] = vc4m_hat[4 * 6 + 6] + cc1[k] * n[1 * 3 + a] * n[1 * 3 + b] * (n[2 * 9 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]);
                vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[k] * n[2 * 3 + a] * n[2 * 3 + b] * (n[3 * 9 + a] * n[3 * 3 + b] + n[3 * 3 + b] * n[3 * 3 + a]);
                vc4m_hat[5 * 6 + 6] = vc4m_hat[5 * 6 + 6] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + b] * (n[2 * 9 + a] * n[1 * 3 + b] + n[2 * 3 + b] * n[1 * 3 + a]);
                vc4m_hat[6 * 6 + 6] = vc4m_hat[6 * 6 + 6] + cc1[k] * n[3 * 3 + a] * n[3 * 3 + b] *  n[1 * 9 + a] * n[1 * 3 + b] * load1(2);
            }
        }
        for (int k = 6, i = 0, a, b; k < 9; ++k) {
            switch (k) {
            case 4: a = 0; b = 1; break;
            case 5: a = 1; b = 2; break;
            case 6: a = 0; b = 2; break;
            }
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] *  n[1 * 3 + a] * n[1 * 3 + b] * load1(4);
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] *  n[2 * 3 + a] * n[2 * 3 + b] * load1(4);
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] *  n[3 * 3 + a] * n[3 * 3 + b] * load1(4);
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2);
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] * (n[2 * 3 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[1 * 6 + 6] = vc4m_hat[1 * 6 + 6] + cc3[i] * n[1 * 3 + a] * n[1 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc3[i] * n[2 * 3 + a] * n[2 * 3 + b] *  n[2 * 3 + a] * n[2 * 3 + b] * load1(4);
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc3[i] * n[2 * 3 + a] * n[2 * 3 + b] *  n[3 * 3 + a] * n[3 * 3 + b] * load1(4);
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc3[i] * n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2);
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc3[i] * n[2 * 3 + a] * n[2 * 3 + b] * (n[2 * 3 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[2 * 6 + 6] = vc4m_hat[2 * 6 + 6] + cc3[i] * n[2 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc3[i] * n[3 * 3 + a] * n[3 * 3 + b] *  n[3 * 3 + a] * n[3 * 3 + b] * load1(4);
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc3[i] * n[3 * 3 + a] * n[3 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) * load1(2);
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc3[i] * n[3 * 3 + a] * n[3 * 3 + b] * (n[2 * 3 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[3 * 6 + 6] = vc4m_hat[3 * 6 + 6] + cc3[i] * n[3 * 3 + a] * n[3 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) * load1(2);
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc3[i] * n[1 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[2 * 3 + b] + n[1 * 3 + b] * n[2 * 3 + a]) + n[1 * 3 + b] * n[2 * 3 + a] * (n[1 * 3 + b] * n[2 * 3 + a] + n[1 * 3 + a] * n[2 * 3 + b]);
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc3[i] * n[1 * 3 + a] * n[2 * 3 + b] * (n[2 * 3 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]) + n[1 * 3 + b] * n[2 * 3 + a] * (n[2 * 3 + b] * n[3 * 3 + a] + n[2 * 3 + a] * n[3 * 3 + b]);
            vc4m_hat[4 * 6 + 6] = vc4m_hat[4 * 6 + 6] + cc3[i] * n[1 * 3 + a] * n[2 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) + n[1 * 3 + b] * n[2 * 3 + a] * (n[1 * 3 + b] * n[3 * 3 + a] + n[1 * 3 + a] * n[3 * 3 + b]);
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc3[i] * n[2 * 3 + a] * n[3 * 3 + b] * (n[2 * 3 + a] * n[3 * 3 + b] + n[2 * 3 + b] * n[3 * 3 + a]) + n[2 * 3 + b] * n[3 * 3 + a] * (n[2 * 3 + b] * n[3 * 3 + a] + n[2 * 3 + a] * n[3 * 3 + b]);
            vc4m_hat[5 * 6 + 6] = vc4m_hat[5 * 6 + 6] + cc3[i] * n[2 * 3 + a] * n[3 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) + n[2 * 3 + b] * n[3 * 3 + a] * (n[1 * 3 + b] * n[3 * 3 + a] + n[1 * 3 + a] * n[3 * 3 + b]);
            vc4m_hat[6 * 6 + 6] = vc4m_hat[6 * 6 + 6] + cc3[i] * n[1 * 3 + a] * n[3 * 3 + b] * (n[1 * 3 + a] * n[3 * 3 + b] + n[1 * 3 + b] * n[3 * 3 + a]) + n[1 * 3 + b] * n[3 * 3 + a] * (n[1 * 3 + b] * n[3 * 3 + a] + n[1 * 3 + a] * n[3 * 3 + b]);
        }
        vc4m_hat[2 * 6 + 1] = vc4m_hat[1 * 6 + 2];
        vc4m_hat[3 * 6 + 1] = vc4m_hat[1 * 6 + 3];
        vc4m_hat[3 * 6 + 2] = vc4m_hat[2 * 6 + 3];
        vc4m_hat[4 * 6 + 1] = vc4m_hat[1 * 6 + 4];
        vc4m_hat[4 * 6 + 2] = vc4m_hat[2 * 6 + 4];
        vc4m_hat[4 * 6 + 3] = vc4m_hat[3 * 6 + 4];
        vc4m_hat[5 * 6 + 1] = vc4m_hat[1 * 6 + 5];
        vc4m_hat[5 * 6 + 2] = vc4m_hat[2 * 6 + 5];
        vc4m_hat[5 * 6 + 3] = vc4m_hat[3 * 6 + 5];
        vc4m_hat[5 * 6 + 4] = vc4m_hat[4 * 6 + 5];
        vc4m_hat[6 * 6 + 1] = vc4m_hat[1 * 6 + 6];
        vc4m_hat[6 * 6 + 2] = vc4m_hat[2 * 6 + 6];
        vc4m_hat[6 * 6 + 3] = vc4m_hat[3 * 6 + 6];
        vc4m_hat[6 * 6 + 4] = vc4m_hat[4 * 6 + 6];
        vc4m_hat[6 * 6 + 5] = vc4m_hat[5 * 6 + 6];
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD mu = load1(1), sigmaY0 = load1(1), Hisotropic = load1(1);

        SIMD Jx[9];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
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

        SIMD detF, invF[9];
        inv33(element.F, detF, invF);

        SIMD alpha = load(this->alpha);
        SIMD invCp[9] = {
                load(this->invCp + 0 * SIMD::size), load(this->invCp + 3 * SIMD::size), load(this->invCp + 5 * SIMD::size),
                load(this->invCp + 3 * SIMD::size), load(this->invCp + 1 * SIMD::size), load(this->invCp + 4 * SIMD::size),
                load(this->invCp + 5 * SIMD::size), load(this->invCp + 4 * SIMD::size), load(this->invCp + 2 * SIMD::size)
        };
        // move to the next GP
        this->alpha += SIMD::size;
        this->invCp += 6 * SIMD::size;

        SIMD be_trial[9]; multABAt<3, 3>(be_trial, element.F, invCp);
        SIMD eigVal[3], eigVec[9];
        eigSym33(be_trial, eigVal, eigVec);
        SIMD lambda_trial[3] = { sqrt(eigVal[0]), sqrt(eigVal[1]), sqrt(eigVal[2]) };

        SIMD tau_trial_devdiag[3] = {
                load1(2) * mu * log(lambda_trial[0]) - load1(2/3.) * mu * log(detF),
                load1(2) * mu * log(lambda_trial[1]) - load1(2/3.) * mu * log(detF),
                load1(2) * mu * log(lambda_trial[2]) - load1(2/3.) * mu * log(detF)
        };
        SIMD norm_tau_trial_dev = sqrt(tau_trial_devdiag[0] * tau_trial_devdiag[0] + tau_trial_devdiag[1] * tau_trial_devdiag[1] + tau_trial_devdiag[2] * tau_trial_devdiag[2]);

        SIMD f_trial = load1(std::sqrt(3/2.)) * norm_tau_trial_dev - (sigmaY0 * detF + Hisotropic * alpha);
        SIMD plastized = ispositive(f_trial);
        SIMD elastic = load1(1) - plastized;

        SIMD rnorm = load1(1) / (load1(std::sqrt(3/2.)) * norm_tau_trial_dev);
        SIMD nu[3] = {
                tau_trial_devdiag[0] * rnorm,
                tau_trial_devdiag[1] * rnorm,
                tau_trial_devdiag[2] * rnorm
        };

        SIMD Dgamma = plastized * f_trial / (load1(3) * mu + Hisotropic);
        SIMD tau_devdiag[3] = { elastic * tau_trial_devdiag[0], elastic * tau_trial_devdiag[1], elastic * tau_trial_devdiag[2] };
        SIMD tau_devdiag1 = plastized * (load1(1) - (load1(2) * mu * Dgamma) * rnorm);
        tau_devdiag[0] = tau_devdiag[0] + tau_devdiag1 * tau_trial_devdiag[0];
        tau_devdiag[1] = tau_devdiag[1] + tau_devdiag1 * tau_trial_devdiag[1];
        tau_devdiag[2] = tau_devdiag[2] + tau_devdiag1 * tau_trial_devdiag[2];

        SIMD rDetF = load1(1) / detF;
        SIMD sigma_devdiag[3] = { tau_devdiag[0] * rDetF, tau_devdiag[1] * rDetF, tau_devdiag[2] * rDetF };
        SIMD sigma_diag[3] = { sigma_devdiag[0] + Jbar, sigma_devdiag[1] + Jbar, sigma_devdiag[2] + Jbar };
        SIMD sigma[9];
        multABt<3, 1, 3>(sigma, eigVec + 0, eigVec + 0, sigma_diag[0]);
        multABt<3, 1, 3>(sigma, eigVec + 3, eigVec + 3, sigma_diag[1]);
        multABt<3, 1, 3>(sigma, eigVec + 6, eigVec + 6, sigma_diag[2]);

        SIMD cc1_albe_pl[9] {
            load1( 2/3.), load1(-1/3.), load1(-1/3.),
            load1(-1/3.), load1( 2/3.), load1(-1/3.),
            load1(-1/3.), load1(-1/3.), load1( 2/3.)
        };
        SIMD cc1_albe_el[9] {
            load1( 2/3.), load1(-1/3.), load1(-1/3.),
            load1(-1/3.), load1( 2/3.), load1(-1/3.),
            load1(-1/3.), load1(-1/3.), load1( 2/3.)
        };
        SIMD nunut[9]; multABt<3, 1, 3>(nunut, nu, nu, load1(4) * mu * mu);
        for (size_t i = 0; i < 9; ++i) {
            cc1_albe_el[i] = cc1_albe_el[i] * rDetF * load1(2) * mu;
            cc1_albe_pl[i] = cc1_albe_pl[i] * rDetF * load1(2) * mu * (load1(1) - load1(2) * mu * Dgamma * rnorm);
            cc1_albe_pl[i] = cc1_albe_pl[i] - nunut[i] * (load1(1) / (load1(3) * mu + Hisotropic) - load1(sqrt(2/3.)) * Dgamma / norm_tau_trial_dev);
        }
        for (size_t i = 0; i < 3; ++i) {
            cc1_albe_el[4 * i] = cc1_albe_el[4 * i] - load1(2) * sigma_devdiag[i];
            cc1_albe_pl[4 * i] = cc1_albe_pl[4 * i] - load1(2) * sigma_devdiag[i];
        }
        for (size_t i = 0; i < 9; ++i) {
            cc1_albe_pl[i] = cc1_albe_pl[i] * plastized + cc1_albe_el[i] * elastic;
        }

        SIMD cc3_albe[3];
        for (size_t s = 0; s < SIMD::size; ++s) {
            for (int i = 0, k = 0; i < 2; ++i) {
                for (int j = i + 1; j < 2; ++j, ++k) {
                    if (lambda_trial[i][s] - lambda_trial[j][s] <= 1e-5 * lambda_trial[j][s]) {
                        cc3_albe[k][s] = mu[s] * rDetF[s] - sigma_devdiag[i][s];
                    } else {
                        cc3_albe[k][s] = (sigma_devdiag[i][s] * lambda_trial[j][s] * lambda_trial[j][s] - sigma_devdiag[j][s] * lambda_trial[i][s] * lambda_trial[i][s]) / (lambda_trial[i][s] * lambda_trial[i][s] - lambda_trial[j][s] * lambda_trial[j][s]);
                    }
                }
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_PLASTICITY_MULTIPLICATIVE_H_ */
