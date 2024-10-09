
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixElasticity: SubKernel {
    const char* name() const { return "StructuralMechanicsMatrixElasticity"; }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
    LinearElasticPropertiesConfiguration::MODEL model;
    bool rotated;

    MatrixElasticity()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN),
      model(LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC),
      rotated(false)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, LinearElasticPropertiesConfiguration::MODEL model, bool rotated)
    {
        this->behaviour = behaviour;
        this->model = model;
        this->rotated = rotated;
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct MatrixElasticityKernel;

template <size_t nodes>
struct MatrixElasticityKernel<nodes, 2>: MatrixElasticity {
    MatrixElasticityKernel(const MatrixElasticity &base): MatrixElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (behaviour) {
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
        {
            // C
            // 0 1 _
            //   0 _
            //     2

            // B * C * Bt
            //
            // C = 3x3
            // B = dX  0 dY
            //      0 dY dX

            SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
            SIMD c0 = element.elasticity[0];
            SIMD c1 = element.elasticity[1];
            SIMD c2 = element.elasticity[8];
            for (size_t n = 0, _n, _m; n < nodes; ++n) {
                _n = n; _m = n + 1;
                SIMD a = element.dND[n * 2 + 0] * c0;
                SIMD b;
                SIMD c = element.dND[n * 2 + 1] * c2;
                element.K[_n * 2 * nodes + _n] = element.K[_n * 2 * nodes + _n] + scale * (a * element.dND[n * 2 + 0] + c * element.dND[n * 2 + 1]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD xx = scale * (a * element.dND[m * 2 + 0] + c * element.dND[m * 2 + 1]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + xx;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + xx;
                }

                _n = n; _m = nodes;
                b = element.dND[n * 2 + 0] * c1;
                c = element.dND[n * 2 + 1] * c2;
                for (size_t m = 0; m < nodes; ++m, ++_m) {
                    SIMD xy = scale * (b * element.dND[m * 2 + 1] + c * element.dND[m * 2 + 0]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + xy;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + xy;
                }

                _n = n + nodes; _m = n + nodes + 1;
                b = element.dND[n * 2 + 1] * c0;
                c = element.dND[n * 2 + 0] * c2;
                element.K[_n * 2 * nodes + _n] = element.K[_n * 2 * nodes + _n] + scale * (b * element.dND[n * 2 + 1] + c * element.dND[n * 2 + 0]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD yy = scale * (b * element.dND[m * 2 + 1] + c * element.dND[m * 2 + 0]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + yy;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + yy;
                }
            }
        } break;
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
        {
            // C
            // 0 1 1 _
            //   0 1 _
            //     0 _
            //       2

            // B * C * Bt
            //
            // C = 4x4
            // B = dX  0  C dY
            //      0 dY  0 dX
            SIMD scale = element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0];
            SIMD c0 = element.elasticity[0];
            SIMD c1 = element.elasticity[1];
            SIMD c2 = element.elasticity[15];
            for (size_t n = 0, _n, _m; n < nodes; ++n) {
                SIMD coo = load1(element.N[gp][n]) / element.coords.gp[0];

                _n = n; _m = n + 1;
                SIMD a = element.dND[n * 2 + 0] * c0 + coo * c1;
                SIMD b;
                SIMD c = element.dND[n * 2 + 0] * c1 + coo * c0;
                SIMD d = element.dND[n * 2 + 1] * c2;
                element.K[_n * 2 * nodes + _n] = element.K[_n * 2 * nodes + _n] + scale * (a * element.dND[n * 2 + 0] + c * coo + d * element.dND[n * 2 + 1]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD xx = scale * (a * element.dND[m * 2 + 0] + c * load1(element.N[gp][m]) / element.coords.gp[0] + d * element.dND[m * 2 + 1]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + xx;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + xx;
                }

                _n = n; _m = nodes;
                b = element.dND[n * 2 + 0] * c1 + coo * c1;
                d = element.dND[n * 2 + 1] * c2;
                for (size_t m = 0; m < nodes; ++m, ++_m) {
                    SIMD xy = scale * (b * element.dND[m * 2 + 1] + d * element.dND[m * 2 + 0]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + xy;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + xy;
                }

                _n = n + nodes; _m = n + nodes + 1;
                b = element.dND[n * 2 + 1] * c0;
                d = element.dND[n * 2 + 0] * c2;
                element.K[_n * 2 * nodes + _n] = element.K[_n * 2 * nodes + _n] + scale * (b * element.dND[n * 2 + 1] + d * element.dND[n * 2 + 0]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD yy = scale * (b * element.dND[m * 2 + 1] + d * element.dND[m * 2 + 0]);
                    element.K[_n * 2 * nodes + _m] = element.K[_n * 2 * nodes + _m] + yy;
                    element.K[_m * 2 * nodes + _n] = element.K[_m * 2 * nodes + _n] + yy;
                }
            }
        } break;
        }
    }
};

template <size_t nodes>
struct MatrixElasticityKernel<nodes, 3>: MatrixElasticity {
    MatrixElasticityKernel(const MatrixElasticity &base): MatrixElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        if (rotated && model != LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC) {
            SIMD scale = element.det * load1(element.w[gp]);
            SIMD c00 = element.elasticity[0];
            SIMD c01 = element.elasticity[1], c11 = element.elasticity[ 7];
            SIMD c02 = element.elasticity[2], c12 = element.elasticity[ 8], c22 = element.elasticity[14];
            SIMD c03 = element.elasticity[3], c13 = element.elasticity[ 9], c23 = element.elasticity[15], c33 = element.elasticity[21];
            SIMD c04 = element.elasticity[4], c14 = element.elasticity[10], c24 = element.elasticity[16], c34 = element.elasticity[22], c44 = element.elasticity[28];
            SIMD c05 = element.elasticity[5], c15 = element.elasticity[11], c25 = element.elasticity[17], c35 = element.elasticity[23], c45 = element.elasticity[29], c55 = element.elasticity[35];
            for (size_t n = 0, _n, _m, _p; n < nodes; ++n) {
                SIMD a = element.dND[n * 3 + 0] * c00 + element.dND[n * 3 + 1] * c03 + element.dND[n * 3 + 2] * c05;
                SIMD b;
                SIMD c;
                SIMD d = element.dND[n * 3 + 0] * c03 + element.dND[n * 3 + 1] * c33 + element.dND[n * 3 + 2] * c35;
                SIMD e;
                SIMD f = element.dND[n * 3 + 0] * c05 + element.dND[n * 3 + 1] * c35 + element.dND[n * 3 + 2] * c55;
                _n = n; _m = n + 1;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (a * element.dND[n * 3 + 0] + d * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD xx = scale * (a * element.dND[m * 3 + 0] + d * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xx;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xx;
                }

                _n = n; _m = nodes, _p = 2 * nodes;
                b = element.dND[n * 3 + 0] * c01 + element.dND[n * 3 + 1] * c13 + element.dND[n * 3 + 2] * c15;
                c = element.dND[n * 3 + 0] * c02 + element.dND[n * 3 + 1] * c23 + element.dND[n * 3 + 2] * c25;
                d = element.dND[n * 3 + 0] * c03 + element.dND[n * 3 + 1] * c33 + element.dND[n * 3 + 2] * c35;
                e = element.dND[n * 3 + 0] * c04 + element.dND[n * 3 + 1] * c34 + element.dND[n * 3 + 2] * c45;
                f = element.dND[n * 3 + 0] * c05 + element.dND[n * 3 + 1] * c35 + element.dND[n * 3 + 2] * c55;
                for (size_t m = 0; m < nodes; ++m, ++_m, ++_p) {
                    SIMD xy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0] + e * element.dND[m * 3 + 2]);
                    SIMD xz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xy;
                    element.K[_n * 3 * nodes + _p] = element.K[_n * 3 * nodes + _p] + xz;
                    element.K[_p * 3 * nodes + _n] = element.K[_p * 3 * nodes + _n] + xz;
                }

                _n = n + nodes; _m = n + nodes + 1;
                b = element.dND[n * 3 + 1] * c11 + element.dND[n * 3 + 0] * c13 + element.dND[n * 3 + 2] * c14;
                d = element.dND[n * 3 + 1] * c13 + element.dND[n * 3 + 0] * c33 + element.dND[n * 3 + 2] * c34;
                e = element.dND[n * 3 + 1] * c14 + element.dND[n * 3 + 0] * c34 + element.dND[n * 3 + 2] * c44;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (b * element.dND[n * 3 + 1] + d * element.dND[n * 3 + 0] + e * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD yy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0] + e * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yy;
                }

                _n = n + nodes; _m = nodes * 2;
                c = element.dND[n * 3 + 1] * c12 + element.dND[n * 3 + 0] * c23 + element.dND[n * 3 + 2] * c24;
                e = element.dND[n * 3 + 1] * c14 + element.dND[n * 3 + 0] * c34 + element.dND[n * 3 + 2] * c44;
                f = element.dND[n * 3 + 1] * c15 + element.dND[n * 3 + 0] * c35 + element.dND[n * 3 + 2] * c45;
                for (size_t m = 0; m < nodes; ++m, ++_m) {
                    SIMD yz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yz;
                }

                _n = n + nodes * 2; _m = n + nodes * 2 + 1;
                c = element.dND[n * 3 + 2] * c22 + element.dND[n * 3 + 1] * c24 + element.dND[n * 3 + 0] * c25;
                e = element.dND[n * 3 + 2] * c24 + element.dND[n * 3 + 1] * c44 + element.dND[n * 3 + 0] * c45;
                f = element.dND[n * 3 + 2] * c25 + element.dND[n * 3 + 1] * c45 + element.dND[n * 3 + 0] * c55;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (c * element.dND[n * 3 + 2] + e * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 0]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD zz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + zz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + zz;
                }
            }
            return;
        }

        switch (model) {
        case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
        {
            // C
            // 0 1 1 _ _ _
            //   0 1 _ _ _
            //     0 _ _ _
            //       2 _ _
            //         2 _
            //           2

            // B * C * Bt
            //
            // C = 6x6
            // B = dX  0  0 dY  0 dZ
            //      0 dY  0 dX dZ  0
            //      0  0 dZ  0 dY dX
            //     - - - - - - - - -
            //      0  6 12 18 24 30
            //      a  b  c  d  e  f
            SIMD scale = element.det * load1(element.w[gp]);
            SIMD c0 = element.elasticity[0];
            SIMD c1 = element.elasticity[1];
            SIMD c2 = element.elasticity[21];
            for (size_t n = 0, _n, _m, _p; n < nodes; ++n) {
                SIMD a = element.dND[n * 3 + 0] * c0;
                SIMD b;
                SIMD c;
                SIMD d = element.dND[n * 3 + 1] * c2;
                SIMD e;
                SIMD f = element.dND[n * 3 + 2] * c2;
                _n = n; _m = n + 1;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (a * element.dND[n * 3 + 0] + d * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD xx = scale * (a * element.dND[m * 3 + 0] + d * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xx;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xx;
                }

                _n = n; _m = nodes, _p = 2 * nodes;
                b = element.dND[n * 3 + 0] * c1;
                c = element.dND[n * 3 + 0] * c1;
                d = element.dND[n * 3 + 1] * c2;
                f = element.dND[n * 3 + 2] * c2;
                for (size_t m = 0; m < nodes; ++m, ++_m, ++_p) {
                    SIMD xy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0]);
                    SIMD xz = scale * (c * element.dND[m * 3 + 2] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xy;
                    element.K[_n * 3 * nodes + _p] = element.K[_n * 3 * nodes + _p] + xz;
                    element.K[_p * 3 * nodes + _n] = element.K[_p * 3 * nodes + _n] + xz;
                }

                _n = n + nodes; _m = n + nodes + 1;
                b = element.dND[n * 3 + 1] * c0;
                d = element.dND[n * 3 + 0] * c2;
                e = element.dND[n * 3 + 2] * c2;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (b * element.dND[n * 3 + 1] + d * element.dND[n * 3 + 0] + e * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD yy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0] + e * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yy;
                }

                _n = n + nodes; _m = nodes * 2;
                c = element.dND[n * 3 + 1] * c1;
                e = element.dND[n * 3 + 2] * c2;
                for (size_t m = 0; m < nodes; ++m, ++_m) {
                    SIMD yz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yz;
                }

                _n = n + nodes * 2; _m = n + nodes * 2 + 1;
                c = element.dND[n * 3 + 2] * c0;
                e = element.dND[n * 3 + 1] * c2;
                f = element.dND[n * 3 + 0] * c2;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (c * element.dND[n * 3 + 2] + e * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 0]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD zz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + zz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + zz;
                }
            }

        } break;
        case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
        {
            // C
            // 0 1 2 _ _ _
            //   3 4 _ _ _
            //     5 _ _ _
            //       6 _ _
            //         7 _
            //           8

            // B * C * Bt
            //
            // C = 6x6
            // B = dX  0  0 dY  0 dZ
            //      0 dY  0 dX dZ  0
            //      0  0 dZ  0 dY dX
            //     - - - - - - - - -
            //      0  6 12 18 24 30
            //      a  b  c  d  e  f
            SIMD scale = element.det * load1(element.w[gp]);
            SIMD c00 = element.elasticity[0];
            SIMD c01 = element.elasticity[1], c11 = element.elasticity[7];
            SIMD c02 = element.elasticity[2], c12 = element.elasticity[8], c22 = element.elasticity[14];
            SIMD c33 = element.elasticity[21];
            SIMD c44 = element.elasticity[28];
            SIMD c55 = element.elasticity[35];
            for (size_t n = 0, _n, _m, _p; n < nodes; ++n) {
                SIMD a = element.dND[n * 3 + 0] * c00;
                SIMD b;
                SIMD c;
                SIMD d = element.dND[n * 3 + 1] * c33;
                SIMD e;
                SIMD f = element.dND[n * 3 + 2] * c55;
                _n = n; _m = n + 1;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (a * element.dND[n * 3 + 0] + d * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD xx = scale * (a * element.dND[m * 3 + 0] + d * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xx;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xx;
                }

                b = element.dND[n * 3 + 0] * c01;
                c = element.dND[n * 3 + 0] * c02;
                d = element.dND[n * 3 + 1] * c33;
                f = element.dND[n * 3 + 2] * c55;
                _n = n; _m = nodes; _p = nodes * 2;
                for (size_t m = 0; m < nodes; ++m, ++_m, ++_p) {
                    SIMD xy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0]);
                    SIMD xz = scale * (c * element.dND[m * 3 + 2] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + xy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + xy;
                    element.K[_n * 3 * nodes + _p] = element.K[_n * 3 * nodes + _p] + xz;
                    element.K[_p * 3 * nodes + _n] = element.K[_p * 3 * nodes + _n] + xz;
                }

                b = element.dND[n * 3 + 1] * c11;
                d = element.dND[n * 3 + 0] * c33;
                e = element.dND[n * 3 + 2] * c44;
                _n = n + nodes; _m = n + nodes + 1;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (b * element.dND[n * 3 + 1] + d * element.dND[n * 3 + 0] + e * element.dND[n * 3 + 2]);
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD yy = scale * (b * element.dND[m * 3 + 1] + d * element.dND[m * 3 + 0] + e * element.dND[m * 3 + 2]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yy;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yy;
                }

                c = element.dND[n * 3 + 1] * c12;
                e = element.dND[n * 3 + 2] * c44;
                _n = n + nodes; _m = nodes * 2;
                for (size_t m = 0; m < nodes; ++m, ++_m) {
                    SIMD yz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + yz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + yz;
                }

                c = element.dND[n * 3 + 2] * c22;
                e = element.dND[n * 3 + 1] * c44;
                f = element.dND[n * 3 + 0] * c55;
                _n = n + nodes * 2; _m = n + nodes * 2 + 1;
                element.K[_n * 3 * nodes + _n] = element.K[_n * 3 * nodes + _n] + scale * (c * element.dND[n * 3 + 2] + e * element.dND[n * 3 + 1] + f * element.dND[n * 3 + 0]);;
                for (size_t m = n + 1; m < nodes; ++m, ++_m) {
                    SIMD zz = scale * (c * element.dND[m * 3 + 2] + e * element.dND[m * 3 + 1] + f * element.dND[m * 3 + 0]);
                    element.K[_n * 3 * nodes + _m] = element.K[_n * 3 * nodes + _m] + zz;
                    element.K[_m * 3 * nodes + _n] = element.K[_m * 3 * nodes + _n] + zz;
                }
            }
        } break;
        case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
        {

        } break;
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_ */
