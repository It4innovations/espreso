
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATERIAL_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATERIAL_H_

#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/material.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "config/ecf/material/elasticityproperties.h"

#include <memory>
#include <cmath>

namespace espreso {

struct MaterialStructuralMechanics: Material {
    const char* name() const { return "MaterialStructuralMechanics"; }

    MaterialStructuralMechanics()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN), nonlinear(true)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::SOLUTION;
    }

    void activate(const MaterialConfiguration *configuration, StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, bool nonlinear)
    {
        Material::activate(configuration);
        this->behaviour = behaviour;
        this->nonlinear = nonlinear;
        this->isconst = false;
        this->isactive = true;
    }

    void init(size_t chunks, size_t gps)
    {
        switch (configuration->elasticity_properties.material_model) {
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
            break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
            break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
            if (invCp.empty()) {
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
            break;
        }
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
    bool nonlinear;

    static std::vector<double> invCp, alpha;
};

template <size_t nodes, size_t gps, size_t ndim> struct MaterialStructuralMechanicsKernel;

template <size_t nodes, size_t gps> struct MaterialStructuralMechanicsKernel<nodes, gps, 2>: MaterialStructuralMechanics {
    MaterialStructuralMechanicsKernel(const MaterialStructuralMechanics &base): MaterialStructuralMechanics(base) {}

    template <typename Element>
    void init(Element &element)
    {

    }

    template <typename Element>
    void simd(Element &element)
    {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                element.F[i * 2 + j] = zeros();
                for (size_t n = 0; n < nodes; ++n) {
                    element.F[i * 2 + j] = element.F[i * 2 + j] + element.displacement[n][i] * element.dND[n * 2 + j];
                }
            }
        }
        element.F[0] = element.F[0] + load1(1.);
        element.F[3] = element.F[3] + load1(1.);

        // Voigt   0  1  2
        //        11 22 12
        set<2, 2>(element.C2, zeros());
        multAtB<2, 2, 2>(element.C2, element.F, element.F); // Right Cauchy-Green tensor

        switch (configuration->elasticity_properties.material_model) {
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
            switch (behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C1 = load1(1.);
                SIMD C2 = load1(2.);
                SIMD k = ex * (C1 - mi) / ((C1 + mi) * (C1 - C2 * mi));
                element.vC4[0] = element.vC4[4] = k;
                element.vC4[3] = element.vC4[1] = k * (mi / (C1 - mi));
                element.vC4[8] = k * ((C1 - C2 * mi) / (C2 * (C1 - mi)));
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C1 = load1(1.);
                SIMD C2 = load1(.5);
                SIMD k = ex / (C1 - mi * mi);
                element.vC4[0] = element.vC4[4] = k;
                element.vC4[3] = element.vC4[1] = k * mi;
                element.vC4[8] = k * ((C1 -  mi) * C2);
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC: {
                SIMD ex = element.ecf.youngModulus[0];
                SIMD mi = element.ecf.poissonRatio[0];
                SIMD C05 = load1(.5);
                SIMD C1 = load1(1.);
                SIMD C2 = load1(2.);
                SIMD k = ex / ((C1 + mi) * (C1 - C2 * mi));
                element.vC4[ 0] = element.vC4[5] = element.vC4[10] = k * (C1 - mi);
                element.vC4[ 4] = element.vC4[1] = element.vC4[ 8] = element.vC4[2] = element.vC4[ 9] = element.vC4[6] = k * mi;
                element.vC4[15] = k * (C1 - C2 * mi) * C05;
            } break;
            }
            switch (behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
            {
                SIMD eVec[3] = { load1(.5) * (element.C2[0] - load1(1)), load1(.5) * (element.C2[3] - load1(1)), element.C2[1] };
                set<1, 3>(element.vS, load1(0));
                multAB<3, 3, 1>(element.vS, element.vC4, eVec, load1(1.0));
            } break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
            {
                eslog::error("not implemented material behaviour\n");
            } break;
            }
            break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
        default:
            eslog::error("not implemented material model\n");
            break;
        }
    }
};

template <size_t nodes, size_t gps> struct MaterialStructuralMechanicsKernel<nodes, gps, 3>: MaterialStructuralMechanics {
    MaterialStructuralMechanicsKernel(MaterialStructuralMechanics &base)
    : MaterialStructuralMechanics(base),
      invCp(base.invCp.data()),
      alpha(base.alpha.data()),
      save(false)
    {
        switch (configuration->elasticity_properties.material_model) {
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
            break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
            break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
            invCp = utils::getAligned(SIMD::size, base.invCp);
            alpha = utils::getAligned(SIMD::size, base.alpha);
            break;
        }
    }

    SIMD logJbarDivJbar;

    double *invCp, *alpha;
    bool save;

    void setActiveness(Action action)
    {
        save = action == SubKernel::Action::SOLUTION;
        SubKernel::setActiveness(action);
    }

    template <typename Element>
    void init(Element &element)
    {
        switch (configuration->elasticity_properties.material_model) {
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
        {

        } break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
        {

        } break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
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
                Ve = Ve + load1(element.w[gp]) * determinant33(JX);

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
                SIMD detJxW = detJx * load1(element.w[gp]);
                ve = ve + detJxW;

                for (size_t n = 0; n < nodes; ++n) {
                    SIMD dNX = load1(element.dN[gp][n][0]);
                    SIMD dNY = load1(element.dN[gp][n][1]);
                    SIMD dNZ = load1(element.dN[gp][n][2]);
                    mean_dNdx[0 * nodes + n] = mean_dNdx[0 * nodes + n] + detJxW * (invJx[0] * dNX + invJx[3] * dNY + invJx[6] * dNZ);
                    mean_dNdx[1 * nodes + n] = mean_dNdx[1 * nodes + n] + detJxW * (invJx[1] * dNX + invJx[4] * dNY + invJx[7] * dNZ);
                    mean_dNdx[2 * nodes + n] = mean_dNdx[2 * nodes + n] + detJxW * (invJx[2] * dNX + invJx[5] * dNY + invJx[8] * dNZ);
                }
            }
//            for (size_t s = 0; s < SIMD::size; ++s) {
//                if (ve[s] < 0) {
//                    printf("NEGATIVE\n");
//                    std::ofstream os("negative.vtk");
//                    os << "# vtk DataFile Version 2.0\n";
//                    os << "EXAMPLE\n";
//                    os << "ASCII\n";
//                    os << "DATASET UNSTRUCTURED_GRID\n";
//                    os << "POINTS 8 float\n";
//                    for (size_t n = 0; n < nodes; ++n) {
//                        SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
//                        SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
//                        SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
//                        os << coordsX[s] << " " << coordsY[s] << " " << coordsZ[s] << "\n";
//                    }
//                    os << "CELLS 1 9\n";
//                    os << "8 0 1 2 3 4 5 6 7\n";;
//                    os << "CELL_TYPES 1\n";
//                    os << "12\n";
//                }
//            }
//            if (step.iteration == 0) {
//                for (size_t s = 0; s < SIMD::size; ++s) {
//                    std::ofstream os("alpha_" + std::to_string(chunk) + "_" + std::to_string(s) + ".txt", std::ofstream::app);
//                    os << std::showpos;
//                    os << step.loadstep << " X " << step.substep << " X " << step.iteration << "\n";
//                    for (size_t gp = 0; gp < gps; ++gp) {
//                        os << std::setw(5) << std::scientific << this->alpha[gp * SIMD::size + s] << " ";
//                    }
//                    os << "\n";
//                }
//                for (size_t s = 0; s < SIMD::size; ++s) {
//                    std::ofstream os("invCp_" + std::to_string(chunk) + "_" + std::to_string(s) + ".txt", std::ofstream::app);
//                    os << std::showpos;
//                    os << step.loadstep << " X " << step.substep << " X " << step.iteration << "\n";
//                    for (int r = 0; r < 3; ++r) {
//                        for (size_t gp = 0; gp < gps; ++gp) {
//    //                        os << std::setw(5) << std::scientific << this->alpha + gp * SIMD::size + s << "\n";
//                            SIMD invCp[9] = {
//                                    load(this->invCp + 0 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 3 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 5 * SIMD::size + gp * 6 * SIMD::size),
//                                    load(this->invCp + 3 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 1 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 4 * SIMD::size + gp * 6 * SIMD::size),
//                                    load(this->invCp + 5 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 4 * SIMD::size + gp * 6 * SIMD::size), load(this->invCp + 2 * SIMD::size + gp * 6 * SIMD::size)
//                            };
//                            os << std::setw(5) << std::scientific << invCp[r * 3 + 0][s] << " " << invCp[r * 3 + 1][s] << " " << invCp[r * 3 + 2][s] << " | ";
//                        }
//                        os << "\n";
//                    }
//                    os << "\n";
//                }
//            }
            SIMD Jbar = ve / Ve;
            divJbar = load1(1) / Jbar;
            logJbarDivJbar = log(Jbar) * divJbar;

            SIMD bulk_modulus = element.ecf.youngModulus[0] / (load1(3) - load1(6) * element.ecf.poissonRatio[0]);
            SIMD kappabar = bulk_modulus * divJbar - bulk_modulus * logJbarDivJbar;

            multABt<3 * nodes, 1, 3 * nodes>(element.K, mean_dNdx, mean_dNdx, kappabar / ve);
        } break;
        }
    }

    template <typename Element>
    void simd(Element &element)
    {
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

        // Voigt   0  1  2  3  4  5
        //        11 22 33 12 13 23
        set<3, 3>(element.C2, zeros());
        multAtB<3, 3, 3>(element.C2, element.F, element.F); // Right Cauchy-Green tensor

        switch (configuration->elasticity_properties.material_model) {
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
        {
            SIMD E = element.ecf.youngModulus[0];
            SIMD nu = element.ecf.poissonRatio[0]; // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
            SIMD lambda = E * nu / ((load1(1.) + nu) * (load1(1.) - load1(2.) * nu));
            SIMD mu = E / (load1(2.) + load1(2.) * nu);
            element.vC4[ 0] = element.vC4[ 7] = element.vC4[14] = lambda + load1(2.) * mu;
            element.vC4[ 1] = element.vC4[ 2] = element.vC4[ 8] = lambda;
            element.vC4[21] = element.vC4[28] = element.vC4[35] = mu;

            SIMD al0 = load1(.5) * (element.C2[0] + element.C2[4] + element.C2[8]) * lambda - mu - load1(.5) * load1(3) * lambda;
            element.vS[0] = al0 + mu * element.C2[0];
            element.vS[1] = al0 + mu * element.C2[4];
            element.vS[2] = al0 + mu * element.C2[8];
            element.vS[3] = mu * element.C2[1];
            element.vS[4] = mu * element.C2[5];
            element.vS[5] = mu * element.C2[2];
        }
        break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
        {
            SIMD E = element.ecf.youngModulus[0];
            SIMD nu = element.ecf.poissonRatio[0]; // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
            SIMD lambda = E * nu / ((load1(1.) + nu) * (load1(1.) - load1(2.) * nu));
            SIMD mu = E / (load1(2.) + load1(2.) * nu);

            SIMD J2 = element.C2[0] * element.C2[4] * element.C2[8] + load1(2) * element.C2[1] * element.C2[5] * element.C2[2] - element.C2[0] * element.C2[2] * element.C2[2] - element.C2[4] * element.C2[5] * element.C2[5] - element.C2[8] * element.C2[1] * element.C2[1];
            SIMD rJ2 = load1(1) / J2;
            SIMD vCinv[6] = {
                (element.C2[4] *element.C2[8] - element.C2[2] * element.C2[2] ) * rJ2,
                (element.C2[0] *element.C2[8] - element.C2[5] * element.C2[5] ) * rJ2,
                (element.C2[0] *element.C2[4] - element.C2[1] * element.C2[1] ) * rJ2,
                (element.C2[5] *element.C2[2] - element.C2[8] * element.C2[1] ) * rJ2,
                (element.C2[1] *element.C2[2] - element.C2[4] * element.C2[5] ) * rJ2,
                (element.C2[1] *element.C2[5] - element.C2[0] * element.C2[2] ) * rJ2,
            };

            SIMD C05 = load1(.5);
            SIMD logJ2 = log(J2);
            SIMD al3 = C05 * (lambda * logJ2) - mu;
            element.vS[0] = al3 * vCinv[0] + mu;
            element.vS[1] = al3 * vCinv[1] + mu;
            element.vS[2] = al3 * vCinv[2] + mu;
            element.vS[3] = al3 * vCinv[3];
            element.vS[4] = al3 * vCinv[4];
            element.vS[5] = al3 * vCinv[5];

            auto ij = [] (int i, int j) { return (i - 1) * 6 + j - 1; };
            SIMD symodot_vCinv[36];
            symodot_vCinv[ij(1,1)] =  vCinv[0] * vCinv[0];
            symodot_vCinv[ij(1,2)] =  vCinv[3] * vCinv[3];                                      symodot_vCinv[ij(2,1)] = symodot_vCinv[ij(1,2)];
            symodot_vCinv[ij(1,3)] =  vCinv[5] * vCinv[5];                                      symodot_vCinv[ij(3,1)] = symodot_vCinv[ij(1,3)];
            symodot_vCinv[ij(1,4)] =  vCinv[0] * vCinv[3];                                      symodot_vCinv[ij(4,1)] = symodot_vCinv[ij(1,4)];
            symodot_vCinv[ij(1,5)] =  vCinv[3] * vCinv[5];                                      symodot_vCinv[ij(5,1)] = symodot_vCinv[ij(1,5)];
            symodot_vCinv[ij(1,6)] =  vCinv[0] * vCinv[5];                                      symodot_vCinv[ij(6,1)] = symodot_vCinv[ij(1,6)];
            symodot_vCinv[ij(2,2)] =  vCinv[1] * vCinv[1];
            symodot_vCinv[ij(2,3)] =  vCinv[4] * vCinv[4];                                      symodot_vCinv[ij(3,2)] = symodot_vCinv[ij(2,3)];
            symodot_vCinv[ij(2,4)] =  vCinv[1] * vCinv[3];                                      symodot_vCinv[ij(4,2)] = symodot_vCinv[ij(2,4)];
            symodot_vCinv[ij(2,5)] =  vCinv[1] * vCinv[4];                                      symodot_vCinv[ij(5,2)] = symodot_vCinv[ij(2,5)];
            symodot_vCinv[ij(2,6)] =  vCinv[3] * vCinv[4];                                      symodot_vCinv[ij(6,2)] = symodot_vCinv[ij(2,6)];
            symodot_vCinv[ij(3,3)] =  vCinv[2] * vCinv[2];
            symodot_vCinv[ij(3,4)] =  vCinv[4] * vCinv[5];                                      symodot_vCinv[ij(4,3)] = symodot_vCinv[ij(3,4)];
            symodot_vCinv[ij(3,5)] =  vCinv[2] * vCinv[4];                                      symodot_vCinv[ij(5,3)] = symodot_vCinv[ij(3,5)];
            symodot_vCinv[ij(3,6)] =  vCinv[2] * vCinv[5];                                      symodot_vCinv[ij(6,3)] = symodot_vCinv[ij(3,6)];
            symodot_vCinv[ij(4,4)] = (vCinv[0] * vCinv[1] + vCinv[3] * vCinv[3]) * load1(.5);
            symodot_vCinv[ij(4,5)] = (vCinv[3] * vCinv[4] + vCinv[1] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(5,4)] = symodot_vCinv[ij(4,5)];
            symodot_vCinv[ij(4,6)] = (vCinv[0] * vCinv[4] + vCinv[3] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(6,4)] = symodot_vCinv[ij(4,6)];
            symodot_vCinv[ij(5,5)] = (vCinv[1] * vCinv[2] + vCinv[4] * vCinv[4]) * load1(.5);
            symodot_vCinv[ij(5,6)] = (vCinv[2] * vCinv[3] + vCinv[4] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(6,5)] = symodot_vCinv[ij(5,6)];
            symodot_vCinv[ij(6,6)] = (vCinv[0] * vCinv[2] + vCinv[5] * vCinv[5]) * load1(.5);

            SIMD be6 = lambda;
            SIMD be7 = load1(2) * mu - lambda * logJ2;
            for (int i = 0; i < 6; i++) {
                for (int j = i; j < 6; j++) {
                    element.vC4[i * 6 + j] = be6 * vCinv[i] * vCinv[j] + be7 * symodot_vCinv[i * 6 + j];
                }
            }
        }
        break;
        case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
        {
            SIMD sigmaY0 = element.ecf.sigma;
            SIMD Hisotropic = element.ecf.isotropicHardening;
            SIMD mu = element.ecf.youngModulus[0] / (load1(2.) + load1(2.) * element.ecf.poissonRatio[0]);
            SIMD bulk_modulus = element.ecf.youngModulus[0] / (load1(3) - load1(6) * element.ecf.poissonRatio[0]);

            SIMD detF, invF[9];
            inv33(element.F, detF, invF);

            SIMD alpha = load(this->alpha);
            SIMD invCp[9] = {
                    load(this->invCp + 0 * SIMD::size), load(this->invCp + 3 * SIMD::size), load(this->invCp + 5 * SIMD::size),
                    load(this->invCp + 3 * SIMD::size), load(this->invCp + 1 * SIMD::size), load(this->invCp + 4 * SIMD::size),
                    load(this->invCp + 5 * SIMD::size), load(this->invCp + 4 * SIMD::size), load(this->invCp + 2 * SIMD::size)
            };

            SIMD C2[9]; multAtB<3, 3, 3>(C2, element.F, element.F);
            SIMD detC2, invC2[9];
            inv33(C2, detC2, invC2);

            SIMD pbar = bulk_modulus * logJbarDivJbar;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    element.vC4[i * 6 + j] = zeros();
                    if (i < 3 && j < 3) {
                        element.vC4[i * 6 + j] = pbar;
                    }
                    if (i == j) {
                        element.vC4[i * 6 + j] = -pbar;
                    }
                }
            }

            SIMD be_trial[9]; multABAt<3, 3>(be_trial, element.F, invCp);
            SIMD eigVal[3], n_trial[9];
            eigSym33Asc(be_trial, eigVal, n_trial);
            SIMD lambda_trial[3] = { sqrt(eigVal[0]), sqrt(eigVal[1]), sqrt(eigVal[2]) };

            SIMD tau_trial_devdiag[3] = {
                    load1(2) * mu * log(lambda_trial[0]) - load1(2/3.) * mu * log(detF),
                    load1(2) * mu * log(lambda_trial[1]) - load1(2/3.) * mu * log(detF),
                    load1(2) * mu * log(lambda_trial[2]) - load1(2/3.) * mu * log(detF)
            };

            SIMD norm_tau_trial_dev = sqrt(tau_trial_devdiag[0] * tau_trial_devdiag[0] + tau_trial_devdiag[1] * tau_trial_devdiag[1] + tau_trial_devdiag[2] * tau_trial_devdiag[2]);
            SIMD f_trial = load1(std::sqrt(3/2.)) * norm_tau_trial_dev - (sigmaY0 * detF + Hisotropic * alpha);
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
                        cc3[k - 3][s] = mu[s] * (1 / detF[s]) - sigma_devdiag[i][s];
                    } else {
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

            SIMD sigma_diag[3] = { sigma_devdiag[0] + pbar, sigma_devdiag[1] + pbar, sigma_devdiag[2] + pbar };
            for (int i = 0; i < 6; ++i) {
                element.vS[i] = zeros();
            }
            for (int a = 0; a < 3; ++a) {  // element.vS is now element.vs
                element.vS[0] = element.vS[0] + sigma_diag[a] * n_trial[a * 3 + 0] * n_trial[a * 3 + 0];
                element.vS[1] = element.vS[1] + sigma_diag[a] * n_trial[a * 3 + 1] * n_trial[a * 3 + 1];
                element.vS[2] = element.vS[2] + sigma_diag[a] * n_trial[a * 3 + 2] * n_trial[a * 3 + 2];
                element.vS[3] = element.vS[3] + sigma_diag[a] * n_trial[a * 3 + 0] * n_trial[a * 3 + 1];
                element.vS[4] = element.vS[4] + sigma_diag[a] * n_trial[a * 3 + 1] * n_trial[a * 3 + 2];
                element.vS[5] = element.vS[5] + sigma_diag[a] * n_trial[a * 3 + 0] * n_trial[a * 3 + 2];
            }
            SIMD SA[3 * 3], SB[3 * 3]; // element.vs --> element.s
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    SA[i * 3 + j] = element.vS[matrix33ToVoigh6(i * 3 + j)];
                }
            }
            // pullback(invF, s, 'sharp')
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    SIMD sum;
                    for (int n = 0; n < 3; ++n) {
                        sum = sum + invF[i * 3 + n] * SA[n * 3 + j];
                    }
                    SB[i * 3 + j] = sum;
                }
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    SIMD sum;
                    for (int n = 0; n < 3; ++n) {
                        sum = sum + invF[j * 3 + n] * SB[i * 3 + n];
                    }
                    SA[i * 3 + j] = sum;
                }
            }
            for (int i = 0; i < 6; ++i) {
                element.vS[i] = detF * SA[voigt6ToMatrix33(i)];
            }

            n4_otimes_symallcomb(cc1, cc3, n_trial, element.vC4); // element.vC4 is now element.vc4
            SIMD C4A[3 * 3 * 3 * 3], C4B[3 * 3 * 3 * 3]; // element.vc4 --> element.c4

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            C4A[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l] = element.vC4[matrix3333ToVoigh36(i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l)];
                        }
                    }
                }
            }

            // pullback(invF, c4, 'sharp')
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            SIMD sum;
                            for (int n = 0; n < 3; ++n) {
                                sum = sum + invF[i * 3 + n] * C4A[n * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l];
                            }
                            C4B[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l] = sum;
                        }
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            SIMD sum;
                            for (int n = 0; n < 3; ++n) {
                                sum = sum + invF[j * 3 + n] * C4B[i * 3 * 3 * 3 + n * 3 * 3 + k * 3 + l];
                            }
                            C4A[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l] = sum;
                        }
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            SIMD sum;
                            for (int n = 0; n < 3; ++n) {
                                sum = sum + invF[k * 3 + n] * C4A[i * 3 * 3 * 3 + j * 3 * 3 + n * 3 + l];
                            }
                            C4B[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l] = sum;
                        }
                    }
                }
            }
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < 3; ++l) {
                            SIMD sum;
                            for (int n = 0; n < 3; ++n) {
                                sum = sum + invF[l * 3 + n] * C4B[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + n];
                            }
                            C4A[i * 3 * 3 * 3 + j * 3 * 3 + k * 3 + l] = sum;
                        }
                    }
                }
            }

            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    element.vC4[i * 6 + j] = detF * C4A[voigt36ToMatrix3333(i * 6 + j)];
                }
            }

            // move to the next GP
            this->alpha += SIMD::size;
            this->invCp += 6 * SIMD::size;
        }
        break;
        default:
            break;
        }

        // make the elasticity symmetric
        element.vC4[ 6] = element.vC4[ 1];
        element.vC4[12] = element.vC4[ 2]; element.vC4[13] = element.vC4[ 8];
        element.vC4[18] = element.vC4[ 3]; element.vC4[19] = element.vC4[ 9]; element.vC4[20] = element.vC4[15];
        element.vC4[24] = element.vC4[ 4]; element.vC4[25] = element.vC4[10]; element.vC4[26] = element.vC4[16]; element.vC4[27] = element.vC4[22];
        element.vC4[30] = element.vC4[ 5]; element.vC4[31] = element.vC4[11]; element.vC4[32] = element.vC4[17]; element.vC4[33] = element.vC4[23]; element.vC4[34] = element.vC4[29];
    }

    void n4_otimes_symallcomb(SIMD cc1[3], SIMD cc3[3], SIMD n[9], SIMD vc4m_hat[6*6])
    {
        for (int a = 0; a < 3; ++a) {
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[0 + 3 * a] * n[0 + 3 * a];
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc1[a] * n[0 + 3 * a] * n[0 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[a] * n[1 + 3 * a] * n[1 + 3 * a] * n[1 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[a] * n[1 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[a] * n[1 + 3 * a] * n[1 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[a] * n[1 + 3 * a] * n[1 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[a] * n[1 + 3 * a] * n[1 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[a] * n[2 + 3 * a] * n[2 + 3 * a] * n[2 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[a] * n[2 + 3 * a] * n[2 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[a] * n[2 + 3 * a] * n[2 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[a] * n[2 + 3 * a] * n[2 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[a] * n[0 + 3 * a] * n[1 + 3 * a] * n[0 + 3 * a] * n[1 + 3 * a];
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[a] * n[0 + 3 * a] * n[1 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[a] * n[0 + 3 * a] * n[1 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[a] * n[1 + 3 * a] * n[2 + 3 * a] * n[1 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[a] * n[1 + 3 * a] * n[2 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[a] * n[0 + 3 * a] * n[2 + 3 * a] * n[0 + 3 * a] * n[2 + 3 * a];
        }
        for (int k = 3, a, b; k < 6; ++k) {
            switch (k) {
            case 3: a = 0; b = 1; break;
            case 4: a = 1; b = 2; break;
            case 5: a = 0; b = 2; break;
            }
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * a] *  n[0 + 3 * b] * n[0 + 3 * b] * load1(2));
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * a] *  n[1 + 3 * b] * n[1 + 3 * b] + n[0 + 3 * b] * n[0 + 3 * b] * n[1 + 3 * a] * n[1 + 3 * a]);
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * a] *  n[2 + 3 * b] * n[2 + 3 * b] + n[0 + 3 * b] * n[0 + 3 * b] * n[2 + 3 * a] * n[2 + 3 * a]);
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[0 + 3 * a] * n[1 + 3 * b] + n[0 + 3 * b] * n[1 + 3 * a]));
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * a] *  n[1 + 3 * b] * n[2 + 3 * b] + n[0 + 3 * b] * n[0 + 3 * b] * n[1 + 3 * a] * n[2 + 3 * a]);
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]));
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * a] *  n[1 + 3 * b] * n[1 + 3 * b] * load1(2));
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * a] *  n[2 + 3 * b] * n[2 + 3 * b] + n[1 + 3 * b] * n[1 + 3 * b] * n[2 + 3 * a] * n[2 + 3 * a]);
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * b] * n[1 + 3 * a] + n[1 + 3 * b] * n[0 + 3 * a]));
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]));
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * a] *  n[0 + 3 * b] * n[2 + 3 * b] + n[1 + 3 * b] * n[1 + 3 * b] * n[0 + 3 * a] * n[2 + 3 * a]);
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * a] *  n[2 + 3 * b] * n[2 + 3 * b] * load1(2));
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * a] *  n[0 + 3 * b] * n[1 + 3 * b] + n[2 + 3 * b] * n[2 + 3 * b] * n[0 + 3 * a] * n[1 + 3 * a]);
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[1 + 3 * b] * n[2 + 3 * a] + n[2 + 3 * b] * n[1 + 3 * a]));
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[0 + 3 * b] * n[2 + 3 * a] + n[2 + 3 * b] * n[0 + 3 * a]));
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * b] *  n[1 + 3 * a] * n[1 + 3 * b] * load1(2));
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]));
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc1[k] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]));
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc1[k] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[2 + 3 * a] * n[2 + 3 * b] + n[2 + 3 * b] * n[2 + 3 * a]));
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[1 + 3 * a] * n[0 + 3 * b] + n[1 + 3 * b] * n[0 + 3 * a]));
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc1[k] * (n[2 + 3 * a] * n[2 + 3 * b] *  n[0 + 3 * a] * n[0 + 3 * b] * load1(2));
        }
        for (int k = 6, a, b; k < 9; ++k) {
            switch (k) {
            case 6: a = 0; b = 1; break;
            case 7: a = 1; b = 2; break;
            case 8: a = 0; b = 2; break;
            }
            vc4m_hat[0 * 6 + 0] = vc4m_hat[0 * 6 + 0] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] *  n[0 + 3 * a] * n[0 + 3 * b] * load1(4));
            vc4m_hat[0 * 6 + 1] = vc4m_hat[0 * 6 + 1] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] *  n[1 + 3 * a] * n[1 + 3 * b] * load1(4));
            vc4m_hat[0 * 6 + 2] = vc4m_hat[0 * 6 + 2] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] *  n[2 + 3 * a] * n[2 + 3 * b] * load1(4));
            vc4m_hat[0 * 6 + 3] = vc4m_hat[0 * 6 + 3] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[0 + 3 * a] * n[1 + 3 * b] + n[0 + 3 * b] * n[1 + 3 * a]) * load1(2));
            vc4m_hat[0 * 6 + 4] = vc4m_hat[0 * 6 + 4] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[0 * 6 + 5] = vc4m_hat[0 * 6 + 5] + cc3[k - 6] * (n[0 + 3 * a] * n[0 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[1 * 6 + 1] = vc4m_hat[1 * 6 + 1] + cc3[k - 6] * (n[1 + 3 * a] * n[1 + 3 * b] *  n[1 + 3 * a] * n[1 + 3 * b] * load1(4));
            vc4m_hat[1 * 6 + 2] = vc4m_hat[1 * 6 + 2] + cc3[k - 6] * (n[1 + 3 * a] * n[1 + 3 * b] *  n[2 + 3 * a] * n[2 + 3 * b] * load1(4));
            vc4m_hat[1 * 6 + 3] = vc4m_hat[1 * 6 + 3] + cc3[k - 6] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * a] * n[1 + 3 * b] + n[0 + 3 * b] * n[1 + 3 * a]) * load1(2));
            vc4m_hat[1 * 6 + 4] = vc4m_hat[1 * 6 + 4] + cc3[k - 6] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[1 * 6 + 5] = vc4m_hat[1 * 6 + 5] + cc3[k - 6] * (n[1 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[2 * 6 + 2] = vc4m_hat[2 * 6 + 2] + cc3[k - 6] * (n[2 + 3 * a] * n[2 + 3 * b] *  n[2 + 3 * a] * n[2 + 3 * b] * load1(4));
            vc4m_hat[2 * 6 + 3] = vc4m_hat[2 * 6 + 3] + cc3[k - 6] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[0 + 3 * a] * n[1 + 3 * b] + n[0 + 3 * b] * n[1 + 3 * a]) * load1(2));
            vc4m_hat[2 * 6 + 4] = vc4m_hat[2 * 6 + 4] + cc3[k - 6] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[2 * 6 + 5] = vc4m_hat[2 * 6 + 5] + cc3[k - 6] * (n[2 + 3 * a] * n[2 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) * load1(2));
            vc4m_hat[3 * 6 + 3] = vc4m_hat[3 * 6 + 3] + cc3[k - 6] * (n[0 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * a] * n[1 + 3 * b] + n[0 + 3 * b] * n[1 + 3 * a]) + n[0 + 3 * b] * n[1 + 3 * a] * (n[0 + 3 * b] * n[1 + 3 * a] + n[0 + 3 * a] * n[1 + 3 * b]));
            vc4m_hat[3 * 6 + 4] = vc4m_hat[3 * 6 + 4] + cc3[k - 6] * (n[0 + 3 * a] * n[1 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]) + n[0 + 3 * b] * n[1 + 3 * a] * (n[1 + 3 * b] * n[2 + 3 * a] + n[1 + 3 * a] * n[2 + 3 * b]));
            vc4m_hat[3 * 6 + 5] = vc4m_hat[3 * 6 + 5] + cc3[k - 6] * (n[0 + 3 * a] * n[1 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) + n[0 + 3 * b] * n[1 + 3 * a] * (n[0 + 3 * b] * n[2 + 3 * a] + n[0 + 3 * a] * n[2 + 3 * b]));
            vc4m_hat[4 * 6 + 4] = vc4m_hat[4 * 6 + 4] + cc3[k - 6] * (n[1 + 3 * a] * n[2 + 3 * b] * (n[1 + 3 * a] * n[2 + 3 * b] + n[1 + 3 * b] * n[2 + 3 * a]) + n[1 + 3 * b] * n[2 + 3 * a] * (n[1 + 3 * b] * n[2 + 3 * a] + n[1 + 3 * a] * n[2 + 3 * b]));
            vc4m_hat[4 * 6 + 5] = vc4m_hat[4 * 6 + 5] + cc3[k - 6] * (n[1 + 3 * a] * n[2 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) + n[1 + 3 * b] * n[2 + 3 * a] * (n[0 + 3 * b] * n[2 + 3 * a] + n[0 + 3 * a] * n[2 + 3 * b]));
            vc4m_hat[5 * 6 + 5] = vc4m_hat[5 * 6 + 5] + cc3[k - 6] * (n[0 + 3 * a] * n[2 + 3 * b] * (n[0 + 3 * a] * n[2 + 3 * b] + n[0 + 3 * b] * n[2 + 3 * a]) + n[0 + 3 * b] * n[2 + 3 * a] * (n[0 + 3 * b] * n[2 + 3 * a] + n[0 + 3 * a] * n[2 + 3 * b]));
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
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATERIAL_H_ */
