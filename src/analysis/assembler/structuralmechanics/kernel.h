
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_

#include "element.h"
#include "operators.h"
#include "mesh/element.h"

#include <iostream>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setElementKernel(StructuralMechanicsOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

    if constexpr(ndim == 2) {
        if (subkernels.thickness.expression) {
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.thickness.expression->evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.thickness.node[n][s] = value; }));
        }

        switch (subkernels.elasticity.coordinateSystem->type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (subkernels.elasticity.coordinateSystem->rotation.z.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (subkernels.elasticity.coordinateSystem->center.x.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->center.y.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            break;
        }

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));

        if (subkernels.acceleration.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.acceleration.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
        }
        if (subkernels.acceleration.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.acceleration.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
        }

        if (subkernels.angularVelocity.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.angularVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
        }
        if (subkernels.angularVelocity.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.angularVelocity.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
        }
    }

    if constexpr(ndim == 3) {
        switch (subkernels.elasticity.coordinateSystem->type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (subkernels.elasticity.coordinateSystem->rotation.x.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->rotation.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->rotation.y.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->rotation.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->rotation.z.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (subkernels.elasticity.coordinateSystem->center.x.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->center.y.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            if (subkernels.elasticity.coordinateSystem->center.x.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->center.y.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (subkernels.elasticity.coordinateSystem->center.z.isset) {
                subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.elasticity.coordinateSystem->center.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        }

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->young_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->young_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->young_modulus.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[2][s] = value; }));

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->poisson_ratio.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->poisson_ratio.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->poisson_ratio.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[2][s] = value; }));

        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->shear_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->shear_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.elasticity.configuration->shear_modulus.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[2][s] = value; }));

        if (subkernels.plasticity.isactive) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.plasticity.configuration->initial_yield_stress.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.initialYieldStress[s] = value; }));
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.plasticity.configuration->isotropic_hardening.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.isotropicHardening[s] = value; }));
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.plasticity.configuration->kinematic_hardening.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.kinematicHardening[s] = value; }));
            subkernels.plasticity.smallStrainTensorPlastic.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6);
            subkernels.plasticity.xi.resize(gps * (subkernels.chunks + 1) * SIMD::size * 6 + SIMD::size);
        }

        if (subkernels.acceleration.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.acceleration.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
        }
        if (subkernels.acceleration.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.acceleration.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
        }
        if (subkernels.acceleration.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.acceleration.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[2][s] = value; }));
        }
        if (subkernels.angularVelocity.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.angularVelocity.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[0][s] = value; }));
        }
        if (subkernels.angularVelocity.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.angularVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
        }
        if (subkernels.angularVelocity.expressionVector) {
            subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    subkernels.angularVelocity.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
        }
    }

    if (subkernels.material.configuration->density.isset) {
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.material.configuration->density.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.density[s] = value; }));
    }
    if (subkernels.material.configuration->heat_capacity.isset) {
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.material.configuration->heat_capacity.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatCapacity[s] = value; }));
    }

    BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
    CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
    IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
    ThicknessToNodes<nodes, ndim> thickness(subkernels.thickness);

    SIMD volume;
    basis.simd(element);
    for (size_t c = 0; c < subkernels.chunks; ++c) {
        coordinates.simd(element);
        thickness.simd(element);

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);
            volume = volume + element.det * load1(element.w[gp]);
        }
    }

    subkernels.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        subkernels.volume += volume[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runElementKernel(const step::Step &step, StructuralMechanicsOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

    BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
    CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
    ThicknessToGp<nodes, ndim> thickness(subkernels.thickness);
    TemperatureKernel<nodes> temperature(subkernels.temperature);
    TemperatureToGPsKernel<nodes> temperatureToGPs(subkernels.temperature);
    IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
    IntegrationDisplacedKernel<nodes, ndim> integrationDisplaced(subkernels.integrationDisplaced);
    DisplacementKernel<nodes, ndim> displacement(subkernels.displacement);
    SmallStrainTensorKernel<nodes, ndim> smallStrainTensor(subkernels.smallStrainTensor);
    ElasticityKernel<ndim> elasticity(subkernels.elasticity);
    ElasticityLargeDisplacementKernel<ndim> elasticityLargeDisplacement(subkernels.elasticityLargeDisplacement);
    PlasticityKernel<nodes, ndim> plasticity(subkernels.plasticity, action);
    MatrixElasticityKernel<nodes, ndim> K(subkernels.K);
    MatrixLargeDisplacementKernel<nodes, ndim> KLD(subkernels.KLD);
    MatrixMassKernel<nodes, ndim> M(subkernels.M);
    AccelerationKernel<nodes, ndim> acceleration(subkernels.acceleration);
    AngularVelocityKernel<nodes, ndim> angularVelocity(subkernels.angularVelocity);
    SigmaKernel<nodes, ndim> sigma(subkernels.sigma);
    StressKernel<nodes, gps, ndim> stress(subkernels.stress);
    MatrixFillerKernel<nodes> outK(subkernels.Kfiller);
    MatrixFillerKernel<nodes> outM(subkernels.Mfiller);
    MatrixFillerKernel<nodes> outC(subkernels.Cfiller);
    RHSFillerKernel<nodes> outReRHS(subkernels.reRHSfiller);
    RHSFillerKernel<nodes> outReNRHS(subkernels.reNRHSfiller);
    RHSFillerKernel<nodes> outImRHS(subkernels.imRHSfiller);
    RHSFillerKernel<nodes> outImNRHS(subkernels.imRHSfiller);

    struct {
        std::vector<ExternalNodeExpression<ndim, Element>*> node;
        std::vector<ExternalGPsExpression<ndim, Element>*> gp;
    } nonconst;

    for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
        ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
        if (subkernels.expressions.node[i]->evaluator->isConst()) {
            for (size_t n = 0; n < nodes; ++n) {
                exp->simd(element, n);
            }
        } else {
            nonconst.node.push_back(exp);
        }
    }

    for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
        ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
        if (subkernels.expressions.gp[i]->evaluator->isConst()) {
            for (size_t gp = 0; gp < gps; ++gp) {
                exp->simd(element, gp);
            }
        } else {
            nonconst.gp.push_back(exp);
        }
    }

    // pre-processing of possible constant parameters from ecf
    basis.simd(element);
    elasticity.simd(element);
    thickness.simd(element, 0);

    coordinatesToGPs.setActiveness(action);
    thickness.setActiveness(action);
    temperature.setActiveness(action);
    elasticity.setActiveness(action);
    elasticityLargeDisplacement.setActiveness(action, step.iteration);
    integrationDisplaced.setActiveness(action, step.iteration);
    plasticity.setActiveness(action);
    displacement.setActiveness(action);
    smallStrainTensor.setActiveness(action);
    KLD.setActiveness(action, step.iteration);
    K.setActiveness(action, !KLD.isactive);
    M.setActiveness(action);
//    C.setActiveness(action);
    acceleration.setActiveness(action);
    angularVelocity.setActiveness(action);
    sigma.setActiveness(action);
    stress.setActiveness(action);

    outK.setActiveness(action);
    outM.setActiveness(action);
    outC.setActiveness(action);
    outReRHS.setActiveness(action);
    outReNRHS.setActiveness(action, step.iteration);
    outImRHS.setActiveness(action);
    outImNRHS.setActiveness(action);

    for (size_t c = 0; c < subkernels.chunks; ++c) {
        if (sigma.isactive) {
            sigma.reset(element);
        }
        coordinates.simd(element);
        if (temperature.isactive) {
            temperature.simd(element);
        }
        if (displacement.isactive) {
            displacement.simd(element);
        }
        for (size_t i = 0; i < nonconst.node.size(); ++i) {
            for (size_t n = 0; n < nodes; ++n) {
                nonconst.node[i]->simd(element, n);
            }
        }

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);

            if (integrationDisplaced.isactive) {
                integrationDisplaced.simd(element, gp);
            }

            if (coordinatesToGPs.isactive) {
                coordinatesToGPs.simd(element, gp);
            }
            if (temperatureToGPs.isactive) {
                temperatureToGPs.simd(element, gp);
            }
            if (thickness.isactive) {
                thickness.simd(element, gp);
            }
            for (size_t i = 0; i < nonconst.gp.size(); ++i) {
                for (size_t n = 0; n < nodes; ++n) {
                    nonconst.gp[i]->simd(element, n);
                }
            }
            if (elasticity.isactive) {
                elasticity.simd(element, gp);
            }
            if (elasticityLargeDisplacement.isactive) {
                elasticityLargeDisplacement.simd(element, gp);
            }
            if (smallStrainTensor.isactive) {
                smallStrainTensor.simd(element, gp);
            }
            if (plasticity.isactive) {
                plasticity.simd(element, gp);
            }
            if (K.isactive) {
                K.simd(element, gp);
            }
            if (KLD.isactive) {
                KLD.simd(element, gp);
            }
            if (M.isactive) {
                M.simd(element, gp);
            }
//            if (C.isactive) {
//                C.simd(element, gp);
//            }

            if (acceleration.isactive) {
                acceleration.simd(element, gp);
            }
            if (angularVelocity.isactive) {
                angularVelocity.simd(element, gp);
            }
            if (sigma.isactive) {
                sigma.simd(element, gp);
            }
        }

        if (outK.isactive) {
            outK.simd(element.K);
        }
        if (outM.isactive) {
            outM.simd(element.M);
        }
        if (outC.isactive) {
            outC.simd(element.C);
        }
        if (outReRHS.isactive) {
            outReRHS.simd(element.f);
        }
        if (outReNRHS.isactive) {
            outReNRHS.simd(element.nf);
        }
        if (outImRHS.isactive) {
            outImRHS.simd(element.f);
        }
        if (outImNRHS.isactive) {
            outImNRHS.simd(element.nf);
        }
        if (stress.isactive) {
            stress.simd(element);
        }
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setBoundaryKernel(StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

    if (subkernels.normalPressure.expression) {
        subkernels.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                subkernels.normalPressure.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.normalPressure[s] = value; }));
    }

    BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
    CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
    IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
    StoreNormalKernel<nodes, ndim> storeNornal(subkernels.normal);

    storeNornal.setActiveness(action);

    SIMD surface;
    basis.simd(element);
    for (size_t c = 0; c < subkernels.chunks; ++c) {
        coordinates.simd(element);
        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);
            surface = surface + element.det * load1(element.w[gp]);
        }
        if (storeNornal.isactive) {
            storeNornal.simd(element);
        }
    }

    subkernels.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        subkernels.surface += surface[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runBoundaryKernel(const StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

    BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
    CoordinatesKernel<nodes, ndim> coordinates(subkernels.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(subkernels.coordinates);
    ThicknessFromNodes<nodes, ndim> thickness(subkernels.thickness);
    ThicknessToGp<nodes, ndim> thicknessToGPs(subkernels.thickness);
    IntegrationKernel<nodes, ndim, edim> integration(subkernels.integration);
    NormalPressureKernel<nodes, ndim> normalPressure(subkernels.normalPressure);
    RHSFillerKernel<nodes> outReRHS(subkernels.reRHSfiller);
    RHSFillerKernel<nodes> outImRHS(subkernels.imRHSfiller);

    std::vector<ExternalGPsExpression<ndim, Element>*> nonconst;
    for (size_t i = 0; i < subkernels.expressions.gp.size(); ++i) {
        ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(subkernels.expressions.gp[i]);
        if (subkernels.expressions.gp[i]->evaluator->isConst()) {
            exp->simd(element, 0);
        } else {
            nonconst.push_back(exp);
        }
    }

    basis.simd(element);
    thickness.setActiveness(action);

    outReRHS.setActiveness(action);
    outImRHS.setActiveness(action);

    for (size_t c = 0; c < subkernels.chunks; ++c) {
        coordinates.simd(element);

        if (thickness.isactive) {
            thickness.simd(element);
        }

        for (size_t gp = 0; gp < gps; ++gp) {
            if (coordinatesToGPs.isactive) {
                coordinatesToGPs.simd(element, gp);
            }

            for (size_t i = 0; i < nonconst.size(); ++i) {
                nonconst[i]->simd(element, gp);
            }

            if (thicknessToGPs.isactive) {
                thicknessToGPs.simd(element, gp);
            }

            integration.simd(element, gp);

            if (normalPressure.isactive) {
                normalPressure.simd(element, gp);
            }
        }

        if (outReRHS.isactive) {
            outReRHS.simd(element.f);
        }
        if (outImRHS.isactive) {
            outImRHS.simd(element.f);
        }
    }
}

template <size_t ndim>
void setNodeKernel(StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsDirichlet<ndim> Element;

    if constexpr(ndim == 2) {
        if (subkernels.displacement.expression) {
            if (subkernels.displacement.expression->x.isset) {
                subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        subkernels.displacement.expression->x.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[0][s] = value; }));
            }
            if (subkernels.displacement.expression->y.isset) {
                subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        subkernels.displacement.expression->y.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[1][s] = value; }));
            }
        }
        if (subkernels.harmonicForce.magnitude.expressionVector) {
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.magnitude.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[0][s] = value; }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.magnitude.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[1][s] = value; }));
        }
        if (subkernels.harmonicForce.phase.expressionVector) {
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.phase.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[0][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[0][s] = std::sin(value * M_PI / 180); }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.phase.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[1][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[1][s] = std::sin(value * M_PI / 180); }));
        }
    }

    if constexpr(ndim == 3) {
        if (subkernels.displacement.expression) {
            if (subkernels.displacement.expression->x.isset) {
                subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        subkernels.displacement.expression->x.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[0][s] = value; }));
            }
            if (subkernels.displacement.expression->y.isset) {
                subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        subkernels.displacement.expression->y.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[1][s] = value; }));
            }
            if (subkernels.displacement.expression->z.isset) {
                subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        subkernels.displacement.expression->z.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[2][s] = value; }));
            }
        }
        if (subkernels.harmonicForce.magnitude.expressionVector) {
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.magnitude.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[0][s] = value; }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.magnitude.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[1][s] = value; }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.magnitude.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[2][s] = value; }));
        }
        if (subkernels.harmonicForce.phase.expressionVector) {
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.phase.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[0][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[0][s] = std::sin(value * M_PI / 180); }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.phase.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[1][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[1][s] = std::sin(value * M_PI / 180); }));
            subkernels.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    subkernels.harmonicForce.phase.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[2][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[2][s] = std::sin(value * M_PI / 180); }));
        }
    }
}

template <size_t ndim>
void runNodeKernel(const StructuralMechanicsBoundaryOperators &subkernels, SubKernel::Action action)
{
    typedef StructuralMechanicsDirichlet<ndim> Element; Element element;

    CoordinatesKernel<1, ndim> coordinates(subkernels.coordinates);
    HarmonicForceKernel<1, ndim> harmonicForce(subkernels.harmonicForce);
    RHSFillerKernel<1> outReRHS(subkernels.reRHSfiller);
    RHSFillerKernel<1> outImRHS(subkernels.imRHSfiller);
    VectorSetterKernel<1, Element> set(subkernels.reDirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.displacement.node[d][s]; });

    harmonicForce.setActiveness(action);
    outReRHS.setActiveness(action);
    outImRHS.setActiveness(action);
    set.setActiveness(action);

    std::vector<ExternalNodeExpression<ndim, Element>*> nonconst;
    for (size_t i = 0; i < subkernels.expressions.node.size(); ++i) {
        ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(subkernels.expressions.node[i]);
        if (subkernels.expressions.node[i]->evaluator->isConst()) {
            exp->simd(element, 0);
        } else {
            nonconst.push_back(exp);
        }
    }

    for (size_t c = 0; c < subkernels.chunks; ++c) {
        coordinates.simd(element);
        for (size_t i = 0; i < nonconst.size(); ++i) {
            nonconst[i]->simd(element, 0);
        }

        if (harmonicForce.isactive) {
            harmonicForce.simd(element);
        }

        if (set.isactive) {
            set.simd(element);
        }

        if (outReRHS.isactive) {
            outReRHS.simd(element.f);
        }
        if (outImRHS.isactive) {
            outImRHS.simd(element.f);
        }
    }
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_ */
