
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_KERNEL_H_

#include "element.h"
#include "operators.h"
#include "mesh/element.h"

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setElementKernel(StructuralMechanicsElementOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

    if constexpr(ndim == 2) {
        if (operators.thickness.expression) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.thickness.expression->evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.thickness.node[n][s] = value; }));
        }

        switch (operators.material.configuration->coordinate_system.type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (operators.material.configuration->coordinate_system.rotation.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (operators.material.configuration->coordinate_system.center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            break;
        }

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.young_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
//        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
//                operators.material.configuration->elasticity_properties.young_modulus.get(1, 1).evaluator,
//                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.poisson_ratio.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
//        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
//                operators.material.configuration->elasticity_properties.poisson_ratio.get(1, 1).evaluator,
//                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.shear_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
//        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
//                operators.material.configuration->elasticity_properties.shear_modulus.get(1, 1).evaluator,
//                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));

        if (operators.acceleration.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.acceleration.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
        }
        if (operators.acceleration.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.acceleration.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
        }

        if (operators.angularVelocity.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.angularVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
        }
        if (operators.angularVelocity.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.angularVelocity.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
        }

        if (operators.initVelocity.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.initVelocity.expressionVector->x.evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.velocity[n][0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.initVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.velocity[n][1][s] = value; }));
        }
    }

    if constexpr(ndim == 3) {
        switch (operators.material.configuration->coordinate_system.type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (operators.material.configuration->coordinate_system.rotation.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.rotation.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.rotation.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.rotation.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.rotation.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (operators.material.configuration->coordinate_system.center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            if (operators.material.configuration->coordinate_system.center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (operators.material.configuration->coordinate_system.center.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.material.configuration->coordinate_system.center.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        }

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.young_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[0][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.young_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[1][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.young_modulus.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.youngModulus[2][s] = value; }));

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.poisson_ratio.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[0][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.poisson_ratio.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[1][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.poisson_ratio.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.poissonRatio[2][s] = value; }));

        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.shear_modulus.get(0, 0).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[0][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.shear_modulus.get(1, 1).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[1][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.shear_modulus.get(2, 2).evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.shearModulus[2][s] = value; }));

        operators.material.init(operators.chunks, gps);

        if (operators.acceleration.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.acceleration.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[0][s] = value; }));
        }
        if (operators.acceleration.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.acceleration.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[1][s] = value; }));
        }
        if (operators.acceleration.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.acceleration.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.acceleration[2][s] = value; }));
        }
        if (operators.angularVelocity.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.angularVelocity.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[0][s] = value; }));
        }
        if (operators.angularVelocity.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.angularVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[1][s] = value; }));
        }
        if (operators.angularVelocity.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.angularVelocity.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.angularVelocity[2][s] = value; }));
        }

        if (operators.initVelocity.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.initVelocity.expressionVector->x.evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.velocity[n][0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.initVelocity.expressionVector->y.evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.velocity[n][1][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.initVelocity.expressionVector->z.evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.velocity[n][2][s] = value; }));
        }
    }

    if (operators.material.configuration->density.isset) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->density.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.density[s] = value; }));
    }
    if (operators.material.configuration->heat_capacity.isset) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->heat_capacity.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatCapacity[s] = value; }));
    }
    if (operators.material.configuration->elasticity_properties.sigma.isset) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.sigma.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.sigma[s] = value; }));
    }
    if (operators.material.configuration->elasticity_properties.Hisotropic.isset) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.material.configuration->elasticity_properties.Hisotropic.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.isotropicHardening[s] = value; }));
    }

    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    InitialVelocityKernel<nodes, ndim> velocity(operators.velocity);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    ThicknessToNodes<nodes, ndim> thickness(operators.thickness);

    struct {
        std::vector<ExternalNodeExpression<ndim, Element>*> node;
        std::vector<ExternalGPsExpression<ndim, Element>*> gp;
    } nonconst;

    for (size_t i = 0; i < operators.expressions.node.size(); ++i) {
        ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(operators.expressions.node[i]);
        if (operators.expressions.node[i]->evaluator->isConst()) {
            for (size_t n = 0; n < nodes; ++n) {
                exp->simd(element, n);
            }
        } else {
            nonconst.node.push_back(exp);
        }
    }

    SIMD volume;
    BaseFunctions<code, gps>::simd(element);
    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        thickness.simd(element);

        for (size_t i = 0; i < nonconst.node.size(); ++i) {
            for (size_t n = 0; n < nodes; ++n) {
                nonconst.node[i]->simd(element, n);
            }
        }

        if (velocity.isactive) {
            velocity.simd(element);
        }

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.coords(element, gp);
            volume = volume + element.det * load1(element.w[gp]);
        }
    }

    operators.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        operators.volume += volume[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runElementKernel(const step::Step &step, StructuralMechanicsElementOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsElement<nodes, gps, ndim, edim> Element; Element element;

    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(operators.coordinates);
    ThicknessToGp<nodes, ndim> thickness(operators.thickness);
    TemperatureKernel<nodes> temperature(operators.temperature);
    TemperatureToGPsKernel<nodes> temperatureToGPs(operators.temperature);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    DisplacementKernel<nodes, ndim> displacement(operators.displacement);
    MaterialStructuralMechanicsKernel<nodes, gps, ndim> material(operators.material);
    MatrixElasticityKernel<nodes, ndim> matrixElasticity(operators.matrixElasticity);
    MatrixCorotationKernel<code, nodes, gps, ndim> matrixCorotation(operators.matrixCorotation);
    MatrixMassKernel<nodes, ndim> M(operators.M);
    MatrixMassKernel<nodes, 1> postM(operators.postM);
    AccelerationKernel<nodes, ndim> acceleration(operators.acceleration);
    AngularVelocityKernel<nodes, ndim> angularVelocity(operators.angularVelocity);
    StressKernel<nodes, gps, ndim> stress(operators.stress);
    MatrixFillerKernel<nodes> outK(operators.Kfiller);
    MatrixFillerKernel<nodes> outM(operators.Mfiller);
    MatrixFillerKernel<nodes> outC(operators.Cfiller);
    MatrixFillerKernel<nodes> outPostM(operators.postMfiller);
    RHSFillerKernel<nodes> outReRHS(operators.reRHSfiller);
    RHSFillerKernel<nodes> outReNRHS(operators.reNRHSfiller);
    RHSFillerKernel<nodes> outImRHS(operators.imRHSfiller);
    RHSFillerKernel<nodes> outImNRHS(operators.imRHSfiller);
    RHSFillerKernel<nodes> outPostB(operators.postBfiller);
    PrintEigenValuesKernel eigvals(operators.print);

    struct {
        std::vector<ExternalNodeExpression<ndim, Element>*> node;
        std::vector<ExternalGPsExpression<ndim, Element>*> gp;
    } nonconst;

    for (size_t i = 0; i < operators.expressions.node.size(); ++i) {
        ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(operators.expressions.node[i]);
        if (operators.expressions.node[i]->evaluator->isConst()) {
            for (size_t n = 0; n < nodes; ++n) {
                exp->simd(element, n);
            }
        } else {
            nonconst.node.push_back(exp);
        }
    }

    for (size_t i = 0; i < operators.expressions.gp.size(); ++i) {
        ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(operators.expressions.gp[i]);
        if (operators.expressions.gp[i]->evaluator->isConst()) {
            for (size_t gp = 0; gp < gps; ++gp) {
                exp->simd(element, gp);
            }
        } else {
            nonconst.gp.push_back(exp);
        }
    }

    // pre-processing of possible constant parameters from ecf
    BaseFunctions<code, gps>::simd(element);
    thickness.simd(element, 0);

    coordinatesToGPs.setActiveness(action);
    thickness.setActiveness(action);
    temperature.setActiveness(action);
    material.setActiveness(action);
    displacement.setActiveness(action);
    matrixCorotation.setActiveness(action, step.loadstep || step.substep || step.iteration);
    matrixElasticity.setActiveness(action);
    M.setActiveness(action);
    postM.setActiveness(action);
//    C.setActiveness(action);
    acceleration.setActiveness(action);
    angularVelocity.setActiveness(action);
    stress.setActiveness(action);

    outK.setActiveness(action);
    outM.setActiveness(action);
    outC.setActiveness(action);
    outPostM.setActiveness(action);
    outReRHS.setActiveness(action);
    outReNRHS.setActiveness(action, step.loadstep || step.substep || step.iteration);
    outImRHS.setActiveness(action);
    outImNRHS.setActiveness(action);
    outPostB.setActiveness(action);

    for (size_t c = 0; c < operators.chunks; ++c) {
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

        if (material.isactive) {
            material.init(element);
        }

        if (stress.isactive) {
            stress.init(element);
        }

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.coords(element, gp);

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
            if (material.isactive) {
                material.simd(element);
            }
            if (stress.isactive) {
                stress.element(element);
            }

            if (matrixElasticity.isactive) {
                matrixElasticity.simd(element, gp);
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

            if (postM.isactive) {
                integration.displacement(element, gp);
                postM.simd(element, gp, element.det * load1(element.w[gp]));
            }
        }

        if (matrixCorotation.isactive) {
            matrixCorotation.simd(element);
        }

        if (outK.isactive) {
            if (eigvals.isactive) eigvals.simd("K[ELM]", ndim * nodes, element.K, 8);
            outK.simd(element.K);
        }
        if (outM.isactive) {
            if (eigvals.isactive) eigvals.simd("M[ELM]", ndim * nodes, element.M, 8);
            outM.simd(element.M);
        }
        if (outC.isactive) {
            outC.simd(element.C);
        }
        if (outPostM.isactive) {
            outPostM.simd(element.M);
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

        for (size_t n = 0; n < nodes; ++n) {
            // TODO: re-evaluate all parameters dependent on coordinates
            if (stress.isactive) {
                integration.coordsInNodes(element, n);
                material.simd(element);
                stress.node(element, n);
            }
        }

        if (outPostB.isactive) {
            outPostB.simd(element.stress, 13);
        }

        if (stress.isactive) {
            stress.move(element);
        }
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setBoundaryKernel(StructuralMechanicsFaceOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

    if (operators.normalPressure.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.normalPressure.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.normalPressure[s] = value; }));
    }

    if (operators.pressure.pressure.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.pressure.pressure.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.pressure.pressure[s] = value; }));
    }
    if (operators.pressure.direction.expressionVector) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.pressure.direction.expressionVector->x.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.pressure.direction[0][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.pressure.direction.expressionVector->y.evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.pressure.direction[1][s] = value; }));

        if constexpr(ndim == 3) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.pressure.direction.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.pressure.direction[2][s] = value; }));
        }
    }

    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    NormalKernel<nodes, ndim, edim> normal(operators.normal);

    normal.setActiveness(action);

    SIMD surface;
    BaseFunctions<code, gps>::simd(element);
    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        for (size_t gp = 0; gp < gps; ++gp) {
            integration.coords(element, gp);
            surface = surface + element.det * load1(element.w[gp]);

            if (normal.isactive) {
                normal.simd(element, gp);
            }
        }
        if (normal.isactive) {
            normal.simd(element);
        }
    }

    operators.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        operators.surface += surface[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runBoundaryKernel(const StructuralMechanicsFaceOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsBoundary<nodes, gps, ndim, edim> Element; Element element;

    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(operators.coordinates);
    DisplacementKernel<nodes, ndim> displacement(operators.displacement);
    ThicknessFromNodes<nodes, ndim> thickness(operators.thickness);
    ThicknessToGp<nodes, ndim> thicknessToGPs(operators.thickness);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    NormalKernel<nodes, ndim, edim> normal(operators.normal);
    NormalPressureKernel<nodes, ndim> normalPressure(operators.normalPressure);
    FluidPressureGatherKernel<nodes, ndim> fluidPressureGather(operators.fluidPressure);
    FluidStressGatherKernel<nodes, ndim> fluidStressGather(operators.fluidStress);
    FluidPressureKernel<nodes, ndim> fluidPressure(operators.fluidPressure);
    FluidStressKernel<nodes, ndim> fluidStress(operators.fluidStress);
    PressureKernel<nodes, ndim> pressure(operators.pressure);
    RHSFillerKernel<nodes> outReRHS(operators.reRHSfiller);
    RHSFillerKernel<nodes> outImRHS(operators.imRHSfiller);

    std::vector<ExternalGPsExpression<ndim, Element>*> nonconst;
    for (size_t i = 0; i < operators.expressions.gp.size(); ++i) {
        ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(operators.expressions.gp[i]);
        if (operators.expressions.gp[i]->evaluator->isConst()) {
            exp->simd(element, 0);
        } else {
            nonconst.push_back(exp);
        }
    }

    BaseFunctions<code, gps>::simd(element);

    thickness.setActiveness(action);
    displacement.setActiveness(action);
    fluidPressureGather.setActiveness(action);
    fluidStressGather.setActiveness(action);
    fluidPressure.setActiveness(action);
    fluidStress.setActiveness(action);
    normal.setActiveness(action);

    outReRHS.setActiveness(action);
    outImRHS.setActiveness(action);

    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);

        if (displacement.isactive) {
            displacement.simd(element);
        }

        if (fluidPressureGather.isactive) {
            fluidPressureGather.simd(element);
        }
        if (fluidStressGather.isactive) {
            fluidStressGather.simd(element);
        }

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

            integration.coords(element, gp);

            if (normal.isactive) {
                normal.simd(element, gp);
            }

            if (normalPressure.isactive) {
                normalPressure.simd(element, gp);
            }
            if (pressure.isactive) {
                pressure.simd(element, gp);
            }
            if (fluidPressure.isactive) {
                fluidPressure.simd(element, gp);
            }
            if (fluidStress.isactive) {
                fluidStress.simd(element, gp);
            }
        }

        if (normal.isactive) {
            normal.simd(element);
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
void setNodeKernel(StructuralMechanicsNodeOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsDirichlet<ndim> Element;

    if constexpr(ndim == 2) {
        if (operators.displacement.expression) {
            if (operators.displacement.expression->x.isset) {
                operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        operators.displacement.expression->x.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[0][s] = value; }));
            }
            if (operators.displacement.expression->y.isset) {
                operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        operators.displacement.expression->y.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[1][s] = value; }));
            }
        }
        if (operators.force.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.force.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.force[0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.force.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.force[1][s] = value; }));
        }
        if (operators.harmonicForce.magnitude.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.magnitude.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.magnitude.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[1][s] = value; }));
        }
        if (operators.harmonicForce.phase.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.phase.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[0][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[0][s] = std::sin(value * M_PI / 180); }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.phase.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[1][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[1][s] = std::sin(value * M_PI / 180); }));
        }
    }

    if constexpr(ndim == 3) {
        if (operators.displacement.expression) {
            if (operators.displacement.expression->x.isset) {
                operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        operators.displacement.expression->x.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[0][s] = value; }));
            }
            if (operators.displacement.expression->y.isset) {
                operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        operators.displacement.expression->y.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[1][s] = value; }));
            }
            if (operators.displacement.expression->z.isset) {
                operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                        operators.displacement.expression->z.evaluator,
                        [] (Element &element, size_t &n, size_t &s, double value) { element.displacement.node[2][s] = value; }));
            }
        }
        if (operators.force.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.force.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.force[0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.force.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.force[1][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.force.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.force[2][s] = value; }));
        }
        if (operators.harmonicForce.magnitude.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.magnitude.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[0][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.magnitude.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[1][s] = value; }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.magnitude.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceMag[2][s] = value; }));
        }
        if (operators.harmonicForce.phase.expressionVector) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.phase.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[0][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[0][s] = std::sin(value * M_PI / 180); }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.phase.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[1][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[1][s] = std::sin(value * M_PI / 180); }));
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.harmonicForce.phase.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.harmonicForceCos[2][s] = std::cos(value * M_PI / 180); element.ecf.harmonicForceSin[2][s] = std::sin(value * M_PI / 180); }));
        }
    }
}

template <size_t ndim>
void runNodeKernel(const StructuralMechanicsNodeOperators &operators, SubKernel::Action action)
{
    typedef StructuralMechanicsDirichlet<ndim> Element; Element element;

    CoordinatesKernel<1, ndim> coordinates(operators.coordinates);
    ForceKernel<ndim> force(operators.force);
    HarmonicForceKernel<ndim> harmonicForce(operators.harmonicForce);
    FluidForceKernel<ndim> fluidForce(operators.fluidForce);
    RHSFillerKernel<1> outReRHS(operators.reRHSfiller);
    RHSFillerKernel<1> outImRHS(operators.imRHSfiller);
    VectorSetterKernel<1, Element> set(operators.reDirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.displacement.node[d][s]; });

    force.setActiveness(action);
    harmonicForce.setActiveness(action);
    fluidForce.setActiveness(action);
    outReRHS.setActiveness(action);
    outImRHS.setActiveness(action);
    set.setActiveness(action);

    std::vector<ExternalNodeExpression<ndim, Element>*> nonconst;
    for (size_t i = 0; i < operators.expressions.node.size(); ++i) {
        ExternalNodeExpression<ndim, Element>* exp = dynamic_cast<ExternalNodeExpression<ndim, Element>*>(operators.expressions.node[i]);
        if (operators.expressions.node[i]->evaluator->isConst()) {
            exp->simd(element, 0);
        } else {
            nonconst.push_back(exp);
        }
    }

    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        for (size_t i = 0; i < nonconst.size(); ++i) {
            nonconst[i]->simd(element, 0);
        }

        if (force.isactive) {
            force.simd(element);
        }

        if (harmonicForce.isactive) {
            harmonicForce.simd(element);
        }

        if (fluidForce.isactive) {
            fluidForce.simd(element);
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
