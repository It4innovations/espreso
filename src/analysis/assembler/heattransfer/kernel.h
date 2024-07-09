
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_

#include "element.h"
#include "operators.h"
#include "mesh/element.h"

#include <iostream>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setElementKernel(HeatTransferElementOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferElement<nodes, gps, ndim, edim> Element; Element element;

    if constexpr(ndim == 2) {
        if (operators.thickness.expression) {
            operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                    operators.thickness.expression->evaluator,
                    [] (Element &element, size_t &n, size_t &s, double value) { element.thickness.node[n][s] = value; }));
        }

        switch (operators.conductivity.coordinateSystem->type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (operators.conductivity.coordinateSystem->rotation.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (operators.conductivity.coordinateSystem->center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            break;
        }

        Evaluator *kxx = operators.conductivity.conductivity->values.get(0, 0).evaluator;
        Evaluator *kxy = operators.conductivity.conductivity->values.get(0, 1).evaluator;
        Evaluator *kyx = operators.conductivity.conductivity->values.get(1, 0).evaluator;
        Evaluator *kyy = operators.conductivity.conductivity->values.get(1, 1).evaluator;
        if (operators.conductivity.conductivity->model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = element.ecf.conductivity[3][s] = value; }));
        } else {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = value; }));
        }
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kxy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[1][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kyx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[2][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kyy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[3][s] = value; }));

        if (operators.advection.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.advection.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[0][s] = value; }));
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.advection.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[1][s] = value; }));
        }
    }

    if constexpr(ndim == 3) {
        switch (operators.conductivity.coordinateSystem->type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:
            if (operators.conductivity.coordinateSystem->rotation.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->rotation.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->rotation.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->rotation.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->rotation.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->rotation.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
            if (operators.conductivity.coordinateSystem->center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:
            if (operators.conductivity.coordinateSystem->center.x.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[0][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->center.y.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[1][s] = value; }));
            }
            if (operators.conductivity.coordinateSystem->center.z.isset) {
                operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.conductivity.coordinateSystem->center.z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.rotation.center[2][s] = value; }));
            }
            break;
        }

        Evaluator *kxx = operators.conductivity.conductivity->values.get(0, 0).evaluator;
        Evaluator *kxy = operators.conductivity.conductivity->values.get(0, 1).evaluator;
        Evaluator *kxz = operators.conductivity.conductivity->values.get(0, 2).evaluator;
        Evaluator *kyx = operators.conductivity.conductivity->values.get(1, 0).evaluator;
        Evaluator *kyy = operators.conductivity.conductivity->values.get(1, 1).evaluator;
        Evaluator *kyz = operators.conductivity.conductivity->values.get(1, 2).evaluator;
        Evaluator *kzx = operators.conductivity.conductivity->values.get(2, 0).evaluator;
        Evaluator *kzy = operators.conductivity.conductivity->values.get(2, 1).evaluator;
        Evaluator *kzz = operators.conductivity.conductivity->values.get(2, 2).evaluator;

        if (operators.conductivity.conductivity->model == ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = element.ecf.conductivity[4][s] = element.ecf.conductivity[9][s] = value; }));
        } else {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                kxx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[0][s] = value; }));
        }
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kxy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[1][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kxz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[2][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kyx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[3][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kyy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[4][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kyz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[5][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kzx, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[6][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kzy, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[7][s] = value; }));
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
            kzz, [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[8][s] = value; }));

        if (operators.advection.expressionVector) {
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.advection.expressionVector->x.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[0][s] = value; }));
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.advection.expressionVector->y.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[1][s] = value; }));
            operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                    operators.advection.expressionVector->z.evaluator,
                    [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[2][s] = value; }));
        }
    }

    if (operators.initTemperature.expression) {
        operators.expressions.node.push_back(new ExternalNodeExpression<ndim, Element>(
                operators.initTemperature.expression->evaluator,
                [] (Element &element, size_t &n, size_t &s, double value) { element.temperature.initial[n][s] = value; }));
    }

    if (operators.heatSource.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.heatSource.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatSource[s] = value; }));
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

    BasisKernel<code, nodes, gps, edim> basis(operators.basis);
    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    InitialTemperatureKernel<nodes> temperature(operators.temperature);
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
    basis.simd(element);
    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        thickness.simd(element);

        for (size_t i = 0; i < nonconst.node.size(); ++i) {
            for (size_t n = 0; n < nodes; ++n) {
                nonconst.node[i]->simd(element, n);
            }
        }

        temperature.simd(element);

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);
            volume = volume + element.det * load1(element.w[gp]);
        }
    }

    operators.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        operators.volume += volume[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runElementKernel(const step::Step &step, const step::Time &time, const HeatTransferElementOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferElement<nodes, gps, ndim, edim> Element; Element element;

    BasisKernel<code, nodes, gps, edim> basis(operators.basis);
    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(operators.coordinates);
    ThicknessToGp<nodes, ndim> thickness(operators.thickness);
    TemperatureKernel<nodes> temperature(operators.temperature);
    TemperatureToGPsKernel<nodes> temperatureToGPs(operators.temperature);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    ConductivityKernel<ndim> conductivity(operators.conductivity);
    HeatSourceKernel<nodes> heatSource(operators.heatSource);
    AdvectionKernel<nodes, ndim> advection(operators.advection);
    MatrixConductivityKernel<nodes, ndim> K(operators.K);
    MatrixMassKernel<nodes, 1> M(operators.M);
    TemperatureGradientKernel<nodes, gps, ndim> gradient(operators.gradient);
    TemperatureFluxKernel<nodes, gps, ndim> flux(operators.flux);
    MatrixApplyKernel<nodes> tempResidual(operators.temperatureResidual);
    MatrixFillerKernel<nodes> outK(operators.Kfiller);
    MatrixFillerKernel<nodes> outM(operators.Mfiller);
    RHSFillerKernel<nodes> outRHS(operators.RHSfiller);
    RHSFillerKernel<nodes> outnRHS(operators.nRHSfiller);

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
    basis.simd(element);
    conductivity.simd(element);
    thickness.simd(element, 0);

    coordinatesToGPs.setActiveness(action);
    thickness.setActiveness(action);
    temperature.setActiveness(action);
    conductivity.setActiveness(action);
    heatSource.setActiveness(action);
    advection.setActiveness(action);
    K.setActiveness(action);
    M.setActiveness(action);
    gradient.setActiveness(action);
    flux.setActiveness(action);
    tempResidual.setActiveness(action);

    outK.setActiveness(action);
    outRHS.setActiveness(action);
    outnRHS.setActiveness(action);

    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        if (temperature.isactive) {
            temperature.simd(element);
        }
        for (size_t i = 0; i < nonconst.node.size(); ++i) {
            for (size_t n = 0; n < nodes; ++n) {
                nonconst.node[i]->simd(element, n);
            }
        }

        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);

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

            if (conductivity.isactive) {
                conductivity.simd(element, gp);
            }
            if (heatSource.isactive) {
                heatSource.simd(element, gp);
            }
            if (advection.isactive) {
                advection.simd(element, gp);
            }
            if (K.isactive) {
                K.simd(element, gp);
            }
            if (M.isactive) {
                M.simd(element, gp);
            }
            if (gradient.isactive) {
                gradient.simd(element, gp);
            }
            if (flux.isactive) {
                flux.simd(element, gp);
            }
        }

        if (gradient.isactive) {
            gradient.store(element);
        }

        if (flux.isactive) {
            flux.store(element);
        }

        if (outK.isactive) {
            tempResidual.simd(element.nf, element.K, element.temperature.node, time.timeIntegrationConstantK);
            outK.simd(element.K);
        }
        if (outM.isactive) {
            tempResidual.simd(element.nf, element.M, element.temperature.node, time.timeIntegrationConstantM);
            outM.simd(element.M);
        }
        if (outRHS.isactive) {
            outRHS.simd(element.f);
        }
        if (outnRHS.isactive) {
            outnRHS.simd(element.nf);
        }
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void setBoundaryKernel(HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferBoundary<nodes, gps, ndim, edim> Element; Element element;

    if (operators.heatFlow.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.heatFlow.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlow[s] = value; }));
    }
    if (operators.heatFlux.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.heatFlux.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatFlux[s] = value; }));
    }

    if (operators.htc.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.htc.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.htc[s] = value; }));
    }
    if (operators.externalTemperature.expression) {
        operators.expressions.gp.push_back(new ExternalGPsExpression<ndim, Element>(
                operators.externalTemperature.expression->evaluator,
                [] (Element &element, size_t &gp, size_t &s, double value) { element.ecf.extTemp[s] = value; }));
    }

    BasisKernel<code, nodes, gps, edim> basis(operators.basis);
    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);

    SIMD surface;
    basis.simd(element);
    for (size_t c = 0; c < operators.chunks; ++c) {
        coordinates.simd(element);
        for (size_t gp = 0; gp < gps; ++gp) {
            integration.simd(element, gp);
            surface = surface + element.det * load1(element.w[gp]);
        }
    }

    operators.esize = sizeof(Element);
    for (size_t s = 0; s < SIMD::size; ++s) {
        operators.surface += surface[s];
    }
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runBoundaryKernel(const HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferBoundary<nodes, gps, ndim, edim> Element; Element element;

    BasisKernel<code, nodes, gps, edim> basis(operators.basis);
    CoordinatesKernel<nodes, ndim> coordinates(operators.coordinates);
    CoordinatesToGPsKernel<nodes, ndim> coordinatesToGPs(operators.coordinates);
    ThicknessFromNodes<nodes, ndim> thickness(operators.thickness);
    ThicknessToGp<nodes, ndim> thicknessToGPs(operators.thickness);
    IntegrationKernel<nodes, ndim, edim> integration(operators.integration);
    ExternalHeatKernel<nodes> externalHeat(operators.externalHeat);
    RHSFillerKernel<nodes> outRHS(operators.RHSfiller);

    std::vector<ExternalGPsExpression<ndim, Element>*> nonconst;
    for (size_t i = 0; i < operators.expressions.gp.size(); ++i) {
        ExternalGPsExpression<ndim, Element>* exp = dynamic_cast<ExternalGPsExpression<ndim, Element>*>(operators.expressions.gp[i]);
        if (operators.expressions.gp[i]->evaluator->isConst()) {
            exp->simd(element, 0);
        } else {
            nonconst.push_back(exp);
        }
    }

    basis.simd(element);
    thickness.setActiveness(action);

    outRHS.setActiveness(action);

    for (size_t c = 0; c < operators.chunks; ++c) {
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

            if (externalHeat.isactive) {
                externalHeat.simd(element, gp);
            }
        }

        if (outRHS.isactive) {
            outRHS.simd(element.f);
        }
    }
}

template <size_t ndim>
void setNodeKernel(HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferNode<ndim> Element; Element element;
    if (operators.temperature.expression) {
        auto setter = [] (Element &element, size_t &n, size_t &s, double value) { element.temperature.node[0][s] = element.temperature.initial[0][s] = value; };
        switch (info::mesh->dimension) {
        case 2: operators.expressions.node.push_back(new ExternalNodeExpression<2, Element>(operators.temperature.expression->evaluator, setter)); break;
        case 3: operators.expressions.node.push_back(new ExternalNodeExpression<3, Element>(operators.temperature.expression->evaluator, setter)); break;
        }
    }

    CoordinatesKernel<1, ndim> coordinates(operators.coordinates);
    InitialTemperatureKernel<1> initTemperature(operators.initialTemperature);

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
        if (initTemperature.isactive) {
            initTemperature.simd(element);
        }
    }
}

template <size_t ndim>
void runNodeKernel(const HeatTransferBoundaryOperators &operators, SubKernel::Action action)
{
    typedef HeatTransferNode<ndim> Element; Element element;

    CoordinatesKernel<1, ndim> coordinates(operators.coordinates);
    VectorSetterKernel<1, Element> set(operators.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.temperature.node[0][s]; });

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
        set.simd(element);
    }
}

}



#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_KERNEL_H_ */
