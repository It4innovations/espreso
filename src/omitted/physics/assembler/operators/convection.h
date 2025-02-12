
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONVECTION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONVECTION_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

#include <cmath>

namespace espreso {

struct ConvectionMaterialParameters: public Operator {
    ConvectionMaterialParameters(ParametersConvection &convection, ParameterData &temp, int region, int interval, int update)
    : Operator(interval, convection.heatTransferCoeficient.gp.regions[region].isconst[interval], update),
      extTemp(convection.externalTemperature.gp.regions[region], interval),
      temp(temp, interval),
      rho(convection.rho.gp.regions[region], interval),
      heatCapacity(convection.heatCapacity.gp.regions[region], interval),
      thermalConductivity(convection.thermalConductivity.gp.regions[region], interval),
      dynamicViscosity(convection.dynamicViscosity.gp.regions[region], interval),
      dynamicViscosityTemp(convection.dynamicViscosityTemp.gp.regions[region], interval)
    {

    }

    InputParameterIterator extTemp, temp;
    OutputParameterIterator rho, heatCapacity, thermalConductivity, dynamicViscosity, dynamicViscosityTemp;
};

struct ConvectionFluidAir: public ConvectionMaterialParameters {
    ConvectionFluidAir(ParametersConvection &convection, ParameterData &temp, int region, int interval)
    : ConvectionMaterialParameters(convection, temp, region, interval,
            convection.rho.gp.regions[region].update[interval] ||
            convection.heatCapacity.gp.regions[region].update[interval] ||
            convection.thermalConductivity.gp.regions[region].update[interval] ||
            convection.dynamicViscosity.gp.regions[region].update[interval] ||
            convection.dynamicViscosityTemp.gp.regions[region].update[interval]),
      absolutePressure(convection.absolutePressure.gp.regions[region], interval)
    {

    }

    InputParameterIterator absolutePressure;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double gas_constant = 286.9;
        rho[gpindex] = absolutePressure[gpindex];
        rho[gpindex] /= (gas_constant * extTemp[gpindex]);

        if ((extTemp[gpindex] >= 200) && (extTemp[gpindex] <= 1600)) {
            heatCapacity[gpindex] = 1047.63657 - 0.372589265 * extTemp[gpindex] + 9.45304214E-4 * std::pow(extTemp[gpindex], 2.0) - 6.02409443E-7 * std::pow(extTemp[gpindex], 3.0) + 1.2858961E-10 * std::pow(extTemp[gpindex], 4.0);
        } else if (extTemp[gpindex] < 200) {
            heatCapacity[gpindex] = 1047.63657 - 0.372589265 * 200.0 + 9.45304214E-4 * std::pow(200.0, 2.0) - 6.02409443E-7 * std::pow(200.0, 3.0) + 1.2858961E-10 * std::pow(200.0, 4.0);
        } else if (extTemp[gpindex] > 1600) {
            heatCapacity[gpindex] = 1047.63657 - 0.372589265 * 1600.0 + 9.45304214E-4 * std::pow(1600.0, 2.0) - 6.02409443E-7 * std::pow(1600.0, 3.0) + 1.2858961E-10 * std::pow(1600.0, 4.0);
        }

        if ((extTemp[gpindex] >= 200) && (extTemp[gpindex] <= 1600)) {
            thermalConductivity[gpindex] = -0.00227583562 + 1.15480022E-4 * extTemp[gpindex] - 7.90252856E-8 * std::pow(extTemp[gpindex], 2.0) + 4.11702505E-11 * std::pow(extTemp[gpindex], 3.0) - 7.43864331E-15 * std::pow(extTemp[gpindex], 4.0);
        } else if (extTemp[gpindex] < 200) {
            thermalConductivity[gpindex] = -0.00227583562 + 1.15480022E-4 * 200.0 - 7.90252856E-8 * std::pow(200.0, 2.0) + 4.11702505E-11 * std::pow(200.0, 3.0) - 7.43864331E-15 * std::pow(200.0, 4.0);
        } else if (extTemp[gpindex] > 1600) {
            thermalConductivity[gpindex] = -0.00227583562 + 1.15480022E-4 * 1600.0 - 7.90252856E-8 * std::pow(1600.0, 2.0) + 4.11702505E-11 * std::pow(1600.0, 3.0) - 7.43864331E-15 * std::pow(1600.0, 4.0);
        }

        if ((extTemp[gpindex] >= 200) && (extTemp[gpindex] <= 1600)) {
            dynamicViscosity[gpindex] = -8.38278E-7 + 8.35717342E-8 * extTemp[gpindex] - 7.69429583E-11 * std::pow(extTemp[gpindex], 2.0) + 4.6437266E-14 * std::pow(extTemp[gpindex], 3.0) - 1.06585607E-17 * std::pow(extTemp[gpindex], 4.0);
        } else if (extTemp[gpindex] < 200) {
            dynamicViscosity[gpindex] = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * std::pow(200.0, 2.0) + 4.6437266E-14 * std::pow(200.0, 3.0) - 1.06585607E-17 * std::pow(200.0, 4.0);
        } else if (extTemp[gpindex] > 1600) {
            dynamicViscosity[gpindex] = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * std::pow(1600.0, 2.0) + 4.6437266E-14 * std::pow(1600.0, 3.0) - 1.06585607E-17 * std::pow(1600.0, 4.0);
        }


        if ((temp[gpindex] >= 200) && (temp[gpindex] <= 1600)) {
            dynamicViscosityTemp[gpindex] = -8.38278E-7 + 8.35717342E-8 * temp[gpindex] - 7.69429583E-11 * std::pow(temp[gpindex], 2.0) + 4.6437266E-14 * std::pow(temp[gpindex], 3.0) - 1.06585607E-17 * std::pow(temp[gpindex], 4.0);
        } else if (temp[gpindex] < 200) {
            dynamicViscosityTemp[gpindex] = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * std::pow(200.0, 2.0) + 4.6437266E-14 * std::pow(200.0, 3.0) - 1.06585607E-17 * std::pow(200.0, 4.0);
        } else if (temp[gpindex] > 1600) {
            dynamicViscosityTemp[gpindex] = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * std::pow(1600.0, 2.0) + 4.6437266E-14 * std::pow(1600.0, 3.0) - 1.06585607E-17 * std::pow(1600.0, 4.0);
        }
    }

    void operator++()
    {
        ++extTemp; ++temp;
        ++absolutePressure;
        ++rho; ++heatCapacity; ++thermalConductivity; ++dynamicViscosity; ++dynamicViscosityTemp;
    }
};

struct ConvectionOperator: public Operator {
    ConvectionOperator(ParametersConvection &convection, ParameterData &temp, int region, int interval, int update)
    : Operator(interval, convection.heatTransferCoeficient.gp.regions[region].isconst[interval], update),
      extTemp(convection.externalTemperature.gp.regions[region], interval),
      temp(temp, interval),
      rho(convection.rho.gp.regions[region], interval),
      heatCapacity(convection.heatCapacity.gp.regions[region], interval),
      thermalConductivity(convection.thermalConductivity.gp.regions[region], interval),
      dynamicViscosity(convection.dynamicViscosity.gp.regions[region], interval),
      dynamicViscosityTemp(convection.dynamicViscosityTemp.gp.regions[region], interval),
      htc(convection.heatTransferCoeficient.gp.regions[region], interval)
    {

    }

    InputParameterIterator extTemp, temp;
    InputParameterIterator rho, heatCapacity, thermalConductivity, dynamicViscosity, dynamicViscosityTemp;
    OutputParameterIterator htc;
};

struct ConvectionExternalNaturalInclinedWall: public ConvectionOperator {
    ConvectionExternalNaturalInclinedWall(ParametersConvection &convection, ParameterData &temp, int region, int interval)
    : ConvectionOperator(convection, temp, region, interval, convection.heatTransferCoeficient.gp.regions[region].update[interval]),
      wallHeight(convection.wallHeight.gp.regions[region], interval),
      tiltAngle(convection.tiltAngle.gp.regions[region], interval)
    {

    }

    InputParameterIterator wallHeight, tiltAngle;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double g = 9.81;
        double avgTemp = (extTemp[gpindex] + temp[gpindex]) / 2;
        double RaL = rho[gpindex] * rho[gpindex] * g * (1 / avgTemp) * heatCapacity[gpindex] * std::fabs(temp[gpindex] - extTemp[gpindex]) * wallHeight[gpindex] * wallHeight[gpindex] * wallHeight[gpindex] / ( thermalConductivity[gpindex] * dynamicViscosity[gpindex]);
        double angle = tiltAngle[gpindex] * M_PI / 180.0;
        if (RaL <= 1e9) {
            htc[gpindex] = (thermalConductivity[gpindex] / wallHeight[gpindex]) * (0.68 + (0.67 * std::cos(angle) * std::pow(RaL, 0.25)) / (std::pow(1 + std::pow((0.492 * thermalConductivity[gpindex]) / (dynamicViscosity[gpindex] * heatCapacity[gpindex]), 9.0 / 16.0), 4.0 / 9.0)));
        } else {
            htc[gpindex] = (thermalConductivity[gpindex] / wallHeight[gpindex]) * std::pow(0.825 + (0.387 * std::pow(RaL, 1.0 / 6.0)) / (std::pow(1 + std::pow((0.492 * thermalConductivity[gpindex]) / (dynamicViscosity[gpindex] * heatCapacity[gpindex]), 9.0 / 16.0), 8.0 / 27.0)), 2);
        }
    }

    void operator++()
    {
        ++extTemp; ++temp;
        ++rho; ++heatCapacity; ++thermalConductivity; ++dynamicViscosity;
        ++wallHeight; ++tiltAngle;
    }
};

struct ConvectionBuilder: public BoundaryOperatorBuilder {
    HeatTransferModuleOpt &kernel;
    ParametersConvection &convection;

    ConvectionBuilder(HeatTransferModuleOpt &kernel, ParametersConvection &convection)
    : BoundaryOperatorBuilder("BOUNDARY CONVECTION"), kernel(kernel), convection(convection)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        for (size_t region = 0; region < convection.configuration.regions.size(); ++region) {
            if (convection.configuration.regions[region].isset) {
                switch (convection.configuration.regions[region].settings.front()->type) {
                case ConvectionConfiguration::TYPE::USER:
                    // automatically from the expression
                    break;
                case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL:
                case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:
                case ConvectionConfiguration::TYPE::EXTERNAL_FORCED:
                case ConvectionConfiguration::TYPE::INTERNAL_FORCED:
                    switch (convection.configuration.regions[region].settings.front()->fluid) {
                    case ConvectionConfiguration::FLUID::AIR:
                        convection.absolutePressure.gp.builder->replace("TEMPERATURE", convection.externalTemperature.gp);
                        convection.rho.gp.regions[region].addInput(convection.absolutePressure.gp.regions[region]);
                        convection.heatCapacity.gp.regions[region].addInput(convection.externalTemperature.gp.regions[region]);
                        convection.thermalConductivity.gp.regions[region].addInput(convection.externalTemperature.gp.regions[region]);
                        convection.dynamicViscosity.gp.regions[region].addInput(convection.externalTemperature.gp.regions[region]);
                        convection.dynamicViscosityTemp.gp.regions[region].addInput(kernel.temp.boundary.gp.regions[region]);

                        convection.rho.gp.regions[region].resize();
                        convection.heatCapacity.gp.regions[region].resize();
                        convection.thermalConductivity.gp.regions[region].resize();
                        convection.dynamicViscosity.gp.regions[region].resize();
                        convection.dynamicViscosityTemp.gp.regions[region].resize();

                        kernel.addParameter(convection.rho.gp.regions[region]);
                        kernel.addParameter(convection.heatCapacity.gp.regions[region]);
                        kernel.addParameter(convection.thermalConductivity.gp.regions[region]);
                        kernel.addParameter(convection.dynamicViscosity.gp.regions[region]);
                        kernel.addParameter(convection.dynamicViscosityTemp.gp.regions[region]);
                        break;
                    default:
                        break;
                    }
                    break;
                }
                switch (convection.configuration.regions[region].settings.front()->type) {
                case ConvectionConfiguration::TYPE::USER: break;
                case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL:
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.rho.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.heatCapacity.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.thermalConductivity.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.dynamicViscosity.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.wallHeight.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].addInput(convection.tiltAngle.gp.regions[region]);
                    convection.heatTransferCoeficient.gp.regions[region].resize();

                    kernel.addParameter(convection.heatTransferCoeficient.gp.regions[region]);
                    break;
                case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:
                case ConvectionConfiguration::TYPE::EXTERNAL_FORCED:
                case ConvectionConfiguration::TYPE::INTERNAL_FORCED:
                    break;
                }
            }
        }
        return true;
    }

    void apply(int region, int interval)
    {
        if (convection.configuration.regions[region].isset) {
            switch (convection.configuration.regions[region].settings.front()->type) {
            case ConvectionConfiguration::TYPE::USER:
                // automatically by from the expression
                break;
            case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL:
                switch (convection.configuration.regions[region].settings.front()->fluid) {
                case ConvectionConfiguration::FLUID::AIR:
                    iterate_boundary_gps<HeatTransferModuleOpt::NGP>(ConvectionFluidAir(convection, kernel.temp.boundary.gp.regions[region], region, interval), region);
                    iterate_boundary_gps<HeatTransferModuleOpt::NGP>(ConvectionExternalNaturalInclinedWall(convection, kernel.temp.boundary.gp.regions[region], region, interval), region);
                    break;
                default:
                    break;
                }
                break;
            case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:
            case ConvectionConfiguration::TYPE::EXTERNAL_FORCED:
            case ConvectionConfiguration::TYPE::INTERNAL_FORCED:
                break;
            }
        }
    }
};


}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONVECTION_H_ */
