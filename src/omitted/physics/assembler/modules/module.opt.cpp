
#include "module.opt.hpp"

#include "basis/utilities/parser.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementsregionstore.h"
#include "physics/assembler/operators/expression.h"

#include <algorithm>

using namespace espreso;

void ModuleOpt::updateVersions()
{
    for (size_t p = 0; p < parameters.size(); ++p) {
        for (size_t ii = 0; ii < parameters[p]->version.size(); ++ii) {
            if (parameters[p]->version[ii] == -1) {
                parameters[p]->version[ii] = 0;
                parameters[p]->update[ii] = 1;
            } else {
                parameters[p]->update[ii] = 0;
                for (size_t i = 0; i < parameters[p]->inputs.size(); ++i) {
                    if (parameters[p]->version[ii] < parameters[p]->inputs[i]->version(ii)) {
                        parameters[p]->version[ii] = parameters[p]->inputs[i]->version(ii);
                        parameters[p]->update[ii] = 1;
                    }
                }
            }
        }
    }

    if (Operator::print) {
        printVersions();
    }
}

void ModuleOpt::printParamtereStats(const char* name, ParameterData &parameter)
{
    printf("parameter[update/version/isconst]:  ");
    for (size_t i = 0; i < parameter.isconst.size(); ++i) {
        if (parameter.version[i] == -1) {
            printf(" [ / / ]");
        } else {
            printf(" [%c/%d/%c]", parameter.update[i] ? 'U' : ' ', parameter.version[i], parameter.isconst[i] ? 'C' : ' ');
        }
    }
    printf(" %s\n", name);
}

void ModuleOpt::printParamtereStats(const char* name, NamedData *data)
{
    printf("nameddata[update/version/isconst]:   [-/%d/ ] %s\n", data->updated, name);
}

void ModuleOpt::setMaterials(const std::map<std::string, std::string> &settings)
{
    for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei) {
        int region = ei->region;
        if (region == -1) { // intersected regions
            for (auto rindex = ei->regions.begin(); rindex != ei->regions.end(); ++rindex) {
                if (settings.find(info::mesh->elementsRegions[*rindex]->name) != settings.end()) {
                    region = *rindex;
                }
            }
        }

        if (region == -1 || settings.find(info::mesh->elementsRegions[region]->name) == settings.end()) {
            region = 0;
        }
        auto mat = settings.find(info::mesh->elementsRegions[region]->name);
        if (mat == settings.end()) {
            eslog::error("Invalid material configuration: a region without a material settings found.\n");
        }

        for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
            if (StringCompare::caseInsensitiveEq(info::mesh->materials[i]->name, mat->second)) {
                ei->material = i;
            }
        }
    }
}

void ModuleOpt::printMaterials(const std::map<std::string, std::string> &settings)
{
    for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg) {
        auto ms = settings.find((*reg)->name);
        if (ms != settings.end()) {
            eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
        } else {
            ms = settings.find("ALL_ELEMENTS");
            if (ms != settings.end()) {
                eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
            } else {
                eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), info::mesh->materials.front()->name.c_str());
            }
        }
    }
}

void ModuleOpt::examineMaterialParameter(const std::string &material, const std::string &name, const ECFExpression &settings, ExpressionsToElements &builder, int dimension)
{
    if (settings.evaluator->variables.size()) {
        std::string params = Parser::join(", ", settings.evaluator->variables);
        eslog::info("   %18s:  %*s       FNC( %s )\n", name.c_str(), 55 - params.size(), "", params.c_str());
    } else {
        eslog::info("   %18s:  %69g \n", name.c_str(), settings.evaluator->eval(Evaluator::Params()));
    }

    for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        if (StringCompare::caseInsensitiveEq(info::mesh->materials[info::mesh->elements->eintervals[i].material]->name, material)) {
            builder.evaluators[builder.dimension * i + dimension] = settings.evaluator;
        }
    }
}

void ModuleOpt::examineElementParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExpressionsToElements &builder)
{
    examineElementParameter<ECFExpression>(name, settings, builder, 0, [] (const ECFExpression &expr) { return expr.evaluator; });
}
void ModuleOpt::examineElementParameter(const std::string &name, const std::map<std::string, ECFExpressionVector> &settings, ExpressionsToElements &builder, int dimension)
{
    examineElementParameter<ECFExpressionVector>(name, settings, builder, dimension, [&] (const ECFExpressionVector &expr) { return expr.data[dimension].evaluator; });
}

void ModuleOpt::examineBoundaryParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExpressionsToBoundary &builder)
{
    if (settings.size()) {
        eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

        int rindex = 0;
        for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
            auto ms = settings.find((*reg)->name);
            if (ms != settings.end()) {
                const Evaluator *evaluator = ms->second.evaluator;
                builder.insert(rindex, 0, evaluator);
                if (evaluator->variables.size()) {
                    std::string params = Parser::join(", ", evaluator->variables);
                    eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
                } else {
                    eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
                }
            }
        }
        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    }
}

void ModuleOpt::examineBoundaryParameter(const std::string &name, const std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection)
{
    if (settings.size()) {
        eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

        int rindex = 0;
        for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
            auto ms = settings.find((*reg)->name);
            if (ms != settings.end()) {

                auto name = [] (const std::string &name) {
                    eslog::info("   %30s:  %*s -- %s \n", "", 53 - name.size(), "", name.c_str());
                };

                auto addParam = [&] (ParametersConvection::ExternalParameter &param, const Evaluator *evaluator, const std::string &name) {
                    param.gp.builder->insert(rindex, 0, evaluator);
                    if (evaluator->variables.size()) {
                        std::string params = Parser::join(", ", evaluator->variables);
                        eslog::info("   %30s: %s %*s       FNC( %s )\n", "", name.c_str(), 43 - params.size() - name.size(), "", params.c_str());
                    } else {
                        eslog::info("   %30s: %s %*g \n", "", name.c_str(), 57 - name.size(), evaluator->eval(Evaluator::Params()));
                    }
                };

                auto setMaterial = [&] () {
                    switch (ms->second.fluid) {
                    case ConvectionConfiguration::FLUID::AIR:
                    case ConvectionConfiguration::FLUID::STEAM:
                        name("AIR");
                        addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
                        break;

                    case ConvectionConfiguration::FLUID::WATER:
                    case ConvectionConfiguration::FLUID::ENGINE_OIL:
                    case ConvectionConfiguration::FLUID::TRANSFORMER_OIL:
                        break;
                    }
                };

                convection.configuration.regions[rindex].isset = true;
                convection.configuration.regions[rindex].settings.front() = &ms->second;
                switch (ms->second.type) {

                case ConvectionConfiguration::TYPE::USER: {
                    eslog::info("   %30s:  %57s \n", ms->first.c_str(), "USER");
                    addParam(convection.heatTransferCoeficient, ms->second.heat_transfer_coefficient.evaluator, "HTC");
                    addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                } break;

                case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL:
                    eslog::info("   %30s:  %57s \n", ms->first.c_str(), "EXTERNAL NATURAL");
                    addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                    setMaterial();

                    switch (ms->second.variant) {
                    case ConvectionConfiguration::VARIANT::INCLINED_WALL:
                        name("INCLINED WALL");
                        addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
                        addParam(convection.wallHeight, ms->second.tilt_angle.evaluator, "TILT ANGLE");
                        break;
                    case ConvectionConfiguration::VARIANT::VERTICAL_WALL:
                        name("VERTICAL WALL");
                        addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
                        break;
                    case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_UP:
                        name("HORIZONTAL PLATE UP");
                        addParam(convection.length, ms->second.length.evaluator, "LENGTH");
                        break;
                    case ConvectionConfiguration::VARIANT::HORIZONTAL_PLATE_DOWN:
                        name("HORIZONTAL PLATE DOWN");
                        addParam(convection.length, ms->second.length.evaluator, "LENGTH");
                        break;
                    case ConvectionConfiguration::VARIANT::HORIZONTAL_CYLINDER:
                        name("HORIZONTAL CYLINDER");
                        addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
                        break;
                    case ConvectionConfiguration::VARIANT::SPHERE:
                        name("SPHERE");
                        addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
                        break;
                    default:
                        break;
                    }
                    break;

                case ConvectionConfiguration::TYPE::INTERNAL_NATURAL:
                    eslog::info("   %30s:  %57s \n", ms->first.c_str(), "INTERNAL NATURAL");
                    addParam(convection.length, ms->second.length.evaluator, "LENGTH");
                    setMaterial();
                    switch (ms->second.variant) {
                    case ConvectionConfiguration::VARIANT::PARALLEL_PLATES:
                        name("PARALLEL PLATES");
                        addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
                        addParam(convection.length, ms->second.length.evaluator, "LENGTH");
                        break;
                    case ConvectionConfiguration::VARIANT::CIRCULAR_TUBE:
                        name("CIRCULAR TUBE");
                        addParam(convection.wallHeight, ms->second.wall_height.evaluator, "WALL HEIGHT");
                        addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
                        break;
                    default:
                        break;
                    }
                    break;

                case ConvectionConfiguration::TYPE::EXTERNAL_FORCED:
                    setMaterial();
                    switch (ms->second.variant) {
                    case ConvectionConfiguration::VARIANT::AVERAGE_PLATE:
                        name("AVERAGE PLATE");
                        addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                        addParam(convection.length, ms->second.length.evaluator, "LENGTH");
                        break;
                    default:
                        break;
                    }
                    break;

                case ConvectionConfiguration::TYPE::INTERNAL_FORCED:
                    switch (ms->second.variant) {
                    case ConvectionConfiguration::VARIANT::TUBE:
                        setMaterial();
                        name("TUBE");
                        addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                        addParam(convection.fluidVelocity, ms->second.fluid_velocity.evaluator, "FLUID VELOCITY");
                        addParam(convection.diameter, ms->second.diameter.evaluator, "DIAMETER");
                        break;
                    case ConvectionConfiguration::VARIANT::QUENCH:
                        name("QUENCH");
                        addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
                        addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                        addParam(convection.experimentalConstant, ms->second.experimental_constant.evaluator, "EXPERIMENTAL CONSTANT");
                        break;
                    case ConvectionConfiguration::VARIANT::QUENCH_PARALLEL:
                        name("QUENCH PARALLEL");
                        addParam(convection.absolutePressure, ms->second.absolute_pressure.evaluator, "ABSOLUTE PRESSURE");
                        addParam(convection.externalTemperature, ms->second.external_temperature.evaluator, "EXTERNAL TEMPERATURE");
                        addParam(convection.experimentalConstant, ms->second.experimental_constant.evaluator, "EXPERIMENTAL CONSTANT");
                        break;
                    default:
                        break;
                    }
                    break;
                }
            }
        }
        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    }
}


