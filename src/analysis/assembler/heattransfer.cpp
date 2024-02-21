
#include "heattransfer.h"
#include "assembler.hpp"

#include "analysis/assembler/heattransfer/element.h"
#include "analysis/math/matrix_feti.h"
#include "analysis/math/vector_feti.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/domainsurfacestore.h"
#include "wrappers/bem/w.bem.h"

#include "math/math.h"

#include <numeric>
#include <algorithm>

#include <iostream>
namespace espreso {

NodeData* HeatTransfer::Results::temperature = nullptr;
NodeData* HeatTransfer::Results::initialTemperature = nullptr;
NodeData* HeatTransfer::Results::thickness = nullptr;
ElementData* HeatTransfer::Results::translationMotion = nullptr;
ElementData* HeatTransfer::Results::gradient = nullptr;
ElementData* HeatTransfer::Results::flux = nullptr;

HeatTransfer::HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
    threaded = configuration.solver == HeatTransferLoadStepConfiguration::SOLVER::FETI;
    subkernels.resize(info::mesh->elements->eintervals.size());
    boundary.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
        } else {
            boundary[r].resize(info::env::threads);
        }
    }

    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                subkernels[i].code = info::mesh->elements->eintervals[i].code;
                subkernels[i].elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
                subkernels[i].chunks = subkernels[i].elements / SIMD::size + (subkernels[i].elements % SIMD::size ? 1 : 0);
            }

            for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
                if (info::mesh->boundaryRegions[r]->dimension) {
                    for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
                        boundary[r][i].code = info::mesh->boundaryRegions[r]->eintervals[i].code;
                        boundary[r][i].elements = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;
                        boundary[r][i].chunks = boundary[r][i].elements / SIMD::size + (boundary[r][i].elements % SIMD::size ? 1 : 0);
                    }
                } else {
                    boundary[r][t].code = static_cast<int>(Element::CODE::POINT1);
                    boundary[r][t].elements = info::mesh->boundaryRegions[r]->nodes->datatarray().size(t);
                    boundary[r][t].chunks = boundary[r][t].elements / SIMD::size + (boundary[r][t].elements % SIMD::size ? 1 : 0);
                }
            }
        }
    }


    GaussPoints<Element::CODE::LINE2    ,  2, HeatTransferGPC::LINE2    , 1>::set();
    GaussPoints<Element::CODE::TRIANGLE3,  3, HeatTransferGPC::TRIANGLE3, 2>::set();
    GaussPoints<Element::CODE::SQUARE4  ,  4, HeatTransferGPC::SQUARE4  , 2>::set();
    GaussPoints<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4   , 3>::set();
    GaussPoints<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5 , 3>::set();
    GaussPoints<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6  , 3>::set();
    GaussPoints<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8    , 3>::set();
    GaussPoints<Element::CODE::LINE3    ,  3, HeatTransferGPC::LINE3    , 1>::set();
    GaussPoints<Element::CODE::TRIANGLE6,  6, HeatTransferGPC::TRIANGLE6, 2>::set();
    GaussPoints<Element::CODE::SQUARE8  ,  8, HeatTransferGPC::SQUARE8  , 2>::set();
    GaussPoints<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10  , 3>::set();
    GaussPoints<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13, 3>::set();
    GaussPoints<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15 , 3>::set();
    GaussPoints<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20   , 3>::set();
}

void HeatTransfer::analyze()
{
    double start = eslog::time();
    eslog::info("\n ============================================================================================= \n");

    validateRegionSettings("MATERIAL", settings.material_set);
    validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
    validateRegionSettings("THICKNESS", settings.thickness);

    if (Results::temperature == nullptr) {
        Results::temperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
    }
    if (Results::initialTemperature == nullptr) {
        Results::initialTemperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
    }
    if (Results::thickness == nullptr && info::mesh->dimension == 2) {
        Results::thickness = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "THICKNESS");
    }
    if (info::ecf->output.results_selection.translation_motions && Results::translationMotion == nullptr) {
        Results::translationMotion = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
    }
    if (info::ecf->output.results_selection.gradient && Results::gradient == nullptr) {
        Results::gradient = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
    }
    if (info::ecf->output.results_selection.flux && Results::flux == nullptr) {
        Results::flux = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
    }

    eslog::info(" ============================================================================================= \n");
    bool correct = true;

    if (configuration.temperature.size()) {
        correct &= checkBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature);
    }

    if (settings.initial_temperature.size()) {
        correct &= checkElementParameter("INITIAL TEMPERATURE", settings.initial_temperature);
    }

    if (true) { // add check
        if (info::mesh->dimension == 2) {
            correct &= checkElementParameter("THICKNESS", settings.thickness);
        }
        ///////////////////////////////////// Set materials and check if there is not any incorrect region intersection
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        eslog::info("\n  MATERIALS                                                                                    \n");
        eslog::info(" --------------------------------------------------------------------------------------------- \n");
        for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
            eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
            MaterialConfiguration *mat = info::mesh->materials[i];

            switch (mat->coordinate_system.type) {
            case CoordinateSystemConfiguration::TYPE::CARTESIAN:
                eslog::info("    COORDINATE SYSTEM:                                                              CARTESIAN \n");
                if (info::mesh->dimension == 2) {
                    correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
                }
                if (info::mesh->dimension == 3) {
                    correct &= checkExpression("ROTATION.X", mat->coordinate_system.rotation.x);
                    correct &= checkExpression("ROTATION.Y", mat->coordinate_system.rotation.y);
                    correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
                }
                break;
            case CoordinateSystemConfiguration::TYPE::SPHERICAL:
                if (info::mesh->dimension == 2) {
                    eslog::error("SPHERICAL coordinate system is not supported in 2D.\n");
                }
                if (info::mesh->dimension == 3) {
                    eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
                    correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
                    correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
                    correct &= checkExpression("CENTER.Z", mat->coordinate_system.center.z);
                }
                break;
            case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
                eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
                if (info::mesh->dimension == 2) {
                    correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
                    correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
                }
                if (info::mesh->dimension == 3) {
                    correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
                    correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
                }
                break;
            }
            eslog::info("                                                                                               \n");

            correct &= checkExpression("DENSITY", mat->density);
            correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
            eslog::info("                                                                                               \n");

        switch (mat->thermal_conductivity.model) {
            case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
                eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
                correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                break;
            case ThermalConductivityConfiguration::MODEL::DIAGONAL:
                eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
                if (info::mesh->dimension == 2) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                }
                if (info::mesh->dimension == 3) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                    correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
                }
                break;
            case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
                eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
                if (info::mesh->dimension == 2) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                    correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
                }
                if (info::mesh->dimension == 3) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                    correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
                    correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
                    correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
                    correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
                }
                break;
            case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
                eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
                if (info::mesh->dimension == 2) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                    correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
                    correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(1, 0));
                }
                if (info::mesh->dimension == 3) {
                    correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
                    correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
                    correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
                    correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
                    correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
                    correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
                    correct &= checkExpression("KYX", mat->thermal_conductivity.values.get(1, 0));
                    correct &= checkExpression("KZY", mat->thermal_conductivity.values.get(2, 1));
                    correct &= checkExpression("KZX", mat->thermal_conductivity.values.get(2, 0));
                }
                break;
            }
            eslog::info("                                                                                               \n");
        }

        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
        printMaterials(settings.material_set);
        eslog::info(" ============================================================================================= \n");
    }

    for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        bem[info::mesh->elements->eintervals[i].domain - info::mesh->domains->offset] = isBEM(i);

        const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
        bool cartesian = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
        bool gptemp = mat->thermal_conductivity.needTemperature();
        esint ebegin = info::mesh->elements->eintervals[i].begin, eend = info::mesh->elements->eintervals[i].end;

        if (info::mesh->dimension == 2) {
            subkernels[i].thickness.activate(getExpression(i, settings.thickness), info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::thickness->data.data());
        }
        subkernels[i].material.activate(mat);
        subkernels[i].conductivity.activate(&mat->thermal_conductivity, &mat->coordinate_system);
        subkernels[i].heatSource.activate(getExpression(i, configuration.heat_source));
        subkernels[i].advection.activate(getExpression(i, configuration.translation_motions), settings.sigma);

        bool gpcoo =
                subkernels[i].conductivity.needCoordinates ||
                subkernels[i].material.needCoordinates ||
                subkernels[i].thickness.needCoordinates ||
                subkernels[i].heatSource.needCoordinates ||
                subkernels[i].advection.needCoordinates;

        subkernels[i].coordinates.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, !cartesian || gpcoo);

        auto model = mat->thermal_conductivity.model;
        if (subkernels[i].conductivity.rotated) {
            switch (mat->thermal_conductivity.model) {
            case ThermalConductivityConfiguration::MODEL::ISOTROPIC:   model = ThermalConductivityConfiguration::MODEL::ISOTROPIC  ; break;
            case ThermalConductivityConfiguration::MODEL::DIAGONAL:    model = ThermalConductivityConfiguration::MODEL::SYMMETRIC  ; break;
            case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   model = ThermalConductivityConfiguration::MODEL::SYMMETRIC  ; break;
            case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: model = ThermalConductivityConfiguration::MODEL::ANISOTROPIC; break;
            }
        }
        subkernels[i].K.activate(model);
        if (configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
            subkernels[i].M.activate();
        }

        if (Results::gradient != nullptr) {
            subkernels[i].gradient.activate(i, Results::gradient);
        }
        if (Results::flux != nullptr) {
            subkernels[i].flux.activate(i, Results::flux, model);
            subkernels[i].conductivity.action |= SubKernel::SOLUTION;
//            subkernels[i].material.action |= SubKernel::SOLUTION;
//            subkernels[i].advection.action |= SubKernel::SOLUTION;
        }

        subkernels[i].initTemperature.activate(getExpression(i, settings.initial_temperature));
        subkernels[i].temperature.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::temperature->data.data(), gptemp);
    }

    for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (info::mesh->boundaryRegions[r]->dimension) {
            for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
                if (info::mesh->dimension == 2) {
                    boundary[r][i].thickness.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, Results::thickness->data.data());
                }
                boundary[r][i].heatFlow.activate(getExpression(info::mesh->boundaryRegions[r]->name, configuration.heat_flow));
                boundary[r][i].heatFlux.activate(getExpression(info::mesh->boundaryRegions[r]->name, configuration.heat_flux));

                auto convection = configuration.convection.find(info::mesh->boundaryRegions[r]->name);
                if (convection != configuration.convection.end()) {
                    boundary[r][i].htc.activate(&convection->second.heat_transfer_coefficient);
                    boundary[r][i].externalTemperature.activate(&convection->second.external_temperature);
                }
                if (boundary[r][i].heatFlow.isactive | boundary[r][i].heatFlux.isactive | boundary[r][i].htc.isactive) {
                    boundary[r][i].externalHeat.activate();
                }

                bool gpcoords =
                        boundary[r][i].thickness.needCoordinates ||
                        boundary[r][i].heatFlow.needCoordinates ||
                        boundary[r][i].heatFlux.needCoordinates ||
                        boundary[r][i].htc.needCoordinates ||
                        boundary[r][i].externalTemperature.needCoordinates;

                boundary[r][i].coordinates.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, gpcoords);
            }
        } else {
            for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
                boundary[r][t].coordinates.activate(region->nodes->cbegin(t), region->nodes->cend(), false);
            }
        }
    }

    for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        for (size_t t = 0; t < region->nodes->threads(); ++t) {
            boundary[r][t].temperature.activate(it->second);
            if (settings.init_temp_respect_bc) {
                boundary[r][t].initialTemperature.activate(region->nodes->cbegin(t), region->nodes->cend(), Results::temperature->data.data(), false);
            }
        }
    }

    assemble(SubKernel::Action::PREPROCESS);
    size_t esize = 0;
    std::vector<double> volume(subkernels.size()), surface(boundary.size());
    for (size_t i = 0; i < subkernels.size(); ++i) {
        esize = std::max(subkernels[i].esize, esize);
        volume[i] = subkernels[i].volume;
    }
    for (size_t r = 1; r < boundary.size(); ++r) {
        for (size_t i = 0; i < boundary[r].size(); ++i) {
            surface[r] += boundary[r][i].surface;
        }
    }
    printElementVolume(volume);
    printBoundarySurface(surface);
    for (size_t r = 1; r < boundary.size(); ++r) {
        for (size_t i = 0; i < boundary[r].size(); ++i) {
            boundary[r][i].surface = boundary[r][i].externalHeat.area = surface[r];
            boundary[r][i].surface = surface[r];
        }
    }
    Results::initialTemperature->data = Results::temperature->data;

    constant.K = true; // currently everything or nothing
    for (size_t i = 0; i < subkernels.size(); ++i) {
        for (size_t e = 0; e < subkernels[i].expressions.node.size(); ++e) {
            constant.K &= subkernels[i].expressions.node[e]->evaluator->isConst();
        }
        for (size_t e = 0; e < subkernels[i].expressions.gp.size(); ++e) {
            constant.K &= subkernels[i].expressions.gp[e]->evaluator->isConst();
        }
    }
    constant.M = constant.f = constant.K;

    constant.f = true;
    for (size_t r = 1; r < boundary.size(); ++r) {
        for (size_t i = 0; i < boundary[r].size(); ++i) {
            for (size_t e = 0; e < boundary[r][i].expressions.node.size(); ++e) {
                constant.f = boundary[r][i].expressions.node[e]->evaluator->isConst();
            }
            for (size_t e = 0; e < boundary[r][i].expressions.gp.size(); ++e) {
                constant.f = boundary[r][i].expressions.gp[e]->evaluator->isConst();
            }
        }
    }
    constant.dirichlet = constant.f;

    std::string constants;
    if (constant.K)         { constants += "K"; }
    if (constant.M)         { if (constants.size()) constants += ", "; constants += "M"; }
    if (constant.f)         { if (constants.size()) constants += ", "; constants += "f"; }
    if (constant.dirichlet) { if (constants.size()) constants += ", "; constants += "dirichlet"; }

    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
    eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize);
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    if (correct) {
        eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
        eslog::info("  CONSTANT MATRICES, VECTORS         %*s  [ %s ] \n", 50 - constants.size(), " ", constants.c_str());
    } else {
        eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
    }
    eslog::info(" ============================================================================================= \n");
}

void HeatTransfer::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
    Matrix_FETI<double> *KBEM = dynamic_cast<Matrix_FETI<double>*>(K);
    for (size_t i = 0; i < bem.size(); ++i) { // when BEM, K is FETI matrix
        if (bem[i]) {
            BETI[i] = KBEM->domains[i].vals;
        }
    }

    for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        subkernels[i].Kfiller.activate(i, 1, subkernels[i].elements, K);
        subkernels[i].Mfiller.activate(i, 1, subkernels[i].elements, M);
        subkernels[i].RHSfiller.activate(i, 1, subkernels[i].elements, f);
        subkernels[i].nRHSfiller.activate(i, 1, subkernels[i].elements, nf);
    }

    for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
                boundary[r][i].RHSfiller.activate(r, i, 1, boundary[r][i].elements, f);
            }
        }
    }
    for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            boundary[r][t].dirichlet.activate(r, t, 1, boundary[r][t].elements, dirichlet);
        }
    }
}

void HeatTransfer::run(SubKernel::Action action, size_t interval)
{
    switch (subkernels[interval].code) {
    case static_cast<size_t>(Element::CODE::TRIANGLE3): runElement<Element::CODE::TRIANGLE3>(action, interval); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE6): runElement<Element::CODE::TRIANGLE6>(action, interval); break;
    case static_cast<size_t>(Element::CODE::SQUARE4  ): runElement<Element::CODE::SQUARE4  >(action, interval); break;
    case static_cast<size_t>(Element::CODE::SQUARE8  ): runElement<Element::CODE::SQUARE8  >(action, interval); break;
    case static_cast<size_t>(Element::CODE::TETRA4   ): runElement<Element::CODE::TETRA4   >(action, interval); break;
    case static_cast<size_t>(Element::CODE::TETRA10  ): runElement<Element::CODE::TETRA10  >(action, interval); break;
    case static_cast<size_t>(Element::CODE::PYRAMID5 ): runElement<Element::CODE::PYRAMID5 >(action, interval); break;
    case static_cast<size_t>(Element::CODE::PYRAMID13): runElement<Element::CODE::PYRAMID13>(action, interval); break;
    case static_cast<size_t>(Element::CODE::PRISMA6  ): runElement<Element::CODE::PRISMA6  >(action, interval); break;
    case static_cast<size_t>(Element::CODE::PRISMA15 ): runElement<Element::CODE::PRISMA15 >(action, interval); break;
    case static_cast<size_t>(Element::CODE::HEXA8    ): runElement<Element::CODE::HEXA8    >(action, interval); break;
    case static_cast<size_t>(Element::CODE::HEXA20   ): runElement<Element::CODE::HEXA20   >(action, interval); break;
    }
}

void HeatTransfer::run(SubKernel::Action action, size_t region, size_t interval)
{
    switch (boundary[region][interval].code) {
    case static_cast<size_t>(Element::CODE::POINT1   ): runBoundary<Element::CODE::POINT1   >(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::LINE2    ): runBoundary<Element::CODE::LINE2    >(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::LINE3    ): runBoundary<Element::CODE::LINE3    >(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE3): runBoundary<Element::CODE::TRIANGLE3>(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE6): runBoundary<Element::CODE::TRIANGLE6>(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::SQUARE4  ): runBoundary<Element::CODE::SQUARE4  >(action, region, interval); break;
    case static_cast<size_t>(Element::CODE::SQUARE8  ): runBoundary<Element::CODE::SQUARE8  >(action, region, interval); break;
    }
}

void HeatTransfer::runBEM(SubKernel::Action action, size_t domain, double *BETI)
{
    if (action == SubKernel::Action::ASSEMBLE) {
        int np = info::mesh->domainsSurface->dnodes[domain].size();
        double *points = &(info::mesh->domainsSurface->coordinates[domain][0].x);
        int ne = info::mesh->domainsSurface->edistribution[domain + 1] - info::mesh->domainsSurface->edistribution[domain];
        int *elemNodes = info::mesh->domainsSurface->denodes[domain].data();

        Matrix_Dense<double> K; K.resize(np, np);
        BEM3DLaplace(K.vals, np, points, ne, elemNodes, 1);
//        math::store(K, "K");

        for (int r = 0, cc = 0; r < K.nrows; ++r) {
            for (int c = r; c < K.ncols; ++c) {
                BETI[cc++] = K.vals[r * K.ncols + c];
            }
        }

        // V: ne * ne;
        // K: ne * np;
        // D: np * np;
        // M: ne * np;
//        Matrix_Dense<double> V; V.resize(ne, ne);
//        Matrix_Dense<double> K; K.resize(ne, np);
//        Matrix_Dense<double> D; D.resize(np, np);
//        Matrix_Dense<double> M; M.resize(ne, np);
//        int order = 7;

//        BEM3DLaplace(np, points, ne, elemNodes, order, V.vals, K.vals, D.vals, M.vals);
//        // make matrix V symmetric
//        for (int r = 0; r < V.nrows; ++r) {
//            for (int c = r + 1; c < V.ncols; ++c) {
//                double avg = (V.vals[r * V.ncols + c] + V.vals[c * V.ncols + r]) / 2;
//                V.vals[r * V.ncols + c] = V.vals[c * V.ncols + r] = avg;
//            }
//        }
//
//        math::store(K, "TB_K");
//        math::store(M, "TB_M");
//        math::store(V, "TB_V");
//        math::store(D, "TB_D");

        // S = D + (.5 * M + K)^T * (V)^-1 * (.5 * M + K)
//        math::add(K, .5, M);
//        math::copy(M, K);
//        math::lapack::solve(V, K);
//        math::blas::multiply(1., M, K, 1., D, true, false);
//        math::store(D, "HT_S");

//        for (int r = 0, cc = 0; r < D.nrows; ++r) {
//            for (int c = r; c < D.ncols; ++c) {
//                BETI[cc++] = D.vals[r * D.ncols + c];
//            }
//        }
    }
}

void HeatTransfer::evaluate(step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
    for (size_t i = 0; i < subkernels.size(); ++i) {
        for (size_t e = 0; e < subkernels[i].expressions.node.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                subkernels[i].expressions.node[e]->evaluator->getTime(t) = time.current;
            }
        }
        for (size_t e = 0; e < subkernels[i].expressions.gp.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                subkernels[i].expressions.gp[e]->evaluator->getTime(t) = time.current;
            }
        }
    }
    for (size_t i = 0; i < boundary.size(); ++i) {
        for (size_t j = 0; j < boundary[i].size(); ++j) {
            for (size_t e = 0; e < boundary[i][j].expressions.node.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    boundary[i][j].expressions.node[e]->evaluator->getTime(t) = time.current;
                }
            }
            for (size_t e = 0; e < boundary[i][j].expressions.gp.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    boundary[i][j].expressions.gp[e]->evaluator->getTime(t) = time.current;
                }
            }
        }
    }

    bool run = reset(K, constant.K) | reset(M, constant.M) | reset(f, constant.f) | reset(dirichlet, constant.dirichlet);
    if (run) { assemble(SubKernel::ASSEMBLE); }
    update(K, constant.K);
    update(M, constant.M);
    update(f, constant.f);
    update(dirichlet, constant.dirichlet);
}

void HeatTransfer::getInitialTemperature(Vector_Base<double> *x)
{
    x->setFrom(Results::initialTemperature->data);
}

void HeatTransfer::updateSolution(Vector_Base<double> *x)
{
    Vector_FETI<Vector_Dense, double> *xBEM = dynamic_cast<Vector_FETI<Vector_Dense, double>*>(x);
    #pragma omp parallel for
    for (size_t i = 0; i < bem.size(); ++i) {
        if (bem[i]) {
            int np = info::mesh->domainsSurface->dnodes[i].size();
            double *points = &(info::mesh->domainsSurface->coordinates[i][0].x);
            int ne = info::mesh->domainsSurface->edistribution[i + 1] - info::mesh->domainsSurface->edistribution[i];
            int *elemNodes = info::mesh->domainsSurface->denodes[i].data();
            int ni = info::mesh->domainsSurface->coordinates[i].size() - info::mesh->domainsSurface->dnodes[i].size();
            double *inner = points + 3 * np;


            BEM3DLaplaceEval(xBEM->domains[i].vals + np, np, points, ne, elemNodes, ni, inner, 1, xBEM->domains[i].vals);
        }
    }
    x->storeTo(Results::temperature->data);
    assemble(SubKernel::SOLUTION);
}

}
