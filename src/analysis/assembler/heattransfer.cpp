
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

HeatTransfer::HeatTransfer(HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
    threaded = configuration.solver == HeatTransferLoadStepConfiguration::SOLVER::FETI;
    elementKernels.resize(info::mesh->elements->eintervals.size());
    faceKernels.resize(info::mesh->boundary.size());
    nodeKernels.resize(info::mesh->boundary.size());
    for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
        if (info::mesh->boundary[r]->dimension) {
            faceKernels[r].resize(info::mesh->boundary[r]->eintervals.size());
        }
        nodeKernels[r].resize(info::env::threads);
    }

    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                elementKernels[i].code = info::mesh->elements->eintervals[i].code;
                elementKernels[i].elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
                elementKernels[i].chunks = elementKernels[i].elements / SIMD::size + (elementKernels[i].elements % SIMD::size ? 1 : 0);
            }

            for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
                if (info::mesh->boundary[r]->dimension) {
                    for (esint i = info::mesh->boundary[r]->eintervalsDistribution[d]; i < info::mesh->boundary[r]->eintervalsDistribution[d + 1]; ++i) {
                        faceKernels[r][i].code = info::mesh->boundary[r]->eintervals[i].code;
                        faceKernels[r][i].elements = info::mesh->boundary[r]->eintervals[i].end - info::mesh->boundary[r]->eintervals[i].begin;
                        faceKernels[r][i].chunks = faceKernels[r][i].elements / SIMD::size + (faceKernels[r][i].elements % SIMD::size ? 1 : 0);
                    }
                }
            }
        }
        for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
            nodeKernels[r][t].code = static_cast<int>(Element::CODE::POINT1);
            nodeKernels[r][t].elements = info::mesh->boundary[r]->nodes->datatarray().size(t);
            nodeKernels[r][t].chunks = nodeKernels[r][t].elements / SIMD::size + (nodeKernels[r][t].elements % SIMD::size ? 1 : 0);
        }
    }
}

bool HeatTransfer::analyze(const step::Step &step)
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
        BEM[info::mesh->elements->eintervals[i].domain - info::mesh->domains->offset] = isBEM(i);

        const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
        bool cartesian = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
        bool gptemp = mat->thermal_conductivity.needTemperature();
        esint ebegin = info::mesh->elements->eintervals[i].begin, eend = info::mesh->elements->eintervals[i].end;

        if (info::mesh->dimension == 2) {
            elementKernels[i].thickness.activate(getExpression(i, settings.thickness), info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::thickness->data.data());
        }
        elementKernels[i].material.activate(mat);
        elementKernels[i].conductivity.activate(&mat->thermal_conductivity, &mat->coordinate_system);
        elementKernels[i].heatSource.activate(getExpression(i, configuration.heat_source));
        elementKernels[i].advection.activate(getExpression(i, configuration.translation_motions), settings.sigma);

        bool gpcoo =
                elementKernels[i].conductivity.needCoordinates ||
                elementKernels[i].material.needCoordinates ||
                elementKernels[i].thickness.needCoordinates ||
                elementKernels[i].heatSource.needCoordinates ||
                elementKernels[i].advection.needCoordinates;

        elementKernels[i].coordinates.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, !cartesian || gpcoo);

        auto model = mat->thermal_conductivity.model;
        if (elementKernels[i].conductivity.rotated) {
            switch (mat->thermal_conductivity.model) {
            case ThermalConductivityConfiguration::MODEL::ISOTROPIC:   model = ThermalConductivityConfiguration::MODEL::ISOTROPIC  ; break;
            case ThermalConductivityConfiguration::MODEL::DIAGONAL:    model = ThermalConductivityConfiguration::MODEL::SYMMETRIC  ; break;
            case ThermalConductivityConfiguration::MODEL::SYMMETRIC:   model = ThermalConductivityConfiguration::MODEL::SYMMETRIC  ; break;
            case ThermalConductivityConfiguration::MODEL::ANISOTROPIC: model = ThermalConductivityConfiguration::MODEL::ANISOTROPIC; break;
            }
        }
        elementKernels[i].K.activate(model);
        if (configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
            elementKernels[i].M.activate();
        }

        if (Results::gradient != nullptr) {
            elementKernels[i].gradient.activate(i, Results::gradient);
        }
        if (Results::flux != nullptr) {
            elementKernels[i].flux.activate(i, Results::flux, model);
            elementKernels[i].conductivity.action |= SubKernel::SOLUTION;
//            elementKernels[i].material.action |= SubKernel::SOLUTION;
//            elementKernels[i].advection.action |= SubKernel::SOLUTION;
        }

        elementKernels[i].initTemperature.activate(getExpression(i, settings.initial_temperature));
        elementKernels[i].temperature.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::temperature->data.data(), gptemp);
    }

    for(size_t r = 1; r < info::mesh->boundary.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundary[r];
        if (region->dimension) {
            for(size_t i = 0; i < region->eintervals.size(); ++i) {
                if (info::mesh->dimension == 2) {
                    faceKernels[r][i].thickness.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, Results::thickness->data.data());
                }
                faceKernels[r][i].heatFlow.activate(getExpression(region->name, configuration.heat_flow));
                faceKernels[r][i].heatFlux.activate(getExpression(region->name, configuration.heat_flux));

                auto convection = configuration.convection.find(region->name);
                if (convection != configuration.convection.end()) {
                    faceKernels[r][i].htc.activate(&convection->second.heat_transfer_coefficient);
                    faceKernels[r][i].externalTemperature.activate(&convection->second.external_temperature);
                }
                if (faceKernels[r][i].heatFlow.isactive || faceKernels[r][i].heatFlux.isactive || faceKernels[r][i].htc.isactive) {
                    faceKernels[r][i].externalHeat.activate();
                }

                bool gpcoords =
                        faceKernels[r][i].thickness.needCoordinates ||
                        faceKernels[r][i].heatFlow.needCoordinates ||
                        faceKernels[r][i].heatFlux.needCoordinates ||
                        faceKernels[r][i].htc.needCoordinates ||
                        faceKernels[r][i].externalTemperature.needCoordinates;

                faceKernels[r][i].coordinates.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, gpcoords);
            }
        }
        for(size_t t = 0; t < region->nodes->threads(); ++t) {
            nodeKernels[r][t].coordinates.activate(region->nodes->cbegin(t), region->nodes->cend(), false);
        }
    }

    for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        const BoundaryRegionStore *region = info::mesh->boundary[r];
        for (size_t t = 0; t < region->nodes->threads(); ++t) {
            nodeKernels[r][t].temperature.activate(it->second);
            if (settings.init_temp_respect_bc) {
                nodeKernels[r][t].initialTemperature.activate(region->nodes->cbegin(t), region->nodes->cend(), Results::temperature->data.data(), false);
            }
        }
    }

    assemble(SubKernel::Action::PREPROCESS, step);
    size_t esize = 0;
    std::vector<double> volume(elementKernels.size()), surface(faceKernels.size());
    for (size_t i = 0; i < elementKernels.size(); ++i) {
        esize = std::max(elementKernels[i].esize, esize);
        volume[i] = elementKernels[i].volume;
    }
    for (size_t r = 1; r < faceKernels.size(); ++r) {
        for (size_t i = 0; i < faceKernels[r].size(); ++i) {
            surface[r] += faceKernels[r][i].surface;
        }
    }
    printElementVolume(volume);
    printBoundarySurface(surface);
    for (size_t r = 1; r < faceKernels.size(); ++r) {
        for (size_t i = 0; i < faceKernels[r].size(); ++i) {
            faceKernels[r][i].surface = faceKernels[r][i].externalHeat.area = surface[r];
            faceKernels[r][i].surface = surface[r];
        }
    }
    Results::initialTemperature->data = Results::temperature->data;

    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    if (BEM.front()) {
        eslog::info("  ASSEMBLER                                                                               BEM \n");
    } else {
        eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
        eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize);
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    eslog::info("  PHYSICS ANALYZED                                                                 %8.3f s \n", eslog::time() - start);

    eslog::info(" ============================================================================================= \n");
    return correct;
}

void HeatTransfer::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
    Matrix_FETI<double> *KBEM = dynamic_cast<Matrix_FETI<double>*>(K);
    for (size_t i = 0; i < BEM.size(); ++i) { // when BEM, K is FETI matrix
        if (BEM[i]) {
            BETI[i] = KBEM->domains[i].vals;
            withBEM = true;
        }
    }

    if (withBEM) {
        xBEM.decomposition = KBEM->decomposition;
        xBEM.domains.resize(KBEM->domains.size());
        for (size_t di = 0; di < KBEM->domains.size(); ++di) {
            xBEM.domains[di].resize(KBEM->decomposition->dsize[di]);
        }
    }

    for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        elementKernels[i].Kfiller.activate(i, 1, elementKernels[i].elements, K);
        elementKernels[i].Mfiller.activate(i, 1, elementKernels[i].elements, M);
        elementKernels[i].RHSfiller.activate(i, 1, elementKernels[i].elements, f);
        elementKernels[i].nRHSfiller.activate(i, 1, elementKernels[i].elements, nf);
    }

    for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
                faceKernels[r][i].RHSfiller.activate(r, i, 1, faceKernels[r][i].elements, f);
            }
        }
    }
    for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].dirichlet.activate(r, t, 1, nodeKernels[r][t].elements, dirichlet);
        }
    }

    std::string constants;

    if (K) {
        K->constant = configuration.mode == LoadStepSolverConfiguration::MODE::LINEAR;

        for (size_t i = 0; i < elementKernels.size(); ++i) {
            for (size_t e = 0; e < elementKernels[i].expressions.node.size(); ++e) {
                K->constant &= elementKernels[i].expressions.node[e]->evaluator->isConst();
            }
            for (size_t e = 0; e < elementKernels[i].expressions.gp.size(); ++e) {
                K->constant &= elementKernels[i].expressions.gp[e]->evaluator->isConst();
            }
        }
        if (K->constant) { constants += "K "; }
    }
    if (M) {
        M->constant = true;
        if (M->constant) { constants += "M "; }
    }
    if (f) {
        for (size_t i = 0; i < elementKernels.size(); ++i) {
            for (size_t e = 0; e < elementKernels[i].expressions.node.size(); ++e) {
                f->constant &= elementKernels[i].expressions.node[e]->evaluator->isConst();
            }
            for (size_t e = 0; e < elementKernels[i].expressions.gp.size(); ++e) {
                f->constant &= elementKernels[i].expressions.gp[e]->evaluator->isConst();
            }
        }
        for (size_t i = 0; i < faceKernels.size(); ++i) {
            for (size_t j = 0; j < faceKernels[i].size(); ++j) {
                for (size_t e = 0; e < faceKernels[i][j].expressions.node.size(); ++e) {
                    f->constant &= faceKernels[i][j].expressions.node[e]->evaluator->isConst();
                }
                for (size_t e = 0; e < faceKernels[i][j].expressions.gp.size(); ++e) {
                    f->constant &= faceKernels[i][j].expressions.gp[e]->evaluator->isConst();
                }
            }
        }
        for (size_t i = 0; i < nodeKernels.size(); ++i) {
            for (size_t j = 0; j < nodeKernels[i].size(); ++j) {
                for (size_t e = 0; e < nodeKernels[i][j].expressions.node.size(); ++e) {
                    f->constant &= nodeKernels[i][j].expressions.node[e]->evaluator->isConst();
                }
            }
        }
        if (f->constant) { constants += "f "; }
    }
    if (nf) {
        nf->constant = false;
        if (nf->constant) { constants += "nf "; }
    }
    if (dirichlet) {
        dirichlet->constant = true;
        for (size_t i = 0; i < nodeKernels.size(); ++i) {
            for (size_t j = 0; j < nodeKernels[i].size(); ++j) {
                for (size_t e = 0; e < nodeKernels[i][j].expressions.node.size(); ++e) {
                    dirichlet->constant &= nodeKernels[i][j].expressions.node[e]->evaluator->isConst();
                }
            }
        }
        if (dirichlet->constant) { constants += "dirichlet "; }
    }

    eslog::info("  CONSTANT MATRICES, VECTORS          %*s  [ %s] \n", 50 - constants.size(), " ", constants.c_str());
}

void HeatTransfer::elements(SubKernel::Action action, const step::Step &step, size_t interval)
{
    switch (elementKernels[interval].code) {
    case static_cast<size_t>(Element::CODE::TRIANGLE3): runElement<Element::CODE::TRIANGLE3>(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE6): runElement<Element::CODE::TRIANGLE6>(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::SQUARE4  ): runElement<Element::CODE::SQUARE4  >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::SQUARE8  ): runElement<Element::CODE::SQUARE8  >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::TETRA4   ): runElement<Element::CODE::TETRA4   >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::TETRA10  ): runElement<Element::CODE::TETRA10  >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::PYRAMID5 ): runElement<Element::CODE::PYRAMID5 >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::PYRAMID13): runElement<Element::CODE::PYRAMID13>(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::PRISMA6  ): runElement<Element::CODE::PRISMA6  >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::PRISMA15 ): runElement<Element::CODE::PRISMA15 >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::HEXA8    ): runElement<Element::CODE::HEXA8    >(step, elementKernels[interval], action); break;
    case static_cast<size_t>(Element::CODE::HEXA20   ): runElement<Element::CODE::HEXA20   >(step, elementKernels[interval], action); break;
    }
}

void HeatTransfer::boundary(SubKernel::Action action, const step::Step &step, size_t region, size_t interval)
{
    switch (faceKernels[region][interval].code) {
    case static_cast<size_t>(Element::CODE::LINE2    ): runBoundary<Element::CODE::LINE2    >(step, faceKernels[region][interval], action); break;
    case static_cast<size_t>(Element::CODE::LINE3    ): runBoundary<Element::CODE::LINE3    >(step, faceKernels[region][interval], action); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE3): runBoundary<Element::CODE::TRIANGLE3>(step, faceKernels[region][interval], action); break;
    case static_cast<size_t>(Element::CODE::TRIANGLE6): runBoundary<Element::CODE::TRIANGLE6>(step, faceKernels[region][interval], action); break;
    case static_cast<size_t>(Element::CODE::SQUARE4  ): runBoundary<Element::CODE::SQUARE4  >(step, faceKernels[region][interval], action); break;
    case static_cast<size_t>(Element::CODE::SQUARE8  ): runBoundary<Element::CODE::SQUARE8  >(step, faceKernels[region][interval], action); break;
    }
}

void HeatTransfer::nodes(SubKernel::Action action, const step::Step &step, size_t region, size_t interval)
{
    switch (nodeKernels[region][interval].code) {
    case static_cast<size_t>(Element::CODE::POINT1   ): runNode<Element::CODE::POINT1   >(step, nodeKernels[region][interval], action); break;
    }
}

void HeatTransfer::bem(SubKernel::Action action, size_t domain, double *BETI)
{
    if (action == SubKernel::Action::ASSEMBLE) {
        esint np = info::mesh->domainsSurface->dnodes[domain].size();
        double *points = &(info::mesh->domainsSurface->coordinates[domain][0].x);
        esint ne = info::mesh->domainsSurface->edistribution[domain + 1] - info::mesh->domainsSurface->edistribution[domain];
        esint *elemNodes = info::mesh->domainsSurface->denodes[domain].data();

        double c = elementKernels[domain].conductivity.conductivity->values.get(0, 0).evaluator->evaluate();

        Matrix_Dense<double> K; K.resize(np, np);
        BEM3DLaplace(K.vals, np, points, ne, elemNodes, c);
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

void HeatTransfer::evaluate(const step::Step &step, step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
    for (size_t i = 0; i < elementKernels.size(); ++i) {
        for (size_t e = 0; e < elementKernels[i].expressions.node.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                elementKernels[i].expressions.node[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
                elementKernels[i].expressions.node[e]->evaluator->getTime(t) = time.current;
            }
        }
        for (size_t e = 0; e < elementKernels[i].expressions.gp.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                elementKernels[i].expressions.gp[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
                elementKernels[i].expressions.gp[e]->evaluator->getTime(t) = time.current;
            }
        }

        elementKernels[i].K.isactive = isactive(K);
        elementKernels[i].Kfiller.isactive = isactive(K);
        elementKernels[i].M.isactive = isactive(M);
        elementKernels[i].Mfiller.isactive = isactive(M);
        elementKernels[i].RHSfiller.isactive = isactive(f);
        elementKernels[i].nRHSfiller.isactive = isactive(nf);
        elementKernels[i].timeIntegrationConstantK = time.timeIntegrationConstantK;
        elementKernels[i].timeIntegrationConstantM = time.timeIntegrationConstantM;
    }
    for (size_t i = 0; i < faceKernels.size(); ++i) {
        for (size_t j = 0; j < faceKernels[i].size(); ++j) {
            for (size_t e = 0; e < faceKernels[i][j].expressions.node.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    faceKernels[i][j].expressions.node[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
                    faceKernels[i][j].expressions.node[e]->evaluator->getTime(t) = time.current;
                }
            }
            for (size_t e = 0; e < faceKernels[i][j].expressions.gp.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    faceKernels[i][j].expressions.gp[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
                    faceKernels[i][j].expressions.gp[e]->evaluator->getTime(t) = time.current;
                }
            }
            faceKernels[i][j].RHSfiller.isactive = isactive(f);
        }
    }
    for (size_t i = 0; i < nodeKernels.size(); ++i) {
        for (size_t j = 0; j < nodeKernels[i].size(); ++j) {
            for (size_t e = 0; e < nodeKernels[i][j].expressions.node.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    nodeKernels[i][j].expressions.node[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
                    nodeKernels[i][j].expressions.node[e]->evaluator->getTime(t) = time.current;
                }
            }
            nodeKernels[i][j].dirichlet.isactive = isactive(dirichlet);
        }
    }

    reset(K, M, f, nf, dirichlet);
    assemble(SubKernel::ASSEMBLE, step);
    update(K, M, f, nf, dirichlet);
}

void HeatTransfer::getInitialTemperature(Vector_Base<double> *x)
{
    x->setFrom(Results::initialTemperature->data);
}

void HeatTransfer::updateSolution(const step::Step &step, Vector_Distributed<Vector_Dense, double> *x)
{
    if (withBEM) {
        x->copyTo(&xBEM);
    }
    #pragma omp parallel for
    for (size_t i = 0; i < BEM.size(); ++i) {
        if (BEM[i]) {
            esint np = info::mesh->domainsSurface->dnodes[i].size();
            double *points = &(info::mesh->domainsSurface->coordinates[i][0].x);
            esint ne = info::mesh->domainsSurface->edistribution[i + 1] - info::mesh->domainsSurface->edistribution[i];
            esint *elemNodes = info::mesh->domainsSurface->denodes[i].data();
            esint ni = info::mesh->domainsSurface->coordinates[i].size() - info::mesh->domainsSurface->dnodes[i].size();
            double *inner = points + 3 * np;

            double c = elementKernels[i].conductivity.conductivity->values.get(0, 0).evaluator->evaluate();
            BEM3DLaplaceEval(xBEM.domains[i].vals + np, np, points, ne, elemNodes, ni, inner, c, xBEM.domains[i].vals);
        }
    }
    if (withBEM) {
        xBEM.copyTo(x);
    }
    x->storeTo(Results::temperature->data);
    assemble(SubKernel::SOLUTION, step);
}

}
