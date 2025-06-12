
#include "structuralmechanics.h"
#include "assembler.hpp"

#include "analysis/assembler/structuralmechanics/element.h"
#include "analysis/math/matrix_feti.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "math/math.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/domainsurfacestore.h"
#include "wrappers/bem/w.bem.h"

#include <numeric>
#include <algorithm>

#include <iostream>

namespace espreso {

NodeData* StructuralMechanics::Results::thickness = nullptr;
NodeData* StructuralMechanics::Results::normal = nullptr;
NodeData* StructuralMechanics::Results::initialVelocity = nullptr;

ElementData* StructuralMechanics::Results::principalStress = nullptr;
ElementData* StructuralMechanics::Results::principalStrain = nullptr;
ElementData* StructuralMechanics::Results::componentStress = nullptr;
ElementData* StructuralMechanics::Results::componentStrain = nullptr;
ElementData* StructuralMechanics::Results::vonMisesStress = nullptr;
ElementData* StructuralMechanics::Results::vonMisesStrain = nullptr;

NodeData* StructuralMechanics::Results::principalStressAvg = nullptr;
NodeData* StructuralMechanics::Results::principalStrainAvg = nullptr;
NodeData* StructuralMechanics::Results::componentStressAvg = nullptr;
NodeData* StructuralMechanics::Results::componentStrainAvg = nullptr;
NodeData* StructuralMechanics::Results::vonMisesStressAvg = nullptr;
NodeData* StructuralMechanics::Results::vonMisesStrainAvg = nullptr;
//ElementData* StructuralMechanics::Results::isPlastized = nullptr;

NodeData* StructuralMechanics::Results::displacement = nullptr;

NodeData* StructuralMechanics::Results::cosDisplacement = nullptr;
NodeData* StructuralMechanics::Results::sinDisplacement = nullptr;
NodeData* StructuralMechanics::Results::displacementAmplitude = nullptr;
NodeData* StructuralMechanics::Results::phase = nullptr;
NodeData* StructuralMechanics::Results::velocity = nullptr;
NodeData* StructuralMechanics::Results::velocityAmplitude = nullptr;
NodeData* StructuralMechanics::Results::acceleration = nullptr;
NodeData* StructuralMechanics::Results::accelerationAmplitude = nullptr;

NodeData* StructuralMechanics::Results::force = nullptr;
NodeData* StructuralMechanics::Results::reactionForce = nullptr;
NodeData* StructuralMechanics::Results::contact_force = nullptr;

std::vector<double> MaterialStructuralMechanics::invCp;
std::vector<double> MaterialStructuralMechanics::alpha;

StructuralMechanics::StructuralMechanics(StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
    threaded = configuration.solver == StructuralMechanicsLoadStepConfiguration::SOLVER::FETI;
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

int StructuralMechanics::postProcessSolverSize()
{
    return info::ecf->output.results_selection.stress && info::ecf->output.results_selection.global_average ? 13 : 0;
}

bool StructuralMechanics::analyze(const step::Step &step)
{
    double start = eslog::time();
    eslog::info("\n ============================================================================================= \n");

    validateRegionSettings("MATERIAL", settings.material_set);
//    validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
    validateRegionSettings("INITIAL VELOCITY", settings.initial_velocity);
    validateRegionSettings("THICKNESS", settings.thickness);

    if (Results::thickness == nullptr && info::mesh->dimension == 2) {
        Results::thickness = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "THICKNESS", step::TYPE::TIME, info::ecf->output.results_selection.thickness);
    }
    if (Results::initialVelocity== nullptr) {
        Results::initialVelocity = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "INITIAL_VELOCITY", step::TYPE::TIME, info::ecf->output.results_selection.initial_values);
    }

    bool correct = true;

    if (configuration.type == StructuralMechanicsLoadStepConfiguration::TYPE::HARMONIC) {
        if (Results::cosDisplacement == nullptr) {
            Results::cosDisplacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT_COS", step::TYPE::FREQUENCY, info::ecf->output.results_selection.displacement);
        }
        if (Results::sinDisplacement == nullptr) {
            Results::sinDisplacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT_SIN", step::TYPE::FREQUENCY, info::ecf->output.results_selection.displacement);
        }
        if (Results::displacementAmplitude == nullptr) {
            Results::displacementAmplitude = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::SCALAR, "DISPLACEMENT_AMPLITUDE", step::TYPE::FREQUENCY, info::ecf->output.results_selection.displacement);
        }
        if (Results::phase == nullptr) {
            Results::phase = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::SCALAR, "PHASE", step::TYPE::FREQUENCY, info::ecf->output.results_selection.displacement);
        }

        if (Results::velocity == nullptr) {
            Results::velocity = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "VELOCITY", step::TYPE::FREQUENCY, info::ecf->output.results_selection.velocity);
        }
        if (Results::velocityAmplitude == nullptr) {
            Results::velocityAmplitude = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::SCALAR, "VELOCITY_AMPLITUDE", step::TYPE::FREQUENCY, info::ecf->output.results_selection.velocity);
        }
        if (Results::acceleration == nullptr) {
            Results::acceleration = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "ACCELERATION", step::TYPE::FREQUENCY, info::ecf->output.results_selection.acceleration);
        }
        if (Results::accelerationAmplitude == nullptr) {
            Results::accelerationAmplitude = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::SCALAR, "ACCELERATION_AMPLITUDE", step::TYPE::FREQUENCY, info::ecf->output.results_selection.acceleration);
        }
    } else {
        if (Results::displacement == nullptr) {
            Results::displacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT", step::TYPE::TIME, info::ecf->output.results_selection.displacement);
        }
        if (Results::velocity == nullptr) {
            Results::velocity = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "VELOCITY", step::TYPE::TIME, info::ecf->output.results_selection.velocity);
        }
        if (Results::acceleration == nullptr) {
            Results::acceleration = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "ACCELERATION", step::TYPE::TIME, info::ecf->output.results_selection.acceleration);
        }
        int tensor = info::mesh->dimension == 3 ? 6 : 3;
        if (info::ecf->output.results_selection.stress && Results::principalStress == nullptr) {
            Results::principalStress = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::NUMBERED   , "PRINCIPAL_STRESS");
            Results::principalStrain = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::NUMBERED   , "PRINCIPAL_STRAIN");
            Results::componentStress = info::mesh->elements->appendData(               tensor, NamedData::DataType::TENSOR_SYMM, "COMPONENT_STRESS");
            Results::componentStrain = info::mesh->elements->appendData(               tensor, NamedData::DataType::TENSOR_SYMM, "COMPONENT_STRAIN");
            Results::vonMisesStress  = info::mesh->elements->appendData(                    1, NamedData::DataType::SCALAR     , "VON_MISES_STRESS");
            Results::vonMisesStrain  = info::mesh->elements->appendData(                    1, NamedData::DataType::SCALAR     , "VON_MISES_STRAIN");
        }
        if (info::ecf->output.results_selection.stress && Results::principalStressAvg == nullptr) {
            Results::principalStressAvg = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::NUMBERED   , "AVG_PRINCIPAL_STRESS");
            Results::principalStrainAvg = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::NUMBERED   , "AVG_PRINCIPAL_STRAIN");
            Results::componentStressAvg = info::mesh->nodes->appendData(               tensor, NamedData::DataType::TENSOR_SYMM, "AVG_COMPONENT_STRESS");
            Results::componentStrainAvg = info::mesh->nodes->appendData(               tensor, NamedData::DataType::TENSOR_SYMM, "AVG_COMPONENT_STRAIN");
            Results::vonMisesStressAvg  = info::mesh->nodes->appendData(                    1, NamedData::DataType::SCALAR     , "AVG_VON_MISES_STRESS");
            Results::vonMisesStrainAvg  = info::mesh->nodes->appendData(                    1, NamedData::DataType::SCALAR     , "AVG_VON_MISES_STRAIN");
        }
        if (Results::reactionForce == nullptr) {
            Results::reactionForce = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "REACTION_FORCES", step::TYPE::TIME, info::ecf->output.results_selection.reactions);
        }
        if (Results::force == nullptr) {
            Results::force = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FORCE", step::TYPE::TIME, info::ecf->output.results_selection.force);
        }
        if (Results::contact_force == nullptr) {
            Results::contact_force = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "CONTACT_FORCE");
        }

//        for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
//            if (info::mesh->materials[i]->material_model == MaterialConfiguration::MATERIAL_MODEL::PLASTICITY) {
//                if (Results::isPlastized == nullptr) {
//                    Results::isPlastized = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "IS_PLASTIZED");
//                }
//            }
//        }
    }

    bool withNormals = false;
    for(size_t r = 1; r < info::mesh->boundary.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundary[r];
        if (StringCompare::caseInsensitivePreffix("CONTACT", region->name)) {
            withNormals = true;
        }
    }

    if (withNormals) {
        if (Results::normal == nullptr) {
            Results::normal = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "NORMAL", step::TYPE::TIME, info::ecf->output.results_selection.normal);
        }
        faceMultiplicity.resize(info::mesh->nodes->size);
        for(size_t r = 1; r < info::mesh->boundary.size(); ++r) {
            const BoundaryRegionStore *region = info::mesh->boundary[r];
            if (StringCompare::caseInsensitivePreffix("CONTACT", region->name)) {
                for (auto face = region->elements->cbegin(); face != region->elements->cend(); ++face) {
                    for (auto n = face->begin(); n != face->end(); ++n) {
                        faceMultiplicity[*n] += 1;
                    }
                }
            }
        }
        NodeData::synchronize(faceMultiplicity, 1);
        for (size_t i = 0; i < faceMultiplicity.size(); ++i) {
            if (faceMultiplicity[i] != 0) {
                faceMultiplicity[i] = 1 / faceMultiplicity[i];
            }
        }
    }

    if (info::ecf->output.results_selection.stress) {
        nodeMultiplicity.resize(info::mesh->nodes->size);
        for (auto enodes = info::mesh->elements->nodes->cbegin(); enodes != info::mesh->elements->nodes->cend(); ++enodes) {
            for (auto n = enodes->begin(); n != enodes->end(); ++n) {
                nodeMultiplicity[*n] += 1;
            }
        }
        NodeData::synchronize(nodeMultiplicity, 1);
        for (size_t i = 0; i < nodeMultiplicity.size(); ++i) {
            nodeMultiplicity[i] = 1 / nodeMultiplicity[i];
        }
    }

    if (settings.initial_temperature.size()) {
        correct &= checkElementParameter("INITIAL TEMPERATURE", settings.initial_temperature);
    }


    eslog::info(" ============================================================================================= \n");

    if (configuration.displacement.size()) {
        correct &= checkBoundaryParameter("FIXED DISPLACEMENT", configuration.displacement);
//        generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 0, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][0][s] = value; });
//        generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 1, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][1][s] = value; });
//        generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 2, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][2][s] = value; });
    }

    if (true) { // add check
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

            if (info::mesh->dimension == 2) {
                switch (settings.element_behaviour) {
                case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
                    eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRAIN \n");
                    break;
                case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
                    eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRESS \n");
                    break;
                case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
                    eslog::info("     ELEMENT BEHAVIOR:                                            PLANE STRESS WITH THICKNESS \n");
                    correct &= checkElementParameter("THICKNESS", settings.thickness);
                    break;
                case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
                    eslog::info("     ELEMENT BEHAVIOR:                                                           AXISYMMETRIC \n");
                    correct &= checkElementParameter("THICKNESS", settings.thickness);
                    break;
                }
                eslog::info("                                                                                               \n");
            }

            correct &= checkExpression("DENSITY", mat->density);
            correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
            eslog::info("                                                                                               \n");

            switch (mat->elasticity_properties.material_model) {
            case ElasticityPropertiesConfiguration::MATERIAL_MODEL::KIRCHHOFF:
                eslog::info("       MATERIAL MODEL:                                                               KIRCHHOFF \n");
                eslog::info("                                                                                               \n");
                break;
            case ElasticityPropertiesConfiguration::MATERIAL_MODEL::NEO_HOOKEAN_CMP:
                eslog::info("       MATERIAL MODEL:                                                COMPRESSIBLE NEO-HOOKEAN \n");
                eslog::info("                                                                                               \n");
                break;
            case ElasticityPropertiesConfiguration::MATERIAL_MODEL::BONET_WOOD:
                eslog::info("       MATERIAL MODEL:                                                              BONET-WOOD \n");
                eslog::info("                                                                                               \n");
                break;
            }
            correct &= checkExpression("EX", mat->elasticity_properties.young_modulus.get(0, 0));
            correct &= checkExpression("MI", mat->elasticity_properties.poisson_ratio.get(0, 0));
//                correct &= checkExpression("INITIAL_YIELD_STRESS", mat->plasticity_properties.initial_yield_stress);
//                correct &= checkExpression("ISOTROPIC_HARDENING", mat->plasticity_properties.isotropic_hardening);
//                correct &= checkExpression("KINEMATIC_HARDENING", mat->plasticity_properties.kinematic_hardening);
//                eslog::info("                                                                                               \n");
                /* no break */

//            case MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
//                switch (mat->linear_elastic_properties.model) {
//                case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
//                    eslog::info(" LINEAR ELASTIC MODEL:                                                              ISOTROPIC \n");
//                    correct &= checkExpression("EX", mat->linear_elastic_properties.young_modulus.get(0, 0));
//                    correct &= checkExpression("MI", mat->linear_elastic_properties.poisson_ratio.get(0, 0));
//                    break;
//                case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
//                    eslog::info("                MODEL:                                                            ORTHOTROPIC \n");
//                    correct &= checkExpression("EX", mat->linear_elastic_properties.young_modulus.get(0, 0));
//                    correct &= checkExpression("EY", mat->linear_elastic_properties.young_modulus.get(1, 1));
//                    correct &= checkExpression("EZ", mat->linear_elastic_properties.young_modulus.get(2, 2));
//                    correct &= checkExpression("MIXY", mat->linear_elastic_properties.poisson_ratio.get(0, 0));
//                    correct &= checkExpression("MIXZ", mat->linear_elastic_properties.poisson_ratio.get(1, 1));
//                    correct &= checkExpression("MIYZ", mat->linear_elastic_properties.poisson_ratio.get(2, 2));
//                    correct &= checkExpression("GXY", mat->linear_elastic_properties.shear_modulus.get(0, 0));
//                    correct &= checkExpression("GXZ", mat->linear_elastic_properties.shear_modulus.get(1, 1));
//                    correct &= checkExpression("GYZ", mat->linear_elastic_properties.shear_modulus.get(2, 2));
//                    break;
//                case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
//                    eslog::info("                MODEL:                                                            ANISOTROPIC \n");
//                    break;
//                }
//                break;
//            case MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
//                eslog::info("                                                                                HYPER ELASTIC \n");
//                break;
//            }
            eslog::info("                                                                                               \n");
        }

        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
        printMaterials(settings.material_set);
        eslog::info(" ============================================================================================= \n");
    }

    if (configuration.acceleration.size()) {
        correct &= checkElementParameter("ACCELERATION", configuration.acceleration);
    }
    if (configuration.angular_velocity.size()) {
        switch (info::mesh->dimension) {
        case 2:
            switch (settings.element_behaviour) {
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
                correct &= checkElementParameter("ANGULAR_VELOCITY.Z", configuration.angular_velocity, 2); break;
            case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
                correct &= checkElementParameter("ANGULAR_VELOCITY.Y", configuration.angular_velocity, 1); break;
            } break;
        case 3:
            correct &= checkElementParameter("ANGULAR_VELOCITY", configuration.angular_velocity); break;
        }
    }

    for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundary[r];
        for (size_t i = 0; i < region->data.size(); ++i) {
            if (StringCompare::caseInsensitiveEq(region->data[i]->name, "FORCE")) {
                eslog::info("  %s%*s \n", region->name.c_str(), 91 - region->name.size(), "");
                eslog::info("   %25s:  %62s \n", "FORCE", "EXTERNAL");
            }
        }
    }
    for (auto pressure = configuration.pressure.begin(); pressure != configuration.pressure.end(); ++pressure) {
        eslog::info("  %s: %*s\n", pressure->first.c_str(), 90 - pressure->first.size(), " ");
        correct &= checkExpression("PRESSURE", pressure->second.pressure);
        correct &= checkExpression("DIRECTION", pressure->second.direction);
    }
    if (configuration.normal_pressure.size()) {
        correct &= checkBoundaryParameter("NORMAL_PRESSURE", configuration.normal_pressure);
    }
    if (configuration.force.size()) {
        correct &= checkBoundaryParameter("FORCE", configuration.force);
    }
    if (configuration.harmonic_force.size()) {
        correct &= checkBoundaryParameter("HARMONIC_FORCE", configuration.harmonic_force);
    }

    for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        BEM[info::mesh->elements->eintervals[i].domain - info::mesh->domains->offset] = isBEM(i);
        const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
        bool cartesian = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
        bool gpcoo = mat->elasticity_properties.needCoordinates() || getExpression(i, configuration.angular_velocity);
        gpcoo |= settings.element_behaviour == StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC;
        gpcoo |= mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
//        bool gptemp = mat->linear_elastic_properties.needTemperature();
        esint ebegin = info::mesh->elements->eintervals[i].begin, eend = info::mesh->elements->eintervals[i].end;

        if (info::mesh->dimension == 2) {
            elementKernels[i].thickness.activate(getExpression(i, settings.thickness), info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::thickness->data.data());
        }

        elementKernels[i].coordinates.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, !cartesian || gpcoo);
        elementKernels[i].displacement.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::displacement->data.data());
        elementKernels[i].material.activate(mat, settings.element_behaviour, configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR && !configuration.corotation);

        elementKernels[i].matrixElasticity.activate(settings.element_behaviour, configuration.mode == LoadStepSolverConfiguration::MODE::NONLINEAR && !configuration.corotation);
        if (configuration.corotation) {
            elementKernels[i].matrixCorotation.activate();
        }

        if (configuration.type != LoadStepSolverConfiguration::TYPE::STEADY_STATE) {
            elementKernels[i].M.activate();
        }
//        if (Results::avgStress) {
//            elementKernels[i].postM.activate();
//            elementKernels[i].postM.action = SubKernel::SOLUTION;
//        }
        elementKernels[i].acceleration.activate(getExpression(i, configuration.acceleration), settings.element_behaviour);
        elementKernels[i].angularVelocity.activate(getExpression(i, configuration.angular_velocity), settings.element_behaviour);
        if (Results::principalStress) {
            elementKernels[i].displacement.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::displacement->data.data());
            elementKernels[i].stress.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, i, nodeMultiplicity.data());
        }

        elementKernels[i].initVelocity.activate(getExpression(i, settings.initial_velocity));
        elementKernels[i].velocity.activate(info::mesh->elements->nodes->cbegin() + ebegin, info::mesh->elements->nodes->cbegin() + eend, Results::initialVelocity->data.data());

        if (info::ecf->output.print_eigen_values > 2) {
            elementKernels[i].print.activate();
        }
    }

    for (auto wall = configuration.fixed_wall.begin(); wall != configuration.fixed_wall.end(); ++wall) {
        eslog::info("  %s: %*s\n", wall->first.c_str(), 90 - wall->first.size(), " ");
        correct &= checkExpression("POINT", wall->second.point);
        correct &= checkExpression("NORMAL", wall->second.normal);
    }

    for(size_t r = 1; r < info::mesh->boundary.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundary[r];
        BoundaryElementData *force = nullptr;
        for (size_t i = 0; i < region->data.size(); ++i) {
            if (StringCompare::caseInsensitiveEq(region->data[i]->name, "FORCE")) {
                force = region->data[i];
            }
        }
        if (info::mesh->boundary[r]->dimension) {
            for (size_t i = 0; i < info::mesh->boundary[r]->eintervals.size(); ++i) {
                faceKernels[r][i].coordinates.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, settings.element_behaviour == StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC);
                faceKernels[r][i].normalPressure.activate(getExpression(region->name, configuration.normal_pressure), settings.element_behaviour);
                faceKernels[r][i].displacement.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, Results::displacement->data.data());
                if (StringCompare::caseInsensitivePreffix("CONTACT", region->name)) {
                    faceKernels[r][i].normal.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cbegin() + region->eintervals[i].end, Results::normal->data.data(), faceMultiplicity.data());
                }
                auto pressure = configuration.pressure.find(region->name);
                if (pressure != configuration.pressure.end()) {
                    faceKernels[r][i].pressure.activate(pressure->second.pressure, pressure->second.direction, settings.element_behaviour);
                }
                if (force) {
                    faceKernels[r][i].force.activate(force->data.data(), region->eintervals[i].begin, region->eintervals[i].end);
                }
            }
        }
        for (size_t t = 0; t < region->nodes->threads(); ++t) {
            nodeKernels[r][t].force.activate(getExpression(region->name, configuration.force));
            nodeKernels[r][t].harmonicForce.activate(getExpression(region->name, configuration.harmonic_force), settings.element_behaviour);
            nodeKernels[r][t].coordinates.activate(region->nodes->cbegin(t), region->nodes->cend(), false);
            if (force && info::mesh->boundary[r]->dimension == 0) {
                nodeKernels[r][t].force.activate(force->data.data(), 0, 0, region->nodes->cbegin(t), region->nodes->cend(t));
            }
        }
    }

    for (auto it = configuration.displacement.begin(); it != configuration.displacement.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].displacement.activate(it->second);
        }
    }

    assemble(SubKernel::PREPROCESS, step);
    if (Results::normal) {
        Results::normal->synchronize();
    }
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
        for (size_t i = 0; i < faceKernels[r].size(); ++i) {
            faceKernels[r][i].surface = surface[r];
        }
    }
    printElementVolume(volume);
    printBoundarySurface(surface);

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

void StructuralMechanics::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet, Matrix_Base<double> *postM, Vector_Base<double> *postB)
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
        if (!BEM[info::mesh->elements->eintervals[i].domain - info::mesh->domains->offset]) {
            elementKernels[i].Kfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, K);
            elementKernels[i].Mfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, M);
            elementKernels[i].postMfiller.activate(i, 1, elementKernels[i].elements, postM);
            elementKernels[i].postMfiller.action = SubKernel::SOLUTION;
        }

        elementKernels[i].reRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, f);
        elementKernels[i].reNRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, nf);
        elementKernels[i].postBfiller.activate(i, 1, elementKernels[i].elements, postB);
        elementKernels[i].postBfiller.action = SubKernel::SOLUTION;
    }

    for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
                faceKernels[r][i].reRHSfiller.activate(r, i, info::mesh->dimension, faceKernels[r][i].elements, f);
            }
        }
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].reRHSfiller.activate(r, t, info::mesh->dimension, nodeKernels[r][t].elements, f);
        }
    }
    for (auto it = configuration.displacement.begin(); it != configuration.displacement.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        for (size_t t = 0; t < info::mesh->boundary[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].reDirichlet.activate(r, t, info::mesh->dimension, nodeKernels[r][t].elements, dirichlet);
        }
    }
}

void StructuralMechanics::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet)
{
    for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        elementKernels[i].Kfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, K);
        elementKernels[i].Mfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, M);
        elementKernels[i].Cfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, C);

        elementKernels[i].reRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, ref);
        elementKernels[i].imRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, imf);
        elementKernels[i].reNRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, renf);
        elementKernels[i].imNRHSfiller.activate(i, info::mesh->dimension, elementKernels[i].elements, imnf);
    }

    for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
                faceKernels[r][i].reRHSfiller.activate(r, i, info::mesh->dimension, faceKernels[r][i].elements, ref);
                faceKernels[r][i].imRHSfiller.activate(r, i, info::mesh->dimension, faceKernels[r][i].elements, imf);
            }
        }
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].reRHSfiller.activate(r, t, info::mesh->dimension, nodeKernels[r][t].elements, ref);
            nodeKernels[r][t].imRHSfiller.activate(r, t, info::mesh->dimension, nodeKernels[r][t].elements, imf);
        }
    }
    for (auto it = configuration.displacement.begin(); it != configuration.displacement.end(); ++it) {
        size_t r = info::mesh->bregionIndex(it->first);
        for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
            nodeKernels[r][t].reDirichlet.activate(r, t, info::mesh->dimension, nodeKernels[r][t].elements, reDirichlet); // imDirichlet is 0
        }
    }
}

void StructuralMechanics::evaluate(const step::Step &step, const step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
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

        elementKernels[i].Kfiller.isactive = isactive(K);
        elementKernels[i].Mfiller.isactive = isactive(M);
        elementKernels[i].reRHSfiller.isactive = isactive(f);
        elementKernels[i].reNRHSfiller.isactive = isactive(nf);
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
            faceKernels[i][j].reRHSfiller.isactive = isactive(f);
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
            nodeKernels[i][j].reRHSfiller.isactive = isactive(f);
            nodeKernels[i][j].reDirichlet.isactive = isactive(dirichlet);
        }
    }

    for (auto wall = configuration.fixed_wall.begin(); wall != configuration.fixed_wall.end(); ++wall) {
        wall->second.normal.x.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.normal.y.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.normal.z.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.point.x.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.point.y.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.point.z.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        wall->second.normal.x.evaluator->getTime() = time.current;
        wall->second.normal.y.evaluator->getTime() = time.current;
        wall->second.normal.z.evaluator->getTime() = time.current;
        wall->second.point.x.evaluator->getTime() = time.current;
        wall->second.point.y.evaluator->getTime() = time.current;
        wall->second.point.z.evaluator->getTime() = time.current;
    }

    for (auto sphere = configuration.fixed_sphere.begin(); sphere != configuration.fixed_sphere.end(); ++sphere) {
        sphere->second.center.x.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        sphere->second.center.y.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        sphere->second.center.z.evaluator->getSubstep() = (step.substep + 1) / (double)step.substeps;
        sphere->second.center.x.evaluator->getTime() = time.current;
        sphere->second.center.y.evaluator->getTime() = time.current;
        sphere->second.center.z.evaluator->getTime() = time.current;
    }

    if (Results::normal) {
        std::fill(Results::normal->data.begin(), Results::normal->data.end(), 0);
    }
    reset(K, M, f, nf, dirichlet);
    assemble(SubKernel::ASSEMBLE, step);
    update(K, M, f, nf, dirichlet);
    if (Results::normal) {
        Results::normal->synchronize();
    }
    if (Results::reactionForce && nf) {
        nf->storeTo(Results::reactionForce->data);
    }
    if (Results::force && f) {
        f->storeTo(Results::force->data);
    }

    if (info::ecf->output.print_eigen_values > 1) {
        if (M && M->updated) {
            M->printEigenValues("M[ASM]", 8);
        }
        if (K && K->updated) {
            K->printEigenValues("K[ASM]", 8);
        }
    }
}

void StructuralMechanics::evaluate(const step::Step &step, const step::Frequency &freq, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *reDirichlet, Vector_Base<double> *imDirichlet)
{
    for (size_t i = 0; i < elementKernels.size(); ++i) {
        for (size_t e = 0; e < elementKernels[i].expressions.node.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                elementKernels[i].expressions.node[e]->evaluator->getFrequency(t) = freq.current;
            }
        }
        for (size_t e = 0; e < elementKernels[i].expressions.gp.size(); ++e) {
            #pragma omp parallel for
            for (int t = 0; t < info::env::threads; ++t) {
                elementKernels[i].expressions.gp[e]->evaluator->getFrequency(t) = freq.current;
            }
        }
        elementKernels[i].Kfiller.isactive = isactive(K);
        elementKernels[i].Mfiller.isactive = isactive(M);
        elementKernels[i].reRHSfiller.isactive = isactive(ref);
        elementKernels[i].imRHSfiller.isactive = isactive(imf);
        elementKernels[i].reNRHSfiller.isactive = isactive(renf);
        elementKernels[i].imNRHSfiller.isactive = isactive(imnf);
    }
    for (size_t i = 0; i < faceKernels.size(); ++i) {
        for (size_t j = 0; j < faceKernels[i].size(); ++j) {
            for (size_t e = 0; e < faceKernels[i][j].expressions.node.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    faceKernels[i][j].expressions.node[e]->evaluator->getFrequency(t) = freq.current;
                }
            }
            for (size_t e = 0; e < faceKernels[i][j].expressions.gp.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    faceKernels[i][j].expressions.gp[e]->evaluator->getFrequency(t) = freq.current;
                }
            }
            faceKernels[i][j].reRHSfiller.isactive = isactive(ref);
            faceKernels[i][j].imRHSfiller.isactive = isactive(imf);
        }
    }

    for (size_t i = 0; i < nodeKernels.size(); ++i) {
        for (size_t j = 0; j < nodeKernels[i].size(); ++j) {
            for (size_t e = 0; e < nodeKernels[i][j].expressions.node.size(); ++e) {
                #pragma omp parallel for
                for (int t = 0; t < info::env::threads; ++t) {
                    nodeKernels[i][j].expressions.node[e]->evaluator->getFrequency(t) = freq.current;
                }
            }
            nodeKernels[i][j].reRHSfiller.isactive = isactive(ref);
            nodeKernels[i][j].imRHSfiller.isactive = isactive(imf);
            nodeKernels[i][j].reDirichlet.isactive = isactive(reDirichlet);
            nodeKernels[i][j].imDirichlet.isactive = isactive(imDirichlet);
        }
    }

    if (Results::normal) {
        std::fill(Results::normal->data.begin(), Results::normal->data.end(), 0);
    }
    reset(K, M, C, ref, imf, renf, imnf, reDirichlet, imDirichlet);
    assemble(SubKernel::ASSEMBLE, step);
    if (Results::normal) {
        Results::normal->synchronize();
    }
    update(K, M, C, ref, imf, renf, imnf, reDirichlet, imDirichlet);
}

void StructuralMechanics::elements(SubKernel::Action action, const step::Step &step, size_t interval)
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

void StructuralMechanics::boundary(SubKernel::Action action, const step::Step &step, size_t region, size_t interval)
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

void StructuralMechanics::nodes(SubKernel::Action action, const step::Step &step, size_t region, size_t interval)
{
    switch (nodeKernels[region][interval].code) {
    case static_cast<size_t>(Element::CODE::POINT1   ): runNode<Element::CODE::POINT1   >(step, nodeKernels[region][interval], action); break;
    }
}

void StructuralMechanics::bem(SubKernel::Action action, size_t domain, double *BETI)
{
    if (action == SubKernel::Action::ASSEMBLE) {
        esint np = info::mesh->domainsSurface->dnodes[domain].size();
        double *points = &(info::mesh->domainsSurface->coordinates[domain][0].x);
        esint ne = info::mesh->domainsSurface->edistribution[domain + 1] - info::mesh->domainsSurface->edistribution[domain];
        esint *elemNodes = info::mesh->domainsSurface->denodes[domain].data();

        double ex = elementKernels[domain].material.configuration->elasticity_properties.young_modulus.get(0, 0).evaluator->evaluate();
        double mu = elementKernels[domain].material.configuration->elasticity_properties.poisson_ratio.get(0, 0).evaluator->evaluate();

        Matrix_Dense<double> K; K.resize(3 * np, 3 * np);
        BEM3DElasticity(K.vals, np, points, ne, elemNodes, ex, mu);

        for (int r = 0, cc = 0; r < K.nrows; ++r) {
            int br = np * (r % 3) + r / 3;
            for (int c = r; c < K.ncols; ++c) {
                int bc = np * (c % 3) + c / 3;
                BETI[cc++] = K.vals[br * K.ncols + bc];
            }
        }
    }
}

void StructuralMechanics::getInitialVelocity(Vector_Base<double> *x)
{
    x->setFrom(Results::initialVelocity->data);
}

void StructuralMechanics::updateSolution(const step::Step &step, Vector_Distributed<Vector_Dense, double> *x, Matrix_Base<double> *M, Vector_Base<double> *B)
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
            double ex = elementKernels[i].material.configuration->elasticity_properties.young_modulus.get(0, 0).evaluator->evaluate();
            double mu = elementKernels[i].material.configuration->elasticity_properties.poisson_ratio.get(0, 0).evaluator->evaluate();
            std::vector<double> xx(xBEM.domains[i].size);
            for (esint p = 0; p < np; ++p) {
                xx[0 * np + p] = xBEM.domains[i].vals[3 * p + 0];
                xx[1 * np + p] = xBEM.domains[i].vals[3 * p + 1];
                xx[2 * np + p] = xBEM.domains[i].vals[3 * p + 2];
            }
            BEM3DElasticityEval(xx.data() + 3 * np, np, points, ne, elemNodes, ni, inner, ex, mu, xx.data());
            for (esint p = 0; p < ni; ++p) {
                xBEM.domains[i].vals[3 * (p + np) + 0] = xx[3 * np + 0 * ni + p];
                xBEM.domains[i].vals[3 * (p + np) + 1] = xx[3 * np + 1 * ni + p];
                xBEM.domains[i].vals[3 * (p + np) + 2] = xx[3 * np + 2 * ni + p];
            }
        }
    }
    if (withBEM) {
        xBEM.copyTo(x);
    }
    x->storeTo(Results::displacement->data);
    reset(M, B);
    if (info::ecf->output.results_selection.stress) {
        std::fill(Results::principalStressAvg->data.begin(), Results::principalStressAvg->data.end(), 0);
        std::fill(Results::principalStrainAvg->data.begin(), Results::principalStrainAvg->data.end(), 0);
        std::fill(Results::componentStressAvg->data.begin(), Results::componentStressAvg->data.end(), 0);
        std::fill(Results::componentStrainAvg->data.begin(), Results::componentStrainAvg->data.end(), 0);
        std::fill(Results::vonMisesStressAvg->data.begin() , Results::vonMisesStressAvg->data.end(), 0);
        std::fill(Results::vonMisesStrainAvg->data.begin() , Results::vonMisesStrainAvg->data.end(), 0);
    }
    assemble(SubKernel::SOLUTION, step);
    if (info::ecf->output.results_selection.stress) {
        Results::principalStressAvg->synchronize();
        Results::principalStrainAvg->synchronize();
        Results::componentStressAvg->synchronize();
        Results::componentStrainAvg->synchronize();
        Results::vonMisesStressAvg ->synchronize();
        Results::vonMisesStrainAvg ->synchronize();
    }
    update(M, B);
}

void StructuralMechanics::updateStress(const step::Step &step, Vector_Base<double> *x)
{
    // TODO
}

void StructuralMechanics::updateSolution(const step::Step &step, Vector_Base<double> *rex, Vector_Base<double> *imx)
{
    rex->storeTo(Results::cosDisplacement->data);
    imx->storeTo(Results::sinDisplacement->data);
    assemble(SubKernel::SOLUTION, step);
}

void StructuralMechanics::nextIteration(const step::Step &step, Vector_Distributed<Vector_Dense, double> *x)
{
    x->storeTo(Results::displacement->data);
    assemble(SubKernel::ITERATION, step);
}

}
