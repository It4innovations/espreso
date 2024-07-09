
#include "acoustic.h"
#include "assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

using namespace espreso;

NodeData* Acoustic::Results::pressure = nullptr;
NodeData* Acoustic::Results::initialPressure = nullptr;

Acoustic::Acoustic(Acoustic *previous, AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{

}

void Acoustic::initParameters()
{
    if (Results::initialPressure == nullptr) {
        Results::initialPressure = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_ACOUSTIC_PRESSURE");
//        Variable::list.node["INITIAL_ACOUSTIC_PRESSURE"] = new OutputVariable(Results::initialPressure, 0, 1);
    }
    if (Results::pressure == nullptr) {
        Results::pressure = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "ACOUSTIC_PRESSURE");
//        Variable::list.node["ACOUSTIC_PRESSURE"] = new OutputVariable(Results::pressure, 0, 1);
    }
}

void Acoustic::analyze(const step::Step &step)
{
    double start = eslog::time();
    eslog::info("\n ============================================================================================= \n");
    bool correct = true;

//    validateRegionSettings("MATERIAL", settings.material_set);
//    validateRegionSettings("THICKNESS", settings.thickness);
//
//    initParameters();
//
//    baseFunction(*this);
//    elementCoordinates(*this);
//    elementIntegration(*this);
//
//    if (configuration.acoustic_pressure.size()) {
//        correct &= examineBoundaryParameter("FIXED ACOUSTIC PRESSURE ON BOUNDARIES", configuration.acoustic_pressure, pressure.node.externalValues);
//        fromExpression(*this, pressure.node, pressure.node.externalValues);
//    }
//
//    if (step::step.loadstep == 0) {
//        eslog::info("\n  MATERIALS                                                                                    \n");
//        eslog::info(" --------------------------------------------------------------------------------------------- \n");
//
//        for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
//            eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
//            MaterialConfiguration *mat = info::mesh->materials[i];
//            correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValues, 0);
//            correct &= examineMaterialParameter(mat->name, "SPEED_OF_SOUND", mat->speed_of_sound, material.speed_of_sound.externalValues, 0);
//        }
//
//        fromExpression(*this, material.density, material.density.externalValues);
//        fromExpression(*this, material.speed_of_sound, material.speed_of_sound.externalValues);
//
//        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
//
//        printMaterials(settings.material_set);
//
//        eslog::info(" ============================================================================================= \n");
//    }
//
//
//    stiffness(*this);
//    mass(*this);
//    // boundary conditions have to added according to boundary settings below
////    acousticBoundaryMass(*this);
//
//    if (configuration.normal_acceleration.size()) {
//        examineBoundaryParameter("NORMAL ACCELERATION", configuration.normal_acceleration, normalAcceleration.gp.externalValues);
//        fromExpression(*this, normalAcceleration.gp, normalAcceleration.gp.externalValues);
//    }
//    if (configuration.acceleration.size()) {
//        correct &= examineBoundaryParameter("ACCELERATION.X", configuration.acceleration, acceleration.gp.externalValues, 0);
//        correct &= examineBoundaryParameter("ACCELERATION.Y", configuration.acceleration, acceleration.gp.externalValues, 1);
//        if (info::mesh->dimension == 3) {
//            correct &= examineBoundaryParameter("ACCELERATION.Z", configuration.acceleration,  acceleration.gp.externalValues, 2);
//        }
//        fromExpression(*this, acceleration.gp, acceleration.gp.externalValues);
//    }
//    if (configuration.impedance.size()) {
//        examineBoundaryParameter("IMPEDANCE", configuration.impedance, impedance.gp.externalValues);
//        fromExpression(*this, impedance.gp, impedance.gp.externalValues);
//    }
//    if (configuration.monopole_source.size()) {
//        correct &= examineElementParameter("MONOPOLE DOMAIN SOURCE", configuration.monopole_source, monopoleSource.gp.externalValues);
//        fromExpression(*this, monopoleSource.gp, monopoleSource.gp.externalValues);
//    }
//    if (configuration.dipole_source.size()) {
//        correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.X", configuration.dipole_source, dipoleSource.gp.externalValues, 0);
//        correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.Y", configuration.dipole_source, dipoleSource.gp.externalValues, 1);
//        if (info::mesh->dimension == 3) {
//            correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.Z", configuration.dipole_source, dipoleSource.gp.externalValues, 2);
//        }
//        fromExpression(*this, dipoleSource.gp, dipoleSource.gp.externalValues);
//    }
//    if (configuration.point_source.size()) {
//        correct &= examineBoundaryParameter("POINT SOURCES", configuration.point_source, pointSource.node.externalValues);
//        fromExpression(*this, pointSource.node, pointSource.node.externalValues);
//    }
//
//    integration.weight.name = "integration.weight";
//    integration.N.name = "integration.N";
//    integration.dN.name = "integration.dN";
//    integration.dND.name = "integration.dND";
//    integration.jacobiDeterminant.name = "integration.jacobiDeterminant";
//    integration.jacobiInversion.name = "integration.jacobiInversion";
//    elements.monopole.name = "elements.monopole";
//    elements.dipole.name = "elements.dipole";
//    material.density.name = "material.density";
//
//    RHS(*this);

    eslog::info(" ============================================================================================= \n");
    if (correct) {
        eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
    } else {
        eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
    }
    eslog::info(" ============================================================================================= \n");
}

void Acoustic::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *dirichlet)
{

}

void Acoustic::evaluate(const step::Step &step, step::Frequency &frequency, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *ref, Vector_Base<double> *imf, Vector_Base<double> *renf, Vector_Base<double> *imnf, Vector_Base<double> *dirichlet)
{

}

void Acoustic::updateSolution(Vector_Base<double> *rex, Vector_Base<double> *imx)
{

}

//void Acoustic::connect(Harmonic &scheme)
//{
////    addFiller(*this, scheme);
//}
//
//void Acoustic::evaluate(Harmonic &scheme)
//{
////    controller.setUpdate();
//    reset(scheme.K, scheme.M, scheme.C, scheme.re.f, scheme.im.f, scheme.re.dirichlet, scheme.im.dirichlet);
////    iterate();
////    fill();
//    update(scheme.K, scheme.M, scheme.C, scheme.re.f, scheme.im.f, scheme.re.dirichlet, scheme.im.dirichlet);
////    controller.resetUpdate();
//}
//
//void Acoustic::updateSolution(Harmonic &scheme)
//{
//    scheme.re.x->storeTo(Results::pressure->data);
//}
