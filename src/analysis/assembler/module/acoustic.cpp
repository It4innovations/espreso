
#include "acoustic.h"
#include "assembler.hpp"

#include "basis/expression/variable.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "analysis/assembler/operators/operators.h"
#include "analysis/scheme/harmonic.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <numeric>
#include <algorithm>

#include "basis/utilities/print.h"

using namespace espreso;

AX_Acoustic::AX_Acoustic(AX_Acoustic *previous, AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration)
: settings(settings), configuration(configuration)
{

}

void AX_Acoustic::initParameters()
{
	if (ParametersAcousticPressure::Initial::output == nullptr) {
		ParametersAcousticPressure::Initial::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_ACOUSTIC_PRESSURE");
		Variable::list.node["INITIAL_ACOUSTIC_PRESSURE"] = new OutputVariable(ParametersAcousticPressure::Initial::output, 0, 1);
	}
	if (ParametersAcousticPressure::output == nullptr) {
		ParametersAcousticPressure::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "ACOUSTIC_PRESSURE");
		Variable::list.node["ACOUSTIC_PRESSURE"] = new OutputVariable(ParametersAcousticPressure::output, 0, 1);
	}
}

void AX_Acoustic::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");
	bool correct = true;

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("THICKNESS", settings.thickness);

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	if (configuration.acoustic_pressure.size()) {
		correct &= examineBoundaryParameter("FIXED ACOUSTIC PRESSURE ON BOUNDARIES", configuration.acoustic_pressure, pressure.node.externalValues);
		fromExpression(*this, pressure.node, pressure.node.externalValues);
	}

	if (step::step.loadstep == 0) {
		eslog::info("\n  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");

		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			MaterialConfiguration *mat = info::mesh->materials[i];
			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValues, 0);
			correct &= examineMaterialParameter(mat->name, "SPEED_OF_SOUND", mat->speed_of_sound, material.speed_of_sound.externalValues, 0);
		}

		fromExpression(*this, material.density, material.density.externalValues);
		fromExpression(*this, material.speed_of_sound, material.speed_of_sound.externalValues);

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");

		printMaterials(settings.material_set);

		eslog::info(" ============================================================================================= \n");
	}


	acousticStiffness(*this);
	acousticMass(*this);
	// boundary conditions have to added according to boundary settings below
//	acousticBoundaryMass(*this);

	if (configuration.normal_acceleration.size()) {
		examineBoundaryParameter("NORMAL ACCELERATION", configuration.normal_acceleration, normalAcceleration.gp.externalValues);
		fromExpression(*this, normalAcceleration.gp, normalAcceleration.gp.externalValues);
	}
	if (configuration.acceleration.size()) {
		correct &= examineBoundaryParameter("ACCELERATION.X", configuration.acceleration, acceleration.gp.externalValues, 0);
		correct &= examineBoundaryParameter("ACCELERATION.Y", configuration.acceleration, acceleration.gp.externalValues, 1);
		if (info::mesh->dimension == 3) {
			correct &= examineBoundaryParameter("ACCELERATION.Z", configuration.acceleration,  acceleration.gp.externalValues, 2);
		}
		fromExpression(*this, acceleration.gp, acceleration.gp.externalValues);
	}
	if (configuration.impedance.size()) {
		examineBoundaryParameter("IMPEDANCE", configuration.impedance, impedance.gp.externalValues);
		fromExpression(*this, impedance.gp, impedance.gp.externalValues);
	}
	if (configuration.monopole_source.size()) {
		correct &= examineElementParameter("MONOPOLE DOMAIN SOURCE", configuration.monopole_source, monopoleSource.gp.externalValues);
		fromExpression(*this, monopoleSource.gp, monopoleSource.gp.externalValues);
	}
	if (configuration.dipole_source.size()) {
		correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.X", configuration.dipole_source, dipoleSource.gp.externalValues, 0);
		correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.Y", configuration.dipole_source, dipoleSource.gp.externalValues, 1);
		if (info::mesh->dimension == 3) {
			correct &= examineElementParameter("DIPOLE_DOMAIN_SOURCE.Z", configuration.dipole_source, dipoleSource.gp.externalValues, 2);
		}
		fromExpression(*this, dipoleSource.gp, dipoleSource.gp.externalValues);
	}
	
	integration.weight.name = "integration.weight";
	integration.N.name = "integration.N";
	integration.dN.name = "integration.dN";
	integration.dND.name = "integration.dND";
	integration.jacobiDeterminant.name = "integration.jacobiDeterminant";
	integration.jacobiInversion.name = "integration.jacobiInversion";
	elements.monopole.name = "elements.monopole";
	elements.dipole.name = "elements.dipole";
	material.density.name = "material.density";
	
	acousticRHS(*this);

	eslog::info(" ============================================================================================= \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void AX_Acoustic::connect(AX_Harmonic &scheme)
{
	addFiller(*this, scheme);
}

void AX_Acoustic::evaluate(AX_Harmonic &scheme)
{
	controller.setUpdate();
	reset(scheme.K, scheme.M, scheme.C, scheme.re.f, scheme.im.f, scheme.re.dirichlet, scheme.im.dirichlet);
	iterate();
	fill();
	update(scheme.K, scheme.M, scheme.C, scheme.re.f, scheme.im.f, scheme.re.dirichlet, scheme.im.dirichlet);
	controller.resetUpdate();
}

void AX_Acoustic::updateSolution(AX_Harmonic &scheme)
{
	scheme.re.x->store(ParametersAcousticPressure::output->data);
}
