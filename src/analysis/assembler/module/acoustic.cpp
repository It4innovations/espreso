
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
: settings(settings), configuration(configuration), K{}, M{}, C{}, re{}, im{}
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

void AX_Acoustic::init(AX_Harmonic &scheme)
{
	this->K = scheme.K;
	this->M = scheme.M;
	this->C = scheme.C;
	this->re.rhs = scheme.re.f;
	this->im.rhs = scheme.im.f;
	this->re.x = scheme.re.x;
	this->im.x = scheme.im.x;
	this->re.dirichlet = scheme.re.dirichlet;
	this->im.dirichlet = scheme.im.dirichlet;

	analyze();
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


	if (K != nullptr) {
		acousticStiffness(*this);
	}
	if (M != nullptr) {
		acousticMass(*this);
	}
	if (C != nullptr) {
		// boundary conditions have to added according to boundary settings below
//		acousticBoundaryMass(*this);
	}

	if (configuration.normal_acceleration.size()) {
		examineBoundaryParameter("NORMAL ACCELERATION", configuration.normal_acceleration, normalAcceleration.gp.externalValues);
		fromExpression(*this, normalAcceleration.gp, normalAcceleration.gp.externalValues);
	}
	if (configuration.impedance.size()) {
		examineBoundaryParameter("IMPEDANCE", configuration.impedance, impedance.gp.externalValues);
		fromExpression(*this, impedance.gp, impedance.gp.externalValues);
	}
	acousticRHS(*this);

	addFiller(*this);

	eslog::info(" ============================================================================================= \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void AX_Acoustic::evaluate()
{
	controller.setUpdate();
	reset(K, M, C, re.rhs, im.rhs, re.dirichlet, im.dirichlet);
	iterate();
	fill();
	update(K, M, C, re.rhs, im.rhs, re.dirichlet, im.dirichlet);
	controller.resetUpdate();
}

void AX_Acoustic::updateSolution()
{
	re.x->store(ParametersAcousticPressure::output->data);
}
