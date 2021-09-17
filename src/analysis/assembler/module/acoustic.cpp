
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

AX_Acoustic::AX_Acoustic(AX_Acoustic *previous, AcousticGlobalSettings &gsettings, AcousticLoadStepConfiguration &configuration)
: gsettings(gsettings), configuration(configuration), K{}, M{}, C{}, re{}, im{}
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
	eslog::info("\n ============================================================================================= \n");
	eslog::info("  PHYSICS                                                                       ACOUSTIC REAL  \n");
	eslog::info(" ============================================================================================= \n");

	bool correct = true;

	if (info::mesh->dimension == 2) {
		setMaterials(info::ecf->acoustics_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->acoustics_2d.material_set);
		validateRegionSettings("THICKNESS", info::ecf->acoustics_2d.thickness);
	}
	if (info::mesh->dimension == 3) {
		setMaterials(info::ecf->acoustics_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->acoustics_3d.material_set);
	}

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	if (configuration.acoustic_pressure.size()) {
		correct &= examineBoundaryParameter("FIXED ACOUSTIC PRESSURE ON BOUNDARIES", configuration.acoustic_pressure, pressure.node.externalValues);
		fromExpression(*this, pressure.node, pressure.node.externalValues);
	}

	if (step::step.loadstep == 0) {
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");

		eslog::info("\n  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");

		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			MaterialConfiguration *mat = info::mesh->materials[i];
			correct &= examineMaterialParameter(mat->name, "DENSITY", mat->density, material.density.externalValue, 0);
			correct &= examineMaterialParameter(mat->name, "SPEED_OF_SOUND", mat->speed_of_sound, material.speed_of_sound.externalValue, 0);
		}

		fromExpression(*this, material.density, material.density.externalValue);
		fromExpression(*this, material.speed_of_sound, material.speed_of_sound.externalValue);

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");

		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->acoustics_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->acoustics_3d.material_set);
		}

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
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	} else {
		eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       FAIL  \n");
	}
	eslog::info(" ============================================================================================= \n");
	if (!correct) {
		eslog::globalerror("                                                               INVALID CONFIGURATION DETECTED \n");
	}
}

void AX_Acoustic::evaluate()
{
	controller.setUpdate();
//	printVersions();

	reset(K, M, C, re.rhs, im.rhs, re.dirichlet, im.dirichlet);

	iterate();
	fill();
	controller.resetUpdate();
}

void AX_Acoustic::updateSolution()
{
	re.x->store(ParametersAcousticPressure::output->data);
}
