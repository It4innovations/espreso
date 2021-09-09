
#include "acoustic.h"
#include "assembler.hpp"

#include "basis/expression/variable.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "analysis/assembler/operators/operators.h"
#include "analysis/scheme/harmonic.real.h"

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

void AX_Acoustic::init(AX_HarmonicReal &scheme)
{
	K = scheme.K;
	M = scheme.M;
	C = scheme.C;
	re.rhs = scheme.re.f;
	im.rhs = scheme.im.f;
	re.x = scheme.re.x;
	im.x = scheme.im.x;

	analyze();
}

void AX_Acoustic::analyze()
{
	eslog::info("\n ============================================================================================= \n");
	eslog::info("  PHYSICS                                                                       ACOUSTIC REAL  \n");
	eslog::info(" ============================================================================================= \n");

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

	bool correct = true;

	if (step::step.loadstep == 0) {
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

	if (configuration.acoustic_pressure.size()) {
		correct &= examineBoundaryParameter("ACOUSTIC_PRESSURE", configuration.acoustic_pressure, dirichlet.gp.externalValues);
		fromExpression(*this, dirichlet.gp, dirichlet.gp.externalValues);
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

void AX_Acoustic::next()
{
	controller.setUpdate();
//	printVersions();

	if (K != nullptr) {
		K->fill(0);
		K->touched = true;
	}
	if (M != nullptr) {
		M->fill(0);
		M->touched = true;
	}
	if (C != nullptr) {
		C->fill(0);
		C->touched = true;
	}
	if (re.rhs != nullptr) {
		re.rhs->fill(0);
		re.rhs->touched = true;
	}
	if (im.rhs != nullptr) {
		im.rhs->fill(0);
		im.rhs->touched = true;
	}

	iterate();
	controller.resetUpdate();
}

void AX_Acoustic::updateSolution()
{
	re.x->store(ParametersAcousticPressure::output->data);
}
