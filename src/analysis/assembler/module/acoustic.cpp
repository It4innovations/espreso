
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
		Variable::list.node["INITIAL_ACOUSTIC_PRESSURE"] = Variable(0, 1, ParametersAcousticPressure::Initial::output->data.data(), false, true);
	}
	if (ParametersAcousticPressure::output == nullptr) {
		ParametersAcousticPressure::output = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "ACOUSTIC_PRESSURE");
		Variable::list.node["ACOUSTIC_PRESSURE"] = Variable(0, 1, ParametersAcousticPressure::output->data.data(), false, true);
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
		setMaterials(info::ecf->acoustic_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->acoustic_2d.material_set);
		validateRegionSettings("THICKNESS", info::ecf->acoustic_2d.thickness);
	}
	if (info::mesh->dimension == 3) {
		setMaterials(info::ecf->acoustic_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->acoustic_3d.material_set);
	}

	initParameters();

	baseFunction(*this);
	elementCoordinates(*this);
	elementIntegration(*this);

	bool correct = true;

	if (step::step.loadstep == 0) {
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->acoustic_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->acoustic_3d.material_set);
		}

		eslog::info(" ============================================================================================= \n");
	}

	if (K != nullptr) {
		acousticStiffness(*this);
	}
	if (M != nullptr) {
		acousticMass(*this);
	}

	if (configuration.acoustic_pressure.size()) {
		correct &= examineBoundaryParameter("ACOUSTIC_PRESSURE", configuration.acoustic_pressure, dirichlet.gp.externalValues);
		fromExpression(*this, dirichlet.gp, dirichlet.gp.externalValues);
	}
	if (configuration.normal_acceleration.size()) {
		examineBoundaryParameter("NORMAL ACCELERATION", configuration.normal_acceleration, normalAcceleration.gp.externalValues);
		fromExpression(*this, normalAcceleration.gp, normalAcceleration.gp.externalValues);
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
	updateVersions();

	if (K != nullptr) {
		K->fill(0);
		K->touched = true;
	}
	if (M != nullptr) {
		M->fill(0);
		M->touched = true;
	}
	if (C != nullptr) {
//		C->fill(0);
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

	std::cout << "COO[nd]: " << *coords.node.data << "\n";
	std::cout << "stiffness: " << *elements.stiffness.data << "\n";
	std::cout << "mass: " << *elements.mass.data << "\n";
}

void AX_Acoustic::initDirichlet(Vector_Sparse<double> &dirichlet)
{
	size_t dsize = 0;
	for (auto it = configuration.acoustic_pressure.begin(); it != configuration.acoustic_pressure.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dsize += region->nodes->datatarray().size();
	}
	dirichletIndices.reserve(dsize);
	for (auto it = configuration.acoustic_pressure.begin(); it != configuration.acoustic_pressure.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		dirichletIndices.insert(dirichletIndices.end(), region->nodes->datatarray().begin(), region->nodes->datatarray().end());
	}
	dirichletPermutation.resize(dsize);
	std::iota(dirichletPermutation.begin(), dirichletPermutation.end(), 0);
	std::sort(dirichletPermutation.begin(), dirichletPermutation.end(), [&] (const esint &i, const esint &j) { return dirichletIndices[i] < dirichletIndices[j]; });
	dsize = 0;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			++dsize;
		}
	}
	dirichlet.resize(info::mesh->nodes->IDs->datatarray().size(), dsize);
	auto dir = dirichlet.indices;
	for (auto i = dirichletPermutation.begin(); i != dirichletPermutation.end(); ++i) {
		if (i == dirichletPermutation.begin() || dirichletIndices[*i] != dirichletIndices[*(i - 1)]) {
			*dir++ = dirichletIndices[*i];
		}
	}
}

void AX_Acoustic::fillDirichlet(Vector_Sparse<double> &dirichlet)
{
	size_t offset = 0;
	std::vector<double> values(dirichletPermutation.size());
	for (auto it = configuration.acoustic_pressure.begin(); it != configuration.acoustic_pressure.end(); ++it) {
		BoundaryRegionStore *region = info::mesh->bregion(it->first);
		it->second.evaluator->evalSelectedSparse(
				region->nodes->datatarray().size(),
				region->nodes->datatarray().data(),
				it->second.evaluator->params,
				values.data() + offset);
		offset += region->nodes->datatarray().size();
	}

	for (size_t i = 0; i < dirichletPermutation.size(); ++i) {
		dirichlet.vals[i] = values[dirichletPermutation[i]];
	}
	dirichlet.touched = true;
}

void AX_Acoustic::updateSolution()
{
	re.x->store(ParametersAcousticPressure::output->data);
}
