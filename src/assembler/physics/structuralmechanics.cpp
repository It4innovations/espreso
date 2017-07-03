
#include "../../configuration/physics/structuralmechanics.h"
#include "structuralmechanics.h"

#include "../../basis/matrices/denseMatrix.h"
#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/elements/element.h"
#include "../../mesh/settings/property.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/elementtypes.h"

#include "../instance.h"
#include "../solution.h"
#include "../step.h"

using namespace espreso;

size_t StructuralMechanics::offset = -1;

StructuralMechanics::StructuralMechanics(const StructuralMechanicsConfiguration &configuration)
: Physics("", NULL, NULL), // skipped because Physics is inherited virtually
  _configuration(configuration)
{

}

MatrixType StructuralMechanics::getMatrixType(const Step &step, size_t domain) const
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

bool StructuralMechanics::isMatrixTimeDependent(const Step &step) const
{
	return _mesh->isAnyPropertyTimeDependent({
		Property::DISPLACEMENT_X,
		Property::DISPLACEMENT_Y,
		Property::DISPLACEMENT_Z,
		Property::PRESSURE,
		Property::OBSTACLE,
		Property::NORMAL_DIRECTION,
		Property::TEMPERATURE
	}, step.step);
}

bool StructuralMechanics::isMatrixTemperatureDependent(const Step &step) const
{
	return _mesh->isAnyPropertyTemperatureDependent({
		Property::DISPLACEMENT_X,
		Property::DISPLACEMENT_Y,
		Property::DISPLACEMENT_Z,
		Property::PRESSURE,
		Property::OBSTACLE,
		Property::NORMAL_DIRECTION,
		Property::TEMPERATURE
	}, step.step);
}

void StructuralMechanics::prepareTotalFETI()
{
	for (size_t s = 1; s <= _configuration.physics_solver.load_steps; s++) {
		if (
				_configuration.displacement.find(s) == _configuration.displacement.end() &&
				!_mesh->hasProperty(Property::DISPLACEMENT_X, s - 1) &&
				!_mesh->hasProperty(Property::DISPLACEMENT_Y, s - 1) &&
				!_mesh->hasProperty(Property::DISPLACEMENT_Z, s - 1)) {

			ESINFO(GLOBAL_ERROR) << "Invalid boundary conditions for Structural mechanics - missing displacement for LOAD_STEP=" << s;
		}
	}

	_instance->domainDOFCount = _mesh->assignUniformDOFsIndicesToNodes(_instance->domainDOFCount, pointDOFs(), _nodesDOFsOffsets);
	_instance->properties = pointDOFs();
	_mesh->computeNodesDOFsCounters(pointDOFs());

	_mesh->loadNodeProperty(_configuration.temperature     , { }, { Property::TEMPERATURE });
	_mesh->loadNodeProperty(_configuration.obstacle        , { }, { Property::OBSTACLE });
	_mesh->loadNodeProperty(_configuration.normal_direction, { }, { Property::NORMAL_DIRECTION });
	_mesh->loadProperty(_configuration.normal_presure      , { }, { Property::PRESSURE });
	_mesh->loadProperty(_configuration.initial_temperature , { }, { Property::INITIAL_TEMPERATURE });

	_mesh->removeDuplicateRegions();
	_mesh->fillDomainsSettings();
}

void StructuralMechanics::preprocessData(const Step &step)
{
	if (offset != (size_t)-1) {
		return;
	}
	offset = _instance->solutions.size();
	_instance->solutions.resize(offset + SolutionIndex::SIZE, NULL);
	_instance->solutions[offset + SolutionIndex::DISPLACEMENT] = new Solution(*_mesh, "displacement", ElementType::NODES, pointDOFs(), _instance->primalSolution);
}

std::vector<size_t> StructuralMechanics::solutions() const
{
	return { offset + SolutionIndex::DISPLACEMENT };
}










