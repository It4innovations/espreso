
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

void StructuralMechanics::prepare()
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


	size_t clusters = *std::max_element(_mesh->getContinuityPartition().begin(), _mesh->getContinuityPartition().end()) + 1;

	_cCenter = _cNorm = std::vector<Point>(clusters, Point(0, 0, 0));
	_cr44 = _cr45 = _cr46 = _cr55 = _cr56 = std::vector<double>(clusters, 0);
	_cNp = std::vector<size_t>(clusters, 0);

	_dCenter = _dNorm = std::vector<Point>(_mesh->parts(), Point(0, 0, 0));
	_dr44 = _dr45 = _dr46 = _dr55 = _dr56 = std::vector<double>(_mesh->parts(), 0);
	_dNp = std::vector<size_t>(_mesh->parts(), 0);

	std::vector<double> cbuffer1(_mesh->parts(), 0), cbuffer2(_mesh->parts(), 0), cbuffer3(_mesh->parts(), 0);

	// Get center
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		Point center;
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			center += _mesh->coordinates().get(n, p);
		}
		_dCenter[p] = center;
	}

	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cCenter[_mesh->getContinuityPartition()[p]] += _dCenter[p];
		_dNp[p] = _mesh->coordinates().localSize(p);
		_dCenter[p] = _dCenter[p] / _dNp[p];
		_cNp[_mesh->getContinuityPartition()[p]] += _dNp[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cCenter[c] /= _cNp[c];
	}

	// Compute norm of column 4 (norm.x)
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		double pnorm = 0, pcnorm = 0;
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			Point dp = _mesh->coordinates().get(n, p) - _dCenter[p];
			pnorm += dp.x * dp.x + dp.y * dp.y;
			Point cp = _mesh->coordinates().get(n, p) - _cCenter[_mesh->getContinuityPartition()[p]];
			pcnorm += cp.x * cp.x + cp.y * cp.y;
		}
		_dNorm[p].x = std::sqrt(pnorm);
		cbuffer1[p] += pcnorm;
	}
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cNorm[_mesh->getContinuityPartition()[p]].x += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].x = std::sqrt(_cNorm[c].x);
	}

	// Compute coefficient r44, r45
	cbuffer1 = cbuffer2 = std::vector<double>(_mesh->parts(), 0);
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		size_t c = _mesh->getContinuityPartition()[p];
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			Point dp = _mesh->coordinates().get(n, p) - _dCenter[p];
			_dr44[p] += (-dp.y / _dNorm[p].x) * (-dp.y / _dNorm[p].x) + (dp.x / _dNorm[p].x) * (dp.x / _dNorm[p].x);
			_dr45[p] += (-dp.y / _dNorm[p].x) * (-dp.z);

			Point cp = _mesh->coordinates().get(n, p) - _cCenter[c];
			cbuffer1[p] += (-cp.y / _cNorm[c].x) * (-cp.y / _cNorm[c].x) + (cp.x / _cNorm[c].x) * (cp.x / _cNorm[c].x);
			cbuffer2[p] += (-cp.y / _cNorm[c].x) * (-cp.z);
		}
	}
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cr44[_mesh->getContinuityPartition()[p]] += cbuffer1[p];
		_cr45[_mesh->getContinuityPartition()[p]] += cbuffer2[p];
	}

	// Compute norm of column 5 (norm.y)
	cbuffer1 = std::vector<double>(_mesh->parts(), 0);
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		double dnorm = 0, cnorm = 0;
		size_t c = _mesh->getContinuityPartition()[p];
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			Point dp = _mesh->coordinates().get(n, p) - _dCenter[p];
			dnorm += (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) * (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x));
			dnorm += (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) * (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x));
			dnorm += dp.x * dp.x;

			Point cp = _mesh->coordinates().get(n, p) - _cCenter[c];
			cnorm += (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) * (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x));
			cnorm += (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) * (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x));
			cnorm += cp.x * cp.x;
		}
		_dNorm[p].y = std::sqrt(dnorm);
		cbuffer1[p] = cnorm;
	}
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cNorm[_mesh->getContinuityPartition()[p]].y += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].y = std::sqrt(_cNorm[c].y);
	}

	// Compute coefficient r46, r55, r56
	cbuffer1 = cbuffer2 = cbuffer3 = std::vector<double>(_mesh->parts(), 0);
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		double c5;
		size_t c = _mesh->getContinuityPartition()[p];
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			Point dp = _mesh->coordinates().get(n, p) - _dCenter[p];
			_dr46[p] += (dp.x / _dNorm[p].x) * (-dp.z);
			c5 = (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) / _dNorm[p].y;
			_dr55[p] += c5 * c5;
			_dr56[p] += c5 * 0;
			c5 = (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) / _dNorm[p].y;
			_dr55[p] += c5 * c5;
			_dr56[p] += c5 * (-dp.z);
			c5 = ( dp.x -                                           0) / _dNorm[p].y;
			_dr55[p] += c5 * c5;
			_dr56[p] += c5 * dp.y;

			Point cp = _mesh->coordinates().get(n, p) - _cCenter[c];
			cbuffer1[p] += (cp.x / _cNorm[c].x) * (-cp.z);
			c5 = (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) / _cNorm[c].y;
			cbuffer2[p] += c5 * c5;
			cbuffer3[p] += c5 * 0;
			c5 = (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) / _cNorm[c].y;
			cbuffer2[p] += c5 * c5;
			cbuffer3[p] += c5 * (-cp.z);
			c5 = ( cp.x -                                           0) / _cNorm[c].y;
			cbuffer2[p] += c5 * c5;
			cbuffer3[p] += c5 * cp.y;
		}
	}
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cr46[_mesh->getContinuityPartition()[p]] += cbuffer1[p];
		_cr55[_mesh->getContinuityPartition()[p]] += cbuffer2[p];
		_cr56[_mesh->getContinuityPartition()[p]] += cbuffer3[p];
	}

	// Compute norm of column 6 (norm.z)
	cbuffer1 = std::vector<double>(_mesh->parts(), 0);
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		double dnorm = 0, cnorm = 0, c6;
		size_t c = _mesh->getContinuityPartition()[p];
		for (size_t n = 0; n < _mesh->coordinates().localSize(p); n++) {
			Point dp = _mesh->coordinates().get(n, p) - _dCenter[p];
			c6 =     0 - _dr56[p] / _dr55[p] * (-dp.z - _dr45[p] / _dr44[p] * (-dp.y / _dNorm[p].x)) / _dNorm[p].y - _dr46[p] / _dr44[p] * (-dp.y / _dNorm[p].x);
			dnorm += c6 * c6;
			c6 = -dp.z - _dr56[p] / _dr55[p] * (    0 - _dr45[p] / _dr44[p] * ( dp.x / _dNorm[p].x)) / _dNorm[p].y - _dr46[p] / _dr44[p] * ( dp.x / _dNorm[p].x);
			dnorm += c6 * c6;
			c6 =  dp.y - _dr56[p] / _dr55[p] * ( dp.x -                                           0) / _dNorm[p].y - _dr46[p] / _dr44[p] * (    0 / _dNorm[p].x);
			dnorm += c6 * c6;

			Point cp = _mesh->coordinates().get(n, p) - _cCenter[c];
			c6 =     0 - _cr56[c] / _cr55[c] * (-cp.z - _cr45[c] / _cr44[c] * (-cp.y / _cNorm[c].x)) / _cNorm[c].y - _cr46[c] / _cr44[c] * (-cp.y / _cNorm[c].x);
			cnorm += c6 * c6;
			c6 = -cp.z - _cr56[c] / _cr55[c] * (    0 - _cr45[c] / _cr44[c] * ( cp.x / _cNorm[c].x)) / _cNorm[c].y - _cr46[c] / _cr44[c] * ( cp.x / _cNorm[c].x);
			cnorm += c6 * c6;
			c6 =  cp.y - _cr56[c] / _cr55[c] * ( cp.x -                                           0) / _cNorm[c].y - _cr46[c] / _cr44[c] * (    0 / _cNorm[c].x);
			cnorm += c6 * c6;
		}
		_dNorm[p].z = std::sqrt(dnorm);
		cbuffer1[p] = cnorm;
	}
	for (size_t p = 0; p < _mesh->parts(); p++) {
		_cNorm[_mesh->getContinuityPartition()[p]].z += cbuffer1[p];
	}
	for (size_t c = 0; c < clusters; c++) {
		_cNorm[c].z = std::sqrt(_cNorm[c].z);
	}
}

void StructuralMechanics::updateMesh(const std::vector<std::vector<eslocal> > &previousDOFMap, const std::vector<std::vector<eslocal> > &previousDomainMap)
{

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

std::vector<size_t> StructuralMechanics::solutionsIndicesToStore() const
{
	return { offset + SolutionIndex::DISPLACEMENT };
}










