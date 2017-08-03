
#include "lamesteklovpoincare3d.h"
#include "../../configuration/physics/structuralmechanics3d.h"

#include "../../basis/logging/logging.h"
#include "../../basis/matrices/sparseVVPMatrix.h"
#include "../../basis/matrices/denseMatrix.h"

#include "../instance.h"
#include "../step.h"
#include "../solution.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/elementtypes.h"

#include "../../solver/generic/SparseMatrix.h"

#ifdef BEM4I
#include "esbem.h"
#endif

#include "../../basis/utilities/utils.h"

using namespace espreso;

size_t LameSteklovPoincare3D::BEMOffset = -1;

LameSteklovPoincare3D::LameSteklovPoincare3D(Mesh *mesh, Instance *instance, const StructuralMechanics3DConfiguration &configuration)
: Physics("LAME STEKLOV POINCARE 3D", mesh, instance), StructuralMechanics3D(mesh, instance, configuration)
{
#ifndef BEM4I
	ESINFO(GLOBAL_ERROR) << "BEM4I is not linked!. Copy BEM4I library to tools/bem4i and re-configure ESPRESO.";
#endif
}

void LameSteklovPoincare3D::prepare()
{
	StructuralMechanics3D::prepare();
	extractBoundaryNodes();
}

void LameSteklovPoincare3D::prepareHybridTotalFETIWithKernels()
{
	// extraction of boundary nodes compute all faces on domains. Hence, no face computation is needed
	prepare();
}

void LameSteklovPoincare3D::preprocessData(const Step &step)
{
	if (offset == -1) {
		offset = _instance->solutions.size();
		_instance->solutions.resize(offset + SolutionIndex::SIZE, NULL);
		_instance->solutions[offset + SolutionIndex::DISPLACEMENT] = new Solution(*_mesh, "displacement", ElementType::NODES, pointDOFs());
	}

	if (BEMOffset == -1) {
		BEMOffset = _instance->solutions.size();
		_instance->solutions.resize(BEMOffset + BEMSolutionIndex::SIZE, NULL);
		_instance->solutions[BEMOffset + BEMSolutionIndex::BOUNDARY] = new Solution(*_mesh, "temperatureOnBoundary", ElementType::NODES, pointDOFs(), _instance->primalSolution);
	}
}

void LameSteklovPoincare3D::updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution)
{
	if (matrices & Matrices::f) {
		_instance->f[domain].clear();
		_instance->f[domain].resize(_instance->domainDOFCount[domain]);
	}

	std::vector<eslocal> elements;
	std::vector<double> coordinates;
	boundaryTriangularization(elements, coordinates, domain);

	_instance->K[domain].rows = pointDOFs().size() * _boundaryIndices[domain].size();
	_instance->K[domain].cols = pointDOFs().size() * _boundaryIndices[domain].size();
	_instance->K[domain].nnz  = _instance->K[domain].rows * _instance->K[domain].cols;
	_instance->K[domain].type = 'G';
	_instance->K[domain].dense_values.resize(_instance->K[domain].nnz);
	_instance->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

	DenseMatrix permuted(_instance->K[domain].rows, _instance->K[domain].cols);

#ifdef BEM4I
	bem4i::getLameSteklovPoincare(
			permuted.values(),
			(eslocal)_boundaryIndices[domain].size(),
			coordinates.data(),
			(eslocal)(elements.size() / 3),
			elements.data(),
			0.3, 2.1e5, // TODO: mi, E -> load from material
			3, 4, false);
#endif

	eslocal n = _instance->K[domain].rows;
	for (eslocal i = 0; i < n / 3; i++) {
		for (eslocal j = 0; j < n / 3; j++) {
			_instance->K[domain].dense_values[3 * i + 0 + (3 * i + 0) * n] = permuted(0 * (n / 3) + i, 0 * (n / 3) + j);
			_instance->K[domain].dense_values[3 * i + 1 + (3 * i + 1) * n] = permuted(1 * (n / 3) + i, 1 * (n / 3) + j);
			_instance->K[domain].dense_values[3 * i + 2 + (3 * i + 2) * n] = permuted(2 * (n / 3) + i, 2 * (n / 3) + j);
		}
	}

	_instance->K[domain].ConvertDenseToCSR(1);
	SparseVVPMatrix<eslocal> K, M; // never filled
	assembleBoundaryConditions(K, M, step, Matrices::f, domain, solution);
}

void LameSteklovPoincare3D::updateMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution)
{
	ESINFO(ERROR) << "BEM discretization not supports assembling of stiffness matrix K for one element.";
}

void LameSteklovPoincare3D::processSolution(const Step &step)
{
	// TODO: get solution for all nodes from BEM library
	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		std::fill(_instance->solutions[offset + SolutionIndex::DISPLACEMENT]->data[p].begin(), _instance->solutions[offset + SolutionIndex::DISPLACEMENT]->data[p].end(), 0);
		for (size_t i = 0; i < _boundaryIndices[p].size(); i++) {
			for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
				_instance->solutions[offset + SolutionIndex::DISPLACEMENT]->data[p][pointDOFs().size() * _mesh->coordinates().localIndex(_boundaryIndices[p][i], p) + dof] = _instance->primalSolution[p][pointDOFs().size() * i + dof];
			}
		}
	}
}






