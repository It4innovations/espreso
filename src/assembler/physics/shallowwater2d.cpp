
#include "../../configuration/physics/shallowwater2d.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../../output/resultstore.h"

#include "../../mesh/settings/property.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"


#include "../../basis/matrices/denseMatrix.h"
#include "../../solver/generic/SparseMatrix.h"
#include "shallowwater2d.h"

using namespace espreso;

ShallowWater2D::ShallowWater2D(Mesh *mesh, Instance *instance, const ShallowWater2DConfiguration &configuration)
: Physics2D("ADVECTION DIFFUSION 2D", mesh, instance), _configuration(configuration)
{

}

MatrixType ShallowWater2D::getMatrixType(const Step &step, size_t domain) const
{
	return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
}

void ShallowWater2D::prepareTotalFETI()
{
	_instance->DOFs = _mesh->assignUniformDOFsIndicesToNodes(_instance->DOFs, pointDOFs());
	_mesh->computeNodesDOFsCounters(pointDOFs());

	_mesh->loadNodeProperty(_configuration.momentum, { "X", "Y" }, { Property::MOMENTUM_X, Property::MOMENTUM_Y });
	_mesh->loadMaterials(_configuration.materials, _configuration.material_set);
	_mesh->removeDuplicateRegions();
	_mesh->fillDomainsSettings();
}

void ShallowWater2D::analyticRegularization(size_t domain)
{
	ESINFO(GLOBAL_ERROR) << "Implement analytic regularization";
}


void ShallowWater2D::processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];
	DenseMatrix coordinates(e->nodes(), 2), J(2, 2), invJ(2, 2), dND;
	double detJ, temp;

	for (size_t i = 0; i < e->nodes(); i++) {
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
	}

	eslocal Ksize = 3 * e->nodes();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & Matrices::K) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}
	if (matrices & Matrices::M) {
		Me.resize(Ksize, Ksize);
		Me = 0;
	}
	if (matrices & Matrices::R) {
		Re.resize(Ksize, 1);
		Re = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	DenseMatrix u(1, 2), v(1, 2), re(1, e->nodes());
	double normGradN = 0;

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		u.multiply(e->N()[gp], U, 1, 0);

		J.multiply(e->dN()[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpThickness.multiply(e->N()[gp], thickness);
		gpK.multiply(e->N()[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, e->dN()[gp]);
	}
}

void ShallowWater2D::processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "Shallow water 2D cannot process face";
}

void ShallowWater2D::processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{

}

void ShallowWater2D::processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{

}

void ShallowWater2D::processSolution(const Step &step)
{
	_instance->solutions.resize(1, NULL);
	if (_instance->solutions[0] != NULL) {
		delete _instance->solutions[0];
	}

	_instance->solutions[0] = new Solution("temperature", store::ElementType::NODES, pointDOFs(), _instance->primalSolution);
}

