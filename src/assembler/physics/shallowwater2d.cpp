
#include "../../configuration/physics/shallowwater2d.h"
#include "../step.h"
#include "../instance.h"
#include "../solution.h"

#include "../../mesh/settings/property.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/structures/elementtypes.h"


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
	_instance->DOFs = _mesh->assignUniformDOFsIndicesToNodes(_instance->DOFs, pointDOFs(), _nodesDOFsOffsets);
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
	DenseMatrix coordinates(e->nodes(), 2), J(2, 2), invJ(2, 2), Cu(e->nodes(), e->nodes());
	double detJ;
	DenseMatrix U(e->nodes(), 2), u(e->nodes(), 2), h(e->nodes(), 1), p(e->nodes(), 1);
	DenseMatrix gpU(1, 2), gpu(1, 2), gph(1, 1), gpp(1, 1);

	for (size_t i = 0; i < e->nodes(); i++) {
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;

		if (solution.size()) {
			U(i, 0) = solution[0]->get(Property::MOMENTUM_X, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
			U(i, 1) = solution[0]->get(Property::MOMENTUM_Y, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
			p(i, 0) = solution[0]->get(Property::PRESSURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
		} else {
			U(i, 0) = 0; // INITIAL MOMENTUM X
			U(i, 1) = 0; // INITIAL MOMENTUM Y
			p(i, 0) = 1; // INITIAL PRESSURE
		}
	}

	eslocal Ksize = pointDOFs().size() * e->nodes();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & Matrices::M) {
		Me.resize(e->nodes(), e->nodes());
		Me = 0;
	}
	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		J.multiply(e->dN()[gp], coordinates);
		detJ = determinant2x2(J.values());
		inverse2x2(J.values(), invJ.values(), detJ);

		gpU.multiply(e->N()[gp], U);
		gpp.multiply(e->N()[gp], p);

		Me.multiply(e->N()[gp], e->N()[gp], e->weighFactor()[gp] * detJ, 1, true, false);
		Cu.multiply(e->N()[gp], e->N()[gp], e->weighFactor()[gp] * detJ, 1, true, false);
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

	_instance->solutions[0] = new Solution("temperature", ElementType::NODES, pointDOFs(), _instance->primalSolution);
}

