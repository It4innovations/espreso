
#include "assembler.h"

#include "../../../../basis/matrices/denseMatrix.h"
#include "../../../../basis/matrices/sparseVVPMatrix.h"
#include "../../../../basis/matrices/sparseCSRMatrix.h"
#include "../../../../basis/utilities/utils.h"
#include "../../../../configuration/physics/advectiondiffusion3d.h"
#include "../../../../solver/generic/SparseMatrix.h"

#include "../../../../mesh/elements/element.h"
#include "../../../../mesh/settings/evaluator.h"
#include "../../../../mesh/structures/mesh.h"
#include "../../../../mesh/structures/region.h"
#include "../../../../mesh/structures/material.h"
#include "../../../../mesh/structures/elementtypes.h"

#include "../../../../mesh/elements/plane/square4.h"
#include "../../../../mesh/elements/plane/square8.h"
#include "../../../../mesh/elements/plane/triangle3.h"
#include "../../../../mesh/elements/plane/triangle6.h"

#include "../../../../mesh/elements/volume/hexahedron20.h"
#include "../../../../mesh/elements/volume/hexahedron8.h"
#include "../../../../mesh/elements/volume/prisma15.h"
#include "../../../../mesh/elements/volume/prisma6.h"
#include "../../../../mesh/elements/volume/pyramid13.h"
#include "../../../../mesh/elements/volume/pyramid5.h"
#include "../../../../mesh/elements/volume/tetrahedron10.h"
#include "../../../../mesh/elements/volume/tetrahedron4.h"
#include "../../../../output/resultstore.h"

#include "../../../constraints/equalityconstraints.h"

#ifdef ESBEM
#include "../../../../../tools/bem4i/esbem.h"
#endif


namespace espreso {

std::vector<Property> LaplaceSteklovPoincare::elementDOFs;
std::vector<Property> LaplaceSteklovPoincare::faceDOFs;
std::vector<Property> LaplaceSteklovPoincare::edgeDOFs;
std::vector<Property> LaplaceSteklovPoincare::pointDOFs = { Property::TEMPERATURE };
std::vector<Property> LaplaceSteklovPoincare::midPointDOFs = { Property::TEMPERATURE };

LaplaceSteklovPoincare::LaplaceSteklovPoincare(Mesh &mesh, Constraints &constraints, const AdvectionDiffusion3DConfiguration &configuration)
: LinearPhysics(
		mesh, constraints, configuration.espreso,
		MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE,
		elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs),
  _configuration(configuration)
{
	if (_configuration.translation_motions.size()) {
		mtype = MatrixType::REAL_UNSYMMETRIC;
	}
	_boundaryIndices.resize(mesh.parts());
};

void LaplaceSteklovPoincare::prepareMeshStructures()
{
	Hexahedron8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Hexahedron20::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron10::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma15::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid5::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid13::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	std::vector<size_t> DOFsOffsets;
	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs, DOFsOffsets);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	_mesh.computeFacesOnDomainsSurface();

	for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 0; i < _mesh.faces().size(); i++) {
			if (std::find(_mesh.faces()[i]->domains().begin(), _mesh.faces()[i]->domains().end(), p) != _mesh.faces()[i]->domains().end()) {
				for (size_t n = 0; n < _mesh.faces()[i]->nodes(); n++) {
					_boundaryIndices[p].push_back(_mesh.faces()[i]->node(n));
				}
			}
		}

		std::sort(_boundaryIndices[p].begin(), _boundaryIndices[p].end());
		Esutils::removeDuplicity(_boundaryIndices[p]);

		for (size_t i = 0; i < _boundaryIndices[p].size(); i++) {
			Element *node = _mesh.nodes()[_boundaryIndices[p][i]];
			size_t d = std::lower_bound(node->domains().begin(), node->domains().end(), p) - node->domains().begin();
			node->DOFsIndices()[d] = i;
		}
		matrixSize[p] = _boundaryIndices[p].size();
	}

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computeVolumeCorners(1, true, true, false);
			break;
		case B0_TYPE::KERNELS:
			// _mesh.computeFacesSharedByDomains();
			break;
		default:
			break;
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.loadProperty(_configuration.initial_temperature, { }         , { Property::INITIAL_TEMPERATURE });

	_mesh.loadNodeProperty(_configuration.temperature, { }, { Property::TEMPERATURE });

	_mesh.loadProperty(_configuration.translation_motions, { "X", "Y", "Z" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z });
	_mesh.loadProperty(_configuration.heat_source        , { }     , { Property::HEAT_SOURCE });
	_mesh.loadFaceProperty(_configuration.heat_flux      , { }     , { Property::HEAT_FLUX });
	_mesh.loadFaceProperty(_configuration.heat_flow      , { }     , { Property::HEAT_FLOW });

	for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
		std::map<std::string, std::string> values;
		for (auto regions = it->second.begin(); regions != it->second.end(); ++regions) {
			values[regions->first] = regions->second->external_temperature;
			_mesh.loadFaceProperty(values, { }, { Property::EXTERNAL_TEMPERATURE });
			values[regions->first] = regions->second->heat_transfer_coefficient;
			_mesh.loadFaceProperty(values, { }, { Property::HEAT_TRANSFER_COEFFICIENT });
		}
	}

	for (size_t r = 0; r < _mesh.regions().size(); r++) {
		for (size_t i = 0; i < _mesh.regions()[r]->settings.size(); i++) {
			if (_mesh.regions()[r]->settings[i].count(Property::HEAT_FLOW)) {
				_mesh.regions()[r]->computeArea(_mesh.coordinates());
				break;
			}
		}
	}

	_mesh.loadMaterials(_configuration.materials, _configuration.material_set);
	_mesh.removeDuplicateRegions();
	if (Test::report(EXPENSIVE)) {
		_mesh.checkRegions(_mesh.nodes());
	}
}

void LaplaceSteklovPoincare::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
}

void LaplaceSteklovPoincare::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		default:
			break;
		}
	}
}

void LaplaceSteklovPoincare::saveMeshProperties(output::ResultStore &store)
{
//	store.storeProperty("translationMotion", { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z }, ElementType::ELEMENTS);
//	store.storeProperty("headSource", { Property::HEAT_SOURCE }, ElementType::ELEMENTS);
//	store.storeProperty("temperature", { Property::TEMPERATURE }, ElementType::NODES);
}

void LaplaceSteklovPoincare::saveMeshResults(output::ResultStore &store, const std::vector<std::vector<double> > &results)
{
	std::vector<std::vector<double> > full(results.size());
	for (size_t p = 0; p < _mesh.parts(); p++) {
		full[p].resize(_mesh.coordinates().localSize(p));
		for (size_t i = 0; i < _boundaryIndices[p].size(); i++) {
			full[p][_mesh.coordinates().localIndex(_boundaryIndices[p][i], p)] = results[p][i];
		}
	}


	store.storeValues("temperature", 1, full, ElementType::NODES);
}

static double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
	);
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

static void processEdge(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* edge, const AdvectionDiffusion3DConfiguration &configuration, double area = 1)
{
	DenseMatrix coordinates(edge->nodes(), 3), dND(1, 3), q(edge->nodes(), 1), htc(edge->nodes(), 1), thickness(edge->nodes(), 1), flow(edge->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1);

	eslocal Ksize = edge->nodes();
	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize);
	Ke = 0;
	std::fill(fe.begin(), fe.end(), 0);

	const std::vector<DenseMatrix> &dN = edge->dN();
	const std::vector<DenseMatrix> &N = edge->N();
	const std::vector<double> &weighFactor = edge->weighFactor();

	for (size_t n = 0; n < edge->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[edge->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[edge->node(n)].y;
		coordinates(n, 2) = mesh.coordinates()[edge->node(n)].z;
		htc(n, 0) = edge->getProperty(Property::HEAT_TRANSFER_COEFFICIENT, n, 0, 0);
		q(n, 0) += htc(n, 0) * edge->getProperty(Property::EXTERNAL_TEMPERATURE, n, 0, 0);
		q(n, 0) += edge->getProperty(Property::HEAT_FLOW, n, 0, 0) / area;
		q(n, 0) += edge->getProperty(Property::HEAT_FLUX, n, 0, 0);
	}

	for (size_t gp = 0; gp < edge->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (edge->hasProperty(Property::EXTERNAL_TEMPERATURE, 0)) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % edge->nodes()) * gpQ(0, 0);
		}
	}
}

static void processFace(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* face, const AdvectionDiffusion3DConfiguration &configuration, double area = 1)
{
	DenseMatrix coordinates(face->nodes(), 3), dND(1, 3), q(face->nodes(), 1), htc(face->nodes(), 1), flow(face->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1);

	eslocal Ksize = face->nodes();
	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize);
	Ke = 0;
	std::fill(fe.begin(), fe.end(), 0);

	const std::vector<DenseMatrix> &dN = face->dN();
	const std::vector<DenseMatrix> &N = face->N();
	const std::vector<double> &weighFactor = face->weighFactor();

	for (size_t n = 0; n < face->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[face->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[face->node(n)].y;
		coordinates(n, 2) = mesh.coordinates()[face->node(n)].z;
		htc(n, 0) = face->getProperty(Property::HEAT_TRANSFER_COEFFICIENT, n, 0, 0);
		q(n, 0) += htc(n, 0) * face->getProperty(Property::EXTERNAL_TEMPERATURE, n, 0, 0);
		q(n, 0) += face->getProperty(Property::HEAT_FLOW, n, 0, 0) / area;
		q(n, 0) += face->getProperty(Property::HEAT_FLUX, n, 0, 0);
	}

	for (size_t gp = 0; gp < face->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();
		gpQ.multiply(N[gp], q);
		if (face->hasProperty(Property::EXTERNAL_TEMPERATURE, 0)) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % face->nodes()) * gpQ(0, 0);
		}
	}
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element, const AdvectionDiffusion3DConfiguration &configuration)
{
	size_t timeStep = 0;
	bool CAU = configuration.stabilization == AdvectionDiffusion3DConfiguration::STABILIZATION::CAU;
	double sigma = configuration.sigma;

	DenseMatrix Ce(3, 3), coordinates, J, invJ, dND;
	double detJ;
	double temp;
	DenseMatrix f(1, element->nodes());
	DenseMatrix U(element->nodes(), 3);
	DenseMatrix K(element->nodes(), 9), gpK(element->nodes(), 9);

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 3);
	for (size_t i = 0; i < element->nodes(); i++) {
		temp = element->getProperty(Property::INITIAL_TEMPERATURE, i, timeStep, 273.15 + 20);
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
		U(i, 0) =
				element->getProperty(Property::TRANSLATION_MOTION_X, i, timeStep, 0) *
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(element->node(i), timeStep, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(element->node(i), timeStep, temp);
		U(i, 1) =
				element->getProperty(Property::TRANSLATION_MOTION_Y, i, timeStep, 0) *
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(element->node(i), timeStep, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(element->node(i), timeStep, temp);
		U(i, 2) =
				element->getProperty(Property::TRANSLATION_MOTION_Z, i, timeStep, 0) *
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(element->node(i), timeStep, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(element->node(i), timeStep, temp);
		f(0, i) = element->sumProperty(Property::HEAT_SOURCE, i, timeStep, 0);

		switch (material->getModel(PHYSICS::ADVECTION_DIFFUSION_3D)) {
		case MATERIAL_MODEL::ISOTROPIC:
			K(i, 0) = K(i, 1) = K(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i), timeStep, temp);
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case MATERIAL_MODEL::DIAGONAL:
			K(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i), timeStep, temp);
			K(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i), timeStep, temp);
			K(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case MATERIAL_MODEL::SYMMETRIC:
			K(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i), timeStep, temp);
			K(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i), timeStep, temp);
			K(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 3) = K(i, 5) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i), timeStep, temp);
			K(i, 4) = K(i, 7) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 6) = K(i, 8) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i), timeStep, temp);
			break;
		case MATERIAL_MODEL::ANISOTROPIC:
			K(i, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i), timeStep, temp);
			K(i, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i), timeStep, temp);
			K(i, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 3) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i), timeStep, temp);
			K(i, 4) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 5) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX)->evaluate(element->node(i), timeStep, temp);
			K(i, 6) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i), timeStep, temp);
			K(i, 7) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZX)->evaluate(element->node(i), timeStep, temp);
			K(i, 8) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZY)->evaluate(element->node(i), timeStep, temp);
			break;
		default:
			ESINFO(ERROR) << "Advection diffusion 3D not supports set material model";
		}
	}

	eslocal Ksize = element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(1, 3), v(1, 3), Re(1, element->nodes());
	double normGradN = 0;

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);

		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);

		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->nodes()), b_e_c(1, element->nodes());
		b_e.multiply(u, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < Re.columns(); i++) {
					Re(0, i) = b_e(0, i) - f(0, i);
				}
				DenseMatrix ReBt(1, 3);
				ReBt.multiply(Re, dND, 1 / pow(normGradN, 2), 0, false, true);
				for (size_t i = 0; i < ReBt.columns(); i++) {
					v(0, i) = u(0, i) - ReBt(0, i);
				}
			} else {
				v = u;
			}
		}


		double norm_u_e = u.norm();
		double h_e = 0, tau_e = 0, konst = 0;
		double C_e = 0;

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);

			if (CAU) {
				DenseMatrix u_v(1, 3);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
				u_v(0, 2) = u(0, 2) - v(0, 2);
				b_e_c.multiply(u_v, dND, 1, 0);
				double norm_u_v = u_v.norm();
				double h_e_c = 2 * norm_u_v / b_e_c.norm();
				double P_e_c = h_e_c * norm_u_v / (2 * Ce.norm());
				double tau_e_c = std::max(0.0, 1 - 1 / P_e_c);

				double konst1 = Re.norm() / normGradN;
				double konst2 = tau_e * h_e != 0 ? tau_e_c * h_e_c / (tau_e * h_e) : 0;
				if (konst1 / norm_u_e < konst2) {
					C_e = tau_e * h_e * konst1 * (konst2 - konst1 / norm_u_e) / 2;
				} else {
					C_e = 0;
				}
			}
		}

		Ce(0, 0) += sigma * h_e * norm_u_e;
		Ce(1, 1) += sigma * h_e * norm_u_e;
		Ce(2, 2) += sigma * h_e * norm_u_e;

		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
		Ke.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
		if (konst * weighFactor[gp] * detJ != 0) {
			Ke.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
		}
		if (CAU) {
			Ke.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
		}

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i) * f(0, i);
			if (norm_u_e != 0) {
				fe[i] += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(0, i) / (2 * norm_u_e);
			}
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain, size_t matrix_size)
{
	R1.rows = matrix_size;
	R1.cols = 1;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.resize(R1.nnz, 1 / sqrt(matrix_size));
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat)
{
	RegMat.rows = K.rows;
	RegMat.cols = K.cols;
	RegMat.nnz  = 1;
	RegMat.type = K.type;

	RegMat.I_row_indices.push_back(1);
	RegMat.J_col_indices.push_back(1);
	RegMat.V_values.push_back(K.getDiagonalMaximum());
	RegMat.ConvertToCSR(1);
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R1, SparseMatrix &R2, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernels_from_nonsym_K(K, RegMat, R1, R2, norm, defect, subdomain);
}


static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);

}

void LaplaceSteklovPoincare::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e, _configuration);
	dofs.resize(e->nodes());
	for (size_t n = 0; n < e->nodes(); n++) {
		dofs[n] = e->node(n);
	}
}

void LaplaceSteklovPoincare::makeStiffnessMatricesRegular()
{
	#pragma omp parallel for
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		switch (_solverConfiguration.regularization) {
		case REGULARIZATION::FIX_POINTS:
			switch (mtype) {
			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
				if (singularK[subdomain]) {
					analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain, matrixSize[subdomain]);
					analyticsRegMat(K[subdomain], RegMat[subdomain]);
					K[subdomain].RemoveLower();
					RegMat[subdomain].RemoveLower();
					K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
					RegMat[subdomain].ConvertToCOO(1);
				} else {
					R1[subdomain].rows = 0;
					R1[subdomain].cols = 0;
					R1[subdomain].nnz  = 0;
					R1[subdomain].type = 'G';
					RegMat[subdomain].rows = 0;
					RegMat[subdomain].cols = 0;
					RegMat[subdomain].nnz  = 0;
					RegMat[subdomain].type = 'G';
				}
				break;
			case MatrixType::REAL_UNSYMMETRIC:
				ESINFO(ERROR) << "Cannot regularize stiffness matrix from fix point. Set REGULARIZATION = NULL_PIVOTS";
				break;
			default:
				ESINFO(ERROR) << "Unknown matrix type for regularization.";
			}
			break;
		case REGULARIZATION::NULL_PIVOTS:
			switch (mtype) {
			case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
				K[subdomain].RemoveLower();
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], RegMat[subdomain], subdomain);
				break;
			case MatrixType::REAL_UNSYMMETRIC:
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
				break;
			default:
				ESINFO(ERROR) << "Unknown matrix type for regularization.";
			}
			break;
		}
	}
}

void LaplaceSteklovPoincare::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

#ifndef ESBEM

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* e = _mesh.elements()[i];
		processElement(Ke, fe, _mesh, e, _configuration);

		for (size_t i = 0; i < e->nodes(); i++) {
			eslocal row = _mesh.nodes()[e->node(i)]->DOFIndex(subdomain, 0);
			for (size_t j = 0; j < e->nodes(); j++) {
				eslocal column = _mesh.nodes()[e->node(j)]->DOFIndex(subdomain, 0);
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

#else

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &nodes = _mesh.nodes();
//
//	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {
//
//		const Element* e = _mesh.elements()[i];
//		processElement(Ke, fe, _mesh, e, _configuration);
//
//		for (size_t i = 0; i < e->nodes(); i++) {
//			eslocal row = _mesh.nodes()[e->node(i)]->DOFIndex(subdomain, 0);
//			for (size_t j = 0; j < e->nodes(); j++) {
//				eslocal column = _mesh.nodes()[e->node(j)]->DOFIndex(subdomain, 0);
//				_K(row, column) = Ke(i, j);
//			}
//			f[subdomain][row] += fe[i];
//		}
//	}

	std::vector<eslocal> elements;
	std::vector<double> coordinates;

	for (size_t i = 0; i < _mesh.faces().size(); i++) {
		if (std::find(_mesh.faces()[i]->domains().begin(), _mesh.faces()[i]->domains().end(), subdomain) != _mesh.faces()[i]->domains().end()) {
			for (size_t n = 0; n < _mesh.faces()[i]->nodes(); n++) {
				elements.push_back(std::lower_bound(_boundaryIndices[subdomain].begin(), _boundaryIndices[subdomain].end(), _mesh.faces()[i]->node(n)) - _boundaryIndices[subdomain].begin());
			}
		}
	}

	for (size_t i = 0; i < _boundaryIndices[subdomain].size(); i++) {
		coordinates.push_back(_mesh.coordinates()[_boundaryIndices[subdomain][i]].x);
		coordinates.push_back(_mesh.coordinates()[_boundaryIndices[subdomain][i]].y);
		coordinates.push_back(_mesh.coordinates()[_boundaryIndices[subdomain][i]].z);
	}

	K[subdomain].rows = _boundaryIndices[subdomain].size();
	K[subdomain].cols = _boundaryIndices[subdomain].size();
	K[subdomain].nnz  = _boundaryIndices[subdomain].size() * _boundaryIndices[subdomain].size();
	K[subdomain].type = 'G';
	K[subdomain].dense_values.resize(K[subdomain].nnz);
	bem4i::getLaplaceSteklovPoincare(
			K[subdomain].dense_values.data(),
			(eslocal)_boundaryIndices[subdomain].size(),
			coordinates.data(),
			(eslocal)(elements.size() / 3),
			elements.data(),
			3, 3, 0);

	K[subdomain].ConvertDenseToCSR(1);
	K[subdomain].MatTranspose();

#endif

	auto processRegion = [&] (const std::vector<Element*> &elments, bool withK = false, double area = 1) {
		for (size_t i = 0; i < elments.size(); i++) {
			if (elments[i]->inDomain(subdomain)) {
				if (elments[0]->type() == Element::Type::PLANE) {
					processFace(Ke, fe, _mesh, elments[i], _configuration, area);
				}
				if (elments[0]->type() == Element::Type::LINE) {
					processEdge(Ke, fe, _mesh, elments[i], _configuration, area);
				}

				for (size_t nx = 0; nx < elments[i]->nodes(); nx++) {
					for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
						size_t row = nodes[elments[i]->node(nx)]->DOFIndex(subdomain, dx);
						f[subdomain][row] += fe[dx * elments[i]->nodes() + nx];
					}
				}
				if (withK) {
					singularK[subdomain] = false;
					for (size_t nx = 0; nx < elments[i]->nodes(); nx++) {
						for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
							size_t row = nodes[elments[i]->node(nx)]->DOFIndex(subdomain, dx);
							for (size_t ny = 0; ny < elments[i]->nodes(); ny++) {
								for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
									size_t column = nodes[elments[i]->node(ny)]->DOFIndex(subdomain, dy);
									_K(row, column) = Ke(dx * elments[i]->nodes() + nx, dy * elments[i]->nodes() + ny);
								}
							}
						}
					}
				}

			}
		}
	};


	size_t step = 0;
	for (size_t r = 0; r < _mesh.regions().size(); r++) {
		if (step < _mesh.regions()[r]->settings.size() && _mesh.regions()[r]->settings[step].count(Property::HEAT_FLUX)) {
			processRegion(_mesh.regions()[r]->elements());
		}
		if (step < _mesh.regions()[r]->settings.size() && _mesh.regions()[r]->settings[step].count(Property::HEAT_FLOW)) {
			processRegion(_mesh.regions()[r]->elements(), false, _mesh.regions()[r]->area);
		}
		if (step < _mesh.regions()[r]->settings.size() && _mesh.regions()[r]->settings[step].count(Property::EXTERNAL_TEMPERATURE)) {
			processRegion(_mesh.regions()[r]->elements(), true);
		}
	}

#ifndef ESBEM
	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
#endif
}

void LaplaceSteklovPoincare::postProcess(output::ResultStore &store, const std::vector<std::vector<double> > &solution)
{

}


}

