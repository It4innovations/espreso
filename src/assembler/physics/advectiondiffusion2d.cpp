
#include "../../configuration/physics/advectiondiffusion2d.h"
#include "advectiondiffusion2d.h"

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

using namespace espreso;

NewAdvectionDiffusion2D::NewAdvectionDiffusion2D(Mesh *mesh, Instance *instance, const AdvectionDiffusion2DConfiguration &configuration)
: Physics2D(mesh, instance), _configuration(configuration)
{

}

MatrixType NewAdvectionDiffusion2D::getMatrixType(const Step &step, size_t domain) const
{
	if (_mesh->hasProperty(domain, Property::TRANSLATION_MOTION_X, step.load) || _mesh->hasProperty(domain, Property::TRANSLATION_MOTION_Y, step.load)) {
		return MatrixType::REAL_UNSYMMETRIC;
	} else {
		return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	}
}

void NewAdvectionDiffusion2D::prepareTotalFETI()
{
	_instance->DOFs = _mesh->assignUniformDOFsIndicesToNodes(_instance->DOFs, pointDOFs());
	_mesh->computeNodesDOFsCounters(pointDOFs());

	_mesh->loadProperty(_configuration.initial_temperature, { }, { Property::INITIAL_TEMPERATURE });

	_mesh->loadNodeProperty(_configuration.temperature    , { }, { Property::TEMPERATURE });
	_mesh->loadNodeProperty(_configuration.thickness      , { }, { Property::THICKNESS });

	_mesh->loadProperty(_configuration.translation_motions, { "X", "Y" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y });
	_mesh->loadProperty(_configuration.heat_source        , { }         , { Property::HEAT_SOURCE });
	_mesh->loadProperty(_configuration.heat_flux          , { }         , { Property::HEAT_FLUX });
	_mesh->loadProperty(_configuration.heat_flow          , { }         , { Property::HEAT_FLOW });

	for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
		std::map<std::string, std::string> values;
		for (auto regions = it->second.begin(); regions != it->second.end(); ++regions) {
			values[regions->first] = regions->second->external_temperature;
			_mesh->loadProperty(values, { }, { Property::EXTERNAL_TEMPERATURE });
			values[regions->first] = regions->second->heat_transfer_coefficient;
			_mesh->loadProperty(values, { }, { Property::HEAT_TRANSFER_COEFFICIENT });
		}
	}

	for (size_t r = 0; r < _mesh->regions().size(); r++) {
		for (size_t i = 0; i < _mesh->regions()[r]->settings.size(); i++) {
			if (_mesh->regions()[r]->settings[i].count(Property::HEAT_FLOW)) {
				_mesh->regions()[r]->computeArea(_mesh->coordinates());
				break;
			}
		}
	}

	_mesh->loadMaterials(_configuration.materials, _configuration.material_set);
	_mesh->removeDuplicateRegions();
	_mesh->fillDomainsSettings();
}

void NewAdvectionDiffusion2D::analyticRegularization(size_t domain)
{
	if (_instance->K[domain].mtype != MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE) {
		ESINFO(ERROR) << "Cannot compute analytic regularization of not REAL_SYMMETRIC_POSITIVE_DEFINITE matrix. Set REGULARIZATION = NULL_PIVOTS";
	}

	if (_mesh->hasProperty(domain, Property::EXTERNAL_TEMPERATURE, 0)) {
		return;
	}

	_instance->N1[domain].rows = _mesh->coordinates().localSize(domain);
	_instance->N1[domain].cols = 1;
	_instance->N1[domain].nnz = _instance->N1[domain].rows * _instance->N1[domain].cols;
	_instance->N1[domain].type = 'G';

	_instance->N1[domain].dense_values.resize(_instance->N1[domain].nnz, 1 / sqrt(_mesh->coordinates().localSize(domain)));

	_instance->RegMat[domain].rows = _instance->K[domain].rows;
	_instance->RegMat[domain].cols = _instance->K[domain].cols;
	_instance->RegMat[domain].nnz  = 1;
	_instance->RegMat[domain].type = _instance->K[domain].type;

	_instance->RegMat[domain].I_row_indices.push_back(1);
	_instance->RegMat[domain].J_col_indices.push_back(1);
	_instance->RegMat[domain].V_values.push_back(_instance->K[domain].getDiagonalMaximum());
	_instance->RegMat[domain].ConvertToCSR(1);
}

void NewAdvectionDiffusion2D::assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K) const
{
	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];
	switch (material->getModel(PHYSICS::ADVECTION_DIFFUSION_2D)) {
	case MATERIAL_MODEL::ISOTROPIC:
		K(node, 0) = K(node, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.load, temp);
		K(node, 2) = K(node, 3) = 0;
		break;
	case MATERIAL_MODEL::DIAGONAL:
		K(node, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.load, temp);
		K(node, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.load, temp);
		K(node, 2) = K(node, 3) = 0;
		break;
	case MATERIAL_MODEL::SYMMETRIC:
		K(node, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.load, temp);
		K(node, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.load, temp);
		K(node, 2) = K(node, 3) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.load, temp);
		break;
	case MATERIAL_MODEL::ANISOTROPIC:
		K(node, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.load, temp);
		K(node, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.load, temp);
		K(node, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.load, temp);
		K(node, 3) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX)->evaluate(e->node(node), step.load, temp);
		break;
	default:
		ESINFO(ERROR) << "Advection diffusion 2D not supports set material model";
	}
}

void NewAdvectionDiffusion2D::assembleResidualForces(const Step &step, const Element *e, DenseMatrix &Re) const
{
	DenseMatrix Ke, fe, T(e->nodes(), 1);
	processElement(step, e, Ke, fe);
	for (size_t i = 0; i < e->nodes(); i++) {
		T(i, 0) = _instance->solutions[0]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
	}

	Re.multiply(Ke, T);
}

void NewAdvectionDiffusion2D::processElement(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const
{
	bool CAU = _configuration.stabilization == AdvectionDiffusion2DConfiguration::STABILIZATION::CAU;

	DenseMatrix Ce(2, 2), coordinates, J(2, 2), invJ(2, 2), dND;
	double detJ, temp;
	DenseMatrix f(e->nodes(), 1);
	DenseMatrix U(e->nodes(), 2);
	DenseMatrix thickness(e->nodes(), 1), K(e->nodes(), 4);
	DenseMatrix gpThickness(1, 1), gpK(1, 4);

	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	coordinates.resize(e->nodes(), 2);

	for (size_t i = 0; i < e->nodes(); i++) {
		if (_instance->solutions.size()) {
			temp = _instance->solutions[0]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
		} else {
			temp = e->getProperty(Property::INITIAL_TEMPERATURE, i, step.load, 273.15 + 20);
		}
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
		thickness(i, 0) = e->getProperty(Property::THICKNESS, i, step.load, 1);
		U(i, 0) =
				e->getProperty(Property::TRANSLATION_MOTION_X, i, step.load, 0) *
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.load, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(e->node(i), step.load, temp) *
				thickness(i, 0);
		U(i, 1) =
				e->getProperty(Property::TRANSLATION_MOTION_Y, i, step.load, 0) *
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.load, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(e->node(i), step.load, temp) *
				thickness(i, 0);
		f(i, 0) = e->sumProperty(Property::HEAT_SOURCE, i, step.load, 0) * thickness(i, 0);
		assembleMaterialMatrix(step, e, i, temp, K);
	}

	eslocal Ksize = e->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize, 1);
	fe = 0;

	DenseMatrix u(1, 2), v(1, 2), Re(1, e->nodes());
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

		DenseMatrix b_e(1, e->nodes()), b_e_c(1, e->nodes());
		b_e.multiply(u, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < Re.columns(); i++) {
					Re(0, i) = b_e(0, i) - f(i, 0);
				}
				DenseMatrix ReBt(1, 2);
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
				DenseMatrix u_v(1, 2);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
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

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;

		Ke.multiply(dND, Ce * dND, detJ * e->weighFactor()[gp] * gpThickness(0, 0), 1, true);
		Ke.multiply(e->N()[gp], b_e, detJ * e->weighFactor()[gp], 1, true);
		if (konst * e->weighFactor()[gp] * detJ != 0) {
			Ke.multiply(b_e, b_e, konst * e->weighFactor()[gp] * detJ, 1, true);
		}
		if (CAU) {
			Ke.multiply(dND, dND, C_e * e->weighFactor()[gp] * detJ, 1, true);
		}

		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += detJ * e->weighFactor()[gp] * e->N()[gp](0, i) * f(i, 0);
			if (norm_u_e != 0) {
				fe(i, 0) += detJ * e->weighFactor()[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
			}
		}
	}
}

void NewAdvectionDiffusion2D::processFace(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const
{
	ESINFO(ERROR) << "Advection diffusion 2D cannot process face";
}

void NewAdvectionDiffusion2D::processEdge(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const
{
	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.load) ||
		e->hasProperty(Property::HEAT_FLOW, step.load) ||
		e->hasProperty(Property::HEAT_FLUX, step.load))) {

		Ke.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	DenseMatrix coordinates(e->nodes(), 2), dND(1, 2), q(e->nodes(), 1), htc(e->nodes(), 1), thickness(e->nodes(), 1), flow(e->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1);

	double area = 1;
	eslocal Ksize = e->nodes();
	Ke.resize(0, 0);
	fe.resize(Ksize, 1);
	fe = 0;

	for (size_t r = 0; r < e->regions().size(); r++) {
		if (step.load < e->regions()[r]->settings.size() && e->regions()[r]->settings[step.load].count(Property::HEAT_FLOW)) {
			Ke.resize(Ksize, Ksize);
			Ke = 0;
			area = e->regions()[r]->area;
			break;
		}
	}

	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	for (size_t n = 0; n < e->nodes(); n++) {
		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;

		double temp = 0;
		if (_instance->solutions.size()) {
			temp = _instance->solutions[0]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
		}
		htc(n, 0) = e->getProperty(Property::HEAT_TRANSFER_COEFFICIENT, n, step.load, 0);
		q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.load, 0) - temp);
		q(n, 0) += e->getProperty(Property::HEAT_FLOW, n, step.load, 0) / area;
		q(n, 0) += e->getProperty(Property::HEAT_FLUX, n, step.load, 0);

		thickness(n, 0) = e->getProperty(Property::THICKNESS, n, step.load, 1);
		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.load)) {
			gpHtc.multiply(N[gp], htc);
			gpThickness.multiply(N[gp], thickness);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0) * gpThickness(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
		}
	}
}

void NewAdvectionDiffusion2D::processNode(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const
{

}

void NewAdvectionDiffusion2D::storeSolution(const Step &step, std::vector<std::vector<double> > &solution, store::ResultStore *store)
{
	_instance->solutions.resize(1, NULL);
	if (_instance->solutions[0] != NULL) {
		delete _instance->solutions[0];
	}

	_instance->solutions[0] = new Solution(pointDOFs(), solution);

	store->storeValues("temperature", 1, _instance->solutions[0]->data, store::ResultStore::ElementType::NODES);
}

