
#include "../../configuration/physics/advectiondiffusion3d.h"
#include "advectiondiffusion3d.h"

#include "../step.h"
#include "../instance.h"
#include "../solution.h"
#include "../constraints/equalityconstraints.h"

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

using namespace espreso;

AdvectionDiffusion3D::AdvectionDiffusion3D(Mesh *mesh, Instance *instance, const AdvectionDiffusion3DConfiguration &configuration)
: Physics("ADVECTION DIFFUSION 3D", mesh, instance), AdvectionDiffusion(configuration), _configuration(configuration)
{
	_equalityConstraints = new EqualityConstraints(*_instance, *_mesh, _mesh->nodes(), _mesh->faces(), pointDOFs(), pointDOFsOffsets());
}

void AdvectionDiffusion3D::prepare()
{
	_mesh->loadProperty(_configuration.translation_motions, { "X", "Y", "Z" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z });
	_mesh->loadMaterials(_configuration.materials, _configuration.material_set);

	_mesh->addPropertyGroup({ Property::FLUX_X, Property::FLUX_Y, Property::FLUX_Z });
	_mesh->addPropertyGroup({ Property::GRADIENT_X, Property::GRADIENT_Y, Property::GRADIENT_Z });

	AdvectionDiffusion::prepare();
}

std::vector<std::pair<ElementType, Property> > AdvectionDiffusion3D::propertiesToStore() const
{
	for (size_t s = 0; s < _mesh->steps(); s++) {
		if (
				_mesh->isPropertyTimeDependent(Property::TRANSLATION_MOTION_X, s) ||
				_mesh->isPropertyTimeDependent(Property::TRANSLATION_MOTION_Y, s) ||
				_mesh->isPropertyTimeDependent(Property::TRANSLATION_MOTION_Z, s)) {

			return {
				{ ElementType::NODES, Property::TRANSLATION_MOTION_X },
				{ ElementType::NODES, Property::TRANSLATION_MOTION_Y },
				{ ElementType::NODES, Property::TRANSLATION_MOTION_Z }
			};
		}
	}

	return {};
}

void AdvectionDiffusion3D::assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K, DenseMatrix &CD, bool tangentCorrection) const
{
	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Point sin, cos;
	switch (material->coordination().type) {
	case MaterialCoordination::Type::CARTESIAN:

		cos.x = std::cos(d2r(material->coordination().rotation[0]->evaluate(e->node(node))));
		cos.y = std::cos(d2r(material->coordination().rotation[1]->evaluate(e->node(node))));
		cos.z = std::cos(d2r(material->coordination().rotation[2]->evaluate(e->node(node))));

		sin.x = std::sin(d2r(material->coordination().rotation[0]->evaluate(e->node(node))));
		sin.y = std::sin(d2r(material->coordination().rotation[1]->evaluate(e->node(node))));
		sin.z = std::sin(d2r(material->coordination().rotation[2]->evaluate(e->node(node))));

		break;

	case MaterialCoordination::Type::CYLINDRICAL: {

		Point origin(material->coordination().center[0]->evaluate(e->node(node)), material->coordination().center[1]->evaluate(e->node(node)), 0);
		const Point &p = _mesh->coordinates()[e->node(node)];
		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));

		cos.x = 1.0;
		cos.y = 1.0;
		cos.z = std::cos(rotation);

		sin.x = 0.0;
		sin.y = 0.0;
		sin.z = std::sin(rotation);

	} break;

	case MaterialCoordination::Type::SPHERICAL: {

		Point origin(
				material->coordination().center[0]->evaluate(e->node(node)),
				material->coordination().center[1]->evaluate(e->node(node)),
				material->coordination().center[2]->evaluate(e->node(node)));

		const Point &p = _mesh->coordinates()[e->node(node)];

		double azimut = std::atan2((p.y - origin.y), (p.x - origin.x));
		double r = std::sqrt(pow((p.x - origin.x), 2) + pow((p.y - origin.y), 2) + pow((p.z - origin.z), 2));
		double elevation = 0.0;

		if (r < 1e-12) {
			elevation = 0.0;
		} else {
			elevation = std::atan2(std::sqrt(pow((p.z - origin.z), 2) + pow((p.x - origin.x), 2)), (p.y - origin.y));
		}

		cos.x = 1.0;
		cos.y = std::cos(elevation);
		cos.z = std::cos(azimut);

		sin.x = 0.0;
		sin.y = std::sin(elevation);
		sin.z = std::sin(azimut);

	} break;

	}


	DenseMatrix TCT(3, 3), T(3, 3), C(3, 3), _CD, TCDT;

	T(0, 0) = cos.y * cos.z;                         T(0, 1) = cos.y * sin.z;                         T(0, 2) = -sin.y;
	T(1, 0) = cos.z * sin.x * sin.y - cos.x * sin.z; T(1, 1) = cos.x * cos.z + sin.x * sin.y * sin.z; T(1, 2) = cos.y * sin.x;
	T(2, 0) = sin.x*sin.z+cos.x*cos.z*sin.y;         T(2, 1) = cos.x * sin.y*sin.z-cos.z*sin.x;       T(2, 2) = cos.x * cos.y;

	if (tangentCorrection) {
		_CD.resize(3, 3);
		TCDT.resize(3, 3);
	}

	auto derivation = [&] (MATERIAL_PARAMETER p, double h) {
		return (
				material->get(p)->evaluate(e->node(node), step.currentTime, temp + h) -
				material->get(p)->evaluate(e->node(node), step.currentTime, temp - h)
				) / (2 * h);
	};

	switch (material->getModel(PHYSICS::ADVECTION_DIFFUSION_3D)) {
	case MATERIAL_MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = C(2, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = _CD(2, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case MATERIAL_MODEL::DIAGONAL:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.currentTime, temp);
		C(2, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(2, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ, temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case MATERIAL_MODEL::SYMMETRIC:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.currentTime, temp);
		C(2, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 1) = C(1, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 2) = C(2, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 2) = C(2, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ)->evaluate(e->node(node), step.currentTime, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(2, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ, temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY, temp / 1e4);
			_CD(0, 2) = _CD(2, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ, temp / 1e4);
			_CD(1, 2) = _CD(2, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ, temp / 1e4);
		}
		break;
	case MATERIAL_MODEL::ANISOTROPIC:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.currentTime, temp);
		C(2, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.currentTime, temp);
		C(0, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX)->evaluate(e->node(node), step.currentTime, temp);
		C(1, 2) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ)->evaluate(e->node(node), step.currentTime, temp);
		C(2, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZX)->evaluate(e->node(node), step.currentTime, temp);
		C(2, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZY)->evaluate(e->node(node), step.currentTime, temp);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(2, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ, temp / 1e4);
			_CD(0, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY, temp / 1e4);
			_CD(0, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ, temp / 1e4);
			_CD(1, 2) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ, temp / 1e4);
			_CD(1, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX, temp / 1e4);
			_CD(2, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZX, temp / 1e4);
			_CD(2, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZY, temp / 1e4);
		}
		break;
	default:
		ESINFO(ERROR) << "Advection diffusion 3D not supports set material model";
	}

	TCT.multiply(T, C * T, 1, 0, true, false);
	if (tangentCorrection) {
		TCDT.multiply(T, _CD * T, 1, 0, true, false);
		CD(node, 0) = TCDT(0, 0);
		CD(node, 1) = TCDT(1, 1);
		CD(node, 2) = TCDT(2, 2);
		CD(node, 3) = TCDT(0, 1);
		CD(node, 4) = TCDT(0, 2);
		CD(node, 5) = TCDT(1, 0);
		CD(node, 6) = TCDT(1, 2);
		CD(node, 7) = TCDT(2, 0);
		CD(node, 8) = TCDT(2, 1);
	}

	K(node, 0) = TCT(0, 0);
	K(node, 1) = TCT(1, 1);
	K(node, 2) = TCT(2, 2);
	K(node, 3) = TCT(0, 1);
	K(node, 4) = TCT(0, 2);
	K(node, 5) = TCT(1, 0);
	K(node, 6) = TCT(1, 2);
	K(node, 7) = TCT(2, 0);
	K(node, 8) = TCT(2, 1);
}

void AdvectionDiffusion3D::processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	bool CAU = _configuration.stabilization == AdvectionDiffusion3DConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = (matrices & Matrices::K) && step.tangentMatrixCorrection;

	DenseMatrix Ce(3, 3), coordinates(e->nodes(), 3), J(3, 3), invJ(3, 3), dND;
	double detJ, temp;
	DenseMatrix f(e->nodes(), 1);
	DenseMatrix U(e->nodes(), 3);
	DenseMatrix m(e->nodes(), 1);
	DenseMatrix T(e->nodes(), 1);
	DenseMatrix K(e->nodes(), 9);
	DenseMatrix gpK(1, 9), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;

	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	if (tangentCorrection) {
		CD.resize(e->nodes(), 9);
		CDe.resize(3, 3);
	}

	for (size_t i = 0; i < e->nodes(); i++) {
		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
		T(i, 0) = temp;
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
		coordinates(i, 2) = _mesh->coordinates()[e->node(i)].z;
		m(i, 0) =
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.currentTime, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(e->node(i), step.currentTime, temp);

		U(i, 0) = e->getProperty(Property::TRANSLATION_MOTION_X, i, step.step, step.currentTime, temp, 0) * m(i, 0);
		U(i, 1) = e->getProperty(Property::TRANSLATION_MOTION_Y, i, step.step, step.currentTime, temp, 0) * m(i, 0);
		U(i, 2) = e->getProperty(Property::TRANSLATION_MOTION_Z, i, step.step, step.currentTime, temp, 0) * m(i, 0);

		f(i, 0) = e->sumProperty(Property::HEAT_SOURCE, i, step.step, step.currentTime, temp, 0);
		assembleMaterialMatrix(step, e, i, temp, K, CD, tangentCorrection);
	}

	eslocal Ksize = e->nodes();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if ((matrices & Matrices::K) || ((matrices & Matrices::R) && step.timeIntegrationConstantK != 0)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}
	if ((matrices & Matrices::M) || ((matrices & Matrices::R) && step.timeIntegrationConstantM != 0)) {
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

	if (tangentCorrection) {
		tangentK.resize(Ksize, Ksize);
	}

	DenseMatrix u(1, 3), v(1, 3), re(1, e->nodes());
	double normGradN = 0;

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		u.multiply(e->N()[gp], U, 1, 0);

		J.multiply(e->dN()[gp], coordinates);
		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(e->N()[gp], K);
		if (tangentCorrection) {
			gpCD.multiply(e->N()[gp], CD);
			CDe(0, 0) = gpCD(0, 0);
			CDe(1, 1) = gpCD(0, 1);
			CDe(2, 2) = gpCD(0, 2);
			CDe(0, 1) = gpCD(0, 3);
			CDe(0, 2) = gpCD(0, 4);
			CDe(1, 0) = gpCD(0, 5);
			CDe(1, 2) = gpCD(0, 6);
			CDe(2, 0) = gpCD(0, 7);
			CDe(2, 1) = gpCD(0, 8);
		}
		gpM.multiply(e->N()[gp], m);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);
		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);
		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, e->dN()[gp]);

		DenseMatrix b_e(1, e->nodes()), b_e_c(1, e->nodes());
		b_e.multiply(u, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < re.columns(); i++) {
					re(0, i) = b_e(0, i) - f(0, i);
				}
				DenseMatrix ReBt(1, 3);
				ReBt.multiply(re, dND, 1 / pow(normGradN, 2), 0, false, true);
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

				double konst1 = re.norm() / normGradN;
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
		Ce(2, 2) += _configuration.sigma * h_e * norm_u_e;

		if (matrices & (Matrices::M | Matrices::R)) {
			Me.multiply(e->N()[gp], e->N()[gp], detJ * gpM(0, 0) * e->weighFactor()[gp], 1, true);
		}
		if (matrices & (Matrices::K | Matrices::R)) {
			if (tangentCorrection) {
				BT.multiply(dND, T);
				BTN.multiply(BT, e->N()[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * e->weighFactor()[gp], 1, true);
			}
			Ke.multiply(dND, Ce * dND, detJ * e->weighFactor()[gp], 1, true);
			Ke.multiply(e->N()[gp], b_e, detJ * e->weighFactor()[gp], 1, true);
			if (konst * e->weighFactor()[gp] * detJ != 0) {
				Ke.multiply(b_e, b_e, konst * e->weighFactor()[gp] * detJ, 1, true);
			}
			if (CAU) {
				Ke.multiply(dND, dND, C_e * e->weighFactor()[gp] * detJ, 1, true);
			}
		}

		if (matrices & Matrices::f) {
			for (eslocal i = 0; i < Ksize; i++) {
				fe(i, 0) += detJ * e->weighFactor()[gp] * e->N()[gp](0, i) * f(i, 0);
				if (norm_u_e != 0) {
					fe(i, 0) += detJ * e->weighFactor()[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if (matrices & Matrices::R) {
		Re.multiply(Ke, T, step.timeIntegrationConstantK, 0);
		Re.multiply(Me, T, step.timeIntegrationConstantM, 1);
		if (!(matrices & Matrices::K)) {
			Ke.resize(0, 0);
		}
		if (!(matrices & Matrices::M)) {
			Me.resize(0, 0);
		}
	}

	if (tangentCorrection) {
		Ke += tangentK;
	}
}

void AdvectionDiffusion3D::processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step) ||
		e->hasProperty(Property::HEAT_FLOW, step.step) ||
		e->hasProperty(Property::HEAT_FLUX, step.step))) {

		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	DenseMatrix coordinates(e->nodes(), 3), dND(1, 3), q(e->nodes(), 1), htc(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

	double area = 1, temp;
	eslocal Ksize = e->nodes();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t r = 0; r < e->regions().size(); r++) {
		if (step.step < e->regions()[r]->settings.size() && e->regions()[r]->settings[step.step].count(Property::HEAT_FLOW)) {
			area = e->regions()[r]->area;
			break;
		}
	}
	if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}

	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	AdvectionDiffusionConvection *convection = NULL;
	auto stepitc = _configuration.convection.find(step.step + 1);
	if (stepitc != _configuration.convection.end()) {
		for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
			auto regionit = stepitc->second.find(e->regions()[r]->name);
			if (regionit != stepitc->second.end()) {
				convection = regionit->second;
			}
		}
	}

	for (size_t n = 0; n < e->nodes(); n++) {
		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
		coordinates(n, 2) = _mesh->coordinates()[e->node(n)].z;

		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, n, step, temp) : 0;

		if (step.iteration) {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0) - temp);
		} else {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0));
		}

		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, n, step.step, step.currentTime, temp, 0);
		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0), 4) - pow(temp, 4));
		q(n, 0) += e->getProperty(Property::HEAT_FLOW, n, step.step, step.currentTime, temp, 0) / area;
		q(n, 0) += e->getProperty(Property::HEAT_FLUX, n, step.step, step.currentTime, temp, 0);

		emiss(n, 0) *= 4 * temp * temp * temp;
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();

		gpQ.multiply(N[gp], q);
		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
			gpHtc.multiply(N[gp], htc);
			gpEmiss.multiply(N[gp], emiss);

			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
		}
	}
}

void AdvectionDiffusion3D::processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step) ||
		e->hasProperty(Property::HEAT_FLOW, step.step) ||
		e->hasProperty(Property::HEAT_FLUX, step.step))) {

		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}
	if (!(matrices & (Matrices::K | Matrices::f))) {
		Ke.resize(0, 0);
		Me.resize(0, 0);
		Re.resize(0, 0);
		fe.resize(0, 0);
		return;
	}

	DenseMatrix coordinates(e->nodes(), 3), dND(1, 3), q(e->nodes(), 1), htc(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

	double area = 1, temp;
	eslocal Ksize = e->nodes();
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);

	if (matrices & Matrices::f) {
		fe.resize(Ksize, 1);
		fe = 0;
	}

	for (size_t r = 0; r < e->regions().size(); r++) {
		if (step.step < e->regions()[r]->settings.size() && e->regions()[r]->settings[step.step].count(Property::HEAT_FLOW)) {
			area = e->regions()[r]->area;
			break;
		}
	}
	if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
		Ke.resize(Ksize, Ksize);
		Ke = 0;
	}

	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();

	AdvectionDiffusionConvection *convection = NULL;
	auto stepitc = _configuration.convection.find(step.step + 1);
	if (stepitc != _configuration.convection.end()) {
		for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
			auto regionit = stepitc->second.find(e->regions()[r]->name);
			if (regionit != stepitc->second.end()) {
				convection = regionit->second;
			}
		}
	}

	for (size_t n = 0; n < e->nodes(); n++) {
		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
		coordinates(n, 2) = _mesh->coordinates()[e->node(n)].z;

		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, n, step, temp) : 0;

		if (step.iteration) {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0) - temp);
		} else {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0));
		}

		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, n, step.step, step.currentTime, temp, 0);
		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, step.currentTime, temp, 0), 4) - pow(temp, 4));
		q(n, 0) += e->getProperty(Property::HEAT_FLOW, n, step.step, step.currentTime, temp, 0) / area;
		q(n, 0) += e->getProperty(Property::HEAT_FLUX, n, step.step, step.currentTime, temp, 0);

		emiss(n, 0) *= 4 * temp * temp * temp;
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
			gpHtc.multiply(N[gp], htc);
			gpEmiss.multiply(N[gp], emiss);

			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
		}
	}
}

void AdvectionDiffusion3D::processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
}

void AdvectionDiffusion3D::postProcessElement(const Step &step, const Element *e, std::vector<Solution*> &solution)
{
	DenseMatrix Ce(3, 3), coordinates, J(3, 3), invJ(3, 3), dND, temp(e->nodes(), 1);
	double detJ, m, norm_u_e, h_e;
	DenseMatrix U(e->nodes(), 3), K(e->nodes(), 9), gpK(1, 9), CD;
	DenseMatrix u(1, 3), matFlux(3, 1), matGradient(3, 1);

	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	coordinates.resize(e->nodes(), 3);

	for (size_t i = 0; i < e->nodes(); i++) {
		temp(i, 0) = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
		coordinates(i, 2) = _mesh->coordinates()[e->node(i)].z;
		m =
			material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.currentTime, temp(i, 0)) *
			material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(e->node(i), step.currentTime, temp(i, 0));

		U(i, 0) = e->getProperty(Property::TRANSLATION_MOTION_X, i, step.step, step.currentTime, temp(i, 0), 0) * m;
		U(i, 1) = e->getProperty(Property::TRANSLATION_MOTION_Y, i, step.step, step.currentTime, temp(i, 0), 0) * m;
		U(i, 2) = e->getProperty(Property::TRANSLATION_MOTION_Z, i, step.step, step.currentTime, temp(i, 0), 0) * m;

		assembleMaterialMatrix(step, e, i, temp(i, 0), K, CD, false);
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		u.multiply(e->N()[gp], U, 1, 0);

		J.multiply(e->dN()[gp], coordinates);
		detJ = determinant3x3(J.values());
		inverse3x3(J.values(), invJ.values(), detJ);

		gpK.multiply(e->N()[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);
		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);
		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, e->dN()[gp]);

		 norm_u_e = u.norm();
		 h_e = 0;

		if (norm_u_e != 0) {
			DenseMatrix b_e(1, e->nodes());
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += _configuration.sigma * h_e * norm_u_e;
		Ce(1, 1) += _configuration.sigma * h_e * norm_u_e;
		Ce(2, 2) += _configuration.sigma * h_e * norm_u_e;

		matGradient.multiply(dND, temp, 1, 1);
		matFlux.multiply(Ce, dND * temp, 1, 1);
	}

	solution[offset + SolutionIndex::GRADIENT]->data[e->domains().front()].push_back(matGradient(0, 0) / e->gaussePoints());
	solution[offset + SolutionIndex::GRADIENT]->data[e->domains().front()].push_back(matGradient(1, 0) / e->gaussePoints());
	solution[offset + SolutionIndex::GRADIENT]->data[e->domains().front()].push_back(matGradient(2, 0) / e->gaussePoints());

	solution[offset + SolutionIndex::FLUX]->data[e->domains().front()].push_back(matFlux(0, 0) / e->gaussePoints());
	solution[offset + SolutionIndex::FLUX]->data[e->domains().front()].push_back(matFlux(1, 0) / e->gaussePoints());
	solution[offset + SolutionIndex::FLUX]->data[e->domains().front()].push_back(matFlux(2, 0) / e->gaussePoints());
}

void AdvectionDiffusion3D::processSolution(const Step &step)
{
	if (!_configuration.post_process) {
		return;
	}

	if (_instance->solutions[offset + SolutionIndex::GRADIENT] == NULL) {
		_instance->solutions[offset + SolutionIndex::GRADIENT] = new Solution(*_mesh, "gradient", ElementType::ELEMENTS, { Property::GRADIENT_X, Property::GRADIENT_Y, Property::GRADIENT_Z });
	}
	if (_instance->solutions[offset + SolutionIndex::FLUX] == NULL) {
		_instance->solutions[offset + SolutionIndex::FLUX] = new Solution(*_mesh, "flux", ElementType::ELEMENTS, { Property::FLUX_X, Property::FLUX_Y, Property::FLUX_Z });
	}

	for (size_t p = 0; p < _mesh->parts(); p++) {
		_instance->solutions[offset + SolutionIndex::GRADIENT]->data[p].clear();
		_instance->solutions[offset + SolutionIndex::FLUX]->data[p].clear();
	}

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		for (size_t e = _mesh->getPartition()[p]; e < _mesh->getPartition()[p + 1]; e++) {
			postProcessElement(step, _mesh->elements()[e], _instance->solutions);
		}
	}
}
