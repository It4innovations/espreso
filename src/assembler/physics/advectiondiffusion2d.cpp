
#include "../../configuration/physics/advectiondiffusion2d.h"
#include "advectiondiffusion2d.h"

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

using namespace espreso;

NewAdvectionDiffusion2D::NewAdvectionDiffusion2D(Mesh *mesh, Instance *instance, const AdvectionDiffusion2DConfiguration &configuration)
: Physics2D("ADVECTION DIFFUSION 2D", mesh, instance), _configuration(configuration)
{

}

MatrixType NewAdvectionDiffusion2D::getMatrixType(const Step &step, size_t domain) const
{
	if (_configuration.tangent_matrix_correction) {
		return MatrixType::REAL_UNSYMMETRIC;
	}
	if (_mesh->hasProperty(domain, Property::TRANSLATION_MOTION_X, step.step) || _mesh->hasProperty(domain, Property::TRANSLATION_MOTION_Y, step.step)) {
		return MatrixType::REAL_UNSYMMETRIC;
	} else {
		return MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	}
}

void NewAdvectionDiffusion2D::prepareTotalFETI()
{
	_instance->domainDOFCount = _mesh->assignUniformDOFsIndicesToNodes(_instance->domainDOFCount, pointDOFs(), _nodesDOFsOffsets);
	_instance->properties = pointDOFs();
	_mesh->computeNodesDOFsCounters(pointDOFs());

	_mesh->loadProperty(_configuration.initial_temperature, { }, { Property::INITIAL_TEMPERATURE });

	_mesh->loadNodeProperty(_configuration.temperature    , { }, { Property::TEMPERATURE });
	_mesh->loadNodeProperty(_configuration.thickness      , { }, { Property::THICKNESS });

	_mesh->loadProperty(_configuration.translation_motions, { "X", "Y" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y });
	_mesh->loadProperty(_configuration.heat_source        , { }         , { Property::HEAT_SOURCE });
	_mesh->loadProperty(_configuration.heat_flux          , { }         , { Property::HEAT_FLUX });
	_mesh->loadProperty(_configuration.heat_flow          , { }         , { Property::HEAT_FLOW });

	for (auto it = _configuration.convection.begin(); it != _configuration.convection.end(); ++it) {
		for (auto regions = it->second.begin(); regions != it->second.end(); ++regions) {
			std::map<std::string, std::string> values;

			values[regions->first] = regions->second->external_temperature;
			_mesh->loadProperty(values, { }, { Property::EXTERNAL_TEMPERATURE });

			switch (regions->second->type) {
			case CONVECTION_TYPE::USER:
				values[regions->first] = regions->second->heat_transfer_coefficient;
				_mesh->loadProperty(values, { }, { Property::HEAT_TRANSFER_COEFFICIENT });
				break;
			case CONVECTION_TYPE::EXTERNAL_NATURAL:

				switch (regions->second->variant) {
				case CONVECTION_VARIANT::VERTICAL_WALL:
					values[regions->first] = regions->second->wall_height;
					_mesh->loadProperty(values, { }, { Property::WALL_HEIGHT });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				case CONVECTION_VARIANT::INCLINED_WALL:
					values[regions->first] = regions->second->wall_height;
					_mesh->loadProperty(values, { }, { Property::WALL_HEIGHT });
					values[regions->first] = regions->second->tilt_angle;
					_mesh->loadProperty(values, { }, { Property::TILT_ANGLE });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				case CONVECTION_VARIANT::HORIZONTAL_CYLINDER:
					values[regions->first] = regions->second->diameter;
					_mesh->loadProperty(values, { }, { Property::DIAMETER });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				case CONVECTION_VARIANT::SPHERE:
					values[regions->first] = regions->second->diameter;
					_mesh->loadProperty(values, { }, { Property::DIAMETER });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				case CONVECTION_VARIANT::HORIZONTAL_PLATE_UP:
				case CONVECTION_VARIANT::HORIZONTAL_PLATE_DOWN:
					values[regions->first] = regions->second->length;
					_mesh->loadProperty(values, { }, { Property::LENGTH });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				default:
					break;
				}

				break;
			case CONVECTION_TYPE::EXTERNAL_FORCED:

				switch (regions->second->variant) {
				case CONVECTION_VARIANT::AVERAGE_PLATE:
					values[regions->first] = regions->second->length;
					_mesh->loadProperty(values, { }, { Property::LENGTH });
					values[regions->first] = regions->second->fluid_velocity;
					_mesh->loadProperty(values, { }, { Property::FLUID_VELOCITY });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				default:
					break;
				}

				break;
			case CONVECTION_TYPE::INTERNAL_NATURAL:

				switch (regions->second->variant) {
				case CONVECTION_VARIANT::CIRCULAR_TUBE:
					values[regions->first] = regions->second->diameter;
					_mesh->loadProperty(values, { }, { Property::DIAMETER });
					values[regions->first] = regions->second->wall_height;
					_mesh->loadProperty(values, { }, { Property::WALL_HEIGHT });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				default:
					break;
				}

				break;
			case CONVECTION_TYPE::INTERNAL_FORCED:

				switch (regions->second->variant) {
				case CONVECTION_VARIANT::TUBE:
					values[regions->first] = regions->second->diameter;
					_mesh->loadProperty(values, { }, { Property::DIAMETER });
					values[regions->first] = regions->second->fluid_velocity;
					_mesh->loadProperty(values, { }, { Property::FLUID_VELOCITY });
					values[regions->first] = regions->second->absolute_pressure;
					_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
					break;
				default:
					break;
				}

				break;

			}

		}
	}

	for (auto it = _configuration.diffuse_radiation.begin(); it != _configuration.diffuse_radiation.end(); ++it) {
		for (auto regions = it->second.begin(); regions != it->second.end(); ++regions) {
			std::map<std::string, std::string> values;
			values[regions->first] = regions->second->external_temperature;
			_mesh->loadProperty(values, { }, { Property::EXTERNAL_TEMPERATURE });
			values[regions->first] = regions->second->emissivity;
			_mesh->loadProperty(values, { }, { Property::EMISSIVITY });
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

void NewAdvectionDiffusion2D::assembleMaterialMatrix(const Step &step, const Element *e, eslocal node, double temp, DenseMatrix &K, DenseMatrix &CD) const
{
	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	double cos, sin;
	switch (material->coordination().type) {
	case MaterialCoordination::Type::CARTESIAN:
		cos = std::cos(d2r(material->coordination().rotation[2]->evaluate(e->node(node))));
		sin = std::sin(d2r(material->coordination().rotation[2]->evaluate(e->node(node))));
		break;
	case MaterialCoordination::Type::CYLINDRICAL: {
		Point origin(material->coordination().center[0]->evaluate(e->node(node)), material->coordination().center[1]->evaluate(e->node(node)), 0);
		const Point &p = _mesh->coordinates()[e->node(node)];
		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));
		cos = std::cos(rotation);
		sin = std::sin(rotation);
		break;
	}
	}

	DenseMatrix TCT(2, 2), T(2, 2), C(2, 2), _CD, TCDT;
	T(0, 0) =  cos; T(0, 1) = sin;
	T(1, 0) = -sin; T(1, 1) = cos;

	if (_configuration.tangent_matrix_correction) {
		_CD.resize(2, 2);
		TCDT.resize(2, 2);
	}

	auto derivation = [&] (MATERIAL_PARAMETER p, double h) {
		return (
				material->get(p)->evaluate(e->node(node), step.step, temp + h) -
				material->get(p)->evaluate(e->node(node), step.step, temp - h)
				) / (2 * h);


	};

	switch (material->getModel(PHYSICS::ADVECTION_DIFFUSION_2D)) {
	case MATERIAL_MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.step, temp);
		C(0, 1) = C(1, 0) = 0;
		if (_configuration.tangent_matrix_correction) {
			_CD(0, 0) = _CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case MATERIAL_MODEL::DIAGONAL:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.step, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.step, temp);
		C(0, 1) = C(1, 0) = 0;
		if (_configuration.tangent_matrix_correction) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case MATERIAL_MODEL::SYMMETRIC:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.step, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.step, temp);
		C(1, 0) = C(0, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.step, temp);
		if (_configuration.tangent_matrix_correction) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY, temp / 1e4);
		}
		break;
	case MATERIAL_MODEL::ANISOTROPIC:
		C(0, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX)->evaluate(e->node(node), step.step, temp);
		C(1, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY)->evaluate(e->node(node), step.step, temp);
		C(0, 1) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY)->evaluate(e->node(node), step.step, temp);
		C(1, 0) = material->get(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX)->evaluate(e->node(node), step.step, temp);
		if (_configuration.tangent_matrix_correction) {
			_CD(0, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, temp / 1e4);
			_CD(1, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, temp / 1e4);
			_CD(0, 1) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY, temp / 1e4);
			_CD(1, 0) = derivation(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX, temp / 1e4);
		}
		break;
	default:
		ESINFO(ERROR) << "Advection diffusion 2D not supports set material model";
	}

	TCT.multiply(T, C * T, 1, 0, true, false);
	if (_configuration.tangent_matrix_correction) {
		TCDT.multiply(T, _CD * T, 1, 0, true, false);
		CD(node, 0) = TCDT(0, 0);
		CD(node, 1) = TCDT(1, 1);
		CD(node, 2) = TCDT(0, 1);
		CD(node, 3) = TCDT(1, 0);
	}

	K(node, 0) = TCT(0, 0);
	K(node, 1) = TCT(1, 1);
	K(node, 2) = TCT(0, 1);
	K(node, 3) = TCT(1, 0);
}

void NewAdvectionDiffusion2D::processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	bool CAU = _configuration.stabilization == AdvectionDiffusion2DConfiguration::STABILIZATION::CAU;

	DenseMatrix Ce(2, 2), coordinates, J(2, 2), invJ(2, 2), dND;
	double detJ, temp;
	DenseMatrix f(e->nodes(), 1);
	DenseMatrix U(e->nodes(), 2);
	DenseMatrix m(e->nodes(), 1);
	DenseMatrix T(e->nodes(), 1);
	DenseMatrix thickness(e->nodes(), 1), K(e->nodes(), 4);
	DenseMatrix gpThickness(1, 1), gpK(1, 4), gpM(1, 1);
	DenseMatrix tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;

	const Material* material = _mesh->materials()[e->param(Element::MATERIAL)];

	coordinates.resize(e->nodes(), 2);

	if (_configuration.tangent_matrix_correction) {
		CD.resize(e->nodes(), 4);
		CDe.resize(2, 2);
	}

	for (size_t i = 0; i < e->nodes(); i++) {
		if (solution.size()) {
			temp = solution[0]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(i), e->domains().front()));
		} else {
			temp = e->getProperty(Property::INITIAL_TEMPERATURE, i, step.step, 273.15 + 20);
		}
		T(i, 0) = temp;
		coordinates(i, 0) = _mesh->coordinates()[e->node(i)].x;
		coordinates(i, 1) = _mesh->coordinates()[e->node(i)].y;
		thickness(i, 0) = e->getProperty(Property::THICKNESS, i, step.step, 1);
		m(i, 0) =
				material->get(MATERIAL_PARAMETER::DENSITY)->evaluate(e->node(i), step.step, temp) *
				material->get(MATERIAL_PARAMETER::HEAT_CAPACITY)->evaluate(e->node(i), step.step, temp) *
				thickness(i, 0);

		U(i, 0) = e->getProperty(Property::TRANSLATION_MOTION_X, i, step.step, 0) * m(i, 0);
		U(i, 1) = e->getProperty(Property::TRANSLATION_MOTION_Y, i, step.step, 0) * m(i, 0);
		f(i, 0) = e->sumProperty(Property::HEAT_SOURCE, i, step.step, 0) * thickness(i, 0);
		assembleMaterialMatrix(step, e, i, temp, K, CD);
	}

	eslocal Ksize = e->nodes();

	Ke.resize(0, 0);
	Me.resize(0, 0);
	Re.resize(0, 0);
	fe.resize(0, 0);
	if (matrices & (Matrices::K | Matrices::R)) {
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

	if (_configuration.tangent_matrix_correction) {
		tangentK.resize(Ksize, Ksize);
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
		if (_configuration.tangent_matrix_correction) {
			gpCD.multiply(e->N()[gp], CD);
			CDe(0, 0) = gpCD(0, 0);
			CDe(1, 1) = gpCD(0, 1);
			CDe(0, 1) = gpCD(0, 2);
			CDe(1, 0) = gpCD(0, 3);
		}
		gpM.multiply(e->N()[gp], m);

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
				for (size_t i = 0; i < re.columns(); i++) {
					re(0, i) = b_e(0, i) - f(i, 0);
				}
				DenseMatrix ReBt(1, 2);
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
				DenseMatrix u_v(1, 2);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
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

		if (matrices & Matrices::M) {
			Me.multiply(e->N()[gp], e->N()[gp], detJ * gpM(0, 0) * e->weighFactor()[gp], 1, true);
		}
		if (matrices & (Matrices::K | Matrices::R)) {
			if (_configuration.tangent_matrix_correction) {
				BT.multiply(dND, T);
				BTN.multiply(BT, e->N()[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * e->weighFactor()[gp] * gpThickness(0, 0), 1, true);
			}
			Ke.multiply(dND, Ce * dND, detJ * e->weighFactor()[gp] * gpThickness(0, 0), 1, true);
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
		Re.multiply(Ke, T);
		if (!(matrices & Matrices::K)) {
			Ke.resize(0, 0);
		}
	}

	if ((matrices & Matrices::K) && _configuration.tangent_matrix_correction) {
		Ke += tangentK;
	}
}

void NewAdvectionDiffusion2D::processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{
	ESINFO(ERROR) << "Advection diffusion 2D cannot process face";
}

static void convectionMatParameters(const AdvectionDiffusionConvection &convection, const Element *e, size_t node, size_t step, double temp, double T_EXT, double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity )
{


	double  gas_constant;

	switch (convection.fluid) {
	case espreso::CONVECTION_FLUID::AIR:{


		gas_constant = 286.9;
		rho = (e->getProperty(Property::ABSOLUTE_PRESSURE, node, step, 0)) / (gas_constant * T_EXT);


		if ((T_EXT >=200) && (T_EXT <= 1600)){
			heat_capacity = 1047.63657-0.372589265*T_EXT+9.45304214E-4*pow(T_EXT,2.0)-6.02409443E-7*pow(T_EXT,3.0)+1.2858961E-10*pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			heat_capacity = 1047.63657-0.372589265*200.0+9.45304214E-4*pow(200.0,2.0)-6.02409443E-7*pow(200.0,3.0)+1.2858961E-10*pow(200.0,4.0);
		}else if (T_EXT > 1600){
			heat_capacity = 1047.63657-0.372589265*1600.0+9.45304214E-4*pow(1600.0,2.0)-6.02409443E-7*pow(1600.0,3.0)+1.2858961E-10*pow(1600.0,4.0);
		}

		if ((T_EXT >=200) && (T_EXT <= 1600)){
			thermal_conductivity = -0.00227583562+1.15480022E-4*T_EXT-7.90252856E-8*pow(T_EXT,2.0)+4.11702505E-11*pow(T_EXT,3.0)-7.43864331E-15*pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			thermal_conductivity = -0.00227583562+1.15480022E-4*200.0-7.90252856E-8*pow(200.0,2.0)+4.11702505E-11*pow(200.0,3.0)-7.43864331E-15*pow(200.0,4.0);
		}else if (T_EXT > 1600){
			thermal_conductivity =  -0.00227583562+1.15480022E-4*1600.0-7.90252856E-8*pow(1600.0,2.0)+4.11702505E-11*pow(1600.0,3.0)-7.43864331E-15*pow(1600.0,4.0);
		}

		if ((T_EXT >=200) && (T_EXT <= 1600)){
		dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * T_EXT - 7.69429583E-11 * pow(T_EXT,2.0) + 4.6437266E-14 * pow(T_EXT,3.0) - 1.06585607E-17 * pow(T_EXT,4.0);
		}else if (T_EXT < 200){
			dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (T_EXT > 1600){
			dynamic_viscosity = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}


		if ((temp >=200) && (temp <= 1600)){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * temp - 7.69429583E-11 * pow(temp,2.0) + 4.6437266E-14 * pow(temp,3.0) - 1.06585607E-17 * pow(temp,4.0);
		}else if (temp < 200){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * 200.0 - 7.69429583E-11 * pow(200.0,2.0) + 4.6437266E-14 * pow(200.0,3.0) - 1.06585607E-17 * pow(200.0,4.0);
		}else if (temp > 1600){
			dynamic_viscosity_T = -8.38278E-7 + 8.35717342E-8 * 1600.0 - 7.69429583E-11 * pow(1600.0,2.0) + 4.6437266E-14 * pow(1600.0,3.0) - 1.06585607E-17 * pow(1600.0,4.0);
		}

	}break;
	case espreso::CONVECTION_FLUID::WATER:{

		if ((T_EXT >=273) && (T_EXT <= 283)){
			rho = 972.7584 + 0.2084 *T_EXT - 4.0E-4 * pow(T_EXT,2.0);
		}else if((T_EXT >283) && (T_EXT <= 373)){
			rho = 345.28 + 5.749816 * T_EXT - 0.0157244 * pow(T_EXT,2.0) + 1.264375E-5 * pow(T_EXT,3.0);
		}else if (T_EXT < 273){
			rho = 972.7584 + 0.2084 *273.0 - 4.0E-4 * pow(273.0,2.0);
		}else if (T_EXT > 373){
			rho = 345.28 + 5.749816 * 373.0 - 0.0157244 * pow(373.0,2.0) + 1.264375E-5 * pow(373.0,3.0);
		}


		if ((T_EXT >= 265) && (T_EXT <= 293)){
			dynamic_viscosity = 5.948859 - 0.08236196 * T_EXT + 4.287142E-4 * pow(T_EXT,2.0) - 9.938045E-7 * pow(T_EXT,3.0) + 8.65316E-10 * pow(T_EXT,4.0);
		}else if((T_EXT >293) && (T_EXT <= 353)){
			dynamic_viscosity = 	0.410191 - 0.004753985 * T_EXT + 2.079795E-5 * pow(T_EXT,2.0) - 4.061698E-8 *  pow(T_EXT,3.0) + 2.983925E-11 * pow(T_EXT,4.0);
		}else if((T_EXT >353) && (T_EXT <= 423)){
			dynamic_viscosity = 0.03625638 - 3.265463E-4 * T_EXT + 1.127139E-6 * pow(T_EXT,2.0) - 1.75363E-9 * pow(T_EXT,3.0) + 1.033976E-12 * pow(T_EXT,4.0);
		}else if (T_EXT < 265){
			dynamic_viscosity = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (T_EXT > 423){
			dynamic_viscosity = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}


		if ((T_EXT >= 275) && (T_EXT <= 370)){
			thermal_conductivity = -0.9003748 + 0.008387698 * T_EXT - 1.118205E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 275){
			thermal_conductivity = -0.9003748 + 0.008387698 * 275.0 - 1.118205E-5 * pow(275.0,2.0);
		}else if (T_EXT > 370){
			thermal_conductivity = -0.9003748 + 0.008387698 * 370.0 - 1.118205E-5 * pow(370.0,2.0);
		}

		if ((T_EXT >= 293) && (T_EXT <= 373)){
			heat_capacity = 4035.841 + 0.492312 * T_EXT;
		}else if (T_EXT < 293){
			heat_capacity = 4035.841 + 0.492312 * 293.0;
		}else if (T_EXT > 373){
			heat_capacity = 4035.841 + 0.492312 * 373.0;
		}


		if ((temp >= 265) && (temp <= 293)){
			dynamic_viscosity_T = 5.948859 - 0.08236196 * temp + 4.287142E-4 * pow(temp,2.0) - 9.938045E-7 * pow(temp,3.0) + 8.65316E-10 * pow(temp,4.0);
		}else if((temp >293) && (temp <= 353)){
			dynamic_viscosity_T = 	0.410191 - 0.004753985 * temp + 2.079795E-5 * pow(temp,2.0) - 4.061698E-8 *  pow(temp,3.0) + 2.983925E-11 * pow(temp,4.0);
		}else if((temp >353) && (temp <= 423)){
			dynamic_viscosity_T = 0.03625638 - 3.265463E-4 * temp + 1.127139E-6 * pow(temp,2.0) - 1.75363E-9 * pow(temp,3.0) + 1.033976E-12 * pow(temp,4.0);
		}else if (temp < 265){
			dynamic_viscosity_T = 5.948859 - 0.08236196 * 265.0 + 4.287142E-4 * pow(265.0,2.0) - 9.938045E-7 * pow(265.0,3.0) + 8.65316E-10 * pow(265.0,4.0);
		}else if (temp > 423){
			dynamic_viscosity_T = 0.03625638 - 3.265463E-4 * 423.0 + 1.127139E-6 * pow(423.0,2.0) - 1.75363E-9 * pow(423.0,3.0) + 1.033976E-12 * pow(423.0,4.0);
		}

	}break;
	case espreso::CONVECTION_FLUID::ENGINE_OIL:{


		if ((T_EXT >=273) && (T_EXT <= 433)){
			rho = 1068.70404 - 0.6393421 * T_EXT + 7.34307359E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			rho = 1068.70404 - 0.6393421 * 273.0 + 7.34307359E-5 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			rho = 1068.70404 - 0.6393421 * 433.0 + 7.34307359E-5 * pow(433.0,2.0);
		}

		if ((T_EXT >= 273) && (T_EXT <= 433)){
			thermal_conductivity = 0.192223542 - 2.0637987E-4 * T_EXT + 1.54220779E-7 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			thermal_conductivity = 0.192223542 - 2.0637987E-4 * 273.0 + 1.54220779E-7 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			thermal_conductivity =0.192223542 - 2.0637987E-4 * 433.0 + 1.54220779E-7 * pow(433.0,2.0);
		}


		if ((T_EXT >= 273) && (T_EXT <= 433)){
			heat_capacity = 761.405625 + 3.47685606 * T_EXT + 0.00115530303 * pow(T_EXT,2.0);
		}else if (T_EXT < 273){
			heat_capacity = 761.405625 + 3.47685606 * 273.0 + 0.00115530303 * pow(273.0,2.0);
		}else if (T_EXT > 433){
			heat_capacity = 761.405625 + 3.47685606 * 433.0 + 0.00115530303 * pow(433.0,2.0);
		}


		if ((T_EXT >= 273) && (T_EXT <= 353)){
			dynamic_viscosity = 42669.28688622 - 741.1718801282 * T_EXT + 5.360521287088 * pow(T_EXT,2.0) - 0.02066027676164 * pow(T_EXT,3.0) + 4.47491538052E-5 * pow(T_EXT,4.0) - 5.164053479202E-8 * pow(T_EXT,5.0) + 2.48033770504E-11 * pow(T_EXT,6.0);
		}else if ((T_EXT > 353) && (T_EXT <= 433 )){
			dynamic_viscosity = 4.94593941 - 0.0351869631 * T_EXT + 8.37935977E-5 * pow(T_EXT,2.0) - 6.67125E-8 * pow(T_EXT,3.0);

		}else if (T_EXT < 273){
			dynamic_viscosity = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (T_EXT > 433){
			dynamic_viscosity = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}

		if ((temp >= 273) && (temp <= 353)){
			dynamic_viscosity_T = 42669.28688622 - 741.1718801282 * temp + 5.360521287088 * pow(temp,2.0) - 0.02066027676164 * pow(temp,3.0) + 4.47491538052E-5 * pow(temp,4.0) - 5.164053479202E-8 * pow(temp,5.0) + 2.48033770504E-11 * pow(temp,6.0);
		}else if ((temp > 353) && (temp <= 433 )){
			dynamic_viscosity_T = 4.94593941 - 0.0351869631 * temp + 8.37935977E-5 * pow(temp,2.0) - 6.67125E-8 * pow(temp,3.0);

		}else if (temp < 273){
			dynamic_viscosity_T = 42669.28688622 - 741.1718801282 * 273.0 + 5.360521287088 * pow(273.0,2.0) - 0.02066027676164 * pow(273.0,3.0) + 4.47491538052E-5 * pow(273.0,4.0) - 5.164053479202E-8 * pow(273.0,5.0) + 2.48033770504E-11 * pow(273.0,6.0);
		}else if (temp > 433){
			dynamic_viscosity_T = 4.94593941 - 0.0351869631 * 433.0 + 8.37935977E-5 * pow(433.0,2.0) - 6.67125E-8 * pow(433.0,3.0);
		}


	}break;
	case espreso::CONVECTION_FLUID::TRANSFORMER_OIL:{

		if ((T_EXT >=223) && (T_EXT <= 373)){
			rho = 1055.04607 - 0.581753034 * T_EXT - 6.40531689E-5 * pow(T_EXT,2.0);
		}else if (T_EXT < 223){
			rho =  1055.04607 - 0.581753034 * 223.0 - 6.40531689E-5 * pow(223.0,2.0);
		}else if (T_EXT > 373){
			rho = 1055.04607 - 0.581753034 * 373.0 - 6.40531689E-5 * pow(373.0,2.0);
		}

		if ((T_EXT >= 273) && (T_EXT <= 433)){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * T_EXT;
		}else if (T_EXT < 273){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * 223.0;
		}else if (T_EXT > 433){
			thermal_conductivity = 0.134299084 - 8.04973822E-5 * 373.0;
		}


		if ((T_EXT >= 223) && (T_EXT <= 293)){
			heat_capacity = -117056.38 + 1816.76208 * T_EXT - 10.305786 * pow(T_EXT,2.0) + 0.0256691919 * pow(T_EXT,3.0) - 2.36742424E-5 * pow(T_EXT,4.0);
		}else if ((T_EXT > 293) && (T_EXT <= 373 )){
			heat_capacity = -13408.1491 + 123.044152 * T_EXT - 0.335401786 * pow(T_EXT,2.0) + 3.125E-4 * pow(T_EXT,3.0);
		}else if (T_EXT < 223){
			heat_capacity = heat_capacity = -117056.38 + 1816.76208 * 223.0 - 10.305786 * pow(223.0,2.0) + 0.0256691919 * pow(223.0,3.0) - 2.36742424E-5 * pow(223.0,4.0);
		}else if (T_EXT > 373){
			heat_capacity = heat_capacity = -13408.1491 + 123.044152 * 373.0 - 0.335401786 * pow(373.0,2.0) + 3.125E-4 * pow(373.0,3.0);
		}


		if ((T_EXT >= 243) && (T_EXT <= 273)){
			dynamic_viscosity = 4492.20229 - 64.7408879 * T_EXT + 0.349900959 * pow(T_EXT,2.0) - 8.40477E-4 * pow(T_EXT,3.0) + 7.57041667E-7 * pow(T_EXT,4.0);
		}else if ((T_EXT > 273) && (T_EXT <= 373 )){
			dynamic_viscosity = 91.4524999 - 1.33227058 * T_EXT + 0.00777680216 * pow(T_EXT,2.0) - 2.27271368E-5 *  pow(T_EXT,3.0) + 3.32419673E-8 * pow(T_EXT,4.0) - 1.94631023E-11 * pow(T_EXT,5.0);
		}else if (T_EXT < 243){
			dynamic_viscosity = dynamic_viscosity = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (T_EXT > 373){
			dynamic_viscosity = dynamic_viscosity = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

		if ((temp >= 243) && (temp <= 273)){
			dynamic_viscosity_T = 4492.20229 - 64.7408879 * temp + 0.349900959 * pow(temp,2.0) - 8.40477E-4 * pow(temp,3.0) + 7.57041667E-7 * pow(temp,4.0);
		}else if ((temp > 273) && (temp <= 373 )){
			dynamic_viscosity_T = 91.4524999 - 1.33227058 * temp + 0.00777680216 * pow(temp,2.0) - 2.27271368E-5 *  pow(temp,3.0) + 3.32419673E-8 * pow(temp,4.0) - 1.94631023E-11 * pow(temp,5.0);
		}else if (temp < 243){
			dynamic_viscosity_T = dynamic_viscosity = 4492.20229 - 64.7408879 * 243.0 + 0.349900959 * pow(243.0,2.0) - 8.40477E-4 * pow(243.0,3.0) + 7.57041667E-7 * pow(243.0,4.0);
		}else if (temp > 373){
			dynamic_viscosity_T = dynamic_viscosity = 91.4524999 - 1.33227058 * 373.0 + 0.00777680216 * pow(373.0,2.0) - 2.27271368E-5 *  pow(373.0,3.0) + 3.32419673E-8 * pow(373.0,4.0) - 1.94631023E-11 * pow(373.0,5.0);

		}

	}break;
	default:
		ESINFO(ERROR) << "Invalid convection fluid type.";
	}


}


static double computeHTC(const AdvectionDiffusionConvection &convection, const Element *e, size_t node, size_t step, double temp)
{
//	e->getProperty(Property::HEAT_TRANSFER_COEFFICIENT, node, step, 0);
//	e->getProperty(Property::WALL_HEIGHT, node, step, 0);
//	e->getProperty(Property::TILT_ANGLE, node, step, 0);
//	e->getProperty(Property::DIAMETER, node, step, 0);
//	e->getProperty(Property::PLATE_LENGTH, node, step, 0);
//	e->getProperty(Property::FLUID_VELOCITY, node, step, 0);
//	e->getProperty(Property::PLATE_DISTANCE, node, step, 0);
//	e->getProperty(Property::LENGTH, node, step, 0);

	double htc = 0;
	switch (convection.type) {
	case espreso::CONVECTION_TYPE::USER:{
		htc = e->getProperty(Property::HEAT_TRANSFER_COEFFICIENT, node, step, 0);
	}break;
	case espreso::CONVECTION_TYPE::EXTERNAL_NATURAL:{

		double T_AVG, g, rho, dynamic_viscosity, heat_capacity, thermal_conductivity, dynamic_viscosity_T;

		T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
		g = 9.81;

		convectionMatParameters(convection, e, node, step, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

		switch (convection.variant) {
		case espreso::CONVECTION_VARIANT::INCLINED_WALL: {

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::WALL_HEIGHT, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double tilt_angle = e->getProperty(Property::TILT_ANGLE, node, step,	0) * M_PI / 180.0;
			if (RaL <= 1e9) {
				htc = (thermal_conductivity	/ e->getProperty(Property::WALL_HEIGHT, node, step, 0)) * (0.68 + (0.67 * cos(tilt_angle) * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity	/ e->getProperty(Property::WALL_HEIGHT, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		}break;
		case espreso::CONVECTION_VARIANT::VERTICAL_WALL: {

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::WALL_HEIGHT, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (RaL <= 1e9) {
				htc = (thermal_conductivity	/ e->getProperty(Property::WALL_HEIGHT, node, step, 0)) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity	/ e->getProperty(Property::WALL_HEIGHT, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		}break;
        case espreso::CONVECTION_VARIANT::HORIZONTAL_PLATE_UP:{

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (temp > e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step, 0)){

				if (RaL <= 1e7) {
					htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.54 * pow(RaL,0.25);
				}else{
					htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.15 * pow(RaL,1.0/3.0);
				}
			}else{
				htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.27 * pow(RaL,0.25);
			}

		}break;
		case espreso::CONVECTION_VARIANT::HORIZONTAL_PLATE_DOWN:{

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (temp <= e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step, 0)){

				if (RaL <= 1e7) {
					htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.54 * pow(RaL,0.25);
				}else{
					htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.15 * pow(RaL,1.0/3.0);
				}
			}else{
				htc = thermal_conductivity / e->getProperty(Property::LENGTH, node, step, 0) * 0.27 * pow(RaL,0.25);
			}
		}break;
		case espreso::CONVECTION_VARIANT::HORIZONTAL_CYLINDER:{

			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

			if ( RaD > 10e12 ){
				// warning!!!!
				ESINFO(ERROR) << "Validated only for RaD <= 10e12 ";
			}

			htc = thermal_conductivity / e->getProperty(Property::DIAMETER, node, step, 0) * pow( 0.6 + ( 0.387*pow(RaD,1.0/6.0)/ pow( 1 + pow( 0.559/Pr, 9.0/16.0), 8.0/27.0) ) ,2.0);

		}break;
		case espreso::CONVECTION_VARIANT::SPHERE:{

			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

			if ( RaD > 10e11 || Pr < 0.7 ){
				// warning!!!!
				ESINFO(ERROR) << "Validated only for RaD <= 10e11 and Pr >= 0.7 ";
			}

			htc = thermal_conductivity / e->getProperty(Property::DIAMETER, node, step, 0) * pow( 2.0 + ( 0.589*pow(RaD,0.25)/ pow( 1 + pow( 0.469/Pr, 9.0/16.0), 4.0/9.0) ) ,2.0);

		}break;
		default:
			ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_NATURAL.";
		}
	}break;

	case espreso::CONVECTION_TYPE::INTERNAL_NATURAL:{

		double T_AVG, g, rho, dynamic_viscosity, heat_capacity, thermal_conductivity,dynamic_viscosity_T;

		T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
		g = 9.81;

		convectionMatParameters(convection, e, node, step, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

		switch (convection.variant) {
		case espreso::CONVECTION_VARIANT::PARALLEL_PLATES: {

			double H_L = e->getProperty(Property::WALL_HEIGHT, node, step, 0) / e->getProperty(Property::LENGTH, node, step, 0);
			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (( RaL < H_L ) && (temp >  e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  )){

				htc = thermal_conductivity / e->getProperty(Property::WALL_HEIGHT, node, step, 0) * ( 1.0 / 24.0 ) * RaL;

			}else{

				double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

				if (RaL <= 1e9) {
					htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		}break;
		case espreso::CONVECTION_VARIANT::CIRCULAR_TUBE: {
			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double H_D = e->getProperty(Property::WALL_HEIGHT, node, step, 0) / e->getProperty(Property::DIAMETER, node, step, 0);

			if ( RaD < H_D ){
				htc = thermal_conductivity / e->getProperty(Property::WALL_HEIGHT, node, step, 0) * ( 1.0 / 128.0 ) * RaD;
			}else{

				double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * std::fabs( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
				if (RaD <= 1e9) {
					htc = (thermal_conductivity	/ e->getProperty(Property::DIAMETER, node, step, 0)) * (0.68 + (0.67 * pow(RaD,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermal_conductivity	/ e->getProperty(Property::DIAMETER, node, step, 0)) * pow(0.825 + (0.387 * pow(RaD,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		}break;
	   	default:
		   ESINFO(ERROR) << "Invalid convection variant for INTERNAL_NATURAL.";
		}
	}break;

	case espreso::CONVECTION_TYPE::EXTERNAL_FORCED:{

			switch (convection.variant) {
			case espreso::CONVECTION_VARIANT::AVERAGE_PLATE: {

				double T_AVG, g, rho, dynamic_viscosity, heat_capacity, thermal_conductivity,dynamic_viscosity_T;

				T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
				g = 9.81;

				convectionMatParameters(convection, e, node, step, temp, T_AVG, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

				double Re = rho	* e->getProperty(Property::FLUID_VELOCITY, node, step, 0) * e->getProperty(Property::LENGTH, node, step, 0)	/ dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;
				if (Re <= 5e5) {
					htc = 2	* (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * ((0.3387 * pow(Pr, 1.0 / 3.0) * pow(Re, 0.5)) / (pow(1 + pow(0.0468 / Pr, 2.0 / 3.0), 0.25)));
				} else {
					htc = 2	* (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * pow(Pr, 1.0 / 3.0)	* (0.037 * pow(Re, 0.8) - 871);
				}

			}break;
			default:
				ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_FORCED.";
			}
	}break;


	case espreso::CONVECTION_TYPE::INTERNAL_FORCED:{

			switch (convection.variant) {
			case espreso::CONVECTION_VARIANT::TUBE: {

				double T_EXT, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity;

				convectionMatParameters(convection, e, node, step, temp, T_EXT, rho, dynamic_viscosity, dynamic_viscosity_T, heat_capacity, thermal_conductivity );

				double Re = rho	* e->getProperty(Property::FLUID_VELOCITY, node, step, 0) * e->getProperty(Property::DIAMETER, node, step, 0) / dynamic_viscosity;
				double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;
				double n = temp	< e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step, 0) ? 0.3 : 0.4;
				htc = thermal_conductivity / e->getProperty(Property::DIAMETER, node, step, 0);
				if (Re <= 2500) {
					htc *= 3.66;
				} else {
					htc *= 0.027 * pow(Re, .8) * pow(Pr, n)	* pow(dynamic_viscosity / dynamic_viscosity_T, 0.14);
				}
			}break;
			default:
				ESINFO(ERROR) << "Invalid convection variant for EXTERNAL_FORCED.";
			}
	}break;

    default:
	ESINFO(ERROR) << "Invalid convection TYPE.";
	}

	return htc;
}

void NewAdvectionDiffusion2D::processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
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

	DenseMatrix coordinates(e->nodes(), 2), dND(1, 2), q(e->nodes(), 1), htc(e->nodes(), 1), thickness(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

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
			Ke.resize(Ksize, Ksize);
			Ke = 0;
			area = e->regions()[r]->area;
			break;
		}
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

		if (solution.size()) {
			temp = solution[0]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
		} else {
			temp = e->getProperty(Property::INITIAL_TEMPERATURE, n, step.step, 273.15 + 20);
		}
		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, n, step.step, temp) : 0;


		if (solution.size()) {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, 0) - temp);
		} else {
			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, 0));
		}

		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, n, step.step, 0);
		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, 0), 4) - pow(temp, 4));
		q(n, 0) += e->getProperty(Property::HEAT_FLOW, n, step.step, 0) / area;
		q(n, 0) += e->getProperty(Property::HEAT_FLUX, n, step.step, 0);

		emiss(n, 0) *= 4 * temp * temp * temp;

		thickness(n, 0) = e->getProperty(Property::THICKNESS, n, step.step, 1);
		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
			gpHtc.multiply(N[gp], htc);
			gpEmiss.multiply(N[gp], emiss);
			gpThickness.multiply(N[gp], thickness);

			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0) * gpThickness(0, 0), 1, true);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0) * gpThickness(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
		}
	}
}

void NewAdvectionDiffusion2D::processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const
{

}

void NewAdvectionDiffusion2D::processSolution(const Step &step)
{
	_instance->solutions.resize(1, NULL);
	if (_instance->solutions[0] != NULL) {
		delete _instance->solutions[0];
	}

	_instance->solutions[0] = new Solution(*_mesh, "temperature", ElementType::NODES, pointDOFs(), _instance->primalSolution);
}

