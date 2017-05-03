
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
		std::map<std::string, std::string> values;
		for (auto regions = it->second.begin(); regions != it->second.end(); ++regions) {
			values[regions->first] = regions->second->external_temperature;
			_mesh->loadProperty(values, { }, { Property::EXTERNAL_TEMPERATURE });
			values[regions->first] = regions->second->heat_transfer_coefficient;
			_mesh->loadProperty(values, { }, { Property::HEAT_TRANSFER_COEFFICIENT });

			values[regions->first] = regions->second->wall_height;
			_mesh->loadProperty(values, { }, { Property::WALL_HEIGHT });
			values[regions->first] = regions->second->tilt_angle;
			_mesh->loadProperty(values, { }, { Property::TILT_ANGLE });
			values[regions->first] = regions->second->diameter;
			_mesh->loadProperty(values, { }, { Property::DIAMETER });
			values[regions->first] = regions->second->plate_length;
			_mesh->loadProperty(values, { }, { Property::PLATE_LENGTH });
			values[regions->first] = regions->second->fluid_velocity;
			_mesh->loadProperty(values, { }, { Property::FLUID_VELOCITY });
			values[regions->first] = regions->second->plate_distance;
			_mesh->loadProperty(values, { }, { Property::PLATE_DISTANCE });
			values[regions->first] = regions->second->length;
			_mesh->loadProperty(values, { }, { Property::LENGTH });
			values[regions->first] = regions->second->absolute_pressure;
			_mesh->loadProperty(values, { }, { Property::ABSOLUTE_PRESSURE });
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
				BT.multiply(e->dN()[gp], T);
				BTN.multiply(BT, e->N()[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(e->dN()[gp], CDBTN, 1, 1, true);
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

		double T_AVG, g, gas_constant, rho, dynamic_viscosity, heat_capacity, thermal_conductivity;

		T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
		g = 9.81;
		switch (convection.fluid) {
		case espreso::CONVECTION_FLUID::AIR:{

			gas_constant = 286.9;
			rho = (e->getProperty(Property::ABSOLUTE_PRESSURE, node, step, 0)) / (gas_constant * T_AVG);
			heat_capacity = 1002.5 + 275e-6 * pow((T_AVG - 200), 2.0);
			thermal_conductivity = 0.02626 * pow(T_AVG / 300, 0.8646);
			dynamic_viscosity = (1.458e-6 * pow(T_AVG, 1.5)) / (T_AVG + 110.4);
		}break;
		case espreso::CONVECTION_FLUID::WATER:{

			rho = 765.33 + 1.8142 * T_AVG - 0.0035 * pow(T_AVG, 2);
			heat_capacity = 28.07 - 0.2817 * T_AVG + 1.25e-3 * pow(T_AVG, 2.0) - 2.48e-6 * pow(T_AVG, 3.0) + 1.857e-9 * pow(T_AVG, 4.0);
			thermal_conductivity = -0.5752 + 6.397e-3 * T_AVG - 8.151e-6 * pow(T_AVG, 2.0);
			dynamic_viscosity = 9.67e-2 - 8.207e-4 * T_AVG + 2.344e-6 * pow(T_AVG, 2.0) - 2.244e-9 * pow(T_AVG, 3.0);
		}break;
		default:
			ESINFO(ERROR) << "Invalid convection fluid type.";
		}

		switch (convection.variant) {
		case espreso::CONVECTION_VARIANT::INCLINED_WALL: {

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double tilt_angle = e->getProperty(Property::TILT_ANGLE, node, step,	0) * M_PI / 180.0;
			if (RaL <= 1e9) {
				htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * (0.68 + (0.67 * cos(tilt_angle) * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		}break;
		case espreso::CONVECTION_VARIANT::VERTICAL_WALL: {

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (RaL <= 1e9) {
				htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
			} else {
				htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
			}

		}break;
        case espreso::CONVECTION_VARIANT::HORIZONTAL_PLATE_UP:{

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

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

			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

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

			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double Pr = dynamic_viscosity * heat_capacity / thermal_conductivity;

			if ( RaD > 10e12 ){
				// warning!!!!
				ESINFO(ERROR) << "Validated only for RaD <= 10e12 ";
			}

			htc = thermal_conductivity / e->getProperty(Property::DIAMETER, node, step, 0) * pow( 0.6 + ( 0.387*pow(RaD,1.0/6.0)/ pow( 1 + pow( 0.559/Pr, 9.0/16.0), 8.0/27.0) ) ,2.0);

		}break;
		case espreso::CONVECTION_VARIANT::SPHERE:{

			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
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

		double T_AVG, g, gas_constant, rho, dynamic_viscosity, heat_capacity, thermal_conductivity;

		T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
		g = 9.81;
		switch (convection.fluid) {
		case espreso::CONVECTION_FLUID::AIR:{
			gas_constant = 286.9;
			rho = (e->getProperty(Property::ABSOLUTE_PRESSURE, node, step, 0)) / (gas_constant * T_AVG);
			heat_capacity = 1002.5 + 275e-6 * pow((T_AVG - 200), 2.0);
			thermal_conductivity = 0.02626 * pow(T_AVG / 300, 0.8646);
			dynamic_viscosity = (1.458e-6 * pow(T_AVG, 1.5)) / (T_AVG + 110.4);
		}break;
		case espreso::CONVECTION_FLUID::WATER:{
			rho = 765.33 + 1.8142 * T_AVG - 0.0035 * pow(T_AVG, 2);
			heat_capacity = 28.07 - 0.2817 * T_AVG + 1.25e-3 * pow(T_AVG, 2.0) - 2.48e-6 * pow(T_AVG, 3.0) + 1.857e-9 * pow(T_AVG, 4.0);
			thermal_conductivity = -0.5752 + 6.397e-3 * T_AVG - 8.151e-6 * pow(T_AVG, 2.0);
			dynamic_viscosity = 9.67e-2 - 8.207e-4 * T_AVG + 2.344e-6 * pow(T_AVG, 2.0) - 2.244e-9 * pow(T_AVG, 3.0);
		}break;
		default:
			ESINFO(ERROR) << "Invalid convection fluid type.";
		}

		switch (convection.variant) {
		case espreso::CONVECTION_VARIANT::PARALLEL_PLATES: {

			double H_L = e->getProperty(Property::WALL_HEIGHT, node, step, 0) / e->getProperty(Property::LENGTH, node, step, 0);
			double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

			if (( RaL < H_L ) && (temp >  e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  )){

				htc = thermal_conductivity / e->getProperty(Property::WALL_HEIGHT, node, step, 0) * ( 1.0 / 24.0 ) * RaL;

			}else{

				double RaL = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::LENGTH, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);

				if (RaL <= 1e9) {
					htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * (0.68 + (0.67 * pow(RaL,0.25))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),4.0/9.0)) );
				} else {
					htc = (thermal_conductivity	/ e->getProperty(Property::LENGTH, node, step, 0)) * pow(0.825 + (0.387 * pow(RaL,1.0/6.0))/(pow( 1+ pow((0.492 * thermal_conductivity)/(dynamic_viscosity * heat_capacity),9.0/16.0),8.0/27.0)),2 );
				}
			}

		}break;
		case espreso::CONVECTION_VARIANT::CIRCULAR_TUBE: {
			double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
			double H_D = e->getProperty(Property::WALL_HEIGHT, node, step, 0) / e->getProperty(Property::DIAMETER, node, step, 0);

			if ( RaD < H_D ){
				htc = thermal_conductivity / e->getProperty(Property::WALL_HEIGHT, node, step, 0) * ( 1.0 / 128.0 ) * RaD;
			}else{

				double RaD = pow(rho,2)	* g * (1/T_AVG) * heat_capacity * ( temp - e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0)  ) *pow(e->getProperty(Property::DIAMETER, node, step, 0),3.0)/ ( thermal_conductivity * dynamic_viscosity);
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

				double T_AVG, g, gas_constant, rho, dynamic_viscosity, heat_capacity, thermal_conductivity;

				T_AVG = (e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step,	0) + temp) / 2.0;
				g = 9.81;
				switch (convection.fluid) {
				case espreso::CONVECTION_FLUID::AIR:{
					gas_constant = 286.9;
					rho = (e->getProperty(Property::ABSOLUTE_PRESSURE, node, step, 0)) / (gas_constant * T_AVG);
					heat_capacity = 1002.5 + 275e-6 * pow((T_AVG - 200), 2.0);
					thermal_conductivity = 0.02626 * pow(T_AVG / 300, 0.8646);
					dynamic_viscosity = (1.458e-6 * pow(T_AVG, 1.5)) / (T_AVG + 110.4);
				}break;
				case espreso::CONVECTION_FLUID::WATER:{
					rho = 765.33 + 1.8142 * T_AVG - 0.0035 * pow(T_AVG, 2.0);
					heat_capacity = 28.07 - 0.2817 * T_AVG + 1.25e-3 * pow(T_AVG, 2.0) - 2.48e-6 * pow(T_AVG, 3.0) + 1.857e-9 * pow(T_AVG, 4.0);
					thermal_conductivity = -0.5752 + 6.397e-3 * T_AVG - 8.151e-6 * pow(T_AVG, 2.0);
					dynamic_viscosity = 9.67e-2 - 8.207e-4 * T_AVG + 2.344e-6 * pow(T_AVG, 2.0) - 2.244e-9 * pow(T_AVG, 3.0);
				}break;
				default:
					ESINFO(ERROR) << "Invalid convection fluid type.";
				}

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

				double T_EXT, gas_constant, rho, dynamic_viscosity,	dynamic_viscosity_T, heat_capacity, thermal_conductivity;
				switch (convection.fluid) {
				case espreso::CONVECTION_FLUID::AIR:{

					T_EXT = e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step, 0);
					gas_constant = 286.9;
					rho = (e->getProperty(Property::ABSOLUTE_PRESSURE, node, step, 0)) / (gas_constant * T_EXT);
					heat_capacity = 1002.5 + 275e-6 * pow((T_EXT - 200), 2.0);
					thermal_conductivity = 0.02626 * pow(T_EXT / 300, 0.8646);
					dynamic_viscosity = (1.458e-6 * pow(T_EXT, 1.5)) / (T_EXT + 110.4);
					dynamic_viscosity_T = (1.458e-6 * pow(temp, 1.5)) / (temp + 110.4);

				}break;
				case espreso::CONVECTION_FLUID::WATER:{

					T_EXT = e->getProperty(Property::EXTERNAL_TEMPERATURE, node, step, 0);
					rho = 765.33 + 1.8142 * T_EXT - 0.0035 * pow(T_EXT, 2.0);
					heat_capacity = 28.07 - 0.2817 * T_EXT + 1.25e-3 * pow(T_EXT, 2.0) - 2.48e-6 * pow(T_EXT, 3.0) + 1.857e-9 * pow(T_EXT, 4.0);
					thermal_conductivity = -0.5752 + 6.397e-3 * T_EXT - 8.151e-6 * pow(T_EXT, 2.0);
					dynamic_viscosity = 9.67e-2 - 8.207e-4 * T_EXT + 2.344e-6 * pow(T_EXT, 2.0) - 2.244e-9 * pow(T_EXT, 3.0);
					dynamic_viscosity_T = 9.67e-2 - 8.207e-4 * temp + 2.344e-6 * pow(temp, 2.0) - 2.244e-9 * pow(temp, 3.0);

				}break;
				default:
					ESINFO(ERROR) << "Invalid convection fluid type.";
				}

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

	DenseMatrix coordinates(e->nodes(), 2), dND(1, 2), q(e->nodes(), 1), htc(e->nodes(), 1), thickness(e->nodes(), 1), flow(e->nodes(), 1);
	DenseMatrix T(e->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1);

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
	auto stepit = _configuration.convection.find(step.step + 1);
	if (stepit != _configuration.convection.end()) {
		for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
			auto regionit = stepit->second.find(e->regions()[r]->name);
			if (regionit != stepit->second.end()) {
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
			temp = 0;
		}
		T(n, 0) = temp;
		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, n, step.step, temp) : 0;
		q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, n, step.step, 0) - temp);
		q(n, 0) += e->getProperty(Property::HEAT_FLOW, n, step.step, 0) / area;
		q(n, 0) += e->getProperty(Property::HEAT_FLUX, n, step.step, 0);

		thickness(n, 0) = e->getProperty(Property::THICKNESS, n, step.step, 1);
		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, step.step)) {
			gpHtc.multiply(N[gp], htc);
			gpThickness.multiply(N[gp], thickness);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0) * gpThickness(0, 0), 1, true);
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

