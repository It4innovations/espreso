
#include "heattransfer3d.kernel.h"
#include "solverdataprovider/heattransfer.provider.h"
#include "utils/geometry.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/mesh.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "physics/system/builder/builder.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

// TODO: create file with constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

HeatTransfer3DKernel::HeatTransfer3DKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: KernelExecutor(new HeatTransferSolverDataProvider(configuration)),
  gsettings(gsettings),
  iterator(previous ? &previous->iterator : NULL, physics, gsettings, configuration, 3)
{
	geometry::computeBoundaryRegionsArea();

	solutions.push_back(VectorDense(iterator.temperature.output.data->data.size(), iterator.temperature.output.data->data.data()));

	boundaries.reserve(info::mesh->boundaryRegions.size());
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 3);
	}
}

HeatTransfer3DKernel::~HeatTransfer3DKernel()
{
	iterator.clear();
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries[i].clear();
	}
	delete solverDataProvider;
}

void HeatTransfer3DKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Point p(coordinates[0], coordinates[1], coordinates[2]);
	Evaluator::Params params;
	params.coords(3, coordinates);
	params.temp(&temp);

	Point sin, cos;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:

		cos.x = std::cos(d2r(mat->coordinate_system.rotation.x.evaluator->eval(params)));
		cos.y = std::cos(d2r(mat->coordinate_system.rotation.y.evaluator->eval(params)));
		cos.z = std::cos(d2r(mat->coordinate_system.rotation.z.evaluator->eval(params)));

		sin.x = std::sin(d2r(mat->coordinate_system.rotation.x.evaluator->eval(params)));
		sin.y = std::sin(d2r(mat->coordinate_system.rotation.y.evaluator->eval(params)));
		sin.z = std::sin(d2r(mat->coordinate_system.rotation.z.evaluator->eval(params)));

		break;

	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {

		Point origin(
				mat->coordinate_system.center.x.evaluator->eval(params),
				mat->coordinate_system.center.y.evaluator->eval(params),
				0);

		double rotation = std::atan2((p.y - origin.y), (p.x - origin.x));

		cos.x = 1.0;
		cos.y = 1.0;
		cos.z = std::cos(rotation);

		sin.x = 0.0;
		sin.y = 0.0;
		sin.z = std::sin(rotation);

	} break;

	case CoordinateSystemConfiguration::TYPE::SPHERICAL: {

		Point origin(
				mat->coordinate_system.center.x.evaluator->eval(params),
				mat->coordinate_system.center.y.evaluator->eval(params),
				mat->coordinate_system.center.z.evaluator->eval(params));

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

	MatrixDense TCT(3, 3), T(3, 3), C(3, 3), _CD, CT, _CDT, TCDT;

	T(0, 0) = cos.y * cos.z;                         T(0, 1) = cos.y * sin.z;                         T(0, 2) = -sin.y;
	T(1, 0) = cos.z * sin.x * sin.y - cos.x * sin.z; T(1, 1) = cos.x * cos.z + sin.x * sin.y * sin.z; T(1, 2) = cos.y * sin.x;
	T(2, 0) = sin.x * sin.z + cos.x * cos.z * sin.y; T(2, 1) = cos.x * sin.y * sin.z - cos.z * sin.x; T(2, 2) = cos.x * cos.y;

	if (tangentCorrection) {
		_CD.resize(3, 3);
		TCDT.resize(3, 3);
	}

	auto derivation = [&] (const ECFExpression &expression, double h) {
		double temp1 = temp + h, temp2 = temp - h;
		return (
				expression.evaluator->eval(Evaluator::Params().coords(3, coordinates).temp(&temp1)) -
				expression.evaluator->eval(Evaluator::Params().coords(3, coordinates).temp(&temp2))
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = C(2, 2) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = _CD(2, 2) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->eval(params);
		C(0, 1) = C(0, 2) = C(1, 0) = C(1, 2) = C(2, 0) = C(2, 1) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = _CD(0, 2) = _CD(1, 0) = _CD(1, 2) = _CD(2, 0) = _CD(2, 1) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->eval(params);
		C(0, 1) = C(1, 0) = mat->thermal_conductivity.values.get(0, 1).evaluator->eval(params);
		C(0, 2) = C(2, 0) = mat->thermal_conductivity.values.get(0, 2).evaluator->eval(params);
		C(1, 2) = C(2, 1) = mat->thermal_conductivity.values.get(1, 2).evaluator->eval(params);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(0, 2) = _CD(2, 0) = derivation(mat->thermal_conductivity.values.get(0, 2), temp / 1e4);
			_CD(1, 2) = _CD(2, 1) = derivation(mat->thermal_conductivity.values.get(1, 2), temp / 1e4);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(2, 2) = mat->thermal_conductivity.values.get(2, 2).evaluator->eval(params);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->eval(params);
		C(0, 2) = mat->thermal_conductivity.values.get(0, 2).evaluator->eval(params);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->eval(params);
		C(1, 2) = mat->thermal_conductivity.values.get(1, 2).evaluator->eval(params);
		C(2, 0) = mat->thermal_conductivity.values.get(2, 0).evaluator->eval(params);
		C(2, 1) = mat->thermal_conductivity.values.get(2, 1).evaluator->eval(params);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(2, 2) = derivation(mat->thermal_conductivity.values.get(2, 2), temp / 1e4);
			_CD(0, 1) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(0, 2) = derivation(mat->thermal_conductivity.values.get(0, 2), temp / 1e4);
			_CD(1, 2) = derivation(mat->thermal_conductivity.values.get(1, 2), temp / 1e4);
			_CD(1, 0) = derivation(mat->thermal_conductivity.values.get(1, 0), temp / 1e4);
			_CD(2, 0) = derivation(mat->thermal_conductivity.values.get(2, 0), temp / 1e4);
			_CD(2, 1) = derivation(mat->thermal_conductivity.values.get(2, 1), temp / 1e4);
		}
		break;
	default:
		eslog::error("Advection diffusion 3D not supports set material model.\n");
	}

	CT.multiply(C, T);
	TCT.multiply(T, CT, 1, 0, true, false);
	if (tangentCorrection) {
		_CDT.multiply(_CD, T);
		TCDT.multiply(T, _CDT, 1, 0, true, false);
		CD(node, 0) += phase * TCDT(0, 0);
		CD(node, 1) += phase * TCDT(1, 1);
		CD(node, 2) += phase * TCDT(2, 2);
		CD(node, 3) += phase * TCDT(0, 1);
		CD(node, 4) += phase * TCDT(0, 2);
		CD(node, 5) += phase * TCDT(1, 0);
		CD(node, 6) += phase * TCDT(1, 2);
		CD(node, 7) += phase * TCDT(2, 0);
		CD(node, 8) += phase * TCDT(2, 1);
	}

	K(node, 0) += phase * TCT(0, 0);
	K(node, 1) += phase * TCT(1, 1);
	K(node, 2) += phase * TCT(2, 2);
	K(node, 3) += phase * TCT(0, 1);
	K(node, 4) += phase * TCT(0, 2);
	K(node, 5) += phase * TCT(1, 0);
	K(node, 6) += phase * TCT(1, 2);
	K(node, 7) += phase * TCT(2, 0);
	K(node, 8) += phase * TCT(2, 1);
}

void HeatTransfer3DKernel::processElement(const Builder &builder, const HeatTransferElementIterator &iterator, InstanceFiller &filler) const
{
	filler.DOFs = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	bool CAU = gsettings.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = builder.tangentMatrixCorrection;

	MatrixDense Ce(3, 3), coordinates(filler.DOFs, 3), J(3, 3), invJ(3, 3), dND, CedND;
	double detJ, tauK, xi = 1, C1 = 1, C2 = 6;
	MatrixDense f(filler.DOFs, 1);
	MatrixDense U(filler.DOFs, 3);
	MatrixDense m(filler.DOFs, 1);
	MatrixDense T(filler.DOFs, 1);
	MatrixDense K(filler.DOFs, 9);
	MatrixDense gpK(1, 9), gpM(1, 1), gpF(1, 1), gpDF(1, 1);
	std::vector<MatrixDense> gpT(N.size());
	MatrixDense tangentK, BT, BTN, gpCD, CD, CDBTN, CDe, bhK(filler.DOFs, filler.DOFs);
	MatrixDense gKe(filler.DOFs, filler.DOFs);

	if (tangentCorrection) {
		CD.resize(filler.DOFs, 9); CD.fill(0);
		CDe.resize(3, 3); CDe.fill(0);
	}

	const MaterialBaseConfiguration *phase1 = NULL, *phase2 = NULL;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < filler.DOFs; n++) {
		Evaluator::Params params;
		params.coords(3, iterator.coordinates.data);
		params.temp(iterator.temperature.data);
		T(n, 0) = iterator.temperature.data[n];
		coordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, phase1, phase, step::time::current, T(n, 0), K, CD, tangentCorrection);
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, phase2, (1 - phase), step::time::current, T(n, 0), K, CD, tangentCorrection);
			double dens1 = phase1->density.evaluator->eval(params);
			double dens2 = phase2->density.evaluator->eval(params);
			double hc1 = phase1->heat_capacity.evaluator->eval(params);
			double hc2 = phase2->heat_capacity.evaluator->eval(params);

			m(n, 0) = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation);
		} else {
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, iterator.material, 1, step::time::current, T(n, 0), K, CD, tangentCorrection);
			double dens = iterator.material->density.evaluator->eval(params);
			double hc = iterator.material->heat_capacity.evaluator->eval(params);
			m(n, 0) = dens * hc;
		}

		U(n, 0) = iterator.motion.data[3 * n + 0] * m(n, 0);
		U(n, 1) = iterator.motion.data[3 * n + 1] * m(n, 0);
		U(n, 2) = iterator.motion.data[3 * n + 2] * m(n, 0);
		f(n, 0) = iterator.heatSource.data[n];
	}

	if ((builder.matrices & Builder::Request::K) || ((builder.matrices & Builder::Request::R) && builder.timeIntegrationConstantK != 0)) {
		filler.Ke.resize(filler.DOFs, filler.DOFs);
		filler.Ke.fill(0);
	}

	if ((builder.matrices & Builder::Request::M) || ((builder.matrices & Builder::Request::R) && builder.timeIntegrationConstantM != 0)) {
		filler.Me.resize(filler.DOFs, filler.DOFs);
		filler.Me.fill(0);
	}
	if (builder.matrices & Builder::Request::R) {
		filler.Re.resize(filler.DOFs);
		filler.Re.fill(0);
	}
	if (builder.matrices & Builder::Request::f) {
		filler.Fe.resize(filler.DOFs);
		filler.Fe.fill(0);
	}

	if (tangentCorrection) {
		tangentK.resize(filler.DOFs, filler.DOFs);
		tangentK.fill(0);
	}

	MatrixDense g(1, 3), u(1, 3), v(1, 3), re(1, filler.DOFs);
	double normGradN = 0;

	if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split) {
		g(0, 0) = iterator.gradient.data[0];
		g(0, 1) = iterator.gradient.data[1];
		g(0, 2) = iterator.gradient.data[2];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		gpT[gp].resize(1, 1);
		gpT[gp].multiply(N[gp], T);
	}

	Evaluator::Params params;
	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = MATH::determinant3x3(J.vals);
		if (detJ < 0) { ++filler.invalid; detJ = -detJ; }
		MATH::Dense3x3inverse(J.vals, invJ.vals, detJ);

		gpK.multiply(N[gp], K);
		if (tangentCorrection) {
			gpCD.multiply(N[gp], CD);
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
		gpM.multiply(N[gp], m);
		gpF.multiply(N[gp], f);
		gpDF(0, 0) = 0;

//		 TODO: optimize
//		if (gsettings.load_steps_settings.at(step::loadstep + 1).bio_heat.size()) {
//			for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
//				int maskSize = info::mesh->elements->regionMaskSize;
//				esint maskOffset = r / (8 * sizeof(esint));
//				esint bit = 1 << (r % (8 * sizeof(esint)));
//				if (info::mesh->elements->regions->datatarray()[iterator.offset * maskSize + maskOffset] & bit) {
//					auto bioheat = gsettings.load_steps_settings.at(step::loadstep + 1).bio_heat.find(info::mesh->elementsRegions[r]->name);
//					if (bioheat != gsettings.load_steps_settings.at(step::loadstep + 1).bio_heat.end()) {
//						params.coords(3, iterator.coordinates.data);
//						params.inittemp(iterator.initialTemperature.data);
//						params.temp(gpT[gp].vals);
//						double arte = bioheat->second.arteriar_blood_temperature.evaluator->eval(params);
//						double heat = bioheat->second.blood_specific_heat.evaluator->eval(params);
//						double dens = bioheat->second.blood_density.evaluator->eval(params);
//						double meta = bioheat->second.metabolic_heat_source.evaluator->eval(params);
//						double perf = bioheat->second.blood_perfusion.evaluator->eval(params);
//						double ref =  bioheat->second.reference_temperature.evaluator->eval(params);
//						double awm =  bioheat->second.physical_activity_scatter_factor.evaluator->eval(params);
//						double mu =   bioheat->second.mu.evaluator->eval(params);
//						if (std::fabs(mu) < 1e-15) {
//							gpF(0, 0) += perf * dens * heat  * (arte - gpT[gp](0, 0)) + meta;
//							gpDF(0, 0) -= perf * dens * heat;
//						} else {
//							double act = 5.7;
//							double ny = 0.0;
//
//							if (act < 1.6){
//								ny = 0.0;
//							} else if (act >= 1.6 && act <= 5) {
//								ny = 0.2 * (1 - tanh(0.39*act-0.6));
//							} else {
//								ny = 0.2;
//							}
//
//							gpF(0, 0) += perf * (dens * heat + mu * (awm*((87.1*act*(ny - 1.0))/0.8 - 87.1 ) + pow(2, (gpT[gp](0, 0) - ref) / 10) - 1)) * (arte - gpT[gp](0, 0))  + awm*((87.1*act*(ny - 1.0))/0.8 - 87.1 ) + meta + perf * (pow(2, (gpT[gp](0, 0) -ref) / 10) - 1);
//							gpDF(0, 0) += (pow(2,(gpT[gp](0, 0)/10 - ref/10))*perf*log(2))/10 - perf*(mu*(pow(2,(gpT[gp](0, 0)/10 - ref/10)) + awm*((871*act*(ny - 1))/8 - 871/10) - 1) + dens*heat) - (pow(2,(gpT[gp](0, 0)/10 - ref/10))*mu*perf*log(2)*(gpT[gp](0, 0) - arte))/10;
//							//	gpF(0, 0) += perf * (dens * heat + mu * (pow(2, (gpT[gp](0, 0) - ref) / 10) - 1)) * (arte - gpT[gp](0, 0)) + meta + perf * (pow(2, (gpT[gp](0, 0) -ref) / 10) - 1);
//							//	gpDF(0, 0) += perf * (pow(2, (gpT[gp](0, 0) - ref) / 10) * mu * log(2) * (arte - gpT[gp](0, 0)) / 10 - (dens * heat + mu * (pow(2, (gpT[gp](0, 0) - ref) / 10) - 1))) + perf * pow(2, (gpT[gp](0, 0) - ref) / 10) * log(2) / 10;
//						}
//					}
//				}
//			}
//		}


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

		MatrixDense b_e(1, filler.DOFs), b_e_c(1, filler.DOFs), g_e(1, filler.DOFs);
		b_e.multiply(u, dND, 1, 0);
		g_e.multiply(g, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (esint i = 0; i < re.ncols; i++) {
					re(0, i) = b_e(0, i) - f(0, i);
				}
				MatrixDense ReBt(1, 3);
				ReBt.multiply(re, dND, 1 / pow(normGradN, 2), 0, false, true);
				for (esint i = 0; i < ReBt.ncols; i++) {
					v(0, i) = u(0, i) - ReBt(0, i);
				}
			} else {
				v = u;
			}
		}

		double norm_u_e = u.norm();
		double h_e = 0, tau_e = 0, konst = 0, gh_e = 0;
		double C_e = 0;

		if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split && g.norm() != 0) {
			gh_e = 2 * g.norm() / g_e.norm();
			tauK = (C1 * gh_e * gh_e) / (Ce(0, 0) * C2 + gh_e * gh_e * (gpM(0, 0) / step::time::shift));
			xi = std::max(1., 1 / (1 - tauK * gpM(0, 0) / step::time::shift));
		}

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);

			if (CAU) {
				MatrixDense u_v(1, 3);
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

		Ce(0, 0) += gsettings.sigma * h_e * norm_u_e;
		Ce(1, 1) += gsettings.sigma * h_e * norm_u_e;
		Ce(2, 2) += gsettings.sigma * h_e * norm_u_e;

		if (builder.matrices & (Builder::Request::M | Builder::Request::R)) {
			filler.Me.multiply(N[gp], N[gp], detJ * gpM(0, 0) * weighFactor[gp], 1, true);
		}
		if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
			if (tangentCorrection) {
				BT.multiply(dND, T);
				BTN.multiply(BT, N[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * weighFactor[gp], 1, true);
			}
			bhK.multiply(N[gp], N[gp], detJ * gpDF(0, 0) * weighFactor[gp], 1, true);
			if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split) {
				CedND.multiply(Ce, dND);
				gKe.multiply(dND, CedND, detJ * weighFactor[gp], 1, true);
				gKe.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
				if (konst * weighFactor[gp] * detJ != 0) {
					gKe.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
				}
				if (CAU) {
					gKe.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
				}
			}
			CedND.multiply(Ce, dND);
			filler.Ke.multiply(dND, CedND, xi * detJ * weighFactor[gp], 1, true);
			filler.Ke.multiply(N[gp], b_e, xi * detJ * weighFactor[gp], 1, true);
			if (konst * weighFactor[gp] * detJ != 0) {
				filler.Ke.multiply(b_e, b_e, xi * konst * weighFactor[gp] * detJ, 1, true);
			}
			if (CAU) {
				filler.Ke.multiply(dND, dND, xi * C_e * weighFactor[gp] * detJ, 1, true);
			}
		}

		if (builder.matrices & Builder::Request::f) {
			for (esint i = 0; i < filler.DOFs; i++) {
				filler.Fe[0][i] += detJ * weighFactor[gp] * gpF(0, 0) * N[gp](0, i);
				if (norm_u_e != 0) {
					filler.Fe[0][i] += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * gpF(0, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split) {
		MatrixDense T1, T2;
		T1.multiply(filler.Ke, T, 1, 0);
		T2.multiply(gKe, T, 1, 0);
		for (esint i = 0; i < filler.DOFs; i++) {
			filler.Fe[0][i] += T1(i, 0) - T2(i, 0);
		}
	}

	if (builder.matrices & Builder::Request::R) {
		filler.Re.multiply(filler.Ke, T, builder.timeIntegrationConstantK, 0);
		if (builder.matrices & Builder::Request::M) {
			filler.Re.multiply(filler.Me, T, builder.timeIntegrationConstantM, 1);
		}
	}

	if (tangentCorrection) {
		filler.Ke.add(1, &tangentK);
	}
	filler.Ke.add(-1, &bhK);

	filler.insertK = builder.matrices & Builder::Request::K;
	filler.insertM = builder.matrices & Builder::Request::M;
	filler.insertC = false;
	filler.insertR = builder.matrices & Builder::Request::R;
	filler.insertF = builder.matrices & Builder::Request::f;
}

void HeatTransfer3DKernel::processFace(const Builder &builder, const HeatTransferBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = iterator.element->nodes;
	if (!iterator.convection && !iterator.heatflow.data && !iterator.heatflux.data && !iterator.radiation) {
		return;
	}
	if (!(builder.matrices & (Builder::Request::K | Builder::Request::f))) {
		return;
	}

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense coordinates(filler.DOFs, 3), dND(1, 3), q(filler.DOFs, 1), htc(filler.DOFs, 1), flow(filler.DOFs, 1), emiss(filler.DOFs, 1);
	MatrixDense gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);

	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(filler.DOFs);
		filler.Fe.fill(0);
	}
	if ((filler.insertK = (iterator.convection || iterator.radiation))) {
		filler.Ke.resize(filler.DOFs, filler.DOFs);
		filler.Ke.fill(0);
	}

	for (esint n = 0; n < filler.DOFs; n++) {
		double temp = iterator.temperature.data[n];
		coordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		if (iterator.convection) {
			double text = iterator.extemperature.data[n];
			htc(n, 0) = iterator.htc.data[n];

			if (step::iteration) {
				q(n, 0) += htc(n, 0) * (text - temp);
			} else {
				q(n, 0) += htc(n, 0) * (text);
			}
		}

		if (iterator.radiation) {
			emiss(n, 0) = CONST_Stefan_Boltzmann * iterator.emissivity.data[n];
			q(n, 0) += emiss(n, 0) * (pow(iterator.extemperature.data[n], 4) - pow(temp, 4));
			emiss(n, 0) *= 4 * temp * temp * temp;
		}
		if (iterator.heatflow.data) {
			q(n, 0) += iterator.heatflow.data[n] / iterator.regionArea;
		}
		if (iterator.heatflux.data) {
			q(n, 0) += iterator.heatflux.data[n];
		}
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();

		gpQ.multiply(N[gp], q);
		if (iterator.convection) {
			gpHtc.multiply(N[gp], htc);
			filler.Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		if (iterator.radiation) {
			gpEmiss.multiply(N[gp], emiss);
			filler.Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
		}
		for (esint i = 0; i < filler.DOFs; i++) {
			filler.Fe[0][i] += J * weighFactor[gp] * N[gp](0, i % filler.DOFs) * gpQ(0, 0);
		}
	}
}

void HeatTransfer3DKernel::processEdge(const Builder &builder, const HeatTransferBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = iterator.element->nodes;
//	if (!(e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step) ||
//		e->hasProperty(Property::HEAT_FLOW, _step->step) ||
//		e->hasProperty(Property::HEAT_FLUX, _step->step))) {
//
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(0, 0);
//		return;
//	}
//	if (!(builder.matrices & (Builder::Request::K | Builder::Request::f))) {
//		Ke.resize(0, 0);
//		Me.resize(0, 0);
//		Re.resize(0, 0);
//		fe.resize(0, 0);
//		return;
//	}
//
//	MatrixDense coordinates(e->nodes(), 3), dND(1, 3), q(e->nodes(), 1), htc(e->nodes(), 1), flow(e->nodes(), 1), emiss(e->nodes(), 1);
//	MatrixDense gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1), gpEmiss(1, 1);
//
//	double area = 1, temp;
//	esint Ksize = e->nodes();
//	Ke.resize(0, 0);
//	Me.resize(0, 0);
//	Re.resize(0, 0);
//	fe.resize(0, 0);
//
//	if (builder.matrices & Builder::Request::f) {
//		fe.resize(Ksize, 1);
//		fe = 0;
//	}
//
//	for (size_t r = 0; r < e->regions().size(); r++) {
//		if (_step->step < e->regions()[r]->settings.size() && e->regions()[r]->settings[_step->step].count(Property::HEAT_FLOW)) {
//			area = e->regions()[r]->area;
//			break;
//		}
//	}
//	if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step)) {
//		Ke.resize(Ksize, Ksize);
//		Ke = 0;
//	}
//
//	const std::vector<MatrixDense> &dN = e->dN();
//	const std::vector<MatrixDense> &N = e->N();
//	const std::vector<double> &weighFactor = e->weighFactor();
//
//	const ConvectionConfiguration *convection = NULL;
//	for (size_t r = 0; convection == NULL && r < e->regions().size(); r++) {
//		auto regionit = _configuration.load_stepsgsettings.at(_step->step + 1).convection.find(e->regions()[r]->name);
//		if (regionit != _configuration.load_stepsgsettings.at(_step->step + 1).convection.end()) {
//			convection = &regionit->second;
//		}
//	}

//	for (size_t n = 0; n < e->nodes(); n++) {
//		coordinates(n, 0) = _mesh->coordinates()[e->node(n)].x;
//		coordinates(n, 1) = _mesh->coordinates()[e->node(n)].y;
//		coordinates(n, 2) = _mesh->coordinates()[e->node(n)].z;
//
//		temp = solution[offset + SolutionIndex::TEMPERATURE]->get(Property::TEMPERATURE, e->domains().front(), _mesh->coordinates().localIndex(e->node(n), e->domains().front()));
//		htc(n, 0) = convection != NULL ? computeHTC(*convection, e, _mesh->coordinates()[e->node(n)], step, temp) : 0;
//
//		if (_step->iteration) {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0) - temp);
//		} else {
//			q(n, 0) += htc(n, 0) * (e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0));
//		}
//
//		emiss(n, 0) = CONST_Stefan_Boltzmann * e->getProperty(Property::EMISSIVITY, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0);
//		q(n, 0) += emiss(n, 0) * (pow(e->getProperty(Property::EXTERNAL_TEMPERATURE, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0), 4) - pow(temp, 4));
//		q(n, 0) += e->getProperty(Property::HEAT_FLOW, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0) / area;
//		q(n, 0) += e->getProperty(Property::HEAT_FLUX, _step->step, _mesh->coordinates()[e->node(n)], _step->time::currentTime, temp, 0);
//
//		emiss(n, 0) *= 4 * temp * temp * temp;
//	}
//
//	for (size_t gp = 0; gp < e->gaussePoints(); gp++) {
//		dND.multiply(dN[gp], coordinates);
//		double J = dND.norm();
//		gpQ.multiply(N[gp], q);
//		if (e->hasProperty(Property::EXTERNAL_TEMPERATURE, _step->step)) {
//			gpHtc.multiply(N[gp], htc);
//			gpEmiss.multiply(N[gp], emiss);
//
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
//			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpEmiss(0, 0), 1, true);
//		}
//		for (esint i = 0; i < Ksize; i++) {
//			fe(i, 0) += J * weighFactor[gp] * N[gp](0, i % e->nodes()) * gpQ(0, 0);
//		}
//	}
}

void HeatTransfer3DKernel::elementSolution(const HeatTransferElementIterator &iterator)
{
	esint size = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);

	MatrixDense Ce(3, 3), coordinates(size, 3), J(3, 3), invJ(3, 3), dND, T(size, 1), dNDT;
	double detJ, norm_u_e, h_e;
	MatrixDense U(size, 3), K(size, 9), gpK(1, 9), CD;
	MatrixDense u(1, 3), matFlux(3, 1), matGradient(3, 1);

	const MaterialBaseConfiguration *phase1 = NULL, *phase2 = NULL;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
		iterator.phase.data[0] = 0;
		iterator.latentHeat.data[0] = 0;
	}

	for (int n = 0; n < size; n++) {
		T(n, 0) = iterator.temperature.data[n];
		coordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, phase1, phase, step::time::current, T(n, 0), K, CD, false);
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, phase2, (1 - phase), step::time::current, T(n, 0), K, CD, false);
			iterator.phase.data[0] += phase;
			iterator.latentHeat.data[0] += iterator.material->latent_heat * derivation;
		} else {
			assembleMaterialMatrix(n, iterator.coordinates.data + 3 * n, iterator.material, 1, step::time::current, T(n, 0), K, CD, false);
		}

		U(n, 0) = iterator.motion.data[3 * n + 0];
		U(n, 1) = iterator.motion.data[3 * n + 1];
		U(n, 2) = iterator.motion.data[3 * n + 2];
	}
	if (iterator.material->phase_change) {
		iterator.phase.data[0] /= size;
		iterator.latentHeat.data[0] /= size;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = MATH::determinant3x3(J.vals);
		if (detJ < 0) { detJ = -detJ; }
		MATH::Dense3x3inverse(J.vals, invJ.vals, detJ);

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

		norm_u_e = u.norm();
		h_e = 0;

		if (norm_u_e != 0) {
			MatrixDense b_e(1, size);
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += gsettings.sigma * h_e * norm_u_e;
		Ce(1, 1) += gsettings.sigma * h_e * norm_u_e;
		Ce(2, 2) += gsettings.sigma * h_e * norm_u_e;

		if (info::ecf->output.results_selection.gradient) {
			matGradient.multiply(dND, T, 1, 1);
		}
		if (info::ecf->output.results_selection.flux) {
			dNDT.multiply(dND, T);
			matFlux.multiply(Ce, dNDT, 1, 1);
		}
	}

	if (info::ecf->output.results_selection.gradient) {
		iterator.gradient.data[0] = matGradient(0, 0) / N.size();
		iterator.gradient.data[1] = matGradient(1, 0) / N.size();
		iterator.gradient.data[2] = matGradient(2, 0) / N.size();
	}

	if (info::ecf->output.results_selection.flux) {
		iterator.flux.data[0] = matFlux(0, 0) / N.size();
		iterator.flux.data[1] = matFlux(1, 0) / N.size();
		iterator.flux.data[2] = matFlux(2, 0) / N.size();
	}
}

