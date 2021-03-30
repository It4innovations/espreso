
#include "heattransfer2d.kernel.h"
#include "solverdataprovider/heattransfer.provider.h"
#include "utils/geometry.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "physics/system/builder/builder.h"
#include "mesh/mesh.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

// TODO: create file with constants
#define CONST_Stefan_Boltzmann 5.6703e-8

using namespace espreso;

HeatTransfer2DKernel::HeatTransfer2DKernel(HeatTransfer2DKernel *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration)
: KernelExecutor(new HeatTransferSolverDataProvider(configuration)),
  gsettings(gsettings),
  iterator(previous ? &previous->iterator : NULL, physics, gsettings, configuration, 2)
{
	geometry::computeBoundaryRegionsArea();

	solutions.push_back(VectorDense(iterator.temperature.output.data->data.size(), iterator.temperature.output.data->data.data()));

	boundaries.reserve(info::mesh->boundaryRegions.size());
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 2);
	}
}

HeatTransfer2DKernel::~HeatTransfer2DKernel()
{
	iterator.clear();
	for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
		boundaries[i].clear();
	}
	delete solverDataProvider;
}

void boundaryWithSettings(size_t rindex)
{

}

void HeatTransfer2DKernel::assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const
{
	auto d2r = [] (double degree) -> double {
		return M_PI * degree / 180;
	};

	Evaluator::Params params;
	params.coords(2, coordinates);
	params.temp(&temp);

	double cos = 0, sin = 0;
	switch (mat->coordinate_system.type) {
	case CoordinateSystemConfiguration::TYPE::CARTESIAN:
		cos = std::cos(d2r(mat->coordinate_system.rotation.z.evaluator->eval(params)));
		sin = std::sin(d2r(mat->coordinate_system.rotation.z.evaluator->eval(params)));
		break;
	case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: {
		Point origin(
				mat->coordinate_system.center.x.evaluator->eval(params),
				mat->coordinate_system.center.y.evaluator->eval(params),
				0);
		double rotation = std::atan2((coordinates[1] - origin.y), (coordinates[0] - origin.x));
		cos = std::cos(rotation);
		sin = std::sin(rotation);
		break;
	}
	default:
		eslog::error("Invalid material type (SPHERICAL for 2D).\n");
	}

	MatrixDense TCT(2, 2), T(2, 2), C(2, 2), CT, _CD, _CDT, TCDT;
	T(0, 0) =  cos; T(0, 1) = sin;
	T(1, 0) = -sin; T(1, 1) = cos;

	if (tangentCorrection) {
		_CD.resize(2, 2);
		TCDT.resize(2, 2);
	}

	auto derivation = [&] (const ECFExpression &expression, double h) {
		double temp1 = temp + h, temp2 = temp - h;
		return (
				expression.evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp1)) -
				expression.evaluator->eval(Evaluator::Params().coords(2, coordinates).temp(&temp2))
				) / (2 * h);
	};

	switch (mat->thermal_conductivity.model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		C(0, 0) = C(1, 1) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = _CD(1, 1) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(0, 1) = C(1, 0) = 0;
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = 0;
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(1, 0) = C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->eval(params);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = _CD(1, 0) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		C(0, 0) = mat->thermal_conductivity.values.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->thermal_conductivity.values.get(1, 1).evaluator->eval(params);
		C(0, 1) = mat->thermal_conductivity.values.get(0, 1).evaluator->eval(params);
		C(1, 0) = mat->thermal_conductivity.values.get(1, 0).evaluator->eval(params);
		if (tangentCorrection) {
			_CD(0, 0) = derivation(mat->thermal_conductivity.values.get(0, 0), temp / 1e4);
			_CD(1, 1) = derivation(mat->thermal_conductivity.values.get(1, 1), temp / 1e4);
			_CD(0, 1) = derivation(mat->thermal_conductivity.values.get(0, 1), temp / 1e4);
			_CD(1, 0) = derivation(mat->thermal_conductivity.values.get(1, 0), temp / 1e4);
		}
		break;
	default:
		eslog::error("Advection diffusion 2D not supports set material model.\n");
	}

	if (tangentCorrection) {
		_CDT.multiply(_CD, T);
		TCDT.multiply(T, _CDT, 1, 0, true, false);
		CD(node, 0) += phase * TCDT(0, 0);
		CD(node, 1) += phase * TCDT(1, 1);
		CD(node, 2) += phase * TCDT(0, 1);
		CD(node, 3) += phase * TCDT(1, 0);
	}

	CT.multiply(C, T);
	TCT.multiply(T, CT, 1, 0, true, false);
	K(node, 0) += phase * TCT(0, 0);
	K(node, 1) += phase * TCT(1, 1);
	K(node, 2) += phase * TCT(0, 1);
	K(node, 3) += phase * TCT(1, 0);
}

void HeatTransfer2DKernel::processElement(const Builder &builder, HeatTransferElementIterator &iterator, InstanceFiller &filler) const
{
	filler.DOFs = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	bool CAU = gsettings.stabilization == HeatTransferConfiguration::STABILIZATION::CAU;
	bool tangentCorrection = builder.tangentMatrixCorrection;

	MatrixDense Ce(2, 2), J(2, 2), invJ(2, 2), dND, CedND;
	double detJ, tauK, xi = 1, C1 = 1, C2 = 6;
	MatrixDense f(filler.DOFs, 1);
	MatrixDense U(filler.DOFs, 2);
	MatrixDense m(filler.DOFs, 1);
	MatrixDense thickness(filler.DOFs, 1), K(filler.DOFs, 4);
	MatrixDense gpThickness(1, 1), gpK(1, 4), gpM(1, 1);
	MatrixDense tangentK, BT, BTN, gpCD, CD, CDBTN, CDe;
	MatrixDense gKe(filler.DOFs, filler.DOFs);

	if (tangentCorrection) {
		CD.resize(filler.DOFs, 4);
		CDe.resize(2, 2);
	}

	const MaterialBaseConfiguration *phase1 = NULL, *phase2 = NULL;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
	}

	for (int n = 0; n < filler.DOFs; n++) {
		Evaluator::Params params;
		params.coords(2, iterator.coordinates.data + 2 * n);
		params.temp(iterator.temperature.data + n);
		thickness(n, 0) = iterator.thickness.data[n];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, iterator.temperature.data[n], iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, phase1, phase, step::time.current, iterator.temperature.data[n], K, CD, tangentCorrection);
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, phase2, (1 - phase), step::time.current, iterator.temperature.data[n], K, CD, tangentCorrection);
			double dens1 = phase1->density.evaluator->eval(params);
			double dens2 = phase2->density.evaluator->eval(params);
			double hc1 = phase1->heat_capacity.evaluator->eval(params);
			double hc2 = phase2->heat_capacity.evaluator->eval(params);

			m(n, 0) = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation) * iterator.thickness.data[0];
		} else {
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, iterator.material, 1, step::time.current, iterator.temperature.data[n], K, CD, tangentCorrection);
			double dens = iterator.material->density.evaluator->eval(params);
			double hc = iterator.material->heat_capacity.evaluator->eval(params);
			m(n, 0) = dens * hc * thickness(n, 0);
		}
		U(n, 0) = iterator.motion.data[2 * n + 0] * m(n, 0);
		U(n, 1) = iterator.motion.data[2 * n + 1] * m(n, 0);
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
	}

	MatrixDense g(1, 2), u(1, 2), v(1, 2), re(1, filler.DOFs);
	double normGradN = 0;

	if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split) {
		g(0, 0) = iterator.gradient.data[0];
		g(0, 1) = iterator.gradient.data[1];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], filler.DOFs, 2, iterator.coordinates.data);
		detJ = MATH::determinant2x2(J.vals);
		if (detJ < 0) { ++filler.invalid; detJ = -detJ; }
		MATH::Dense2x2inverse(J.vals, invJ.vals, detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);
		if (tangentCorrection) {
			gpCD.multiply(N[gp], CD);
			CDe(0, 0) = gpCD(0, 0);
			CDe(1, 1) = gpCD(0, 1);
			CDe(0, 1) = gpCD(0, 2);
			CDe(1, 0) = gpCD(0, 3);
		}
		gpM.multiply(N[gp], m);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		MatrixDense b_e(1, filler.DOFs), b_e_c(1, filler.DOFs), g_e(1, filler.DOFs);
		b_e.multiply(u, dND, 1, 0);
		g_e.multiply(g, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (esint i = 0; i < re.ncols; i++) {
					re(0, i) = b_e(0, i) - f(i, 0);
				}
				MatrixDense ReBt(1, 2);
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

		if (gsettings.diffusion_split && g.norm() != 0) {
			gh_e = 2 * g.norm() / g_e.norm();
			tauK = (C1 * gh_e * gh_e) / (Ce(0, 0) * C2 + gh_e * gh_e * (gpM(0, 0) / step::time.shift));
			xi = std::max(1., 1 / (1 - tauK * gpM(0, 0) / step::time.shift));
		}

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);

			if (CAU) {
				MatrixDense u_v(1, 2);
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

		Ce(0, 0) += gsettings.sigma * h_e * norm_u_e;
		Ce(1, 1) += gsettings.sigma * h_e * norm_u_e;

		if (builder.matrices & (Builder::Request::M | Builder::Request::R)) {
			filler.Me.multiply(N[gp], N[gp], detJ * gpM(0, 0) * weighFactor[gp], 1, true);
		}
		if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
			if (tangentCorrection) {
				BT.multiply(dND, filler.DOFs, 1, iterator.temperature.data);
				BTN.multiply(BT, N[gp]);
				CDBTN.multiply(CDe, BTN);
				tangentK.multiply(dND, CDBTN,  detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
			}
			if (gsettings.diffusion_split) {
				CedND.multiply(Ce, dND);
				gKe.multiply(dND, CedND, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
				gKe.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
				if (konst * weighFactor[gp] * detJ != 0) {
					gKe.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
				}
				if (CAU) {
					gKe.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
				}
			}
			CedND.multiply(Ce, dND);
			filler.Ke.multiply(dND, CedND, xi * detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
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
				filler.Fe[0][i] += detJ * weighFactor[gp] * N[gp](0, i) * f(i, 0);
				if (norm_u_e != 0) {
					filler.Fe[0][i] += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(i, 0) / (2 * norm_u_e);
				}
			}
		}
	}

	if ((builder.matrices & Builder::Request::M) && gsettings.diffusion_split) {
		MatrixDense T1, T2;
		T1.multiply(filler.Ke, filler.DOFs, 1, iterator.temperature.data, 1, 0);
		T2.multiply(gKe, filler.DOFs, 1, iterator.temperature.data, 1, 0);
		for (esint i = 0; i < filler.DOFs; i++) {
			f(i, 0) += T1(i, 0) - T2(i, 0);
		}
	}

	if (builder.matrices & Builder::Request::R) {
		filler.Re.multiply(filler.Ke, filler.DOFs, 1, iterator.temperature.data, builder.timeIntegrationConstantK, 0);
		if (builder.matrices & Builder::Request::M) {
			filler.Re.multiply(filler.Me, filler.DOFs, 1, iterator.temperature.data, builder.timeIntegrationConstantM, 1);
		}
	}

	if (tangentCorrection) {
		filler.Ke.add(1, &tangentK);
	}

	filler.insertK = builder.matrices & Builder::Request::K;
	filler.insertM = builder.matrices & Builder::Request::M;
	filler.insertC = false;
	filler.insertR = builder.matrices & Builder::Request::R;
	filler.insertF = builder.matrices & Builder::Request::f;
}

void HeatTransfer2DKernel::processEdge(const Builder &builder, HeatTransferBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = iterator.element->nodes;
	if (!iterator.convection && !iterator.heatflow.data && !iterator.heatflux.data && !iterator.radiation) {
		return;
	}
	if (!(builder.matrices & (Builder::Request::K | Builder::Request::f))) {
		return;
	}

	filler.DOFs = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense coordinates(filler.DOFs, 2), dND(1, 2), q(filler.DOFs, 1), htc(filler.DOFs, 1), thickness(filler.DOFs, 1), emiss(filler.DOFs, 1);
	MatrixDense gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpEmiss(1, 1);

	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(filler.DOFs);
		filler.Fe.fill(0);
	}

	if ((filler.insertK = (iterator.convection || iterator.radiation))) {
		filler.Ke.resize(filler.DOFs, filler.DOFs);
		filler.Ke.fill(0);
	}

	for (int n = 0; n < filler.DOFs; n++) {
		double temp = iterator.temperature.data[n];
		coordinates(n, 0) = iterator.coordinates.data[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[2 * n + 1];
		thickness(n, 0) = 1; // FIXME: iterator.thickness.data[n];

		if (iterator.convection) {
			double text = iterator.extemperature.data[n];
			htc(n, 0) = iterator.htc.data[n];

			if (step::outstep.iteration) {
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

		q(n, 0) *= thickness(n, 0);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);

		if (iterator.convection || iterator.radiation) {
			gpThickness.multiply(N[gp], thickness);
		}
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

void HeatTransfer2DKernel::elementSolution(HeatTransferElementIterator &iterator)
{
	int DOFs = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);

	MatrixDense Ce(2, 2), coordinates(DOFs, 2), J(2, 2), invJ(2, 2), dND, dNDT, T(DOFs, 1);
	double detJ, m, norm_u_e, h_e;
	MatrixDense thickness(DOFs, 1), U(DOFs, 2), K(DOFs, 4), gpK(1, 4), CD;
	MatrixDense u(1, 2), matFlux(2, 1), matGradient(2, 1);

	const MaterialBaseConfiguration *phase1 = NULL, *phase2 = NULL;
	if (iterator.material->phase_change) {
		phase1 = &iterator.material->phases.find(1)->second;
		phase2 = &iterator.material->phases.find(2)->second;
		iterator.phase.data[0] = 0;
		iterator.latentHeat.data[0] = 0;
	}

	for (int n = 0; n < DOFs; n++) {
		Evaluator::Params params;
		params.coords(2, iterator.coordinates.data + 2 * n);
		params.temp(iterator.temperature.data + n);
		T(n, 0) = iterator.temperature.data[n];
		coordinates(n, 0) = iterator.coordinates.data[2 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[2 * n + 1];
		thickness(n, 0) = iterator.thickness.data[n];
		if (iterator.material->phase_change) {
			double phase, derivation;
			smoothstep(phase, derivation, iterator.material->phase_change_temperature - iterator.material->transition_interval / 2, iterator.material->phase_change_temperature + iterator.material->transition_interval / 2, T(n, 0), iterator.material->smooth_step_order);
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, phase1, phase, step::time.current, T(n, 0), K, CD, false);
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, phase2, (1 - phase), step::time.current, T(n, 0), K, CD, false);
			double dens1 = phase1->density.evaluator->eval(params);
			double dens2 = phase2->density.evaluator->eval(params);
			double hc1 = phase1->heat_capacity.evaluator->eval(params);
			double hc2 = phase2->heat_capacity.evaluator->eval(params);
			iterator.phase.data[0] += phase;
			iterator.latentHeat.data[0] += iterator.material->latent_heat * derivation;
			m = (phase * dens1 + (1 - phase) * dens2) * (phase * hc1 + (1 - phase) * hc2 + iterator.material->latent_heat * derivation) * iterator.thickness.data[0];
		} else {
			assembleMaterialMatrix(n, iterator.coordinates.data + 2 * n, iterator.material, 1, step::time.current, T(n, 0), K, CD, false);
			double dens = iterator.material->density.evaluator->eval(params);
			double hc = iterator.material->heat_capacity.evaluator->eval(params);
			m = dens * hc * thickness(n, 0);
		}

		U(n, 0) = iterator.motion.data[2 * n + 0] * m;
		U(n, 1) = iterator.motion.data[2 * n + 1] * m;
	}
	if (iterator.material->phase_change) {
		iterator.phase.data[0] /= DOFs;
		iterator.latentHeat.data[0] /= DOFs;
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = MATH::determinant2x2(J.vals);
		if (detJ < 0) { detJ = -detJ; }
		MATH::Dense2x2inverse(J.vals, invJ.vals, detJ);
		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		norm_u_e = u.norm();
		h_e = 0;

		if (norm_u_e != 0) {
			MatrixDense b_e(1, DOFs);
			b_e.multiply(u, dND, 1, 0);
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += gsettings.sigma * h_e * norm_u_e;
		Ce(1, 1) += gsettings.sigma * h_e * norm_u_e;

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
	}

	if (info::ecf->output.results_selection.flux) {
		iterator.flux.data[0] = matFlux(0, 0) / N.size();
		iterator.flux.data[1] = matFlux(1, 0) / N.size();
	}
}





