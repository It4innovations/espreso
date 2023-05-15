
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_

#include "subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"

namespace espreso {

struct Plasticity: SubKernel {
	Plasticity()
	: behaviour(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN), configuration(nullptr),
	  isPlastized(nullptr), isPlastizedEnd(nullptr), nrhs(nullptr)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::ITERATION | Assembler::SOLUTION;
	}

	void activate(size_t interval, StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour, const PlasticityPropertiesConfiguration *configuration, NamedData *isPlastized, double *nrhs)
	{
		this->behaviour = behaviour;
		this->configuration = configuration;
		this->isconst = false;
		this->isactive = true;
		this->isPlastized = isPlastized->data.data() + info::mesh->elements->eintervals[interval].begin;
		this->isPlastizedEnd = isPlastized->data.data() + isPlastized->data.size();
		this->nrhs = nrhs;
	}

	StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour;
	const PlasticityPropertiesConfiguration *configuration;
	double *isPlastized, *isPlastizedEnd;
	double *nrhs;

	std::vector<double> smallStrainTensorPlastic, xi;
};

struct PlasticityStorage: Plasticity {

	PlasticityStorage(Plasticity &base, Assembler::Action action)
	: Plasticity(base),
	  smallStrainTensorPlastic(base.smallStrainTensorPlastic.data()),
	  xi(base.xi.data()),
	  save(action == Assembler::SOLUTION)
	{
		if (isactive) {
			size_t elements = base.smallStrainTensorPlastic.size() / 6;

			size_t size = sizeof(double) * (elements + 1);
			void *_smallStrainTensorPlastic = static_cast<void*>(smallStrainTensorPlastic);
			smallStrainTensorPlastic = static_cast<double*>(std::align(SIMD::size * sizeof(double), elements * sizeof(double), _smallStrainTensorPlastic, size));

			size = sizeof(double) * (elements + 1);
			void *_xi = static_cast<void*>(xi);
			xi = static_cast<double*>(std::align(SIMD::size * sizeof(double), elements * sizeof(double), _xi, size));
		}
	}

	double *smallStrainTensorPlastic, *xi;
	bool save;
};

template <size_t nodes, size_t gps, size_t ndim, enum ElasticityModel model, class Physics> struct PlasticityKernel: PlasticityStorage, Physics {
	PlasticityKernel(Plasticity &base, Assembler::Action action): PlasticityStorage(base, action) {}

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t gps, class Physics> struct PlasticityKernel<nodes, gps, 3, ElasticityModel::ISOTROPIC, Physics>: PlasticityStorage, Physics {
	PlasticityKernel(Plasticity &base, Assembler::Action action): PlasticityStorage(base, action), csqr32(std::sqrt(3./2)), crsqr32(1/csqr32), csqrt6(std::sqrt(6.)) {}

	double csqr32, crsqr32, csqrt6;

//	function [C_el,C_pla,sigma,plas] = radial_return_small(eps_, plasticity, material)
//	nqp  = size(eps_,2);
//	% von Mises linear mixed hardening
//	IxI  = [ones(3) zeros(3);zeros(3) zeros(3)];                                                       % 6 . 6
	// 1 1 1 0 0 0
	// 1 1 1 0 0 0
	// 1 1 1 0 0 0
	// 0 0 0 0 0 0
	// 0 0 0 0 0 0
	// 0 0 0 0 0 0

//	Idev = diag([1 1 1 1/2 1/2 1/2])-IxI/3;                                                            % 6 . 6
	//  2/3 -1/3 -1/3   0   0   0
	// -1/3  2/3 -1/3   0   0   0
	// -1/3 -1/3  2/3   0   0   0
	//    0    0    0 1/2   0   0
	//    0    0    0   0 1/2   0
	//    0    0    0   0   0 1/2

//	Idev_voi_operator = eye(6)      -IxI/3;
	//  2/3 -1/3 -1/3   0   0   0
	// -1/3  2/3 -1/3   0   0   0
	// -1/3 -1/3  2/3   0   0   0
	//    0    0    0   1   0   0
	//    0    0    0   0   1   0
	//    0    0    0   0   0   1
//	C_el = material.bulk_modulus*IxI + 2*material.shear_modulus*Idev;                                  % 6 . 6

//	% elastic predictor
//	eps_tr      = eps_-plasticity.smallStrainTensorPlastic;                                            % 6 . q
//	sigma_tr    = C_el*eps_tr;                                                                         % 6 . q
//	sigmadev_tr = Idev_voi_operator*sigma_tr;                                                          % 6 . q
//	eta_tr      = sigmadev_tr - plasticity.kinematic.xi(2:end,:);                                      % 6 . q
//	norm_eta_tr = sqrt(sum(eta_tr(1:3,:).^2 + 2*eta_tr(4:6,:).^2,1));                                  % 1 . q

//	fyield_tr   = sqrt(3/2)*(norm_eta_tr - material.Y);                                                % 1 . q
//	denom       = 3*material.shear_modulus+material.Hkinematic+material.Hisotropic;
//	Dgamma      = fyield_tr/denom;
//	q_tr        = sqrt(3/2)*norm_eta_tr;

//	ind_p       = fyield_tr > 0;
//	sigma       = sigma_tr;                                                                            % 6 . q
//	C_pla       = C_el;
//	% plastic corrector
//	if np >0
//	  n_tr      = eta_tr(:,ind_p)./norm_eta_tr(1,ind_p);                                               % 6 . p
//	  nn_tr     = pagemtimes(reshape(n_tr,6,1,np),reshape(n_tr,1,6,np));                               % 6 . 6 . p
//	  tmp_const = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p));                       % 1 . p
//	  tmp_vec   = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p)-1/denom);               % 1 . p
//	  C_pla(:,:,ind_p) = C_pla(:,:,ind_p) ...
//	    - reshape(reshape( Idev,[], 1) *tmp_const,6,6,np)     ...
//	    + reshape(reshape(nn_tr,36,np).*tmp_vec  ,6,6,np);
//	  sigma(  :,ind_p) = sigma(  :,ind_p) - (sqrt(6)*material.shear_modulus*Dgamma(:,ind_p)).*n_tr;
//	end
//	if nargout > 1
//	  plas = plasticity;
//	  if np > 0
//	    plas.isplastized(:   ,ind_p) = true;
//	    plas.smallStrainTensorPlastic(      :   ,ind_p) = plas.smallStrainTensorPlastic(:   ,ind_p) + sqrt(3/2)*[1;1;1;2;2;2].*(     n_tr.*Dgamma(:,ind_p));
//	    plas.xi(        2:end,ind_p) = plas.xi(  2:end,ind_p) + sqrt(2/3)*material.Hkinematic*(n_tr.*Dgamma(:,ind_p));
//	    plas.xi(        1    ,ind_p) = plas.xi(  1    ,ind_p) +                                      Dgamma(:,ind_p);
//	  end
//	end
	void simd(typename Physics::Element &element)
	{
		SIMD c1  = load1(1);
		SIMD c2  = load1(2);
		SIMD c3  = load1(3);
		SIMD c6  = load1(6);
		SIMD c12 = load1(1/2.);
		SIMD c13 = load1(1/3.);
		SIMD c23 = load1(2/3.);
		SIMD csqr32 = load1(this->csqr32);
		SIMD crsqr32 = load1(this->crsqr32);
		SIMD csqrt6 = load1(this->csqrt6);
		size_t size = std::min((size_t)SIMD::size, (size_t)(isPlastizedEnd - isPlastized));

		double * __restrict__ smallStrainTensorPlastic = this->smallStrainTensorPlastic;
		double * __restrict__ xi = this->xi;
		double * __restrict__ nrhs = this->nrhs;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD shearModulus = element.ecf.elasticity[gp][2];

			SIMD eps0 = element.smallStrainTensor[gp][0] - load(smallStrainTensorPlastic + (gp * 6 + 0) * SIMD::size);
			SIMD eps1 = element.smallStrainTensor[gp][1] - load(smallStrainTensorPlastic + (gp * 6 + 1) * SIMD::size);
			SIMD eps2 = element.smallStrainTensor[gp][2] - load(smallStrainTensorPlastic + (gp * 6 + 2) * SIMD::size);
			SIMD eps3 = element.smallStrainTensor[gp][3] - load(smallStrainTensorPlastic + (gp * 6 + 3) * SIMD::size);
			SIMD eps4 = element.smallStrainTensor[gp][4] - load(smallStrainTensorPlastic + (gp * 6 + 4) * SIMD::size);
			SIMD eps5 = element.smallStrainTensor[gp][5] - load(smallStrainTensorPlastic + (gp * 6 + 5) * SIMD::size);
			//	sigma_tr    = C_el*eps_tr;
			SIMD sigma0 = element.ecf.elasticity[gp][0] * eps0 + element.ecf.elasticity[gp][1] * (eps1 + eps2);
			SIMD sigma1 = element.ecf.elasticity[gp][0] * eps1 + element.ecf.elasticity[gp][1] * (eps0 + eps2);
			SIMD sigma2 = element.ecf.elasticity[gp][0] * eps2 + element.ecf.elasticity[gp][1] * (eps0 + eps1);
			SIMD sigma3 = element.ecf.elasticity[gp][2] * eps3;
			SIMD sigma4 = element.ecf.elasticity[gp][2] * eps4;
			SIMD sigma5 = element.ecf.elasticity[gp][2] * eps5;

			//	sigmadev_tr = Idev_voi_operator*sigma_tr;
			//	eta_tr      = sigmadev_tr - plasticity.kinematic.xi(2:end,:);
			SIMD eta_tr0 =  c23 * sigma0 - c13 * sigma1 - c13 * sigma2 - load(xi + (gp * 6 + 0) * SIMD::size);
			SIMD eta_tr1 = -c13 * sigma0 + c23 * sigma1 - c13 * sigma2 - load(xi + (gp * 6 + 1) * SIMD::size);
			SIMD eta_tr2 = -c13 * sigma0 - c13 * sigma1 + c23 * sigma2 - load(xi + (gp * 6 + 2) * SIMD::size);
			SIMD eta_tr3 =  sigma3 - load(xi + (gp * 6 + 3) * SIMD::size);
			SIMD eta_tr4 =  sigma4 - load(xi + (gp * 6 + 4) * SIMD::size);
			SIMD eta_tr5 =  sigma5 - load(xi + (gp * 6 + 5) * SIMD::size);

			//	norm_eta_tr = sqrt(sum(eta_tr(1:3,:).^2 + 2*eta_tr(4:6,:).^2,1));
			SIMD eta_tr_c = eta_tr0 * eta_tr0 + eta_tr1 * eta_tr1 + eta_tr2 * eta_tr2 + c2 * (eta_tr3 * eta_tr3 + eta_tr4 * eta_tr4 + eta_tr5 * eta_tr5);
			SIMD eta_tr_norm = sqrt(eta_tr_c);
			SIMD eta_tr_rnorm = rsqrt14(eta_tr_c);

			//	fyield_tr   = sqrt(3/2)*norm_eta_tr - material.initialYieldStress;
			SIMD fyield_tr = csqr32 * eta_tr_norm - element.ecf.initialYieldStress[gp];
			//	denom       = 3*material.shear_modulus+material.Hkinematic+material.Hisotropic;
			SIMD rdenom = c1 / (c3 * shearModulus + element.ecf.isotropicHardening[gp] + element.ecf.kinematicHardening[gp]);
			//	Dgamma      = fyield_tr/denom;
			SIMD Dgamma = fyield_tr * rdenom;
			//	q_tr        = sqrt(3/2)*norm_eta_tr;
			SIMD rq_tr = crsqr32 * eta_tr_rnorm;

			eta_tr0 = eta_tr0 * eta_tr_rnorm;
			eta_tr1 = eta_tr1 * eta_tr_rnorm;
			eta_tr2 = eta_tr2 * eta_tr_rnorm;
			eta_tr3 = eta_tr3 * eta_tr_rnorm;
			eta_tr4 = eta_tr4 * eta_tr_rnorm;
			eta_tr5 = eta_tr5 * eta_tr_rnorm;

			//	  tmp_const = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p));                       % 1 . p
			//	  tmp_vec   = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p)-1/denom);               % 1 . p
			SIMD shear26 = c6 * shearModulus * shearModulus;
			SIMD c = shear26 * Dgamma * rq_tr;
			SIMD v = shear26 * (Dgamma * rq_tr - rdenom);

			// radial return begin
			SIMD pl = ispositive(fyield_tr);
			// C
			//  0  1  2  3  4  5       0 1 1 _ _ _
			//     6  7  8  9 10         0 1 _ _ _
			//       11 12 13 14           0 _ _ _
			//          15 16 17             2 _ _
			//             18 19               2 _
			//                20                 2
			element.elasticity[gp][ 0] = pl * (v * eta_tr0 * eta_tr0 - c23 * c) + element.ecf.elasticity[gp][0];
			element.elasticity[gp][ 1] = pl * (v * eta_tr0 * eta_tr1 + c13 * c) + element.ecf.elasticity[gp][1];
			element.elasticity[gp][ 2] = pl * (v * eta_tr0 * eta_tr2 + c13 * c) + element.ecf.elasticity[gp][1];
			element.elasticity[gp][ 3] = pl * (v * eta_tr0 * eta_tr3);
			element.elasticity[gp][ 4] = pl * (v * eta_tr0 * eta_tr4);
			element.elasticity[gp][ 5] = pl * (v * eta_tr0 * eta_tr5);
			element.elasticity[gp][ 6] = pl * (v * eta_tr1 * eta_tr1 - c23 * c) + element.ecf.elasticity[gp][0];
			element.elasticity[gp][ 7] = pl * (v * eta_tr1 * eta_tr2 + c13 * c) + element.ecf.elasticity[gp][1];
			element.elasticity[gp][ 8] = pl * (v * eta_tr1 * eta_tr3);
			element.elasticity[gp][ 9] = pl * (v * eta_tr1 * eta_tr4);
			element.elasticity[gp][10] = pl * (v * eta_tr1 * eta_tr5);
			element.elasticity[gp][11] = pl * (v * eta_tr2 * eta_tr2 - c23 * c) + element.ecf.elasticity[gp][0];
			element.elasticity[gp][12] = pl * (v * eta_tr2 * eta_tr3);
			element.elasticity[gp][13] = pl * (v * eta_tr2 * eta_tr4);
			element.elasticity[gp][14] = pl * (v * eta_tr2 * eta_tr5);
			element.elasticity[gp][15] = pl * (v * eta_tr3 * eta_tr3 - c12 * c) + element.ecf.elasticity[gp][2];
			element.elasticity[gp][16] = pl * (v * eta_tr3 * eta_tr4);
			element.elasticity[gp][17] = pl * (v * eta_tr3 * eta_tr5);
			element.elasticity[gp][18] = pl * (v * eta_tr4 * eta_tr4 - c12 * c) + element.ecf.elasticity[gp][2];
			element.elasticity[gp][19] = pl * (v * eta_tr4 * eta_tr5);
			element.elasticity[gp][20] = pl * (v * eta_tr5 * eta_tr5 - c12 * c) + element.ecf.elasticity[gp][2];

			SIMD plSigma = pl * csqrt6 * shearModulus * Dgamma;
			sigma0 = sigma0 - plSigma * eta_tr0;
			sigma1 = sigma1 - plSigma * eta_tr1;
			sigma2 = sigma2 - plSigma * eta_tr2;
			sigma3 = sigma3 - plSigma * eta_tr3;
			sigma4 = sigma4 - plSigma * eta_tr4;
			sigma5 = sigma5 - plSigma * eta_tr5;
			// radial return end

			//      0  1  2  3  4  5
			// B = dX  0  0 dY  0 dZ
			//      0 dY  0 dX dZ  0
			//      0  0 dZ  0 dY dX
			SIMD scale = element.det[gp] * load1(element.w[gp]);
			for (size_t n = 0; n < nodes; ++n) {
				SIMD fx = load(nrhs + (0 * nodes + n) * SIMD::size);
				SIMD fy = load(nrhs + (1 * nodes + n) * SIMD::size);
				SIMD fz = load(nrhs + (2 * nodes + n) * SIMD::size);
				fx = fx + scale * (element.dND[gp][n][0] * sigma0 + element.dND[gp][n][1] * sigma3 + element.dND[gp][n][2] * sigma5);
				fy = fy + scale * (element.dND[gp][n][1] * sigma1 + element.dND[gp][n][0] * sigma3 + element.dND[gp][n][2] * sigma4);
				fz = fz + scale * (element.dND[gp][n][2] * sigma2 + element.dND[gp][n][1] * sigma4 + element.dND[gp][n][0] * sigma5);
				store(nrhs + (0 * nodes + n) * SIMD::size, fx);
				store(nrhs + (1 * nodes + n) * SIMD::size, fy);
				store(nrhs + (2 * nodes + n) * SIMD::size, fz);
			}

			if (this->save) {
				double * __restrict__ isPlastized = this->isPlastized;
				for (size_t s = 0; s < size; ++s) {
					isPlastized[s] += std::max(isPlastized[s], pl[s]);
				}
				store(smallStrainTensorPlastic + (6 * gp + 0) * SIMD::size, csqr32 * eta_tr0 * Dgamma);
				store(smallStrainTensorPlastic + (6 * gp + 1) * SIMD::size, csqr32 * eta_tr1 * Dgamma);
				store(smallStrainTensorPlastic + (6 * gp + 2) * SIMD::size, csqr32 * eta_tr2 * Dgamma);
				store(smallStrainTensorPlastic + (6 * gp + 3) * SIMD::size, csqr32 * eta_tr3 * Dgamma * c2);
				store(smallStrainTensorPlastic + (6 * gp + 4) * SIMD::size, csqr32 * eta_tr4 * Dgamma * c2);
				store(smallStrainTensorPlastic + (6 * gp + 5) * SIMD::size, csqr32 * eta_tr5 * Dgamma * c2);
				store(xi + (6 * gp + 0) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr0 * Dgamma);
				store(xi + (6 * gp + 1) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr1 * Dgamma);
				store(xi + (6 * gp + 2) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr2 * Dgamma);
				store(xi + (6 * gp + 3) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr3 * Dgamma * c2);
				store(xi + (6 * gp + 4) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr4 * Dgamma * c2);
				store(xi + (6 * gp + 5) * SIMD::size, crsqr32 * element.ecf.kinematicHardening[gp] * eta_tr5 * Dgamma * c2);
			}
		}
		this->nrhs += SIMD::size * 3 * nodes;
		this->smallStrainTensorPlastic += SIMD::size * 6 * gps;
		this->xi += SIMD::size * 6 * gps;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_ */

