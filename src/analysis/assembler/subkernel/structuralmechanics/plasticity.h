
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_

#include "subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct Plasticity: SubKernel {
	Plasticity()
	: behaviour(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN), configuration(nullptr)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE | Assembler::SOLUTION;
	}

	void activate(StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour, const PlasticityPropertiesConfiguration *configuration)
	{
		this->behaviour = behaviour;
		this->configuration = configuration;
		this->isconst = false;
	}

	StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR behaviour;
	const PlasticityPropertiesConfiguration *configuration;

	std::vector<double> scale, eps, xi;
};

template <size_t gps, size_t ndim, enum ElasticityModel model, class Physics> struct PlasticityKernel: Plasticity, Physics {
	PlasticityKernel(Plasticity &base): Plasticity(base), scale(base.scale.data()), eps(base.eps.data()), xi(base.xi.data()) {}

	double *scale, *eps, *xi;

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t gps, class Physics> struct PlasticityKernel<gps, 3, ElasticityModel::ISOTROPIC, Physics>: Plasticity, Physics {
	PlasticityKernel(Plasticity &base): Plasticity(base), scale(base.scale.data()), eps(base.eps.data()), xi(base.xi.data()), csqr32(std::sqrt(3./2)), crsqr32(1/csqr32) {}

	double *scale, *eps, *xi;
	double csqr32, crsqr32;

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
	// eps = B * displacement
	// plasticity.eps_p from memory
//	eps_tr      = eps_-plasticity.eps_p;                                                               % 6 . q
//	sigma_tr    = C_el*eps_tr;                                                                         % 6 . q
//	sigmadev_tr = Idev_voi_operator*A;                                                          % 6 . q

//	plasticity.xi from memory
//	eta_tr      = sigmadev_tr - plasticity.kinematic.xi(2:end,:);                                      % 6 . q
//	norm_eta_tr = sqrt(sum(eta_tr(1:3,:).^2 + 2*eta_tr(4:6,:).^2,1));                                  % 1 . q

//	ecf parameter
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
//	    plas.eps_p(      :   ,ind_p) = plas.eps_p(:   ,ind_p) + sqrt(3/2)*[1;1;1;2;2;2].*(     n_tr.*Dgamma(:,ind_p));
//	    plas.xi(        2:end,ind_p) = plas.xi(  2:end,ind_p) + sqrt(2/3)*material.Hkinematic*(n_tr.*Dgamma(:,ind_p));
//	    plas.xi(        1    ,ind_p) = plas.xi(  1    ,ind_p) +                                      Dgamma(:,ind_p);
//	  end
//	end

	// C
	// 0 1 1 _ _ _
	//   0 1 _ _ _
	//     0 _ _ _
	//       2 _ _
	//         2 _
	//           2
	void simd(typename Physics::Element &element)
	{
		printf("plasticity\n");
		SIMD c1  = load1(1);
		SIMD c2  = load1(2);
		SIMD c3  = load1(3);
		SIMD c6  = load1(6);
		SIMD c12 = load1(1/2.);
		SIMD c13 = load1(1/3.);
		SIMD c23 = load1(2/3.);
		SIMD csqr32 = load1(this->csqr32);
		SIMD crsqr32 = load1(this->crsqr32);
		double * __restrict__ scale = this->scale;
		double * __restrict__ eps = this->eps;
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD eps0 = element.dispTensor[gp][0] - load(eps + (gp * 6 + 0) * SIMD::size);
			SIMD eps1 = element.dispTensor[gp][1] - load(eps + (gp * 6 + 1) * SIMD::size);
			SIMD eps2 = element.dispTensor[gp][2] - load(eps + (gp * 6 + 2) * SIMD::size);
			SIMD eps3 = element.dispTensor[gp][3] - load(eps + (gp * 6 + 3) * SIMD::size);
			SIMD eps4 = element.dispTensor[gp][4] - load(eps + (gp * 6 + 4) * SIMD::size);
			SIMD eps5 = element.dispTensor[gp][5] - load(eps + (gp * 6 + 5) * SIMD::size);
			SIMD sigma0 = element.ecf.elasticity[gp][0] * eps0 + element.ecf.elasticity[gp][1] * (eps1 + eps2);
			SIMD sigma1 = element.ecf.elasticity[gp][0] * eps1 + element.ecf.elasticity[gp][1] * (eps0 + eps2);
			SIMD sigma2 = element.ecf.elasticity[gp][0] * eps2 + element.ecf.elasticity[gp][1] * (eps0 + eps1);
			SIMD sigma3 = element.ecf.elasticity[gp][2] * eps3;
			SIMD sigma4 = element.ecf.elasticity[gp][2] * eps4;
			SIMD sigma5 = element.ecf.elasticity[gp][2] * eps5;

			SIMD eta_tr0 =  c23 * sigma0 - c13 * sigma1 - c13 * sigma2 - element.ecf.kinematic[gp][0];
			SIMD eta_tr1 = -c13 * sigma0 + c23 * sigma1 - c13 * sigma2 - element.ecf.kinematic[gp][1];
			SIMD eta_tr2 = -c13 * sigma0 - c13 * sigma1 + c23 * sigma2 - element.ecf.kinematic[gp][2];
			SIMD eta_tr3 =  sigma3 - element.ecf.kinematic[gp][3];
			SIMD eta_tr4 =  sigma4 - element.ecf.kinematic[gp][4];
			SIMD eta_tr5 =  sigma5 - element.ecf.kinematic[gp][5];

			SIMD eta_tr_c = eta_tr0 * eta_tr0 + eta_tr1 * eta_tr1 + eta_tr2 * eta_tr2 + c2 * (eta_tr3 * eta_tr3 + eta_tr4 * eta_tr4 + eta_tr5 * eta_tr5);
			SIMD eta_tr_norm = sqrt(eta_tr_c);
			SIMD eta_tr_rnorm = rsqrt14(eta_tr_c);
			SIMD fyield_tr = csqr32 * (eta_tr_norm - element.ecf.Y[gp]);
			SIMD rdenom = c1 / (c3 * element.ecf.shearModulus[gp]); // + material.Hkinematic+material.Hisotropic;
			SIMD Dgamma = fyield_tr * rdenom;
			SIMD q_tr = crsqr32 * eta_tr_rnorm;

			SIMD eta_tr[6] = {
				eta_tr0 * eta_tr_rnorm,
				eta_tr1 * eta_tr_rnorm,
				eta_tr2 * eta_tr_rnorm,
				eta_tr3 * eta_tr_rnorm,
				eta_tr4 * eta_tr_rnorm,
				eta_tr5 * eta_tr_rnorm,
			};

			SIMD shear26 = c6 * element.ecf.shearModulus[gp] * element.ecf.shearModulus[gp];
			for (size_t n = 0; n < 6; ++n) {
				for (size_t m = 0; m < 6; ++m) {
					SIMD nn_tr = eta_tr[n] * eta_tr[m];
					SIMD c = shear26 * Dgamma * q_tr;
					SIMD v = shear26 * Dgamma * q_tr - rdenom;
//					element.elasticity[n * 6 + m] = (c1 - scale)
				}
			}
			//	  n_tr      = eta_tr(:,ind_p)./norm_eta_tr(1,ind_p);                                               % 6 . p
			//	  nn_tr     = pagemtimes(reshape(n_tr,6,1,np),reshape(n_tr,1,6,np));                               % 6 . 6 . p
			//	  tmp_const = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p));                       % 1 . p
			//	  tmp_vec   = (6*material.shear_modulus^2)*(Dgamma(:,ind_p)./q_tr(:,ind_p)-1/denom);               % 1 . p
			//	  C_pla(:,:,ind_p) = C_pla(:,:,ind_p) ...
			//	    - reshape(reshape( Idev,[], 1) *tmp_const,6,6,np)     ...
			//	    + reshape(reshape(nn_tr,36,np).*tmp_vec  ,6,6,np);
			//	  sigma(  :,ind_p) = sigma(  :,ind_p) - (sqrt(6)*material.shear_modulus*Dgamma(:,ind_p)).*n_tr;

		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_PLASTICITY_H_ */

