
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETERS_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETERS_H_

#include "analysis/assembler/parameter.h"
#include "math/simd/simd.h"

namespace espreso {

class ElementData;
class ConvectionConfiguration;

template <int dimension>
struct ParametersElements {
	struct Boundary {
		BoundaryParameter<dimension * dimension * enodes * enodes> stiffness;
		BoundaryParameter<dimension * enodes> rhs;
	};

	ElementParameter<dimension * dimension * enodes * enodes> stiffness;
	ElementParameter<dimension * dimension * enodes * enodes> mass;
	ElementParameter<dimension * dimension * enodes * enodes> damping;
	ElementParameter<dimension * enodes> rhs;

	Boundary boundary;
};

struct ParametersTemperature {
	static NodeData *output;

	struct Boundary {
		BoundaryParameter<enodes> node;
		BoundaryParameter<egps> gp;
	};

	struct Initial {
		static NodeData *output;

		ElementNodesExternalParameter<enodes> node;
		ElementGPsExternalParameter<egps> gp;
		Boundary boundary;
	};

	Initial initial;
	ElementNodesExternalParameter<enodes> node;
	ElementGPsExternalParameter<egps> gp;
	Boundary boundary;
};

struct ParametersAcousticPressure {
	static NodeData *output;

	struct Boundary {
		BoundaryParameter<enodes> node;
		BoundaryParameter<egps> gp;
	};

	struct Initial {
		static NodeData *output;

		ElementGPsExternalParameter<enodes> node;
		ElementGPsExternalParameter<egps> gp;
		Boundary boundary;
	};

	Initial initial;
	ElementParameter<enodes> node;
	ElementGPsExternalParameter<egps> gp;
	Boundary boundary;
};


struct ParametersCoordinates {
	struct Boundary {
		BoundaryParameter<ndim * enodes> node;
		BoundaryParameter<ndim * egps> gp;
	};

	ElementParameter<ndim * enodes> node;
	ElementGPsExternalParameter<ndim * egps> gp;
	Boundary boundary;
};

struct ParametersIntegration {
	struct Boundary {
		BoundaryParameter<egps> weight;
		BoundaryParameter<enodes * egps> N;
		BoundaryParameter<edim * enodes * egps> dN;

		BoundaryParameter<egps> jacobian;
	};

	ElementParameter<egps> weight;
	ElementParameter<enodes * egps> N;
	ElementParameter<edim * enodes * egps> dN;

	ElementParameter<egps> jacobiDeterminant;
	ElementParameter<ndim * ndim * egps> jacobiInversion;
	ElementParameter<edim * enodes * egps> dND;

	Boundary boundary;
};

struct ParametersIntegrationSimd {
	struct Boundary {
		BoundaryParameter<egps> weight;
		BoundaryParameter<enodes * egps> N;
		BoundaryParameter<edim * enodes * egps> dN;

		BoundaryParameter<egps> jacobian;
	};

	ElementParameter<egps> weight;
	ElementParameter<enodes * egps + SIMD::size> N;
	ElementParameter<edim * enodes * egps + SIMD::size> dN;

	ElementParameter<egps> jacobiDeterminant;
	ElementParameter<ndim * ndim * egps> jacobiInversion;
	ElementParameter<edim * enodes * egps> dND;

	Boundary boundary;
};

struct ParametersCoordinateSystem {
	ElementGPsExternalParameter<egps> cartesian2D;
	ElementGPsExternalParameter<ndim * egps> cartesian3D, spherical;
	ElementGPsExternalParameter<2 * egps> cylindric;
};

struct ParametersMaterial {
	struct Model {
		ElementGPsExternalParameter<egps> isotropic;
		ElementGPsExternalParameter<ndim * egps> diagonal;
		ElementGPsExternalParameter<3 * egps> symmetric2D;
		ElementGPsExternalParameter<6 * egps> symmetric3D;
		ElementGPsExternalParameter<ndim * ndim * egps> anisotropic;
	};

	ElementGPsExternalParameter<egps> density, heatCapacity;
	Model model;

	ElementParameter<egps> mass;
	ElementParameter<egps> conductivityIsotropic;
	ElementParameter<ndim * ndim * egps> conductivity;
};

struct ParametersThickness {
	struct Boundary {
		BoundaryExternalParameter<egps> gp;
	};

	ElementGPsExternalParameter<egps> gp;
	Boundary boundary;
};

struct ParametersTranslationMotions {
	static ElementData *output;
	double sigma;
	bool CAU, SUPG;

	ElementGPsExternalParameter<ndim * egps> gp;
	ElementParameter<enodes * enodes> stiffness;
	ElementParameter<enodes> rhs;
};

struct ParametersElementNodeFunction {
	ElementGPsExternalParameter<enodes> node;
};

struct ParametersBoundaryFunction {
	BoundaryExternalParameter<egps> gp;
};

struct ParametersConvection {
	struct ExternalParameter {
		BoundaryExternalParameter<egps> gp;
	};

	struct InternalParameter {
		BoundaryParameter<egps> gp;
	};

	// input
	ExternalParameter wallHeight, tiltAngle, diameter, plateLength, fluidVelocity, plateDistance, length, experimentalConstant, volumeFraction, absolutePressure;
	// output
	ExternalParameter heatTransferCoeficient, externalTemperature;

	InternalParameter rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity;

	BoundarySettings<ConvectionConfiguration> configuration;
};

struct ParametersGradient {
	static ElementData *output;

	ElementParameter<egps> xi;
};

struct ParametersFlux {
	static ElementData *output;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_PARAMETERS_H_ */
