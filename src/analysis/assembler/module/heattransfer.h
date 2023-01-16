
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_

#include "assembler.h"
#include "config/ecf/physics/heattransfer.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include <cstddef>
#include <map>

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct SteadyState;

class HeatTransfer: public Assembler
{
public:
	struct NGP {
		static const size_t POINT1 = 1;

		static const size_t LINE2 = 2;
		static const size_t LINE3 = 3;

		static const size_t TRIANGLE3 = 6;
		static const size_t TRIANGLE6 = 6;
		static const size_t SQUARE4   = 4;
		static const size_t SQUARE8   = 9;

		static const size_t TETRA4    = 4;
		static const size_t TETRA10   = 15;
		static const size_t PYRAMID5  = 8;
		static const size_t PYRAMID13 = 14;
		static const size_t PRISMA6   = 9;
		static const size_t PRISMA15  = 9;
		static const size_t HEXA8     = 8;
		static const size_t HEXA20    = 8;
	};

	HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme);
	void updateSolution(SteadyState &scheme);

	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	struct {
		ElementParameter<egps> weight;
		ElementParameter<enodes * egps> N;
		ElementParameter<edim * enodes * egps> dN;

		ElementParameter<egps> jacobiDeterminant;
		ElementParameter<ndim * ndim * egps> jacobiInversion;
		ElementParameter<edim * enodes * egps> dND;

		struct Boundary {
			BoundaryParameter<egps> weight;
			BoundaryParameter<enodes * egps> N;
			BoundaryParameter<edim * enodes * egps> dN;

			BoundaryParameter<egps> jacobian;
		} boundary;
	} integration;

	struct {
		ElementParameter<ndim * enodes> node;
		ElementParameter<ndim * egps> gp;

		struct Boundary {
			BoundaryParameter<ndim * enodes> node;
			BoundaryParameter<ndim * egps> gp;
		} boundary;
	} coords;

	struct ParametersThickness {
		ElementGPsExternalParameter<egps> gp;

		struct {
			BoundaryExternalParameter<egps> gp;
		} boundary;
	} thickness;

	struct ParametersCoordinateSystem {
		ElementGPsExternalParameter<egps> cartesian2D;
		ElementGPsExternalParameter<ndim * egps> cartesian3D, spherical;
		ElementGPsExternalParameter<2 * egps> cylindric;

		ElementParameter<2 * egps> angle2D;
		ElementParameter<6 * egps> angle3D;
	} cooSystem;

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
	} material;

	struct ParametersTemperature {
		ElementNodesExternalParameter<enodes> node;
		ElementGPsExternalParameter<egps> gp;

		struct Initial {
			ElementNodesExternalParameter<enodes> node;
			ElementGPsExternalParameter<egps> gp;
			struct Boundary {
				BoundaryParameter<enodes> node;
				BoundaryParameter<egps> gp;
			} boundary;
		} initial;

		struct Boundary {
			BoundaryParameter<enodes> node;
			BoundaryParameter<egps> gp;
		} boundary;
	} temp;

	struct ParametersTranslationMotions {
		double sigma;
		bool CAU, SUPG;

		ElementGPsExternalParameter<ndim * egps> gp;
		ElementParameter<enodes * enodes> stiffness;
		ElementParameter<enodes> rhs;
	} translationMotions;

	struct ParametersElementNodeFunction {
		ElementGPsExternalParameter<egps> gp;
	} heatSource;

	struct ParametersConvection {
		// input
		struct {
			BoundaryExternalParameter<egps> gp;
		} wallHeight, tiltAngle, diameter, plateLength, fluidVelocity, plateDistance, length, experimentalConstant, volumeFraction, absolutePressure;
		// output
		struct {
			BoundaryExternalParameter<egps> gp;
		} heatTransferCoeficient, externalTemperature;

		struct {
			BoundaryParameter<egps> gp;
		} rho, dynamicViscosity, dynamicViscosityTemp, heatCapacity, thermalConductivity;

		BoundarySettings<ConvectionConfiguration> configuration;
	} convection;

	struct {
		BoundaryExternalParameter<enodes> node;
	} temperature;

	struct {
		BoundaryExternalParameter<egps> gp;
	} heatFlow, heatFlux, q;

	struct ParametersElements {
		ElementParameter<enodes * enodes> stiffness;
		ElementParameter<enodes * enodes> mass;
		ElementParameter<enodes> rhs;

		struct {
			BoundaryParameter<enodes * enodes> stiffness;
			BoundaryParameter<enodes> rhs;
		} boundary;
	} elements;

	struct  {
		ElementParameter<egps> xi;
	} gradient;

	struct Results {
		static NodeData *temperature, *initialTemperature;
		static ElementData *translationMotion, *gradient, *flux;
	};

protected:
	bool initTemperature();
	void initParameters();
	void initNames();
	void printVersions();

	void _evaluate();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_ */
