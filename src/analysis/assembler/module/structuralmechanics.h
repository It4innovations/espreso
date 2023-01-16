
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_

#include "assembler.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "config/holders/expression.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include <cstddef>
#include <map>

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;
struct SteadyState;

class StructuralMechanics: public Assembler
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

	StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme);
	void updateSolution(SteadyState &scheme);

	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	struct {
		ElementParameter<ndim * enodes> node;
		ElementParameter<ndim * egps> gp;

		struct {
			BoundaryParameter<ndim * enodes> node;
			BoundaryParameter<ndim * egps> gp;
		} boundary;
	} coords;

	struct {
		ElementGPsExternalParameter<egps> gp;
		struct {
			BoundaryExternalParameter<egps> gp;
		} boundary;
	} thickness;

	struct {
		ElementParameter<egps> weight;
		ElementParameter<enodes * egps> N;
		ElementParameter<edim * enodes * egps> dN;

		ElementParameter<egps> jacobiDeterminant;
		ElementParameter<ndim * ndim * egps> jacobiInversion;
		ElementParameter<edim * enodes * egps> dND;

		struct {
			BoundaryParameter<egps> weight;
			BoundaryParameter<enodes * egps> N;
			BoundaryParameter<edim * enodes * egps> dN;
			BoundaryParameter<ndim * egps> normal;

			BoundaryParameter<egps> jacobian;
		} boundary;
	} integration;

	struct {
		struct {
			ElementGPsExternalParameter<egps> isoPoissonRatio, isoYoungModulus;
			ElementGPsExternalParameter<ndim * egps> poissonRatio, youngModulus, shearModulus;
			ElementGPsExternalParameter<21 * egps> anisotropic3D;
		} model;

		ElementGPsExternalParameter<egps> density;

		ElementParameter<egps> mass;
		ElementParameter<3 * 3 * egps> elasticityPlane;
		ElementParameter<4 * 4 * egps> elasticityAxisymm;
		ElementParameter<6 * 6 * egps> elasticity3D;
	} material;

	struct {
		ElementGPsExternalParameter<egps> cartesian2D;
		ElementGPsExternalParameter<ndim * egps> cartesian3D;

		ElementParameter<2 * egps> angle2D;
		ElementParameter<6 * egps> angle3D;
	} cooSystem;

	struct ParametersElements {
		ElementParameter<ndim * enodes * ndim * enodes> stiffness;
//		ElementParameter<ndim * enodes * ndim * enodes> mass;
		ElementParameter<ndim * enodes> rhs;

		struct {
			BoundaryParameter<ndim * enodes * ndim * enodes> stiffness;
//			BoundaryParameter<ndim * enodes * ndim * enodes> mass;
			BoundaryParameter<ndim * enodes> rhs;
		} boundary;
	} elements;

	struct {
		ElementGPsExternalParameter<ndim * egps> gp;
	} acceleration;

	struct {
		ElementGPsExternalParameter<3 * egps> gp;
	} angularVevocity;

	struct {
		BoundaryExternalParameter<ndim * enodes> node;
	} displacement;

	struct {
		BoundaryExternalParameter<egps> gp;
	} normalPressure;

	struct Results {
		static NodeData *displacement;
	};

protected:
	bool initTemperature();
	void initParameters();
	void initNames();
	void printVersions();

	void _evaluate();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_ */
