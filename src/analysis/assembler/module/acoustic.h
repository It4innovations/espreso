
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_

#include "assembler.h"
#include "analysis/assembler/parameter.h"
#include "config/ecf/physics/acoustic.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include <vector>

namespace espreso {

struct AcousticConfiguration;
struct AcousticLoadStepConfiguration;
struct Harmonic;

class Acoustic: public Assembler
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

	Acoustic(Acoustic *previous, AcousticConfiguration &settings, AcousticLoadStepConfiguration &configuration);

	void analyze();

	void connect(Harmonic &scheme);
	void evaluate(Harmonic &scheme);
	void updateSolution(Harmonic &scheme);

	AcousticConfiguration &settings;
	AcousticLoadStepConfiguration &configuration;

//	struct ParametersAcousticPressure {
//		ElementParameter<enodes> node;
//		ElementGPsExternalParameter<egps> gp;
//
//		struct {
//			ElementGPsExternalParameter<enodes> node;
//			ElementGPsExternalParameter<egps> gp;
//		} initial;
//	} acoustic_pressure;
//
//	struct {
//		ElementParameter<egps> weight;
//		ElementParameter<enodes * egps> N;
//		ElementParameter<edim * enodes * egps> dN;
//
//		ElementParameter<egps> jacobiDeterminant;
//		ElementParameter<ndim * ndim * egps> jacobiInversion;
//		ElementParameter<edim * enodes * egps> dND;
//
//		struct {
//			BoundaryParameter<egps> weight;
//			BoundaryParameter<enodes * egps> N;
//			BoundaryParameter<edim * enodes * egps> dN;
//
//			BoundaryParameter<egps> jacobian;
//		} boundary;
//	} integration;
//
//	struct {
//		ElementParameter<ndim * enodes> node;
//		ElementParameter<ndim * egps> gp;
//		struct {
//			BoundaryParameter<ndim * enodes> node;
//			BoundaryParameter<ndim * egps> gp;
//		} boundary;
//	} coords;

//	struct {
//		BoundaryExternalParameter<enodes> node;
//	} pressure;
//
//	struct {
//		BoundaryExternalParameter<egps> gp;
//	} normalAcceleration, impedance, q, proj_acceleration;
//
//	struct {
//		BoundaryExternalParameter<ndim * egps> gp;
//	} acceleration, normals;
//
//	struct {
//		ElementGPsExternalParameter<egps> density, speed_of_sound;
//	} material;
//
//	struct {
//		BoundaryExternalParameter<enodes> node;
//	} pointSource;

//	struct {
//		ElementGPsExternalParameter<egps> gp;
//	} monopoleSource;
//
//	struct {
//		ElementGPsExternalParameter<ndim * egps> gp;
//	} dipoleSource;

	struct {
		ElementParameter<enodes * enodes> stiffness;
		ElementParameter<enodes * enodes> mass;
		ElementParameter<enodes * enodes> damping;
		ElementParameter<enodes> rhs;

		ElementParameter<enodes> monopole;
		ElementParameter<enodes> dipole;

		struct {
			BoundaryParameter<enodes * enodes> stiffness;
			BoundaryParameter<enodes * enodes> mass;
			BoundaryParameter<enodes> rhs;
		} boundary;
	} elements;

	struct Results {
		static NodeData *pressure, *initialPressure;
	};
protected:
	void initParameters();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_ */
