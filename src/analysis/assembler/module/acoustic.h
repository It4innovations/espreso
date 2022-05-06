
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_

#include "assembler.h"
#include "config/ecf/physics/acoustic.h"
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
		static const size_t POINT1 = 0;

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

	ParametersAcousticPressure acoustic_pressure;
	ParametersIntegration integration;
	ParametersCoordinates coords;

	ParametersBoundaryNodeFunction pressure;
	ParametersBoundaryFunction normalAcceleration, impedance, q, proj_acceleration;
	ParametersBoundaryVectorFunction acceleration;
	ParametersMaterial material;

	ParameterMonopoleSource monopoleSource;
	ParameterDipoleSource dipoleSource;

	ParametersBoundaryVectorFunction normals;

	ParametersElements<1> elements;
protected:
	void initParameters();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_ */
