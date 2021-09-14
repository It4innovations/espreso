
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_

#include "assembler.h"
#include "config/ecf/physics/acoustic.h"
#include "math2/primitives/vector_sparse.h"
#include "math2/primitives/matrix_info.h"
#include "math2/generalization/matrix_base.h"

#include <vector>

namespace espreso {

struct AcousticGlobalSettings;
struct AcousticLoadStepConfiguration;
struct AX_HarmonicReal;

class AX_Acoustic: public Assembler
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

	AX_Acoustic(AX_Acoustic *previous, AcousticGlobalSettings &gsettings, AcousticLoadStepConfiguration &configuration);

	void init(AX_HarmonicReal &scheme, Vector_Base<double> *dirichlet);
	void analyze();
	void next();

	void updateSolution();

	Matrix_Type matrixType() { return Matrix_Type::REAL_SYMMETRIC_INDEFINITE; }
	bool hasKernel(int domain) { return true; }

	AcousticGlobalSettings &gsettings;
	AcousticLoadStepConfiguration &configuration;

	ParametersAcousticPressure acoustic_pressure;

	ParametersIntegration integration;
	ParametersCoordinates coords;

	ParametersBoundaryNodeFunction pressure;
	ParametersBoundaryFunction normalAcceleration, impedance, q;

	ParametersElements<1> elements;

	std::vector<esint> dirichletIndices, dirichletPermutation;

	struct Fragment {
		Vector_Base<double> *rhs, *x;
	};

	Matrix_Base<double> *K, *M, *C;
	Fragment re, im;
	Vector_Base<double> *dirichlet;

protected:
	void initParameters();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ACOUSTIC_H_ */
