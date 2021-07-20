
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_

#include "assembler.h"
#include "math2/primitives/vector_sparse.h"
#include "math2/primitives/matrix_info.h"
#include "math2/generalization/matrix_base.h"

#include <cstddef>

namespace espreso {

struct HeatTransferGlobalSettings;
struct HeatTransferLoadStepConfiguration;
struct AX_SteadyState;

class AX_HeatTransfer: public Assembler
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

	AX_HeatTransfer(AX_HeatTransfer *previous, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);

	void init(AX_SteadyState &scheme);
	void analyze();
	void next();

	void initDirichlet(Vector_Sparse<double> &dirichlet);
	void fillDirichlet(Vector_Sparse<double> &dirichlet);

	void updateSolution();

	Matrix_Type matrixType() { return Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }
	bool hasKernel(int domain) { return true; }

	HeatTransferGlobalSettings &gsettings;
	HeatTransferLoadStepConfiguration &configuration;

	ParametersIntegration integration;
	ParametersIntegrationSimd integrationSimd;
	ParametersCoordinates coords;
	ParametersCoordinates coordsSimd;
	ParametersThickness thickness;
	ParametersThickness thicknessSimd;

	ParametersCoordinateSystem cooSystem;
	ParametersMaterial material;
	ParametersMaterial materialSimd;

	ParametersTemperature temp;
	ParametersTranslationMotions translationMotions;
	ParametersElementNodeFunction heatSource;

	ParametersConvection convection;
	ParametersBoundaryFunction temperature, heatFlow, heatFlux, q;

	ParametersElements<1> elements;
	ParametersElements<1> elementsSimd;

	ParametersGradient gradient;
	ParametersGradient gradientSimd;
	ParametersFlux flux;

	std::vector<esint> dirichletIndices, dirichletPermutation;

	Matrix_Base<double> *K, *M;
	Vector_Base<double> *rhs, *x;

protected:
	void initTemperature();
	void printVersions();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_ */
