
#ifndef SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_
#define SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_

#include "module.opt.h"
#include "module.parameters.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "config/ecf/physics/heattransfer.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct Builder;

class HeatTransferModuleOpt: public ModuleOpt
{
public:
	struct NGP {
		static const int LINE2 = 2;
		static const int LINE3 = 3;

		static const int TRIANGLE3 = 6;
		static const int TRIANGLE6 = 6;
		static const int SQUARE4   = 4;
		static const int SQUARE8   = 9;

		static const int TETRA4    = 4;
		static const int TETRA10   = 15;
		static const int PYRAMID5  = 8;
		static const int PYRAMID13 = 14;
		static const int PRISMA6   = 9;
		static const int PRISMA15  = 9;
		static const int HEXA8     = 8;
		static const int HEXA20    = 8;
	};

	static ElementData *phase, *latentHeat;
	static void createParameters();
	static void insertParameters(Evaluator *evaluator);

	HeatTransferModuleOpt(HeatTransferModuleOpt *previous, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);
	~HeatTransferModuleOpt() {};

	void nextSubstep();
	void solutionChanged(Vectors *solution);

	void updateStiffness(double *K, esint *perm, int interval);

	void updateStiffness(double *K, esint *perm, int region, int interval);
	void updateRHS(double *RHS, esint *perm, int region, int interval);

	void fillElementsInterval(int interval);
//	void processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler);
	void processSolution();

	const HeatTransferGlobalSettings &gsettings;
	const HeatTransferLoadStepConfiguration &configuration;

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
	ParametersBoundaryFunction dirichlet, heatFlow, heatFlux, q;

	ParametersElements<1> elements;
	ParametersElements<1> elementsSimd;

	ParametersGradient gradient;
	ParametersGradient gradientSimd;
	ParametersFlux flux;

protected:
	void initTemperature();
	void printVersions();
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_ */
