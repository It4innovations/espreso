
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
