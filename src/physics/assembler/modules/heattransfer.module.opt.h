
#ifndef SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_
#define SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_

#include "module.opt.h"
#include "module.parameters.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct HeatTransferLoadStepConfiguration;
struct Builder;

class HeatTransferModuleOpt: public ModuleOpt
{
public:
	static const int GP_LINE2 = 2;
	static const int GP_LINE3 = 3;

	static const int GP_TRIANGLE3 = 6;
	static const int GP_TRIANGLE6 = 6;
	static const int GP_SQUARE4   = 4;
	static const int GP_SQUARE8   = 9;

	static const int GP_TETRA4    = 4;
	static const int GP_TETRA10   = 15;
	static const int GP_PYRAMID5  = 8;
	static const int GP_PYRAMID13 = 14;
	static const int GP_PRISMA6   = 9;
	static const int GP_PRISMA15  = 9;
	static const int GP_HEXA8     = 8;
	static const int GP_HEXA20    = 8;

	static ElementData *phase, *latentHeat;
	static void createParameters();
	static void insertParameters(Evaluator *evaluator);

	HeatTransferModuleOpt(HeatTransferModuleOpt *previous, HeatTransferLoadStepConfiguration &configuration);
	~HeatTransferModuleOpt() {};

	void nextSubstep();
	void solutionChanged(Vectors *solution);

	void updateStiffness(double *K, esint *perm, int interval);

	void updateStiffness(double *K, esint *perm, int region, int interval);
	void updateRHS(double *RHS, esint *perm, int region, int interval);

	void fillElementsInterval(int interval);
//	void processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler);
	void processSolution();

	const HeatTransferLoadStepConfiguration &configuration;

	ParametersIntegration integration;
	ParametersCoordinates coords;
	ParametersThickness thickness;

	ParametersCoordinateSystem cooSystem;
	ParametersMaterial material;

	ParametersTemperature temp;
	ParametersTranslationMotions translationMotions;

	ParametersConvection convection;
	ParametersBoundaryFunction dirichlet, heatFlow, heatFlux, q;

	ParametersLinearSystem<1> linearSystem;

	ParametersGradient gradient;
	ParametersFlux flux;

protected:
	void initTemperature();
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_MODULES_HEATTRANSFER_MODULE_OPT_H_ */
