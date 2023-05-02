
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_

#include "assembler.h"
#include "analysis/assembler/parameter.h"
#include "config/ecf/physics/heattransfer.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include "analysis/assembler/subkernel/basis.h"
#include "analysis/assembler/subkernel/thickness.h"
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
#include "analysis/assembler/subkernel/filler.h"
#include "analysis/assembler/subkernel/heattransfer/conductivity.h"
#include "analysis/assembler/subkernel/heattransfer/coordinatesystem.h"
#include "analysis/assembler/subkernel/heattransfer/advection.h"
#include "analysis/assembler/subkernel/heattransfer/matrix.h"
#include "analysis/assembler/subkernel/heattransfer/flux.h"
#include "analysis/assembler/subkernel/heattransfer/gradient.h"

#include <cstddef>
#include <map>

#include <type_traits>

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct SteadyState;

class HeatTransfer: public Assembler
{
public:
	struct SubKernels {
		int code;
		size_t elements, chunks;

		size_t esize;
		double volume;

		Basis basis;
		Thickness thickness;
		Coordinates coordinates;
		Temperature temperature;
		Integration integration;
		Conductivity conductivity;
		HeatTransferCoordinateSystem coosystem;
		Advection advection;
		HeatTransferMatrix K;

		TemperatureGradient gradient;
		TemperatureFlux flux;

		DataFiller Kfiller, RHSfiller;

		std::vector<ExternalEvaluator*> expressions;
	};

	struct BoundarySubKernels {
		int code;
		size_t elements, chunks;

		size_t esize;
		double surface;

		Basis basis;
		Thickness thickness;
		Coordinates coordinates;
		Integration integration;

		ExternalExpression temperature;

		DataFiller RHSfiller, dirichlet;

		std::vector<ExternalEvaluator*> expressions;
	};

	HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme, step::Time &time);
	void updateSolution(SteadyState &scheme);

	HeatTransferConfiguration &settings;
	HeatTransferLoadStepConfiguration &configuration;

	struct ParametersElements {
		ElementParameter<enodes * enodes> stiffness;
		ElementParameter<enodes * enodes> mass;
		ElementParameter<enodes> rhs;

		struct {
			BoundaryParameter<enodes * enodes> stiffness;
			BoundaryParameter<enodes> rhs, dirichlet;
		} boundary;
	} elements;

	struct Results {
		static NodeData *temperature, *initialTemperature, *thickness;
		static ElementData *translationMotion, *gradient, *flux;
	};

protected:
	void run(Action action, size_t interval);
	void run(Action action, size_t region, size_t interval);

	void runPlane(Action action, size_t interval);
	void runAxisymmetric(Action action, size_t interval);
	void runVolume(Action action, size_t interval);
	void runPreprocess(Action action, size_t interval);
	void runBoundary(Action action, size_t region, size_t interval);

	std::vector<SubKernels> subkernels;
	std::vector<std::vector<BoundarySubKernels> > boundary;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_ */

