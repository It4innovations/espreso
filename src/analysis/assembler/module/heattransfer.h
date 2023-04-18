
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
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
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
	struct TransferElementType {
		enum: int {
			SYMMETRIC_ISOTROPIC  = 0,
			SYMMETRIC_GENERAL    = 1,
			ASYMMETRIC_ISOTROPIC = 2,
			ASYMMETRIC_GENERAL   = 3,
			FACE                 = 4,
			EDGE                 = 5,
			NODE                 = 6
		};
	};

public:
	struct SubKernels {
		int code;
		size_t etype;
		esint elements;

		Basis basis;
		Coordinates coordinates;
		Temperature temperature;
		Integration integration;
		Conductivity conductivity;
		HeatTransferCoordinateSystem translation;
		Advection advection;
		HeatTransferMatrix K;

		TemperatureGradient gradient;
		TemperatureFlux flux;

		std::vector<ExternalExpression*> expressions;
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

	template <int etype> void instantiate2D(Action action, size_t interval);
	template <int etype> void instantiate3D(Action action, size_t interval);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	void hybridloop(Action action, size_t interval);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	void hybridpreprocess(size_t interval);

	void initParameters();

	void initTemperatureAndThickness();
	void volume();
	size_t esize();

	std::vector<std::vector<double> > cossin_conditions;
	std::vector<SubKernels> subkernels;
private:
	void generateConductivity();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_ */

