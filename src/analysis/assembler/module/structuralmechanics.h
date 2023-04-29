
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_

#include "assembler.h"
#include "analysis/assembler/parameter.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include "analysis/assembler/subkernel/basis.h"
#include "analysis/assembler/subkernel/thickness.h"
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/material.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
#include "analysis/assembler/subkernel/filler.h"
#include "analysis/assembler/subkernel/structuralmechanics/coordinatesystem.h"
#include "analysis/assembler/subkernel/structuralmechanics/displacement.h"
#include "analysis/assembler/subkernel/structuralmechanics/elasticity.h"
#include "analysis/assembler/subkernel/structuralmechanics/plasticity.h"
#include "analysis/assembler/subkernel/structuralmechanics/matrix.h"
#include "analysis/assembler/subkernel/structuralmechanics/rhs.h"
#include "analysis/assembler/subkernel/structuralmechanics/stress.h"

#include <cstddef>
#include <map>

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;
struct SteadyState;

class StructuralMechanics: public Assembler
{
public:
	struct SubKernels {
		int code;
		size_t elements, chunks;

		size_t esize;
		double volume;

		Basis basis;
		Thickness thickness;
		Material material;
		Coordinates coordinates;
		Displacement displacement;
		Temperature temperature;
		Integration integration;
		Elasticity elasticity;
		Plasticity plasticity;
		StructuralMechanicsCoordinateSystem coosystem;
		StructuralMechanicsMatrix K;
		StructuralMechanicsRHS acceleration, angularVelocity;
		Stress stress;

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

		ExternalExpressionVector displacement;
		StructuralMechanicsRHS normalPressure;

		DataFiller RHSfiller, dirichlet;

		std::vector<ExternalEvaluator*> expressions;
	};

	StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme, step::Time &time);
	void updateSolution(SteadyState &scheme);

	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

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

	struct Results {
		static NodeData *displacement, *thickness;
		static ElementData *principalStress, *componentStress, *vonMisesStress;
	};

protected:
	void run(Action action, size_t interval);
	void run(Action action, size_t region, size_t interval);

	std::vector<SubKernels> subkernels;
	std::vector<std::vector<BoundarySubKernels> > boundary;
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_ */
