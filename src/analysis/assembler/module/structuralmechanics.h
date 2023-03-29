
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

#include <cstddef>
#include <map>

namespace espreso {

struct StructuralMechanicsConfiguration;
struct StructuralMechanicsLoadStepConfiguration;
struct SteadyState;

class StructuralMechanics: public Assembler
{
public:
	struct StructuralElementType {
		enum: int {
			SYMMETRIC_PLANE              = 0,
			SYMMETRIC_PLANE_AXISYMMETRIC = 1,
			SYMMETRIC_VOLUME             = 2,
			FACE                         = 4,
			EDGE                         = 5,
			EDGE_AXISYMMETRIC            = 6,
			NODE                         = 7
		};
	};
	StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme);
	void updateSolution(SteadyState &scheme);

	void dryrun();

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
		static NodeData *displacement;
		static ElementData *principalStress, *componentStress, *vonMisesStress;
	};

protected:
	int axisymmetric;

	template <int etype> Assembler::measurements instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	template <int etype> Assembler::measurements instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	Assembler::measurements instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	template <int etype> Assembler::measurements instantiateConditions2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	template <int etype> Assembler::measurements instantiateConditions3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	Assembler::measurements instantiateConditions(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	template <int etype> Assembler::measurements instantiateManual2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	template <int etype> Assembler::measurements instantiateManual3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	Assembler::measurements instantiateManual(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	Assembler::measurements conditionsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
	Assembler::measurements manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<!
			(ndim == 2 &&
			edim == 2 &&
			ETYPE == StructuralElementType::SYMMETRIC_PLANE), int>::type* = 0);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t ETYPE>
	Assembler::measurements manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements,
	typename std::enable_if<
			ndim == 2 &&
			edim == 2 &&
			ETYPE == StructuralElementType::SYMMETRIC_PLANE, int>::type* = 0);

	bool initDisplacement();
	void initParameters();

	void volume();
	size_t esize();

	std::vector<std::vector<double> > cossin_conditions;
private:
	void generateElasticity();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_ */
