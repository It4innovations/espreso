
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
	StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme);
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
		static NodeData *displacement;
	};

protected:
	int axisymmetric;

	template <int etype> double instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	template <int etype> double instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	double instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	double operatorsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);
	template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
	double manualloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements);

	bool initDisplacement();
	void initParameters();

	void volume();
	size_t esize();

private:
	void generateElasticity();
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_ */
