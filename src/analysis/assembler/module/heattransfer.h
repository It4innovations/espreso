
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_

#include "assembler.h"
#include "parameter.h"
#include "config/ecf/physics/heattransfer.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "math/primitives/vector_sparse.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"

#include <cstddef>
#include <map>

namespace espreso {

struct HeatTransferConfiguration;
struct HeatTransferLoadStepConfiguration;
struct SteadyState;

class HeatTransfer: public Assembler
{
public:
	HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

	void analyze();

	void connect(SteadyState &scheme);
	void evaluate(SteadyState &scheme);
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
		static NodeData *temperature, *initialTemperature;
		static ElementData *translationMotion, *gradient, *flux;
	};

protected:
	template <int etype> double instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, esint elements);
	template <int etype> double instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, esint elements);
	double instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, esint elements);

	bool initTemperature();
	void initParameters();

	void volume();
	size_t esize();

private:
	void gatherInputs();
	void generateConductivity();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEATTRANSFER_H_ */

