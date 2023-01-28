
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_HEAT_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_HEAT_H_

#include "assembler.h"
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

class Heat: public Assembler
{
public:
	struct NGP {
		static const size_t POINT1 = 1;

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

	Heat(Heat *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration);

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
			BoundaryParameter<enodes> rhs;
		} boundary;
	} elements;

	struct Results {
		static NodeData *temperature, *initialTemperature;
		static ElementData *translationMotion, *gradient, *flux;
	};

protected:
	double instantiate(size_t interval, int code, const std::vector<ActionOperator*> &ops, esint elements);

	bool initTemperature();
	void initParameters();
	void initNames();
	void printVersions();

	void volume();
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_HEAT_H_ */

