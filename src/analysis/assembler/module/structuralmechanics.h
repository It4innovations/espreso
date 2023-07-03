
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

#include "analysis/assembler/kernel/structuralmechanics/subkernellist.h"

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

	void connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
	void evaluate(const step::Step &step, const step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet);
	void updateSolution(Vector_Base<double> *x);
	void nextIteration(Vector_Base<double> *x);

	StructuralMechanicsConfiguration &settings;
	StructuralMechanicsLoadStepConfiguration &configuration;

	struct ParametersElements {
		ElementParameter<ndim * enodes * ndim * enodes> stiffness;
//		ElementParameter<ndim * enodes * ndim * enodes> mass;
		ElementParameter<ndim * enodes> rhs, nrhs;

		struct {
			BoundaryParameter<ndim * enodes * ndim * enodes> stiffness;
//			BoundaryParameter<ndim * enodes * ndim * enodes> mass;
			BoundaryParameter<ndim * enodes> rhs;
		} boundary;
	} elements;

	struct Results {
		static NodeData *displacement, *thickness;
		static ElementData *principalStress, *componentStress, *vonMisesStress;
		static ElementData *isPlastized;
	};

protected:
	void run(Action action, size_t interval);
	void run(Action action, size_t region, size_t interval);

	void runPlane(Action action, size_t interval);
	void runAxisymmetric(Action action, size_t interval);
	void runVolume(Action action, size_t interval);
	void runPreprocess(Action action, size_t interval);
	void runBoundary(Action action, size_t region, size_t interval);

	std::vector<StructuralMechanicsSubKernelsList> subkernels;
	std::vector<std::vector<StructuralMechanicsBoundarySubKernelsList> > boundary;
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_STRUCTURALMECHANICS_H_ */
