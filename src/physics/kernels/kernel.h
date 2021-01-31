
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_KERNEL_H_

#include "basis/evaluator/evaluator.h"
#include "basefunctions/basefunctions.h"
#include "math/math.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cstddef>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <functional>

namespace espreso {

class NodeData;
class ECFExpression;
class Builder;
struct SolverDataProvider;

class Kernel {
public:
	struct InstanceFiller {
		esint begin, end;
		esint interval;

		int invalid;
		int DOFs;
		int insertK, insertM, insertC, insertR, insertF;

		MatrixDense Ke, Me, Ce, CMe;
		VectorsDense Re, Fe;

		std::function<void()> insert;

		double reduction;
		double *K, *M, *C, *R, *F;
		int *offset;

		InstanceFiller(esint nvectors);
	};

	virtual void nextSubstep() =0;
	virtual void solutionChanged() =0;

	virtual void processElements(const Builder &builder, InstanceFiller &filler) =0;
	virtual void processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler) =0;
	virtual void processSolution() =0;

	Kernel(SolverDataProvider *provider): solverDataProvider(provider)
	{
		BaseFunctions::setBaseFunctions();
	}
	virtual ~Kernel() {}

	SolverDataProvider *solverDataProvider;
	std::vector<VectorDense> solutions;
protected:
	void smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order) const;
};


template<class Instance>
class KernelExecutor: public Kernel {
public:
	void nextSubstep();
	void solutionChanged();

	bool boundaryWithSettings(size_t rindex);

	void processElements(const Builder &builder, InstanceFiller &filler);
	void processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler);
	void processSolution();

	KernelExecutor(SolverDataProvider *provider): Kernel(provider) {}
	virtual ~KernelExecutor() {}

private:
	template<typename T> typename std::enable_if<std::is_member_function_pointer<decltype(&T::processNode)>::value>::type
	processNodes(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler);
	template<typename T> void
	processNodes(T * t, const Builder &builder, size_t rindex, const InstanceFiller &filler) {}

	template<typename T> typename std::enable_if<std::is_member_function_pointer<decltype(&T::processEdge)>::value>::type
	processEdges(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler);
	template<typename T> void
	processEdges(T * t, const Builder &builder, size_t rindex, const InstanceFiller &filler) {}

	template<typename T> typename std::enable_if<std::is_member_function_pointer<decltype(&T::processFace)>::value>::type
	processFaces(T * t, const Builder &builder, size_t rindex, InstanceFiller &filler);
	template<typename T> void
	processFaces(T * t, const Builder &builder, size_t rindex, const InstanceFiller &filler) {}

	template<typename T> typename std::enable_if<std::is_member_function_pointer<decltype(&T::nodeSolution)>::value>::type
	nodesSolution(T * t);
	template<typename T> void
	nodesSolution(const T * t) {}

	template<typename T> typename std::enable_if<std::is_member_function_pointer<decltype(&T::elementSolution)>::value>::type
	elementsSolution(T * t);
	template<typename T> void
	elementsSolution(const T * t) {}
};

}

#include "kernel.hpp"

#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_KERNEL_H_ */
