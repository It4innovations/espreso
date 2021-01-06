
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
		size_t interval;

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

	virtual bool boundaryWithSettings(size_t rindex) =0;

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

	struct ParameterInfo;
	struct ParameterInfoLink {
		ParameterInfo *parameter;
		int offset;

		ParameterInfoLink(ParameterInfo &parameter): parameter(&parameter), offset(0) {}
		ParameterInfoLink(ParameterInfo &parameter, int offset, int increment): parameter(&parameter), offset(offset) {}
	};

	struct ParameterInfo {
		enum class Range {
			AGREGATED,
			PER_NODE,
			PER_GP,
			PER_NODE_GP,
			PER_GP_GP,
			EACH_NODE,
			EACH_GP,
			EACH_NODE_GP
		};

		enum class Status {
			GLOBAL,
			PER_INTERVAL,
			EMPTY,
			EXPRESSION
		};

		Range range;
		Status constness;
		int dimensions;

		serializededata<esint, double>* data;

		std::vector<int> isset, isconstant, nodes, gps;
		std::vector<double> ivalues; // in the case of constant parameter (per interval)
		std::vector<const Evaluator*> evaluator;

		std::vector<std::vector<ParameterInfoLink> > inputs;

		ParameterInfo(size_t intervals, int dimensions, ParameterInfo::Range range, double initial = 0);

		void forceResize();
		void smartResize();
		void eval();
		void swapInput(ParameterInfo &toreplace, ParameterInfo &replacement);

		bool constantInterval(esint interval)
		{
			if (isconstant.size()) {
				for (int d = 0; d < dimensions; ++d) {
					if (!isconstant[dimensions * interval + d]) {
						return false;
					}
				}
			}
			return true;
		}

		void addGeneralInput(const ParameterInfoLink &input)
		{
			for (size_t i = 0; i < inputs.size(); ++i) {
				inputs[i].push_back(input);
			}
		}

		void addGeneralInputs() {}

		template <typename...Rest>
		void addGeneralInputs(const ParameterInfoLink &input, Rest&&... rest) { addGeneralInput(input); addGeneralInputs(rest...); }

		void addInput(const ParameterInfoLink &input, esint interval, int indimension = 0, int outdimension = 0)
		{
			inputs[dimensions * interval + indimension].push_back(input);
		}

		void insertParameterValues(Evaluator::Params &params, esint interval, int dimension = 0);
	};

	struct ElementParameterInfo: public ParameterInfo {
		ElementParameterInfo(int dimensions, ParameterInfo::Range range, double initial = 0);
	};

	struct BoundaryParameterInfo: public ParameterInfo {
		BoundaryParameterInfo(int dimensions, ParameterInfo::Range range);
	};

	struct InputParameterIterator {
		const int inc;
		const double * __restrict data;

		InputParameterIterator(const double * data, int increment): inc(increment), data(data) {}
		InputParameterIterator(ParameterInfo &info, esint interval, int increment);

		inline InputParameterIterator& operator++() { data += inc; return *this; }
		inline const double& operator[](esint i) const { return data[i]; }
	};

	struct OuputParameterIterator {
		const int inc;
		double * __restrict data;

		OuputParameterIterator(double * data, int increment): inc(increment), data(data) {}
		OuputParameterIterator(ParameterInfo &info, esint interval, int increment);

		inline OuputParameterIterator& operator++() { data += inc; return *this; }
		inline double& operator[](esint i) { return data[i]; }
		inline const double& operator[](esint i) const { return data[i]; }
	};

protected:
	std::vector<ParameterInfo> eparam;

	void move(ElementParameterInfo &info1, ElementParameterInfo &info2);
	void move(NodeData *data, ElementParameterInfo &info);
	void move(ElementParameterInfo &info, NodeData *data);
	void move(BoundaryParameterInfo &info, NodeData *data);

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
