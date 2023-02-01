
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_INFO_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_INFO_H_

#include "analysis/assembler/operator.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct DataDescriptorElementSize: ActionOperator, Physics {
	DataDescriptorElementSize(size_t interval, std::vector<size_t> &esize)
	{
		esize[interval] = (sizeof(typename Physics::Element));
	}

	void sisd(typename Physics::Element &element) {}
	void simd(typename Physics::Element &element) {}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_INFO_H_ */
