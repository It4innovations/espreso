
#ifndef SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_

#include "math2/primitives/matrix_info.h"

namespace espreso {

class Assembler {

public:
	Matrix_Type matrixType() { return Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE; }

	void init() {}

};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_ASSEMBLER_H_ */
