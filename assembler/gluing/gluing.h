
#ifndef ASSEMBLER_GLUING_GLUING_H_
#define ASSEMBLER_GLUING_GLUING_H_

#include "../assembler.h"

namespace assembler {

// rename to equality constrains
template <class TInput>
class Gluing: public Assembler<TInput> {

public:
	virtual ~Gluing() {};

protected:
	Gluing(TInput &input);

	// computeDirichlet
	// computeGluing
	// computeRigidBodyModes

	void computeSubdomainGluing();
	void computeClusterGluing(std::vector<size_t> &rows);

	virtual size_t DOFs() = 0;

	// Matrices for Linear Solver
	// B0 B1
	std::vector<SparseMatrix> _localB, _globalB;

	// Description ??
	std::vector<std::vector<eslocal> > _lambda_map_sub_B1;
	std::vector<std::vector<eslocal> > _lambda_map_sub_B0;
	std::vector<std::vector<eslocal> > _lambda_map_sub_clst;
	std::vector<std::vector<double> > _B1_duplicity;
	std::vector<std::vector<double> > _vec_c;
	std::vector<eslocal> _neighClusters;

private:
	// used for computation
	std::vector<SparseIJVMatrix<eslocal> > _B0, _B1;
};

}

#include "gluing.hpp"


#endif /* ASSEMBLER_GLUING_GLUING_H_ */
