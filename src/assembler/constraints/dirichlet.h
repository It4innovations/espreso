
#ifndef SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_
#define SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_

#include "constraints.h"

namespace espreso {

//class Dirichlet: public ConstraintsBase
//{
//public:
//	Dirichlet(const Mesh &mesh, std::vector<Property>& DOFs, size_t offset, const std::vector<eslocal> &indices, const std::vector<double> &values)
//	:Constraints(mesh, DOFs, offset),
//	 _dirichletSize(indices.size()), _dirichletOffset(0),
//	 _dirichletIndices(indices.data()), _dirichletValues(values.data()) { };
//
//	Dirichlet(const Mesh &mesh, std::vector<Property> &DOFs, size_t firstIndex, size_t size, eslocal offset, const eslocal *indices, const double *values)
//	:Constraints(mesh, DOFs, firstIndex),
//	 _dirichletSize(size), _dirichletOffset(offset),
//	_dirichletIndices(indices), _dirichletValues(values) { };
//
//	size_t assemble(
//			std::vector<SparseMatrix> &B1,
//			std::vector<std::vector<esglobal> > &B1clustersMap,
//			std::vector<std::vector<double> > &values
//	);
//
//protected:
//	const size_t _dirichletSize;
//	const eslocal _dirichletOffset;
//	const eslocal *_dirichletIndices;
//	const double *_dirichletValues;
//};

}

#endif /* SRC_ASSEMBLER_CONSTRAINTS_DIRICHLET_H_ */
