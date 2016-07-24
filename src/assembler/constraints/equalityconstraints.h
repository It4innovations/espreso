
#ifndef ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_
#define ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_

#include "../assembler.h"

namespace espreso {

class Constraints
{
protected:
	Constraints(const Mesh &mesh, std::vector<DOFType> &DOFs, size_t firstIndex);

	const Mesh &_mesh;
	std::vector<int> _neighbours;

	size_t _subdomains;
	size_t _firstIndex;
	std::vector<DOFType> _DOFs;
};

class Dirichlet: public Constraints
{
public:
	Dirichlet(const Mesh &mesh, std::vector<DOFType>& DOFs, size_t offset, const std::vector<eslocal> &indices, const std::vector<double> &values)
	:Constraints(mesh, DOFs, offset),
	 _dirichletSize(indices.size()), _dirichletOffset(0),
	 _dirichletIndices(indices.data()), _dirichletValues(values.data()) { };

	Dirichlet(const Mesh &mesh, std::vector<DOFType> &DOFs, size_t firstIndex, size_t size, eslocal offset, const eslocal *indices, const double *values)
	:Constraints(mesh, DOFs, firstIndex),
	 _dirichletSize(size), _dirichletOffset(offset),
	_dirichletIndices(indices), _dirichletValues(values) { };

	size_t assemble(
			std::vector<SparseMatrix> &B1,
			std::vector<std::vector<esglobal> > &B1clustersMap,
			std::vector<std::vector<double> > &values
	);

protected:
	const size_t _dirichletSize;
	const eslocal _dirichletOffset;
	const eslocal *_dirichletIndices;
	const double *_dirichletValues;
};

class Gluing: public Constraints
{
public:
	Gluing(const Mesh &mesh, std::vector<DOFType> &DOFs, size_t firstIndex, const std::vector<eslocal> &ignoredDOFs)
	:Constraints(mesh, DOFs, firstIndex),
	  _ignoredDOFsSize(ignoredDOFs.size()), _ignoredDOFsOffset(0), _ignoredDOFs(ignoredDOFs.data()),
	  _sBoundary(mesh.subdomainBoundaries()), _cBoundary(mesh.clusterBoundaries()),
	  _c2g(mesh.coordinates().clusterToGlobal()) { };

	Gluing(const Mesh &mesh, std::vector<DOFType> &DOFs, size_t firstIndex, size_t ignoredDOFsSize, eslocal ignoredDOFsOffset, const eslocal *ignoredDOFs)
	:Constraints(mesh, DOFs, firstIndex),
	 _ignoredDOFsSize(ignoredDOFsSize), _ignoredDOFsOffset(ignoredDOFsOffset), _ignoredDOFs(ignoredDOFs),
	 _sBoundary(mesh.subdomainBoundaries()), _cBoundary(mesh.clusterBoundaries()),
	 _c2g(mesh.coordinates().clusterToGlobal()) { };

	Gluing(const Mesh &mesh, std::vector<DOFType> &DOFs)
	:Constraints(mesh, DOFs, 0),
	  _ignoredDOFsSize(0), _ignoredDOFsOffset(0), _ignoredDOFs(NULL),
	  _sBoundary(mesh.subdomainBoundaries()), _cBoundary(mesh.clusterBoundaries()),
	  _c2g(mesh.coordinates().clusterToGlobal()) { };

	size_t assembleB1(
			std::vector<SparseMatrix> &B1,
			std::vector<std::vector<esglobal> > &B1clustersMap,
			std::vector<std::vector<double> > &B1duplicity,
			const std::vector<eslocal> &excludes
	);

	size_t assembleB0(
			std::vector<SparseMatrix> &B0
	);

	size_t assembleB0fromKernels(
			std::vector<SparseMatrix> &B0
	);

protected:
	const size_t _ignoredDOFsSize;
	const eslocal _ignoredDOFsOffset;
	const eslocal *_ignoredDOFs;

	const Boundaries &_sBoundary;
	const Boundaries &_cBoundary;
	const std::vector<esglobal> &_c2g;

private:
	// calculating lambdas IDs
	size_t subdomainsLambdaCounters(std::vector<esglobal> &ids, std::vector<std::vector<esglobal> > &skippedNodes, const std::vector<eslocal> &excludes);
	size_t clustersLambdaCounters(std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity, size_t localGluingSize);
	size_t lambdaCountersToIds(std::vector<esglobal> &ids, std::vector<size_t> &distribution, std::vector<size_t> &offsets, size_t offset);

	// exchanges IDs among clusters
	void computeNodesMultiplicity(const std::vector<std::vector<esglobal> > &skippedNodes, std::vector<std::vector<eslocal> > &multiplicity);
	void exchangeGlobalIds(std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity);

	// compose gluing matrices
	void composeSubdomainGluing(const std::vector<esglobal> &ids, std::vector<SparseIJVMatrix<esglobal> > &gluing, std::vector<std::vector<double> > &duplicity, std::vector<std::vector<std::vector<esglobal> > > &clusterMap);
	void composeClustersGluing(const std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity,
			std::vector<SparseIJVMatrix<esglobal> > &gluing, std::vector<std::vector<double> > &duplicity, std::vector<std::vector<std::vector<esglobal> > > &clusterMap);

	size_t composeCornersGluing(const std::vector<eslocal> &corners, std::vector<SparseIJVMatrix<esglobal> > &gluing);
};

template <class TInput>
class EqualityConstraints: public Assembler<TInput> {

	EqualityConstraints(const Mesh &mesh, std::vector<DOFType> &DOFs);
	EqualityConstraints(const Mesh &mesh, std::vector<DOFType> &DOFs, eslocal dirichlet_size, eslocal* dirichlet_indices, double* dirichlet_values, eslocal indexBase);

protected:
	EqualityConstraints(TInput &input);

	void assembleConstraints(std::vector<size_t> columns);

	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> _B0;
	std::vector<std::vector<esglobal> > _B0subdomainsMap; // TODO: not needed

protected:
	const Mesh &_mesh;
	std::vector<DOFType> _DOFs;

	std::vector<std::vector<double> > _B1c;
	std::vector<std::vector<double> > _B1duplicity;
};

}

#include "../constraints/equalityconstraints.hpp"


#endif /* ASSEMBLER_CONSTRAINTS_EQUALITYCONSTRAINTS_H_ */
