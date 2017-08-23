
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <vector>

#include "../../mesh/settings/property.h"

namespace espreso {

struct Step;
enum Matrices : int;
enum class MatrixType;
enum class ElementType;
template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Element;
class Mesh;
class Instance;
class EqualityConstraints;
class Solution;
class SparseMatrix;
namespace store { class ResultStore; }

enum class REGULARIZATION;

enum class SumOperation {
	SUM,
	AVERAGE
};

enum class SumRestriction {
	NONE,
	DIRICHLET,
	NON_DIRICHLET
};

struct Physics {

	Physics(const std::string &name, Mesh *mesh, Instance *instance);
	const std::string& name() const { return _name; }

	virtual std::vector<size_t> solutionsIndicesToStore() const =0;
	virtual std::vector<std::pair<ElementType, Property> > propertiesToStore() const =0;

	virtual void prepare() =0;
	virtual void prepareHybridTotalFETIWithCorners() =0;
	virtual void prepareHybridTotalFETIWithKernels() =0;
	virtual void updateMesh(const std::vector<std::vector<eslocal> > &previousDOFMap, const std::vector<std::vector<eslocal> > &previousDomainMap);

	virtual void preprocessData(const Step &step) =0;

	virtual void updateMatrix(const Step &step, Matrices matrices, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution);

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const =0;
	virtual bool isMatrixTimeDependent(const Step &step) const =0;
	virtual bool isMatrixTemperatureDependent(const Step &step) const =0;

	virtual void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processSolution(const Step &step) =0;

	virtual void makeStiffnessMatricesRegular(REGULARIZATION regularization, size_t scSize, bool ortogonalCluster);
	virtual void makeStiffnessMatrixRegular(REGULARIZATION regularization, size_t scSize, size_t domains, bool ortogonalCluster);
	virtual void analyticRegularization(size_t domain, bool ortogonalCluster) =0;

	virtual void assembleB1(const Step &step, bool withRedundantMultipliers, bool withGluing, bool withScaling);
	virtual void updateDirichletInB1(const Step &step, bool withRedundantMultipliers);
	virtual void assembleB0FromCorners();
	virtual void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels);

	virtual double sumSquares(const std::vector<std::vector<double> > &data, SumOperation operation, SumRestriction restriction = SumRestriction::NONE, size_t loadStep = 0) const;

	virtual ~Physics();

	virtual const std::vector<Property>& pointDOFs() const =0;
	virtual const std::vector<Property>& midPointDOFs() const =0;
	virtual const std::vector<Property>& edgeDOFs() const =0;
	virtual const std::vector<Property>& faceDOFs() const =0;
	virtual const std::vector<Property>& elementDOFs() const =0;

	inline const std::vector<size_t>& pointDOFsOffsets() const
	{
		return _nodesDOFsOffsets;
	}

	inline const std::vector<size_t>& midPointDOFsOffsets() const
	{
		return _midNodesDOFsOffsets;
	}

	inline const std::vector<size_t>& edgeDOFsOffsets() const
	{
		return _edgesDOFsOffsets;
	}

	inline const std::vector<size_t>& faceDOFsOffsets() const
	{
		return _facesDOFsOffsets;
	}

	inline const std::vector<size_t>& elementDOFsOffsets() const
	{
		return _elementsDOFsOffsets;
	}

protected:
	virtual void fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs) const;
	virtual void insertElementToDomain(
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M,
			const std::vector<eslocal> &DOFs,
			const DenseMatrix &Ke, const DenseMatrix &Me, const DenseMatrix &Re, const DenseMatrix &fe,
			const Step &step, size_t domain, bool isBoundaryCondition);

	virtual void assembleBoundaryConditions(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);

	std::string _name;
	Mesh *_mesh;
	Instance *_instance;
	EqualityConstraints *_equalityConstraints;

	std::vector<size_t> _nodesDOFsOffsets;
	std::vector<size_t> _midNodesDOFsOffsets;
	std::vector<size_t> _edgesDOFsOffsets;
	std::vector<size_t> _facesDOFsOffsets;
	std::vector<size_t> _elementsDOFsOffsets;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
