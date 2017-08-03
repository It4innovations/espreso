
#ifndef SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_
#define SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_

#include "physics.h"

namespace espreso {

class APIMesh;

struct Precomputed: public virtual Physics
{
	enum SolutionIndex {
		DECOMPOSED = 0,
		MERGED     = 1,

		SIZE       = 2
	};

	Precomputed(Mesh *mesh, Instance *instance, MatrixType type, double *rhs, size_t rhsSize);

	std::vector<size_t> solutionsIndicesToStore() const;
	std::vector<std::pair<ElementType, Property> > propertiesToStore() const;

	MatrixType getMatrixType(const Step &step, size_t domain) const;
	bool isMatrixTimeDependent(const Step &step) const;
	bool isMatrixTemperatureDependent(const Step &step) const;

	void prepare();
	void prepareHybridTotalFETIWithCorners();
	void prepareHybridTotalFETIWithKernels();
	void preprocessData(const Step &step);

	void analyticRegularization(size_t domain);

	void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	void assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling);
	void assembleB0FromCorners();
	void assembleB0FromKernels(const std::vector<SparseMatrix> &kernels);

	void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const;
	void processSolution(const Step &step);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::UNKNOWN };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::UNKNOWN };
		return midPointDOFs;
	}
	const std::vector<Property>& edgeDOFs() const
	{
		static std::vector<Property> edgeDOFs = { };
		return edgeDOFs;
	}
	const std::vector<Property>& faceDOFs() const
	{
		static std::vector<Property> faceDOFs = { };
		return faceDOFs;
	}
	const std::vector<Property>& elementDOFs() const
	{
		static std::vector<Property> elementDOFs = { };
		return elementDOFs;
	}

	virtual ~Precomputed() {}

protected:
	MatrixType _mtype;
	double *_rhs;
	size_t _rhsSize;
	std::vector<std::vector<double> > _mergedSolution;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_PRECOMPUTED_H_ */
