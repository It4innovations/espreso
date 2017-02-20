
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <vector>

#include "../../mesh/settings/property.h"

namespace espreso {

struct Step;
enum Matrices : int;
enum class MatrixType;
template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Element;
class Mesh;
class Instance;
class Solution;
namespace store { class ResultStore; }

enum class REGULARIZATION;

struct Physics {

	enum class SumOperation {
		SUM,
		AVERAGE
	};

	enum class SumRestriction {
		NONE,
		DIRICHLET
	};

	Physics(const std::string &name, Mesh *mesh, Instance *instance);
	const std::string& name() const { return _name; }

	virtual void prepareTotalFETI() =0;
	virtual void prepareHybridTotalFETIWithCorners() =0;
	virtual void prepareHybridTotalFETIWithKernels() =0;

	virtual void assembleMatrix(const Step &step, Matrices matrices);
	virtual void assembleMatrix(const Step &step, Matrices matrices, size_t domain);
	virtual void assembleMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe);

	virtual void updateMatrix(const Step &step, Matrices matrices, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, size_t domain, const std::vector<Solution*> &solution);
	virtual void updateMatrix(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution);

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const =0;

	virtual void processElement(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processFace(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processEdge(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processNode(const Step &step, Matrices matrices, const Element *e, DenseMatrix &Ke, DenseMatrix &Me, DenseMatrix &Re, DenseMatrix &fe, const std::vector<Solution*> &solution) const =0;
	virtual void processSolution(const Step &step) =0;

	virtual void fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs) const;
	virtual void insertElementToDomain(SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &Me, const DenseMatrix &Re, const DenseMatrix &fe, size_t domain);

	virtual void makeStiffnessMatricesRegular(REGULARIZATION regularization);
	virtual void analyticRegularization(size_t domain) =0;

	virtual void assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling);
	virtual void assembleB0FromCorners(const Step &step) =0;
	virtual void assembleB0FromKernels(const Step &step) =0;

	virtual double sumSquares(const std::vector<std::vector<double> > &data, SumOperation operation, SumRestriction restriction = SumRestriction::NONE, size_t loadStep = 0) const;

	virtual ~Physics() {}

	virtual const std::vector<Property>& pointDOFs() const =0;
	virtual const std::vector<Property>& midPointDOFs() const =0;
	virtual const std::vector<Property>& edgeDOFs() const =0;
	virtual const std::vector<Property>& faceDOFs() const =0;
	virtual const std::vector<Property>& elementDOFs() const =0;

	Instance* instance() { return _instance; }

//protected:
	std::string _name;
	Mesh *_mesh;
	Instance *_instance;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
