
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <vector>

#include "../../mesh/settings/property.h"

namespace espreso {

struct Step;
enum class MatrixType;
template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Element;
class Mesh;
class Instance;
namespace store { class ResultStore; }

enum class REGULARIZATION;

struct NewPhysics {

	NewPhysics(Mesh *mesh, Instance *instance);

	virtual void prepareTotalFETI() =0;
	virtual void prepareHybridTotalFETIWithCorners() =0;
	virtual void prepareHybridTotalFETIWithKernels() =0;

	virtual void assembleStiffnessMatrices(const Step &step);
	virtual void assembleStiffnessMatrix(const Step &step, size_t domain);
	virtual void assembleStiffnessMatrix(const Step &step, Element *e, std::vector<eslocal> &DOFs, DenseMatrix &Ke, DenseMatrix &fe);

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const =0;

	virtual void processElement(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processFace(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processEdge(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processNode(const Step &step, const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;

	virtual void fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs);
	virtual void insertElementToDomain(SparseVVPMatrix<eslocal> &K, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &fe, size_t domain);

	virtual void makeStiffnessMatricesRegular(REGULARIZATION regularization);
	virtual void analyticRegularization(size_t domain) =0;

	virtual void assembleB1(const Step &step, bool withRedundantMultipliers, bool withScaling);
	virtual void assembleB0FromCorners(const Step &step) =0;
	virtual void assembleB0FromKernels(const Step &step) =0;

	virtual void storeSolution(const Step &step, std::vector<std::vector<double> > &solution, store::ResultStore *store) =0;

	virtual ~NewPhysics() {}

	virtual const std::vector<Property>& pointDOFs() const =0;
	virtual const std::vector<Property>& midPointDOFs() const =0;
	virtual const std::vector<Property>& edgeDOFs() const =0;
	virtual const std::vector<Property>& faceDOFs() const =0;
	virtual const std::vector<Property>& elementDOFs() const =0;

protected:
	Mesh *_mesh;
	Instance *_instance;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
