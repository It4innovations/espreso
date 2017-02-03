
#ifndef SRC_ASSEMBLER_PHYSICS_PHYSICS_H_
#define SRC_ASSEMBLER_PHYSICS_PHYSICS_H_

#include <cstddef>
#include <vector>

#include "../../mesh/settings/property.h"

namespace espreso {

template<typename TIndices> class SparseVVPMatrix;
class DenseMatrix;
class Element;
class Mesh;
class NewInstance;
namespace store { class ResultStore; }

enum class REGULARIZATION;

struct NewPhysics {

	NewPhysics(Mesh *mesh, NewInstance *instance);

	virtual void prepareTotalFETI() =0;
	virtual void prepareHybridTotalFETIWithCorners() =0;
	virtual void prepareHybridTotalFETIWithKernels() =0;

	void assembleStiffnessMatrices();
	virtual void assembleStiffnessMatrix(size_t domain) =0;
	virtual void assembleStiffnessMatrix(Element *e, std::vector<eslocal> &DOFs, DenseMatrix &Ke, DenseMatrix &fe) =0;

	virtual void processElement(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processFace(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processEdge(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;
	virtual void processNode(const Element *e, DenseMatrix &Ke, DenseMatrix &fe) const =0;

	virtual void fillDOFsIndices(const Element *e, eslocal domain, std::vector<eslocal> &DOFs);
	virtual void insertElementToDomain(SparseVVPMatrix<eslocal> &K, const std::vector<eslocal> &DOFs, const DenseMatrix &Ke, const DenseMatrix &fe, size_t domain);

	void makeStiffnessMatricesRegular(REGULARIZATION regularization);
	virtual void analyticRegularization(size_t domain) =0;

	virtual void assembleB1(bool withRedundantMultipliers, bool withScaling);
	virtual void assembleB0FromCorners() =0;
	virtual void assembleB0FromKernels() =0;

	virtual void storeSolution(std::vector<std::vector<double> > &solution, store::ResultStore *store) =0;

	virtual ~NewPhysics() {}

	virtual const std::vector<Property>& pointDOFs() const =0;
	virtual const std::vector<Property>& midPointDOFs() const =0;
	virtual const std::vector<Property>& edgeDOFs() const =0;
	virtual const std::vector<Property>& faceDOFs() const =0;
	virtual const std::vector<Property>& elementDOFs() const =0;

protected:
	Mesh *_mesh;
	NewInstance *_instance;
};

}

#endif /* SRC_ASSEMBLER_PHYSICS_PHYSICS_H_ */
