
#ifndef SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_

#include "esmesh.h"
#include "../../solver/generic/SparseMatrix.h"

namespace espreso {

class Constraints
{
public:
	// matrices for Hybrid FETI constraints
	std::vector<SparseMatrix> B0;
	std::vector<std::vector<esglobal> > B0subdomainsMap; // TODO: not needed

	// matrices for FETI constraints
	std::vector<SparseMatrix> B1;
	std::vector<std::vector<esglobal> > B1subdomainsMap; // TODO: not needed
	std::vector<std::vector<esglobal> > B1clustersMap; // TODO: get it directly

	std::vector<std::vector<double> > B1c;
	std::vector<std::vector<double> > B1duplicity;


	void initMatrices(const std::vector<size_t> &columns);
	void save();

	virtual void insertDirichletToB1(const std::vector<Element*> &nodes, const std::vector<Property> &DOFs) =0;
	virtual void insertElementGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;
	virtual void insertMortarGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;

	virtual void insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;
	virtual void insertKernelsToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel) = 0;

	virtual ~Constraints() {};

protected:
	Constraints(Mesh &mesh): _mesh(mesh) {};

	size_t synchronizeOffsets(size_t &offset);

	Mesh &_mesh;
};

}


#endif /* SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_ */
