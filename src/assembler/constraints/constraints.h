
#ifndef SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_
#define SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_

#include "esmesh.h"
#include "essolver.h"

namespace espreso {

class ConstraintsBase
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


	virtual void insertDirichletToB1(const std::vector<Element*> &nodes, const Coordinates &coordinates, const std::vector<Property> &DOFs) =0;

	virtual void insertDomainGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;
	virtual void insertClusterGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;

	virtual void insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs) =0;

	virtual ~ConstraintsBase() {};

protected:
	ConstraintsBase(Mesh &mesh, Physics &physics);

	size_t synchronizeOffsets(size_t &offset);

	Mesh &_mesh;
	Physics &_physics;
};

}


#endif /* SRC_ASSEMBLER_CONSTRAINTS_CONSTRAINTS_H_ */
