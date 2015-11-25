
#include "wrapper.h"

using namespace assembler;

std::list<FETI4IStructRHS*> DataHolder::RHSs;
std::list<FETI4IStructMatrix*> DataHolder::matrices;
std::list<FETI4IStructIntance*> DataHolder::instances;

using namespace assembler;

int FETI4ICreateStiffnessMatrix(
		FETI4IMatrix *stiffnessMatrix,
		FETI4IInt n,
		FETI4IInt nelt,
		FETI4IInt* eltptr,
		FETI4IInt* eltvar,
		FETI4IReal* values)
{
	FETI4IInt indexing = eltptr[0];

	SparseVVPMatrix<eslocal> matrix(n, n);

	FETI4IInt value = 0;
	for (FETI4IInt e = 0; e < nelt; e++) {
		for (FETI4IInt i = eltptr[e] - indexing; i < eltptr[e + 1] - indexing; i++) {
			FETI4IInt size = eltptr[e + 1] - eltptr[e];
			for (FETI4IInt j = 0; j < size; j++) {
				matrix(eltvar[eltptr[e] - indexing + j] - indexing, eltvar[i] - indexing) = values[value++];
			}
		}
	}

	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexing));
	DataHolder::matrices.back()->data = matrix;
	*stiffnessMatrix = DataHolder::matrices.back();
	return 0;
}

int FETI4ICreateStiffnessMatrixAndRHS(
		FETI4IMatrix 	*stiffnessMatrix,
		FETI4IRHS 		*rhs,
		FETI4IInt		indexBase)
{
	DataHolder::RHSs.push_back(new FETI4IStructRHS(indexBase));
	DataHolder::matrices.push_back(new FETI4IStructMatrix(indexBase));
	*stiffnessMatrix = DataHolder::matrices.back();
	*rhs = DataHolder::RHSs.back();
	return 0;
}

int FETI4IAddElement(
		FETI4IMatrix 	stiffnessMatrix,
		FETI4IRHS 		rhs,
		FETI4IInt 		size,
		FETI4IInt* 		indices,
		FETI4IReal* 	eMatrix,
		FETI4IReal*		eRHS)
{
	eslocal offset = stiffnessMatrix->offset;
	for (eslocal i = 0; i < size; i++) {
		for (eslocal j = 0; j < size; j++) {
			stiffnessMatrix->data(indices[i] - offset,  indices[j] - offset) = eMatrix[i * size + j];
		}
	}
	offset = rhs->offset;
	for (eslocal i = 0; i < size; i++) {
		rhs->data(0, indices[i] - offset) = eRHS[i];
	}
	return 0;
}

int FETI4ICreateInstance(
		FETI4IInstance 	*instance,
		FETI4IMatrix 	stiffnessMatrix,
		FETI4IRHS 		rhs,
		FETI4IInt 		dirichlet_size,
		FETI4IInt* 		dirichlet_indices,
		FETI4IReal* 	dirichlet_values,
		FETI4IInt* 		l2g,
		FETI4IInt 		neighbours_size,
		FETI4IInt* 		neighbours)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);
  
  std::cout.setstate(std::ios_base::failbit);
	API2 api;
	DataHolder::instances.push_back(new FETI4IStructIntance(api));
	DataHolder::instances.back()->K = stiffnessMatrix->data;
	api.K = &(DataHolder::instances.back()->K);
	SparseIJVMatrix<eslocal> ijv = rhs->data;
	DenseMatrix d = ijv;
	DataHolder::instances.back()->rhs = std::vector<double>(d.values(), d.values() + d.columns());
	api.rhs = &(DataHolder::instances.back()->rhs);
	api.dirichlet_size = dirichlet_size;
	api.dirichlet_indices = dirichlet_indices;
	api.dirichlet_values = dirichlet_values;
	api.l2g_size = api.K->rows();
	api.l2g = l2g;
	api.neighbours_size = neighbours_size;
	api.neighbours = neighbours;
	DataHolder::instances.back()->data = assembler::LinearElasticity<assembler::API2>(api);

	DataHolder::instances.back()->data.init();
	*instance = DataHolder::instances.back();
	return 0;
}

int FETI4ISolve(
	FETI4IInstance instance,
	FETI4IInt solution_size,
	FETI4IReal* solution)
{
	std::vector<std::vector<double> > solutions(1);
	solutions[0] = std::vector<double>(solution, solution + solution_size);
	instance->data.solve(solutions);
	memcpy(solution, &solutions[0][0], solution_size);
	return 0;
}

template <typename TFETI4I>
static void destroy(std::list<TFETI4I*> &list, void *value)
{
	for (typename std::list<TFETI4I*>::iterator it = list.begin(); it != list.end(); ++it) {
		if (*it == value) {
			delete *it;
			list.erase(it);
			return;
		}
	}
}

int FETI4IDestroy(void *data)
{
	destroy(DataHolder::RHSs, data);
	destroy(DataHolder::matrices, data);
	destroy(DataHolder::instances, data);
	return 0;
}


