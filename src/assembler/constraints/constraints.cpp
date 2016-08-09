
#include "constraints.h"

using namespace espreso;


ConstraintsBase::ConstraintsBase(Mesh &mesh, Physics &physics)
: _mesh(mesh), _physics(physics)
{
	B0.resize(_mesh.parts());
	B0subdomainsMap.resize(_mesh.parts());

	B1.resize(_mesh.parts());
	B1subdomainsMap.resize(_mesh.parts());
	B1duplicity.resize(_mesh.parts());
	B1c.resize(_mesh.parts());

	for (size_t p = 0; p < _mesh.parts(); p++) {
		B0[p].rows = 0;
		B0[p].cols = _mesh.coordinates().localSize(p) * _physics.pointDOFs.size();
		B0[p].nnz = 0;
		B0[p].type = 'G';

		B1[p].rows = 0;
		B1[p].cols = _mesh.coordinates().localSize(p) * _physics.pointDOFs.size();
		B1[p].nnz = 0;
		B1[p].type = 'G';
	}
}

static void offsetSum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<size_t*>(out)) += *(static_cast<size_t*>(in));
}

size_t ConstraintsBase::synchronizeOffsets(size_t &offset)
{
	size_t size = offset;
	if (config::env::MPIsize == 1) {
		offset = 0;
		return size;
	}

	MPI_Op op;
	MPI_Op_create(offsetSum, 1, &op);
	MPI_Exscan(&size, &offset, sizeof(size_t), MPI_BYTE, op, MPI_COMM_WORLD);

	size = offset + size;
	MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, config::env::MPIsize - 1, MPI_COMM_WORLD);
	if (config::env::MPIrank == 0) {
		offset = 0;
	}

	return size;
}


