
#include "constraints.h"

using namespace espreso;


void Constraints::initMatrices(const std::vector<size_t> &columns)
{
	B0.resize(_mesh.parts());
	B0subdomainsMap.resize(_mesh.parts());

	B1.resize(_mesh.parts());
	B1subdomainsMap.resize(_mesh.parts());
	B1duplicity.resize(_mesh.parts());
	B1c.resize(_mesh.parts());
	LB.resize(_mesh.parts());

	for (size_t p = 0; p < _mesh.parts(); p++) {
		B0[p].rows = 0;
		B0[p].cols = columns[p];
		B0[p].nnz = 0;
		B0[p].type = 'G';

		B1[p].rows = 0;
		B1[p].cols = columns[p];
		B1[p].nnz = 0;
		B1[p].type = 'G';
	}
}

void Constraints::save()
{
	ESINFO(PROGRESS2) << "Save matrices B0 and B1";
	for (size_t p = 0; p < _mesh.parts(); p++) {
		std::ofstream osK(Logging::prepareFile(p, "B0").c_str());
		osK << B0[p];
		osK.close();

		std::ofstream osF(Logging::prepareFile(p, "B1").c_str());
		osF << B1[p];
		osF.close();

		std::ofstream osd(Logging::prepareFile(p, "weight").c_str());
		osd << B1duplicity[p];
		osd.close();

		std::ofstream osdi(Logging::prepareFile(p, "loc_ind_weight").c_str());
		osdi << B1subdomainsMap[p];
		osdi.close();

		std::ofstream osc(Logging::prepareFile(p, "c1").c_str());
		osc << B1c[p];
		osc.close();
	}
}

static void offsetSum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<size_t*>(out)) += *(static_cast<size_t*>(in));
}

size_t Constraints::synchronizeOffsets(size_t &offset)
{
	size_t size = offset;
	if (environment->MPIsize == 1) {
		offset = 0;
		return size;
	}

	MPI_Op op;
	MPI_Op_create(offsetSum, 1, &op);
	MPI_Exscan(&size, &offset, sizeof(size_t), MPI_BYTE, op, MPI_COMM_WORLD);

	size = offset + size;
	MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, environment->MPIsize - 1, MPI_COMM_WORLD);
	if (environment->MPIrank == 0) {
		offset = 0;
	}

	return size;
}


