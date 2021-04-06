
#include <cstring>
#include "data.synchronization.h"
#include "matrix.csr.distributed.h"
#include "vector.dense.distributed.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "basis/utilities/communication.h"

using namespace espreso;

void DataSynchronization::uniformCombination(const DataSynchronization *first, const DataSynchronization *second, int nfirst, int nsecond)
{
	esint sum = nfirst + nsecond;

	sIndices.resize(first->sIndices.size());
	rIndices.resize(first->rIndices.size());
	sBuffer.resize(first->sBuffer.size());
	rBuffer.resize(first->rBuffer.size());
	for (size_t n = 0; n < first->neighbors.size(); ++n) {
		for (size_t i = 0; i < first->sIndices[n].size(); ++i) {
			sIndices[n].push_back(sum * (first->sIndices[n][i] / nfirst) + first->sIndices[n][i] % nfirst);
		}
		for (size_t i = 0; i < first->rIndices[n].size(); ++i) {
			rIndices[n].push_back(sum * (first->rIndices[n][i] / nfirst) + first->rIndices[n][i] % nfirst);
		}
		for (size_t i = 0; i < second->sIndices[n].size(); ++i) {
			sIndices[n].push_back(sum * (second->sIndices[n][i] / nsecond) + second->sIndices[n][i] % nsecond + nfirst);
		}
		for (size_t i = 0; i < second->rIndices[n].size(); ++i) {
			rIndices[n].push_back(sum * (second->rIndices[n][i] / nsecond) + second->rIndices[n][i] % nsecond + nfirst);
		}

		std::sort(sIndices[n].begin(), sIndices[n].end());
		std::sort(rIndices[n].begin(), rIndices[n].end());

		sBuffer[n].resize(first->sBuffer[n].size() + second->sBuffer[n].size());
		rBuffer[n].resize(first->rBuffer[n].size() + second->rBuffer[n].size());
	}

	neighbors = first->neighbors;
}

void DataSynchronization::init(const MatrixCSRDistributed *matrix)
{
	sBuffer.resize(neighbors.size());
	rBuffer.resize(neighbors.size());
	rIndices.resize(neighbors.size());
	sIndices.resize(neighbors.size());

	std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

	for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
		for (esint r = matrix->nintervals[2 * n]; r < matrix->nintervals[2 * n + 1]; ++r) {
			sbuffer[n].push_back(matrix->halo[r]);
			sbuffer[n].push_back(matrix->rows[r + 1] - matrix->rows[r]);
			for (esint c = matrix->rows[r]; c < matrix->rows[r + 1]; ++c) {
				sbuffer[n].push_back(matrix->cols[c - matrix->rows[0]]);
			}
		}
		sBuffer[n].resize(matrix->rows[matrix->nintervals[2 * n + 1]] - matrix->rows[matrix->nintervals[2 * n]]);
	}

	if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, neighbors)) {
		eslog::error("ESPRESO internal error: receive MatrixCSRDistribution pattern.\n");
	}

	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rbuffer[n].size(); ) {
			esint *c = matrix->cols + matrix->rows[rbuffer[n][i++] - matrix->distribution[info::mpi::rank] + matrix->nhalo] - 1;
			esint columns = rbuffer[n][i++];
			for (esint cc = 0; cc < columns; ++cc, ++i) {
				while (*c != rbuffer[n][i]) { ++c; }
				rIndices[n].push_back(c - matrix->cols);
			}
		}
		rBuffer[n].resize(rIndices[n].size());
	}
}

void DataSynchronization::init(const VectorDenseDistributed *vector)
{
	sBuffer.resize(neighbors.size());
	rBuffer.resize(neighbors.size());
	rIndices.resize(neighbors.size());
	sIndices.resize(neighbors.size());

	std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

	for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
		for (esint r = vector->nintervals[2 * n]; r < vector->nintervals[2 * n + 1]; ++r) {
			sbuffer[n].push_back(vector->halo[r]);
		}
		sBuffer[n].resize(vector->nintervals[2 * n + 1] - vector->nintervals[2 * n]);
	}

	if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, neighbors)) {
		eslog::error("ESPRESO internal error: receive MatrixDenseDistributed pattern.\n");
	}

	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rbuffer[n].size(); ++i) {
			rIndices[n].push_back(rbuffer[n][i] - vector->distribution[info::mpi::rank] + vector->nhalo);
		}
		rBuffer[n].resize(rIndices[n].size());
	}
}

void DataSynchronization::init(const VectorSparseDistributed *vector)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
}

void DataSynchronization::gatherFromUpper(const MatrixCSRDistributed *matrix)
{
	for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
		memcpy(sBuffer[n].data(), matrix->vals + matrix->rows[matrix->nintervals[2 * n]] - matrix->rows[0], sizeof(double) * sBuffer[n].size());
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, neighbors)) {
		eslog::error("ESPRESO internal error: receive MatrixCSRDistribution data.\n");
	}

	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rIndices[n].size(); ++i) {
			matrix->vals[rIndices[n][i]] += rBuffer[n][i];
		}
	}
}

void DataSynchronization::gatherFromUpper(const VectorDenseDistributed *vector)
{
	for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
		memcpy(sBuffer[n].data(), vector->vals + vector->nintervals[2 * n], sizeof(double) * sBuffer[n].size());
	}

	if (!Communication::receiveUpperKnownSize(sBuffer, rBuffer, neighbors)) {
		eslog::error("ESPRESO internal error: receive VectorDenseDistributed data.\n");
	}

	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rIndices[n].size(); ++i) {
			vector->vals[rIndices[n][i]] += rBuffer[n][i];
		}
	}
}

void DataSynchronization::gatherFromUpper(const VectorSparseDistributed *vector)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
}

void DataSynchronization::scatterToUpper(const MatrixCSRDistributed *matrix)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
}

void DataSynchronization::scatterToUpper(const VectorDenseDistributed *vector)
{
	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rIndices[n].size(); ++i) {
			rBuffer[n][i] = vector->vals[rIndices[n][i]];
		}
	}

	if (!Communication::receiveLowerKnownSize(rBuffer, sBuffer, neighbors)) {
		eslog::error("ESPRESO internal error: scatter VectorDenseDistributed data.\n");
	}

	for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
		memcpy(vector->vals + vector->nintervals[2 * n], sBuffer[n].data(), sizeof(double) * sBuffer[n].size());
	}
}

void DataSynchronization::scatterToUpper(const VectorSparseDistributed *vector)
{
	eslog::error("ESPRESO internal error: call empty function.\n");
}
