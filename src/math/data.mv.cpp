
#include <cstring>
#include "data.mv.h"
#include "matrix.indices.h"
#include "matrix.csr.distributed.h"
#include "vector.dense.distributed.h"
#include "math.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"

#include <algorithm>

using namespace espreso;

void DataMV::uniformCombination(const DataMV *first, const DataMV *second, int nfirst, int nsecond)
{
	esint sum = nfirst + nsecond;

	send.resize(first->send.size());
	recv.resize(first->recv.size());
	sBuffer.resize(first->sBuffer.size());
	rBuffer.resize(first->rBuffer.size());
	for (size_t n = 0; n < first->neighbors.size(); ++n) {
		for (size_t i = 0; i < first->send[n].size(); ++i) {
			send[n].push_back(sum * (first->send[n][i] / nfirst) + first->send[n][i] % nfirst);
		}
		for (size_t i = 0; i < first->recv[n].size(); ++i) {
			recv[n].push_back(sum * (first->recv[n][i] / nfirst) + first->recv[n][i] % nfirst);
		}
		for (size_t i = 0; i < second->send[n].size(); ++i) {
			send[n].push_back(sum * (second->send[n][i] / nsecond) + second->send[n][i] % nsecond + nfirst);
		}
		for (size_t i = 0; i < second->recv[n].size(); ++i) {
			recv[n].push_back(sum * (second->recv[n][i] / nsecond) + second->recv[n][i] % nsecond + nfirst);
		}

		std::sort(send[n].begin(), send[n].end());
		std::sort(recv[n].begin(), recv[n].end());

		sBuffer[n].resize(first->sBuffer[n].size() + second->sBuffer[n].size());
		rBuffer[n].resize(first->rBuffer[n].size() + second->rBuffer[n].size());
	}

	neighbors = first->neighbors;

	if (first->minCol <= second->minCol) {
		minCol = sum * (first->minCol / nfirst) + first->minCol % nfirst;
	} else {
		minCol = sum * (second->minCol / nsecond) + second->minCol % nsecond + nfirst;
	}

	m.uniformCombination(&first->m, &second->m, nfirst, nsecond);
	v.uniformCombination(&first->v, &second->v, nfirst, nsecond);
}

void DataMV::init(const MatrixCSRDistributed *matrix)
{
	std::vector<IJ> synchronization;

	auto minmax = std::minmax_element(matrix->cols, matrix->cols + matrix->nnz);

	minCol = *minmax.first;
	v.resize(*minmax.second - *minmax.first + 1);
	m.resize(matrix->nrows - matrix->nhalo, v.size, matrix->rows[matrix->nrows] - matrix->rows[matrix->nhalo]);

	m.rows[0] = 1;
	for (esint r = matrix->nhalo; r < matrix->nrows; ++r) {
		for (esint c = matrix->rows[r]; c < matrix->rows[r + 1]; ++c) {
			m.cols[c - matrix->rows[matrix->nhalo]] = matrix->cols[c - 1] - *minmax.first + 1;
			esint cindex = matrix->cols[c - 1] - 1;
			if (cindex < matrix->distribution[info::mpi::rank] || matrix->distribution[info::mpi::rank + 1] <= cindex) {
				synchronization.push_back({r, matrix->cols[c - 1] - 1});
			}
		}
		m.rows[r - matrix->nhalo + 1] = matrix->rows[r + 1] - matrix->rows[matrix->nhalo] + 1;
	}

	std::sort(synchronization.begin(), synchronization.end(), [] (const IJ &i, const IJ &j) {
		if (i.column == j.column) {
			return i.row < j.row;
		}
		return i.column < j.column;
	});

	neighbors.assign(matrix->neighbors, matrix->neighbors + matrix->nneighbors);
	auto sbegin = synchronization.begin();
	while (sbegin != synchronization.end()) {
		int neighbor = std::lower_bound(matrix->distribution, matrix->distribution + info::mpi::size + 1, sbegin->column + 1) - matrix->distribution - 1;
		if (neighbor != info::mpi::rank) {
			neighbors.push_back(neighbor);
		}
		sbegin = std::lower_bound(sbegin, synchronization.end(), matrix->distribution[neighbor + 1], [] (const IJ &index, esint value) {
			return index.column < value;
		});
	}
	utils::sortAndRemoveDuplicates(neighbors);

	sBuffer.resize(neighbors.size());
	rBuffer.resize(neighbors.size());
	send.resize(neighbors.size());
	recv.resize(neighbors.size());

	for (size_t n = 0, i = 0; n < neighbors.size(); n++) {
		while (i < synchronization.size() && synchronization[i].column < matrix->distribution[neighbors[n] + 1]) {
			send[n].push_back(synchronization[i].row);
			recv[n].push_back(synchronization[i].column - *minmax.first + 1);
			++i;
		}
	}

	for (size_t n = 0; n < neighbors.size(); n++) {
		utils::sortAndRemoveDuplicates(send[n]);
		utils::sortAndRemoveDuplicates(recv[n]);
		sBuffer[n].resize(send[n].size());
		rBuffer[n].resize(recv[n].size());
	}

	m.structureUpdated();
}

void DataMV::apply(const MatrixCSRDistributed *matrix, const VectorDenseDistributed *in, VectorDenseDistributed *out)
{
	memcpy(m.vals, matrix->vals + matrix->rows[matrix->nhalo] - 1, sizeof(double) * m.nnz);
	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < send[n].size(); ++i) {
			sBuffer[n][i] = in->vals[send[n][i]];
		}
	}

	if (!Communication::exchangeKnownSize(sBuffer, rBuffer, neighbors)) {
		eslog::internalFailure("exchange MV data.\n");
	}

	std::copy(in->vals + in->nhalo, in->vals + in->size, v.vals + in->distribution[info::mpi::rank] - minCol + 1);
	for (size_t n = 0; n < neighbors.size(); ++n) {
		for (size_t i = 0; i < rBuffer[n].size(); ++i) {
			v.vals[recv[n][i]] = rBuffer[n][i];
		}
	}

	out->vals += out->nhalo;
	m.apply(&v, out);
	out->vals -= out->nhalo;
	out->scatterToUpper();
}
