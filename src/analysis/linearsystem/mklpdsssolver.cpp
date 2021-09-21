
#include "mklpdsssystem.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

using namespace espreso;

struct __Dirichlet__ {
	double value;
	esint row, column;
};

void espreso::setDirichlet(Matrix_Distributed<Matrix_CSR, double> &A, Vector_Distributed<Vector_Dense, double> &b, const Vector_Sparse<double> &dirichlet, const DOFsDistribution &distribution)
{
	std::vector<__Dirichlet__> tosend;

	esint nhalo = distribution.halo.size();
	auto getrow = [&] (esint dof) {
		esint row = -1;
		if (distribution.begin <= dof && dof < distribution.end) {
			row = dof - distribution.begin + nhalo;
		} else {
			auto inhalo = std::lower_bound(distribution.halo.begin(), distribution.halo.end(), dof);
			if (inhalo != distribution.halo.end() && *inhalo == dof) {
				row = inhalo - distribution.halo.begin();
			}
		}
		return row;
	};

	for (esint i = 0; i < dirichlet.nnz; ++i) {
		b.cluster.vals[dirichlet.indices[i]] = dirichlet.vals[i];
		esint col = 0;
		if (dirichlet.indices[i] < nhalo) {
			col = distribution.halo[dirichlet.indices[i]] + 1;
		} else {
			col = distribution.begin + dirichlet.indices[i] - nhalo + 1;
		}
		for (esint j = A.cluster.rows[dirichlet.indices[i]]; j < A.cluster.rows[dirichlet.indices[i] + 1]; j++) {
			if (A.cluster.cols[j - 1] == col) {
				A.cluster.vals[j - 1] = 1;
			} else {
				A.cluster.vals[j - 1] = 0;
				esint row = getrow(A.cluster.cols[j - 1] - 1);
				if (row != -1 && !std::binary_search(dirichlet.indices, dirichlet.indices + dirichlet.nnz, row)) {
					for (esint c = A.cluster.rows[row]; c < A.cluster.rows[row + 1]; c++) {
						if (A.cluster.cols[c - 1] == col) {
							if (row < nhalo) {
								tosend.push_back({A.cluster.vals[c - 1] * b.cluster.vals[dirichlet.indices[i]], A.cluster.cols[j - 1] - 1, A.cluster.cols[c - 1] - 1});
							}
							b.cluster.vals[row] -= A.cluster.vals[c - 1] * b.cluster.vals[dirichlet.indices[i]];
							A.cluster.vals[c - 1] = 0;
						}
					}
				} else {
					// the column will be updated by a higher process
				}
			}
		}
	}

	std::sort(tosend.begin(), tosend.end(), [] (const __Dirichlet__ &i, const __Dirichlet__ &j) {
		if (i.row == j.row) {
			return i.column < j.column;
		}
		return i.row < j.row;
	});

	std::vector<std::vector<__Dirichlet__> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());
	for (size_t n = 0, i = 0; n < info::mesh->neighbors.size(); n++) {
		while (i < tosend.size() && tosend[i].row < distribution.neighDOF[n + 1]) {
			sBuffer[n].push_back(tosend[i++]);
		}
	}

	if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("synchronize dirichlet data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
		for (size_t i = 0; i < rBuffer[n].size(); i++) {
			esint row = getrow(rBuffer[n][i].column);
			if (row == -1) {
				row = getrow(rBuffer[n][i].row);
				b.cluster.vals[row] -= rBuffer[n][i].value;
				esint c = std::lower_bound(A.cluster.cols + A.cluster.rows[row] - 1, A.cluster.cols + A.cluster.rows[row + 1] - 1, rBuffer[n][i].column + 1) - A.cluster.cols;
				A.cluster.vals[c] = 0;
			}
		}
	}
}


void espreso::setDirichlet(Matrix_Distributed<Matrix_CSR, std::complex<double>> &A, Vector_Distributed<Vector_Dense, std::complex<double>> &b, const Vector_Sparse<std::complex<double>> &dirichlet, const DOFsDistribution &distribution)
{
}
