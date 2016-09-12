
#include "assembler.h"

using namespace espreso;

void Physics::assembleScalingMatrices()
{
	D.resize(K.size());
	cilk_for (size_t p = 0; p < K.size(); p++) {
		D[p] = K[p].getDiagonal();
	}

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _mesh.nodes().size());

	// threads x neighbours x (nodes x DOFs)
	std::vector<std::vector<std::vector<double> > > dofBuffer(threads, std::vector<std::vector<double> >(_mesh.neighbours().size()));
	std::vector<std::vector<std::vector<esglobal> > > nBuffer(threads, std::vector<std::vector<esglobal> >(_mesh.neighbours().size()));

	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t c = 0; _mesh.nodes()[i]->clusters().size() > 1 && c < _mesh.nodes()[i]->clusters().size(); c++) {
				if (_mesh.nodes()[i]->clusters()[c] != config::env::MPIrank) {
					nBuffer[t][n2i(_mesh.nodes()[i]->clusters()[c])].push_back(_mesh.coordinates().globalIndex(i));
					for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
						double k = 0;
						for (size_t d = 0; d < _mesh.nodes()[i]->domains().size(); d++) {
							if (_mesh.nodes()[i]->DOFIndex(_mesh.nodes()[i]->domains()[d], dof) != -1) {
								k += D[_mesh.nodes()[i]->domains()[d]][_mesh.nodes()[i]->DOFIndex(_mesh.nodes()[i]->domains()[d], dof)];
							}
						}
						dofBuffer[t][n2i(_mesh.nodes()[i]->clusters()[c])].push_back(k);
					}
				}
			}

		}
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
			dofBuffer[0][n].insert(dofBuffer[0][n].end(), dofBuffer[t][n].begin(), dofBuffer[t][n].end());
			nBuffer[0][n].insert(nBuffer[0][n].end(), nBuffer[t][n].begin(), nBuffer[t][n].end());
		}
	}

	std::vector<std::vector<double> > sBuffer(_mesh.neighbours().size());
	std::vector<std::vector<double> > rBuffer(_mesh.neighbours().size());
	std::vector<std::vector<esglobal> > commonNodes(_mesh.neighbours().size());

	cilk_for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		std::vector<eslocal> permutation(nBuffer[0][n].size());
		std::iota(permutation.begin(), permutation.end(), 0);
		std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) { return nBuffer[0][n][i] < nBuffer[0][n][j]; });
		commonNodes[n].reserve(permutation.size());
		for (size_t i = 0; i < permutation.size(); i++) {
			commonNodes[n].push_back(nBuffer[0][n][permutation[i]]);
		}
		rBuffer[n].resize(permutation.size() * pointDOFs.size());
		sBuffer[n].reserve(permutation.size() * pointDOFs.size());
		for (size_t i = 0; i < permutation.size(); i++) {
			sBuffer[n].insert(sBuffer[n].end(), dofBuffer[0][n].begin() + pointDOFs.size() * permutation[i], dofBuffer[0][n].begin() + pointDOFs.size() * (permutation[i] + 1));
		}
	}

	std::vector<MPI_Request> req(2 * _mesh.neighbours().size());
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		MPI_Isend(sBuffer[n].data(), sBuffer[n].size() * sizeof(double), MPI_BYTE, _mesh.neighbours()[n], 0, MPI_COMM_WORLD, req.data() + 2 * n);
		MPI_Irecv(rBuffer[n].data(), rBuffer[n].size() * sizeof(double), MPI_BYTE, _mesh.neighbours()[n], 0, MPI_COMM_WORLD, req.data() + 2 * n + 1);
	}
	MPI_Waitall(2 * _mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

	cilk_for (size_t p = 0; p < K.size(); p++) {
		for (size_t i = 0; i < _mesh.coordinates().localSize(p); i++) {
			size_t n = _mesh.coordinates().clusterIndex(i, p);

			if (_mesh.nodes()[n]->clusters().size() == 1 && _mesh.nodes()[n]->domains().size() == 1) {
				continue;
			}

			for (size_t dof = 0; dof < pointDOFs.size(); dof++) {
				double k = 0;
				for (size_t d = 0; d < _mesh.nodes()[n]->domains().size(); d++) {
					k += D[_mesh.nodes()[n]->domains()[d]][_mesh.nodes()[n]->DOFIndex(_mesh.nodes()[n]->domains()[d], dof)];
				}
				for (size_t c = 0; c < _mesh.nodes()[n]->clusters().size(); c++) {
					if (_mesh.nodes()[i]->clusters()[c] != config::env::MPIrank) {
						size_t neigh = n2i(_mesh.nodes()[n]->clusters()[c]);
						size_t index = std::lower_bound(commonNodes[neigh].begin(), commonNodes[neigh].end(), _mesh.coordinates().globalIndex(n)) - commonNodes[neigh].begin();
						k += rBuffer[neigh][pointDOFs.size() * index + dof];
					}
				}
				D[p][_mesh.nodes()[n]->DOFIndex(p, dof)] /= k;
			}
		}
	}
}




