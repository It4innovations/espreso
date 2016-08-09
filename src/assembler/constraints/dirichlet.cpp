
#include "dirichlet.h"

using namespace espreso;


void Dirichlet::insertDirichletToB1(const std::vector<Element*> &nodes, const Coordinates &coordinates, const std::vector<Property> &DOFs)
{
	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

	// part x thread x Dirichlet
	std::vector<std::vector<std::vector<esglobal> > > dirichlet(_mesh.parts(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<double> > > dirichletValues(_mesh.parts(), std::vector<std::vector<double> >(threads));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (nodes[i]->settings().isSet(DOFs[dof])) {
					const Point &p = coordinates[i];
					const std::vector<eslocal>& indices = nodes[i]->DOFsIndices();
					double value = nodes[i]->settings(DOFs[dof]).back()->evaluate(p.x, p.y, p.z);
					for(size_t d = 0; d < nodes[i]->domains().size(); d++) {
						if (indices[d * DOFs.size() + dof] != -1) {
							dirichlet[nodes[i]->domains()[d]][t].push_back(indices[d * DOFs.size() + dof] + IJVMatrixIndexing);
							dirichletValues[nodes[i]->domains()[d]][t].push_back(value);
						}
					}
				}
			}

		}
	}

	std::vector<size_t> dirichletSizes(_mesh.parts());
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[p][t].size();
		}
		dirichletSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t p = 0; p < _mesh.parts(); p++) {
		clusterOffset += dirichletSizes[p];
		if (dirichletSizes[p]) {
			subdomainsWithDirichlet.push_back(p);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t globalDirichletSize = synchronizeOffsets(clusterOffset);

	clusterOffset += B1[0].rows;
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		B1[p].rows += globalDirichletSize;
	}

	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		B1[s].nnz += dirichletSizes[s];
		B1[s].I_row_indices.reserve(B1[s].nnz);
		B1[s].J_col_indices.reserve(B1[s].nnz);
		B1[s].V_values.resize(B1[s].nnz, 1);
	}


	Esutils::sizesToOffsets(dirichletSizes);
	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (size_t i = 0; i < B1[s].nnz; i++) {
			B1[s].I_row_indices.push_back(clusterOffset + dirichletSizes[s] + i + IJVMatrixIndexing);
		}
		for (size_t t = 0; t < threads; t++) {
			B1[s].J_col_indices.insert(B1[s].J_col_indices.end(), dirichlet[s][t].begin(), dirichlet[s][t].end());
			B1c[s].insert(B1c[s].end(), dirichletValues[s][t].begin(), dirichletValues[s][t].end());
		}
		B1duplicity[s].resize(B1[s].I_row_indices.size(), 1);
	}

	B1clustersMap.reserve(B1clustersMap.size() + clusterDirichletSize);
	for (esglobal i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		B1clustersMap.push_back({ i, config::env::MPIrank });
	}

	if (globalDirichletSize == 0) {
		ESINFO(ERROR) << "ESPRESO requires some nodes with Dirichlet condition.";
	}

	ESINFO(DETAILS) << "Lambdas with Dirichlet in B1: " << B1[0].rows;
}
