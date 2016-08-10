
#include "equalitygluing.h"

using namespace espreso;

void EqualityGluing::insertDomainGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
	std::vector<size_t> offsets(threads);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->domains().size() > 1 && elements[e]->clusters().size() < 2) {
				for (size_t i = 0; i < DOFs.size(); i++) {
					if (!elements[e]->settings().isSet(DOFs[i])) { // Dirichlet is not glued
						size_t n = elements[e]->numberOfLocalDomainsWithDOF(i);
						lambdasSize += n * (n - 1) / 2;
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t clusterOffset = Esutils::sizesToOffsets(offsets);
	size_t clusterGluingSize = clusterOffset;
	size_t globalSize = synchronizeOffsets(clusterOffset);
	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + B1[0].rows;
	}

	// threads x domains x values
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<double> > > duplicity(threads, std::vector<std::vector<double> >(_mesh.parts()));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			if (elements[e]->domains().size() > 1 && elements[e]->clusters().size() < 2) {
				const std::vector<eslocal> &DOFIndices = elements[e]->DOFsIndices();
				for (size_t i = 0; i < DOFs.size(); i++) {
					if (elements[e]->settings().isSet(DOFs[i])) { // Dirichlet is not glued
						continue;
					}

					size_t dup = 0, start_offset = offset;
					for (size_t d1 = 0; d1 < elements[e]->domains().size(); d1++) {
						if (DOFIndices[d1 * DOFs.size() + i] != -1) {
							for (size_t d2 = d1 + 1; d2 < elements[e]->domains().size(); d2++) {
								if (DOFIndices[d2 * DOFs.size() + i] != -1) {
									rows[t][elements[e]->domains()[d1]].push_back(offset + IJVMatrixIndexing);
									cols[t][elements[e]->domains()[d1]].push_back(DOFIndices[d1 * DOFs.size() + i] + IJVMatrixIndexing);
									vals[t][elements[e]->domains()[d1]].push_back(1);

									rows[t][elements[e]->domains()[d2]].push_back(offset + IJVMatrixIndexing);
									cols[t][elements[e]->domains()[d2]].push_back(DOFIndices[d2 * DOFs.size() + i] + IJVMatrixIndexing);
									vals[t][elements[e]->domains()[d2]].push_back(-1);
									offset++;
								}
							}
							if (dup == 0) { // set duplicity only once
								dup = offset - start_offset + 1;
							}
							if (dup) { // if there was some lambda, fill dup and map
								duplicity[t][elements[e]->domains()[d1]].insert(duplicity[t][elements[e]->domains()[d1]].end(), dup - 1, 1. / dup);
							}
						}
					}

				}
			}
		}
	}

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t t = 0; t < threads; t++) {
			B1[p].I_row_indices.insert(B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			B1[p].J_col_indices.insert(B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			B1[p].V_values.insert(B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			B1duplicity[p].insert(B1duplicity[p].end(), duplicity[t][p].begin(), duplicity[t][p].end());
		}
		B1[p].rows += globalSize;
		B1[p].nnz = B1[p].I_row_indices.size();
	}

	B1clustersMap.reserve(B1clustersMap.size() + clusterGluingSize);
	for (esglobal i = clusterOffset; i < clusterOffset + clusterGluingSize; i++) {
		B1clustersMap.push_back({ i, config::env::MPIrank });
	}

	ESINFO(DETAILS) << "Lambdas gluing domains in B1: " << B1[0].rows;
}

std::vector<esglobal> EqualityGluing::computeLambdasID(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	std::vector<esglobal> lambdasID(elements.size() * DOFs.size(), -1);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
	std::vector<size_t> offsets(threads);

	// neighbours x threads x data
	std::vector<std::vector<std::vector<esglobal> > > sLambdas(_mesh.neighbours().size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<esglobal> > > rLambdas(_mesh.neighbours().size(), std::vector<std::vector<esglobal> >(threads));


	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				size_t n = elements[e]->numberOfGlobalDomainsWithDOF(dof);
				if (n > 1 && !elements[e]->settings().isSet(DOFs[dof])) { // Dirichlet and inner nodes are not glued
					if (elements[e]->clusters()[0] == config::env::MPIrank) { // set lambda ID
						lambdasID[e * DOFs.size() + dof] = n * (n - 1) / 2;
						lambdasSize += n * (n - 1) / 2;
						for (size_t c = 1; c < elements[e]->clusters().size(); c++) { // send to higher clusters
							sLambdas[n2i(elements[e]->clusters()[c])][t].push_back(e * DOFs.size() + dof);
						}
					} else { // pick ID from lower cluster
						rLambdas[n2i(elements[e]->clusters()[0])][t].push_back(e * DOFs.size() + dof);
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t numberOfClusterLambdas = Esutils::sizesToOffsets(offsets);
	size_t clusterOffset = numberOfClusterLambdas;
	size_t totalNumberOfLambdas = synchronizeOffsets(clusterOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + B1[0].rows;
	}

	cilk_for (size_t t = 0; t < threads; t++) {
		esglobal offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (lambdasID[e * DOFs.size() + dof] > 0) {
					offset += lambdasID[e * DOFs.size() + dof];
					lambdasID[e * DOFs.size() + dof] = offset - lambdasID[e * DOFs.size() + dof];
				}
			}

		}

		for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
			for (size_t i = 0; i < sLambdas[n][t].size(); i++) {
				sLambdas[n][t][i] = lambdasID[sLambdas[n][t][i]];
			}
		}
	}

	std::vector<std::vector<esglobal> > rBuffer(_mesh.neighbours().size());

	cilk_for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		size_t size = rLambdas[n][0].size();
		for (size_t t = 1; t < threads; t++) {
			sLambdas[n][0].insert(sLambdas[n][0].end(), sLambdas[n][t].begin(), sLambdas[n][t].end());
			size += rLambdas[n][t].size();
		}
		rBuffer[n].resize(size);
	}


	std::vector<MPI_Request> req(_mesh.neighbours().size());
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		if (_mesh.neighbours()[n] > config::env::MPIrank) {
			MPI_Isend(sLambdas[n][0].data(), sizeof(esglobal) * sLambdas[n][0].size(), MPI_BYTE, _mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
		if (_mesh.neighbours()[n] < config::env::MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(esglobal) * rBuffer[n].size(), MPI_BYTE, _mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
	}

	MPI_Waitall(_mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
				size_t offset = 0;
				for (size_t i = 0; i < t; i++) {
					offset += rLambdas[n][i].size();
				}
				for (size_t i = 0; i < rLambdas[n][t].size(); i++) {
					lambdasID[rLambdas[n][t][i]] = rBuffer[n][offset + i];
				}
			}

		}
	}

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		B1[p].rows += totalNumberOfLambdas;
	}

	return lambdasID;
}

void EqualityGluing::insertClusterGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	std::vector<esglobal> lambdasID = computeLambdasID(elements, DOFs);

	std::vector<eslocal> permutation(lambdasID.size());
	std::iota(permutation.begin(), permutation.end(), 0);

	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return (lambdasID[i] < 0) ? false : (lambdasID[j] < 0) ? true : lambdasID[i] < lambdasID[j];
	});

	auto it = std::find_if(permutation.begin(), permutation.end(), [&] (eslocal i) { return lambdasID[i] == -1; });
	permutation.resize(it - permutation.begin());


	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, permutation.size());

	// threads x domains x data
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(_mesh.parts()));

	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);

	auto findDomain = [&] (const Element *e, size_t d, size_t dof) {
		auto &DOFIndices = e->DOFsIndices();
		size_t c = 0, DOFs = DOFIndices.size() / e->domains().size();
		for (size_t i = 0; i < e->domains().size(); i++) {
			if (DOFIndices[i * DOFs + dof] != -1) {
				if (d == c++) {
					return e->domains()[i];
				}
			}
		}
		return 0;
	};

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();
			esglobal offset = 0;
			double duplicity = 1. / e->numberOfGlobalDomainsWithDOF(dof);

			for (auto c1 = e->clusters().begin(); c1 != e->clusters().end(); ++c1) {
				for (size_t d1 = 0; d1 < e->DOFCounter(*c1, dof); d1++) {

					for (auto c2 = c1; c2 != e->clusters().end(); ++c2) {
						for (size_t d2 = (*c1 == *c2 ? d1 + 1 : 0); d2 < e->DOFCounter(*c2, dof); d2++) {

							if (*c1 == config::env::MPIrank) {
								eslocal d = findDomain(e, d1, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + IJVMatrixIndexing);
								cols[t][d].push_back(e->DOFIndex(d, dof) + IJVMatrixIndexing);
								vals[t][d].push_back(1);
								dup[t][d].push_back(duplicity);
							}

							if (*c2 == config::env::MPIrank) {
								eslocal d = findDomain(e, d2, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + IJVMatrixIndexing);
								cols[t][d].push_back(e->DOFIndex(d, dof) + IJVMatrixIndexing);
								vals[t][d].push_back(-1);
								dup[t][d].push_back(duplicity);
							}

							if (*c1 == config::env::MPIrank || *c2 == config::env::MPIrank) {
								cMap[t].push_back({ lambdasID[permutation[i]] + offset });
								if (*c1 == *c2) {
									cMap[t].back().push_back(*c1);
								} else {
									cMap[t].back().push_back(*c1);
									cMap[t].back().push_back(*c2);
								}
							}

							offset++;

						}
					}

				}

			}

		}
	}



	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t t = 0; t < threads; t++) {
			B1[p].I_row_indices.insert(B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			B1[p].J_col_indices.insert(B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			B1[p].V_values.insert(B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			B1duplicity[p].insert(B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
		}
		B1[p].nnz = B1[p].I_row_indices.size();
	}

	for (size_t t = 0; t < threads; t++) {
		B1clustersMap.insert(B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
	}

	ESINFO(DETAILS) << "Lambdas gluing domains in B1: " << B1[0].rows;
}

void EqualityGluing::insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{

}



