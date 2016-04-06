
#include "../constraints/equalityconstraints.h"

using namespace espreso;

#define LAMBDA_FIRST 0
#define LAMBDA_INNER -1
#define LAMBDA_REMOVED -2

static void offsetSum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<size_t*>(out)) += *(static_cast<size_t*>(in));
}

static size_t scanOffsets(size_t &offset)
{
	size_t size = offset;
	if (config::MPIsize == 1) { // because MPI hates all users
		offset = 0;
		return size;
	}

	MPI_Op op;
	MPI_Op_create(offsetSum, 1, &op);
	MPI_Exscan(&size, &offset, sizeof(size_t), MPI_BYTE, op, MPI_COMM_WORLD);

	size = offset + size;
	MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, config::MPIsize - 1, MPI_COMM_WORLD);
	if (config::MPIrank == 0) {
		// MPI really hates all -> it set offset on all processes to zero except on process 0
		offset = 0;
	}

	return size;
}

static size_t countersToOffset(std::vector<size_t> &offsets)
{
	size_t sum = 0;
	for (size_t i = 0; i < offsets.size(); i++) {
		size_t tmp = offsets[i];
		offsets[i] = sum;
		sum += tmp;
	}
	return sum;
}

static size_t gluingPairs(size_t n) {
	size_t sum = 1;
	for(size_t i = 2; i < n; i++) {
		sum += i;
	}
	return sum;
};

static esglobal pairOffset(size_t i, size_t j, size_t size) {
	esglobal index = 0;
	for (esglobal s = 0; s < i; s++) {
		index += size - s - 1;
	}
	index += j - i - 1;
	return index;
};

static bool compareAccordingFirst(const std::pair<esglobal, esglobal> &p1, const std::pair<esglobal, esglobal> &p2)
{
	return p1.first < p2.first;
}

static bool clusterMappingCompare(const std::vector<esglobal> &v1, const std::vector<esglobal> &v2)
{
	return v1[0] < v2[0];
}

Constraints::Constraints(const Mesh &mesh, size_t firstIndex)
: _mesh(mesh), _subdomains(mesh.parts()), _firstIndex(firstIndex), _DOFs(mesh.DOFs())
{
	_neighbours = mesh.neighbours();
	_neighbours.push_back(config::MPIrank);
	std::sort(_neighbours.begin(), _neighbours.end());
}

size_t Dirichlet::assemble(
		std::vector<SparseMatrix> &B1,
		std::vector<std::vector<esglobal> > &B1clustersMap,
		std::vector<std::vector<double> > &values)
{
	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _dirichletSize);

	std::vector<std::vector<std::vector<esglobal> > > dirichlet(threads, std::vector<std::vector<esglobal> >(_subdomains));
	std::vector<std::vector<std::vector<double> > > dirichletValues(threads, std::vector<std::vector<double> >(_subdomains));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			const eslocal index = (_dirichletIndices[i] - _dirichletOffset) / _DOFs;
			const std::vector<eslocal> &neighs = _mesh.subdomainBoundaries()[index];
			for(size_t s = 0; s < neighs.size(); s++) {
				const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(neighs[s]);
				auto it = lower_bound(l2c.begin(), l2c.end(), index);
				if (it != l2c.end() && *it == index) {
					dirichlet[t][neighs[s]].push_back(_DOFs * (it - l2c.begin()) + (_dirichletIndices[i] - _dirichletOffset) % _DOFs + IJVMatrixIndexing);
					dirichletValues[t][neighs[s]].push_back(_dirichletValues[i]);
				}
			}

		}
	}

	ESINFO(PROGRESS3) << "Dirichlet conditions assigned to subdomains";

	std::vector<size_t> subdomainsCounters(_subdomains, 0);
	cilk_for (size_t s = 0; s < _subdomains; s++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[t][s].size();
		}
		subdomainsCounters[s] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t s = 0; s < _subdomains; s++) {
		clusterOffset += subdomainsCounters[s];
		if (subdomainsCounters[s]) {
			subdomainsWithDirichlet.push_back(s);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t dirichletSize = scanOffsets(clusterOffset);

	clusterOffset += _firstIndex;
	cilk_for (size_t s = 0; s < _subdomains; s++) {
		B1[s].rows += dirichletSize;
	}

	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		B1[s].nnz += subdomainsCounters[s];
		B1[s].I_row_indices.reserve(B1[s].nnz);
		B1[s].J_col_indices.reserve(B1[s].nnz);
		B1[s].V_values.resize(B1[s].nnz, 1);
	}


	countersToOffset(subdomainsCounters);
	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (size_t i = 0; i < B1[s].nnz; i++) {
			B1[s].I_row_indices.push_back(clusterOffset + subdomainsCounters[s] + i + IJVMatrixIndexing);
		}
		for (size_t t = 0; t < threads; t++) {
			B1[s].J_col_indices.insert(B1[s].J_col_indices.end(), dirichlet[t][s].begin(), dirichlet[t][s].end());
			values[s].insert(values[s].end(), dirichletValues[t][s].begin(), dirichletValues[t][s].end());
		}
	}

	B1clustersMap.reserve(clusterDirichletSize);
	for (esglobal i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		B1clustersMap.push_back({ i, config::MPIrank });
	}

	if (dirichletSize == 0) {
		ESINFO(ERROR) << "ESPRESO requires some nodes with Dirichlet condition.";
	}

	ESINFO(DETAILS) << "Number of DOFs with Dirichlet in B1 is " << dirichletSize;

	return dirichletSize;
}

size_t Gluing::assembleB1(
		std::vector<SparseMatrix> &B1,
		std::vector<std::vector<esglobal> > &B1clustersMap,
		std::vector<std::vector<double> > &B1duplicity,
		const std::vector<eslocal> &excludes)
{
	std::vector<esglobal> ids;
	std::vector<std::vector<esglobal> > skippedNodes;

	size_t subdomainsGluingSize = subdomainsLambdaCounters(ids, skippedNodes, excludes);
	ESINFO(PROGRESS3) << "B1 lambdas IDs for subdomains gluing computed";

	std::vector<SparseIJVMatrix<esglobal> > gluing;
	std::vector<std::vector<double> > duplicity;
	std::vector<std::vector<std::vector<esglobal> > > clusterMap;

	composeSubdomainGluing(ids, gluing, duplicity, clusterMap);
	ESINFO(PROGRESS3) << "Composed B1 for subdomains gluing";

	std::vector<std::vector<eslocal> > multiplicity;

	computeNodesMultiplicity(skippedNodes, multiplicity);
	ESINFO(PROGRESS3) << "Computed lambdas multiplicity for cluster gluing";
	size_t clustersGluingSize = clustersLambdaCounters(ids, multiplicity, subdomainsGluingSize);
	ESINFO(PROGRESS3) << "Computed cluster lambdas size";
	exchangeGlobalIds(ids, multiplicity);
	ESINFO(PROGRESS3) << "Lambdas IDs exchanged";
	composeClustersGluing(ids, multiplicity, gluing, duplicity, clusterMap);
	ESINFO(PROGRESS3) << "Composed B1 for clusters gluing";

	cilk_for (size_t s = 0; s < _subdomains; s++) {
		B1[s].rows += subdomainsGluingSize + clustersGluingSize;
		B1[s].nnz += gluing[s].rowIndices().size();
		B1[s].I_row_indices.reserve(B1[s].nnz);
		B1[s].J_col_indices.reserve(B1[s].nnz);
		B1[s].V_values.reserve(B1[s].nnz);
		B1duplicity[s].reserve(B1[s].nnz);

		B1[s].I_row_indices.insert(B1[s].I_row_indices.end(), gluing[s].rowIndices().begin(), gluing[s].rowIndices().end());
		B1[s].J_col_indices.insert(B1[s].J_col_indices.end(), gluing[s].columnIndices().begin(), gluing[s].columnIndices().end());
		B1[s].V_values.insert(B1[s].V_values.end(), gluing[s].values().begin(), gluing[s].values().end());
		B1duplicity[s].insert(B1duplicity[s].end(), duplicity[s].begin(), duplicity[s].end());
	}
	ESINFO(PROGRESS3) << "B1 post-process finished";

	size_t previousSize = B1clustersMap.size(), clusterGluingSize = 0;
	for (size_t s = 0; s < _subdomains; s++) {
		clusterGluingSize += clusterMap[s].size();
	}
	B1clustersMap.reserve(previousSize + clusterGluingSize);
	for (size_t s = 0; s < _subdomains; s++) {
		B1clustersMap.insert(B1clustersMap.end(), clusterMap[s].begin(), clusterMap[s].end());
	}
	std::sort(B1clustersMap.begin() + previousSize, B1clustersMap.end(), clusterMappingCompare);

	ESINFO(DETAILS) << "Total number of lambdas in B1(including Dirichlet) is " << subdomainsGluingSize + clustersGluingSize;

	return subdomainsGluingSize + clusterGluingSize;
}

size_t Gluing::assembleB0(std::vector<SparseMatrix> &B0)
{
	std::vector<SparseIJVMatrix<esglobal> > gluing;

	std::vector<eslocal> corners;
	for (size_t i = 0; i < _mesh.subdomainBoundaries().size(); i++) {
		if (_mesh.subdomainBoundaries().isCorner(i)) {
			corners.push_back(i);
		}
	}

	size_t cornersGluingSize = composeCornersGluing(corners, gluing);

	cilk_for (size_t s = 0; s < _subdomains; s++) {
		B0[s].rows += cornersGluingSize;
		B0[s].nnz += gluing[s].rowIndices().size();
		B0[s].I_row_indices.reserve(B0[s].nnz);
		B0[s].J_col_indices.reserve(B0[s].nnz);
		B0[s].V_values.reserve(B0[s].nnz);

		B0[s].I_row_indices.insert(B0[s].I_row_indices.end(), gluing[s].rowIndices().begin(), gluing[s].rowIndices().end());
		B0[s].J_col_indices.insert(B0[s].J_col_indices.end(), gluing[s].columnIndices().begin(), gluing[s].columnIndices().end());
		B0[s].V_values.insert(B0[s].V_values.end(), gluing[s].values().begin(), gluing[s].values().end());
	}

	ESINFO(DETAILS) << "Total number of lambdas in B0 is " << Info::averageValue(cornersGluingSize);
	return cornersGluingSize;
}

size_t Gluing::subdomainsLambdaCounters(std::vector<esglobal> &ids, std::vector<std::vector<esglobal> > &skippedNodes, const std::vector<eslocal> &excludes)
{
	ids.clear();
	ids.resize(_DOFs * _sBoundary.size(), 0);

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _sBoundary.size());

	skippedNodes.resize(threads);
	std::vector<size_t> offsets(threads, 0);

	// compute nodes on sub-domains boundaries and chunks offsets
	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t offset = 0, exclude_pointer = std::lower_bound(excludes.begin(), excludes.end(), distribution[t]) - excludes.begin();
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {


			if (_sBoundary[n].size() > 1) { // node 'n' belongs to more sub-domains
				if (excludes.size() && excludes[exclude_pointer] == n) {
					exclude_pointer++;
					continue;
				}
				if (_cBoundary[n].size() > 1) {
					skippedNodes[t].push_back(n);
					continue; // skip nodes belong to more clusters -> it will be glued later
				}

				size_t pairs = gluingPairs(_sBoundary[n].size());
				for (size_t i = 0; i < _DOFs; i++) {
					ids[n * _DOFs + i] = pairs;
					offset += pairs;
				}
			}


		}
		offsets[t] = offset;
	}

	return lambdaCountersToIds(ids, distribution, offsets, 0);
}

size_t Gluing::clustersLambdaCounters(std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity, size_t localGluingSize)
{
	ids.clear();
	ids.resize(_cBoundary.size() * _DOFs, 0);

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _cBoundary.size());

	std::vector<size_t> offsets(threads, 0);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t offset = 0;
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {


			if (_cBoundary[n].size() > 1 && config::MPIrank == _cBoundary[n][0]) {
				size_t pairs = gluingPairs(multiplicity[n].back());
				for (size_t i = 0; i < _DOFs; i++) {
					ids[n * _DOFs + i] = pairs;
					offset += pairs;
				}
			}


		}
		offsets[t] = offset;
	}

	return lambdaCountersToIds(ids, distribution, offsets, localGluingSize);
}

size_t Gluing::lambdaCountersToIds(std::vector<esglobal> &ids, std::vector<size_t> &distribution, std::vector<size_t> &offsets, size_t offset)
{
	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> ignoredDistribution = Esutils::getDistribution(threads, _ignoredDOFsSize);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = ignoredDistribution[t]; i < ignoredDistribution[t + 1]; i++) {

			if (ids[_ignoredDOFs[i] - _ignoredDOFsOffset]) {
				ids[_ignoredDOFs[i] - _ignoredDOFsOffset] *= -1;
			}

		}
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = _DOFs * distribution[t]; i < _DOFs * distribution[t + 1]; i++) {

			if (ids[i] < 0) {
				offsets[t] += ids[i];
				ids[i] = LAMBDA_REMOVED;
			}

		}
	}

	size_t subdomainGluingOffset = countersToOffset(offsets);
	size_t subdomainsGluingSize = scanOffsets(subdomainGluingOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += subdomainGluingOffset + offset + _firstIndex;
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t offset = offsets[t];
		for (size_t i = distribution[t] * _DOFs; i < distribution[t + 1] * _DOFs; i++) {
			if (ids[i] > 0) {
				offset += ids[i];
				ids[i] = offset - ids[i];
			} else if (ids[i] != LAMBDA_REMOVED) {
				ids[i] = LAMBDA_INNER;
			}
		}
	}

	return subdomainsGluingSize;
}

void Gluing::composeSubdomainGluing(const std::vector<esglobal> &ids, std::vector<SparseIJVMatrix<esglobal> > &gluing, std::vector<std::vector<double> > &duplicity, std::vector<std::vector<std::vector<esglobal> > > &clusterMap)
{
	gluing.resize(_subdomains);
	duplicity.resize(_subdomains);
	clusterMap.resize(_subdomains);

	cilk_for (size_t s = 0; s < _subdomains; s++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(s);
		for (size_t n = 0; n < l2c.size(); n++) {
			for (size_t d = 0; d < _DOFs; d++) {
				if (ids[l2c[n] * _DOFs + d] >= 0) {
					// for each DOF in sub-domain s

					const std::vector<eslocal> &neighs = _sBoundary[l2c[n]];
					size_t index = std::lower_bound(neighs.begin(), neighs.end(), s) - neighs.begin();
					for (size_t i = 0; i < neighs.size(); i++) {
						if (i != index) {
							esglobal offset = i < index ? pairOffset(i, index, neighs.size()) : pairOffset(index, i, neighs.size());
							gluing[s].push(ids[l2c[n] * _DOFs + d] + offset, n * _DOFs + d, index < i ? 1 : -1);
							duplicity[s].push_back(1. / neighs.size());
							if (index < i) {
								clusterMap[s].push_back({ ids[l2c[n] * _DOFs + d] + offset, config::MPIrank });
							}
						}
					}

				}
			}
		}
	}
}

void Gluing::composeClustersGluing(
		const std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity,
		std::vector<SparseIJVMatrix<esglobal> > &gluing, std::vector<std::vector<double> > &duplicity, std::vector<std::vector<std::vector<esglobal> > > &clusterMap)
{
	clusterMap.resize(_subdomains);
	cilk_for (size_t s = 0; s < _subdomains; s++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(s);
		for (size_t n = 0; n < l2c.size(); n++) {

			if (std::all_of(ids.begin() + l2c[n] * _DOFs, ids.begin() + (l2c[n] + 1) * _DOFs, [] (const esglobal &v) { return v < 0; })) {
				continue;
			}
			size_t index = 0, myFirst;
			const std::vector<eslocal> &neighs = _cBoundary[l2c[n]];
			for (size_t i = 0; i < neighs.size(); i++) {
				if (neighs[i] != config::MPIrank) {
					index += multiplicity[l2c[n]][i];
				} else {
					myFirst = index;
					index += std::lower_bound(_sBoundary[l2c[n]].begin(), _sBoundary[l2c[n]].end(), s) - _sBoundary[l2c[n]].begin();
					break;
				}
			}

			for (size_t d = 0; d < _DOFs; d++) {
				if (ids[l2c[n] * _DOFs + d] >= 0) {
					// for each DOF in sub-domain s

					size_t i = 0;
					for (size_t r = 0; r < multiplicity[l2c[n]].size() - 1; r++) {
						for (size_t m = 0; m < multiplicity[l2c[n]][r]; m++, i++) {
							if (i != index) {
								esglobal offset = i < index ? pairOffset(i, index, multiplicity[l2c[n]].back()) : pairOffset(index, i, multiplicity[l2c[n]].back());
								gluing[s].push(ids[l2c[n] * _DOFs + d] + offset, n * _DOFs + d, index < i ? 1 : -1);
								duplicity[s].push_back(1. / multiplicity[l2c[n]].back());
								if (_cBoundary[l2c[n]][r] == config::MPIrank) {
									if (index < i) {
										clusterMap[s].push_back({ ids[l2c[n] * _DOFs + d] + offset, config::MPIrank });
									}
								} else {
									clusterMap[s].push_back({ ids[l2c[n] * _DOFs + d] + offset, config::MPIrank, _cBoundary[l2c[n]][r] });
								}
							}
						}
					}
				}
			}
		}
	}
}


void Gluing::computeNodesMultiplicity(const std::vector<std::vector<esglobal> > &skippedNodes, std::vector<std::vector<eslocal> > &multiplicity)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_neighbours.begin(), _neighbours.end(), neighbour) - _neighbours.begin();
	};

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _cBoundary.size());

	std::vector<std::vector<std::pair<esglobal, esglobal> > > sBuffer(_neighbours.size());
	std::vector<std::vector<std::pair<esglobal, esglobal> > > rBuffer(_neighbours.size());

	for (size_t t = 0; t < skippedNodes.size(); t++) {
		for (size_t i = 0; i < skippedNodes[t].size(); i++) {
			for (size_t c = 0; c < _cBoundary[skippedNodes[t][i]].size(); c++) {
				sBuffer[n2i(_cBoundary[skippedNodes[t][i]][c])].push_back({ _c2g[skippedNodes[t][i]], _sBoundary[skippedNodes[t][i]].size() });
			}
		}
	}

	cilk_for (size_t n = 0; n < _neighbours.size(); n++) {
		std::sort(sBuffer[n].begin(), sBuffer[n].end(), compareAccordingFirst);
	}

	std::vector<MPI_Request> req(_neighbours.size());

	std::sort(_neighbours.begin(), _neighbours.end());
	for (size_t n = 0; n < _neighbours.size(); n++) {
		MPI_Isend(sBuffer[n].data(), 2 * sizeof(esglobal) * sBuffer[n].size(), MPI_BYTE, _neighbours[n], 0, MPI_COMM_WORLD, req.data() + n);
	}

	int flag, counter = 0;
	MPI_Status status;
	while (counter < _neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / (2 * sizeof(esglobal)));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(_neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	multiplicity.resize(_cBoundary.size());
	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			multiplicity[n].reserve(_cBoundary[n].size());
			size_t size = 0;
			for (size_t i = 0; i < _cBoundary[n].size(); i++) {

				std::pair<esglobal, esglobal> node(_c2g[n], 0);
				auto it = std::lower_bound(rBuffer[n2i(_cBoundary[n][i])].begin(), rBuffer[n2i(_cBoundary[n][i])].end(), node, compareAccordingFirst);
				if (it != rBuffer[n2i(_cBoundary[n][i])].end() && it->first == _c2g[n]) {
					multiplicity[n].push_back(it->second);
				} else {
					multiplicity[n].push_back(1);
				}
				size += multiplicity[n].back();
			}
			multiplicity[n].push_back(size);
		}
	}
}

void Gluing::exchangeGlobalIds(std::vector<esglobal> &ids, const std::vector<std::vector<eslocal> > &multiplicity)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_neighbours.begin(), _neighbours.end(), neighbour) - _neighbours.begin();
	};

	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _cBoundary.size());

	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > nodesIds(threads, std::vector<std::vector<std::pair<esglobal, esglobal> > >(_neighbours.size()));

	for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			if (_cBoundary[n].size() > 1) {
				if (config::MPIrank == _cBoundary[n][0]) {
					// send my ids to higher ranks
					for (size_t i = 1; i < _cBoundary[n].size(); i++) {
						nodesIds[t][n2i(_cBoundary[n][i])].push_back({ _c2g[n], ids[_DOFs * n] });
					}
				} else {
					// get id from the smallest rank
					nodesIds[t][n2i(_cBoundary[n][0])].push_back({ 0, 0 }); // buffer for ids from neighbours
				}
			}
		}
	}

	cilk_for (size_t n = 0; n < _neighbours.size(); n++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += nodesIds[t][n].size();
		}
		if (_neighbours[n] < config::MPIrank) {
			nodesIds[0][n].resize(size);
		} else {
			nodesIds[0][n].reserve(size);
			for (size_t t = 1; t < threads; t++) {
				nodesIds[0][n].insert(nodesIds[0][n].end(), nodesIds[t][n].begin(), nodesIds[t][n].end());
			}
			std::sort(nodesIds[0][n].begin(), nodesIds[0][n].end(), compareAccordingFirst);
		}
	}

	std::vector<MPI_Request> req(_neighbours.size() - 1); // not send to my id
	size_t nCounter = 0;
	for (size_t n = 0; n < _neighbours.size(); n++) {
		if (_neighbours[n] > config::MPIrank) {
			MPI_Isend(nodesIds[0][n].data(), 2 * sizeof(esglobal) * nodesIds[0][n].size(), MPI_BYTE, _neighbours[n], 1, MPI_COMM_WORLD, req.data() + nCounter++);
		}
		if (_neighbours[n] < config::MPIrank) {
			MPI_Irecv(nodesIds[0][n].data(), 2 * sizeof(esglobal) * nodesIds[0][n].size(), MPI_BYTE, _neighbours[n], 1, MPI_COMM_WORLD, req.data() + nCounter++);
		}
	}

	MPI_Waitall(_neighbours.size() - 1, req.data(), MPI_STATUSES_IGNORE);

	std::vector<size_t> sizes(threads, 0);
	for (size_t t = 0; t < threads; t++) {
		size_t size = 0;
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
			if (_cBoundary[n].size() > 1 && _cBoundary[n][0] < config::MPIrank) {
				size_t pairs = gluingPairs(multiplicity[n].back());
				std::pair<esglobal, esglobal> node(_c2g[n], 0);
				auto it = std::lower_bound(nodesIds[0][n2i(_cBoundary[n][0])].begin(), nodesIds[0][n2i(_cBoundary[n][0])].end(), node, compareAccordingFirst);
				if (it->second == LAMBDA_REMOVED) {
					continue; // skip dirichlet
				}
				for (size_t d = 0; d < _DOFs; d++) {
					size += pairs;
					ids[_DOFs * n + d] = it->second + d * pairs;
				}
			}
		}
		sizes[t] += size;
	}
}

size_t Gluing::composeCornersGluing(const std::vector<eslocal> &corners, std::vector<SparseIJVMatrix<esglobal> > &gluing)
{
	gluing.resize(_subdomains);
	if (!corners.size()) {
		return 0;
	}

	auto s2i = [&] (eslocal n, eslocal s) {
		return std::lower_bound(_sBoundary[n].begin(), _sBoundary[n].end(), s) - _sBoundary[n].begin();
	};

	std::vector<eslocal> lambdas;
	lambdas.reserve(corners.size());

	lambdas.push_back(0);
	for (size_t l = 1; l < corners.size(); l++) {
		lambdas.push_back(lambdas.back() + _DOFs * (_sBoundary[corners[l - 1]].size() - 1));
	}

	cilk_for (size_t s = 0; s < _subdomains; s++) {
		const std::vector<eslocal> &l2c = _mesh.coordinates().localToCluster(s);
		for (size_t l = 0; l < lambdas.size(); l++) {

			if (!std::binary_search(_sBoundary[corners[l]].begin(), _sBoundary[corners[l]].end(), s)) {
				continue;
			}

			auto index = std::lower_bound(l2c.begin(), l2c.end(), corners[l]) - l2c.begin();
			if (s > _sBoundary[corners[l]].front()) {
				for (size_t d = 0; d < _DOFs; d++) {
					gluing[s].push(lambdas[l] + _DOFs * (s2i(corners[l], s) - 1) + d, _DOFs * index + d, -1);
				}
			}
			if (s < _sBoundary[corners[l]].back()) {
				for (size_t d = 0; d < _DOFs; d++) {
					gluing[s].push(lambdas[l] + _DOFs * s2i(corners[l], s) + d, _DOFs * index + d, 1);
				}
			}
		}
	}

	return lambdas.back() + _DOFs * (_sBoundary[corners.back()].size() - 1);
}





























