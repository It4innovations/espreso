
#include "equalityconstraints.h"

#include "../instance.h"

#include "../../basis/utilities/utils.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/settings/property.h"

#include "../../config/configuration.h"
#include "../../config/environment.h"

#include <numeric>

using namespace espreso;

static std::vector<std::vector<Region*> > getRegionsWithDOFs(const std::vector<Region*> &regions, size_t loadStep, const std::vector<Property> &DOFs)
{
	std::vector<std::vector<Region*> > result(DOFs.size());
	for (size_t r = 0; r < regions.size(); r++) {
		Region* region = regions[r];
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (loadStep < region->settings.size() && region->settings[loadStep].count(DOFs[dof])) {
				result[dof].push_back(region);
			}
		}
	}
	for (size_t dof = 0; dof < DOFs.size(); dof++) {
		std::sort(result[dof].begin(), result[dof].end());
	}
	return result;
}



static bool commonRegion(const std::vector<Region*> &v1, const std::vector<Region*> &v2) {
	for (size_t i = 0, j = 0; i < v1.size() && j < v2.size(); v1[i] < v2[j] ? i++ : j++) {
		if (v1[i] == v2[j]) {
			return true;
		}
	}
	return false;
};

void EqualityConstraints::insertDirichletToB1(Constraints &constraints, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs)
{
	size_t loadStep = 0;

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

	// part x thread x Dirichlet
	std::vector<std::vector<std::vector<esglobal> > > dirichlet(constraints._mesh.parts(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<double> > > dirichletValues(constraints._mesh.parts(), std::vector<std::vector<double> >(threads));

	std::vector<std::vector<Region*> > fixedRegions = getRegionsWithDOFs(constraints._mesh.regions(), loadStep, DOFs);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (commonRegion(fixedRegions[dof], nodes[i]->regions())) {
					if (!constraints._configuration.redundant_lagrange && nodes[i]->clusters()[0] != environment->MPIrank) {
						continue;
					}
					const std::vector<eslocal>& indices = nodes[i]->DOFsIndices();
					double value = nodes[i]->getProperty(DOFs[dof], 0, loadStep, 0);
					for(size_t d = 0; d < nodes[i]->domains().size(); d++) {
						if (indices[d * DOFs.size() + dof] != -1) {
							dirichlet[nodes[i]->domains()[d]][t].push_back(indices[d * DOFs.size() + dof] + 1);
							dirichletValues[nodes[i]->domains()[d]][t].push_back(value);
						}
						if (!constraints._configuration.redundant_lagrange) {
							break;
						}
					}
				}
			}

		}
	}


	std::vector<size_t> dirichletSizes(constraints._mesh.parts());
	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[p][t].size();
		}
		dirichletSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		clusterOffset += dirichletSizes[p];
		if (dirichletSizes[p]) {
			subdomainsWithDirichlet.push_back(p);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t globalDirichletSize = constraints.synchronizeOffsets(clusterOffset);
	constraints.block[Constraints::BLOCK::DIRICHLET] += globalDirichletSize;

	clusterOffset += constraints.B1[0].rows;
	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		constraints.B1[p].rows += globalDirichletSize;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		constraints.B1[s].nnz += dirichletSizes[s];
		constraints.B1[s].I_row_indices.reserve(constraints.B1[s].nnz);
		constraints.B1[s].J_col_indices.reserve(constraints.B1[s].nnz);
		constraints.B1[s].V_values.resize(constraints.B1[s].nnz, 1);
	}

	Esutils::sizesToOffsets(dirichletSizes);
	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (eslocal i = 0; i < constraints.B1[s].nnz; i++) {
			constraints.B1[s].I_row_indices.push_back(clusterOffset + dirichletSizes[s] + i + 1);
		}
		for (size_t t = 0; t < threads; t++) {
			constraints.B1[s].J_col_indices.insert(constraints.B1[s].J_col_indices.end(), dirichlet[s][t].begin(), dirichlet[s][t].end());
			constraints.B1c[s].insert(constraints.B1c[s].end(), dirichletValues[s][t].begin(), dirichletValues[s][t].end());
		}
		constraints.B1duplicity[s].resize(constraints.B1[s].I_row_indices.size(), 1);
		for (eslocal r = constraints.B1subdomainsMap[s].size(); r < constraints.B1[s].nnz; r++) {
			constraints.B1subdomainsMap[s].push_back(constraints.B1[s].I_row_indices[r] - 1);
		}
		constraints.LB[s].resize(constraints.B1[s].nnz, -std::numeric_limits<double>::infinity());
	}

	constraints.B1clustersMap.reserve(constraints.B1clustersMap.size() + clusterDirichletSize);
	for (size_t i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		constraints.B1clustersMap.push_back({ (esglobal)i, environment->MPIrank });
	}

	ESINFO(DETAILS) << "Lambdas with Dirichlet in B1: " << constraints.B1[0].rows;
	ESTEST(MANDATORY) << "ESPRESO requires some nodes with Dirichlet condition." << (globalDirichletSize == 0 ? TEST_FAILED : TEST_PASSED);
}

std::vector<esglobal> EqualityConstraints::computeLambdasID(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	size_t loadStep = 0;

	std::vector<esglobal> lambdasID(elements.size() * DOFs.size(), -1);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(constraints._mesh.neighbours().begin(), constraints._mesh.neighbours().end(), neighbour) - constraints._mesh.neighbours().begin();
	};

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
	std::vector<size_t> offsets(threads);

	// neighbours x threads x data
	std::vector<std::vector<std::vector<esglobal> > > sLambdas(constraints._mesh.neighbours().size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > rLambdas(constraints._mesh.neighbours().size(), std::vector<std::vector<std::pair<esglobal,esglobal> > >(threads));

	std::vector<std::vector<Region*> > skippedRegions = getRegionsWithDOFs(constraints._mesh.regions(), loadStep, DOFs);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				size_t n = elements[e]->numberOfGlobalDomainsWithDOF(dof);
				if (n > 1 && (!constraints._configuration.redundant_lagrange || !commonRegion(skippedRegions[dof], elements[e]->regions()))) {
					if (elements[e]->clusters()[0] == environment->MPIrank) { // set lambda ID
						if (constraints._configuration.redundant_lagrange) {
							lambdasID[e * DOFs.size() + dof] = n * (n - 1) / 2;
							lambdasSize += n * (n - 1) / 2;
						} else {
							lambdasID[e * DOFs.size() + dof] = n - 1;
							lambdasSize += n - 1;
						}
						for (size_t c = 1; c < elements[e]->clusters().size(); c++) { // send to higher clusters
							sLambdas[n2i(elements[e]->clusters()[c])][t].push_back(e * DOFs.size() + dof);
						}
					} else { // pick ID from lower cluster
						rLambdas[n2i(elements[e]->clusters()[0])][t].push_back( // offset + lambda
								std::make_pair(elements[e]->clusterOffset(elements[e]->clusters()[0]), e * DOFs.size() + dof)
						);
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t numberOfClusterLambdas = Esutils::sizesToOffsets(offsets);
	size_t clusterOffset = numberOfClusterLambdas;
	size_t totalNumberOfLambdas = constraints.synchronizeOffsets(clusterOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + constraints.B1[0].rows;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esglobal offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (lambdasID[e * DOFs.size() + dof] > 0) {
					offset += lambdasID[e * DOFs.size() + dof];
					lambdasID[e * DOFs.size() + dof] = offset - lambdasID[e * DOFs.size() + dof];
				}
			}

		}

		for (size_t n = 0; n < constraints._mesh.neighbours().size(); n++) {
			for (size_t i = 0; i < sLambdas[n][t].size(); i++) {
				sLambdas[n][t][i] = lambdasID[sLambdas[n][t][i]];
			}
		}
	}

	std::vector<std::vector<esglobal> > rBuffer(constraints._mesh.neighbours().size());

	#pragma omp parallel for
	for (size_t n = 0; n < constraints._mesh.neighbours().size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sLambdas[n][0].insert(sLambdas[n][0].end(), sLambdas[n][t].begin(), sLambdas[n][t].end());
			rLambdas[n][0].insert(rLambdas[n][0].end(), rLambdas[n][t].begin(), rLambdas[n][t].end());
		}
		std::sort(rLambdas[n][0].begin(), rLambdas[n][0].end());
		rBuffer[n].resize(rLambdas[n][0].size());
	}

	std::vector<MPI_Request> req(constraints._mesh.neighbours().size());
	for (size_t n = 0; n < constraints._mesh.neighbours().size(); n++) {
		if (constraints._mesh.neighbours()[n] > environment->MPIrank) {
			MPI_Isend(sLambdas[n][0].data(), sizeof(esglobal) * sLambdas[n][0].size(), MPI_BYTE, constraints._mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
		if (constraints._mesh.neighbours()[n] < environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(esglobal) * rBuffer[n].size(), MPI_BYTE, constraints._mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
	}

	MPI_Waitall(constraints._mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

	for (size_t n = 0; n < constraints._mesh.neighbours().size(); n++) {
		for (size_t i = 0; i < rLambdas[n][0].size(); i++) {
			lambdasID[rLambdas[n][0][i].second] = rBuffer[n][i];
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		constraints.B1[p].rows += totalNumberOfLambdas;
	}
	constraints.block[Constraints::BLOCK::EQUALITY_CONSTRAINTS] += totalNumberOfLambdas;

	return lambdasID;
}

void EqualityConstraints::insertElementGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &K)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(constraints._mesh.neighbours().begin(), constraints._mesh.neighbours().end(), neighbour) - constraints._mesh.neighbours().begin();
	};

	std::vector<esglobal> lambdasID = computeLambdasID(constraints, elements, DOFs);

	std::vector<eslocal> permutation(lambdasID.size());
	std::iota(permutation.begin(), permutation.end(), 0);

	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return (lambdasID[i] < 0) ? false : (lambdasID[j] < 0) ? true : lambdasID[i] < lambdasID[j];
	});

	auto it = std::find_if(permutation.begin(), permutation.end(), [&] (eslocal i) { return lambdasID[i] == -1; });
	permutation.resize(it - permutation.begin());


	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, permutation.size());

	// threads x domains x data
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(constraints._mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(constraints._mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(constraints._mesh.parts()));
	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(constraints._mesh.parts()));

	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);

	auto findDomain = [&] (const Element *e, size_t d, size_t dof) -> eslocal {
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

	std::vector<std::vector<double> > diagonals;
	if (constraints._configuration.scaling) {
		diagonals.resize(permutation.size());
		std::vector<std::vector<double> > D(constraints._mesh.parts());

		#pragma omp parallel for
		for  (size_t p = 0; p < K.size(); p++) {
			D[p] = K[p].getDiagonal();
		}

		std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(constraints._mesh.neighbours().size()));
		std::vector<std::vector<double> > rBuffer(constraints._mesh.neighbours().size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

				const Element *e = elements[permutation[i] / DOFs.size()];
				size_t dof = permutation[i] % DOFs.size();

				for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
					if (*c != environment->MPIrank) {
						for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
							sBuffer[t][n2i(*c)].push_back(D[*d][e->DOFIndex(*d, dof)]);
						}
					}
				}

			}
		}

		for (size_t t = 1; t < threads; t++) {
			for (size_t n = 0; n < sBuffer[t].size(); n++) {
				sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
			}
		}

		std::vector<MPI_Request> req(constraints._mesh.neighbours().size());
		for (size_t n = 0; n < constraints._mesh.neighbours().size(); n++) {
			MPI_Isend(sBuffer[0][n].data(), sizeof(double) * sBuffer[0][n].size(), MPI_BYTE, constraints._mesh.neighbours()[n], 0, MPI_COMM_WORLD, req.data() + n);
		}

		int flag;
		size_t counter = 0;
		MPI_Status status;
		while (counter < constraints._mesh.neighbours().size()) {
			MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
			if (flag) {
				int count;
				MPI_Get_count(&status, MPI_BYTE, &count);
				rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(double));
				MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				counter++;
			}
		}

		MPI_Waitall(constraints._mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

		std::vector<eslocal> nPointer(constraints._mesh.neighbours().size());
		for (size_t i = 0; i < diagonals.size(); i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();

			for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
				if (*c == environment->MPIrank) {
					for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
						diagonals[i].push_back(D[*d][e->DOFIndex(*d, dof)]);
					}
				} else {
					for (eslocal d = 0; d < e->DOFCounter(*c, dof); d++) {
						diagonals[i].push_back(rBuffer[n2i(*c)][nPointer[n2i(*c)]++]);
					}
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();
			esglobal offset = 0;
			double duplicity = 0;
			if (constraints._configuration.scaling) {
				std::for_each(diagonals[i].begin(), diagonals[i].end(), [&] (double v) { duplicity += v; });
			} else {
				duplicity = e->numberOfGlobalDomainsWithDOF(dof);
			}

			eslocal diag1 = 0;
			for (auto c1 = e->clusters().begin(); c1 != e->clusters().end(); ++c1) {
				for (eslocal d1 = 0; d1 < e->DOFCounter(*c1, dof); d1++, diag1++) {

					eslocal diag2 = diag1 + 1;
					for (auto c2 = c1; c2 != e->clusters().end(); ++c2) {
						for (eslocal d2 = (*c1 == *c2 ? d1 + 1 : 0); d2 < e->DOFCounter(*c2, dof); d2++, diag2++) {

							if (*c1 == environment->MPIrank) {
								eslocal d = findDomain(e, d1, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
								vals[t][d].push_back(1);
								if (constraints._configuration.scaling) {
									dup[t][d].push_back(diagonals[i][diag2] / duplicity);
								} else {
									dup[t][d].push_back(1 / duplicity);
								}
							}

							if (*c2 == environment->MPIrank) {
								eslocal d = findDomain(e, d2, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
								vals[t][d].push_back(-1);
								if (constraints._configuration.scaling) {
									dup[t][d].push_back(diagonals[i][diag1] / duplicity);
								} else {
									dup[t][d].push_back(1 / duplicity);
								}
							}

							if (*c1 == environment->MPIrank || *c2 == environment->MPIrank) {
								cMap[t].push_back({ lambdasID[permutation[i]] + offset });
								if (*c1 == *c2) {
									cMap[t].back().push_back(*c1);
								} else if (*c1 == environment->MPIrank) {
									cMap[t].back().push_back(*c1);
									cMap[t].back().push_back(*c2);
								} else {
									cMap[t].back().push_back(*c2);
									cMap[t].back().push_back(*c1);
								}
							}

							offset++;
						}
					}
					if (!constraints._configuration.redundant_lagrange) {
						break;
					}
				}
				if (!constraints._configuration.redundant_lagrange) {
					break;
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		for (size_t t = 0; t < threads; t++) {
			constraints.B1[p].I_row_indices.insert(constraints.B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			constraints.B1[p].J_col_indices.insert(constraints.B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			constraints.B1[p].V_values.insert(constraints.B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			constraints.B1duplicity[p].insert(constraints.B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
		}
		constraints.B1[p].nnz = constraints.B1[p].I_row_indices.size();
		constraints.B1c[p].resize(constraints.B1[p].nnz, 0);
		constraints.LB[p].resize(constraints.B1[p].nnz, -std::numeric_limits<double>::infinity());
		for (eslocal r = constraints.B1subdomainsMap[p].size(); r < constraints.B1[p].nnz; r++) {
			constraints.B1subdomainsMap[p].push_back(constraints.B1[p].I_row_indices[r] - 1);
		}
	}

	for (size_t t = 0; t < threads; t++) {
		constraints.B1clustersMap.insert(constraints.B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
	}

	ESINFO(DETAILS) << "Lambdas in B1: " << constraints.B1[0].rows;
}

#ifdef HAVE_MORTAR

#include "mortar.h"

void EqualityConstraints::insertMortarGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	size_t loadStep = 0;

	std::vector<int> rows;
	std::vector<int> columns;
	std::vector<double> values;

	std::vector<std::vector<int> > masterElements;
	std::vector<Point_3D> masterCoordinates;
	std::vector<std::vector<int> > slaveElements;
	std::vector<Point_3D> slaveCoordinates;
	std::vector<int> nodes;

	for (size_t r = 0; r < constraints._mesh.regions().size(); r++) {
		if (loadStep < constraints._mesh.regions()[r]->settings.size() && constraints._mesh.regions()[r]->settings[loadStep].count(Property::NONMATCHING_ELEMENT)) {
			for (size_t i = 0; i < elements.size(); i++) {
					if (environment->MPIrank) {
						masterElements.push_back(std::vector<int>());
						for (size_t n = 0; n < elements[i]->nodes(); n++) {
							masterElements.back().push_back(elements[i]->node(n));
							nodes.push_back(elements[i]->node(n));
						}
					} else {
						slaveElements.push_back(std::vector<int>());
						for (size_t n = 0; n < elements[i]->nodes(); n++) {
							slaveElements.back().push_back(elements[i]->node(n));
							nodes.push_back(elements[i]->node(n));
						}
					}

			}
		}
	}

	if (!masterElements.size() && !slaveElements.size()) {
		// no MORTAR interface founded
		return;
	}

	std::sort(nodes.begin(), nodes.end());
	Esutils::removeDuplicity(nodes);

	for (size_t n = 0; n < nodes.size(); n++) {
		if (environment->MPIrank) {
			masterCoordinates.push_back(Point_3D());
			masterCoordinates.back().x = constraints._mesh.coordinates()[nodes[n]].x;
			masterCoordinates.back().y = constraints._mesh.coordinates()[nodes[n]].y;
			masterCoordinates.back().z = constraints._mesh.coordinates()[nodes[n]].z;
		} else {
			slaveCoordinates.push_back(Point_3D());
			slaveCoordinates.back().x = constraints._mesh.coordinates()[nodes[n]].x;
			slaveCoordinates.back().y = constraints._mesh.coordinates()[nodes[n]].y;
			slaveCoordinates.back().z = constraints._mesh.coordinates()[nodes[n]].z;
		}
	}

	std::vector<int> buffer;

	for (size_t e = 0; e < masterElements.size(); e++) {
		buffer.push_back(masterElements[e].size());
		for (size_t n = 0; n < masterElements[e].size(); n++) {
			masterElements[e][n] = std::lower_bound(nodes.begin(), nodes.end(), masterElements[e][n]) - nodes.begin();
			buffer.push_back(masterElements[e][n]);
		}
	}

	for (size_t e = 0; e < slaveElements.size(); e++) {
		for (size_t n = 0; n < slaveElements[e].size(); n++) {
			slaveElements[e][n] = std::lower_bound(nodes.begin(), nodes.end(), slaveElements[e][n]) - nodes.begin();
		}
	}

	if (environment->MPIrank) {
		MPI_Send(buffer.data(), buffer.size() * sizeof(int), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(masterCoordinates.data(), masterCoordinates.size() * sizeof(Point_3D), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
	} else {
		MPI_Status status;
		int size;

		// ELEMENTS
		MPI_Probe(1, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		buffer.resize(size / sizeof(int));
		MPI_Recv(buffer.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		// COORDINATES
		MPI_Probe(1, 1, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		masterCoordinates.resize(size / sizeof(Point_3D));
		MPI_Recv(masterCoordinates.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		for (size_t i = 0; i < buffer.size(); i++) {
			masterElements.push_back(std::vector<int>());
			for (size_t n = 0; n < buffer[i - n]; n++, i++) {
				masterElements.back().push_back(buffer[i + 1]);
			}
		}
	}

	if (!environment->MPIrank && (masterElements.size() || slaveElements.size())) {
		computeMortarEqualityConstraints(rows, columns, values, masterElements, masterCoordinates, slaveElements, slaveCoordinates);
	}
}

#else

void EqualityConstraints::insertMortarGluingToB1(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	size_t loadStep = 0;
	for (size_t r = 0; r < constraints._mesh.regions().size(); r++) {
		if (loadStep < constraints._mesh.regions()[r]->settings.size() && constraints._mesh.regions()[r]->settings[loadStep].count(Property::NONMATCHING_ELEMENT)) {
			ESINFO(GLOBAL_ERROR) << "Gluing of non-matching grids is not supported. Link MORTAR library!";
		}
	}
}

#endif


void EqualityConstraints::insertDomainGluingToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	if (!elements.size()) {
		return;
	}

	size_t lambdas = constraints.B0[0].rows;

	for (size_t e = 0; e < elements.size(); e++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (elements[e]->numberOfLocalDomainsWithDOF(dof) > 1) { // inner nodes are not glued
				const std::vector<eslocal> &DOFIndices = elements[e]->DOFsIndices();

				for (size_t d1 = 0, d2 = 1; d2 < elements[e]->domains().size(); d1++, d2++) {

					constraints.B0[elements[e]->domains()[d1]].I_row_indices.push_back(lambdas + 1);
					constraints.B0[elements[e]->domains()[d1]].J_col_indices.push_back(DOFIndices[d1 * DOFs.size() + dof] + 1);
					constraints.B0[elements[e]->domains()[d1]].V_values.push_back(1);

					constraints.B0[elements[e]->domains()[d2]].I_row_indices.push_back(lambdas + 1);
					constraints.B0[elements[e]->domains()[d2]].J_col_indices.push_back(DOFIndices[d2 * DOFs.size() + dof] + 1);
					constraints.B0[elements[e]->domains()[d2]].V_values.push_back(-1);

					lambdas++;
				}

			}
		}
	}


	#pragma omp parallel for
	for  (size_t p = 0; p < constraints._mesh.parts(); p++) {
		constraints.B0[p].rows = lambdas;
		constraints.B0[p].nnz = constraints.B0[p].I_row_indices.size();

		constraints.B0subdomainsMap[p].reserve(constraints.B0[p].nnz);
		for (eslocal i = constraints.B0subdomainsMap[p].size(); i < constraints.B0[p].nnz; i++) {
			constraints.B0subdomainsMap[p].push_back(constraints.B0[p].I_row_indices[i] - 1);
		}
	}

	ESINFO(DETAILS) << "Average number of lambdas in B0 is " << Info::averageValue(lambdas);
}

void EqualityConstraints::insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel)
{
	std::vector<Element*> el(elements);

	std::sort(el.begin(), el.end(), [] (Element* e1, Element* e2) {
		if (e1->domains().size() != e2->domains().size()) {
			return e1->domains().size() < e2->domains().size();
		}
		return e1->domains() < e2->domains();
	});

	std::vector<size_t> part;
	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (Element *e, size_t size) { return e->domains().size() < size; }) - el.begin());
	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((elements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
	for (size_t i = part[0] + 1; i < el.size(); i++) {
		if (i && el[i - 1]->domains() != el[i]->domains()) {
			part.push_back(i);
		}
	}
	part.push_back(el.size());

	#pragma omp parallel for
	for  (size_t p = 0; p < constraints._mesh.parts(); p++) {
		for (size_t i = 0; i < part.size() - 1; i++) {
			const std::vector<eslocal> &domains = el[part[i]]->domains();
			int sign = domains[0] == (eslocal)p ? 1 : domains[1] == (eslocal)p ? -1 : 0;
			if (sign == 0) {
				continue;
			}


			std::vector<Element*> nodes;
			for (size_t e = part[i]; e < part[i + 1]; e++) {
				for (size_t n = 0; n < el[e]->nodes(); n++) {
					nodes.push_back(constraints._mesh.nodes()[el[e]->node(n)]);
				}
			}
			std::sort(nodes.begin(), nodes.end());
			Esutils::removeDuplicity(nodes);

			for (eslocal col = 0; col < kernel[domains[0]].cols; col++) {
				for (size_t n = 0; n < nodes.size(); n++) {
					for (size_t dof = 0; dof < DOFs.size(); dof++) {
						constraints.B0[p].I_row_indices.push_back(i * kernel[0].cols + col + 1);
						constraints.B0[p].J_col_indices.push_back(nodes[n]->DOFIndex(p, dof) + 1);
						constraints.B0[p].V_values.push_back(sign * kernel[domains[0]].dense_values[kernel[domains[0]].rows * col + nodes[n]->DOFIndex(domains[0], dof)]);
					}
				}
			}
		}


		constraints.B0[p].rows = kernel[0].cols * (part.size() - 1);
		constraints.B0[p].nnz = constraints.B0[p].I_row_indices.size();
		constraints.B0subdomainsMap[p].reserve(constraints.B0[p].nnz);
		for (eslocal i = constraints.B0subdomainsMap[p].size(); i < constraints.B0[p].nnz; i++) {
			constraints.B0subdomainsMap[p].push_back(constraints.B0[p].I_row_indices[i] - 1);
		}
	}
}

void EqualityConstraints::insertKernelsToB0(Constraints &constraints, const std::vector<Element*> &elements, const std::vector<Element*> &DOFs, const std::vector<SparseMatrix> &kernel)
{
	std::vector<Element*> el(elements);

	std::sort(el.begin(), el.end(), [] (Element* e1, Element* e2) {
		if (e1->domains().size() != e2->domains().size()) {
			return e1->domains().size() < e2->domains().size();
		}
		return e1->domains() < e2->domains();
	});

	std::vector<size_t> part;
	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (Element *e, size_t size) { return e->domains().size() < size; }) - el.begin());
	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((elements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
	for (size_t i = part[0] + 1; i < el.size(); i++) {
		if (i && el[i - 1]->domains() != el[i]->domains()) {
			part.push_back(i);
		}
	}
	part.push_back(el.size());

	#pragma omp parallel for
	for  (size_t p = 0; p < constraints._mesh.parts(); p++) {
		for (size_t i = 0; i < part.size() - 1; i++) {
			const std::vector<eslocal> &domains = el[part[i]]->domains();
			int sign = domains[0] == (eslocal)p ? 1 : domains[1] == (eslocal)p ? -1 : 0;
			if (sign == 0) {
				continue;
			}

			std::vector<eslocal> interfaceDOFs;
			for (size_t e = part[i]; e < part[i + 1]; e++) {
				interfaceDOFs.insert(interfaceDOFs.end(), el[e]->DOFsIndices().begin(), el[e]->DOFsIndices().end());
			}

			std::sort(interfaceDOFs.begin(), interfaceDOFs.end());
			Esutils::removeDuplicity(interfaceDOFs);

			size_t DOFIndex = 0;
			for (eslocal col = 0; col < kernel[domains[0]].cols; col++) {
				for (size_t n = 0; n < interfaceDOFs.size(); n++) {
					constraints.B0[p].I_row_indices.push_back(i * kernel[0].cols + col + 1);
					constraints.B0[p].J_col_indices.push_back(DOFs[interfaceDOFs[n]]->DOFIndex(p, DOFIndex) + 1);
					constraints.B0[p].V_values.push_back(sign * kernel[domains[0]].dense_values[kernel[domains[0]].rows * col + DOFs[interfaceDOFs[n]]->DOFIndex(domains[0], DOFIndex)]);
				}
			}
		}


		constraints.B0[p].rows = kernel[0].cols * (part.size() - 1);
		constraints.B0[p].nnz = constraints.B0[p].I_row_indices.size();
		constraints.B0subdomainsMap[p].reserve(constraints.B0[p].nnz);
		for (eslocal i = constraints.B0subdomainsMap[p].size(); i < constraints.B0[p].nnz; i++) {
			constraints.B0subdomainsMap[p].push_back(constraints.B0[p].I_row_indices[i] - 1);
		}
	}
}

void EqualityConstraints::insertDirichletToB1(Instance &instance, const std::vector<Region*> &regions, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs, bool withRedundantMultiplier)
{
	size_t loadStep = 0;

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

	// part x thread x Dirichlet
	std::vector<std::vector<std::vector<esglobal> > > dirichlet(instance.domains, std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<double> > > dirichletValues(instance.domains, std::vector<std::vector<double> >(threads));

	std::vector<std::vector<Region*> > fixedRegions = getRegionsWithDOFs(regions, loadStep, DOFs);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (commonRegion(fixedRegions[dof], nodes[i]->regions())) {
					if (!withRedundantMultiplier && nodes[i]->clusters()[0] != environment->MPIrank) {
						continue;
					}
					const std::vector<eslocal>& indices = nodes[i]->DOFsIndices();
					double value = nodes[i]->getProperty(DOFs[dof], 0, loadStep, 0);
					for(size_t d = 0; d < nodes[i]->domains().size(); d++) {
						if (indices[d * DOFs.size() + dof] != -1) {
							dirichlet[nodes[i]->domains()[d]][t].push_back(indices[d * DOFs.size() + dof] + 1);
							dirichletValues[nodes[i]->domains()[d]][t].push_back(value);
						}
						if (!withRedundantMultiplier) {
							break;
						}
					}
				}
			}

		}
	}


	std::vector<size_t> dirichletSizes(instance.domains);
	#pragma omp parallel for
	for (size_t p = 0; p < instance.domains; p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[p][t].size();
		}
		dirichletSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t p = 0; p < instance.domains; p++) {
		clusterOffset += dirichletSizes[p];
		if (dirichletSizes[p]) {
			subdomainsWithDirichlet.push_back(p);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t globalDirichletSize = Constraints::synchronizeOffsets(clusterOffset);
	instance.block[Constraints::BLOCK::DIRICHLET] += globalDirichletSize;

	clusterOffset += instance.B1[0].rows;
	#pragma omp parallel for
	for (size_t p = 0; p < instance.domains; p++) {
		instance.B1[p].rows += globalDirichletSize;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		instance.B1[s].nnz += dirichletSizes[s];
		instance.B1[s].I_row_indices.reserve(instance.B1[s].nnz);
		instance.B1[s].J_col_indices.reserve(instance.B1[s].nnz);
		instance.B1[s].V_values.resize(instance.B1[s].nnz, 1);
	}

	Esutils::sizesToOffsets(dirichletSizes);
	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (eslocal i = 0; i < instance.B1[s].nnz; i++) {
			instance.B1[s].I_row_indices.push_back(clusterOffset + dirichletSizes[s] + i + 1);
		}
		for (size_t t = 0; t < threads; t++) {
			instance.B1[s].J_col_indices.insert(instance.B1[s].J_col_indices.end(), dirichlet[s][t].begin(), dirichlet[s][t].end());
			instance.B1c[s].insert(instance.B1c[s].end(), dirichletValues[s][t].begin(), dirichletValues[s][t].end());
		}
		instance.B1duplicity[s].resize(instance.B1[s].I_row_indices.size(), 1);
		for (eslocal r = instance.B1subdomainsMap[s].size(); r < instance.B1[s].nnz; r++) {
			instance.B1subdomainsMap[s].push_back(instance.B1[s].I_row_indices[r] - 1);
		}
		instance.LB[s].resize(instance.B1[s].nnz, -std::numeric_limits<double>::infinity());
	}

	instance.B1clustersMap.reserve(instance.B1clustersMap.size() + clusterDirichletSize);
	for (size_t i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		instance.B1clustersMap.push_back({ (esglobal)i, environment->MPIrank });
	}

	ESINFO(DETAILS) << "Lambdas with Dirichlet in B1: " << instance.B1[0].rows;
	ESTEST(MANDATORY) << "ESPRESO requires some nodes with Dirichlet condition." << (globalDirichletSize == 0 ? TEST_FAILED : TEST_PASSED);
}


std::vector<esglobal> EqualityConstraints::computeLambdasID(Instance &instance, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, bool withRedundantMultiplier)
{
	size_t loadStep = 0;

	std::vector<esglobal> lambdasID(elements.size() * DOFs.size(), -1);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
	std::vector<size_t> offsets(threads);

	// neighbours x threads x data
	std::vector<std::vector<std::vector<esglobal> > > sLambdas(neighbours.size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > rLambdas(neighbours.size(), std::vector<std::vector<std::pair<esglobal,esglobal> > >(threads));

	std::vector<std::vector<Region*> > skippedRegions = getRegionsWithDOFs(regions, loadStep, DOFs);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				size_t n = elements[e]->numberOfGlobalDomainsWithDOF(dof);
				if (n > 1 && (!withRedundantMultiplier || !commonRegion(skippedRegions[dof], elements[e]->regions()))) {
					if (elements[e]->clusters()[0] == environment->MPIrank) { // set lambda ID
						if (withRedundantMultiplier) {
							lambdasID[e * DOFs.size() + dof] = n * (n - 1) / 2;
							lambdasSize += n * (n - 1) / 2;
						} else {
							lambdasID[e * DOFs.size() + dof] = n - 1;
							lambdasSize += n - 1;
						}
						for (size_t c = 1; c < elements[e]->clusters().size(); c++) { // send to higher clusters
							sLambdas[n2i(elements[e]->clusters()[c])][t].push_back(e * DOFs.size() + dof);
						}
					} else { // pick ID from lower cluster
						rLambdas[n2i(elements[e]->clusters()[0])][t].push_back( // offset + lambda
								std::make_pair(elements[e]->clusterOffset(elements[e]->clusters()[0]), e * DOFs.size() + dof)
						);
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t numberOfClusterLambdas = Esutils::sizesToOffsets(offsets);
	size_t clusterOffset = numberOfClusterLambdas;
	size_t totalNumberOfLambdas = Constraints::synchronizeOffsets(clusterOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + instance.B1[0].rows;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esglobal offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (lambdasID[e * DOFs.size() + dof] > 0) {
					offset += lambdasID[e * DOFs.size() + dof];
					lambdasID[e * DOFs.size() + dof] = offset - lambdasID[e * DOFs.size() + dof];
				}
			}

		}

		for (size_t n = 0; n < neighbours.size(); n++) {
			for (size_t i = 0; i < sLambdas[n][t].size(); i++) {
				sLambdas[n][t][i] = lambdasID[sLambdas[n][t][i]];
			}
		}
	}

	std::vector<std::vector<esglobal> > rBuffer(neighbours.size());

	#pragma omp parallel for
	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sLambdas[n][0].insert(sLambdas[n][0].end(), sLambdas[n][t].begin(), sLambdas[n][t].end());
			rLambdas[n][0].insert(rLambdas[n][0].end(), rLambdas[n][t].begin(), rLambdas[n][t].end());
		}
		std::sort(rLambdas[n][0].begin(), rLambdas[n][0].end());
		rBuffer[n].resize(rLambdas[n][0].size());
	}

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] > environment->MPIrank) {
			MPI_Isend(sLambdas[n][0].data(), sizeof(esglobal) * sLambdas[n][0].size(), MPI_BYTE, neighbours[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
		if (neighbours[n] < environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(esglobal) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

	for (size_t n = 0; n < neighbours.size(); n++) {
		for (size_t i = 0; i < rLambdas[n][0].size(); i++) {
			lambdasID[rLambdas[n][0][i].second] = rBuffer[n][i];
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < instance.domains; p++) {
		instance.B1[p].rows += totalNumberOfLambdas;
	}
	instance.block[Constraints::BLOCK::EQUALITY_CONSTRAINTS] += totalNumberOfLambdas;

	return lambdasID;
}


void EqualityConstraints::insertElementGluingToB1(Instance &instance, const std::vector<int> &neighbours, const std::vector<Region*> &regions, const std::vector<Element*> &elements, const std::vector<Property> &DOFs, bool withRedundantMultiplier, bool withScaling)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	std::vector<esglobal> lambdasID = computeLambdasID(instance, neighbours, regions, elements, DOFs, withRedundantMultiplier);

	std::vector<eslocal> permutation(lambdasID.size());
	std::iota(permutation.begin(), permutation.end(), 0);

	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return (lambdasID[i] < 0) ? false : (lambdasID[j] < 0) ? true : lambdasID[i] < lambdasID[j];
	});

	auto it = std::find_if(permutation.begin(), permutation.end(), [&] (eslocal i) { return lambdasID[i] == -1; });
	permutation.resize(it - permutation.begin());


	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, permutation.size());

	// threads x domains x data
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(instance.domains));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(instance.domains));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(instance.domains));
	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(instance.domains));

	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);

	auto findDomain = [&] (const Element *e, size_t d, size_t dof) -> eslocal {
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

	std::vector<std::vector<double> > diagonals;
	if (withScaling) {
		diagonals.resize(permutation.size());
		std::vector<std::vector<double> > D(instance.domains);

		#pragma omp parallel for
		for  (size_t d = 0; d < instance.domains; d++) {
			D[d] = instance.K[d].getDiagonal();
		}

		std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(neighbours.size()));
		std::vector<std::vector<double> > rBuffer(neighbours.size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

				const Element *e = elements[permutation[i] / DOFs.size()];
				size_t dof = permutation[i] % DOFs.size();

				for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
					if (*c != environment->MPIrank) {
						for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
							sBuffer[t][n2i(*c)].push_back(D[*d][e->DOFIndex(*d, dof)]);
						}
					}
				}

			}
		}

		for (size_t t = 1; t < threads; t++) {
			for (size_t n = 0; n < sBuffer[t].size(); n++) {
				sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
			}
		}

		std::vector<MPI_Request> req(neighbours.size());
		for (size_t n = 0; n < neighbours.size(); n++) {
			MPI_Isend(sBuffer[0][n].data(), sizeof(double) * sBuffer[0][n].size(), MPI_BYTE, neighbours[n], 0, MPI_COMM_WORLD, req.data() + n);
		}

		int flag;
		size_t counter = 0;
		MPI_Status status;
		while (counter < neighbours.size()) {
			MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
			if (flag) {
				int count;
				MPI_Get_count(&status, MPI_BYTE, &count);
				rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(double));
				MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				counter++;
			}
		}

		MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);

		std::vector<eslocal> nPointer(neighbours.size());
		for (size_t i = 0; i < diagonals.size(); i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();

			for (auto c = e->clusters().begin(); c != e->clusters().end(); ++c) {
				if (*c == environment->MPIrank) {
					for (auto d = e->domains().begin(); d != e->domains().end(); d++) {
						diagonals[i].push_back(D[*d][e->DOFIndex(*d, dof)]);
					}
				} else {
					for (eslocal d = 0; d < e->DOFCounter(*c, dof); d++) {
						diagonals[i].push_back(rBuffer[n2i(*c)][nPointer[n2i(*c)]++]);
					}
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();
			esglobal offset = 0;
			double duplicity = 0;
			if (withScaling) {
				std::for_each(diagonals[i].begin(), diagonals[i].end(), [&] (double v) { duplicity += v; });
			} else {
				duplicity = e->numberOfGlobalDomainsWithDOF(dof);
			}

			eslocal diag1 = 0;
			for (auto c1 = e->clusters().begin(); c1 != e->clusters().end(); ++c1) {
				for (eslocal d1 = 0; d1 < e->DOFCounter(*c1, dof); d1++, diag1++) {

					eslocal diag2 = diag1 + 1;
					for (auto c2 = c1; c2 != e->clusters().end(); ++c2) {
						for (eslocal d2 = (*c1 == *c2 ? d1 + 1 : 0); d2 < e->DOFCounter(*c2, dof); d2++, diag2++) {

							if (*c1 == environment->MPIrank) {
								eslocal d = findDomain(e, d1, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
								vals[t][d].push_back(1);
								if (withScaling) {
									dup[t][d].push_back(diagonals[i][diag2] / duplicity);
								} else {
									dup[t][d].push_back(1 / duplicity);
								}
							}

							if (*c2 == environment->MPIrank) {
								eslocal d = findDomain(e, d2, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + 1);
								cols[t][d].push_back(e->DOFIndex(d, dof) + 1);
								vals[t][d].push_back(-1);
								if (withScaling) {
									dup[t][d].push_back(diagonals[i][diag1] / duplicity);
								} else {
									dup[t][d].push_back(1 / duplicity);
								}
							}

							if (*c1 == environment->MPIrank || *c2 == environment->MPIrank) {
								cMap[t].push_back({ lambdasID[permutation[i]] + offset });
								if (*c1 == *c2) {
									cMap[t].back().push_back(*c1);
								} else if (*c1 == environment->MPIrank) {
									cMap[t].back().push_back(*c1);
									cMap[t].back().push_back(*c2);
								} else {
									cMap[t].back().push_back(*c2);
									cMap[t].back().push_back(*c1);
								}
							}

							offset++;
						}
					}
					if (!withRedundantMultiplier) {
						break;
					}
				}
				if (!withRedundantMultiplier) {
					break;
				}
			}

		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < instance.domains; p++) {
		for (size_t t = 0; t < threads; t++) {
			instance.B1[p].I_row_indices.insert(instance.B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			instance.B1[p].J_col_indices.insert(instance.B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			instance.B1[p].V_values.insert(instance.B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			instance.B1duplicity[p].insert(instance.B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
		}
		instance.B1[p].nnz = instance.B1[p].I_row_indices.size();
		instance.B1c[p].resize(instance.B1[p].nnz, 0);
		instance.LB[p].resize(instance.B1[p].nnz, -std::numeric_limits<double>::infinity());
		for (eslocal r = instance.B1subdomainsMap[p].size(); r < instance.B1[p].nnz; r++) {
			instance.B1subdomainsMap[p].push_back(instance.B1[p].I_row_indices[r] - 1);
		}
	}

	for (size_t t = 0; t < threads; t++) {
		instance.B1clustersMap.insert(instance.B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
	}

	ESINFO(DETAILS) << "Lambdas in B1: " << instance.B1[0].rows;
}


void EqualityConstraints::insertCornersGluingToB0(Instance &instance, const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	if (!elements.size()) {
		return;
	}

	for (size_t d = 0; d < instance.domains; d++) {
		instance.B0[d].cols = instance.K[d].cols;
	}

	size_t lambdas = instance.B0[0].rows;

	for (size_t e = 0; e < elements.size(); e++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (elements[e]->numberOfLocalDomainsWithDOF(dof) > 1) { // inner nodes are not glued
				const std::vector<eslocal> &DOFIndices = elements[e]->DOFsIndices();

				for (size_t d1 = 0, d2 = 1; d2 < elements[e]->domains().size(); d1++, d2++) {

					instance.B0[elements[e]->domains()[d1]].I_row_indices.push_back(lambdas + 1);
					instance.B0[elements[e]->domains()[d1]].J_col_indices.push_back(DOFIndices[d1 * DOFs.size() + dof] + 1);
					instance.B0[elements[e]->domains()[d1]].V_values.push_back(1);

					instance.B0[elements[e]->domains()[d2]].I_row_indices.push_back(lambdas + 1);
					instance.B0[elements[e]->domains()[d2]].J_col_indices.push_back(DOFIndices[d2 * DOFs.size() + dof] + 1);
					instance.B0[elements[e]->domains()[d2]].V_values.push_back(-1);

					lambdas++;
				}

			}
		}
	}


	#pragma omp parallel for
	for  (size_t p = 0; p < instance.domains; p++) {
		instance.B0[p].rows = lambdas;
		instance.B0[p].nnz = instance.B0[p].I_row_indices.size();

		instance.B0subdomainsMap[p].reserve(instance.B0[p].nnz);
		for (eslocal i = instance.B0subdomainsMap[p].size(); i < instance.B0[p].nnz; i++) {
			instance.B0subdomainsMap[p].push_back(instance.B0[p].I_row_indices[i] - 1);
		}
	}

	ESINFO(DETAILS) << "Average number of lambdas in B0 is " << Info::averageValue(lambdas);
}

void EqualityConstraints::insertKernelsGluingToB0(Instance &instance, const std::vector<Element*> &elements, const std::vector<Element*> &nodes, const std::vector<Property> &DOFs)
{
	std::vector<Element*> el(elements);

	std::sort(el.begin(), el.end(), [] (Element* e1, Element* e2) {
		if (e1->domains().size() != e2->domains().size()) {
			return e1->domains().size() < e2->domains().size();
		}
		return e1->domains() < e2->domains();
	});

	std::vector<size_t> part;
	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (Element *e, size_t size) { return e->domains().size() < size; }) - el.begin());
	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((elements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
	for (size_t i = part[0] + 1; i < el.size(); i++) {
		if (i && el[i - 1]->domains() != el[i]->domains()) {
			part.push_back(i);
		}
	}
	part.push_back(el.size());

	#pragma omp parallel for
	for  (size_t p = 0; p < instance.domains; p++) {
		for (size_t i = 0; i < part.size() - 1; i++) {
			const std::vector<eslocal> &domains = el[part[i]]->domains();
			int sign = domains[0] == (eslocal)p ? 1 : domains[1] == (eslocal)p ? -1 : 0;
			if (sign == 0) {
				continue;
			}


			std::vector<Element*> nodesOnInterface;
			for (size_t e = part[i]; e < part[i + 1]; e++) {
				for (size_t n = 0; n < el[e]->nodes(); n++) {
					nodesOnInterface.push_back(nodes[el[e]->node(n)]);
				}
			}
			std::sort(nodesOnInterface.begin(), nodesOnInterface.end());
			Esutils::removeDuplicity(nodesOnInterface);

			for (eslocal col = 0; col < instance.R1[domains[0]].cols; col++) {
				for (size_t n = 0; n < nodesOnInterface.size(); n++) {
					for (size_t dof = 0; dof < DOFs.size(); dof++) {
						instance.B0[p].I_row_indices.push_back(i * instance.R1[0].cols + col + 1);
						instance.B0[p].J_col_indices.push_back(nodesOnInterface[n]->DOFIndex(p, dof) + 1);
						instance.B0[p].V_values.push_back(sign * instance.R1[domains[0]].dense_values[instance.R1[domains[0]].rows * col + nodesOnInterface[n]->DOFIndex(domains[0], dof)]);
					}
				}
			}
		}


		instance.B0[p].rows = instance.R1[0].cols * (part.size() - 1);
		instance.B0[p].nnz = instance.B0[p].I_row_indices.size();
		instance.B0subdomainsMap[p].reserve(instance.B0[p].nnz);
		for (eslocal i = instance.B0subdomainsMap[p].size(); i < instance.B0[p].nnz; i++) {
			instance.B0subdomainsMap[p].push_back(instance.B0[p].I_row_indices[i] - 1);
		}
	}
}

