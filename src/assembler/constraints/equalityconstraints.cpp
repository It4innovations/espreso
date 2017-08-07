
#include "equalityconstraints.h"

#include "../instance.h"
#include "../step.h"

#include "../../basis/utilities/utils.h"
#include "../../basis/utilities/communication.h"

#include "../../solver/generic/SparseMatrix.h"

#include "../../mesh/elements/element.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/settings/property.h"

#include "../../configuration/configuration.h"
#include "../../configuration/environment.h"

#include <numeric>

using namespace espreso;

EqualityConstraints::EqualityConstraints(
		Instance &instance,
		const Mesh &mesh,
		const std::vector<Element*> &gluedElements,
		const std::vector<Element*> &gluedInterfaceElements,
		const std::vector<Property> &gluedDOFs,
		const std::vector<size_t> &gluedDOFsMeshOffsets,
		bool interfaceElementContainsGluedDOFs)
: _instance(instance),
  _mesh(mesh),
  _gluedElements(gluedElements),
  _gluedInterfaceElements(gluedInterfaceElements),
  _gluedDOFs(gluedDOFs),
  _gluedDOFsMeshOffsets(gluedDOFsMeshOffsets),
  _interfaceElementContainsGluedDOFs(interfaceElementContainsGluedDOFs)
{

}

void EqualityConstraints::insertDirichletToB1(const Step &step, bool withRedundantMultiplier)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _gluedElements.size());

	// part x thread x Dirichlet
	std::vector<std::vector<std::vector<esglobal> > > dirichlet(_instance.domains, std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<double> > > dirichletValues(_instance.domains, std::vector<std::vector<double> >(threads));

	std::vector<std::vector<Region*> > fixedRegions = _mesh.getRegionsWithProperties(step.step, _gluedDOFs);

	double temp = 0; // irrelevant -> set to zero

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
				if (Mesh::commonRegion(fixedRegions[dof], _gluedElements[i]->regions())) {
					if (!withRedundantMultiplier && _gluedElements[i]->clusters()[0] != environment->MPIrank) {
						continue;
					}
					double value = _gluedElements[i]->getProperty(_gluedDOFs[dof], 0, step.step, step.currentTime, temp, 0);
					for(size_t d = 0; d < _gluedElements[i]->domains().size(); d++) {
						if (_gluedElements[i]->DOFIndex(_gluedElements[i]->domains()[d], _gluedDOFsMeshOffsets[dof]) != -1) {
							dirichlet[_gluedElements[i]->domains()[d]][t].push_back(_gluedElements[i]->DOFIndex(_gluedElements[i]->domains()[d], _gluedDOFsMeshOffsets[dof]) + 1);
							dirichletValues[_gluedElements[i]->domains()[d]][t].push_back(value);
						}
						if (!withRedundantMultiplier) {
							break;
						}
					}
				}
			}

		}
	}


	std::vector<size_t> dirichletSizes(_instance.domains);
	#pragma omp parallel for
	for (size_t p = 0; p < _instance.domains; p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[p][t].size();
		}
		dirichletSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t p = 0; p < _instance.domains; p++) {
		clusterOffset += dirichletSizes[p];
		if (dirichletSizes[p]) {
			subdomainsWithDirichlet.push_back(p);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t globalDirichletSize = Communication::exscan(clusterOffset);
	_instance.block[Instance::CONSTRAINT::DIRICHLET] += globalDirichletSize;

	clusterOffset += _instance.B1[0].rows;
	#pragma omp parallel for
	for (size_t p = 0; p < _instance.domains; p++) {
		_instance.B1[p].rows += globalDirichletSize;
		_instance.B1[p].cols = _instance.domainDOFCount[p];
	}

	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		_instance.B1[s].nnz += dirichletSizes[s];
		_instance.B1[s].I_row_indices.reserve(_instance.B1[s].nnz);
		_instance.B1[s].J_col_indices.reserve(_instance.B1[s].nnz);
		_instance.B1[s].V_values.resize(_instance.B1[s].nnz, 1);
	}

	Esutils::sizesToOffsets(dirichletSizes);
	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (eslocal i = 0; i < _instance.B1[s].nnz; i++) {
			_instance.B1[s].I_row_indices.push_back(clusterOffset + dirichletSizes[s] + i + 1);
		}
		for (size_t t = 0; t < threads; t++) {
			_instance.B1[s].J_col_indices.insert(_instance.B1[s].J_col_indices.end(), dirichlet[s][t].begin(), dirichlet[s][t].end());
			_instance.B1c[s].insert(_instance.B1c[s].end(), dirichletValues[s][t].begin(), dirichletValues[s][t].end());
		}
		_instance.B1duplicity[s].resize(_instance.B1[s].I_row_indices.size(), 1);
		for (eslocal r = _instance.B1subdomainsMap[s].size(); r < _instance.B1[s].nnz; r++) {
			_instance.B1subdomainsMap[s].push_back(_instance.B1[s].I_row_indices[r] - 1);
		}
		_instance.LB[s].resize(_instance.B1[s].nnz, -std::numeric_limits<double>::infinity());
	}

	_instance.B1clustersMap.reserve(_instance.B1clustersMap.size() + clusterDirichletSize);
	for (size_t i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		_instance.B1clustersMap.push_back({ (esglobal)i, environment->MPIrank });
	}

	ESINFO(EXHAUSTIVE) << "Lambdas with Dirichlet in B1: " << _instance.B1[0].rows;
}

std::vector<esglobal> EqualityConstraints::computeLambdasID(const Step &step, bool withRedundantMultiplier)
{
	std::vector<esglobal> lambdasID(_gluedElements.size() * _gluedDOFs.size(), -1);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, _gluedElements.size());
	std::vector<size_t> offsets(threads);

	// neighbours x threads x data
	std::vector<std::vector<std::vector<esglobal> > > sLambdas(threads, std::vector<std::vector<esglobal> >(_mesh.neighbours().size()));
	std::vector<std::vector<std::vector<std::pair<esglobal, esglobal> > > > rLambdas(threads, std::vector<std::vector<std::pair<esglobal,esglobal> > >(_mesh.neighbours().size()));

	std::vector<std::vector<Region*> > skippedRegions = _mesh.getRegionsWithProperties(step.step, _gluedDOFs);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
				size_t n = _gluedElements[e]->numberOfGlobalDomainsWithDOF(_gluedDOFsMeshOffsets[dof]);
				if (n > 1 && (!withRedundantMultiplier || !Mesh::commonRegion(skippedRegions[dof], _gluedElements[e]->regions()))) {
					if (_gluedElements[e]->clusters()[0] == environment->MPIrank) { // set lambda ID
						if (withRedundantMultiplier) {
							lambdasID[e * _gluedDOFs.size() + dof] = n * (n - 1) / 2;
							lambdasSize += n * (n - 1) / 2;
						} else {
							lambdasID[e * _gluedDOFs.size() + dof] = n - 1;
							lambdasSize += n - 1;
						}
						for (size_t c = 1; c < _gluedElements[e]->clusters().size(); c++) { // send to higher clusters
							sLambdas[t][n2i(_gluedElements[e]->clusters()[c])].push_back(e * _gluedDOFs.size() + dof);
						}
					} else { // pick ID from lower cluster
						rLambdas[t][n2i(_gluedElements[e]->clusters()[0])].push_back( // offset + lambda
								std::make_pair(_gluedElements[e]->clusterOffset(_gluedElements[e]->clusters()[0]), e * _gluedDOFs.size() + dof)
						);
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t numberOfClusterLambdas = Esutils::sizesToOffsets(offsets);
	size_t clusterOffset = numberOfClusterLambdas;
	size_t totalNumberOfLambdas = Communication::exscan(clusterOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + _instance.B1[0].rows;
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esglobal offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
				if (lambdasID[e * _gluedDOFs.size() + dof] > 0) {
					offset += lambdasID[e * _gluedDOFs.size() + dof];
					lambdasID[e * _gluedDOFs.size() + dof] = offset - lambdasID[e * _gluedDOFs.size() + dof];
				}
			}

		}

		for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
			for (size_t i = 0; i < sLambdas[t][n].size(); i++) {
				sLambdas[t][n][i] = lambdasID[sLambdas[t][n][i]];
			}
		}
	}

	std::vector<std::vector<esglobal> > rBuffer(_mesh.neighbours().size());

	#pragma omp parallel for
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		for (size_t t = 1; t < threads; t++) {
			sLambdas[0][n].insert(sLambdas[0][n].end(), sLambdas[t][n].begin(), sLambdas[t][n].end());
			rLambdas[0][n].insert(rLambdas[0][n].end(), rLambdas[t][n].begin(), rLambdas[t][n].end());
		}
		std::sort(rLambdas[0][n].begin(), rLambdas[0][n].end());
		rBuffer[n].resize(rLambdas[0][n].size());
	}

	if (!Communication::receiveLowerKnownSize(sLambdas[0], rBuffer, _mesh.neighbours())) {
		ESINFO(ERROR) << "problem while synchronization of lambdas.";
	}

	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		for (size_t i = 0; i < rLambdas[0][n].size(); i++) {
			lambdasID[rLambdas[0][n][i].second] = rBuffer[n][i];
		}
	}

	#pragma omp parallel for
	for (size_t p = 0; p < _instance.domains; p++) {
		_instance.B1[p].rows += totalNumberOfLambdas;
		_instance.B1[p].cols = _instance.domainDOFCount[p];
	}
	_instance.block[Instance::CONSTRAINT::EQUALITY_CONSTRAINTS] += totalNumberOfLambdas;

	return lambdasID;
}

void EqualityConstraints::insertElementGluingToB1(const Step &step, bool withRedundantMultiplier, bool withScaling)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	std::vector<esglobal> lambdasID = computeLambdasID(step, withRedundantMultiplier);

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
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(_instance.domains));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(_instance.domains));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(_instance.domains));
	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(_instance.domains));

	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);

	auto findDomain = [&] (const Element *e, size_t d, size_t dof) -> eslocal {
		auto &DOFIndices = e->DOFsIndices();
		size_t c = 0, DOFsSize = DOFIndices.size() / e->domains().size();
		for (size_t i = 0; i < e->domains().size(); i++) {
			if (DOFIndices[i * DOFsSize + dof] != -1) {
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
		std::vector<std::vector<double> > D(_instance.domains);

		#pragma omp parallel for
		for  (size_t d = 0; d < _instance.domains; d++) {
			D[d] = _instance.K[d].getDiagonal();
		}

		std::vector<std::vector<std::vector<double> > > sBuffer(threads, std::vector<std::vector<double> >(_mesh.neighbours().size()));
		std::vector<std::vector<double> > rBuffer(_mesh.neighbours().size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

				const Element *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
				size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];

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

		if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, _mesh.neighbours())) {
			ESINFO(ERROR) << "problem while exchange K diagonal in B1 scaling.";
		}

		std::vector<eslocal> nPointer(_mesh.neighbours().size());
		for (size_t i = 0; i < diagonals.size(); i++) {

			const Element *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
			size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];

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

			const Element *e = _gluedElements[permutation[i] / _gluedDOFs.size()];
			size_t dof = _gluedDOFsMeshOffsets[permutation[i] % _gluedDOFs.size()];
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
	for (size_t p = 0; p < _instance.domains; p++) {
		for (size_t t = 0; t < threads; t++) {
			_instance.B1[p].I_row_indices.insert(_instance.B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			_instance.B1[p].J_col_indices.insert(_instance.B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			_instance.B1[p].V_values.insert(_instance.B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			_instance.B1duplicity[p].insert(_instance.B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
		}
		_instance.B1[p].cols = _instance.domainDOFCount[p];
		_instance.B1[p].nnz = _instance.B1[p].I_row_indices.size();
		_instance.B1c[p].resize(_instance.B1[p].nnz, 0);
		_instance.LB[p].resize(_instance.B1[p].nnz, -std::numeric_limits<double>::infinity());
		for (eslocal r = _instance.B1subdomainsMap[p].size(); r < _instance.B1[p].nnz; r++) {
			_instance.B1subdomainsMap[p].push_back(_instance.B1[p].I_row_indices[r] - 1);
		}
	}

	for (size_t t = 0; t < threads; t++) {
		_instance.B1clustersMap.insert(_instance.B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
	}

	ESINFO(EXHAUSTIVE) << "Lambdas in B1: " << _instance.B1[0].rows;
}

void EqualityConstraints::insertCornersGluingToB0()
{
	const std::vector<Element*> &corners = _mesh.corners();
	if (!corners.size()) {
		ESINFO(ERROR) << "Mesh not contains corners.";
		return;
	}

	for (size_t d = 0; d < _instance.domains; d++) {
		_instance.B0[d].cols = _instance.K[d].cols;
	}

	size_t lambdas = _instance.B0[0].rows;

	for (size_t e = 0; e < corners.size(); e++) {
		for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
			if (corners[e]->numberOfLocalDomainsWithDOF(_gluedDOFsMeshOffsets[dof]) > 1) { // inner nodes are not glued

				for (size_t d1 = 0, d2 = 1; d2 < corners[e]->domains().size(); d1++, d2++) {

					_instance.B0[corners[e]->domains()[d1]].I_row_indices.push_back(lambdas + 1);
					_instance.B0[corners[e]->domains()[d1]].J_col_indices.push_back(corners[e]->DOFIndex(corners[e]->domains()[d1], _gluedDOFsMeshOffsets[dof]) + 1);
					_instance.B0[corners[e]->domains()[d1]].V_values.push_back(1);

					_instance.B0[corners[e]->domains()[d2]].I_row_indices.push_back(lambdas + 1);
					_instance.B0[corners[e]->domains()[d2]].J_col_indices.push_back(corners[e]->DOFIndex(corners[e]->domains()[d2], _gluedDOFsMeshOffsets[dof]) + 1);
					_instance.B0[corners[e]->domains()[d2]].V_values.push_back(-1);

					lambdas++;
				}

			}
		}
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < _instance.domains; p++) {
		_instance.B0[p].rows = lambdas;
		_instance.B0[p].cols = _instance.domainDOFCount[p];
		_instance.B0[p].nnz = _instance.B0[p].I_row_indices.size();

		_instance.B0subdomainsMap[p].reserve(_instance.B0[p].nnz);
		for (eslocal i = _instance.B0subdomainsMap[p].size(); i < _instance.B0[p].nnz; i++) {
			_instance.B0subdomainsMap[p].push_back(_instance.B0[p].I_row_indices[i] - 1);
		}
	}

	ESINFO(EXHAUSTIVE) << "Average number of lambdas in B0 is " << Info::averageValue(lambdas);
}

void EqualityConstraints::insertKernelsGluingToB0(const std::vector<SparseMatrix> &kernels)
{
	std::vector<Element*> el(_gluedInterfaceElements);

	std::sort(el.begin(), el.end(), [] (Element* e1, Element* e2) {
		if (e1->domains().size() != e2->domains().size()) {
			return e1->domains().size() < e2->domains().size();
		}
		return e1->domains() < e2->domains();
	});

	std::vector<size_t> part;
	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (Element *e, size_t size) { return e->domains().size() < size; }) - el.begin());
	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((_gluedInterfaceElements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
	for (size_t i = part[0] + 1; i < el.size(); i++) {
		if (i && el[i - 1]->domains() != el[i]->domains()) {
			part.push_back(i);
		}
	}
	part.push_back(el.size());

	std::vector<eslocal> rowIndex;
	std::vector<eslocal> clusterRowIndex(_instance.clustersMap.size(), 1);
	for (size_t i = 0; i < part.size() - 1; i++) {
		const std::vector<eslocal> &domains = el[part[i]]->domains();
		if (_instance.clustersMap[domains[0]] == _instance.clustersMap[domains[1]]) {
			eslocal master = kernels[domains[0]].cols > kernels[domains[1]].cols ? domains[0] : domains[1];
			eslocal rows = kernels[master].cols > 0 ? kernels[master].cols : 1;
			rowIndex.push_back(clusterRowIndex[_instance.clustersMap[domains[0]]]);
			clusterRowIndex[_instance.clustersMap[domains[0]]] += rows;
		} else {
			rowIndex.push_back(-1);
		}
	}

	#pragma omp parallel for
	for  (size_t p = 0; p < _instance.domains; p++) {
		for (size_t i = 0; i < part.size() - 1; i++) {
			const std::vector<eslocal> &domains = el[part[i]]->domains();
			if (_instance.clustersMap[domains[0]] != _instance.clustersMap[domains[1]]) {
				continue;
			}

			int sign = domains[0] == (eslocal)p ? 1 : domains[1] == (eslocal)p ? -1 : 0;
			if (sign == 0) {
				continue;
			}

			std::vector<Element*> DOFsOnInterface;
			for (size_t e = part[i]; e < part[i + 1]; e++) {
				if (_interfaceElementContainsGluedDOFs) {
					for (size_t n = 0; n < el[e]->DOFsIndices().size(); n++) {
						DOFsOnInterface.push_back(_gluedElements[el[e]->DOFsIndices()[n]]);
					}
				} else {
					for (size_t n = 0; n < el[e]->nodes(); n++) {
						DOFsOnInterface.push_back(_gluedElements[el[e]->node(n)]);
					}
				}
			}
			std::sort(DOFsOnInterface.begin(), DOFsOnInterface.end());
			Esutils::removeDuplicity(DOFsOnInterface);

			eslocal master = kernels[domains[0]].cols > kernels[domains[1]].cols ? domains[0] : domains[1];
			if (kernels[master].cols == 0) {
				for (size_t n = 0; n < DOFsOnInterface.size(); n++) {
					for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
						_instance.B0[p].I_row_indices.push_back(rowIndex[i]);
						_instance.B0[p].J_col_indices.push_back(DOFsOnInterface[n]->DOFIndex(p, _gluedDOFsMeshOffsets[dof]) + 1);
						_instance.B0[p].V_values.push_back(sign);
					}
				}
			} else {
				for (eslocal col = 0; col < kernels[master].cols; col++) {
					for (size_t n = 0; n < DOFsOnInterface.size(); n++) {
						for (size_t dof = 0; dof < _gluedDOFs.size(); dof++) {
							_instance.B0[p].I_row_indices.push_back(rowIndex[i] + col);
							_instance.B0[p].J_col_indices.push_back(DOFsOnInterface[n]->DOFIndex(p, _gluedDOFsMeshOffsets[dof]) + 1);
							_instance.B0[p].V_values.push_back(sign * kernels[master].dense_values[kernels[master].rows * col + DOFsOnInterface[n]->DOFIndex(master, _gluedDOFsMeshOffsets[dof])]);
						}
					}
				}
			}
		}


		_instance.B0[p].rows = clusterRowIndex[_instance.clustersMap[p]] - 1;
		_instance.B0[p].cols = _instance.domainDOFCount[p];
		_instance.B0[p].nnz = _instance.B0[p].I_row_indices.size();
		_instance.B0subdomainsMap[p].reserve(_instance.B0[p].nnz);
		for (eslocal i = _instance.B0subdomainsMap[p].size(); i < _instance.B0[p].nnz; i++) {
			_instance.B0subdomainsMap[p].push_back(_instance.B0[p].I_row_indices[i] - 1);
		}
	}
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
		MPI_Send(buffer.data(), buffer.size() * sizeof(int), MPI_BYTE, 0, 0, environment->MPICommunicator);
		MPI_Send(masterCoordinates.data(), masterCoordinates.size() * sizeof(Point_3D), MPI_BYTE, 0, 1, environment->MPICommunicator);
	} else {
		MPI_Status status;
		int size;

		// ELEMENTS
		MPI_Probe(1, 0, environment->MPICommunicator, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		buffer.resize(size / sizeof(int));
		MPI_Recv(buffer.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, environment->MPICommunicator, MPI_STATUSES_IGNORE);

		// COORDINATES
		MPI_Probe(1, 1, environment->MPICommunicator, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		masterCoordinates.resize(size / sizeof(Point_3D));
		MPI_Recv(masterCoordinates.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, environment->MPICommunicator, MPI_STATUSES_IGNORE);

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


#endif

