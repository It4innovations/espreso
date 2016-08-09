
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
				const std::vector<eslocal> &DOFIndices = elements[e]->DOFsIndices();
				for (size_t i = 0; i < DOFs.size(); i++) {
					if (!elements[e]->settings().isSet(DOFs[i])) { // Dirichlet is not glued
						size_t n = 0;
						for (size_t d = 0; d < elements[e]->domains().size(); d++) {
							if (DOFIndices[d * DOFs.size() + i] != -1) {
								n++;
							}
						}
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

void EqualityGluing::insertClusterGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{

}

void EqualityGluing::insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{

}



