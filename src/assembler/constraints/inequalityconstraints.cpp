
#include "inequalityconstraints.h"

using namespace espreso;

void InequalityConstraints::insertLowerBoundToB1(Constraints &constraints, const std::vector<Property> &eDOFs, const std::vector<Property> &boundDOFs)
{
	size_t threads = environment->OMP_NUM_THREADS;
	std::vector<size_t> distribution;

	// part x thread x indices
	std::vector<std::vector<std::vector<eslocal> > > indices(constraints._mesh.parts(), std::vector<std::vector<eslocal> >(threads));
	std::vector<std::vector<std::vector<double> > > values(constraints._mesh.parts(), std::vector<std::vector<double> >(threads));
	std::vector<std::vector<std::vector<double> > > normal(constraints._mesh.parts(), std::vector<std::vector<double> >(threads));

	size_t loadStep = 0;
	for (size_t r = 0; r < constraints._mesh.regions().size(); r++) {
		const Region *region = constraints._mesh.regions()[r];
		if (loadStep < region->settings.size()) {

			for (size_t dof = 0; dof < boundDOFs.size(); dof++) {
				auto settings = region->settings[loadStep].find(boundDOFs[dof]);
				auto normalDirection = region->settings[loadStep].find(Property::NORMAL_DIRECTION);
				if (settings != region->settings[loadStep].end()) {
					if (normalDirection == region->settings[loadStep].end()) {
						ESINFO(GLOBAL_ERROR) << "You have to set normal direction for elements with obstacle";
					}

					distribution = Esutils::getDistribution(threads, region->elements.size());

					#pragma omp parallel for
					for (size_t t = 0; t < threads; t++) {
						for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
							if (region->elements[n]->clusters()[0] != environment->MPIrank) {
								continue;
							}
							Point direction;
							direction.x = normalDirection->second.back()->evaluate(Point(1, 0, 0));
							direction.y = normalDirection->second.back()->evaluate(Point(0, 1, 0));
							direction.z = normalDirection->second.back()->evaluate(Point(0, 0, 1));
							direction.normalize();

							double value[3] { direction.x, direction.y, direction.z };

							eslocal domain = region->elements[n]->domains().front();
							for (size_t dof = 0; dof < eDOFs.size(); dof++) {
								if (value[dof]) {
									indices[domain][t].push_back(region->elements[n]->DOFIndex(domain, dof) + IJVMatrixIndexing);
									values[domain][t].push_back(settings->second.back()->evaluate(n));
									normal[domain][t].push_back(std::abs(value[dof]) * values[domain][t].back());
								}
							}
						}
					}

				}
			}

		}
	}

	std::vector<size_t> indicesSizes(constraints._mesh.parts());
	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += indices[p][t].size();
		}
		indicesSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithLowerBounds;
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		clusterOffset += indicesSizes[p];
		if (indicesSizes[p]) {
			subdomainsWithLowerBounds.push_back(p);
		}
	}

	size_t clusterIndicesSize = clusterOffset;
	size_t globalIndicesSize = constraints.synchronizeOffsets(clusterOffset);

	clusterOffset += constraints.B1[0].rows;
	#pragma omp parallel for
	for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		constraints.B1[p].rows += globalIndicesSize;
	}

	Esutils::sizesToOffsets(indicesSizes);
	#pragma omp parallel for
	for (size_t i = 0; i < subdomainsWithLowerBounds.size(); i++) {
		size_t s = subdomainsWithLowerBounds[i];
		for (size_t t = 0, row = clusterOffset + indicesSizes[s]; t < threads; row += indices[s][t++].size()) {
			constraints.B1[s].I_row_indices.resize(constraints.B1[s].I_row_indices.size() + indices[s][t].size());
			std::iota(constraints.B1[s].I_row_indices.end() - indices[s][t].size(), constraints.B1[s].I_row_indices.end(), row + IJVMatrixIndexing);
			constraints.B1[s].J_col_indices.insert(constraints.B1[s].J_col_indices.end(), indices[s][t].begin(), indices[s][t].end());
			constraints.B1[s].V_values.insert(constraints.B1[s].V_values.end(), values[s][t].begin(), values[s][t].end());
			constraints.B1c[s].insert(constraints.B1c[s].end(), normal[s][t].begin(), normal[s][t].end());
		}
		constraints.B1[s].nnz = constraints.B1[s].I_row_indices.size();
		constraints.B1duplicity[s].resize(constraints.B1[s].I_row_indices.size(), 1);
		for (eslocal r = constraints.B1subdomainsMap[s].size(); r < constraints.B1[s].nnz; r++) {
			constraints.B1subdomainsMap[s].push_back(constraints.B1[s].I_row_indices[r] - 1);
		}
		constraints.LB[s].resize(constraints.B1[s].nnz, 0);
	}

	constraints.B1clustersMap.reserve(constraints.B1clustersMap.size() + clusterIndicesSize);
	for (size_t i = clusterOffset; i < clusterOffset + clusterIndicesSize; i++) {
		constraints.B1clustersMap.push_back({ (eslocal)i, environment->MPIrank });
	}

	constraints.block[Constraints::BLOCK::INEQUALITY_CONSTRAINTS] += globalIndicesSize;
	ESINFO(DETAILS) << "Lambdas with lower bounds in B1: " << globalIndicesSize;
}



