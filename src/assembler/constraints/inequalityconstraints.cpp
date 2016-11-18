
#include "inequalityconstraints.h"

using namespace espreso;

void InequalityConstraints::insertLowerBoundToB1(Constraints &constraints, const std::vector<Element*> &nodes, const std::vector<Property> &eDOFs, const std::vector<Property> &boundDOFs)
{
	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

	// part x thread x indices
	std::vector<std::vector<std::vector<esglobal> > > indices(constraints._mesh.parts(), std::vector<std::vector<eslocal> >(threads));
	std::vector<std::vector<std::vector<double> > > values(constraints._mesh.parts(), std::vector<std::vector<double> >(threads));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {

			for (size_t dof = 0; dof < boundDOFs.size(); dof++) {
				if (nodes[n]->clusters()[0] != config::env::MPIrank) {
					continue;
				}
				if (nodes[n]->settings().isSet(boundDOFs[dof])) {
					double value = nodes[n]->settings(boundDOFs[dof]).back()->evaluate(n);
					indices[nodes[n]->domains().front()][t].push_back(n);
					values[nodes[n]->domains().front()][t].push_back(value);
				}
			}

		}
	}

	std::vector<size_t> indicesSizes(constraints._mesh.parts());
	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < constraints._mesh.parts(); p++) {
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
	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < constraints._mesh.parts(); p++) {
		constraints.B1[p].rows += globalIndicesSize;
	}

	Esutils::sizesToOffsets(indicesSizes);
	#pragma cilk grainsize = 1
	cilk_for (size_t i = 0; i < subdomainsWithLowerBounds.size(); i++) {
		size_t s = subdomainsWithLowerBounds[i];
		for (size_t t = 0, row = clusterOffset + indicesSizes[s]; t < threads; t++) {
			for (size_t n = 0; n < indices[s][t].size(); n++, row++) {
				if (nodes[indices[s][t][n]]->settings().isSet(Property::NORMAL_DIRECTION)) {
					Point normal;
					normal.x = nodes[indices[s][t][n]]->settings(Property::NORMAL_DIRECTION).back()->evaluate(Point(1, 0, 0));
					normal.y = nodes[indices[s][t][n]]->settings(Property::NORMAL_DIRECTION).back()->evaluate(Point(0, 1, 0));
					normal.z = nodes[indices[s][t][n]]->settings(Property::NORMAL_DIRECTION).back()->evaluate(Point(0, 0, 1));
					normal.normalize();
					double value[3] { normal.x, normal.y, normal.z };

					for (size_t dof = 0; dof < eDOFs.size(); dof++) {
						if (value[dof]) {
							constraints.B1[s].I_row_indices.push_back(row + IJVMatrixIndexing);
							constraints.B1[s].J_col_indices.push_back(nodes[indices[s][t][n]]->DOFIndex(s, dof) + IJVMatrixIndexing);
							constraints.B1[s].V_values.push_back(value[dof]);
							constraints.B1c[s].push_back(std::abs(value[dof]) * values[s][t][n]);
						}
					}
				} else {
					ESINFO(GLOBAL_ERROR) << "You have to set normal direction for elements with obstacle";
				}

			}
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
		constraints.B1clustersMap.push_back({ (eslocal)i, config::env::MPIrank });
	}

	constraints.block[Constraints::BLOCK::INEQUALITY_CONSTRAINTS] += globalIndicesSize;
	ESINFO(DETAILS) << "Lambdas with lower bounds in B1: " << globalIndicesSize;
}



