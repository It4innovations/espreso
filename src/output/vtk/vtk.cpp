
#include "vtk.h"

#include "../../basis/utilities/utils.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/elements/element.h"

#include "../../assembler/constraints/constraints.h"
#include "../../assembler/instance.h"
#include "../../solver/generic/SparseMatrix.h"

#include "../../configuration/environment.h"
#include "../../configuration/output.h"

namespace espreso {
namespace store {

void VTK::computeCenters()
{
	_sCenters.resize(_mesh->parts());
	for (size_t p = 0; p < _mesh->coordinates().parts(); p++) {
		for (size_t i = 0; i < _mesh->coordinates().localSize(p); i++) {
			_sCenters[p] += _mesh->coordinates().get(i, p);
		}
		_sCenters[p] /= _mesh->coordinates().localSize(p);
	}

	for (size_t i = 0; i < _mesh->coordinates().clusterSize(); i++) {
		_cCenter += _mesh->coordinates()[i];
	}
	_cCenter /= _mesh->coordinates().clusterSize();
}

Point VTK::shrink(const Point &p, size_t part) const
{
	Point x = p;
	x = _sCenters[part] + (x - _sCenters[part]) * _output.domain_shrink_ratio;
	x = _cCenter + (x - _cCenter) * _output.cluster_shrink_ratio;
	return x;
}

Point VTK::shrink(const Point &p, const Point &sCenter, const Point &cCenter) const
{
	Point x = p;
	x = sCenter + (x - sCenter) * _output.domain_shrink_ratio;
	x = cCenter + (x - cCenter) * _output.cluster_shrink_ratio;
	return x;
}

void VTK::storeGeometry(size_t timeStep)
{
	coordinates();
	cells(ElementType::ELEMENTS);
}

void VTK::storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType)
{
	std::vector<std::vector<int> > selection(_mesh->parts());
	std::vector<std::vector<double> > values(_mesh->parts());

	switch (eType) {
	case ElementType::ELEMENTS:
		for (size_t e = 0; e < _mesh->elements().size(); e++) {
			const Element *element = _mesh->elements()[e];
			for (size_t p = 0; p < properties.size(); p++) {
				size_t domain = element->domains()[0];
				double value = 0;
				for (size_t n = 0; n < element->nodes(); n++) {
					value += element->getProperty(properties[p], n, 0, 0);
				}
				values[domain].push_back(value / element->nodes());
			}
		}
		data(name + "FixedValue", properties.size(), values, eType);
		break;
	case ElementType::FACES:
		ESINFO(GLOBAL_ERROR) << "Implement store properties for faces";
		break;
	case ElementType::EDGES:
		ESINFO(GLOBAL_ERROR) << "Implement store properties for edges";
		break;
	case ElementType::NODES:
		for (size_t n = 0; n < _mesh->nodes().size(); n++) {
			const Element *node = _mesh->nodes()[n];
			for (size_t p = 0; p < properties.size(); p++) {
				for (size_t d = 0; d < node->domains().size(); d++) {
					size_t domain = node->domains()[d];
					selection[domain].push_back(node->hasProperty(properties[p], 0) ? 1 : 0);
					values[domain].push_back(node->getProperty(properties[p], 0, 0, 0));
				}
			}
		}
		data(name + "IsSet", properties.size(), selection, eType);
		data(name + "FixedValue", properties.size(), values, eType);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown element type";
	}
}

void VTK::storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType)
{
	data(name, dimension, values, eType);
}

void VTK::mesh(const OutputConfiguration &output, const Mesh &mesh, const std::string &path, ElementType eType)
{
	VTK vtk(output, mesh, path);

	vtk.coordinates();
	vtk.cells(eType);
	vtk.finalize();
}

void VTK::fixPoints(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	std::vector<Element*> fixPoints;
	for (size_t p = 0; p < mesh.parts(); p++) {
		fixPoints.insert(fixPoints.end(), mesh.fixPoints(p).begin(), mesh.fixPoints(p).end());
	}

	std::sort(fixPoints.begin(), fixPoints.end());
	Esutils::removeDuplicity(fixPoints);

	VTK vtk(output, mesh, path);

	vtk.nodes(fixPoints);
	vtk.finalize();
}

void VTK::corners(const OutputConfiguration &output, const Mesh &mesh, const std::string &path)
{
	VTK vtk(output, mesh, path);

	vtk.nodes(mesh.corners());
	vtk.finalize();
}

void VTK::gluing(const OutputConfiguration &output, const Mesh &mesh, const Constraints &constraints, const std::string &path, size_t dofs)
{
	VTK vtk(output, mesh, path);

	std::vector<Point> cCenter(environment->MPIsize);
	std::vector<Point> sCenters(environment->MPIsize * mesh.parts());

	MPI_Allgather(&vtk._cCenter, sizeof(Point), MPI_BYTE, cCenter.data(), sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);
	MPI_Allgather(vtk._sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);

	std::vector<size_t> sOffset(mesh.parts());
	std::vector<size_t> eOffset(mesh.parts());

	std::vector<std::vector<std::vector<std::pair<esglobal, Element*> > > > DOF2e(mesh.parts(), std::vector<std::vector<std::pair<esglobal, Element*> > >(dofs));

	#pragma omp parallel for
	for (size_t p = 0; p < mesh.parts(); p++) {
		eOffset[p] += std::lower_bound(constraints.B1[p].I_row_indices.begin(), constraints.B1[p].I_row_indices.end(), constraints.block[Constraints::DIRICHLET] + 1) - constraints.B1[p].I_row_indices.begin();
	}
	for (size_t n = 0; n < mesh.nodes().size(); n++) {
		for (size_t d = 0; d < mesh.nodes()[n]->domains().size(); d++) {
			size_t p = mesh.nodes()[n]->domains()[d];
			for (size_t dof = 0; dof < dofs; dof++) {
				DOF2e[p][dof].push_back(std::make_pair(mesh.nodes()[n]->DOFIndex(p, dof) + 1, mesh.nodes()[n]));
			}
		}
	}
	#pragma omp parallel for
	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t dof = 0; dof < dofs; dof++) {
			std::sort(DOF2e[p][dof].begin(), DOF2e[p][dof].end());
		}
	}

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "dirichlet";
		VTK dirichlet(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());

		#pragma omp parallel for
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), constraints.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == constraints.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
				}
			}
		}

		dirichlet.nodes(dnodes);
		dirichlet.finalize();
	}

	size_t maxLambda = constraints.block[Constraints::DIRICHLET] + constraints.block[Constraints::EQUALITY_CONSTRAINTS] + 1;
	for (size_t p = 0; p < mesh.parts(); p++) {
		sOffset[p] = eOffset[p];
		eOffset[p] = std::lower_bound(constraints.B1[p].I_row_indices.begin() + eOffset[p], constraints.B1[p].I_row_indices.end(), maxLambda) - constraints.B1[p].I_row_indices.begin();
	}

	auto getDomain = [&] (esglobal lambda, size_t exclude) -> size_t {
		for (size_t p = 0; p < mesh.parts(); p++) {
			if (p == exclude) {
				continue;
			}
			auto it = std::lower_bound(constraints.B1subdomainsMap[p].begin(), constraints.B1subdomainsMap[p].end(), lambda - 1);
			if (it != constraints.B1subdomainsMap[p].end() && *it == lambda - 1) {
				return p;
			}
		}
		ESINFO(ERROR) << "Internal error: Broken exporting of gluing matrices";
		return 0;
	};

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "gluing";
		VTK gluing(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());
		std::vector<std::vector<size_t> > indices(mesh.parts());

		std::vector<std::vector<std::vector<std::pair<esglobal, eslocal> > > > CDneighbour(mesh.parts(), std::vector<std::vector<std::pair<esglobal, eslocal> > >(environment->MPIsize));
		std::vector<std::vector<std::pair<eslocal, eslocal> > > CDindex(mesh.parts()); // cluster x domain index
		std::vector<std::vector<std::pair<esglobal, eslocal> > > sBuffer(environment->MPIsize);
		std::vector<std::vector<std::pair<esglobal, eslocal> > > rBuffer(environment->MPIsize);

		#pragma omp parallel for
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), constraints.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == constraints.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
					indices[p].push_back(i);

					auto it = std::lower_bound(constraints.B1clustersMap.begin(), constraints.B1clustersMap.end(), constraints.B1[p].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
						return v[0] < i;
					});
					if (it->size() == 2) { // local gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), (eslocal)getDomain(constraints.B1[p].I_row_indices[i], p)));
					} else { // global gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), -1));
						CDneighbour[p][it->back()].push_back(std::make_pair(constraints.B1[p].I_row_indices[i], p));
					}

				}
			}
		}

		#pragma omp parallel for
		for (int c = 0; c < environment->MPIsize; c++) {
			for (size_t p = 0; p < mesh.parts(); p++) {
				sBuffer[c].insert(sBuffer[c].end(), CDneighbour[p][c].begin(), CDneighbour[p][c].end());
			}
			std::sort(sBuffer[c].begin(), sBuffer[c].end(), [] (const std::pair<esglobal, eslocal> &p1, const std::pair<esglobal, eslocal> &p2) {
				return p1.first < p2.first;
			});
		}

		std::vector<MPI_Request> req(2 * environment->MPIsize);
		for (int n = 0; n < environment->MPIsize; n++) {
			rBuffer[n].resize(sBuffer[n].size());
			MPI_Isend(sBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * sBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n);
			MPI_Irecv(rBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * rBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n + 1);
		}

		MPI_Waitall(2 * environment->MPIsize, req.data(), MPI_STATUSES_IGNORE);

		for (size_t p = 0; p < mesh.parts(); p++) {
			std::vector<int> offsets(environment->MPIsize);
			for (size_t i = 0; i < CDindex[p].size(); i++) {
				if (CDindex[p][i].second == -1) {
					int n = CDindex[p][i].first;
					esglobal lambda = CDneighbour[p][n][offsets[n]++].first;
					auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), lambda, [] (const std::pair<esglobal, eslocal> &p, esglobal l) { return p.first < l; });
					if (it != rBuffer[n].end() && it->first == lambda) {
						CDindex[p][i].second = it->second;
					} else {
						ESINFO(ERROR) << "Internal error: gluing matrices are weird.";
					}
				}
			}
		}

		size_t cSize = 0;
		for (size_t p = 0; p < dnodes.size(); p++) {
			cSize += dnodes[p].size();
		}

		gluing.lambdas(dnodes, [&] (const Point &point, size_t part, size_t index, bool start) -> Point {
			if (start) {
				return vtk.shrink(point, part);
			} else {
				Point p1 = vtk.shrink(point, part);
				Point p2 = vtk.shrink(point, sCenters[CDindex[part][index].first * mesh.parts() + CDindex[part][index].second], cCenter[CDindex[part][index].first]);
				return p1 + (p2 - p1) * constraints.B1duplicity[part][indices[part][index]];
			}
		});

		std::vector<std::vector<double> > values(mesh.parts());
		for (size_t p = 0; p < dnodes.size(); p++) {
			for (size_t n = 0; n < dnodes[p].size(); n++) {
				values[p].push_back(constraints.B1[p].V_values[indices[p][n]]);
			}
		}

		gluing.data("values", 1, values, ElementType::ELEMENTS);
		gluing.finalize();
	}
}

void VTK::gluing(const OutputConfiguration &output, const Mesh &mesh, const Instance &instance, const std::string &path, size_t dofs)
{
	VTK vtk(output, mesh, path);

	std::vector<Point> cCenter(environment->MPIsize);
	std::vector<Point> sCenters(environment->MPIsize * mesh.parts());

	MPI_Allgather(&vtk._cCenter, sizeof(Point), MPI_BYTE, cCenter.data(), sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);
	MPI_Allgather(vtk._sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, sCenters.data(), mesh.parts() * sizeof(Point), MPI_BYTE, MPI_COMM_WORLD);

	std::vector<size_t> sOffset(mesh.parts());
	std::vector<size_t> eOffset(mesh.parts());

	std::vector<std::vector<std::vector<std::pair<esglobal, Element*> > > > DOF2e(mesh.parts(), std::vector<std::vector<std::pair<esglobal, Element*> > >(dofs));

	#pragma omp parallel for
	for (size_t p = 0; p < mesh.parts(); p++) {
		eOffset[p] += std::lower_bound(instance.B1[p].I_row_indices.begin(), instance.B1[p].I_row_indices.end(), instance.block[Constraints::DIRICHLET] + 1) - instance.B1[p].I_row_indices.begin();
	}
	for (size_t n = 0; n < mesh.nodes().size(); n++) {
		for (size_t d = 0; d < mesh.nodes()[n]->domains().size(); d++) {
			size_t p = mesh.nodes()[n]->domains()[d];
			for (size_t dof = 0; dof < dofs; dof++) {
				DOF2e[p][dof].push_back(std::make_pair(mesh.nodes()[n]->DOFIndex(p, dof) + 1, mesh.nodes()[n]));
			}
		}
	}
	#pragma omp parallel for
	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t dof = 0; dof < dofs; dof++) {
			std::sort(DOF2e[p][dof].begin(), DOF2e[p][dof].end());
		}
	}

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "dirichlet";
		VTK dirichlet(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());

		#pragma omp parallel for
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), instance.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == instance.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
				}
			}
		}

		dirichlet.nodes(dnodes);
		dirichlet.finalize();
	}

	size_t maxLambda = instance.block[Constraints::DIRICHLET] + instance.block[Constraints::EQUALITY_CONSTRAINTS] + 1;
	for (size_t p = 0; p < mesh.parts(); p++) {
		sOffset[p] = eOffset[p];
		eOffset[p] = std::lower_bound(instance.B1[p].I_row_indices.begin() + eOffset[p], instance.B1[p].I_row_indices.end(), maxLambda) - instance.B1[p].I_row_indices.begin();
	}

	auto getDomain = [&] (esglobal lambda, size_t exclude) -> size_t {
		for (size_t p = 0; p < mesh.parts(); p++) {
			if (p == exclude) {
				continue;
			}
			auto it = std::lower_bound(instance.B1subdomainsMap[p].begin(), instance.B1subdomainsMap[p].end(), lambda - 1);
			if (it != instance.B1subdomainsMap[p].end() && *it == lambda - 1) {
				return p;
			}
		}
		ESINFO(ERROR) << "Internal error: Broken exporting of gluing matrices";
		return 0;
	};

	for (size_t dof = 0; dof < dofs; dof++) {
		std::stringstream ss;
		ss << dof << "gluing";
		VTK gluing(output, mesh, ss.str());

		std::vector<std::vector<eslocal> > dnodes(mesh.parts());
		std::vector<std::vector<size_t> > indices(mesh.parts());

		std::vector<std::vector<std::vector<std::pair<esglobal, eslocal> > > > CDneighbour(mesh.parts(), std::vector<std::vector<std::pair<esglobal, eslocal> > >(environment->MPIsize));
		std::vector<std::vector<std::pair<eslocal, eslocal> > > CDindex(mesh.parts()); // cluster x domain index
		std::vector<std::vector<std::pair<esglobal, eslocal> > > sBuffer(environment->MPIsize);
		std::vector<std::vector<std::pair<esglobal, eslocal> > > rBuffer(environment->MPIsize);

		#pragma omp parallel for
		for (size_t p = 0; p < mesh.parts(); p++) {
			for (size_t i = sOffset[p]; i < eOffset[p]; i++) {
				auto it = std::lower_bound(DOF2e[p][dof].begin(), DOF2e[p][dof].end(), instance.B1[p].J_col_indices[i], [] (const std::pair<esglobal, Element*> &pair, esglobal index) {
					return pair.first < index;
				});
				if (it != DOF2e[p][dof].end() && it->first == instance.B1[p].J_col_indices[i]) {
					dnodes[p].push_back(it->second->node(0));
					indices[p].push_back(i);

					auto it = std::lower_bound(instance.B1clustersMap.begin(), instance.B1clustersMap.end(), instance.B1[p].I_row_indices[i] - 1, [&] (const std::vector<esglobal> &v, esglobal i) {
						return v[0] < i;
					});
					if (it->size() == 2) { // local gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), (eslocal)getDomain(instance.B1[p].I_row_indices[i], p)));
					} else { // global gluing
						CDindex[p].push_back(std::make_pair((eslocal)it->back(), -1));
						CDneighbour[p][it->back()].push_back(std::make_pair(instance.B1[p].I_row_indices[i], p));
					}

				}
			}
		}

		#pragma omp parallel for
		for (int c = 0; c < environment->MPIsize; c++) {
			for (size_t p = 0; p < mesh.parts(); p++) {
				sBuffer[c].insert(sBuffer[c].end(), CDneighbour[p][c].begin(), CDneighbour[p][c].end());
			}
			std::sort(sBuffer[c].begin(), sBuffer[c].end(), [] (const std::pair<esglobal, eslocal> &p1, const std::pair<esglobal, eslocal> &p2) {
				return p1.first < p2.first;
			});
		}

		std::vector<MPI_Request> req(2 * environment->MPIsize);
		for (int n = 0; n < environment->MPIsize; n++) {
			rBuffer[n].resize(sBuffer[n].size());
			MPI_Isend(sBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * sBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n);
			MPI_Irecv(rBuffer[n].data(), sizeof(std::pair<esglobal, eslocal>) * rBuffer[n].size(), MPI_BYTE, n, 1, MPI_COMM_WORLD, req.data() + 2 * n + 1);
		}

		MPI_Waitall(2 * environment->MPIsize, req.data(), MPI_STATUSES_IGNORE);

		for (size_t p = 0; p < mesh.parts(); p++) {
			std::vector<int> offsets(environment->MPIsize);
			for (size_t i = 0; i < CDindex[p].size(); i++) {
				if (CDindex[p][i].second == -1) {
					int n = CDindex[p][i].first;
					esglobal lambda = CDneighbour[p][n][offsets[n]++].first;
					auto it = std::lower_bound(rBuffer[n].begin(), rBuffer[n].end(), lambda, [] (const std::pair<esglobal, eslocal> &p, esglobal l) { return p.first < l; });
					if (it != rBuffer[n].end() && it->first == lambda) {
						CDindex[p][i].second = it->second;
					} else {
						ESINFO(ERROR) << "Internal error: gluing matrices are weird.";
					}
				}
			}
		}

		size_t cSize = 0;
		for (size_t p = 0; p < dnodes.size(); p++) {
			cSize += dnodes[p].size();
		}

		gluing.lambdas(dnodes, [&] (const Point &point, size_t part, size_t index, bool start) -> Point {
			if (start) {
				return vtk.shrink(point, part);
			} else {
				Point p1 = vtk.shrink(point, part);
				Point p2 = vtk.shrink(point, sCenters[CDindex[part][index].first * mesh.parts() + CDindex[part][index].second], cCenter[CDindex[part][index].first]);
				return p1 + (p2 - p1) * instance.B1duplicity[part][indices[part][index]];
			}
		});

		std::vector<std::vector<double> > values(mesh.parts());
		for (size_t p = 0; p < dnodes.size(); p++) {
			for (size_t n = 0; n < dnodes[p].size(); n++) {
				values[p].push_back(instance.B1[p].V_values[indices[p][n]]);
			}
		}

		gluing.data("values", 1, values, ElementType::ELEMENTS);
		gluing.finalize();
	}
}

}
}


