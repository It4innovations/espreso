
#ifdef HAVE_GMSH
#include "gmsh.h"
#endif

#include "w.gmsh.h"
#include "basis/utilities/parser.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "input/meshbuilder.h"
#include "mesh/element.h"

using namespace espreso;

bool GMSH::islinked()
{
#ifdef HAVE_GMSH
	return true;
#else
	return false;
#endif
}

void GMSH::generate(MeshData &mesh)
{
#ifndef HAVE_GMSH
	eslog::globalerror("ESPRESO run-time error: cannot call GMSH library (the library is not linked).\n");
#else
	if (info::mpi::rank) {
		// GMSH is sequential only
		return;
	}

	const GMSHConfiguration &configuration = info::ecf->input.generation.gmsh_options;
	bool isvalid = false;
	gmsh::initialize();
	gmsh::model::add("model");
	if (
			StringCompare::caseSensitiveSuffix(info::ecf->input.path, "stp") ||
			StringCompare::caseSensitiveSuffix(info::ecf->input.path, "step") ||
			StringCompare::caseSensitiveSuffix(info::ecf->input.path, "igs")) {
		std::vector<std::pair<int, int> > v;
		try {
			gmsh::model::occ::importShapes(info::ecf->input.path, v);
		} catch(...) {
			eslog::error("ESPRESO run-time error: GMSH cannot load the file.\n");
		}
		gmsh::model::occ::synchronize();

		isvalid = true;
	}

	if (StringCompare::caseSensitiveSuffix(info::ecf->input.path, "stl")) {
		try {
			gmsh::merge(info::ecf->input.path);
		} catch(...) {
			eslog::error("ESPRESO run-time error: GMSH cannot load the file.\n");
		}

		double angle = configuration.stl_angle;
		bool includeBoundary = true;
		bool forceParametrizablePatches = false;
		gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches);

		gmsh::model::mesh::createGeometry();
		std::vector<std::pair<int, int> > s;
		gmsh::model::getEntities(s, 2);
		std::vector<int> sl;
		for(auto surf: s) sl.push_back(surf.second);
		int l = gmsh::model::geo::addSurfaceLoop(sl);
		gmsh::model::geo::addVolume({l});
		gmsh::model::geo::synchronize();
		int f = gmsh::model::mesh::field::add("MathEval");
		gmsh::model::mesh::field::setString(f, "F", std::to_string(configuration.stl_precision));
		gmsh::model::mesh::field::setAsBackgroundMesh(f);

		isvalid = true;
	}

	if (!isvalid) {
		eslog::error("MESIO run-time error: GMSH does not support provided format. Supported formats are: 'stp', 'step', 'igs', 'stl'.\n");
	}

	gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", configuration.characteristic_length.extern_from_boundary);
	gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", configuration.characteristic_length.from_points);
	gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", configuration.characteristic_length.from_curvature);
	gmsh::option::setNumber("Mesh.CharacteristicLengthMin", configuration.characteristic_length.min);
	gmsh::option::setNumber("Mesh.CharacteristicLengthMax", configuration.characteristic_length.max);
	gmsh::option::setNumber("Mesh.Algorithm3D", configuration.algorithm3D);
	gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", configuration.subdivisionAlgorithm);
	gmsh::option::setNumber("Mesh.Optimize", configuration.optimize);
	gmsh::option::setNumber("Mesh.OptimizeNetgen", 1);
	gmsh::option::setNumber("Mesh.MaxNumThreads3D", info::env::OMP_NUM_THREADS);
	gmsh::option::setNumber("Mesh.MaxNumThreads2D", info::env::OMP_NUM_THREADS);
	gmsh::option::setNumber("Mesh.MaxNumThreads1D", info::env::OMP_NUM_THREADS);

	gmsh::model::mesh::generate(3);

	{ // insert nodes
		std::vector<size_t> ids;
		std::vector<double> coords, params;
		gmsh::model::mesh::getNodes(ids, coords, params, -1, -1, false, false);

		mesh.nIDs.reserve(ids.size());
		mesh.coordinates.reserve(ids.size());
		for (size_t i = 0; i < ids.size(); ++i) {
			mesh.nIDs.push_back(ids[i] - 1);
			mesh.coordinates.push_back(Point(coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]));
		}
	}

	std::vector<std::pair<int, int> > entities;
	gmsh::model::getEntities(entities);
	for (auto entity = entities.begin(); entity != entities.end(); ++entity) { // insert elements
		int dim = entity->first;
		int tag = entity->second;
		std::vector<int> types;
		std::vector<std::vector<size_t> > ids, enodes;
		gmsh::model::mesh::getElements(types, ids, enodes, dim, tag);

		for (size_t t = 0; t < types.size(); ++t) {
			for (size_t i = 0; i < ids[t].size(); ++i) {
				ids[t][i] -= 1;
			}
			for (size_t i = 0; i < enodes[t].size(); ++i) {
				enodes[t][i] -= 1;
			}
		}

		std::string name;
		gmsh::model::getEntityName(dim, tag, name);

		if (dim == 0) {
			if (name.size() == 0) {
				name = "NODES_" + std::to_string(tag);
			}
			for (size_t t = 0; t < types.size(); ++t) {
				mesh.nregions[name].insert(mesh.nregions[name].end(), ids[t].begin(), ids[t].end());
			}
		} else if (dim > 1) { // there is probably bug in ESPRESO during loading 1D elements
			for (size_t t = 0; t < types.size(); ++t) {
				switch (types[t]) {

				// 1D linear
				case 1:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::LINE2);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 2);
					break;

				// 1D quadratic
				case 8:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::LINE3);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 3);
					break;

				// 2D linear
				case 2:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::TRIANGLE3);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 3);
					break;
				case 3:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::SQUARE4);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 4);
					break;

				// 2D quadratic
				case 9:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::TRIANGLE6);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 6);
					break;
				case 16:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::SQUARE8);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 8);
					break;

				// 3D linear
				case 4:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::TETRA4);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 4);
					break;
				case 5:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::HEXA8);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 8);
					break;
				case 6:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::PRISMA6);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 6);
					break;
				case 7:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::PYRAMID5);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 5);
					break;

				// 3D quadratic
				case 11:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::TETRA10);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 11);
					break;
				case 17:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::HEXA20);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 20);
					break;
				case 18:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::PRISMA15);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 15);
					break;
				case 19:
					mesh.etype.insert(mesh.etype.end(), ids[t].size(), (int)Element::CODE::PYRAMID13);
					mesh.esize.insert(mesh.esize.end(), ids[t].size(), 13);
					break;
				default:
					eslog::error("ESPRESO internal error: GMSH returns unsupported element type: %d\n", types[t]);
				}
				mesh.eIDs.insert(mesh.eIDs.end(), ids[t].begin(), ids[t].end());
				mesh.enodes.insert(mesh.enodes.end(), enodes[t].begin(), enodes[t].end());
			}

			if (name.size() == 0) {
				if (dim == 3) {
					name = "ELEMENTS_" + std::to_string(tag);
				} else {
					name = "BOUNDARY_" + std::to_string(tag);
				}
			}
			for (size_t t = 0; t < types.size(); ++t) {
				mesh.eregions[name].insert(mesh.eregions[name].end(), ids[t].begin(), ids[t].end());
			}
		}
	}
#endif
}
