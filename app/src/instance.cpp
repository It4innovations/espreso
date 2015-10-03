
#include "instance.h"

Instance::Instance(const Configuration &configuration, int rank, int size)
	: _configuration(configuration), _rank(rank), _size(size), _localBoundaries(_mesh), _globalBoundaries(_mesh)
{
	if (_configuration.parameter(Configuration::MESH_FILE)->isSet()) {
		std::stringstream ssm;
		std::stringstream ssb;
		ssm << _configuration.value<std::string>(Configuration::MESH_FILE) << _rank << ".dat";
		ssb << _configuration.value<std::string>(Configuration::BOUNDARIES_FILE) << _rank << ".dat";

		std::cout << ssm.str() << "\n";

		_mesh.loadData(ssm.str().c_str());
		_mesh.partitiate(
				_configuration.value<eslocal>(Configuration::MESH_SUBDOMAINS),
				_configuration.value<eslocal>(Configuration::MESH_FIX_POINTS)
		);
		_localBoundaries.compute();
		_globalBoundaries.loadData(ssb.str().c_str());
		return;
	}

	if (_configuration.parameter(Configuration::ANSYS_DIR)->isSet()) {
		mesh::Ansys ansys(_configuration.value<std::string>(Configuration::ANSYS_DIR));
		ansys.coordinatesProperty(mesh::CP::DIRICHLET_X) = _configuration.value<std::string>(Configuration::ANSYS_DIRICHLET_X);
		ansys.coordinatesProperty(mesh::CP::DIRICHLET_Y) = _configuration.value<std::string>(Configuration::ANSYS_DIRICHLET_Y);
		ansys.coordinatesProperty(mesh::CP::DIRICHLET_Z) = _configuration.value<std::string>(Configuration::ANSYS_DIRICHLET_Z);
		ansys.coordinatesProperty(mesh::CP::FORCES_X) = _configuration.value<std::string>(Configuration::ANSYS_FORCES_X);
		ansys.coordinatesProperty(mesh::CP::FORCES_Y) = _configuration.value<std::string>(Configuration::ANSYS_FORCES_Y);
		ansys.coordinatesProperty(mesh::CP::FORCES_Z) = _configuration.value<std::string>(Configuration::ANSYS_FORCES_Z);

		_mesh.loadAnsys(
			ansys,
			_configuration.value<eslocal>(Configuration::MESH_SUBDOMAINS),
			_configuration.value<eslocal>(Configuration::MESH_FIX_POINTS)
		);
		_localBoundaries.compute();

		_mesh.computeCorners(
				_localBoundaries,
				_configuration.value<bool>(Configuration::MESH_CORNERS_NUMBER),
				_configuration.value<bool>(Configuration::MESH_CORNERS_IN_CORNER),
				_configuration.value<bool>(Configuration::MESH_CORNERS_IN_EDGES),
				_configuration.value<bool>(Configuration::MESH_CORNERS_IN_FACES)
		);
		return;
	}

//	esinput::Settings settings;
//	settings.clusters[0] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_X);
//	settings.clusters[1] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_Y);
//	settings.clusters[2] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_Z);
//	settings.subdomainsInCluster[0] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_X);
//	settings.subdomainsInCluster[1] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_Y);
//	settings.subdomainsInCluster[2] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_Z);
//	settings.elementsInSubdomain[0] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_X);
//	settings.elementsInSubdomain[1] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_Y);
//	settings.elementsInSubdomain[2] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_Z);
//
//	esinput::Generator *generator;
//	size_t cluster[3];
//
//	switch (_configuration.value<eslocal>(Configuration::PMCUBE_SHAPE)) {
//	case Configuration::CUBE: {
//		esinput::CubeSettings cubeSettings(settings);
//		switch (_configuration.value<eslocal>(Configuration::PMCUBE_ELEMENT_TYPE)) {
//		case Configuration::HEXA8: {
//			generator = new esinput::CubeGenerator<esinput::Hexahedron8>(cubeSettings);
//			break;
//		}
//		case Configuration::TETRA10: {
//			generator = new esinput::CubeGenerator<esinput::Tetrahedron10>(cubeSettings);
//			break;
//		}
//		case Configuration::TETRA4: {
//			generator = new esinput::CubeGenerator<esinput::Tetrahedron4>(cubeSettings);
//			break;
//		}
//		case Configuration::HEXA20: {
//			generator = new esinput::CubeGenerator<esinput::Hexahedron20>(cubeSettings);
//			break;
//		}
//		case Configuration::PRISMA6: {
//			generator = new esinput::CubeGenerator<esinput::Prisma6>(cubeSettings);
//			break;
//		}
//		case Configuration::PRISMA15: {
//			generator = new esinput::CubeGenerator<esinput::Prisma15>(cubeSettings);
//			break;
//		}
//		case Configuration::PYRAMID5: {
//			generator = new esinput::CubeGenerator<esinput::Pyramid5>(cubeSettings);
//			break;
//		}
//		case Configuration::PYRAMID13: {
//			generator = new esinput::CubeGenerator<esinput::Pyramid13>(cubeSettings);
//			break;
//		}
//		default:
//			std::cerr << "Unknown PermonCube element type.\n";
//			exit(EXIT_FAILURE);
//		}
//
//		if (generator->assumedProcessCount() != _size) {
//			if (_rank == 0) {
//				std::cerr << "Number of clusters(" << generator->assumedProcessCount();
//				std::cerr << ") does not accord number of MPI processes(" << _size << ").\n";
//			}
//			MPI_Finalize();
//			exit(EXIT_FAILURE);
//		}
//		generator->fillCluster(_rank, cluster);
//		generator->mesh(_mesh, cluster);
//		_localBoundaries.compute();
//
//		if (_configuration.value<bool>(Configuration::PMCUBE_FIX_ZERO_PLANES)) {
//			generator->setDirichlet(_mesh, cluster, esinput::CubeGeneratorOptions::DIRICHLET_ZERO_PLANES);
//		}
//		if (_configuration.value<bool>(Configuration::PMCUBE_FIX_BOTTOM)) {
//			generator->setDirichlet(_mesh, cluster, esinput::CubeGeneratorOptions::DIRICHLET_BOTTOM);
//		}
//
//		break;
//	}
//	case Configuration::SPHERE: {
//		esinput::SphereSettings sphereSettings(settings);
//		switch (_configuration.value<eslocal>(Configuration::PMCUBE_ELEMENT_TYPE)) {
//		case Configuration::HEXA8: {
//			generator = new esinput::SphereGenerator<esinput::Hexahedron8>(sphereSettings);
//			break;
//		}
//		case Configuration::TETRA10: {
//			generator = new esinput::SphereGenerator<esinput::Tetrahedron10>(sphereSettings);
//			break;
//		}
//		case Configuration::TETRA4: {
//			generator = new esinput::SphereGenerator<esinput::Tetrahedron4>(sphereSettings);
//			break;
//		}
//		case Configuration::HEXA20: {
//			generator = new esinput::SphereGenerator<esinput::Hexahedron20>(sphereSettings);
//			break;
//		}
//		case Configuration::PRISMA6: {
//			generator = new esinput::SphereGenerator<esinput::Prisma6>(sphereSettings);
//			break;
//		}
//		case Configuration::PRISMA15: {
//			generator = new esinput::SphereGenerator<esinput::Prisma15>(sphereSettings);
//			break;
//		}
//		case Configuration::PYRAMID5: {
//			generator = new esinput::SphereGenerator<esinput::Pyramid5>(sphereSettings);
//			break;
//		}
//		case Configuration::PYRAMID13: {
//			generator = new esinput::SphereGenerator<esinput::Pyramid13>(sphereSettings);
//			break;
//		}
//		default:
//			std::cerr << "Unknown PermonCube element type.\n";
//			exit(EXIT_FAILURE);
//		}
//
//		if (generator->assumedProcessCount() != _size) {
//			if (_rank == 0) {
//				std::cerr << "Number of clusters(" << generator->assumedProcessCount();
//				std::cerr << ") does not accord number of MPI processes(" << _size << ").\n";
//			}
//			MPI_Finalize();
//			exit(EXIT_FAILURE);
//		}
//		generator->fillCluster(_rank, cluster);
//		generator->mesh(_mesh, cluster);
//		_localBoundaries.compute();
//
//		generator->setDirichlet(_mesh, cluster, esinput::SphereGeneratorOptions::DIRICHLET_INNER_SURFACE);
//
//		break;
//	}
//	}
//
//	generator->fillGlobalBoundaries(_globalBoundaries, cluster);
//	generator->setFixPoints(_mesh);
//
//	size_t number[3] = {
//			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_X),
//			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_Y),
//			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_Z)
//	};
//	generator->setCorners(
//			_localBoundaries,
//			number,
//			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_CORNER),
//			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_EDGES),
//			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_FACES)
//	);
//
//	delete generator;
}



