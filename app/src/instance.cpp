
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

	permoncube::Settings settings;

	settings.clusters[0] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_X);
	settings.clusters[1] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_Y);
	settings.clusters[2] = _configuration.value<eslocal>(Configuration::PMCUBE_CLUSTERS_Z);
	settings.subdomainsInCluster[0] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_X);
	settings.subdomainsInCluster[1] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_Y);
	settings.subdomainsInCluster[2] = _configuration.value<eslocal>(Configuration::PMCUBE_SUBDOMAINS_Z);
	settings.elementsInSubdomain[0] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_X);
	settings.elementsInSubdomain[1] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_Y);
	settings.elementsInSubdomain[2] = _configuration.value<eslocal>(Configuration::PMCUBE_ELEMENTS_Z);

	if (settings.clusters[0] * settings.clusters[1] * settings.clusters[2] != _size) {
		std::cerr << "Number of clusters(";
		std::cerr << settings.clusters[0] * settings.clusters[1] * settings.clusters[2];
		std::cerr << ") does not accord number of MPI processes(";
		std::cerr << _size;
		std::cerr << ").\n";
		exit(EXIT_SUCCESS);
	}

	permoncube::Generator *generator;

	switch (_configuration.value<eslocal>(Configuration::PMCUBE_ELEMENT_TYPE)) {
	case Configuration::HEXA8: {
		generator = new permoncube::ElementGenerator<permoncube::Hexahedron8>(settings);
		break;
	}
	case Configuration::TETRA10: {
		generator = new permoncube::ElementGenerator<permoncube::Tetrahedron10>(settings);
		break;
	}
	case Configuration::TETRA4: {
		generator = new permoncube::ElementGenerator<permoncube::Tetrahedron4>(settings);
		break;
	}
	case Configuration::HEXA20: {
		generator = new permoncube::ElementGenerator<permoncube::Hexahedron20>(settings);
		break;
	}
	case Configuration::PRISMA6: {
		generator = new permoncube::ElementGenerator<permoncube::Prisma6>(settings);
		break;
	}
	case Configuration::PRISMA15: {
		generator = new permoncube::ElementGenerator<permoncube::Prisma15>(settings);
		break;
	}
	case Configuration::PYRAMID5: {
		generator = new permoncube::ElementGenerator<permoncube::Pyramid5>(settings);
		break;
	}
	case Configuration::PYRAMID13: {
		generator = new permoncube::ElementGenerator<permoncube::Pyramid13>(settings);
		break;
	}
	default:
		std::cerr << "Unknown PermonCube element type.\n";
		exit(EXIT_FAILURE);
	}

	size_t cluster[3];
	cluster[0] = _rank % settings.clusters[0];
	cluster[1] = (_rank / settings.clusters[0]) % settings.clusters[1];
	cluster[2] = _rank / (settings.clusters[0] * settings.clusters[1]);

	generator->mesh(_mesh, cluster);
	_localBoundaries.compute();
	if (_configuration.value<bool>(Configuration::PMCUBE_FIX_ZERO_PLANES)) {
		generator->fixZeroPlanes(_mesh, cluster);
	}
	if (_configuration.value<bool>(Configuration::PMCUBE_FIX_BOTTOM)) {
		generator->fixBottom(_mesh, cluster);
	}
	generator->fillGlobalBoundaries(_globalBoundaries, cluster);
	generator->setFixPoints(_mesh);

	size_t number[3] = {
			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_X),
			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_Y),
			_configuration.value<eslocal>(Configuration::PMCUBE_CORNERS_Z)
	};
	generator->setCorners(
			_localBoundaries,
			number,
			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_CORNER),
			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_EDGES),
			_configuration.value<bool>(Configuration::MESH_CORNERS_IN_FACES)
	);

	delete generator;
}



