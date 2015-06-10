
#include "hexahedron20.h"

using namespace permoncube;

eslocal Hexahedron20::subelements = Hexahedron8Subelements;

eslocal Hexahedron20::subnodes[3] = {
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
		Hexahedron20Subnodes,
};

std::vector<eslocal> Hexahedron20::_projection;

void Hexahedron20::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal hexa[20];
	hexa[0] = _projection[indices[2]];
	hexa[1] = _projection[indices[8]];
	hexa[2] = _projection[indices[6]];
	hexa[3] = _projection[indices[0]];
	hexa[4] = _projection[indices[20]];
	hexa[5] = _projection[indices[26]];
	hexa[6] = _projection[indices[24]];
	hexa[7] = _projection[indices[18]];

	hexa[8] = _projection[indices[8]];
	hexa[9] = _projection[indices[7]];
	hexa[10] = _projection[indices[3]];
	hexa[11] = _projection[indices[1]];
	hexa[12] = _projection[indices[23]];
	hexa[13] = _projection[indices[25]];
	hexa[14] = _projection[indices[21]];
	hexa[15] = _projection[indices[19]];
	hexa[16] = _projection[indices[11]];
	hexa[17] = _projection[indices[17]];
	hexa[18] = _projection[indices[15]];
	hexa[19] = _projection[indices[9]];
	mesh.pushElement(new mesh::Hexahedron20(hexa));
}

bool odd(esglobal x)
{
	return x % 2 == 1;
}

void Hexahedron20::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	eslocal cNodes[3];
	eslocal cElems;
	eslocal nCount;
	Utils<Hexahedron20>::clusterNodesCount(settings, cNodes);
	cElems = Utils<Hexahedron20>::clusterElementsCount(settings);

	nCount = cNodes[0] * cNodes[1] * cNodes[2]; // Full mesh
	nCount -= 4 * cElems;						// remove unused nodes in mesh
	for (int i = 0; i < 3; i++) {				// remove unused nodes from the surface
		nCount -=
			settings.subdomainsInCluster[i] * settings.elementsInSubdomain[i] *
			settings.subdomainsInCluster[(i + 1) % 3] * settings.elementsInSubdomain[(i + 1) % 3];
	}

	mesh.coordinates().reserve(nCount);
	_projection.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	esglobal global = 0;
	eslocal local = 0;
	esglobal s[3], e[3];
	double step[3];
	for (eslocal i = 0; i < 3; i++) {
		s[i] = (cNodes[i] - 1) * cluster[i];
		e[i] = (cNodes[i] - 1) * (cluster[i] + 1);
	}
	for (eslocal i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (cNodes[i] - 1);
	}

	esglobal gNodes[3];
	Utils<Hexahedron20>::globalNodesCount(settings, gNodes);

	for (esglobal z = 0; z < gNodes[2]; z++) {
		for (esglobal y = 0; y < gNodes[1]; y++) {
			for (esglobal x = 0; x < gNodes[0]; x++) {
				if (s[2] <= z && z <= e[2] && s[1] <= y && y <= e[1] && s[0] <= x && x <= e[0]) {
					_projection.push_back(local);
					if ((odd(x) && odd(y)) || (odd(y) && odd(z)) || (odd(x) && odd(z))) {
						continue;
					}
					coordinates.add(mesh::Point(x * step[0], y * step[1], z * step[2]), local, global);
					local++;
				}
				global++;
			}
		}
	}

	if (local != nCount) {
		std::cerr << "Permoncube: internal error while adding coordinates of Hexahedron20\n";
		exit(EXIT_FAILURE);
	}
}

void Hexahedron20::fixZeroPlanes(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	eslocal nodes[3];
	Utils<Hexahedron20>::clusterNodesCount(settings, nodes);

	if (cluster[0] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal y = 0; y < nodes[1]; y++) {
				if (!(odd(z) && odd(y))) {
					dirichlet_x[_projection[index]] = 0;
				}
				index += nodes[0];
			}
		}
	}
	if (cluster[1] == 0) {
		eslocal index = 0;
		for (eslocal z = 0; z < nodes[2]; z++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (!(odd(z) && odd(x))) {
					dirichlet_y[_projection[index]] = 0;
				}
				index++;
			}
			index = (z + 1) * nodes[1] * nodes[0];
		}
	}
	if (cluster[2] == 0) {
		eslocal index = 0;
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				if (!(odd(y) && odd(x))) {
					dirichlet_z[_projection[index]] = 0;
				}
				index++;
			}
		}
	}
}

void Hexahedron20::fixBottom(
		const permoncube::Settings &settings,
		std::map<eslocal, double> &dirichlet_x,
		std::map<eslocal, double> &dirichlet_y,
		std::map<eslocal, double> &dirichlet_z,
		const size_t cluster[])
{
	if (cluster[2] > 0) {
		return;
	}
	eslocal nodes[3];
	Utils<Hexahedron20>::clusterNodesCount(settings, nodes);

	eslocal index = 0;
	for (eslocal y = 0; y < nodes[1]; y++) {
		for (eslocal x = 0; x < nodes[0]; x++) {
			if (odd(y) && odd(x)) {
				continue;
			}
			dirichlet_z[index] = 0;
			dirichlet_y[index] = 0;
			dirichlet_x[index] = 0;
			index++;
		}
	}
}

void Hexahedron20::fillGlobalBoundaries(
		const permoncube::Settings &settings,
		mesh::Boundaries &boundaries)
{
	esglobal gNodes[3];
	Utils<Hexahedron20>::globalNodesCount(settings, gNodes);
	eslocal cNodes[3];
	Utils<Hexahedron20>::clusterNodesCount(settings, cNodes);
	boundaries.resize(gNodes[0] * gNodes[1] * gNodes[2]);

	eslocal cluster[3];
	bool border[3];
	eslocal cIndex;
	esglobal gIndex = 0;

	for (esglobal z = 0; z < gNodes[2]; z++) {
		cluster[2] = z ? (z - 1) / ( cNodes[2] - 1) : 0;
		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
		for (esglobal y = 0; y < gNodes[1]; y++) {
			cluster[1] = y ? (y - 1) / ( cNodes[1] - 1) : 0;
			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
			for (esglobal x = 0; x < gNodes[0]; x++) {
				if ((odd(x) && odd(y)) || (odd(y) && odd(z)) || (odd(x) && odd(z))) {
					continue;
				}
				cluster[0] = x ? (x - 1) / ( cNodes[0] - 1) : 0;
				border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
				cIndex = cluster[0] + cluster[1] * settings.clusters[0] + cluster[2] * settings.clusters[0] * settings.clusters[1];
				boundaries[gIndex].insert(cIndex);
				for (int i = 0; i < 8; i++) {
					eslocal tmp = cIndex;
					if (border[0] && (i & 1)) {
						tmp += 1;
					}
					if (border[1] && (i & 2)) {
						tmp += settings.clusters[0];
					}
					if (border[2] && (i & 4)) {
						tmp += settings.clusters[0] * settings.clusters[1];
					}
					boundaries[gIndex].insert(tmp);
				}
				gIndex++;
			}
		}
	}
}


