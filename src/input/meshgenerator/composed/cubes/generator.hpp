
#include "generator.h"

#include "esconfig.h"

namespace espreso {
namespace input {

static void setCluster(size_t cluster[], CubesSettings &settings, size_t &cubeIndex)
{
	size_t clusters0 = settings.cube[0].clusters[0] * settings.cube[0].clusters[1] * settings.cube[0].clusters[2];
	size_t clusters1 = settings.cube[1].clusters[0] * settings.cube[1].clusters[1] * settings.cube[1].clusters[2];

	size_t rank = settings.cube[0].index;
	size_t size = settings.cube[0].size;

	if (clusters0 + clusters1 != size) {
		ESINFO(espreso::GLOBAL_ERROR)
				<< "The number of clusters(" << clusters0 + clusters1
				<< ") does not accord the number of MPI processes("
				<< settings.cube[0].size + settings.cube[1].size << ").";
	}

	cubeIndex = rank < clusters0 ? 0 : 1;

	settings.cube[0].size = clusters0;
	settings.cube[1].size = clusters1;
	settings.cube[0].index = settings.cube[0].index % settings.cube[0].size;
	settings.cube[1].index = settings.cube[1].index % settings.cube[1].size;
	settings.cube[1].clusterOffset = settings.cube[0].size;
}

template<class TElement>
CubesGenerator<TElement>::CubesGenerator(Mesh &mesh, CubesSettings &settings)
: Loader(mesh), _settings(settings)
{
	setCluster(_cluster, _settings, _cubeIndex);
	switch(_cubeIndex) {
	case 0:
		_loader = new CubeGenerator<TElement>(mesh, _settings.cube[0]);
		break;
	case 1:
		_loader = new CubeGenerator<TElement>(mesh, _settings.cube[1]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Wrong cube index";
	}
}

template<class TElement>
void CubesGenerator<TElement>::points(Coordinates &coordinates)
{
	_loader->points(coordinates);
}

template<class TElement>
void CubesGenerator<TElement>::elements(std::vector<Element*> &elements)
{
	_loader->elements(elements);
}

template<class TElement>
void CubesGenerator<TElement>::materials(std::vector<Material> &materials)
{
	_loader->materials(materials);
}

template<class TElement>
void CubesGenerator<TElement>::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	_loader->clusterBoundaries(nodes, neighbours);

	neighbours.push_back(_cubeIndex ? 0 : 1);

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemLength[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Z UP";
		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[1].problemLength[2] == _settings.cube[0].problemOrigin[2]) {

		std::cout << "Z BOTTOM";
		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemLength[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Y RIGHT";


		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[1].problemOrigin[0] == _settings.cube[0].problemLength[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Y LEFT";
		return;
	}

	if (
			_settings.cube[0].problemLength[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "X FRONT";
		return;
	}

	if (
			_settings.cube[1].problemOrigin[0] == _settings.cube[0].problemLength[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "X BACK";
		return;
	}

	ESINFO(GLOBAL_ERROR) << "Cubes do not fit";
}

template<class TElement>
void CubesGenerator<TElement>::settings(
	std::vector<Evaluator*> &evaluators,
	std::vector<Element*> &elements,
	std::vector<Element*> &faces,
	std::vector<Element*> &edges,
	std::vector<Element*> &nodes)
{
	_loader->settings(evaluators, elements, faces, edges, nodes);


	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemLength[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Z UP";
		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[1].problemLength[2] == _settings.cube[0].problemOrigin[2]) {

		std::cout << "Z BOTTOM";
		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemLength[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Y RIGHT";


		if (_cubeIndex) {
			eslocal cNodes[3];
			UniformUtils<TElement>::clusterNodesCount(_settings.cube[_cubeIndex], cNodes);
			std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));
			eslocal subdomain[3], element[3], subdomainOffset[3], elementOffset[3];


			for (subdomain[2] = 0; subdomain[2] < _settings.cube[_cubeIndex].subdomainsInCluster[2]; subdomain[2]++) {
				subdomainOffset[2] = subdomain[2] * (_settings.cube[_cubeIndex].elementsInSubdomain[2] * (1 + TElement::subnodes[2]));
				for (subdomain[1] = 0; subdomain[1] < 1; subdomain[1]++) {
					subdomainOffset[1] = subdomain[1] * (_settings.cube[_cubeIndex].elementsInSubdomain[1] * (1 + TElement::subnodes[1]));
					for (subdomain[0] = 0; subdomain[0] < _settings.cube[_cubeIndex].subdomainsInCluster[0]; subdomain[0]++) {
						subdomainOffset[0] = subdomain[0] * (_settings.cube[_cubeIndex].elementsInSubdomain[0] * (1 + TElement::subnodes[0]));
						// for each sub-domain

						for (element[2] = 0; element[2] < _settings.cube[_cubeIndex].elementsInSubdomain[2]; element[2]++) {
							elementOffset[2] = subdomainOffset[2] + element[2] * (1 + TElement::subnodes[2]);
							for (element[1] = 0; element[1] < 1; element[1]++) {
								elementOffset[1] = subdomainOffset[1] + element[1] * (1 + TElement::subnodes[1]);
								for (element[0] = 0; element[0] < _settings.cube[_cubeIndex].elementsInSubdomain[0]; element[0]++) {
									elementOffset[0] = subdomainOffset[0] + element[0] * (1 + TElement::subnodes[0]);
									// for each element

									eslocal i = 0;
									for (eslocal z = 0; z < 2 + TElement::subnodes[2]; z++) {
										for (eslocal y = 0; y < 2 + TElement::subnodes[1]; y++) {
											for (eslocal x = 0; x < 2 + TElement::subnodes[0]; x++) {
												// fill node indices

												indices[i++] =
														(elementOffset[2] + z) * cNodes[0] * cNodes[1] +
														(elementOffset[1] + y) * cNodes[0] +
														(elementOffset[0] + x);
											}
										}
									}
									_e.addFaces(faces, &indices[0], 3);
								}
							}
						}
					}
				}
			}
		} else {
			eslocal cNodes[3];
			UniformUtils<TElement>::clusterNodesCount(_settings.cube[_cubeIndex], cNodes);
			std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));
			eslocal subdomain[3], element[3], subdomainOffset[3], elementOffset[3];


			for (subdomain[2] = 0; subdomain[2] < _settings.cube[_cubeIndex].subdomainsInCluster[2]; subdomain[2]++) {
				subdomainOffset[2] = subdomain[2] * (_settings.cube[_cubeIndex].elementsInSubdomain[2] * (1 + TElement::subnodes[2]));
				for (subdomain[1] = _settings.cube[_cubeIndex].subdomainsInCluster[1] - 1; subdomain[1] < _settings.cube[_cubeIndex].subdomainsInCluster[1]; subdomain[1]++) {
					subdomainOffset[1] = subdomain[1] * (_settings.cube[_cubeIndex].elementsInSubdomain[1] * (1 + TElement::subnodes[1]));
					for (subdomain[0] = 0; subdomain[0] < _settings.cube[_cubeIndex].subdomainsInCluster[0]; subdomain[0]++) {
						subdomainOffset[0] = subdomain[0] * (_settings.cube[_cubeIndex].elementsInSubdomain[0] * (1 + TElement::subnodes[0]));
						// for each sub-domain

						for (element[2] = 0; element[2] < _settings.cube[_cubeIndex].elementsInSubdomain[2]; element[2]++) {
							elementOffset[2] = subdomainOffset[2] + element[2] * (1 + TElement::subnodes[2]);
							for (element[1] = _settings.cube[_cubeIndex].elementsInSubdomain[2] - 1; element[1] < _settings.cube[_cubeIndex].elementsInSubdomain[2]; element[1]++) {
								elementOffset[1] = subdomainOffset[1] + element[1] * (1 + TElement::subnodes[1]);
								for (element[0] = 0; element[0] < _settings.cube[_cubeIndex].elementsInSubdomain[0]; element[0]++) {
									elementOffset[0] = subdomainOffset[0] + element[0] * (1 + TElement::subnodes[0]);
									// for each element

									eslocal i = 0;
									for (eslocal z = 0; z < 2 + TElement::subnodes[2]; z++) {
										for (eslocal y = 0; y < 2 + TElement::subnodes[1]; y++) {
											for (eslocal x = 0; x < 2 + TElement::subnodes[0]; x++) {
												// fill node indices

												indices[i++] =
														(elementOffset[2] + z) * cNodes[0] * cNodes[1] +
														(elementOffset[1] + y) * cNodes[0] +
														(elementOffset[0] + x);
											}
										}
									}
									_e.addFaces(faces, &indices[0], 1);
								}
							}
						}
					}
				}
			}
		}

		for (size_t i = 0; i < faces.size(); i++) {
			faces[i]->clusters().push_back(0);
			faces[i]->clusters().push_back(1);
			faces[i]->settings(Property::NONMATCHING_ELEMENT);
		}
		return;
	}

	if (
			_settings.cube[0].problemOrigin[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[1].problemOrigin[0] == _settings.cube[0].problemLength[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "Y LEFT";
		return;
	}

	if (
			_settings.cube[0].problemLength[0] == _settings.cube[1].problemOrigin[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "X FRONT";
		return;
	}

	if (
			_settings.cube[1].problemOrigin[0] == _settings.cube[0].problemLength[0] &&
			_settings.cube[0].problemOrigin[1] == _settings.cube[1].problemOrigin[1] &&
			_settings.cube[0].problemOrigin[2] == _settings.cube[1].problemOrigin[2]) {

		std::cout << "X BACK";
		return;
	}
}

template<class TElement>
bool CubesGenerator<TElement>::partitiate(std::vector<eslocal> &parts)
{
	return _loader->partitiate(parts);
}

template<class TElement>
void CubesGenerator<TElement>::fixPoints(std::vector<std::vector<eslocal> > &fixPoints)
{
	_loader->fixPoints(fixPoints);
}

template<class TElement>
void CubesGenerator<TElement>::corners(std::vector<eslocal> &corners)
{
	_loader->corners(corners);
}

}
}

