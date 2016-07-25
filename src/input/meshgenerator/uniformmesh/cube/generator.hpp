
#include "generator.h"

namespace espreso {
namespace input {


//	###################################################
//	#                                                 #
//	#             A z-coord.                          #
//	#             |                                   #
//	#             |            E3                     #
//	#             |_ _ _ _ _ _ _                      #
//	#            /     E5      /|                     #
//	#           /_ _ _ _ _ _  / |                     #
//	#          |      |      |  |                     #
//	#        E4|      |      |E2|                     #
//	#          |_ _ _ |_ _ _ |  |       y-coord.      #
//	#          |    E1|      |  |------->             #
//	#          |      |      | /                      #
//	#          |_ _ _ |_ _ _ |/                       #
//	#         /                                       #
//	#        /       E0                               #
//	#       /                                         #
//	#      v  x-coord.                                #
//	#                                                 #
//	###################################################

// Intel compilator error: do not use generic contructor

//template<class TElement>
//CubeGenerator<TElement>::CubeGenerator(Mesh &mesh, const CubeSettings &settings)
//: UniformGenerator<TElement>(mesh, settings), _settings(settings)
//{
//	ESINFO(GLOBAL_ERROR) << "Cube generator does not support the selected element type.";
//}

template<class TElement>
void CubeGenerator<TElement>::elementsMaterials(std::vector<Element*> &elements)
{
	esglobal cubeElements[3], partSize[3], cOffset[3], offset[3];
	eslocal subdomain[3], element[3], material, counter;

	for (size_t i = 0; i < 3; i++) {
		cubeElements[i] = _settings.clusters[i] * _settings.subdomainsInCluster[i] * _settings.elementsInSubdomain[i];
		cOffset[i] = _cluster[i] * _settings.subdomainsInCluster[i] * _settings.elementsInSubdomain[i];
		partSize[i] = std::ceil(cubeElements[i] / (double)_settings.materialsLayers[i]);
	}

	counter = 0;
	for (subdomain[2] = 0; subdomain[2] < _settings.subdomainsInCluster[2]; subdomain[2]++) {
			for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
				for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {

					for (element[2] = 0; element[2] < _settings.elementsInSubdomain[2]; element[2]++) {
						for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
							for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {

								material = 0;
								for (eslocal i = 0; i < 3; i++) {
									offset[i] = cOffset[i] + subdomain[i] * _settings.elementsInSubdomain[i] + element[i];
									if (offset[i] / partSize[i] % 2 == 1) {
										material = (material + 1) % 2;
									}
								}
								for (size_t e = 0; e < TElement::subelements; e++) {
									elements[counter++]->setParam(Element::MATERIAL, material);
								}
							}
						}
					}

				}
			}
	}

}


template<class TElement>
void CubeGenerator<TElement>::points(Coordinates &coordinates, size_t &DOFs)
{
	DOFs = this->_DOFs;

	eslocal cNodes[3];
	esglobal gNodes[3];

	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1] * cNodes[2]);

	esglobal cs[3], ce[3];
	double step[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}
	for (eslocal i = 0; i < 3; i++) {
		step[i] = _settings.problemLength[i] / ((cNodes[i] - 1) * _settings.clusters[i]);
	}

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				coordinates.add(
					Point(x * step[0], y * step[1], z * step[2]),
					(z - cs[2]) * cNodes[0] * cNodes[1] + (y - cs[1]) * cNodes[0] + (x - cs[0]),
					z * gNodes[0] * gNodes[1] + y * gNodes[0] + x
				);
			}
		}
	}
}

template<class TElement>
void CubeGenerator<TElement>::boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions)
{
	for (auto it = _settings.regions.begin(); it != _settings.regions.end(); ++it) {

		if (_settings.dirichlet.find(it->first) == _settings.dirichlet.end()) {
			ESINFO(ALWAYS) << TextColor::YELLOW << "Warning: Region without dirichlet.";
			continue;
		}

		const Interval &interval = it->second;
		NodeCondition *region = new NodeCondition();
		region->setType(ConditionType::DIRICHLET);

		for (size_t i = 0; i < coordinates.clusterSize(); i++) {
			if (interval.isIn(coordinates[i])) {
				region->nodes.push_back(i);
			}
		}

		std::vector<std::string> expressions = Parser::split(Parser::strip(_settings.dirichlet.find(it->first)->second), ",");

		for (size_t i = 0; i < expressions.size(); i++) {
			std::string param = Parser::getParameter(expressions[i], ":");
			if (StringCompare::caseInsensitiveEq(param, "all")) {
				region->setValue(DOFType::DISPLACEMENT_X, Parser::getValue(expressions[i], ":"));
				region->setValue(DOFType::DISPLACEMENT_Y, Parser::getValue(expressions[i], ":"));
				region->setValue(DOFType::DISPLACEMENT_Z, Parser::getValue(expressions[i], ":"));
				region->setValue(DOFType::TEMPERATURE, Parser::getValue(expressions[i], ":"));
				region->setValue(DOFType::PRESSURE, Parser::getValue(expressions[i], ":"));
				continue;
			}
			if (StringCompare::caseInsensitiveEq(param, "ux")) {
				region->setValue(DOFType::DISPLACEMENT_X, Parser::getValue(expressions[i], ":"));
				continue;
			}
			if (StringCompare::caseInsensitiveEq(param, "uy")) {
				region->setValue(DOFType::DISPLACEMENT_Y, Parser::getValue(expressions[i], ":"));
				continue;
			}
			if (StringCompare::caseInsensitiveEq(param, "uz")) {
				region->setValue(DOFType::DISPLACEMENT_Z, Parser::getValue(expressions[i], ":"));
				continue;
			}
			if (StringCompare::caseInsensitiveEq(param, "t")) {
				region->setValue(DOFType::TEMPERATURE, Parser::getValue(expressions[i], ":"));
				continue;
			}
			if (StringCompare::caseInsensitiveEq(param, "p")) {
				region->setValue(DOFType::PRESSURE, Parser::getValue(expressions[i], ":"));
				continue;
			}
			ESINFO(GLOBAL_ERROR) << "Unknown DOF type " << param;
		}

		conditions.push_back(region);
	}
}

template <class TElement>
void CubeGenerator<TElement>::initialConditions(const Coordinates &coordinates, std::vector<InitialCondition*> &conditions)
{

}


template <class TElement>
void CubeGenerator<TElement>::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	esglobal gNodes[3];
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(cNodes[0] * cNodes[1] * cNodes[2]);

	bool border[3];
	eslocal cIndex = _cluster[0] + _cluster[1] * _settings.clusters[0] + _cluster[2] * _settings.clusters[0] * _settings.clusters[1];
	esglobal index = 0;

	esglobal cs[3], ce[3];
	for (eslocal i = 0; i < 3; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}

	// TODO: optimize this
	std::set<int> neighs;

	for (esglobal z = cs[2]; z <= ce[2]; z++) {
		border[2] = (z == 0 || z == gNodes[2] - 1) ? false : z % ( cNodes[2] - 1) == 0;
		for (esglobal y = cs[1]; y <= ce[1]; y++) {
			border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
			for (esglobal x = cs[0]; x <= ce[0]; x++) {
				border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
				for (int i = 0; i < 8; i++) {
					eslocal tmp = cIndex;
					if (border[0] && (i & 1)) {
						tmp += (x == cs[0]) ? -1 : 1;
					}
					if (border[1] && (i & 2)) {
						tmp += ((y == cs[1]) ? -1 : 1) * _settings.clusters[0];
					}
					if (border[2] && (i & 4)) {
						tmp += ((z == cs[2]) ? -1 : 1) * _settings.clusters[0] * _settings.clusters[1];
					}
					boundaries[index].push_back(tmp);
					neighs.insert(tmp);
				}
				std::sort(boundaries[index].begin(), boundaries[index].end());
				auto end = std::unique(boundaries[index].begin(), boundaries[index].end());
				boundaries[index].resize(end - boundaries[index].begin());
				index++;
			}
		}
	}

	neighs.erase(config::env::MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

}
}

