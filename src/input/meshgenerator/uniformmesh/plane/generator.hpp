
#include "generator.h"

namespace espreso {
namespace input {

// Intel compilator error: do not use generoc constructor

//template<class TElement>
//PlaneGenerator<TElement>::PlaneGenerator(Mesh &mesh, const PlaneSettings &settings)
//	: UniformGenerator<TElement>(mesh, settings), _settings(settings)
//{
//	ESINFO(GLOBAL_ERROR) << "Plane generator does not support the selected element type.";
//}


template<class TElement>
void PlaneGenerator<TElement>::points(Coordinates &coordinates)
{
	eslocal cNodes[3];
	esglobal gNodes[3];

	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);

	coordinates.clear();
	coordinates.reserve(cNodes[0] * cNodes[1]);

	esglobal cs[2], ce[2];
	double step[2];
	for (eslocal i = 0; i < 2; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}
	for (eslocal i = 0; i < 2; i++) {
		step[i] = _settings.problemLength[i] / ((cNodes[i] - 1) * _settings.clusters[i]);
	}

	for (esglobal y = cs[1]; y <= ce[1]; y++) {
		for (esglobal x = cs[0]; x <= ce[0]; x++) {
			coordinates.add(
				Point(x * step[0], y * step[1], 0),
				(y - cs[1]) * cNodes[0] + (x - cs[0]),
				y * gNodes[0] + x
			);
		}
	}
}

template<class TElement>
void PlaneGenerator<TElement>::elementsMesh(std::vector<Element*> &elements)
{
	eslocal cNodes[3];

	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]));

	elements.clear();
	elements.reserve(UniformUtils<TElement>::clusterElementsCount(_settings));


	eslocal subdomain[2];
	eslocal element[2];

	eslocal subdomainOffset[2];
	eslocal elementOffset[2];

	eslocal params[6] = {0, 0, 0, 0, 0, 0};


	for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
		for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {
			// for each sub-domain

			for (eslocal i = 0; i < 2; i++) {
				subdomainOffset[i] = subdomain[i] * (_settings.elementsInSubdomain[i] * (1 + TElement::subnodes[i]));
			}
			for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
				for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {
					// for each element

					for (eslocal i = 0; i < 2; i++) {
						elementOffset[i] = subdomainOffset[i] + element[i] * (1 + TElement::subnodes[i]);
					}
					eslocal i = 0;
					for (eslocal y = 0; y < 2 + TElement::subnodes[1]; y++) {
						for (eslocal x = 0; x < 2 + TElement::subnodes[0]; x++) {
							// fill node indices

							indices[i++] =
									(elementOffset[1] + y) * cNodes[0] +
									(elementOffset[0] + x);
						}
					}
					this->_e.addElements(elements, &indices[0], params);
				}
			}
		}
	}
}

template<class TElement>
void PlaneGenerator<TElement>::elementsMaterials(std::vector<Element*> &elements)
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
void PlaneGenerator<TElement>::boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions)
{
	for (auto it = _settings.regions.begin(); it != _settings.regions.end(); ++it) {
		const Interval &interval = it->second;

		NodeCondition *region = new NodeCondition();
		std::vector<std::string> expressions;

		if (_settings.dirichlet.find(it->first) != _settings.dirichlet.end()) {
			region->setType(ConditionType::DIRICHLET);
			expressions = Parser::split(Parser::strip(_settings.dirichlet.find(it->first)->second), ",");
		}
		if (_settings.forces.find(it->first) != _settings.forces.end()) {
			region->setType(ConditionType::FORCES);
			expressions = Parser::split(Parser::strip(_settings.forces.find(it->first)->second), ",");
		}

		for (size_t i = 0; i < coordinates.clusterSize(); i++) {
			if (interval.isIn(coordinates[i])) {
				region->nodes.push_back(i);
			}
		}

		for (size_t i = 0; i < expressions.size(); i++) {
			std::string param = Parser::getParameter(expressions[i], ":");
			if (StringCompare::caseInsensitiveEq(param, "all")) {
				region->setValue(DOFType::DISPLACEMENT_X, Parser::getValue(expressions[i], ":"));
				region->setValue(DOFType::DISPLACEMENT_Y, Parser::getValue(expressions[i], ":"));
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
void PlaneGenerator<TElement>::initialConditions(const Coordinates &coordinates, std::vector<InitialCondition*> &conditions)
{
	auto getInterval = [] (const std::map<std::string, Interval> &regions, const std::string &name) {
		if (regions.find(name) == regions.end()) {
			ESINFO(GLOBAL_ERROR) << "Not specified region " << name;
		}
		return regions.find(name)->second;
	};

	if (_settings.heat_sources.size()) {
		if (_settings.heat_sources.size() == 1) {
			const Interval &interval = getInterval(_settings.regions, _settings.heat_sources.begin()->first);
			if (!interval.isIn(0, 0, 0) || !interval.isIn(_settings.problemLength[0], _settings.problemLength[1], 0)) {
				ESINFO(GLOBAL_ERROR) << "Heat source has to be set to whole problem.";
			}
			conditions[InitialCondition::index(InitialConditionType::HEAT_SOURCE)] = new ExpressionInitialization(_settings.heat_sources.begin()->second);
		} else {
			IntervalInitialization *condition = new IntervalInitialization();
			for (auto it = _settings.heat_sources.begin(); it != _settings.heat_sources.end(); ++it) {
				condition->add(getInterval(_settings.regions, it->first), it->second);
			}
		}
	}

	if (_settings.u.size()) {
		if (_settings.heat_sources.size() == 1) {
			const Interval &interval = getInterval(_settings.regions, _settings.u.begin()->first);
			if (!interval.isIn(0, 0, 0) || !interval.isIn(_settings.problemLength[0], _settings.problemLength[1], 0)) {
				ESINFO(GLOBAL_ERROR) << "Translation motion has to be set to whole problem.";
			}
			std::vector<std::string> expressions = Parser::split(Parser::strip(_settings.u.begin()->second), ",");
			if (expressions.size() != 2) {
				ESINFO(GLOBAL_ERROR) << "Invalid translation motion: INTERVAL = x: ? , y: ?";
			}
			conditions[InitialCondition::index(InitialConditionType::TRANSLATION_MOTION_X)] = new ExpressionInitialization(Parser::getValue(expressions[0], ":"));
			conditions[InitialCondition::index(InitialConditionType::TRANSLATION_MOTION_Y)] = new ExpressionInitialization(Parser::getValue(expressions[1], ":"));
		} else {
			IntervalInitialization *condition = new IntervalInitialization();
			for (auto it = _settings.u.begin(); it != _settings.u.end(); ++it) {
				std::vector<std::string> expressions = Parser::split(Parser::strip(it->second), ",");
				if (expressions.size() != 2) {
					ESINFO(GLOBAL_ERROR) << "Invalid translation motion: INTERVAL = x: ? , y: ?";
				}
				condition->add(getInterval(_settings.regions, it->first), Parser::getValue(expressions[0], ":"));
				condition->add(getInterval(_settings.regions, it->first), Parser::getValue(expressions[0], ":"));
			}
		}
	}
}

template<class TElement>
void PlaneGenerator<TElement>::fixPoints(std::vector<std::vector<eslocal> > &fixPoints)
{
	if (_settings.useMetis) {
		Loader::fixPoints(fixPoints);
		return;
	}

	fixPoints.reserve(_settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1]);
	eslocal shift_offset[2] = {TElement::subnodes[0] + 1, TElement::subnodes[1] + 1};

	eslocal nodes[2];
	eslocal cNodes[2];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	for (int i = 0; i < 2; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
		if (2 * (shift_offset[i] + 1) > nodes[i] + 1) { // not enough nodes
			shift_offset[i] = (nodes[i] + 1) / 2 - TElement::subnodes[i] - 1;
		}
		if (2 * shift_offset[i] == nodes[i]) { // offset to the same node
			shift_offset[i] -= TElement::subnodes[i] + 1;
		}
	}

	eslocal offset[2];
	eslocal shift[2];
	for (eslocal sy = 0; sy < _settings.subdomainsInCluster[1]; sy++) {
		for (eslocal sx = 0; sx < _settings.subdomainsInCluster[0]; sx++) {
			fixPoints.push_back(std::vector<eslocal>());
			fixPoints.back().reserve(4);
			for (int i = 0; i < 4; i++) {
				offset[0] = (i & 1) ? 1 : 0;
				offset[1] = (i & 2) ? 1 : 0;
				shift[0] = (i & 1) ? -shift_offset[0] : shift_offset[0];
				shift[1] = (i & 2) ? -shift_offset[1] : shift_offset[1];
				fixPoints.back().push_back(
						((sy + offset[1]) * nodes[1] + shift[1]) * cNodes[0] +
						((sx + offset[0]) * nodes[0] + shift[0]));
			}
		}
	}

	for (size_t p = 0; p < fixPoints.size(); p++) {
		for (size_t i = 0; i < fixPoints[p].size(); i++) {
			fixPoints[p][i] = this->mesh.coordinates().localIndex(fixPoints[p][i], p);
		}
		std::sort(fixPoints[p].begin(), fixPoints[p].end());

		// Remove the same points
		auto it = std::unique(fixPoints[p].begin(), fixPoints[p].end());
		fixPoints[p].resize(it - fixPoints[p].begin());
	}
}

template <class TElement>
void PlaneGenerator<TElement>::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	esglobal gNodes[3];
	CubeUtils<TElement>::globalNodesCount(_settings, gNodes);
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	boundaries.resize(cNodes[0] * cNodes[1]);

	bool border[2];
	eslocal cIndex = _cluster[0] + _cluster[1] * _settings.clusters[0];
	esglobal index = 0;

	esglobal cs[2], ce[2];
	for (eslocal i = 0; i < 2; i++) {
		cs[i] = (cNodes[i] - 1) * _cluster[i];
		ce[i] = (cNodes[i] - 1) * (_cluster[i] + 1);
	}

	// TODO: optimize this
	std::set<int> neighs;

	for (esglobal y = cs[1]; y <= ce[1]; y++) {
		border[1] = (y == 0 || y == gNodes[1] - 1) ? false : y % ( cNodes[1] - 1) == 0;
		for (esglobal x = cs[0]; x <= ce[0]; x++) {
			border[0] = (x == 0 || x == gNodes[0] - 1) ? false : x % ( cNodes[0] - 1) == 0;
			for (int i = 0; i < 4; i++) {
				eslocal tmp = cIndex;
				if (border[0] && (i & 1)) {
					tmp += (x == cs[0]) ? -1 : 1;
				}
				if (border[1] && (i & 2)) {
					tmp += ((y == cs[1]) ? -1 : 1) * _settings.clusters[0];
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

	neighs.erase(config::env::MPIrank);
	neighbours.insert(neighbours.end(), neighs.begin(), neighs.end());
}

template <class TElement>
void PlaneGenerator<TElement>::corners(Boundaries &boundaries)
{
	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::TOTAL_FETI) {
		// corners are not used in the case of TOTAL FETI
		return;
	}

	if (_settings.useMetis) {
		Loader::corners(boundaries);
		return;
	}

	if (_settings.corners) {
		ESINFO(DETAILS) << "Set corners to vertices";
	}
	ESINFO(DETAILS) << "Number of corners on each edge is " << (_settings.edges ? _settings.cornerCount : 0) << ".";

	eslocal nodes[2];
	eslocal cNodes[2];
	for (int i = 0; i < 2; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal step[2];
	for (int i = 0; i < 2; i++) {
		step[i] = _settings.elementsInSubdomain[i] / (_settings.cornerCount + 1);
		step[i] *= TElement::subnodes[i] + 1;
	}
	std::vector<std::vector<size_t> > offsets(2);
	std::vector<size_t> mul(2);

	for (int i = 0; i < 2; i++) {
		for (size_t j = 0; j < _settings.subdomainsInCluster[i]; j++) {
			for (size_t k = 0; k <= _settings.cornerCount / 2; k++) {
				offsets[i].push_back(j * nodes[i] + k * step[i]);
				offsets[i].push_back(j * nodes[i] + nodes[i] - k * step[i]);
			}
			if (_settings.cornerCount % 2 == 1) {
				eslocal mid = (_settings.elementsInSubdomain[i] / 2) * (TElement::subnodes[i] + 1);
				offsets[i].push_back(j * nodes[i] + mid);
			}
		}
	}
	mul[0] = 1;
	mul[1] = cNodes[0];

	eslocal index;
	for (size_t d = 0; d < 2; d++) {
		for (eslocal i = 1; i < _settings.subdomainsInCluster[d]; i++) {
			for (size_t j = 0; j < offsets[(d + 1) % 2].size(); j++) {
				if (!_settings.corners && offsets[(d + 1) % 2][j] % nodes[(d + 1) % 2] == 0)
				{
					continue;
				}
				if (!_settings.edges && offsets[(d + 1) % 2][j] % nodes[(d + 1) % 2] != 0)
				{
					continue;
				}

				index = i * nodes[d] * mul[d];
				index += offsets[(d + 1) % 2][j] * mul[(d + 1) % 2];
				boundaries.setCorner(index);
			}
		}
	}

	if (config::mesh::AVERAGE_EDGES || config::mesh::AVERAGE_FACES) {
		// TODO: check correctness
		this->mesh.computeCorners(0, true, false, false, config::mesh::AVERAGE_EDGES, config::mesh::AVERAGE_FACES);
	}

}

}
}


