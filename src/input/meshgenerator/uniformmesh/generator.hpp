
#include "generator.h"

namespace espreso {
namespace input {

template<class TElement>
void UniformGenerator<TElement>::elementsMesh(std::vector<Element*> &elements)
{
	eslocal cNodes[3];

	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	std::vector<eslocal> indices((2 + TElement::subnodes[0]) * (2 + TElement::subnodes[1]) * (2 + TElement::subnodes[2]));

	elements.clear();
	elements.reserve(UniformUtils<TElement>::clusterElementsCount(_settings));


	eslocal subdomain[3];
	eslocal element[3];

	eslocal subdomainOffset[3];
	eslocal elementOffset[3];

	eslocal params[6] = {0, 0, 0, 0, 0, 0};


	for (subdomain[2] = 0; subdomain[2] < _settings.subdomainsInCluster[2]; subdomain[2]++) {
		for (subdomain[1] = 0; subdomain[1] < _settings.subdomainsInCluster[1]; subdomain[1]++) {
			for (subdomain[0] = 0; subdomain[0] < _settings.subdomainsInCluster[0]; subdomain[0]++) {
				// for each sub-domain

				for (eslocal i = 0; i < 3; i++) {
					subdomainOffset[i] = subdomain[i] * (_settings.elementsInSubdomain[i] * (1 + TElement::subnodes[i]));
				}
				for (element[2] = 0; element[2] < _settings.elementsInSubdomain[2]; element[2]++) {
					for (element[1] = 0; element[1] < _settings.elementsInSubdomain[1]; element[1]++) {
						for (element[0] = 0; element[0] < _settings.elementsInSubdomain[0]; element[0]++) {
							// for each element

							for (eslocal i = 0; i < 3; i++) {
								elementOffset[i] = subdomainOffset[i] + element[i] * (1 + TElement::subnodes[i]);
							}
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
							_e.addElements(elements, &indices[0], params);
						}
					}
				}
			}
		}
	}
}

template<class TElement>
void UniformGenerator<TElement>::partitiate(std::vector<eslocal> &parts)
{
	config::mesh::SUBDOMAINS = _settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1] * _settings.subdomainsInCluster[2];
	if (_settings.useMetis) {
		Loader::partitiate(parts);
		return;
	}

	parts.clear();
	parts.reserve(config::mesh::SUBDOMAINS + 1);

	parts.push_back(0);

	for (size_t p = 0; p < config::mesh::SUBDOMAINS; p++) {
		parts.push_back(parts.back() + mesh.getElements().size() / config::mesh::SUBDOMAINS);
	}

	Loader::remapElementsToSubdomains();
	Loader::computeBoundaries();
}

template<class TElement>
void UniformGenerator<TElement>::fixPoints(std::vector<std::vector<eslocal> > &fixPoints)
{
	if (_settings.useMetis) {
		Loader::fixPoints(fixPoints);
		return;
	}

	fixPoints.reserve(_settings.subdomainsInCluster[0] * _settings.subdomainsInCluster[1] * _settings.subdomainsInCluster[2]);
	eslocal SHIFT = 1;
	eslocal shift_offset[3] = {SHIFT, SHIFT, SHIFT};

	eslocal nodes[3];
	eslocal cNodes[3];
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
		if (2 * (shift_offset[i] + 1) > nodes[i] + 1) { // not enough nodes
			shift_offset[i] = (nodes[i] + 1) / 2 - 1;
		}
		if (2 * shift_offset[i] == nodes[i]) { // offset to the same node
			shift_offset[i]--;
		}
	}

	eslocal offset[3];
	eslocal shift[3];
	for (eslocal sz = 0; sz < _settings.subdomainsInCluster[2]; sz++) {
		for (eslocal sy = 0; sy < _settings.subdomainsInCluster[1]; sy++) {
			for (eslocal sx = 0; sx < _settings.subdomainsInCluster[0]; sx++) {
				fixPoints.push_back(std::vector<eslocal>());
				fixPoints.back().reserve(8);
				for (int i = 0; i < 8; i++) {
					offset[0] = (i & 1) ? 1 : 0;
					offset[1] = (i & 2) ? 1 : 0;
					offset[2] = (i & 4) ? 1 : 0;
					shift[0] = (i & 1) ? -shift_offset[0] : shift_offset[0];
					shift[1] = (i & 2) ? -shift_offset[1] : shift_offset[1];
					shift[2] = (i & 4) ? -shift_offset[2] : shift_offset[2];
					fixPoints.back().push_back(
							((sz + offset[2]) * nodes[2] + shift[2]) * cNodes[0] * cNodes[1] +
							((sy + offset[1]) * nodes[1] + shift[1]) * cNodes[0] +
							((sx + offset[0]) * nodes[0] + shift[0]));
				}
			}
		}
	}

	for (size_t p = 0; p < fixPoints.size(); p++) {
		for (size_t i = 0; i < fixPoints[p].size(); i++) {
			fixPoints[p][i] = mesh.coordinates().localIndex(fixPoints[p][i], p);
		}
		std::sort(fixPoints[p].begin(), fixPoints[p].end());

		// Remove the same points
		auto it = std::unique(fixPoints[p].begin(), fixPoints[p].end());
		fixPoints[p].resize(it - fixPoints[p].begin());
	}
}

template <class TElement>
void UniformGenerator<TElement>::corners(Boundaries &boundaries)
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
	ESINFO(DETAILS)
		<< "Number of corners on each edge is " << (_settings.edges ? _settings.cornerCount : 0) << "."
		<< "Number of corners on each face is " << (_settings.edges ? _settings.cornerCount * _settings.cornerCount : 0);

	eslocal nodes[3];
	eslocal cNodes[3];
	for (int i = 0; i < 3; i++) {
		nodes[i] = (TElement::subnodes[i] + 1) * _settings.elementsInSubdomain[i];
	}
	UniformUtils<TElement>::clusterNodesCount(_settings, cNodes);

	eslocal step[3];
	for (int i = 0; i < 3; i++) {
		step[i] = _settings.elementsInSubdomain[i] / (_settings.cornerCount + 1);
		step[i] *= TElement::subnodes[i] + 1;
	}
	std::vector<std::vector<size_t> > offsets(3);
	std::vector<size_t> mul(3);

	for (int i = 0; i < 3; i++) {
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
	mul[2] = cNodes[0] * cNodes[1];

	eslocal index;
	for (size_t d = 0; d < 3; d++) {
		for (eslocal i = 1; i < _settings.subdomainsInCluster[d]; i++) {
			for (size_t j = 0; j < offsets[(d + 1) % 3].size(); j++) {
				for (size_t k = 0; k < offsets[(d + 2) % 3].size(); k++) {
					if (!_settings.corners
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!_settings.edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] == 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					if (!_settings.edges
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] == 0)
					{
						continue;
					}
					if (!_settings.faces
						&& offsets[(d + 1) % 3][j] % nodes[(d + 1) % 3] != 0
						&& offsets[(d + 2) % 3][k] % nodes[(d + 2) % 3] != 0)
					{
						continue;
					}
					index = i * nodes[d] * mul[d];
					index += offsets[(d + 1) % 3][j] * mul[(d + 1) % 3];
					index += offsets[(d + 2) % 3][k] * mul[(d + 2) % 3];
					boundaries.setCorner(index);
				}
			}
		}
	}

	if (config::mesh::AVERAGE_EDGES || config::mesh::AVERAGE_FACES) {
		// TODO: check correctness
		mesh.computeCorners(0, true, false, false, config::mesh::AVERAGE_EDGES, config::mesh::AVERAGE_FACES);
	}

}

}
}

