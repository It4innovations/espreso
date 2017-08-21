
#include "block.h"

#include "../../../mesh/structures/region.h"
#include "../../../basis/utilities/utils.h"
#include "../generator.h"

#include <numeric>

namespace espreso {
namespace input {

template <class TElement>
void Block<TElement>::points(std::vector<Point> &points)
{
	Triple<size_t> nodes = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) - 1) + 1;
	Triple<double> step = (block.end - block.start) / Triple<double>(nodes.x, nodes.y, nodes.z).steps();

	points.reserve(nodes.mul());
	std::vector<double> p;
	for (size_t z = 0; z < nodes.z; z++) {
		for (size_t y = 0; y < nodes.y; y++) {
			for (size_t x = 0; x < nodes.x; x++) {
				p = {
						Generator::precision * (block.start.x + x * step.x),
						Generator::precision * (block.start.y + y * step.y),
						Generator::precision * (block.start.z + z * step.z)};
				points.push_back(Point(block.projection.x(p), block.projection.y(p), block.projection.z(p), block.rotation.x(p), block.rotation.y(p), block.rotation.z(p)));
			}
		}
	}
}

template <class TElement>
void Block<TElement>::forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation)
{
	forEachElement(start, end, operation, [] (Triple<size_t> &offset) {});
}

template <class TElement>
void Block<TElement>::forEachElement(const Triple<size_t> &start, const Triple<size_t> &end, std::function<void(std::vector<eslocal> &indices)> operation, std::function<void(Triple<size_t> &offset)> restriction)
{
	Triple<size_t> enodes = Triple<size_t>(TElement::subnodes) - 1;
	Triple<size_t> size = (block.domains * block.elements * enodes + 1).toSize();

	Triple<size_t> _start = start;
	Triple<size_t> _end = end;

	auto correct = [] (size_t &s, size_t &e) {
		if (s == e) { if (s == 0) { e++; } else { s--; } }
	};

	correct(_start.x, _end.x);
	correct(_start.y, _end.y);
	correct(_start.z, _end.z);

	std::vector<eslocal> indices((enodes + 1).mul());

	Triple<size_t> e;
	for (e.z = _start.z; e.z < _end.z; e.z++) {
		for (e.y = _start.y; e.y < _end.y; e.y++) {
			for (e.x = _start.x; e.x < _end.x; e.x++) {

				size_t index = 0;
				Triple<size_t> eOffset, offset = e * enodes;
				for (eOffset.z = 0; eOffset.z <= enodes.z; eOffset.z++) {
					for (eOffset.y = 0; eOffset.y <= enodes.y; eOffset.y++) {
						for (eOffset.x = 0; eOffset.x <= enodes.x; eOffset.x++) {
							Triple<size_t> coffset = eOffset + offset;
							restriction(coffset);
							indices[index++] = (coffset * size).sum();
						}
					}
				}
				operation(indices);
			}
		}
	}
}


template <class TElement>
void Block<TElement>::elements(std::vector<Element*> &elements, size_t body)
{
	elements.reserve(TElement::subelements * (block.domains * block.elements).mul());
	std::vector<eslocal> params(6);
	params[Element::Params::BODY] = body;

	Triple<size_t> offset;
	for (offset.z = 0; offset.z < block.domains.z; offset.z++) {
		for (offset.y = 0; offset.y < block.domains.y; offset.y++) {
			for (offset.x = 0; offset.x < block.domains.x; offset.x++) {
				forEachElement(offset * block.elements, (offset + 1) * block.elements, [&] (std::vector<eslocal> &indices){
					TElement::addElements(elements, indices.data(), params.data());
				});
			}
		}
	}
}

template<class TElement>
void Block<TElement>::uniformPartition(std::vector<eslocal> &partsPtrs, size_t subdomains)
{
	partsPtrs.clear();
	partsPtrs.reserve(subdomains + 1);

	size_t elements = TElement::subelements * (block.domains * block.elements).mul();
	partsPtrs.push_back(0);
	for (size_t p = 0; p < subdomains; p++) {
		partsPtrs.push_back(partsPtrs.back() + elements / subdomains);
	}
}

template<class TElement>
void Block<TElement>::uniformFixPoints(const std::vector<Element*> &nodes, std::vector<std::vector<Element*> > &fixPoints)
{
	Triple<int>shift = Triple<int>(TElement::subnodes) - 1;
	Triple<eslocal> dnodes = (Triple<int>(TElement::subnodes) - 1) * block.elements;
	Triple<eslocal> size = (dnodes * block.domains + 1).toSize();

	auto shiftCorrection = [] (int &shift, eslocal &nodes, eslocal subnodes) {
		if (2 * (shift + 1) > nodes + 1) { // not enough nodes
			shift = (nodes + 1) / 2 - subnodes - 1;
		}
		if (2 * shift == nodes) { // offset to the same node
			shift -= subnodes + 1;
		}
		if (shift < 0) {
			shift = 0;
		}
	};

	shiftCorrection(shift.x, dnodes.x, TElement::subnodes[0]);
	shiftCorrection(shift.y, dnodes.y, TElement::subnodes[1]);
	shiftCorrection(shift.z, dnodes.z, TElement::subnodes[2]);

	size_t number = 1;
	number *= TElement::subnodes[0] > 1 ? 2 : 1;
	number *= TElement::subnodes[1] > 1 ? 2 : 1;
	number *= TElement::subnodes[2] > 1 ? 2 : 1;

	fixPoints.resize(block.domains.mul(), std::vector<Element*>(number, NULL));

	Triple<size_t> domain;
	for (domain.z = 0; domain.z < block.domains.z; domain.z++) {
		for (domain.y = 0; domain.y < block.domains.y; domain.y++) {
			for (domain.x = 0; domain.x < block.domains.x; domain.x++) {
				for (size_t i = 0; i < number; i++) {
					Triple<eslocal> corner((i & 1) ? domain.x + 1 : domain.x, (i & 2) ? domain.y + 1 : domain.y, (i & 4) ? domain.z + 1 : domain.z);
					Triple<eslocal> offset((i & 1) ? -shift.x : shift.x, (i & 2) ? -shift.y : shift.y, (i & 4) ? -shift.z : shift.z);
					fixPoints[(domain * block.domains.toSize()).sum()][i] = nodes[((corner * dnodes + offset) * size).sum()];
				}
			}
		}
	}
}

template<class TElement>
void Block<TElement>::uniformCorners(const std::vector<Element*> &nodes, std::vector<Element*> &corners, size_t number, bool point, bool edge, bool face)
{
	Triple<size_t> dnodes = (Triple<size_t>(TElement::subnodes) - 1) * block.elements;
	Triple<size_t> cnodes = block.domains * dnodes;
	Triple<size_t> size = (cnodes + 1).toSize();
	Triple<size_t> step = dnodes / (number + 1);

	Triple<std::vector<size_t> > offsets;

	auto addOffset = [&] (std::vector<size_t> &offset, size_t domains, size_t nodes, size_t step) {
		for (size_t j = 0; j < domains; j++) {
			for (size_t k = 0; k <= number / 2; k++) {
				offset.push_back(j * nodes + k * step);
				offset.push_back(j * nodes + nodes - k * step);
			}
			if (number % 2 == 1) {
				offset.push_back(j * nodes + nodes / 2);
			}
		}
		std::sort(offset.begin(), offset.end());
		Esutils::removeDuplicity(offset);
	};

	addOffset(offsets.x, block.domains.x, dnodes.x, step.x);
	addOffset(offsets.y, block.domains.y, dnodes.y, step.y);
	addOffset(offsets.z, block.domains.z, dnodes.z, step.z);

	auto zero = [] (const Triple<size_t> &offset, size_t n) {
		size_t c = 0;
		if (!offset.x) { c++; }
		if (!offset.y) { c++; }
		if (!offset.z) { c++; }
		return c == n;
	};

	dnodes = (dnodes + 1).steps();
	cnodes = (cnodes + 1).steps();

	Triple<size_t> offset;
	for (size_t z = 0; z < offsets.z.size(); z++) {
		offset.z = offsets.z[z];
		for (size_t y = 0; y < offsets.y.size(); y++) {
			offset.y = offsets.y[y];
			for (size_t x = 0; x < offsets.x.size(); x++) {
				offset.x = offsets.x[x];

				if (zero(offset % cnodes, 3)) {
					continue;
				}

				if (zero(offset % cnodes, 2) && !zero(offset % dnodes, 3)) {
					continue;
				}

				if (zero(offset % cnodes, 1) && zero(offset % dnodes, 1)) {
					continue;
				}

				if (!point && zero(offset % dnodes, 3)) {
					continue;
				}

				if (!edge && zero(offset % dnodes, 2)) {
					continue;
				}

				if (!face && zero(offset % dnodes, 1)) {
					continue;
				}

				corners.push_back(nodes[(offset * size).sum()]);
			}
		}
	}
}

template <class TElement>
void Block<TElement>::boundaries(std::vector<Element*> &nodes, const std::vector<int> &neighbours)
{
	std::vector<int> sorted(neighbours.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(), [&] (int i, int j) { return neighbours[i] < neighbours[j]; });

	Triple<size_t> count = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) - 1);
	Triple<size_t> size = (count + 1).toSize();

	for (int i = 0; i < 27; i++) {
		if (neighbours[sorted[i]] < 0) {
			continue;
		}
		size_t  sx = sorted[i] % 3 == 2 ? count.x : 0,
				ex = sorted[i] % 3 != 0 ? count.x : 0,
				sy = sorted[i] % 9 >= 6 ? count.y : 0,
				ey = sorted[i] % 9 >= 3 ? count.y : 0,
				sz = sorted[i] / 9 == 2 ? count.z : 0,
				ez = sorted[i] / 9 >= 1 ? count.z : 0;


		Triple<size_t> offset;
		for (offset.x = sx; offset.x <= ex; offset.x++) {
			for (offset.y = sy; offset.y <= ey; offset.y++) {
				for (offset.z = sz; offset.z <= ez; offset.z++) {
					nodes[(offset * size).sum()]->clusters().push_back(neighbours[sorted[i]]);
				}
			}
		}
	}
}

template <class TElement>
void Block<TElement>::pickElements(const Triple<size_t> &start, const Triple<size_t> &end, const std::vector<Element*> &elements, Region *region)
{
	auto isIn = [&] (const Triple<size_t> &index) {
		return
				start.x <= index.x && index.x < end.x &&
				start.y <= index.y && index.y < end.y &&
				start.z <= index.z && index.z < end.z;
	};

	Triple<size_t> domain, element; size_t eindex = 0;
	for (domain.z = 0; domain.z < block.domains.z; domain.z++) {
		for (domain.y = 0; domain.y < block.domains.y; domain.y++) {
			for (domain.x = 0; domain.x < block.domains.x; domain.x++) {

				for (element.z = 0; element.z < block.elements.z; element.z++) {
					for (element.y = 0; element.y < block.elements.y; element.y++) {
						for (element.x = 0; element.x < block.elements.x; element.x++, eindex += TElement::subelements) {

							if (isIn(domain * block.elements + element)) {
								for (size_t e = 0; e < TElement::subelements; e++) {
									region->elements().push_back(elements[eindex + e]);
								}
							}

						}
					}
				}

			}
		}
	}
}

template <class TElement>
void Block<TElement>::region(const std::vector<Element*> &elements, Region *region, const BlockBorder &border, size_t dimension)
{
	BlockBorder intersection = border.intersect(block);
	if (intersection.dimension() > 3) {
		return;
	}

	Triple<size_t> beginElement, endElement = block.domains * block.elements;
	Triple<size_t> minNode, maxNode = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) - 1);
	Triple<size_t> lastNode = maxNode;

	// project intersection to interval <0, 1>
	Triple<double> istart = (Triple<double>)(intersection.start - block.start) / (block.end - block.start);
	Triple<double> iend   = (Triple<double>)(intersection.end - block.start) / (block.end - block.start);

	// get start, end node
	Triple<double> nStart = istart * (maxNode - minNode);
	Triple<double> nEnd   = iend   * (maxNode - minNode);
	// shrink interval to the closest nodes
	minNode = nStart.ceil();
	maxNode = nEnd.floor();

	// If interval is open, checks whether min, max node should be modev
	Triple<double> cStart = nStart.ceil();
	Triple<double> fEnd   = nEnd.floor();
	nStart = (nStart / Generator::precision).ceil();
	nEnd   = (nEnd   / Generator::precision).floor();

	if (intersection.excludeStart.x && cStart.x / Generator::precision == nStart.x) { minNode.x++; }
	if (intersection.excludeStart.y && cStart.y / Generator::precision == nStart.y) { minNode.y++; }
	if (intersection.excludeStart.z && cStart.z / Generator::precision == nStart.z) { minNode.z++; }
	if (intersection.excludeEnd.x && fEnd.x / Generator::precision == nEnd.x) { maxNode.x--; }
	if (intersection.excludeEnd.y && fEnd.y / Generator::precision == nEnd.y) { maxNode.y--; }
	if (intersection.excludeEnd.z && fEnd.z / Generator::precision == nEnd.z) { maxNode.z--; }

	// No matter what above happen, neighbor nodes have to match! Check it again
	if (intersection.start.x == block.start.x) { minNode.x = intersection.excludeStart.x ? 1 : 0; }
	if (intersection.start.y == block.start.y) { minNode.y = intersection.excludeStart.y ? 1 : 0; }
	if (intersection.start.z == block.start.z) { minNode.z = intersection.excludeStart.z ? 1 : 0; }
	if (intersection.end.x == block.end.x) { maxNode.x = intersection.excludeEnd.x ? lastNode.x - 1 : lastNode.x; }
	if (intersection.end.y == block.end.y) { maxNode.y = intersection.excludeEnd.y ? lastNode.y - 1 : lastNode.y; }
	if (intersection.end.z == block.end.z) { maxNode.z = intersection.excludeEnd.z ? lastNode.z - 1 : lastNode.z; }

	Triple<size_t> enodes = Triple<size_t>(TElement::subnodes);
	if (enodes.z == 1) { enodes.z++; }

	// Add only elements that have included more that one node
	beginElement = minNode / (enodes - 1);
	endElement = ((Triple<double>)(maxNode) / (enodes - 1)).ceil();
	if (TElement::subnodes[2] == 1) {
		endElement.z++;
	}

	// There is still change to miss all nodes by interval.
	if (minNode.x > maxNode.x) { return; }
	if (minNode.y > maxNode.y) { return; }
	if (minNode.z > maxNode.z) { return; }

	switch (intersection.dimension()) {

	case 0:   // Point
	case 1: { // Line

		CubeEdge edge = intersection.getEdge(block);
		switch (dimension) {
		case 0:
			forEachElement(beginElement, endElement,
			[&] (std::vector<eslocal> &indices) {
				TElement::pickNodes(elements, region->elements(), indices.data(), edge);
			},
			[&] (Triple<size_t> &offset) {
				offset.x = offset.x < minNode.x ? minNode.x : maxNode.x < offset.x ? maxNode.x : offset.x;
				offset.y = offset.y < minNode.y ? minNode.y : maxNode.y < offset.y ? maxNode.y : offset.y;
				offset.z = offset.z < minNode.z ? minNode.z : maxNode.z < offset.z ? maxNode.z : offset.z;
			});
			std::sort(region->elements().begin(), region->elements().end());
			Esutils::removeDuplicity(region->elements());
			break;
		case 1:
			forEachElement(beginElement, endElement,
			[&] (std::vector<eslocal> &indices) {
				TElement::addEdges(region->elements(), indices.data(), edge);
			});
			break;
		default:
			if (border.dimension() == intersection.dimension()) {
				ESINFO(ERROR) << "Cannot select element of dimension " << dimension << " on 1D line.";
			}
		}
	} break;

	case 2: { // Face
		CubeFace face = intersection.getFace(block);
		switch (dimension) {
		case 0:
			forEachElement(beginElement, endElement,
			[&] (std::vector<eslocal> &indices) {
				TElement::pickNodes(elements, region->elements(), indices.data(), face);
			},
			[&] (Triple<size_t> &offset) {
				offset.x = offset.x < minNode.x ? minNode.x : maxNode.x < offset.x ? maxNode.x : offset.x;
				offset.y = offset.y < minNode.y ? minNode.y : maxNode.y < offset.y ? maxNode.y : offset.y;
				offset.z = offset.z < minNode.z ? minNode.z : maxNode.z < offset.z ? maxNode.z : offset.z;
			});
			std::sort(region->elements().begin(), region->elements().end());
			Esutils::removeDuplicity(region->elements());
			break;
		case 1:
			ESINFO(ERROR) << "Implement selection of edges on face";
			break;
		case 2:
			forEachElement(beginElement, endElement,
			[&] (std::vector<eslocal> &indices) {
				TElement::addFaces(region->elements(), indices.data(), face);
			});
			break;
		case 3:
			pickElements(beginElement, endElement, elements, region);
			break;
		default:
			ESINFO(ERROR) << "Cannot select element of dimension " << dimension << " on 2D plane.";
		}
	} break;

	case 3: { // Volume

		switch (dimension) {
		case 0:
		case 1:
		case 2:
			ESINFO(ERROR) << "Cannot select the interval.";
			break;
		case 3:
			pickElements(beginElement, endElement, elements, region);
			break;
		}
	} break;

	}
}

template <class TElement>
void Block<TElement>::pattern(const std::vector<Element*> &elements, Region *region, const Triple<size_t> &offset, const Triple<size_t> &size, Pattern pattern, size_t psize)
{
	Triple<size_t> divisor(size * block.elements * block.domains / psize);
	Triple<size_t> modulo;

	divisor.z = divisor.z ? divisor.z : 1;

	Triple<size_t> domain, element, begin, index;
	begin = offset * block.elements * block.domains;

	size_t color;
	switch (pattern) {
	case Pattern::CHESSBOARD_BLACK:
		color = 1;
		break;
	case Pattern::CHESSBOARD_WHITE:
		color = 0;
		break;
	}
	size_t eindex = 0;
	for (domain.z = 0; domain.z < block.domains.z; domain.z++) {
		for (domain.y = 0; domain.y < block.domains.y; domain.y++) {
			for (domain.x = 0; domain.x < block.domains.x; domain.x++) {

				for (element.z = 0; element.z < block.elements.z; element.z++) {
					for (element.y = 0; element.y < block.elements.y; element.y++) {
						for (element.x = 0; element.x < block.elements.x; element.x++, eindex += TElement::subelements) {

							index = begin + domain * block.elements + element;
							if ((index / divisor).sum() % 2 == color) {
								for (size_t e = 0; e < TElement::subelements; e++) {
									region->elements().push_back(elements[eindex + e]);
								}
							}

						}
					}
				}

			}
		}
	}
}

}
}



