
#include "block.h"

namespace espreso {
namespace input {

template <class TElement>
void Block<TElement>::points(std::vector<Point> &points)
{
	Triple<size_t> nodes = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) + 1);
	Triple<double> step = (block.end - block.start) / Triple<double>(nodes.x, nodes.y, nodes.z ? nodes.z : 1);

	points.reserve((nodes + 1).mul());
	for (size_t z = 0; z <= nodes.z; z++) {
		for (size_t y = 0; y <= nodes.y; y++) {
			for (size_t x = 0; x <= nodes.x; x++) {
				points.push_back(Point(block.start.x + x * step.x, block.start.y + y * step.y, block.start.z + z * step.z));
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
	Triple<size_t> elems = block.domains * block.elements;
	Triple<size_t> nodes = Triple<size_t>(TElement::subnodes) + 1;
	Triple<size_t> size = (elems * nodes + 1).toSize();;

	Triple<size_t> _start = start;
	Triple<size_t> _end = end;

	if (_start.x == _end.x) {
		if (_start.x == 0) {
			_end.x++;
		} else {
			_start.x--;
		}
	}
	if (_start.y == _end.y) {
		if (_start.y == 0) {
			_end.y++;
		} else {
			_start.y--;
		}
	}
	if (_start.z == _end.z) {
		if (_start.z == 0) {
			_end.z++;
		} else {
			_start.z--;
		}
	}

	std::vector<eslocal> indices((nodes + 1).mul());

	Triple<size_t> e;
	for (e.z = _start.z; e.z < _end.z; e.z++) {
		for (e.y = _start.y; e.y < _end.y; e.y++) {
			for (e.x = _start.x; e.x < _end.x; e.x++) {

				size_t index = 0;
				Triple<size_t> eOffset, offset = e * nodes;
				for (eOffset.z = 0; eOffset.z <= nodes.z; eOffset.z++) {
					for (eOffset.y = 0; eOffset.y <= nodes.y; eOffset.y++) {
						for (eOffset.x = 0; eOffset.x <= nodes.x; eOffset.x++) {
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
void Block<TElement>::elements(std::vector<Element*> &elements)
{
	elements.reserve(TElement::subelements * (block.domains * block.elements).mul());
	std::vector<eslocal> params(6);

	Triple<size_t> start(0, 0, 0), end = block.domains * block.elements;

	forEachElement(start, end, [&] (std::vector<eslocal> &indices){
		TElement::addElements(elements, indices.data(), params.data());
	});
}


template <class TElement>
void Block<TElement>::boundaries(std::vector<Element*> &nodes, const std::vector<int> &neighbours)
{
	std::vector<int> sorted(neighbours.size());
	std::iota(sorted.begin(), sorted.end(), 0);
	std::sort(sorted.begin(), sorted.end(), [&] (int i, int j) { return neighbours[i] < neighbours[j]; });

	Triple<size_t> count = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) + 1);
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
void Block<TElement>::region(const std::vector<Element*> &nodes, Region &region, const BlockBorder &border, size_t dimension)
{
	if (!border.intersect(block)) {
		return;
	}
	BlockBorder bborder = border;
	if (bborder.start.x < block.start.x) { bborder.start.x = block.start.x; bborder.excludeStart.x = false; }
	if (bborder.start.y < block.start.y) { bborder.start.y = block.start.y; bborder.excludeStart.y = false; }
	if (bborder.start.z < block.start.z) { bborder.start.z = block.start.z; bborder.excludeStart.z = false; }
	if (bborder.end.x > block.end.x) { bborder.end.x = block.end.x; bborder.excludeEnd.x = false; }
	if (bborder.end.y > block.end.y) { bborder.end.y = block.end.y; bborder.excludeEnd.y = false; }
	if (bborder.end.z > block.end.z) { bborder.end.z = block.end.z; bborder.excludeEnd.z = false; }

	Triple<size_t> cnodes = block.domains * block.elements * (Triple<size_t>(TElement::subnodes) + 1);
	Triple<double> step = (block.end - block.start) / Triple<double>(cnodes.x, cnodes.y, cnodes.z ? cnodes.z : 1);
	Triple<size_t> elems = block.domains * block.elements;
	Triple<double> estep = (block.end - block.start) / Triple<double>(elems.x, elems.y, elems.z ? elems.z : 1);

	Triple<double> offsetMin = (bborder.start - block.start) / step;
	Triple<double> offsetMax = (bborder.end - block.start) / step;
	Triple<size_t> start = ((bborder.start - block.start) / estep).floor();
	Triple<size_t> end = ((bborder.end - block.start) / estep).ceil();

	if (offsetMin.x == std::round(offsetMin.x)) { offsetMin.x = std::round(offsetMin.x) + (bborder.excludeStart.x ? 1 : 0); }
	if (offsetMin.y == std::round(offsetMin.y)) { offsetMin.y = std::round(offsetMin.y) + (bborder.excludeStart.y ? 1 : 0); }
	if (offsetMin.z == std::round(offsetMin.z)) { offsetMin.z = std::round(offsetMin.z) + (bborder.excludeStart.z ? 1 : 0); }

	if (offsetMax.x == std::round(offsetMax.x)) { offsetMax.x = std::round(offsetMax.x) - (bborder.excludeEnd.x ? 1 : 0); }
	if (offsetMax.y == std::round(offsetMax.y)) { offsetMax.y = std::round(offsetMax.y) - (bborder.excludeEnd.y ? 1 : 0); }
	if (offsetMax.z == std::round(offsetMax.z)) { offsetMax.z = std::round(offsetMax.z) - (bborder.excludeEnd.z ? 1 : 0); }

	Triple<size_t> minOffset = offsetMin.ceil();
	Triple<size_t> maxOffset = offsetMax.floor();

	switch (bborder.dimension()) {

	case 0: { // Point

		ESINFO(GLOBAL_ERROR) << "Implement selection of points";
	} break;

	case 1: { // Line

		ESINFO(GLOBAL_ERROR) << "Implement selection of edges";
	} break;

	case 2: { // Face
		CubeFace face = bborder.getFace(block);
		switch (dimension) {
		case 0:
			forEachElement(start, end,
			[&] (std::vector<eslocal> &indices) {
				TElement::pickNodes(nodes, region.elements, indices.data(), face);
			},
			[&] (Triple<size_t> &offset) {
				offset.x = offset.x < minOffset.x ? minOffset.x : maxOffset.x < offset.x ? maxOffset.x : offset.x;
				offset.y = offset.y < minOffset.y ? minOffset.y : maxOffset.y < offset.y ? maxOffset.y : offset.y;
				offset.z = offset.z < minOffset.z ? minOffset.z : maxOffset.z < offset.z ? maxOffset.z : offset.z;
			});
			std::sort(region.elements.begin(), region.elements.end());
			Esutils::removeDuplicity(region.elements);
			break;
		case 1:
			ESINFO(GLOBAL_ERROR) << "Implement selection of edges on face";
			break;
		case 2:
			forEachElement(start, end,
			[&] (std::vector<eslocal> &indices) {
				TElement::addFaces(region.elements, indices.data(), face);
			});
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Cannot select element of dimension " << dimension << " on 2D plane.";
		}
	} break;

	case 3: { // Volume

		ESINFO(GLOBAL_ERROR) << "Implement selection of volume";
	} break;

	}
}

}
}



