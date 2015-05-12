#include "triangle.h"

// TODO: Implement base functions
std::vector<std::vector<double> > Triangle::_dN;
std::vector<std::vector<double> > Triangle::_N;
std::vector<double> Triangle::_weighFactor;

bool Triangle::match(idx_t *indices, idx_t n)
{
	if (n != 3) {
		return false;
	}

	for (idx_t i = 0; i < TriangleNodesCount - 1; i++) {
		for (idx_t j = i + 1; j < TriangleNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

void Triangle::fillNeighbour(BoundaryNodes &nodes) const
{
	idx_t r, c;
	for (idx_t i = 0; i < TriangleNodesCount; i++) {
		for (idx_t j = 0; j < TriangleNodesCount; j++) {
			r = _indices[i];
			c = _indices[i];
			if (r < c) {
				nodes[r].insert(c);
			}
		}
	}
}

void Triangle::fillFaces(BoundaryFaces &faces, int part) const
{
	BoundaryFaces::iterator it = faces.find((Element*)this);
	if (it == faces.end()) {
		faces.insert(std::pair<Element*, std::pair<int, int> >(
		                 new Triangle(*this), std::pair<int, int>(part, -1)
		             ));
	} else {
		it->second.second = part;
	}
}

void Triangle::fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const
{
	idx_t positions[3] = { 0, 1, 2 };
	if (isOnBorder(nodes, positions, 3)) {
		fillFaces(faces, part);
	}
}

void Triangle::fillLines(BoundaryLines &lines, int parts[]) const
{
	Line* line;
	idx_t ids[2];
	for (int i = 0; i < TriangleNodesCount; i++) {
		ids[0] = _indices[i];
		ids[1] = _indices[(i + 1) % TriangleNodesCount];
		line = new Line(ids);
		line->fillLines(lines, parts);
		delete line;
	}
}

Triangle::Triangle(idx_t *indices)
{
	memcpy(_indices, indices, sizeof(idx_t) * TriangleNodesCount);
	std::sort(_indices, _indices + 3);
}

