#include "square.h"

// TODO: Implement base functions
std::vector<std::vector<double> > Square::_dN;
std::vector<std::vector<double> > Square::_N;
std::vector<double> Square::_weighFactor;

bool Square::match(idx_t *indices, idx_t n)
{
	if (n != 4) {
		return false;
	}

	for (idx_t i = 0; i < SquareNodesCount - 1; i++) {
		for (idx_t j = i + 1; j < SquareNodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

void Square::fillNeighbour(BoundaryNodes &nodes) const
{
	idx_t r, c;
	for (idx_t i = 0; i < SquareNodesCount; i++) {
		r = _indices[i];
		c = _indices[(i + 1) % 4];
		if (r < c) {
			nodes[r].insert(c);
		}
		c = _indices[(i + 3) % 4];
		if (r < c) {
			nodes[r].insert(c);
		}
	}
}

void Square::fillFaces(BoundaryFaces &faces, int part) const
{
	BoundaryFaces::iterator it = faces.find((Element*)this);
	if (it == faces.end()) {
		faces.insert(std::pair<Element*, std::pair<int, int> >(
		                 new Square(*this), std::pair<int, int>(part, -1)
		             ));
	} else {
		it->second.second = part;
	}
}

void Square::fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const
{
	idx_t positions[4] = { 0, 1, 2, 3 };
	if (isOnBorder(nodes, positions, 4)) {
		fillFaces(faces, part);
	}
}

void Square::fillLines(BoundaryLines &lines, int parts[]) const
{
	Line* line;
	idx_t ids[2];
	for (int i = 0; i < SquareNodesCount; i++) {
		ids[0] = _indices[i];
		ids[1] = _indices[(i + 1) % SquareNodesCount];
		line = new Line(ids);
		line->fillLines(lines, parts);
		delete line;
	}
}

Square::Square(idx_t *indices)
{
	idx_t min = 0;
	for (idx_t i = 1; i < SquareNodesCount; i++) {
		if (indices[min] > indices[i]) {
			min = i;
		}
	}
	_indices[0] = indices[min];
	if (indices[(min + 1) % 4] < indices[(min + 3) % 4]) {
		_indices[1] = indices[(min + 1) % 4];
		_indices[2] = indices[(min + 2) % 4];
		_indices[3] = indices[(min + 3) % 4];
	} else {
		_indices[1] = indices[(min + 3) % 4];
		_indices[2] = indices[(min + 2) % 4];
		_indices[3] = indices[(min + 1) % 4];
	}
}



