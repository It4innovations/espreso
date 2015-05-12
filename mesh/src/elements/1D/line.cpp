#include "line.h"

// TODO: Implement base functions
std::vector<std::vector<double> > Line::_dN;
std::vector<std::vector<double> > Line::_N;
std::vector<double> Line::_weighFactor;

bool Line::match(idx_t *indices, idx_t n)
{
	if (n != 2) {
		return false;
	}

	if (Element::match(indices, 0, 1)) {
		return false;
	}

	return true;
}

void Line::fillNeighbour(BoundaryNodes &nodes) const
{
	if (_indices[0] < _indices[1]) {
		nodes[_indices[0]].insert(_indices[1]);
	} else {
		nodes[_indices[1]].insert(_indices[0]);
	}
}

void Line::fillFaces(BoundaryFaces &faces, int part) const
{
	return;
}

void Line::fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const
{
	return;
}

void Line::fillLines(BoundaryLines &lines, int parts[]) const
{
	BoundaryLines::iterator it = lines.find((Element*)this);
	if (it == lines.end()) {
		lines[new Line(*this)].insert(parts, parts + 2);
	} else {
		it->second.insert(parts, parts + 2);
	}
}

Line::Line(idx_t *indices)
{
	memcpy(_indices, indices, LineNodesCount * sizeof(idx_t));
}




