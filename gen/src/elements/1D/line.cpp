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

void Line::fillNeighbour(BoundaryNodes &nodes, int indexing) const
{
	const idx_t *indices;
	switch (indexing){
		case Element::LOCAL: {
			indices = _indices;
			break;
		}
		case Element::GLOBAL: {
			indices = _localIndices;
			break;
		}
		default:
			fprintf(stderr, "Incorrect indexing.\n");
			exit(-1);
	}

	if (indices[0] < indices[1]) {
		nodes[indices[0]].insert(indices[1]);
	} else {
		nodes[indices[1]].insert(indices[0]);
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
	if (indices[0] < indices[1]) {
		_indices[0] = indices[0];
		_indices[1] = indices[1];
	} else {
		_indices[0] = indices[1];
		_indices[1] = indices[0];
	}
	_localIndices[0] = _indices[0];
	_localIndices[1] = _indices[1];
}




