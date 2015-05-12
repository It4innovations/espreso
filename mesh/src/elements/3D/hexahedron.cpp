#include "hexahedron.h"

std::vector<std::vector<double> > Hexa_dN() {
	std::vector<std::vector<double> > dN(HexahedronGPCount);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < HexahedronGPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		dN[i].resize(Point::size() * HexahedronNodesCount);

		// dNr - derivation of basis function
		dN[i][0] = 0.125 * (-(1 - s) * (1 - t));
		dN[i][1] = 0.125 * ( (1 - s) * (1 - t));
		dN[i][2] = 0.125 * ( (1 + s) * (1 - t));
		dN[i][3] = 0.125 * (-(1 + s) * (1 - t));
		dN[i][4] = 0.125 * (-(1 - s) * (1 + t));
		dN[i][5] = 0.125 * ( (1 - s) * (1 + t));
		dN[i][6] = 0.125 * ( (1 + s) * (1 + t));
		dN[i][7] = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		dN[i][8]  = 0.125 * (-(1 - r) * (1 - t));
		dN[i][9]  = 0.125 * (-(1 + r) * (1 - t));
		dN[i][10] = 0.125 * ( (1 + r) * (1 - t));
		dN[i][11] = 0.125 * ( (1 - r) * (1 - t));
		dN[i][12] = 0.125 * (-(1 - r) * (1 + t));
		dN[i][13] = 0.125 * (-(1 + r) * (1 + t));
		dN[i][14] = 0.125 * ( (1 + r) * (1 + t));
		dN[i][15] = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		dN[i][16] = 0.125 * (-(1 - r) * (1 - s));
		dN[i][17] = 0.125 * (-(1 + r) * (1 - s));
		dN[i][18] = 0.125 * (-(1 + r) * (1 + s));
		dN[i][19] = 0.125 * (-(1 - r) * (1 + s));
		dN[i][20] = 0.125 * ( (1 - r) * (1 - s));
		dN[i][21] = 0.125 * ( (1 + r) * (1 - s));
		dN[i][22] = 0.125 * ( (1 + r) * (1 + s));
		dN[i][23] = 0.125 * ( (1 - r) * (1 + s));
	}

	return dN;
}

std::vector<std::vector<double> > Hexa_N() {
	std::vector<std::vector<double> > N(HexahedronGPCount);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < HexahedronGPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i].resize(HexahedronNodesCount);
		N[i][0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
		N[i][1] = 0.125 * (r + 1) * (1 - s) * (1 - t);
		N[i][2] = 0.125 * (r + 1) * (s + 1) * (1 - t);
		N[i][3] = 0.125 * (1 - r) * (s + 1) * (1 - t);
		N[i][4] = 0.125 * (1 - r) * (1 - s) * (t + 1);
		N[i][5] = 0.125 * (r + 1) * (1 - s) * (t + 1);
		N[i][6] = 0.125 * (r + 1) * (s + 1) * (t + 1);
		N[i][7] = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	return N;
}

std::vector<std::vector<double> > Hexahedron::_dN = Hexa_dN();
std::vector<std::vector<double> > Hexahedron::_N = Hexa_N();
std::vector<double> Hexahedron::_weighFactor(HexahedronNodesCount, 1);

bool Hexahedron::match(idx_t *indices, idx_t n) {

#ifndef D3
	// Hexahedron is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

	for (idx_t i = 0; i < 7; i++) {
		for (idx_t j = i + 1; j < 8; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

void Hexahedron::fillNeighbour(BoundaryNodes &nodes) const
{
	idx_t r, c;
	for (idx_t i = 0; i < 4; i++) {
		// BOTTOM FACE
		r = _indices[i];
		c = _indices[(i + 1) % 4];
		if (r < c) {
			nodes[r].insert(c);
		}
		c = _indices[(i + 3) % 4];
		if (r < c) {
			nodes[r].insert(c);
		}
		c = _indices[i + 4];
		if (r < c) {
			nodes[r].insert(c);
		}

		// TOP FACE
		r = _indices[i + 4];
		c = _indices[(i + 1) % 4 + 4];
		if (r < c) {
			nodes[r].insert(c);
		}
		c = _indices[(i + 3) % 4 + 4];
		if (r < c) {
			nodes[r].insert(c);
		}
		c = _indices[i];
		if (r < c) {
			nodes[r].insert(c);
		}
	}
}

void Hexahedron::setFaceNodes(idx_t nodes[], idx_t face) const
{
	if (face == 4) {
		nodes[0] = _indices[0];
		nodes[1] = _indices[1];
		nodes[2] = _indices[2];
		nodes[3] = _indices[3];
		return;
	}
	if (face == 5) {
		nodes[0] = _indices[4];
		nodes[1] = _indices[5];
		nodes[2] = _indices[6];
		nodes[3] = _indices[7];
		return;
	}
	nodes[0] = _indices[ face              ];
	nodes[1] = _indices[(face + 1) % 4     ];
	nodes[2] = _indices[(face + 1) % 4 + 4 ];
	nodes[3] = _indices[ face + 4          ];
}

void Hexahedron::fillFaces(BoundaryFaces &faces, int part) const
{
	Square *s;
	idx_t ids[4];
	for (idx_t i = 0; i < HexahedronFacesCount; i++) {
		setFaceNodes(ids, i);
		s = new Square(ids);
		s->fillFaces(faces, part);
		delete s;
	}
}

void Hexahedron::fillFacesOnBorder(BoundaryFaces &faces, const BoundaryNodes &nodes, int part) const
{
	Square *s;
	idx_t ids[4];
	for (idx_t i = 0; i < HexahedronFacesCount; i++) {
		setFaceNodes(ids, i);
		s = new Square(ids);
		s->fillFacesOnBorder(faces, nodes, part);
		delete s;
	}
}

void Hexahedron::fillLines(BoundaryLines &lines, int parts[]) const
{
	Line* line;
	idx_t ids[2];
	for (int i = 0; i < 4; i++) {
		// bottom
		ids[0] = _indices[i];
		ids[1] = _indices[(i + 1) % 4];
		line = new Line(ids);
		line->fillLines(lines, parts);
		delete line;

		// middle
		ids[0] = _indices[i];
		ids[1] = _indices[i + 4];
		line = new Line(ids);
		line->fillLines(lines, parts);
		delete line;

		// upper
		ids[0] = _indices[i + 4];
		ids[1] = _indices[(i + 1) % 4 + 4];
		line = new Line(ids);
		line->fillLines(lines, parts);
		delete line;
	}
}

Hexahedron::Hexahedron(idx_t *indices)
{
	idx_t min = 0;
	for (idx_t i = 1; i < HexahedronNodesCount; i++) {
		if (indices[min] > indices[i]) {
			min = i;
		}
	}

	// rotate to move minimal index to position 0 or 4
	idx_t offset = min % 4;
	for (int i = 0; i < 4; i++) {
		_indices[i] = indices[(i + offset) % 4];
		_indices[i + 4] = indices[(i + offset) % 4 + 4];
	}

	// rotate to move minimal index to position 0
	idx_t tmp[8];
	if (min > 3) {
		tmp[0] = _indices[4];
		tmp[1] = _indices[7];
		tmp[2] = _indices[6];
		tmp[3] = _indices[5];
		tmp[4] = _indices[0];
		tmp[5] = _indices[3];
		tmp[6] = _indices[2];
		tmp[7] = _indices[1];
	} else {
		memcpy(tmp, _indices, sizeof(idx_t) * HexahedronNodesCount);
	}


	idx_t right = tmp[1];
	idx_t left = tmp[3];
	idx_t up = tmp[4];
	if (left < up && right < up) {
		memcpy(_indices, tmp, sizeof(idx_t) * HexahedronNodesCount);
	}
	if (left < right && up < right) {
		_indices[0] = tmp[0];
		_indices[1] = tmp[3];
		_indices[2] = tmp[7];
		_indices[3] = tmp[4];
		_indices[4] = tmp[1];
		_indices[5] = tmp[2];
		_indices[6] = tmp[6];
		_indices[7] = tmp[5];
	}
	if (right < left && up < left) {
		_indices[0] = tmp[0];
		_indices[1] = tmp[4];
		_indices[2] = tmp[5];
		_indices[3] = tmp[1];
		_indices[4] = tmp[3];
		_indices[5] = tmp[7];
		_indices[6] = tmp[6];
		_indices[7] = tmp[2];
	}
}


