
#include "hexahedron8.h"
#include "../line/line2.h"
#include "../plane/square4.h"

using namespace espreso;

size_t Hexahedron8::_counter = 0;

std::vector<Property> Hexahedron8::_DOFElement;
std::vector<Property> Hexahedron8::_DOFFace;
std::vector<Property> Hexahedron8::_DOFEdge;
std::vector<Property> Hexahedron8::_DOFPoint;
std::vector<Property> Hexahedron8::_DOFMidPoint;

static std::vector<DenseMatrix> Hexa_dN() {
	std::vector<DenseMatrix> dN(
		Hexahedron8GPCount,
		DenseMatrix(Point::dimension(), Hexahedron8NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNr - derivation of basis function
		m(0, 0) = 0.125 * (-(1 - s) * (1 - t));
		m(0, 1) = 0.125 * ( (1 - s) * (1 - t));
		m(0, 2) = 0.125 * ( (1 + s) * (1 - t));
		m(0, 3) = 0.125 * (-(1 + s) * (1 - t));
		m(0, 4) = 0.125 * (-(1 - s) * (1 + t));
		m(0, 5) = 0.125 * ( (1 - s) * (1 + t));
		m(0, 6) = 0.125 * ( (1 + s) * (1 + t));
		m(0, 7) = 0.125 * (-(1 + s) * (1 + t));

		// dNs - derivation of basis function
		m(1, 0)  = 0.125 * (-(1 - r) * (1 - t));
		m(1, 1)  = 0.125 * (-(1 + r) * (1 - t));
		m(1, 2) = 0.125 * ( (1 + r) * (1 - t));
		m(1, 3) = 0.125 * ( (1 - r) * (1 - t));
		m(1, 4) = 0.125 * (-(1 - r) * (1 + t));
		m(1, 5) = 0.125 * (-(1 + r) * (1 + t));
		m(1, 6) = 0.125 * ( (1 + r) * (1 + t));
		m(1, 7) = 0.125 * ( (1 - r) * (1 + t));

		// dNt - derivation of basis function
		m(2, 0) = 0.125 * (-(1 - r) * (1 - s));
		m(2, 1) = 0.125 * (-(1 + r) * (1 - s));
		m(2, 2) = 0.125 * (-(1 + r) * (1 + s));
		m(2, 3) = 0.125 * (-(1 - r) * (1 + s));
		m(2, 4) = 0.125 * ( (1 - r) * (1 - s));
		m(2, 5) = 0.125 * ( (1 + r) * (1 - s));
		m(2, 6) = 0.125 * ( (1 + r) * (1 + s));
		m(2, 7) = 0.125 * ( (1 - r) * (1 + s));
	}

	return dN;
}

static std::vector<DenseMatrix> Hexa_N() {
	std::vector<DenseMatrix> N(
		Hexahedron8GPCount,
		DenseMatrix(1, Hexahedron8NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Hexahedron8GPCount; i++) {
		double r = (i & 4) ? CsQ_scale : -CsQ_scale;
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i](0, 0) = 0.125 * (1 - r) * (1 - s) * (1 - t);
		N[i](0, 1) = 0.125 * (r + 1) * (1 - s) * (1 - t);
		N[i](0, 2) = 0.125 * (r + 1) * (s + 1) * (1 - t);
		N[i](0, 3) = 0.125 * (1 - r) * (s + 1) * (1 - t);
		N[i](0, 4) = 0.125 * (1 - r) * (1 - s) * (t + 1);
		N[i](0, 5) = 0.125 * (r + 1) * (1 - s) * (t + 1);
		N[i](0, 6) = 0.125 * (r + 1) * (s + 1) * (t + 1);
		N[i](0, 7) = 0.125 * (1 - r) * (s + 1) * (t + 1);
	}

	return N;
}

std::vector<DenseMatrix> Hexahedron8::_dN = Hexa_dN();
std::vector<DenseMatrix> Hexahedron8::_N = Hexa_N();
std::vector<double> Hexahedron8::_weighFactor(Hexahedron8NodesCount, 1);

bool Hexahedron8::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron8 is 3D element
	return false;
#endif

	switch (n) {
	case Hexahedron8NodesCount:
		for (eslocal i = 0; i < Hexahedron8NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Hexahedron8NodesCount; j++) {
				if (Element::match(indices, i, j)) {
					return false;
				}
			}
		}
		return true;
	default:
		return false;
	}
}

std::vector<eslocal> Hexahedron8::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(3);

	if (nodeIndex > 3) {
		result[0] = _indices[nodeIndex - 4];
		result[1] = _indices[(nodeIndex + 1) % 4 + 4];
		result[2] = _indices[(nodeIndex + 3) % 4 + 4];
	} else {
		result[0] = _indices[nodeIndex + 4];
		result[1] = _indices[(nodeIndex + 1) % 4];
		result[2] = _indices[(nodeIndex + 3) % 4];
	}

	return result;
}

void Hexahedron8::fillEdges()
{
	eslocal line[Line2NodesCount];
	_edges.reserve(Hexahedron8EdgeCount);

	for (size_t edge = 0; edge < 4; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 4];
		if (3 * edge < _edges.size()) {
			if (_edges[3 * edge] == NULL) {
				_edges[3 * edge] = new Line2(line);
			}
		} else {
			_edges.push_back(new Line2(line));
		}
		_edges[3 * edge]->parentElements().push_back(this);

		line[0] = _indices[ edge          +  4];
		line[1] = _indices[(edge + 1) % 4 +  4];
		if (3 * edge + 1 < _edges.size()) {
			if (_edges[3 * edge + 1] == NULL) {
				_edges[3 * edge + 1] = new Line2(line);
			}
		} else {
			_edges.push_back(new Line2(line));
		}
		_edges[3 * edge + 1]->parentElements().push_back(this);

		line[0] = _indices[edge     ];
		line[1] = _indices[edge +  4];
		if (3 * edge + 2 < _edges.size()) {
			if (_edges[3 * edge + 2] == NULL) {
				_edges[3 * edge + 2] = new Line2(line);
			}
		} else {
			_edges.push_back(new Line2(line));
		}
		_edges[3 * edge + 2]->parentElements().push_back(this);
	}
}

void Hexahedron8::fillFaces()
{
	eslocal square[Square4NodesCount];
	_faces.reserve(Hexahedron8FacesCount);

	for (size_t face = 0; face < 4; face++) {
		square[0] = _indices[ face             ];
		square[1] = _indices[(face + 1) % 4    ];
		square[2] = _indices[(face + 1) % 4 + 4];
		square[3] = _indices[ face + 4         ];
		_faces.push_back(new Square4(square));
		_faces.back()->parentElements().push_back(this);
	}

	square[0] = _indices[0];
	square[1] = _indices[3];
	square[2] = _indices[2];
	square[3] = _indices[1];
	_faces.push_back(new Square4(square));
	_faces.back()->parentElements().push_back(this);

	square[0] = _indices[4];
	square[1] = _indices[5];
	square[2] = _indices[6];
	square[3] = _indices[7];
	_faces.push_back(new Square4(square));
	_faces.back()->parentElements().push_back(this);
}

void Hexahedron8::setFace(Element* face)
{
	ESINFO(GLOBAL_ERROR) << "Set face";
}

void Hexahedron8::setEdge(Element* edge)
{
	_edges.push_back(edge);
	//ESINFO(GLOBAL_ERROR) << "Set edge";
}

Point Hexahedron8::faceNormal(const Element *face)
{
	ESINFO(GLOBAL_ERROR) << "compute normal";
	return Point();
}

Point Hexahedron8::edgeNormal(const Element *edge, const Coordinates &coordinates)
{
	ESINFO(GLOBAL_ERROR) << "compute normal";
	return Point();
}

Hexahedron8::Hexahedron8(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 8:
		memcpy(_indices, indices, Hexahedron8NodesCount * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Hexahedron8 from " << n << " elements.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Hexahedron8::Hexahedron8(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}


