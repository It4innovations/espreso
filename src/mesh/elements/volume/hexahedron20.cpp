
#include "hexahedron20.h"
#include "../line/line3.h"
#include "../plane/square8.h"

using namespace espreso;

std::vector<Property> Hexahedron20::_DOFElement;
std::vector<Property> Hexahedron20::_DOFFace;
std::vector<Property> Hexahedron20::_DOFEdge;
std::vector<Property> Hexahedron20::_DOFPoint;
std::vector<Property> Hexahedron20::_DOFMidPoint;

static std::vector<std::vector< double> > Hexa20_rst()
{
	std::vector< std::vector<double> > rst(3, std::vector<double>(Hexahedron20GPCount));

	switch (Hexahedron20GPCount) {
	case 8: {
		double v = 0.577350269189625953;
		rst[0] = {  v,  v,  v,  v, -v, -v, -v, -v };
		rst[1] = { -v, -v,  v,  v, -v, -v,  v,  v };
		rst[2] = { -v,  v, -v,  v, -v,  v, -v,  v };
		return rst;
	}
	case 14: {
		double v1 = 0.758786910639329015;
		double v2 = 0.795822425754222018;
		double v3 = 0;
		rst[0] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
		rst[1] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
		rst[2] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };
		return rst;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Hexahedron20 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Hexa20_dN()
{
	std::vector<DenseMatrix> dN(
		Hexahedron20GPCount,
		DenseMatrix(Point::dimension(), Hexahedron20NodesCount)
	);

	std::vector<std::vector< double> > rst = Hexa20_rst();

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = dN[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];


		// dNr - derivation of basis function
		m(0, 0) =  ((s - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(0, 1) =  ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0 - ((s - 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0;
		m(0, 2) = -((s + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 3) = -((s + 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(0, 4) = -((s - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 5) = -((s - 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(0, 6) =  ((s + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(0, 7) =  ((s + 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 + ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(0, 8)  = -(r * (s - 1.0) * (t - 1.0)) / 2.0;
		m(0, 9)  =  ((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 10) =  (r * (s + 1.0) * (t - 1.0)) / 2.0;
		m(0, 11) = -((pow(s, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(0, 12) =  (r * (s - 1.0) * (t + 1.0)) / 2.0;
		m(0, 13) = -((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 14) = -(r * (s + 1.0) * (t + 1.0)) / 2.0;
		m(0, 15) =  ((pow(s, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(0, 16) = -((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 17) =  ((pow(t, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(0, 18) = -((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(0, 19) =  ((pow(t, 2.0) - 1.0) * (s + 1.0)) / 4.0;


		// dNs - derivation of basis function
		m(1, 0) =  ((r - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 1) = -((r + 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(1, 2) = -((r + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(1, 3) =  ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r - 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0;
		m(1, 4) = -((r - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(1, 5) =  ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r + 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0;
		m(1, 6) =  ((r + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(1, 7) =  ((r - 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(1, 8)  = -((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 9)  =  (s * (r + 1.0) * (t - 1.0)) / 2.0;
		m(1, 10) =  ((pow(r, 2.0) - 1.0) * (t - 1.0)) / 4.0;
		m(1, 11) = -(s * (r - 1.0) * (t - 1.0)) / 2.0;
		m(1, 12) =  ((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 13) = -(s * (r + 1.0) * (t + 1.0)) / 2.0;
		m(1, 14) = -((pow(r, 2.0) - 1.0) * (t + 1.0)) / 4.0;
		m(1, 15) =  (s * (r - 1.0) * (t + 1.0)) / 2.0;
		m(1, 16) = -((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(1, 17) =  ((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 18) = -((pow(t, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(1, 19) =  ((pow(t, 2.0) - 1.0) * (r - 1.0)) / 4.0;

		// dNt - derivation of basis function
		m(2, 0) =  ((r - 1.0) * (s - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 1) = -((r + 1.0) * (s - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
		m(2, 2) =  ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r + 1.0) * (s + 1.0) * (r + s - t - 2.0)) / 8.0;
		m(2, 3) = -((r - 1.0) * (s + 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
		m(2, 4) =  ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r - 1.0) * (s - 1.0) * (r + s - t + 2.0)) / 8.0;
		m(2, 5) = -((r + 1.0) * (s - 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
		m(2, 6) =  ((r + 1.0) * (s + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
		m(2, 7) =  ((r - 1.0) * (s + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;

		m(2, 8)  = -((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 9)  =  ((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 10) =  ((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 11) = -((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 12) =  ((pow(r, 2.0) - 1.0) * (s - 1.0)) / 4.0;
		m(2, 13) = -((pow(s, 2.0) - 1.0) * (r + 1.0)) / 4.0;
		m(2, 14) = -((pow(r, 2.0) - 1.0) * (s + 1.0)) / 4.0;
		m(2, 15) =  ((pow(s, 2.0) - 1.0) * (r - 1.0)) / 4.0;
		m(2, 16) = -(t * (r - 1.0) * (s - 1.0)) / 2.0;
		m(2, 17) =  (t * (r + 1.0) * (s - 1.0)) / 2.0;
		m(2, 18) = -(t * (r + 1.0) * (s + 1.0)) / 2.0;
		m(2, 19) =  (t * (r - 1.0) * (s + 1.0)) / 2.0;
	}

	return dN;
}

static std::vector<DenseMatrix> Hexa20_N() {
	std::vector<DenseMatrix> N(
		Hexahedron20GPCount,
		DenseMatrix(1, Hexahedron20NodesCount)
	);

	std::vector< std::vector< double> > rst = Hexa20_rst();

	for (unsigned int i = 0; i < Hexahedron20GPCount; i++) {
		DenseMatrix &m = N[i];

		double r = rst[0][i];
		double s = rst[1][i];
		double t = rst[2][i];

		// basis function
		m(0, 0) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 - t) * (-r - s - t - 2.0));
		m(0, 1) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 - t) * ( r - s - t - 2.0));
		m(0, 2) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 - t) * ( r + s - t - 2.0));
		m(0, 3) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 - t) * (-r + s - t - 2.0));
		m(0, 4) = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 + t) * (-r - s + t - 2.0));
		m(0, 5) = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 + t) * ( r - s + t - 2.0));
		m(0, 6) = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 + t) * ( r + s + t - 2.0));
		m(0, 7) = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 + t) * (-r + s + t - 2.0));

		m(0, 8) =  0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 - t));
		m(0, 9) =  0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 10) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 - t));
		m(0, 11) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 - t));
		m(0, 12) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 - s) * (1.0 + t));
		m(0, 13) = 0.25 * ((1.0 + r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 14) = 0.25 * ((1.0 - pow(r, 2.0)) * (1.0 + s) * (1.0 + t));
		m(0, 15) = 0.25 * ((1.0 - r) * (1.0 - pow(s, 2.0)) * (1.0 + t));
		m(0, 16) = 0.25 * ((1.0 - r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 17) = 0.25 * ((1.0 + r) * (1.0 - s) * (1.0 - pow(t, 2.0)));
		m(0, 18) = 0.25 * ((1.0 + r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
		m(0, 19) = 0.25 * ((1.0 - r) * (1.0 + s) * (1.0 - pow(t, 2.0)));
	}

	return N;
}

static std::vector<double> Hexa20_weight()
{
	switch (Hexahedron20GPCount) {
	case 8: {
		return std::vector<double> (8, 1.0);
	}
	case 14: {
		std::vector<double> w(14, 0.0);
		double WF_scale_1 = 0.335180055401662;
		double WF_scale_2 = 0.886426592797784;
		for (int i = 0; i < 8; i++) {
			w[i] = (i < 8) ? WF_scale_1 : WF_scale_2;
		}
		return w;
	}
	default:
		ESINFO(ERROR) << "Unknown number of Hexahedron20 GP count.";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Hexahedron20::_dN = Hexa20_dN();
std::vector<DenseMatrix> Hexahedron20::_N = Hexa20_N();
std::vector<double> Hexahedron20::_weighFactor = Hexa20_weight();

bool Hexahedron20::match(const eslocal *indices, eslocal n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Hexahedron20 is 3D element
	return false;
#endif

	switch (n) {
	case Hexahedron20NodesCount:
		for (eslocal i = 0; i < Hexahedron20NodesCount - 1; i++) {
			for (eslocal j = i + 1; j < Hexahedron20NodesCount; j++) {
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

std::vector<eslocal> Hexahedron20::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result;
	if (nodeIndex > 8) {
		result.resize(2);
	} else {
		result.resize(3);
	}

	switch (nodeIndex) {
	case 0: {
		result[0] = _indices[11];
		result[1] = _indices[8];
		result[2] = _indices[16];
		return result;
	}
	case 1: {
		result[0] = _indices[8];
		result[1] = _indices[9];
		result[2] = _indices[17];
		return result;
	}
	case 2: {
		result[0] = _indices[9];
		result[1] = _indices[10];
		result[2] = _indices[18];
		return result;
	}
	case 3: {
		result[0] = _indices[10];
		result[1] = _indices[11];
		result[2] = _indices[19];
		return result;
	}
	case 4: {
		result[0] = _indices[15];
		result[1] = _indices[12];
		result[2] = _indices[16];
		return result;
	}
	case 5: {
		result[0] = _indices[12];
		result[1] = _indices[13];
		result[2] = _indices[17];
		return result;
	}
	case 6: {
		result[0] = _indices[13];
		result[1] = _indices[14];
		result[2] = _indices[18];
		return result;
	}
	case 7: {
		result[0] = _indices[14];
		result[1] = _indices[15];
		result[2] = _indices[19];
		return result;
	}
	case 8: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		return result;
	}
	case 9: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		return result;
	}
	case 10: {
		result[0] = _indices[2];
		result[1] = _indices[3];
		return result;
	}
	case 11: {
		result[0] = _indices[3];
		result[1] = _indices[0];
		return result;
	}
	case 12: {
		result[0] = _indices[4];
		result[1] = _indices[5];
		return result;
	}
	case 13: {
		result[0] = _indices[5];
		result[1] = _indices[6];
		return result;
	}
	case 14: {
		result[0] = _indices[6];
		result[1] = _indices[7];
		return result;
	}
	case 15: {
		result[0] = _indices[7];
		result[1] = _indices[4];
		return result;
	}
	case 16: {
		result[0] = _indices[0];
		result[1] = _indices[4];
		return result;
	}
	case 17: {
		result[0] = _indices[1];
		result[1] = _indices[5];
		return result;
	}
	case 18: {
		result[0] = _indices[1];
		result[1] = _indices[6];
		return result;
	}
	case 19: {
		result[0] = _indices[3];
		result[1] = _indices[7];
		return result;
	}
	}
	return result;
}

void Hexahedron20::fillEdges()
{
	eslocal line[Line3NodesCount];
	_edges.reserve(Hexahedron20EdgeCount);

	for (size_t edge = 0; edge < 4; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 4];
		line[2] = _indices[ edge + 8     ];
		if (3 * edge < _edges.size()) {
			if (_edges[3 * edge] == NULL) {
				_edges[3 * edge] = new Line3(line);
			}
		} else {
			_edges.push_back(new Line3(line));
		}
		_edges[3 * edge]->parentElements().push_back(this);

		line[0] = _indices[ edge          +  4];
		line[1] = _indices[(edge + 1) % 4 +  4];
		line[2] = _indices[ edge          + 12];
		if (3 * edge + 1 < _edges.size()) {
			if (_edges[3 * edge + 1] == NULL) {
				_edges[3 * edge + 1] = new Line3(line);
			}
		} else {
			_edges.push_back(new Line3(line));
		}
		_edges[3 * edge + 1]->parentElements().push_back(this);

		line[0] = _indices[edge     ];
		line[1] = _indices[edge +  4];
		line[2] = _indices[edge + 16];
		if (3 * edge + 2 < _edges.size()) {
			if (_edges[3 * edge + 2] == NULL) {
				_edges[3 * edge + 2] = new Line3(line);
			}
		} else {
			_edges.push_back(new Line3(line));
		}
		_edges[3 * edge + 2]->parentElements().push_back(this);
	}
}

void Hexahedron20::fillFaces()
{
	eslocal square[Square8NodesCount];
	_faces.reserve(Hexahedron20FacesCount);

	for (size_t face = 0; face < 4; face++) {
		square[0] = _indices[ face               ];
		square[1] = _indices[(face + 1) % 4      ];
		square[2] = _indices[(face + 1) % 4 + 4  ];
		square[3] = _indices[ face + 4           ];

		square[4] = _indices[ face          + 8  ];
		square[5] = _indices[(face + 1) % 4 + 16 ];
		square[6] = _indices[ face          + 12 ];
		square[7] = _indices[ face          + 16 ];
		_faces.push_back(new Square8(square));
		_faces.back()->parentElements().push_back(this);
	}

	square[0] = _indices[0];
	square[1] = _indices[3];
	square[2] = _indices[2];
	square[3] = _indices[1];

	square[4] = _indices[11];
	square[5] = _indices[10];
	square[6] = _indices[9];
	square[7] = _indices[8];
	_faces.push_back(new Square8(square));
	_faces.back()->parentElements().push_back(this);

	square[0] = _indices[4];
	square[1] = _indices[5];
	square[2] = _indices[6];
	square[3] = _indices[7];

	square[4] = _indices[12];
	square[5] = _indices[13];
	square[6] = _indices[14];
	square[7] = _indices[15];
	_faces.push_back(new Square8(square));
	_faces.back()->parentElements().push_back(this);
}

void Hexahedron20::setFace(Element* face)
{
	ESINFO(GLOBAL_ERROR) << "Set face";
}

void Hexahedron20::setEdge(Element* edge)
{
	ESINFO(GLOBAL_ERROR) << "Set edge";
}

Point Hexahedron20::faceNormal(const Element *face) const
{
	ESINFO(GLOBAL_ERROR) << "compute normal";
	return Point();
}

Point Hexahedron20::edgeNormal(const Element *edge, const Coordinates &coordinates) const
{
	ESINFO(GLOBAL_ERROR) << "compute normal";
	return Point();
}

Hexahedron20::Hexahedron20(const eslocal *indices, eslocal n, const eslocal *params)
{
	switch (n) {
	case 20:
		memcpy(_indices, indices, Hexahedron20NodesCount * sizeof(eslocal));
		break;
	default:
		ESINFO(ERROR) << "It is not possible to create Hexahedron20 from " << n << " indices.";
	}

	memcpy(_params, params, PARAMS_SIZE * sizeof(eslocal));
}

Hexahedron20::Hexahedron20(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	is.read(reinterpret_cast<char *>(_params), sizeof(eslocal) * PARAMS_SIZE);
}


