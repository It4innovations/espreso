
#include "square4.h"
#include "../line/line2.h"

using namespace espreso;

size_t Square4::_counter = 0;

std::vector<Property> Square4::_DOFElement;
std::vector<Property> Square4::_DOFFace;
std::vector<Property> Square4::_DOFEdge;
std::vector<Property> Square4::_DOFPoint;
std::vector<Property> Square4::_DOFMidPoint;

static std::vector<DenseMatrix> get_dN() {
	std::vector<DenseMatrix> dN(
		Square4GPCount,
		DenseMatrix(2, Square4NodesCount)
	);

	double CsQ_scale = 0.577350269189626;

	for (unsigned int i = 0; i < Square4GPCount; i++) {
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		///dN contains [dNr, dNs, dNt]
		DenseMatrix &m = dN[i];

		// dNs - derivation of basis function
		m(0, 0) = 0.25 * ( t - 1);
		m(0, 1) = 0.25 * (-t + 1);
		m(0, 2) = 0.25 * ( t + 1);
		m(0, 3) = 0.25 * (-t - 1);

		// dNt - derivation of basis function
		m(1, 0) = 0.25 * ( s - 1);
		m(1, 1) = 0.25 * (-s - 1);
		m(1, 2) = 0.25 * ( s + 1);
		m(1, 3) = 0.25 * (-s + 1);

//		m(2, 0) = 0;
//		m(2, 1) = 0;
//		m(2, 2) = 0;
//		m(2, 3) = 0;
	}

	return dN;
}

static std::vector<DenseMatrix> get_N() {
	std::vector<DenseMatrix> N(
		Square4GPCount,
		DenseMatrix(1, Square4NodesCount)
	);

	double CsQ_scale = 0.5773502691896258;

	for (unsigned int i = 0; i < Square4GPCount; i++) {
		double s = (i & 2) ? CsQ_scale : -CsQ_scale;
		double t = (i & 1) ? CsQ_scale : -CsQ_scale;

		// basis function
		N[i](0, 0) = 0.25 * (1 - s) * (1 - t);
		N[i](0, 1) = 0.25 * (s + 1) * (1 - t);
		N[i](0, 2) = 0.25 * (s + 1) * (t + 1);
		N[i](0, 3) = 0.25 * (1 - s) * (t + 1);
	}

	return N;
}

std::vector<DenseMatrix> Square4::_dN = get_dN();
std::vector<DenseMatrix> Square4::_N = get_N();
std::vector<double> Square4::_weighFactor(Square4GPCount, 1);

bool Square4::match(const eslocal *indices, eslocal n)
{
	if (n != 4) {
		return false;
	}

	for (eslocal i = 0; i < Square4NodesCount - 1; i++) {
		for (eslocal j = i + 1; j < Square4NodesCount; j++) {
			if (Element::match(indices, i, j)) {
				return false;
			}
		}
	}

	return true;
}

std::vector<eslocal> Square4::getNeighbours(size_t nodeIndex) const
{
	std::vector<eslocal> result(2);

	result[0] = _indices[(nodeIndex + 1) % 4];
	result[1] = _indices[(nodeIndex + 3) % 4];

	return result;
}

void Square4::fillEdges()
{
	eslocal line[Line2NodesCount];
	_edges.reserve(Square4EdgeCount);

	for (size_t edge = 0; edge < 4; edge++) {
		line[0] = _indices[ edge         ];
		line[1] = _indices[(edge + 1) % 4];
		if (edge < _edges.size()) {
			if (_edges[edge] == NULL) {
				_edges[edge] = new Line2(line);
			}
		} else {
			_edges.push_back(new Line2(line));
		}
		_edges[edge]->parentElements().push_back(this);
	}
}

void Square4::setEdge(Element* edge)
{
	eslocal line[Line2NodesCount];
	_edges.resize(Square4EdgeCount, NULL);

	for (size_t e = 0; e < 4; e++) {
		line[0] = _indices[ e         ];
		line[1] = _indices[(e + 1) % 4];
		if (std::is_permutation(line, line + Line2NodesCount, edge->indices())) {
			if (_edges[e] == NULL) {
				_edges[e] = edge;
				edge->parentElements().push_back(this);
			} else {
				ESINFO(GLOBAL_ERROR) << "Merge element";
			}
			return;
		}
	}
	ESINFO(GLOBAL_ERROR) << "Invalid edge";
}

Point Square4::edgeNormal(const Element *edge, const Coordinates &coordinates)
{
	auto it = std::find(_edges.begin(), _edges.end(), edge);
	size_t i = it - _edges.begin();
	Point s = coordinates[_indices[i]];
	Point e = coordinates[_indices[(i + 1) % Square4NodesCount]];

	Point t = e - s;
	Point n(t.y, -t.x, 0);
	n.normalize();
	return n;
}

Square4::Square4(const eslocal *indices)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
}

Square4::Square4(const eslocal *indices, const eslocal *params)
{
	memcpy(_indices, indices, Square4NodesCount * sizeof(eslocal));
	_params.insert(_params.end(), params, params + PARAMS_SIZE);
}

Square4::Square4(std::ifstream &is)
{
	eslocal params;
	is.read(reinterpret_cast<char *>(_indices), sizeof(eslocal) * nodes());
	is.read(reinterpret_cast<char *>(&params), sizeof(eslocal));
	if (params) {
		_params.resize(params);
		is.read(reinterpret_cast<char *>(_params.data()), sizeof(eslocal) * params);
	}
}

