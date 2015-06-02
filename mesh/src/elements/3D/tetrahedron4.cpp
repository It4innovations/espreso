
#include "tetrahedron4.h"

using namespace mesh;

std::vector<DenseMatrix> Tetra4_dN()
{
	// dN contains [dNr, dNs, dNt]
	std::vector<DenseMatrix> dN(
		Tetrahedron4GPCount,
		DenseMatrix(Point::size(), Tetrahedron4NodesCount)
	);


	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		//  N = [ r, s, t,  1 - r - s - t ];
		DenseMatrix &m = dN[i];

		// dNr = [ 1, 0, 0, -1 ];
		m(0, 0) =  1.0;
		m(0, 1) =  0.0;
		m(0, 2) =  0.0;
		m(0, 3) = -1.0;

		// dNs = [ 0, 1, 0, -1 ];
		m(1, 0) =  0.0;
		m(1, 1) =  1.0;
		m(1, 2) =  0.0;
		m(1, 3) = -1.0;

		// dNs = [ 0, 0, 1, -1 ];
		m(2, 0) =  0.0;
		m(2, 1) =  0.0;
		m(2, 2) =  1.0;
		m(2, 3) = -1.0;
	}

	return dN;
}

std::vector<DenseMatrix> Tetra4_N()
{
	std::vector<DenseMatrix> N(
			Tetrahedron4GPCount,
			DenseMatrix(1, Tetrahedron4NodesCount));

	std::vector<double> rv;
	std::vector<double> sv;
	std::vector<double> tv;

	if (Tetrahedron4GPCount == 4) {
		double _rv[] = {0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105};
		double _sv[] = {0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685};
		double _tv[] = {0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		sv.assign(_sv, _sv + Tetrahedron4GPCount);
		tv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 5) {
		double _rv[] = {0.2500000000000000, 0.5000000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667};
		double _sv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000};
		double _tv[] = {0.2500000000000000, 0.1666666666666667, 0.1666666666666667, 0.5000000000000000, 0.1666666666666667};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		sv.assign(_sv, _sv + Tetrahedron4GPCount);
		tv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 11) {
		double _rv[] = {0.2500000000000000, 0.7857142857142857, 0.0714285714285714, 0.0714285714285714,
		0.0714285714285714, 0.1005964238332008, 0.3994035761667992, 0.3994035761667992,
		0.3994035761667992, 0.1005964238332008, 0.1005964238332008};
		double _sv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
		0.7857142857142857, 0.3994035761667992, 0.1005964238332008, 0.3994035761667992,
		0.1005964238332008, 0.3994035761667992, 0.1005964238332008};
		double _tv[] = {0.2500000000000000, 0.0714285714285714, 0.0714285714285714, 0.7857142857142857,
		0.0714285714285714, 0.3994035761667992, 0.3994035761667992, 0.1005964238332008,
		0.1005964238332008, 0.1005964238332008, 0.3994035761667992};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		sv.assign(_sv, _sv + Tetrahedron4GPCount);
		tv.assign(_tv, _tv + Tetrahedron4GPCount);
	}
	else if (Tetrahedron4GPCount == 15) {
		double _rv[] = {0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
		0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
		0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
		0.0665501535736643, 0.4334498464263357, 0.4334498464263357};
		double _sv[] = {0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
		0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
		0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
		0.4334498464263357, 0.0665501535736643, 0.4334498464263357};
		double _tv[] = {0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
		0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
		0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
		0.4334498464263357, 0.4334498464263357, 0.0665501535736643};
		rv.assign(_rv, _rv + Tetrahedron4GPCount);
		sv.assign(_sv, _sv + Tetrahedron4GPCount);
		tv.assign(_tv, _tv + Tetrahedron4GPCount);
	}


	// N = [r s t  1-r-s-t];
	for (unsigned int i = 0; i < Tetrahedron4GPCount; i++) {
		double r = rv[i];
		double s = sv[i];
		double t = tv[i];

		N[i](0, 0) = r;
		N[i](0, 1) = s;
		N[i](0, 2) = t;
		N[i](0, 3) = 1.0 - r - s - t;
	}

	return N;
}

std::vector<double> Tetra4_Weight()
{
	switch (Tetrahedron4GPCount) {
	case 4: {
		return std::vector<double> (4, 1 / 24.0);
	}
	case 5: {
		std::vector<double> w(5, 3 / 40.0);
		w[0] = - 2 / 15.0;
		return w;
	}
	case 11: {
		std::vector<double> w(11);
		w[0] = -0.013155555555556;
		w[1] = w[2] = w[3] = w[4] = 0.007622222222222;
		w[5] = w[6] = w[7] = w[8] = w[9] = w[10] = 0.024888888888889;
		return w;
	}
	case 15: {
		std::vector<double> w(15);
		w[0] = 0.030283678097089;
		w[1] = w[2] = w[3] = w[4] = 0.006026785714286;
		w[5] = w[6] = w[7] = w[8] = 0.011645249086029;
		w[9] = w[10] = w[11] = w[12] = w[13] = w[14] = 0.010949141561386;
		return w;
	}
	default:
		std::cerr << "Unknown number of Tatrahedron4 GP count\n";
		exit(EXIT_FAILURE);
	}
}

std::vector<DenseMatrix> Tetrahedron4::_dN = Tetra4_dN();
std::vector<DenseMatrix> Tetrahedron4::_N = Tetra4_N();
std::vector<double> Tetrahedron4::_weighFactor = Tetra4_Weight();

bool Tetrahedron4::match(esint *indices, esint n) {

#if ESPRESO_POINT_DIMENSION == 2
	// Tetrahedron4 is 3D element
	return false;
#endif

	if (n != 8) {
		return false;
	}

	if (!Element::match(indices, 2, 3)) {
		return false;
	}
	if (!Element::match(indices, 4, 5)) {
		return false;
	}
	if (!Element::match(indices, 5, 6)) {
		return false;
	}
	if (!Element::match(indices, 6, 7)) {
		return false;
	}

	esint various[4] = { 0, 1, 2, 4 };
	for (esint i = 0; i < 3; i++) {
		for (esint j = i + 1; j < 4; j++) {
			if (Element::match(indices, various[i], various[j])) {
				return false;
			}
		}
	}

	return true;
}

std::vector<esint> Tetrahedron4::getNeighbours(size_t nodeIndex) const
{
	std::vector<esint> result;
	result.reserve(3);
	for (size_t i = 0; i < Tetrahedron4NodesCount; i++) {
		if (i != nodeIndex) {
			result.push_back(_indices[i]);
		}
	}
	return result;
}

std::vector<esint> Tetrahedron4::getFace(size_t face) const
{
	std::vector<esint> result(3);
	switch (face){
	case 0: {
		result[0] = _indices[1];
		result[1] = _indices[0];
		result[2] = _indices[2];
		break;
	}
	case 1: {
		result[0] = _indices[0];
		result[1] = _indices[1];
		result[2] = _indices[3];
		break;
	}
	case 2: {
		result[0] = _indices[1];
		result[1] = _indices[2];
		result[2] = _indices[3];
		break;
	}
	case 3: {
		result[0] = _indices[2];
		result[1] = _indices[0];
		result[2] = _indices[3];
		break;
	}
	}
	return result;
}

Tetrahedron4::Tetrahedron4(esint *indices)
{
	memcpy(_indices, indices, 3 * sizeof(esint));
	_indices[3] = indices[4];
}


