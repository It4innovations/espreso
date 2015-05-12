#include "element.h"
#include <iomanip>
#include <iostream>

void Element::fillNodes(idx_t *nodes) const
{
	for (size_t i = 0; i < size(); i++) {
		nodes[i] = node(i);
	}
}

void Element::setLocalIndices(std::vector<idx_t> &mapping)
{
	for (size_t i = 0; i < size(); i++) {
		indices()[i] = mapping[node(i)];
	}
}

void Element::fillBoundaries(BoundaryNodes &nodes, int part) const
{
	for (size_t i = 0; i < size(); i++) {
		nodes[node(i)].insert(part);
	}
}

void Element::coordinatesToVector(
		std::vector<double> &vector,
		const Coordinates &coordinates,
		IndicesType indicesType,
		size_t part) const
{
	for (size_t i = 0; i < size(); i++) {
		if (indicesType == Element::GLOBAL) {
			&vector[i * Point::size()] << coordinates[node(i)];
		}
		if (indicesType == Element::LOCAL) {
			&vector[i * Point::size()] << coordinates.localPoint(part, node(i));
		}
	}
}

void Element::_elaticity(
		std::vector<double> &Ke,
		std::vector<double> &Me,
		std::vector<double> &fe,
		std::vector<double> &coordinates,
		std::vector<double> &inertia,
		double ex,
		double mi,
		bool dynamic
	) const
{
	const std::vector<std::vector<double> > &dN = this->dN();
	const std::vector<std::vector<double> > &N = this->N();
	const std::vector<double> &weighFactor = this->weighFactor();

	int dimension = Point::size();
	int Ksize = dimension * size();
	int Csize = 6;	// TODO: even for D2??
	int nodes = size();
	int gausePoints = gpSize();

	std::vector<double> MatC(Csize * Csize, 0.0);

	double E = ex / ((1 + mi) * (1 - 2 * mi));
	double mi2 = E * (1.0 - mi);
	double mi3 = E * (0.5 - mi);
	mi = E * mi;

	MatC[0]  = mi2; MatC[1]  = mi;  MatC[2]  = mi;
	MatC[6]  = mi;  MatC[7]  = mi2; MatC[8]  = mi;
	MatC[12] = mi;  MatC[13] = mi;  MatC[14] = mi2;

	MatC[21] = mi3; MatC[28] = mi3; MatC[35] = mi3;


	std::vector<double> MatJ (dimension * dimension, 0);
	std::vector<double> invJ (dimension * dimension, 0);

	std::vector<double> dND (Ksize, 0);
	std::vector<double> CB (Ksize * Csize, 0);
	std::vector<double> B (Ksize * Csize, 0);

	Ke.resize(Ksize * Ksize);
	fe.resize(Ksize);
	fill(Ke.begin(), Ke.end(), 0);
	fill(fe.begin(), fe.end(), 0);
	if (dynamic) {
		Me.resize(nodes * nodes);
		fill(Me.begin(), Me.end(), 0);
	}

	for (int gp = 0; gp < gausePoints; gp++) {

		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, dimension, nodes,
			1, &dN[gp][0], nodes, &coordinates[0], dimension,
			0, &MatJ[0], dimension
		);

		double detJ = fabs( MatJ[0] * MatJ[4] * MatJ[8] +
							MatJ[1] * MatJ[5] * MatJ[6] +
							MatJ[2] * MatJ[3] * MatJ[7] -
							MatJ[2] * MatJ[4] * MatJ[6] -
							MatJ[1] * MatJ[3] * MatJ[8] -
							MatJ[0] * MatJ[5] * MatJ[7]);

		double detJx = 1 / detJ;

		invJ[0] = detJx * (  MatJ[8] * MatJ[4] - MatJ[7] * MatJ[5] );
		invJ[1] = detJx * (- MatJ[8] * MatJ[1] + MatJ[7] * MatJ[2] );
		invJ[2] = detJx * (  MatJ[5] * MatJ[1] - MatJ[4] * MatJ[2] );
		invJ[3] = detJx * (- MatJ[8] * MatJ[3] + MatJ[6] * MatJ[5] );
		invJ[4] = detJx * (  MatJ[8] * MatJ[0] - MatJ[6] * MatJ[2] );
		invJ[5] = detJx * (- MatJ[5] * MatJ[0] + MatJ[3] * MatJ[2] );
		invJ[6] = detJx * (  MatJ[7] * MatJ[3] - MatJ[6] * MatJ[4] );
		invJ[7] = detJx * (- MatJ[7] * MatJ[0] + MatJ[6] * MatJ[1] );
		invJ[8] = detJx * (  MatJ[4] * MatJ[0] - MatJ[3] * MatJ[1] );

		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, nodes, dimension,
			1, &invJ[0], dimension, &dN[gp][0], nodes,
			0, &dND.front(), nodes
		);

		// TODO: block ordering inside B
		int columns = Ksize;
		const double *dNDx = &dND[0];
		const double *dNDy = &dND[size()];
		const double *dNDz = &dND[2 * size()];

		// B =
		// dX   0   0
		//  0  dY   0
		//  0   0  dZ
		// dY  dX   0
		//  0  dZ  dY
		// dZ   0  dX

		memcpy(&B[0],                        dNDx, sizeof(double) * size());
		memcpy(&B[3 * columns + size()],     dNDx, sizeof(double) * size());
		memcpy(&B[5 * columns + 2 * size()], dNDx, sizeof(double) * size());

		memcpy(&B[1 * columns + size()],     dNDy, sizeof(double) * size());
		memcpy(&B[3 * columns],              dNDy, sizeof(double) * size());
		memcpy(&B[4 * columns + 2 * size()], dNDy, sizeof(double) * size());

		memcpy(&B[2 * columns + 2 * size()], dNDz, sizeof(double) * size());
		memcpy(&B[4 * columns + size()],     dNDz, sizeof(double) * size());
		memcpy(&B[5 * columns],              dNDz, sizeof(double) * size());


		// C * B
		cblas_dgemm(
			CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Csize, Ksize, Csize,
			1, &MatC[0], Csize, &B[0], Ksize,
			0, &CB.front(), Ksize
		);
		//Ke = Ke + (B' * (C * B)) * dJ * WF;
		cblas_dgemm(
			CblasRowMajor, CblasTrans, CblasNoTrans,
			Ksize, Ksize, Csize,
			detJ * weighFactor[gp], &B[0], Ksize, &CB[0], Ksize,
			1, &Ke[0], Ksize
		);

		for (int i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp][i % nodes] * inertia[i / nodes];
		}

		if (dynamic) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			double dense = 7.85e-9;
			cblas_dgemm(
				CblasRowMajor, CblasTrans, CblasNoTrans,
				nodes, nodes, 1,
				dense * detJ * weighFactor[gp], &N[gp][0], nodes, &N[gp][0], nodes,
				1, &Me[0], nodes
			);
		}
	}
}

void Element::_addLocalValues(
	SparseVVPMatrix &K,
	SparseVVPMatrix &M,
	std::vector<double> &f,
	const std::vector<double> &Ke,
	const std::vector<double> &Me,
	const std::vector<double> &fe,
	int offset,
	bool dynamic) const
{
	// Element ordering: xxxx, yyyy, zzzz,...
	// Global ordering:  xyz, xyz, xyz, xyz, ...
	size_t row, column;
	size_t s = Point::size();
	for (size_t i = 0; i < s * size(); i++) {
		row = s * (node(i % size()) - offset) + i / size();
		for (size_t j = 0; j < s * size(); j++) {
			column = s * (node(j % size()) - offset) + j / size();
			K(row, column) += Ke[i * s * size() + j];
		}
		f[row] += fe[i];
	}
	if (!dynamic) {
		return;
	}
	for (size_t i = 0; i < size(); i++) {
		row = s * (node(i) - offset);
		for (size_t j = 0; j < size(); j++) {
			column = s * (node(i) - offset);
			for (size_t k = 0; k < s; k++) {
				M(row + k, column + k) += Me[i * size() + j];
			}
		}
	}
}

bool Element::isOnBorder(const BoundaryNodes &nodes, const idx_t *positions, idx_t n) const
{
	const idx_t* _indices = indices();
	std::vector<idx_t> result(nodes[_indices[positions[0]]].begin(), nodes[_indices[positions[0]]].end());
	std::vector<idx_t>::iterator it = result.end();

	for (size_t i = 1; i < n; i++) {
		std::set<idx_t> tmp(result.begin(), it);
		it = std::set_intersection(
		         tmp.begin(), tmp.end(),
		         nodes[_indices[positions[i]]].begin(), nodes[_indices[positions[i]]].end(),
		         result.begin());
		if (it - result.begin() == 1) {
			return false;
		}
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, const Element &e)
{
	for (size_t i = 0; i < e.size(); i++) {
		os << e.node(i) << " ";
	}
	return os;
}

inline void operator<<(double *nodeArray, const Element &e)
{
	for (size_t i = 0; i < e.size(); i++) {
		nodeArray[i] = e.node(i);
	}
}

