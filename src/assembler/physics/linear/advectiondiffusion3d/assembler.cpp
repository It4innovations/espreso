
#include "assembler.h"

using namespace espreso;

static double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
   );
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	DenseMatrix Ce(3, 3), coordinates, J, invJ, dND;
	double detJ, inertia;

	const Material &material = mesh.materials()[element->getParam(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	Ce(0, 0) = material.termalConduction.x;
	Ce(1, 1) = material.termalConduction.y;
	Ce(2, 2) = material.termalConduction.z;
	inertia = 0;

	coordinates.resize(element->size(), 3);
	for (size_t i = 0; i < element->size(); i++) {
		coordinates(i, 0) = mesh.coordinates().get(element->node(i), subdomain).x;
		coordinates(i, 1) = mesh.coordinates().get(element->node(i), subdomain).y;
		coordinates(i, 2) = mesh.coordinates().get(element->node(i), subdomain).z;
	}

	eslocal Ksize = element->size();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(3, 1);
	u(0, 0) = 0;
	u(1, 0) = 0;
	u(2, 0) = 0;
	double sigma = 0;

	for (eslocal gp = 0; gp < element->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->size());
		b_e.multiply(u, dND, 1, 0, true);
		double norm_u_e = pow(pow(u(0, 0), 2) + pow(u(1, 0), 2) + pow(u(2, 0), 2), 0.5);
		double h_e = 0, tau_e = 0, konst = 0;
		if (norm_u_e != 0) {
			double nn = 0;
			for (size_t i = 0; i < element->size(); i++) {
				nn += b_e(0, i) * b_e(0, i);
			}
			h_e = 2 * norm_u_e / sqrt(nn);
			double P_e = h_e * norm_u_e / 2 * Ce(0, 0);
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);
		}

		Ce(0, 0) += sigma * h_e * norm_u_e;
		Ce(1, 1) += sigma * h_e * norm_u_e;
		Ce(2, 2) += sigma * h_e * norm_u_e;

		Ke.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
		Ke.multiply(N[0], b_e, detJ * weighFactor[gp], 1, true);
		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i) * inertia;
			if (norm_u_e != 0) {
				fe[i] += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * inertia / (2 * norm_u_e);
			}
		}
	}
}

void AdvectionDiffusion3D::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = _mesh.coordinates().localSize(subdomain);
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.getElements();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* e = _mesh.getElements()[i];
		processElement(Ke, fe, _mesh, subdomain, e);

		for (size_t i = 0; i < e->size(); i++) {
			eslocal row = e->node(i);
			for (size_t j = 0; j < e->size(); j++) {
				eslocal column = e->node(j);
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}



