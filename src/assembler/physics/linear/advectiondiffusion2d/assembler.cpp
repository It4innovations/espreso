
#include "assembler.h"

using namespace espreso;

static double determinant2x2(DenseMatrix &m)
{
	return fabs(
		m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)
	);
}

static void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	inv.resize(m.rows(), m.columns());
	double detJx = 1 / det;
	inv(0, 0) =   detJx * m(1, 1);
	inv(0, 1) = - detJx * m(1, 0);
	inv(1, 0) = - detJx * m(0, 1);
	inv(1, 1) =   detJx * m(0, 0);
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	bool CAU = true;
	DenseMatrix Ce(2, 2), coordinates, J, invJ, dND;
	double detJ;
	DenseMatrix f(1, element->size());
	DenseMatrix U(element->size(), 2);

	const std::vector<Evaluator*> &heat_sources = element->settings(Property::HEAT_SOURCE);
	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);

	const Material &material = mesh.materials()[element->getParam(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	Ce(0, 0) = material.termalConduction.x;
	Ce(1, 1) = material.termalConduction.y;

	coordinates.resize(element->size(), 2);
	for (size_t i = 0; i < element->size(); i++) {
		const Point &p = mesh.coordinates().get(element->node(i), subdomain);
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		U(i, 0) = ux.back()->evaluate(p.x, p.y, p.z) * material.density * material.termalCapacity;
		U(i, 1) = uy.back()->evaluate(p.x, p.y, p.z) * material.density * material.termalCapacity;
		for (size_t j = 0; j < heat_sources.size(); j++) {
			f(0, i) += heat_sources[j]->evaluate(p.x, p.y, p.z);
		}
	}

	eslocal Ksize = element->size();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(1, 2), v(1, 2), Re(1, element->size());
	double sigma = 0;
	double normGradN = 0;

	for (eslocal gp = 0; gp < element->gpSize(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->size()), b_e_c(1, element->size());
		b_e.multiply(u, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < Re.columns(); i++) {
					Re(0, i) = b_e(0, i) - f(0, i);
				}
				DenseMatrix ReBt(1, 2);
				ReBt.multiply(Re, dND, 1 / pow(normGradN, 2), 0, false, true);
				for (size_t i = 0; i < ReBt.columns(); i++) {
					v(0, i) = u(0, i) - ReBt(0, i);
				}
			} else {
				v = u;
			}
		}


		double norm_u_e = u.norm();
		double h_e = 0, tau_e = 0, konst = 0;
		double C_e;

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
			double P_e = h_e * norm_u_e / 2 * Ce(0, 0);
			tau_e = std::max(0.0, 1 - 1 / P_e);
			konst = h_e * tau_e / (2 * norm_u_e);

			if (CAU) {
				DenseMatrix u_v(1, 2);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
				b_e_c.multiply(u_v, dND, 1, 0);
				double norm_u_v = u_v.norm();
				double h_e_c = 2 * norm_u_v / b_e_c.norm();
				double P_e_c = h_e_c * norm_u_v / (2 * Ce.norm());
				double tau_e_c = std::max(0.0, 1 - 1 / P_e_c);

				double konst1 = Re.norm() / normGradN;
				double konst2 = tau_e * h_e != 0 ? tau_e_c * h_e_c / (tau_e * h_e) : 0;
				if (konst1 / norm_u_e < konst2) {
					C_e = tau_e * h_e * konst1 * (konst2 - konst1 / norm_u_e) / 2;
				} else {
					C_e = 0;
				}
			}
		}

		Ce(0, 0) += sigma * h_e * norm_u_e;
		Ce(1, 1) += sigma * h_e * norm_u_e;

		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
		Ke.multiply(N[gp], b_e, detJ * weighFactor[gp], 1, true);
		Ke.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
		if (CAU) {
			Ke.multiply(dND, dND, C_e * weighFactor[gp] * detJ, 1, true);
		}

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i) * f(0, i);
			if (norm_u_e != 0) {
				fe[i] += detJ * weighFactor[gp] * h_e * tau_e * b_e(0, i) * f(0, i) / (2 * norm_u_e);
			}
		}
	}
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R1, SparseMatrix &R2, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernels_from_nonsym_K(K, RegMat, R1, R2, norm, defect, subdomain);
}

void AdvectionDiffusion2D::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = _mesh.coordinates().localSize(subdomain);
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();

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


  if (config::info::PRINT_MATRICES){
			std::ofstream osK(Logging::prepareFile(subdomain, "K").c_str());
			osK << K[subdomain];
			osK.close();
  }
	algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
}



