
#include "assembler.h"

using namespace espreso;

double HypreTemperature::sigma = 0;
HypreTemperature::STABILIZATION HypreTemperature::stabilization = HypreTemperature::STABILIZATION::CAU;

void HypreTemperature::init()
{
	if (Hexahedron8::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Hexahedron8.";
	}
	if (Hexahedron20::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Hexahedron20.";
	}
	if (Tetrahedron4::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Tetrahedron4.";
	}
	if (Tetrahedron10::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Tetrahedron10.";
	}
	if (Prisma6::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Prisma6.";
	}
	if (Prisma15::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Prisma15.";
	}
	if (Pyramid5::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Pyramid5.";
	}
	if (Pyramid13::counter()) {
		ESINFO(GLOBAL_ERROR) << "2D advection diffusion does not support Pyramid13.";
	}

	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	_mesh.prepare(
		false,
		config::solver::B0_TYPE == config::solver::B0_TYPEalternative::KERNELS &&
		config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI
	);
}

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
	inv(0, 1) = - detJx * m(0, 1);
	inv(1, 0) = - detJx * m(1, 0);
	inv(1, 1) =   detJx * m(0, 0);
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	bool CAU = HypreTemperature::stabilization == HypreTemperature::STABILIZATION::CAU;
	double sigma = HypreTemperature::sigma;

	DenseMatrix Ce(2, 2), coordinates, J, invJ, dND;
	double detJ;
	DenseMatrix f(1, element->nodes());
	DenseMatrix U(element->nodes(), 2);

	const std::vector<Evaluator*> &heat_sources = element->settings(Property::HEAT_SOURCE);
	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	Ce(0, 0) = material.termalConduction.x;
	Ce(1, 1) = material.termalConduction.y;

	coordinates.resize(element->nodes(), 2);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates().get(element->node(i), subdomain);
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		U(i, 0) = ux.back()->evaluate(p.x, p.y, p.z) * material.density * material.termalCapacity;
		U(i, 1) = uy.back()->evaluate(p.x, p.y, p.z) * material.density * material.termalCapacity;
		for (size_t j = 0; j < heat_sources.size(); j++) {
			f(0, i) += heat_sources[j]->evaluate(p.x, p.y, p.z);
		}
	}

	eslocal Ksize = element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(1, 2), v(1, 2), Re(1, element->nodes());
	double normGradN = 0;

	for (eslocal gp = 0; gp < element->gaussePoints(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->nodes()), b_e_c(1, element->nodes());
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
			double P_e = h_e * norm_u_e / (2 * Ce(0, 0));
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
		if (konst * weighFactor[gp] * detJ != 0) {
			Ke.multiply(b_e, b_e, konst * weighFactor[gp] * detJ, 1, true);
		}
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


void HypreTemperature::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = _mesh.coordinates().localSize(subdomain);
	DenseMatrix Ke;
	std::vector<double> fe;
	const std::vector<eslocal> &partition = _mesh.getPartition();

	for (eslocal p = partition[subdomain]; p < partition[subdomain + 1]; p++) {

		const Element* e = _mesh.elements()[p];
		processElement(Ke, fe, _mesh, subdomain, e);

		std::cout << "FILL STIFFNESS\n";

		for (size_t i = 0; i < e->nodes(); i++) {
			for (size_t j = 0; j < e->nodes(); j++) {
				// TODO: set to HYPRE STIFFNESS
				//_K(row, column) = Ke(i, j);

				std::cout << Ke(i, j) << " ";

				// ------------------
				//feiPtr->setStiffnes() //TODO : implement me
				// ------------------
			}
			std::cout << "\n";
			// TODO: set to HYPRE RHS
			//f[subdomain][row] += fe[i];

			// ------------------
			//feiPtr->setRHS() //TODO : implement me
			// ------------------

		}
		std::cout << "\n";
	}

}




