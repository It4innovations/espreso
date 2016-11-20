
#include "assembler.h"

namespace espreso {

double AdvectionDiffusion2D::sigma = 0;
AdvectionDiffusion2D::STABILIZATION AdvectionDiffusion2D::stabilization = AdvectionDiffusion2D::STABILIZATION::CAU;

std::vector<Property> AdvectionDiffusion2D::elementDOFs;
std::vector<Property> AdvectionDiffusion2D::faceDOFs;
std::vector<Property> AdvectionDiffusion2D::edgeDOFs;
std::vector<Property> AdvectionDiffusion2D::pointDOFs = { Property::TEMPERATURE };
std::vector<Property> AdvectionDiffusion2D::midPointDOFs = { Property::TEMPERATURE };

void AdvectionDiffusion2D::prepareMeshStructures()
{
	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_mesh.computePlaneCorners(config::mesh::CORNERS, config::mesh::VERTEX_CORNERS, config::mesh::EDGE_CORNERS);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_mesh.computeEdgesSharedByDomains();
			break;
		default:
			break;
		}
	}
}

void AdvectionDiffusion2D::saveMeshProperties(output::Store &store)
{
	store.storeProperty("translationMotion", { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y }, output::Store::ElementType::ELEMENTS);
	store.storeProperty("headSource", { Property::HEAT_SOURCE }, output::Store::ElementType::ELEMENTS);
	store.storeProperty("temperature", { Property::TEMPERATURE }, output::Store::ElementType::NODES);
}

void AdvectionDiffusion2D::saveMeshResults(output::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("temperature", 1, results, output::Store::ElementType::NODES);
	store.finalize();
}

void AdvectionDiffusion2D::assembleGluingMatrices()
{
	_constraints.initMatrices(matrixSize);

	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			std::for_each(R1.begin(), R1.end(), [] (SparseMatrix &m) { m.ConvertCSRToDense(0); });
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.edges(), pointDOFs, R1);
			break;
		default:
			break;
		}
	}
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

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element)
{
	bool CAU = AdvectionDiffusion2D::stabilization == AdvectionDiffusion2D::STABILIZATION::CAU;
	double sigma = AdvectionDiffusion2D::sigma;

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

	Ce(0, 0) = material.termalConductionX(0);
	Ce(1, 1) = material.termalConductionY(0);

	coordinates.resize(element->nodes(), 2);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		U(i, 0) = ux.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i));
		U(i, 1) = uy.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i));
		for (size_t j = 0; j < heat_sources.size(); j++) {
			f(0, i) += heat_sources[j]->evaluate(element->node(i));
		}
	}

	eslocal Ksize = element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(1, 2), v(1, 2), Re(1, element->nodes());
	double normGradN = 0;

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
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
		double C_e = 0;

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

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R1, SparseMatrix &R2, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernels_from_nonsym_K(K, RegMat, R1, R2, norm, defect, subdomain);
}

void AdvectionDiffusion2D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e);
	dofs.resize(e->nodes());
	for (size_t n = 0; n < e->nodes(); n++) {
		dofs[n] = e->node(n);
	}

	for (size_t n = 0; n < e->nodes(); n++) {
		if (_mesh.nodes()[e->node(n)]->settings().isSet(Property::TEMPERATURE)) {
			fe[n] = _mesh.nodes()[e->node(n)]->settings(Property::TEMPERATURE).back()->evaluate(e->node(n)) / _mesh.nodes()[e->node(n)]->domains().size();
		}
	}
}

void AdvectionDiffusion2D::makeStiffnessMatricesRegular()
{
	for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
		algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
	}
}

void AdvectionDiffusion2D::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal e = partition[subdomain]; e < partition[subdomain + 1]; e++) {

		processElement(Ke, fe, _mesh, elements[e]);

		for (size_t nx = 0; nx < elements[e]->nodes(); nx++) {
			for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
				size_t row = nodes[elements[e]->node(nx)]->DOFIndex(subdomain, dx);
				for (size_t ny = 0; ny < elements[e]->nodes(); ny++) {
					for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
						size_t column = nodes[elements[e]->node(ny)]->DOFIndex(subdomain, dy);
						_K(row, column) = Ke(dx * elements[e]->nodes() + nx, dy * elements[e]->nodes() + ny);
					}
				}
				f[subdomain][row] += fe[dx * elements[e]->nodes() + nx];
			}
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

}


