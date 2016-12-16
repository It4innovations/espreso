
#include "assembler.h"

namespace espreso {

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

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computePlaneCorners(1, true, true);
			break;
		case B0_TYPE::KERNELS:
			_mesh.computeEdgesSharedByDomains();
			break;
		default:
			break;
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.loadNodeProperty(_configuration.temperature.values        , { "T" }     , { Property::TEMPERATURE });
	_mesh.loadNodeProperty(_configuration.thickness.values          , { }         , { Property::THICKNESS });

	_mesh.loadProperty(_configuration.translation_motions.values      , { "x", "y" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y });
	_mesh.loadProperty(_configuration.initial_temperature.values      , { }         , { Property::INITIAL_TEMPERATURE });
	_mesh.loadProperty(_configuration.heat_source.values              , { "T" }     , { Property::HEAT_SOURCE });
	switch (_configuration.heat_flux) {
	case AdvectionDiffusion2DConfiguration::HEAT_FLUX::GENERAL:
		_mesh.loadProperty(_configuration.general_heat_flux.values, { "q" }, { Property::GENERAL_HEAT_FLUX });
		break;
	case AdvectionDiffusion2DConfiguration::HEAT_FLUX::CONVECTIVE:
		_mesh.loadProperty(_configuration.heat_transfer_coefficient.values, { "HTC" }  , { Property::HEAT_TRANSFER_COEFFICIENT });
		_mesh.loadProperty(_configuration.external_temperature.values     , { "T_ext" }, { Property::EXTERNAL_TEMPERATURE });
		break;
	}

	_mesh.loadMaterials(_configuration.materials.configurations, _configuration.material_set.values);
}

void AdvectionDiffusion2D::saveMeshProperties(store::Store &store)
{
	store.storeProperty("translationMotion", { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y }, store::Store::ElementType::ELEMENTS);
	store.storeProperty("headSource", { Property::HEAT_SOURCE }, store::Store::ElementType::ELEMENTS);
	store.storeProperty("temperature", { Property::TEMPERATURE }, store::Store::ElementType::NODES);
}

void AdvectionDiffusion2D::saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("temperature", 1, results, store::Store::ElementType::NODES);
	store.finalize();
}

void AdvectionDiffusion2D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
}

void AdvectionDiffusion2D::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
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

static void processEdge(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* edge, const AdvectionDiffusion2DConfiguration &configuration)
{
	DenseMatrix coordinates(edge->nodes(), 2), dND(1, 2), q(edge->nodes(), 1), htc(edge->nodes(), 1), thickness(edge->nodes(), 1);
	DenseMatrix gpQ(1, 2), gpHtc(1, 1), gpThickness(1, 1);

	Ke.resize(2, 2);
	Ke = 0;

	const std::vector<DenseMatrix> &dN = edge->dN();
	const std::vector<DenseMatrix> &N = edge->N();
	const std::vector<double> &weighFactor = edge->weighFactor();

	for (size_t n = 0; n < edge->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[edge->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[edge->node(n)].y;
		switch (configuration.heat_flux) {
		case AdvectionDiffusion2DConfiguration::HEAT_FLUX::GENERAL:
			q(n, 0) = edge->settings(Property::GENERAL_HEAT_FLUX).back()->evaluate(edge->node(n));
			break;
		case AdvectionDiffusion2DConfiguration::HEAT_FLUX::CONVECTIVE:
			q(n, 0) = edge->settings(Property::EXTERNAL_TEMPERATURE).back()->evaluate(edge->node(n));
			q(n, 0) *= edge->settings(Property::HEAT_TRANSFER_COEFFICIENT).back()->evaluate(edge->node(n));
			htc(n, 0) = edge->settings(Property::HEAT_TRANSFER_COEFFICIENT).back()->evaluate(edge->node(n));
			thickness(n, 0) = mesh.nodes()[edge->node(n)]->settings(Property::THICKNESS).back()->evaluate(edge->node(n));
			break;
		}
		q(n, 0) *= mesh.nodes()[edge->node(n)]->settings(Property::THICKNESS).back()->evaluate(edge->node(n));
	}

	eslocal Ksize = 2 * edge->nodes();
	fe.resize(Ksize);
	std::fill(fe.begin(), fe.end(), 0);

	for (size_t gp = 0; gp < edge->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		switch (configuration.heat_flux) {
		case AdvectionDiffusion2DConfiguration::HEAT_FLUX::CONVECTIVE:
			gpHtc.multiply(N[gp], htc);
			gpThickness.multiply(N[gp], thickness);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0) * gpThickness(0, 0), 1, true);
			break;
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % edge->nodes()) * gpQ(0, i / edge->nodes());
		}
	}
}


static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element, const AdvectionDiffusion2DConfiguration &configuration)
{
	bool CAU = configuration.stabilization == AdvectionDiffusion2DConfiguration::STABILIZATION::CAU;
	double sigma = configuration.sigma;

	DenseMatrix Ce(2, 2), coordinates, J, invJ, dND;
	double detJ;
	DenseMatrix f(1, element->nodes());
	DenseMatrix U(element->nodes(), 2);
	DenseMatrix thickness(element->nodes(), 1), dens(element->nodes(), 1), Cp(element->nodes(), 1), K(element->nodes(), 4);
	DenseMatrix gpThickness(1, 1), gpDens(1, 1), gpCp(1, 1), gpK(1, 4);

	const std::vector<Evaluator*> &heat_sources = element->settings(Property::HEAT_SOURCE);
	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 2);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		thickness(i, 0) = mesh.nodes()[element->node(i)]->settings(Property::THICKNESS).back()->evaluate(element->node(i));
		U(i, 0) = ux.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i)) * thickness(i, 0);
		U(i, 1) = uy.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i)) * thickness(i, 0);
		for (size_t j = 0; j < heat_sources.size(); j++) {
			f(0, i) += heat_sources[j]->evaluate(element->node(i)) * thickness(i, 0);
		}

		K(i, 0) = 1; // TODO KXX
		K(i, 1) = 1; // TODO KYY
		K(i, 2) = 0; // TODO KXY
		K(i, 3) = 0; // TODO KYX
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

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

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

		Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp] * gpThickness(0, 0), 1, true);
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

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R1, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R1, norm, defect, subdomain);
}

void AdvectionDiffusion2D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e, _configuration);
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
	switch (mtype) {
		case SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
				K[subdomain].RemoveLower();
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], RegMat[subdomain], subdomain);
			}
			break;
		case SparseMatrix::MatrixType::REAL_UNSYMMETRIC:
			for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
			}
			break;
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

		processElement(Ke, fe, _mesh, elements[e], _configuration);

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

	for (size_t i = 0; i < _mesh.edges().size(); i++) {
		if (_mesh.edges()[i]->inDomain(subdomain) && _mesh.edges()[i]->settings().size()) {
			processEdge(Ke, fe, _mesh, _mesh.edges()[i], _configuration);

			for (size_t nx = 0; nx < _mesh.edges()[i]->nodes(); nx++) {
				for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
					size_t row = nodes[_mesh.edges()[i]->node(nx)]->DOFIndex(subdomain, dx);
					f[subdomain][row] += fe[dx * _mesh.edges()[i]->nodes() + nx];
				}
			}
			switch (_configuration.heat_flux) {
			case AdvectionDiffusion2DConfiguration::HEAT_FLUX::CONVECTIVE:
				for (size_t nx = 0; nx < _mesh.edges()[i]->nodes(); nx++) {
					for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
						size_t row = nodes[_mesh.edges()[i]->node(nx)]->DOFIndex(subdomain, dx);
						for (size_t ny = 0; ny < _mesh.edges()[i]->nodes(); ny++) {
							for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
								size_t column = nodes[_mesh.edges()[i]->node(ny)]->DOFIndex(subdomain, dy);
								_K(row, column) = Ke(dx * _mesh.edges()[i]->nodes() + nx, dy * _mesh.edges()[i]->nodes() + ny);
							}
						}
					}
				}
			}

		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

static void postProcessElement(std::vector<double> &gradient, std::vector<double> &flux, DenseMatrix &solution, const Element* element, const Mesh &mesh, const AdvectionDiffusion2DConfiguration &configuration)
{
	bool CAU = configuration.stabilization == AdvectionDiffusion2DConfiguration::STABILIZATION::CAU;
	double sigma = configuration.sigma;

	DenseMatrix Ce(2, 2), coordinates, J, invJ, dND;
	double detJ;
	DenseMatrix U(element->nodes(), 2);
	DenseMatrix thickness(element->nodes(), 1), dens(element->nodes(), 1), Cp(element->nodes(), 1), K(element->nodes(), 4);
	DenseMatrix gpThickness(1, 1), gpDens(1, 1), gpCp(1, 1), gpK(1, 4);

	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	DenseMatrix matFlux(2, 1), matGradient(2, 1);

	coordinates.resize(element->nodes(), 2);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		thickness(i, 0) = mesh.nodes()[element->node(i)]->settings(Property::THICKNESS).back()->evaluate(element->node(i));
		U(i, 0) = ux.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i)) * thickness(i, 0);
		U(i, 1) = uy.back()->evaluate(element->node(i)) * material.density(element->node(i)) * material.termalCapacity(element->node(i)) * thickness(i, 0);

		K(i, 0) = 1; // TODO KXX
		K(i, 1) = 1; // TODO KYY
		K(i, 2) = 0; // TODO KXY
		K(i, 3) = 0; // TODO KYX
	}

	DenseMatrix u(1, 2), v(1, 2), Re(1, element->nodes());

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant2x2(J);
		inverse(J, invJ, detJ);

		gpThickness.multiply(N[gp], thickness);
		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(0, 1) = gpK(0, 2);
		Ce(1, 0) = gpK(0, 3);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->nodes()), b_e_c(1, element->nodes());
		b_e.multiply(u, dND, 1, 0);


		double norm_u_e = u.norm();
		double h_e = 0;

		if (norm_u_e != 0) {
			h_e = 2 * norm_u_e / b_e.norm();
		}

		Ce(0, 0) += sigma * h_e * norm_u_e;
		Ce(1, 1) += sigma * h_e * norm_u_e;

		matGradient.multiply(dND, solution, gpThickness(0, 0), 1);
		matFlux.multiply(Ce, dND * solution, gpThickness(0, 0), 1);
	}
	gradient.push_back(matGradient(0, 0) / element->gaussePoints());
	gradient.push_back(matGradient(1, 0) / element->gaussePoints());
	flux.push_back(matFlux(0, 0) / element->gaussePoints());
	flux.push_back(matFlux(1, 0) / element->gaussePoints());
}

void AdvectionDiffusion2D::postProcess(store::Store &store, const std::vector<std::vector<double> > &solution)
{
	std::vector<std::vector<double> > termalGradient(_mesh.parts()), termalFlux(_mesh.parts());
	DenseMatrix eSolution;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		termalGradient[p].reserve(matrixSize[p]);
		termalFlux[p].reserve(matrixSize[p]);
		for (size_t e = _mesh.getPartition()[p]; e < _mesh.getPartition()[p + 1]; e++) {
			eSolution.resize(_mesh.elements()[e]->nodes(), 1);
			for (size_t n = 0; n < _mesh.elements()[e]->nodes(); n++) {
				eSolution(n, 0) = solution[p][_mesh.nodes()[_mesh.elements()[e]->node(n)]->DOFIndex(p, 0)];
			}
			postProcessElement(termalGradient[p], termalFlux[p], eSolution, _mesh.elements()[e], _mesh, _configuration);
		}
	}

	store.storeValues("gradient", 2, termalGradient, store::Store::ElementType::ELEMENTS);
	store.storeValues("flux", 2, termalFlux, store::Store::ElementType::ELEMENTS);
	store.finalize();

}

}


