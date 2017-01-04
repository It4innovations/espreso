
#include "assembler.h"

namespace espreso {

std::vector<Property> AdvectionDiffusion3D::elementDOFs;
std::vector<Property> AdvectionDiffusion3D::faceDOFs;
std::vector<Property> AdvectionDiffusion3D::edgeDOFs;
std::vector<Property> AdvectionDiffusion3D::pointDOFs = { Property::TEMPERATURE };
std::vector<Property> AdvectionDiffusion3D::midPointDOFs = { Property::TEMPERATURE };

void AdvectionDiffusion3D::prepareMeshStructures()
{
	Hexahedron8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Hexahedron20::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Tetrahedron10::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Prisma15::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid5::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Pyramid13::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	Square4::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Square8::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle3::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);
	Triangle6::setDOFs(elementDOFs, faceDOFs, edgeDOFs, pointDOFs, midPointDOFs);

	matrixSize = _mesh.assignUniformDOFsIndicesToNodes(matrixSize, pointDOFs);
	_mesh.computeNodesDOFsCounters(pointDOFs);

	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			_mesh.computeVolumeCorners(1, true, true, false);
			break;
		case B0_TYPE::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		default:
			break;
		}
	}

	_constraints.initMatrices(matrixSize);

	_mesh.loadNodeProperty(_configuration.temperature.values, { }, { Property::TEMPERATURE });

	_mesh.loadProperty(_configuration.translation_motions.values      , { "X", "Y", "Z" }, { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z });
	_mesh.loadProperty(_configuration.initial_temperature.values      , { }         , { Property::INITIAL_TEMPERATURE });
	_mesh.loadProperty(_configuration.heat_source.values              , { }     , { Property::HEAT_SOURCE });
	_mesh.loadProperty(_configuration.heat_flux.values                , { }     , { Property::HEAT_FLUX });
	_mesh.loadProperty(_configuration.heat_flow.values                , { }     , { Property::HEAT_FLOW });

	for (auto it = _configuration.convection.configurations.begin(); it != _configuration.convection.configurations.end(); ++it) {
		std::map<std::string, std::string> values;
		values[it->first] = it->second->external_temperature;
		_mesh.loadProperty(values, { }     , { Property::EXTERNAL_TEMPERATURE });
		values[it->first] = it->second->heat_transfer_coefficient;
		_mesh.loadProperty(values, { }     , { Property::HEAT_TRANSFER_COEFFICIENT });
	}

	for (size_t r = 0; r < _mesh.regions().size(); r++) {
		if (_mesh.regions()[r].settings.isSet(Property::HEAT_FLOW)) {
			_mesh.regions()[r].computeArea(_mesh.coordinates());
		}
	}

	_mesh.loadMaterials(_configuration.materials.configurations, _configuration.material_set.values);
}

void AdvectionDiffusion3D::assembleB1()
{
	EqualityConstraints::insertDirichletToB1(_constraints, _mesh.nodes(), pointDOFs);
	EqualityConstraints::insertElementGluingToB1(_constraints, _mesh.nodes(), pointDOFs, K);
}

void AdvectionDiffusion3D::assembleB0()
{
	if (_solverConfiguration.method == ESPRESO_METHOD::HYBRID_FETI) {
		switch (_solverConfiguration.B0_type) {
		case B0_TYPE::CORNERS:
			EqualityConstraints::insertDomainGluingToB0(_constraints, _mesh.corners(), pointDOFs);
			break;
		case B0_TYPE::KERNELS:
			std::for_each(R1.begin(), R1.end(), [] (SparseMatrix &m) { m.ConvertCSRToDense(0); });
			EqualityConstraints::insertKernelsToB0(_constraints, _mesh.faces(), pointDOFs, R1);
			break;
		default:
			break;
		}
	}
}

void AdvectionDiffusion3D::saveMeshProperties(store::Store &store)
{
	store.storeProperty("translationMotion", { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z }, store::Store::ElementType::ELEMENTS);
	store.storeProperty("headSource", { Property::HEAT_SOURCE }, store::Store::ElementType::ELEMENTS);
	store.storeProperty("temperature", { Property::TEMPERATURE }, store::Store::ElementType::NODES);
}

void AdvectionDiffusion3D::saveMeshResults(store::Store &store, const std::vector<std::vector<double> > &results)
{
	store.storeValues("temperature", 1, results, store::Store::ElementType::NODES);
}

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

static void processEdge(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* edge, const AdvectionDiffusion3DConfiguration &configuration, double area = 1)
{
	DenseMatrix coordinates(edge->nodes(), 3), dND(1, 3), q(edge->nodes(), 1), htc(edge->nodes(), 1), thickness(edge->nodes(), 1), flow(edge->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpThickness(1, 1), gpFlow(1, 1);

	eslocal Ksize = edge->nodes();
	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize);
	Ke = 0;
	std::fill(fe.begin(), fe.end(), 0);

	const std::vector<DenseMatrix> &dN = edge->dN();
	const std::vector<DenseMatrix> &N = edge->N();
	const std::vector<double> &weighFactor = edge->weighFactor();

	for (size_t n = 0; n < edge->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[edge->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[edge->node(n)].y;
		coordinates(n, 2) = mesh.coordinates()[edge->node(n)].z;
		if (edge->settings().isSet(Property::HEAT_FLOW)) {
			q(n, 0) = edge->settings(Property::HEAT_FLOW).back()->evaluate(edge->node(n)) / area;
		}

		if (edge->settings().isSet(Property::HEAT_FLUX)) {
			q(n, 0) = edge->settings(Property::HEAT_FLUX).back()->evaluate(edge->node(n));
		}

		if (edge->settings().isSet(Property::EXTERNAL_TEMPERATURE)) {
			htc(n, 0) = edge->settings(Property::HEAT_TRANSFER_COEFFICIENT).back()->evaluate(edge->node(n));
			q(n, 0) = htc(n, 0) * edge->settings(Property::EXTERNAL_TEMPERATURE).back()->evaluate(edge->node(n));
		}
	}

	for (size_t gp = 0; gp < edge->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		double J = dND.norm();
		gpQ.multiply(N[gp], q);
		if (edge->settings().isSet(Property::EXTERNAL_TEMPERATURE)) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % edge->nodes()) * gpQ(0, 0);
		}
	}
}

static void processFace(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* face, const AdvectionDiffusion3DConfiguration &configuration, double area = 1)
{
	DenseMatrix coordinates(face->nodes(), 3), dND(1, 3), q(face->nodes(), 1), htc(face->nodes(), 1), flow(face->nodes(), 1);
	DenseMatrix gpQ(1, 1), gpHtc(1, 1), gpFlow(1, 1);

	eslocal Ksize = face->nodes();
	Ke.resize(Ksize, Ksize);
	fe.resize(Ksize);
	Ke = 0;
	std::fill(fe.begin(), fe.end(), 0);

	const std::vector<DenseMatrix> &dN = face->dN();
	const std::vector<DenseMatrix> &N = face->N();
	const std::vector<double> &weighFactor = face->weighFactor();

	for (size_t n = 0; n < face->nodes(); n++) {
		coordinates(n, 0) = mesh.coordinates()[face->node(n)].x;
		coordinates(n, 1) = mesh.coordinates()[face->node(n)].y;
		coordinates(n, 2) = mesh.coordinates()[face->node(n)].z;
		if (face->settings().isSet(Property::HEAT_FLOW)) {
			q(n, 0) = face->settings(Property::HEAT_FLOW).back()->evaluate(face->node(n)) / area;
		}

		if (face->settings().isSet(Property::HEAT_FLUX)) {
			q(n, 0) = face->settings(Property::HEAT_FLUX).back()->evaluate(face->node(n));
		}

		if (face->settings().isSet(Property::EXTERNAL_TEMPERATURE)) {
			htc(n, 0) = face->settings(Property::HEAT_TRANSFER_COEFFICIENT).back()->evaluate(face->node(n));
			q(n, 0) = htc(n, 0) * face->settings(Property::EXTERNAL_TEMPERATURE).back()->evaluate(face->node(n));
		}
	}

	for (size_t gp = 0; gp < face->gaussePoints(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		double J = va.norm();
		gpQ.multiply(N[gp], q);
		if (face->settings().isSet(Property::EXTERNAL_TEMPERATURE)) {
			gpHtc.multiply(N[gp], htc);
			Ke.multiply(N[gp], N[gp], weighFactor[gp] * J * gpHtc(0, 0), 1, true);
		}
		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += J * weighFactor[gp] * N[gp](0, i % face->nodes()) * gpQ(0, 0);
		}
	}
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, const Element* element, const AdvectionDiffusion3DConfiguration &configuration)
{
	bool CAU = configuration.stabilization == AdvectionDiffusion3DConfiguration::STABILIZATION::CAU;
	double sigma = configuration.sigma;

	DenseMatrix Ce(3, 3), coordinates, J, invJ, dND;
	double detJ;
	DenseMatrix f(1, element->nodes());
	DenseMatrix U(element->nodes(), 3);
	DenseMatrix K(element->nodes(), 9), gpK(element->nodes(), 9);

	const std::vector<Evaluator*> &heat_sources = element->settings(Property::HEAT_SOURCE);
	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);
	const std::vector<Evaluator*> &uz = element->settings(Property::TRANSLATION_MOTION_Z);

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	coordinates.resize(element->nodes(), 3);
	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
		U(i, 0) = ux.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));
		U(i, 1) = uy.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));
		U(i, 2) = uz.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));
		for (size_t j = 0; j < heat_sources.size(); j++) {
			f(0, i) += heat_sources[j]->evaluate(element->node(i));
		}

		switch (material->getModel()) {
		case AdvectionDiffusion3DMaterial::ISOTROPIC:
			K(i, 0) = K(i, 1) = K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case AdvectionDiffusion3DMaterial::DIAGONAL:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case AdvectionDiffusion3DMaterial::SYMMETRIC:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = K(i, 5) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i));
			K(i, 4) = K(i, 7) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i));
			K(i, 6) = K(i, 8) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i));
			break;
		case AdvectionDiffusion3DMaterial::ANISOTROPIC:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i));
			K(i, 4) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i));
			K(i, 5) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YX)->evaluate(element->node(i));
			K(i, 6) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i));
			K(i, 7) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZX)->evaluate(element->node(i));
			K(i, 8) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZY)->evaluate(element->node(i));
			break;
		}
	}

	eslocal Ksize = element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	DenseMatrix u(1, 3), v(1, 3), Re(1, element->nodes());
	double normGradN = 0;

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);

		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);

		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

		dND.multiply(invJ, dN[gp]);

		DenseMatrix b_e(1, element->nodes()), b_e_c(1, element->nodes());
		b_e.multiply(u, dND, 1, 0);

		if (CAU) {
			normGradN = dND.norm();
			if (normGradN >= 1e-12) {
				for (size_t i = 0; i < Re.columns(); i++) {
					Re(0, i) = b_e(0, i) - f(0, i);
				}
				DenseMatrix ReBt(1, 3);
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
				DenseMatrix u_v(1, 3);
				u_v(0, 0) = u(0, 0) - v(0, 0);
				u_v(0, 1) = u(0, 1) - v(0, 1);
				u_v(0, 2) = u(0, 2) - v(0, 2);
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
		Ce(2, 2) += sigma * h_e * norm_u_e;

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


static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &R, SparseMatrix &RegMat, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);

}

void AdvectionDiffusion3D::assembleStiffnessMatrix(const Element* e, DenseMatrix &Ke, std::vector<double> &fe, std::vector<eslocal> &dofs) const
{
	processElement(Ke, fe, _mesh, e, _configuration);
	dofs.resize(e->nodes());
	for (size_t n = 0; n < e->nodes(); n++) {
		dofs[n] = e->node(n);
	}
}

void AdvectionDiffusion3D::makeStiffnessMatricesRegular()
{
	switch (mtype) {
		case SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			#pragma omp parallel for
			for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
				K[subdomain].RemoveLower();
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], RegMat[subdomain], subdomain);
			}
			break;
		case SparseMatrix::MatrixType::REAL_UNSYMMETRIC:
			#pragma omp parallel for
			for (size_t subdomain = 0; subdomain < K.size(); subdomain++) {
				algebraicKernelsAndRegularization(K[subdomain], R1[subdomain], R2[subdomain], RegMat[subdomain], subdomain);
			}
			break;
	}
}

void AdvectionDiffusion3D::composeSubdomain(size_t subdomain)
{
	SparseVVPMatrix<eslocal> _K;
	DenseMatrix Ke;
	std::vector<double> fe;

	_K.resize(matrixSize[subdomain], matrixSize[subdomain]);
	f[subdomain].resize(matrixSize[subdomain]);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &nodes = _mesh.nodes();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* e = _mesh.elements()[i];
		processElement(Ke, fe, _mesh, e, _configuration);

		for (size_t i = 0; i < e->nodes(); i++) {
			eslocal row = _mesh.nodes()[e->node(i)]->DOFIndex(subdomain, 0);
			for (size_t j = 0; j < e->nodes(); j++) {
				eslocal column = _mesh.nodes()[e->node(j)]->DOFIndex(subdomain, 0);
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}
	}

	auto processRegion = [&] (const std::vector<Element*> &elments, bool withK = false, double area = 1) {
		for (size_t i = 0; i < elments.size(); i++) {
			if (elments[i]->inDomain(subdomain)) {
				if (elments[0]->type() == Element::Type::PLANE) {
					processFace(Ke, fe, _mesh, elments[i], _configuration, area);
				}
				if (elments[0]->type() == Element::Type::LINE) {
					processEdge(Ke, fe, _mesh, elments[i], _configuration, area);
				}

				for (size_t nx = 0; nx < elments[i]->nodes(); nx++) {
					for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
						size_t row = nodes[elments[i]->node(nx)]->DOFIndex(subdomain, dx);
						f[subdomain][row] += fe[dx * elments[i]->nodes() + nx];
					}
				}
				if (withK) {
					for (size_t nx = 0; nx < elments[i]->nodes(); nx++) {
						for (size_t dx = 0; dx < pointDOFs.size(); dx++) {
							size_t row = nodes[elments[i]->node(nx)]->DOFIndex(subdomain, dx);
							for (size_t ny = 0; ny < elments[i]->nodes(); ny++) {
								for (size_t dy = 0; dy < pointDOFs.size(); dy++) {
									size_t column = nodes[elments[i]->node(ny)]->DOFIndex(subdomain, dy);
									_K(row, column) = Ke(dx * elments[i]->nodes() + nx, dy * elments[i]->nodes() + ny);
								}
							}
						}
					}
				}

			}
		}
	};


	for (size_t r = 0; r < _mesh.regions().size(); r++) {
		if (_mesh.regions()[r].settings.isSet(Property::HEAT_FLUX)) {
			processRegion(_mesh.regions()[r].elements);
		}
		if (_mesh.regions()[r].settings.isSet(Property::HEAT_FLOW)) {
			processRegion(_mesh.regions()[r].elements, false, _mesh.regions()[r].area);
		}
		if (_mesh.regions()[r].settings.isSet(Property::EXTERNAL_TEMPERATURE)) {
			processRegion(_mesh.regions()[r].elements, true);
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	K[subdomain] = csrK;
}

static void postProcessElement(std::vector<double> &gradient, std::vector<double> &flux, DenseMatrix &solution, const Element* element, const Mesh &mesh, const AdvectionDiffusion3DConfiguration &configuration)
{
	bool CAU = configuration.stabilization == AdvectionDiffusion3DConfiguration::STABILIZATION::CAU;
	double sigma = configuration.sigma;

	DenseMatrix Ce(3, 3), coordinates(element->nodes(), 3), J, invJ, dND;
	double detJ;
	DenseMatrix U(element->nodes(), 3);
	DenseMatrix K(element->nodes(), 9);
	DenseMatrix gpK(1, 9);

	const std::vector<Evaluator*> &ux = element->settings(Property::TRANSLATION_MOTION_X);
	const std::vector<Evaluator*> &uy = element->settings(Property::TRANSLATION_MOTION_Y);
	const std::vector<Evaluator*> &uz = element->settings(Property::TRANSLATION_MOTION_Z);

	const Material* material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN(Element::ElementPointType::GAUSSE_POINT);
	const std::vector<DenseMatrix> &N = element->N(Element::ElementPointType::GAUSSE_POINT);
	const std::vector<double> &weighFactor = element->weighFactor(Element::ElementPointType::GAUSSE_POINT);

	DenseMatrix matFlux(3, 1), matGradient(3, 1);

	for (size_t i = 0; i < element->nodes(); i++) {
		const Point &p = mesh.coordinates()[element->node(i)];
		coordinates(i, 0) = p.x;
		coordinates(i, 1) = p.y;
		coordinates(i, 2) = p.z;
		U(i, 0) = ux.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));
		U(i, 1) = uy.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));
		U(i, 2) = uz.back()->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::DENSITY)->evaluate(element->node(i)) * material->get(AdvectionDiffusion3DMaterial::HEAT_CAPACITY)->evaluate(element->node(i));

		switch (material->getModel()) {
		case AdvectionDiffusion3DMaterial::ISOTROPIC:
			K(i, 0) = K(i, 1) = K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case AdvectionDiffusion3DMaterial::DIAGONAL:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = K(i, 4) = K(i, 5) = K(i, 6) = K(i, 7) = K(i, 8) = 0;
			break;
		case AdvectionDiffusion3DMaterial::SYMMETRIC:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = K(i, 5) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i));
			K(i, 4) = K(i, 7) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i));
			K(i, 6) = K(i, 8) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i));
			break;
		case AdvectionDiffusion3DMaterial::ANISOTROPIC:
			K(i, 0) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XX)->evaluate(element->node(i));
			K(i, 1) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YY)->evaluate(element->node(i));
			K(i, 2) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZZ)->evaluate(element->node(i));
			K(i, 3) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XY)->evaluate(element->node(i));
			K(i, 4) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_XZ)->evaluate(element->node(i));
			K(i, 5) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YX)->evaluate(element->node(i));
			K(i, 6) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_YZ)->evaluate(element->node(i));
			K(i, 7) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZX)->evaluate(element->node(i));
			K(i, 8) = material->get(AdvectionDiffusion3DMaterial::THERMAL_CONDUCTIVITY_ZY)->evaluate(element->node(i));
			break;
		}
	}

	DenseMatrix u(1, 3), v(1, 3), Re(1, element->nodes());

	for (size_t gp = 0; gp < element->gaussePoints(); gp++) {
		u.multiply(N[gp], U, 1, 0);

		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		gpK.multiply(N[gp], K);

		Ce(0, 0) = gpK(0, 0);
		Ce(1, 1) = gpK(0, 1);
		Ce(2, 2) = gpK(0, 2);

		Ce(0, 1) = gpK(0, 3);
		Ce(0, 2) = gpK(0, 4);
		Ce(1, 2) = gpK(0, 6);

		Ce(1, 0) = gpK(0, 5);
		Ce(2, 0) = gpK(0, 7);
		Ce(2, 1) = gpK(0, 8);

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
		Ce(3, 3) += sigma * h_e * norm_u_e;

		matGradient.multiply(dND, solution, 1, 1);
		matFlux.multiply(Ce, dND * solution, 1, 1);
	}
	gradient.push_back(matGradient(0, 0) / element->gaussePoints());
	gradient.push_back(matGradient(1, 0) / element->gaussePoints());
	gradient.push_back(matGradient(2, 0) / element->gaussePoints());
	flux.push_back(matFlux(0, 0) / element->gaussePoints());
	flux.push_back(matFlux(1, 0) / element->gaussePoints());
	flux.push_back(matFlux(2, 0) / element->gaussePoints());
}

void AdvectionDiffusion3D::postProcess(store::Store &store, const std::vector<std::vector<double> > &solution)
{
	if (!_configuration.post_process) {
		return;
	}

	std::vector<std::vector<double> > termalGradient(_mesh.parts()), termalFlux(_mesh.parts());
	DenseMatrix eSolution;

	for (size_t p = 0; p < _mesh.parts(); p++) {
		termalGradient[p].reserve(3 * matrixSize[p]);
		termalFlux[p].reserve(3 * matrixSize[p]);
		for (size_t e = _mesh.getPartition()[p]; e < _mesh.getPartition()[p + 1]; e++) {
			eSolution.resize(_mesh.elements()[e]->nodes(), 1);
			for (size_t n = 0; n < _mesh.elements()[e]->nodes(); n++) {
				eSolution(n, 0) = solution[p][_mesh.nodes()[_mesh.elements()[e]->node(n)]->DOFIndex(p, 0)];
			}
			postProcessElement(termalGradient[p], termalFlux[p], eSolution, _mesh.elements()[e], _mesh, _configuration);
		}
	}

	store.storeValues("gradient", 3, termalGradient, store::Store::ElementType::ELEMENTS);
	store.storeValues("flux", 3, termalFlux, store::Store::ElementType::ELEMENTS);
}


}

