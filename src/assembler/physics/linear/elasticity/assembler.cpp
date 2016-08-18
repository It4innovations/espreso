
#include "assembler.h"

#include "../../../constraints/equalitygluing.h"

using namespace espreso;

std::vector<Property> LinearElasticity::elementDOFs;
std::vector<Property> LinearElasticity::faceDOFs;
std::vector<Property> LinearElasticity::edgeDOFs;
std::vector<Property> LinearElasticity::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> LinearElasticity::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

void LinearElasticity::prepareMeshStructures()
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

	if (config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS) {
		_mesh.computeFixPoints(config::mesh::FIX_POINTS);
	}

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_mesh.computeCorners(config::mesh::CORNERS, config::mesh::VERTEX_CORNERS, config::mesh::EDGE_CORNERS, config::mesh::FACE_CORNERS);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			_mesh.computeFacesSharedByDomains();
			break;
		}
	}
}


static void composeFacesGluing(Mesh &mesh, Constraints &constrains)
{
	std::vector<Element*> faces;
	for (size_t i = 0; i < mesh.faces().size(); i++) {
		if (mesh.faces()[i]->domains().size() == 2) {
			faces.push_back(mesh.faces()[i]);
		}
	}

	std::sort(faces.begin(), faces.end(), [] (Element* e1, Element* e2) { return e1->domains() < e2->domains(); });

	std::vector<size_t> distribution;
	distribution.push_back(0);
	for (size_t i = 1; i < faces.size(); i++) {
		if (faces[i]->domains() != faces[i - 1]->domains()) {
			distribution.push_back(i);
		}
	}
	distribution.push_back(faces.size());

	std::vector<std::vector<std::vector<eslocal> > > rows(mesh.parts(), std::vector<std::vector<eslocal> >(distribution.size() - 1));
	std::vector<std::vector<std::vector<eslocal> > > cols(mesh.parts(), std::vector<std::vector<eslocal> >(distribution.size() - 1));
	std::vector<std::vector<std::vector<double> > > vals(mesh.parts(), std::vector<std::vector<double> >(distribution.size() - 1));

	cilk_for (size_t interface = 0; interface < distribution.size() - 1; interface++) {
		std::vector<eslocal> nodes;
		for (size_t f = distribution[interface]; f < distribution[interface + 1]; f++) {
			for (size_t n = 0; n < faces[f]->nodes(); n++) {
				nodes.push_back(faces[f]->node(n));
			}
			std::sort(nodes.begin(), nodes.end());
			Esutils::removeDuplicity(nodes);
		}

		Point center(0, 0 ,0);
		for (size_t n = 0; n < nodes.size(); n++) {
			center += mesh.coordinates()[nodes[n]];
		}
		center /= nodes.size();

		std::vector<eslocal> domains = { faces[distribution[interface]]->domains()[0], faces[distribution[interface]]->domains()[1] };
		for (size_t d = 0; d < domains.size(); d++) {
			eslocal domain = domains[d];
			eslocal sign = d ? -1 : 1;

			for (size_t r = 0; r < 3; r++) {
				rows[domain][interface].insert(rows[domain][interface].end(), nodes.size(), interface * 6 + r + IJVMatrixIndexing);
				vals[domain][interface].insert(vals[domain][interface].end(), nodes.size(), sign);
				for (size_t n = 0; n < nodes.size(); n++) {
					cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, r) + IJVMatrixIndexing);
				}
			}

			rows[domain][interface].insert(rows[domain][interface].end(), 2 * nodes.size(), interface * 6 + 3 + IJVMatrixIndexing);
			for (size_t n = 0; n < nodes.size(); n++) {
				Point p = mesh.coordinates()[nodes[n]];
				p -= center;
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 0) + IJVMatrixIndexing);
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 1) + IJVMatrixIndexing);
				vals[domain][interface].push_back(-sign * p.y);
				vals[domain][interface].push_back( sign * p.x);
			}

			rows[domain][interface].insert(rows[domain][interface].end(), 2 * nodes.size(), interface * 6 + 4 + IJVMatrixIndexing);
			for (size_t n = 0; n < nodes.size(); n++) {
				Point p = mesh.coordinates()[nodes[n]];
				p -= center;
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 0) + IJVMatrixIndexing);
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 2) + IJVMatrixIndexing);
				vals[domain][interface].push_back(-sign * p.z);
				vals[domain][interface].push_back( sign * p.x);
			}

			rows[domain][interface].insert(rows[domain][interface].end(), 2 * nodes.size(), interface * 6 + 5 + IJVMatrixIndexing);
			for (size_t n = 0; n < nodes.size(); n++) {
				Point p = mesh.coordinates()[nodes[n]];
				p -= center;
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 1) + IJVMatrixIndexing);
				cols[domain][interface].push_back(mesh.nodes()[nodes[n]]->DOFIndex(domain, 2) + IJVMatrixIndexing);
				vals[domain][interface].push_back(-sign * p.z);
				vals[domain][interface].push_back( sign * p.y);
			}
		}
	}

	for (size_t p = 0; p < mesh.parts(); p++) {
		for (size_t interface = 0; interface < distribution.size() - 1; interface++) {
			constrains.B0[p].I_row_indices.insert(constrains.B0[p].I_row_indices.end(), rows[p][interface].begin(), rows[p][interface].end());
			constrains.B0[p].J_col_indices.insert(constrains.B0[p].J_col_indices.end(), cols[p][interface].begin(), cols[p][interface].end());
			constrains.B0[p].V_values     .insert(constrains.B0[p].V_values.end()     , vals[p][interface].begin(), vals[p][interface].end());
		}
		for (size_t i = 0; i < constrains.B0[p].I_row_indices.size(); i++) {
			constrains.B0subdomainsMap[p].push_back(constrains.B0[p].I_row_indices[i] - IJVMatrixIndexing);
		}
		constrains.B0[p].nnz = constrains.B0[p].I_row_indices.size();
		constrains.B0[p].rows = 6 * (distribution.size() - 1);
		constrains.B0[p].cols = constrains.B1[p].cols;
	}
}

void LinearElasticity::assembleGluingMatrices()
{
	_constraints.initMatrices(matrixSize);

	_constraints.insertDirichletToB1(_mesh.nodes(), _mesh.coordinates(), pointDOFs);
	_constraints.insertElementGluingToB1(_mesh.nodes(), pointDOFs);

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		switch (config::solver::B0_TYPE) {
		case config::solver::B0_TYPEalternative::CORNERS:
			_constraints.insertDomainGluingToB0(_mesh.corners(), pointDOFs);
			break;
		case config::solver::B0_TYPEalternative::KERNELS:
			composeFacesGluing(_mesh, _constraints);
			break;
		}
	}
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

// B =
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
static void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();
	const double *dNDz = dND.values() + 2 * dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns + 2 * dND.columns()], dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + 2 * dND.columns()], dNDy, sizeof(double) * dND.columns());

	memcpy(&v[2 * columns + 2 * dND.columns()], dNDz, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + dND.columns()],     dNDz, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns],                     dNDz, sizeof(double) * dND.columns());
}

static void processElement(DenseMatrix &Ke, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	DenseMatrix Ce(6, 6), coordinates, J, invJ, dND, B;
	std::vector<double> inertia(3, 0);
	double detJ;

	const Material &material = mesh.materials()[element->param(Element::MATERIAL)];
	const std::vector<DenseMatrix> &dN = element->dN();
	const std::vector<DenseMatrix> &N = element->N();
	const std::vector<double> &weighFactor = element->weighFactor();

	// TODO: set the omega from example
	Point omega(50, 50, 0);

	double ex = material.youngModulus;
	double mi = material.poissonRatio;
	double E = ex / ((1 + mi) * (1 - 2 * mi));
	Ce(0, 1) = Ce(0, 2) = Ce(1, 0) = Ce(1, 2) = Ce(2, 0) = Ce(2, 1) = E * mi;
	Ce(0, 0) = Ce(1, 1) = Ce(2, 2) = E * (1.0 - mi);
	Ce(3, 3) = Ce(4, 4) = Ce(5, 5) = E * (0.5 - mi);

	inertia[0] = inertia[1] = 0;
	inertia[2] = 9.8066 * material.density;

	coordinates.resize(element->nodes(), 3);

	Point mid;
	for (size_t i = 0; i < element->nodes(); i++) {
		coordinates(i, 0) = mesh.coordinates()[element->node(i)].x;
		coordinates(i, 1) = mesh.coordinates()[element->node(i)].y;
		coordinates(i, 2) = mesh.coordinates()[element->node(i)].z;
		mid += mesh.coordinates()[element->node(i)];
	}
	mid /= element->nodes();

	eslocal Ksize = 3 * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);

	double rotation[3] = { mid.x * omega.x * omega.x, mid.y * omega.y * omega.y, mid.z * omega.z * omega.z };

	for (eslocal gp = 0; gp < element->gaussePoints(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);
		B.resize(Ce.rows(), Ksize);
		distribute(B, dND);
		Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			// TODO: set rotation from example
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
			//fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * 7850 * rotation[i / e->size()];
		}
	}
}

static void analyticsKernels(SparseMatrix &R1, const Coordinates &coordinates, size_t subdomain)
{
	size_t nodes = coordinates.localSize(subdomain);
	R1.rows = 3 * nodes;
	R1.cols = 6;
	R1.nnz = R1.rows * R1.cols;
	R1.type = 'G';

	R1.dense_values.reserve(R1.nnz);

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };
		kernel[c] = 1;
		for (size_t i = 0; i < nodes; i++) {
			R1.dense_values.insert(R1.dense_values.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.y);
		R1.dense_values.push_back( p.x);
		R1.dense_values.push_back(   0);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back( p.x);
	}

	for (size_t i = 0; i < coordinates.localSize(subdomain); i++) {
		const Point &p = coordinates.get(i, subdomain);
		R1.dense_values.push_back(   0);
		R1.dense_values.push_back(-p.z);
		R1.dense_values.push_back( p.y);
	}
}

static void analyticsRegMat(SparseMatrix &K, SparseMatrix &RegMat, const std::vector<Element*> &fixPoints, const Coordinates &coordinates, size_t subdomain)
{
	ESTEST(MANDATORY) << "Too few FIX POINTS: " << fixPoints.size() << (fixPoints.size() > 3 ? TEST_PASSED : TEST_FAILED);

	SparseMatrix Nt; // CSR matice s DOFY
	Nt.rows = 6;
	Nt.cols = K.cols;
	Nt.nnz  = 9 * fixPoints.size();
	Nt.type = 'G';

	std::vector<eslocal> &ROWS = Nt.CSR_I_row_indices;
	std::vector<eslocal> &COLS = Nt.CSR_J_col_indices;
	std::vector<double>  &VALS = Nt.CSR_V_values;

	ROWS.reserve(Nt.rows + 1);
	COLS.reserve(Nt.nnz);
	VALS.reserve(Nt.nnz);

	ROWS.push_back(1);
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());
	ROWS.push_back(ROWS.back() + 2 * fixPoints.size());

	for (size_t c = 0; c < 3; c++) {
		std::vector<double> kernel = { 0, 0, 0 };

		kernel[c] = 1;
		for (size_t i = 0; i < fixPoints.size(); i++) {
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
			COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + IJVMatrixIndexing);
			VALS.insert(VALS.end(), kernel.begin(), kernel.end());
		}
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
		VALS.push_back(-p.y);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 0) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + IJVMatrixIndexing);
		VALS.push_back(-p.z);
		VALS.push_back( p.x);
	}

	for (size_t i = 0; i < fixPoints.size(); i++) {
		const Point &p = coordinates[fixPoints[i]->node(0)];
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 1) + IJVMatrixIndexing);
		COLS.push_back(fixPoints[i]->DOFIndex(subdomain, 2) + IJVMatrixIndexing);
		VALS.push_back(-p.z);
		VALS.push_back( p.y);
	}

	SparseMatrix N;
	Nt.MatTranspose( N );

	RegMat.MatMat(Nt, 'N', N);
	RegMat.MatTranspose();
	RegMat.RemoveLower();

	SparseSolverCPU NtN;
	NtN.ImportMatrix(RegMat);
	RegMat.Clear();

	NtN.Factorization("Create RegMat");
	NtN.SolveMat_Sparse(Nt);
	NtN.Clear();

	RegMat.MatMat(N, 'N', Nt);
	RegMat.MatScale(K.getDiagonalMaximum());
}

static void algebraicKernelsAndRegularization(SparseMatrix &K, SparseMatrix &RegMat, SparseMatrix &R, size_t subdomain)
{
	double norm;
	eslocal defect;

	K.get_kernel_from_K(K, RegMat, R, norm, defect, subdomain);
}

void LinearElasticity::composeSubdomain(size_t subdomain)
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

		processElement(Ke, fe, _mesh, subdomain, elements[e]);

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

	switch (config::solver::REGULARIZATION) {
	case config::solver::REGULARIZATIONalternative::FIX_POINTS:
		if (config::assembler::DOFS_ORDER == config::assembler::DOFS_ORDERalternative::GROUP_DOFS) {
			ESINFO(GLOBAL_ERROR) << "Implement regularization for GROUP_DOFS alternative";
		}

		analyticsKernels(R1[subdomain], _mesh.coordinates(), subdomain);
		analyticsRegMat(K[subdomain], RegMat[subdomain], _mesh.fixPoints(subdomain), _mesh.coordinates(), subdomain);
		K[subdomain].RemoveLower();
		RegMat[subdomain].RemoveLower();
		K[subdomain].MatAddInPlace(RegMat[subdomain], 'N', 1);
		RegMat[subdomain].ConvertToCOO(1);
		break;
	case config::solver::REGULARIZATIONalternative::NULL_PIVOTS:
		K[subdomain].RemoveLower();
		algebraicKernelsAndRegularization(K[subdomain], RegMat[subdomain], R1[subdomain], subdomain);
		break;
	}
}




