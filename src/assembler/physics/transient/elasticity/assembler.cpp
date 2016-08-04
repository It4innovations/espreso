
#include "assembler.h"

using namespace espreso;

std::vector<Property> TransientElasticity::elementDOFs;
std::vector<Property> TransientElasticity::faceDOFs;
std::vector<Property> TransientElasticity::edgeDOFs;
std::vector<Property> TransientElasticity::pointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };
std::vector<Property> TransientElasticity::midPointDOFs = { Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z };

void TransientElasticity::init()
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

static void processElement(DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe, const espreso::Mesh &mesh, size_t subdomain, const Element* element)
{
	size_t DOFs = 3;
	DenseMatrix Ce(6, 6), coordinates, J, invJ, dND, B;
	std::vector<double> inertia(3, 0);
	double detJ, CP;

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
	CP = 1;

	coordinates.resize(element->nodes(), DOFs);

	Point mid;
	for (size_t i = 0; i < element->nodes(); i++) {
		coordinates(i, 0) = mesh.coordinates().get(element->node(i), subdomain).x;
		coordinates(i, 1) = mesh.coordinates().get(element->node(i), subdomain).y;
		coordinates(i, 2) = mesh.coordinates().get(element->node(i), subdomain).z;
		mid += mesh.coordinates().get(element->node(i), subdomain);
	}
	mid /= element->nodes();

	eslocal Ksize = DOFs * element->nodes();

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	Me.resize(element->nodes(), element->nodes());
	Me = 0;
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
		Me.multiply(N[gp], N[gp], material.density * detJ * weighFactor[gp] * CP, 1, true);

		for (eslocal i = 0; i < Ksize; i++) {
			// TODO: set rotation from example
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % element->nodes()) * inertia[i / element->nodes()];
			//fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * 7850 * rotation[i / e->size()];
		}
	}
}

void TransientElasticity::composeSubdomain(size_t subdomain)
{
	eslocal subdomainSize = DOFs.size() * _mesh.coordinates().localSize(subdomain);
	DenseMatrix Ke, Me;
	std::vector<double> fe;

	SparseVVPMatrix<eslocal> _K;
	SparseVVPMatrix<eslocal> _M;

	_K.resize(subdomainSize, subdomainSize);
	_M.resize(subdomainSize, subdomainSize);
	f[subdomain].resize(subdomainSize);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<Element*> &elements = _mesh.elements();

	for (eslocal i = partition[subdomain]; i < partition[subdomain + 1]; i++) {

		const Element* element = elements[i];
		processElement(Ke, Me, fe, _mesh, subdomain, element);

		for (size_t i = 0; i < DOFs.size() * element->nodes(); i++) {
			size_t row = DOFs.size() * (element->node(i % element->nodes())) + i / element->nodes();
			for (size_t j = 0; j < DOFs.size() * element->nodes(); j++) {
				size_t column = DOFs.size() * (element->node(j % element->nodes())) + j / element->nodes();
				_K(row, column) = Ke(i, j);
			}
			f[subdomain][row] += fe[i];
		}

		for (size_t r = 0; r < element->nodes(); r++) {
			for (size_t c = 0; c < element->nodes(); c++) {
				for (size_t k = 0; k < DOFs.size(); k++) {
					_M(DOFs.size() * element->node(r) + k, DOFs.size() * element->node(c) + k) = Me(r, c);
				}
			}
		}
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	SparseCSRMatrix<eslocal> csrM = _M;
	K[subdomain] = csrK;
	M[subdomain] = csrM;

	K[subdomain].MatAddInPlace(M[subdomain], 'N', timeConstant);
}




