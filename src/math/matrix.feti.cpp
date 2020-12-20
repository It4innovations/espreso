
#include "matrix.feti.h"
#include "vector.dense.feti.h"
#include "esinfo/eslog.h"

#include <cstddef>

using namespace espreso;

MatrixFETI::MatrixFETI()
: DataDecomposition(DataDecomposition::DUPLICATION::SPLIT), domains(0), matrices(NULL)
{

}

MatrixFETI::MatrixFETI(const MatrixFETI &other)
: DataDecomposition(DataDecomposition::DUPLICATION::SPLIT), domains(0), matrices(NULL)
{
	deepCopy(&other);
}

MatrixFETI& MatrixFETI::operator=(const MatrixFETI &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

MatrixFETI::~MatrixFETI()
{
	for (esint d = 0; d < domains; ++d) {
		delete at(d);
	}
	delete[] matrices;
}

Matrix* MatrixFETI::copy()
{
	eslog::internalFailure("MatrixFETI cannot be used.\n");
	return NULL;
}

Matrix* MatrixFETI::create()
{
	eslog::internalFailure("cannot create instance of MatrixFETI.\n");
	return NULL;
}

void MatrixFETI::initDomains(esint domains)
{
	for (esint d = 0; d < this->domains; ++d) {
		delete at(d);
	}
	if (matrices) { delete[] matrices; }

	this->domains = domains;
	matrices = new Matrix*[domains];
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		matrices[d] = create();
	}
}

void MatrixFETI::structureUpdated()
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->structureUpdated();
	}
}

void MatrixFETI::swap(Matrix *other)
{
	Matrix::_swap(other);
	MatrixFETI *_other = other->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->swap(_other->at(d));
	}
	DataDecomposition::swap(_other);
}

void MatrixFETI::shallowCopy(const Matrix *other)
{
	Matrix::_assign(other);
	const MatrixFETI *_other = other->downcast<MatrixFETI>();
	initDomains(_other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->shallowCopy(_other->at(d));
	}
	DataDecomposition::shallowCopy(_other);
}

void MatrixFETI::shallowCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	const MatrixFETI *_other = other->downcast<MatrixFETI>();
	initDomains(_other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->shallowCopyStructure(_other->at(d));
	}
	DataDecomposition::shallowCopyStructure(_other);
}

void MatrixFETI::deepCopy(const Matrix *other)
{
	Matrix::_assign(other);
	const MatrixFETI *_other = other->downcast<MatrixFETI>();
	initDomains(_other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->deepCopy(_other->at(d));
	}
	DataDecomposition::deepCopy(_other);
}

void MatrixFETI::deepCopyStructure(const Matrix *other)
{
	Matrix::_assign(other);
	const MatrixFETI *_other = other->downcast<MatrixFETI>();
	initDomains(_other->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->deepCopyStructure(_other->at(d));
	}
	DataDecomposition::deepCopyStructure(_other);
}

void MatrixFETI::uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond)
{
	const MatrixFETI *_first = first->downcast<MatrixFETI>();
	const MatrixFETI *_second = second->downcast<MatrixFETI>();
	initDomains(_first->domains);
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->uniformCombination(_first->at(d), _second->at(d), nfirst, nsecond);
	}
	DataDecomposition::uniformCombination(_first, _second, nfirst, nsecond);
}

void MatrixFETI::fill(double value)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fill(value);
	}
}

void MatrixFETI::fillData(const Matrix *in)
{
	const MatrixFETI *_in = in->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fillData(_in->at(d));
	}
}

void MatrixFETI::fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	const MatrixFETI *_in = in->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fillCombinedData(_in->at(d), roffset, coffset, nsize, sumsize);
	}
}

void MatrixFETI::apply(const Vector *in, Vector *out)
{
	const VectorDenseFETI *_in = in->downcast<VectorDenseFETI>();
	VectorDenseFETI *_out = out->downcast<VectorDenseFETI>();
	if (_in->duplications != _out->duplications) {
		eslog::internalFailure("various FETI vector duplications.\n");
	}
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->apply(_in->at(d), _out->at(d));
	}

	if (_out->duplications == DataDecomposition::DUPLICATION::DUPLICATE) {
		_out->sumDuplications();
	}
}

void MatrixFETI::apply(const Vectors *in, Vectors *out)
{
	for (esint n = 0; n < in->nvectors; ++n) {
		apply(in->at(n), out->at(n));
	}
}

void MatrixFETI::scale(double alpha)
{
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->scale(alpha);
	}
}

void MatrixFETI::add(double alpha, const Matrix *a)
{
	const MatrixFETI *_a = a->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->add(alpha, _a->at(d));
	}
}

void MatrixFETI::sum(double alpha, const Matrix *a, double beta, const Matrix *b)
{
	const MatrixFETI *_a = a->downcast<MatrixFETI>();
	const MatrixFETI *_b = b->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->sum(alpha, _a->at(d), beta, _b);
	}
}

void MatrixFETI::addToCombination(double alpha, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize)
{
	const MatrixFETI *_in = in->downcast<MatrixFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->addToCombination(alpha, _in->at(d), roffset, coffset, nsize, sumsize);
	}
}

void MatrixFETI::fillDiagonal(Vector *diagonal) const
{
	VectorDenseFETI *_diaonal = diagonal->downcast<VectorDenseFETI>();
	#pragma omp parallel for
	for (esint d = 0; d < domains; ++d) {
		at(d)->fillDiagonal(_diaonal->at(d));
	}
}

double MatrixFETI::norm()
{
	eslog::internalFailure("call empty function.\n");
	return 0;
}

