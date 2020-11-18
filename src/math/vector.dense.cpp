
#include "vector.dense.h"
#include "vector.dense.distributed.h"
#include "vector.dense.feti.h"
#include "vector.sparse.h"
#include "matrix.dense.h"
#include "math.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <cstddef>
#include <cstring>

using namespace espreso;

VectorDense::VectorDense()
{

}

VectorDense::VectorDense(esint nvals)
: DataVectorDense(nvals)
{

}

VectorDense::VectorDense(esint nvals, double *vals)
: DataVectorDense(nvals, vals)
{

}

VectorDense::VectorDense(const VectorDense &other)
{
	deepCopy(&other);
}

VectorDense::VectorDense(VectorDense &&other)
{
	shallowCopy(&other);
	std::swap(_allocated, other._allocated);
}

VectorDense& VectorDense::operator=(const VectorDense &other)
{
	if (this != &other) {
		deepCopy(&other);
	}
	return *this;
}

VectorDense::~VectorDense()
{

}

VectorDense* VectorDense::shallowCopy()
{
	VectorDense *copy = this->copy();
	copy->shallowCopy(this);
	return copy;
}

VectorDense* VectorDense::shallowCopyStructure()
{
	VectorDense *copy = this->copy();
	copy->shallowCopyStructure(this);
	return copy;
}

VectorDense* VectorDense::deepCopy()
{
	VectorDense *copy = this->copy();
	copy->deepCopy(this);
	return copy;
}

VectorDense* VectorDense::deepCopyStructure()
{
	VectorDense *copy = this->copy();
	copy->deepCopyStructure(this);
	return copy;
}

VectorDense* VectorDense::copy()
{
	return new VectorDense();
}

void VectorDense::swap(Vector *other)
{
	DataVectorDense::swap(other->downcast<VectorDense>());
}

void VectorDense::shallowCopy(const Vector *other)
{
	DataVectorDense::shallowCopy(other->downcast<VectorDense>());
}

void VectorDense::shallowCopyStructure(const Vector *other)
{
	DataVectorDense::shallowCopyStructure(other->downcast<VectorDense>());
}

void VectorDense::shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors)
{
	const VectorDense *_other = other->downcast<VectorDense>();
	DataVectorDense::shallowCopyFromHolder(_other, offset, nvectors);
}

void VectorDense::deepCopy(const Vector *other)
{
	DataVectorDense::deepCopy(other->downcast<VectorDense>());
}

void VectorDense::deepCopyStructure(const Vector *other)
{
	DataVectorDense::deepCopyStructure(other->downcast<VectorDense>());
}

void VectorDense::uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond)
{
	DataVectorDense::uniformCombination(first->downcast<VectorDense>(), second->downcast<VectorDense>());
}

void VectorDense::fill(double value)
{
	DataVectorDense::fill(value);
}

void VectorDense::fillData(const Vector *in)
{
	if (dynamic_cast<const VectorDense*>(in)) {
		DataVectorDense::fillValues(dynamic_cast<const VectorDense*>(in)->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseDistributed*>(in)) {
		DataVectorDense::fillValues(dynamic_cast<const VectorDenseDistributed*>(in)->vals);
		return;
	}
	if (dynamic_cast<const VectorSparse*>(in)) {
		const VectorSparse * _in = dynamic_cast<const VectorSparse*>(in);
		DataVectorDense::fillSparseValues(_in->nnz, _in->indices, _in->vals);
		return;
	}
	if (dynamic_cast<const VectorDenseFETI*>(in)) {
		const VectorDenseFETI* _in = in->downcast<VectorDenseFETI>();
		esint i = 0;
		for (auto n = _in->dmap->begin(); n != _in->dmap->end(); ++n, ++i) {
			for (auto di = n->begin(); di != n->end(); ++di) {
				if (_in->ismy(di->domain)) {
					vals[i] = _in->at(di->domain - _in->distribution[info::mpi::rank])->vals[di->index];
					break;
				}
			}
		}
		return;
	}
	eslog::error("ESPRESO internal error: unsupported math operation: %s->fillData(%s).\n", name(), in->name());
}

void VectorDense::fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	DataVectorDense::fillCombinedValues(in->downcast<VectorDense>(), offset, nsize, sumsize);
}

void VectorDense::fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize)
{
	DataVectorDense::fillValuesFromCombination(in->downcast<VectorDense>(), offset, nsize, sumsize);
}

void VectorDense::scale(double alpha)
{
	MATH::vecScale(size, alpha, vals);
}

void VectorDense::add(double alpha, const Vector *a)
{
	if (dynamic_cast<const VectorDense*>(a)) {
		MATH::vecAdd(size, vals, alpha, dynamic_cast<const VectorDense*>(a)->vals);
		return;
	}
	if (dynamic_cast<const VectorSparse*>(a)) {
		const VectorSparse * _a = dynamic_cast<const VectorSparse*>(a);
		MATH::vecAddSparse(_a->nnz, vals, alpha, _a->indices, _a->vals);
		return;
	}
	eslog::error("ESPRESO internal error: unsupported math operation.\n");
}

void VectorDense::sum(double alpha, const Vector *a, double beta, const Vector *b)
{
	const VectorDense *_a = a->downcast<VectorDense>();
	const VectorDense *_b = b->downcast<VectorDense>();
	MATH::vecSum(size, vals, alpha, _a->vals, beta, _b->vals);
}

void VectorDense::addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize)
{
	DataVectorDense::addToCombination(alpha, in->downcast<VectorDense>(), offset, nsize, sumsize);
}

//void VectorDense::addToCombination(double alpha, VectorSparse *in, esint offset, esint nsize, esint sumsize)
//{
//	for (esint i = 0; i < in->nvals / nsize; ++i) {
//		for (esint j = 0; j < nsize; ++j) {
//			vals[sumsize * (in->indices[i * nsize + j] / nsize) + offset + in->indices[i * nsize + j] % nsize] += alpha * in->vals[i * nsize + j];
//		}
//	}
//}

double VectorDense::norm()
{
	return MATH::vecNorm(size, vals);
}

double VectorDense::max()
{
	double max = size ? vals[0] : 0;
	for (esint i = 1; i < size; ++i) {
		if (max < vals[i]) {
			max = vals[i];
		}
	}
	return max;
}

double VectorDense::absmax()
{
	double max = size ? std::fabs(vals[0]) : 0;
	for (esint i = 1; i < size; ++i) {
		if (max < std::fabs(vals[i])) {
			max = std::fabs(vals[i]);
		}
	}
	return max;
}

double VectorDense::dot(const Vector *other)
{
	return MATH::vecDot(size, vals, other->downcast<VectorDense>()->vals);
}

void VectorDense::toFETI(VectorFETI *other) const
{
	if (dynamic_cast<VectorDenseFETI*>(other)) {
		VectorDenseFETI *_other = dynamic_cast<VectorDenseFETI*>(other);
		esint i = 0;
		for (auto map = _other->dmap->begin(); map != _other->dmap->end(); ++map, ++i) {
			switch (_other->duplications) {
			case DataDecomposition::DUPLICATION::DUPLICATE: {
				for (auto di = map->begin(); di != map->end(); ++di) {
					if (_other->ismy(di->domain)) {
						_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i];
					}
				}
			} break;
			case DataDecomposition::DUPLICATION::SPLIT: {
				for (auto di = map->begin(); di != map->end(); ++di) {
					if (_other->ismy(di->domain)) {
						_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i] / map->size();
					}
				}
			} break;
			case DataDecomposition::DUPLICATION::SPLIT_DOMAINS: {
				int size = 0;
				for (auto di = map->begin(); di != map->end(); ++di) {
					if (_other->ismy(di->domain)) {
						++size;
					}
				}
				for (auto di = map->begin(); di != map->end(); ++di) {
					if (_other->ismy(di->domain)) {
						_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i] / size;
					}
				}
			} break;
			}
		}
	}
}

void VectorDense::toCombinedFETI(VectorFETI *other, esint offset, esint nsize, esint sumsize) const
{
	if (dynamic_cast<VectorDenseFETI*>(other)) {
		VectorDenseFETI *_other = dynamic_cast<VectorDenseFETI*>(other);
		esint i = 0;
		for (auto map = _other->dmap->begin(); map != _other->dmap->end(); ++map, ++i) {
			if (offset <= i % sumsize && i % sumsize < offset + nsize) {
				switch (_other->duplications) {
				case DataDecomposition::DUPLICATION::DUPLICATE: {
					for (auto di = map->begin(); di != map->end(); ++di) {
						if (_other->ismy(di->domain)) {
							_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i / sumsize + i % nsize + offset];
						}
					}
				} break;
				case DataDecomposition::DUPLICATION::SPLIT: {
					for (auto di = map->begin(); di != map->end(); ++di) {
						if (_other->ismy(di->domain)) {
							_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i / sumsize + i % nsize + offset] / map->size();
						}
					}
				} break;
				case DataDecomposition::DUPLICATION::SPLIT_DOMAINS: {
					int size = 0;
					for (auto di = map->begin(); di != map->end(); ++di) {
						if (_other->ismy(di->domain)) {
							++size;
						}
					}
					for (auto di = map->begin(); di != map->end(); ++di) {
						if (_other->ismy(di->domain)) {
							_other->at(di->domain - _other->doffset)->vals[di->index] = vals[i / sumsize + i % nsize + offset] / size;
						}
					}
				} break;
				}
			}
		}
	}
}

// this = alpha A * B + beta * this;
void VectorDense::multiply(
		const MatrixDense &A, const MatrixDense &B,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{

//	resize(transposeB ? B.nrows : B.ncols);
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, B.nrows, B.ncols, B.vals,
			beta, vals);
}

void VectorDense::multiply(
		const MatrixDense &A, const VectorDense &B,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, B.size, 1, B.vals,
			beta, vals);
}

VectorsDense::VectorsDense()
{
	initVectors(0);
}

VectorsDense::~VectorsDense()
{

}

Vector* VectorsDense::create()
{
	return new VectorDense();
}

VectorsDense* VectorsDense::copy()
{
	return new VectorsDense();
}

void VectorsDense::resize(esint size)
{
	holder()->resize(nvectors * size);
	for (esint n = 0; n < nvectors; ++n) {
		at(n)->shallowCopyFromHolder(holder(), n, nvectors);
	}
}

// this = alpha A * B + beta * this;
void VectorsDense::multiply(
		const MatrixDense &A, const MatrixDense &B,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
//	resize(transposeB ? B.nrows : B.ncols, transposeA ? A.ncols : A.nrows);
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, B.nrows, B.ncols, B.vals,
			beta,
			reinterpret_cast<VectorDense*>(_holder)->vals);
}

void VectorsDense::multiply(
		const MatrixDense &A,
		esint brows, esint bcols, double* bvals,
		double alpha, double beta,
		bool transposeA, bool transposeB)
{
	MATH::DenseMatDenseMatRowMajorProduct(
			alpha, transposeA, A.nrows, A.ncols, A.vals,
			transposeB, brows, bcols, bvals,
			beta,
			reinterpret_cast<VectorDense*>(_holder)->vals);
}
