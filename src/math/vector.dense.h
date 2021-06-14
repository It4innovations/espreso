
#ifndef SRC_WRAPPERS_MATH_VECTORDENSE_H_
#define SRC_WRAPPERS_MATH_VECTORDENSE_H_

#include "vector.h"
#include "data.vector.dense.h"

#include <cstring>

namespace espreso {

class MatrixDense;
class VectorSparse;
class VectorsSparse;
class VectorFETI;

class VectorDense: public Vector, public DataVectorDense
{
public:
	double& operator[](esint index) { return vals[index]; }
	const double& operator[](esint index) const { return vals[index]; }

	VectorDense();
	VectorDense(esint nvals);
	VectorDense(esint nvals, double *vals);

	VectorDense(const VectorDense &other);
	VectorDense(VectorDense &&other);
	VectorDense& operator=(const VectorDense &other);

	~VectorDense();

	VectorDense* shallowCopy();
	VectorDense* shallowCopyStructure();
	VectorDense* deepCopy();
	VectorDense* deepCopyStructure();

	VectorDense* copy();

	void swap(Vector *other);
	void shallowCopy(const Vector *other);
	void shallowCopyStructure(const Vector *other);
	void shallowCopyFromHolder(const Vector *other, esint offset, esint nvectors);
	void deepCopy(const Vector *other);
	void deepCopyStructure(const Vector *other);
	void uniformCombination(const Vector *first, const Vector *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Vector *in);
	void fillCombinedValues(const Vector *in, esint offset, esint nsize, esint sumsize);
	void fillValuesFromCombination(const Vector *in, esint offset, esint nsize, esint sumsize);

	void scale(double alpha);
	void add(double alpha, const Vector *a);
	void sum(double alpha, const Vector *a, double beta, const Vector *b);
	void addToCombination(double alpha, const Vector *in, esint offset, esint nsize, esint sumsize);

	double norm();
	double max();
	double absmax();
	double dot(const Vector *other);

	void toFETI(VectorFETI *other) const;
	void toCombinedFETI(VectorFETI *other, esint offset, esint nsize, esint sumsize) const;

	void fromFETI(VectorFETI *other) const;

	void multiply(
			const MatrixDense &A, const MatrixDense &B,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void multiply(
			const MatrixDense &A, const VectorDense &B,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	const char* name() const { return "VectorDense"; }
};

class VectorsDense: public Vectors
{
public:
	VectorDense* holder() { return reinterpret_cast<VectorDense*>(_holder); }
	const VectorDense* holder() const { return reinterpret_cast<VectorDense*>(_holder); }

	VectorDense& operator[](esint index) { return *reinterpret_cast<VectorDense*>(_vectors[index]); }
	const VectorDense& operator[](esint index) const { return *reinterpret_cast<VectorDense*>(_vectors[index]); }

	VectorDense* at(esint index) { return reinterpret_cast<VectorDense*>(_vectors[index]); }
	const VectorDense* at(esint index) const { return reinterpret_cast<VectorDense*>(_vectors[index]); }

	VectorsDense();
	~VectorsDense();

	VectorsDense* copy();

	void resize(esint size);

	// this = alpha A * B + beta * this;
	void multiply(
			const MatrixDense &A, const MatrixDense &B,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	void multiply(
			const MatrixDense &A,
			esint brows, esint bcols, double* bvals,
			double alpha = 1, double beta = 0,
			bool transposeA = false, bool transposeB = false);

	const char* name() const { return "VectorsDense"; }
protected:
	Vector* create();
};

namespace utils {

inline size_t packedSize(const VectorDense &data)
{
	return sizeof(data.size) + sizeof(double) * data.size;
}

inline void pack(VectorDense &data, char* &p)
{
	std::memcpy(p, &data.size, sizeof(data.size));
	p += sizeof(data.size);
	std::memcpy(p, data.vals, sizeof(double) * data.size);
	p += sizeof(double) * data.size;
}

inline void unpack(VectorDense &data, const char* &p)
{
	std::memcpy(&data.size, p, sizeof(data.size));
	p += packedSize(data.size);
	data.resize(data.size);
	std::memcpy(data.vals, p, sizeof(double) * data.size);
	p += sizeof(double) * data.size;
}

} // namespace utils

}



#endif /* SRC_WRAPPERS_MATH_VECTORDENSE_H_ */
