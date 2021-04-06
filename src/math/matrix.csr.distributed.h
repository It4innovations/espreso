
#ifndef SRC_WRAPPERS_MATH_MATRIXCSRDISTRIBUTED_H_
#define SRC_WRAPPERS_MATH_MATRIXCSRDISTRIBUTED_H_

#include "data.distributed.h"
#include "matrix.h"
#include "data.matrix.csr.h"

namespace espreso {

class MatrixCSR;
class VectorDense;
struct MatrixCSRDistributedData;

class MatrixCSRDistributed: public Matrix, public DataMatrixCSR, public DataDistributed
{
public:
	MatrixCSRDistributed();
	MatrixCSRDistributed(esint nrows, esint ncols, esint nnz, esint nhalo, esint nneighbors);
	MatrixCSRDistributed(const MatrixCSRDistributed &other);
	MatrixCSRDistributed& operator=(const MatrixCSRDistributed &other);
	~MatrixCSRDistributed();

	MatrixCSRDistributed* copy();

	void resize(esint nrows, esint ncols, esint nnz, esint nhalo, esint nneighbors);
	void structureUpdated();
	void swap(Matrix *other);
	void shallowCopy(const Matrix *other);
	void shallowCopyStructure(const Matrix *other);
	void deepCopy(const Matrix *other);
	void deepCopyStructure(const Matrix *other);
	void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Matrix *in);
	void fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void apply(const Vector *in, Vector *out);
	void apply(const Vectors *in, Vectors *out);

	void scale(double alpha);
	void add(double alpha, const Matrix *a);
	void sum(double alpha, const Matrix *a, double beta, const Matrix *b);
	void addToCombination(double scale, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void fillDiagonal(Vector *diagonal) const;

	double norm();

	void gatherFromUpper();

	const char* name() const { return "MatrixCSRDistributed"; }
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXCSRDISTRIBUTED_H_ */
