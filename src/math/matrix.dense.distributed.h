
#ifndef SRC_WRAPPERS_MATH_MATRIXDENSEDISTRIBUTED_H_
#define SRC_WRAPPERS_MATH_MATRIXDENSEDISTRIBUTED_H_

#include "data.distributed.h"
#include "matrix.h"
#include "data.matrix.dense.h"

namespace espreso {

class MatrixCSR;
class VectorDense;
struct MatrixCSRDistributedData;

class MatrixDenseDistributed: public Matrix, public DataMatrixDense, public DataDistributed
{
public:
	class MatrixDenseRow
	{
	public:
		MatrixDenseRow(double *vals): vals(vals) {}

		double& operator[](esint col) { return vals[col]; }
		const double& operator[](esint col) const { return vals[col]; }

		double *vals;
	};

	MatrixDenseRow operator[](esint row) { return MatrixDenseRow(vals + ncols * row); }
	const MatrixDenseRow operator[](esint row) const { return MatrixDenseRow(vals + ncols * row); }

	double& operator()(esint row, esint col) { return vals[ncols * row + col]; }
	const double& operator()(esint row, esint col) const { return vals[ncols * row + col]; }

	double& at(esint row, esint col) { return vals[ncols * row + col]; }
	const double& at(esint row, esint col) const { return vals[ncols * row + col]; }

	MatrixDenseDistributed();
	MatrixDenseDistributed(esint nrows, esint ncols, esint nhalo, esint nneighbors);
	MatrixDenseDistributed(const MatrixDenseDistributed &other);
	MatrixDenseDistributed& operator=(const MatrixDenseDistributed &other);
	~MatrixDenseDistributed();

	MatrixDenseDistributed* copy();

	void resize(esint nrows, esint ncols, esint nhalo, esint nneighbors);
	void structureUpdated();
	void swap(Matrix *other);
	void shallowCopy(const Matrix *other);
	void shallowCopyStructure(const Matrix *other);
	void deepCopy(const Matrix *other);
	void deepCopyStructure(const Matrix *other);
	void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond);

	void fill(double value);
	void fillData(const Matrix *in);
	void fillCombinedData(const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum);

	void apply(const Vector *in, Vector *out);
	void apply(const Vectors *in, Vectors *out);

	void scale(double alpha);
	void add(double alpha, const Matrix *a);
	void sum(double alpha, const Matrix *a, double beta, const Matrix *b);
	void addToCombination(double scale, const Matrix *in,  esint roffset, esint coffset, esint rsize, esint csize, esint rsum, esint csum);

	void fillDiagonal(Vector *diagonal) const;

	double norm();

	const char* name() const { return "MatrixDenseDistributed"; }
};

}

#endif /* SRC_WRAPPERS_MATH_MATRIXDENSEDISTRIBUTED_H_ */
