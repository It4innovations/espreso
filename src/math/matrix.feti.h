
#ifndef SRC_WRAPPERS_MATH_MATRIXFETI_H_
#define SRC_WRAPPERS_MATH_MATRIXFETI_H_

#include "matrix.h"
#include "data.decomposition.h"

namespace espreso {

class MatrixFETI: public Matrix, public DataDecomposition
{
public:
	Matrix& operator[](esint domain) { return *matrices[domain]; }
	const Matrix& operator[](esint domain) const { return *matrices[domain]; }

	Matrix* at(esint domain) { return matrices[domain]; }
	const Matrix* at(esint domain) const { return matrices[domain]; }

	MatrixFETI();
	MatrixFETI(const MatrixFETI &other);
	MatrixFETI& operator=(const MatrixFETI &other);
	virtual ~MatrixFETI();

	virtual Matrix* copy();

	virtual void initDomains(esint domains);

	virtual void structureUpdated();
	virtual void swap(Matrix *other);
	virtual void shallowCopy(const Matrix *other);
	virtual void shallowCopyStructure(const Matrix *other);
	virtual void deepCopy(const Matrix *other);
	virtual void deepCopyStructure(const Matrix *other);
	virtual void uniformCombination(const Matrix *first, const Matrix *second, int nfirst, int nsecond);

	virtual void fill(double value);
	virtual void fillData(const Matrix *in);
	virtual void fillCombinedData(const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	virtual void apply(const Vector *in, Vector *out);
	virtual void apply(const Vectors *in, Vectors *out);

	virtual void scale(double alpha);
	virtual void add(double alpha, const Matrix *a);
	virtual void sum(double alpha, const Matrix *a, double beta, const Matrix *b);
	virtual void addToCombination(double alpha, const Matrix *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	virtual void fillDiagonal(Vector *diagonal) const;

	virtual double norm();

	const char* name() const { return "MatrixFETI"; };

	esint domains;
	Matrix **matrices;
protected:
	virtual Matrix* create();
};

}



#endif /* SRC_WRAPPERS_MATH_MATRIXFETI_H_ */
