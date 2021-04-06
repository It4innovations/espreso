
#ifndef SRC_WRAPPERS_MATH_DATAMATRIXIJV_H_
#define SRC_WRAPPERS_MATH_DATAMATRIXIJV_H_

namespace espreso {

class DataVectorDense;

struct _DataMatrixIJV {
	esint nrows, ncols, nnz;
	esint *rows, *cols;
	double *vals;

	_DataMatrixIJV();
	void alloc(esint nrows, esint ncols, esint nnz);
	void clear();
};

class DataMatrixIJV: public _DataMatrixIJV
{
public:
	static const esint indexing;
	static void combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond);

	DataMatrixIJV();
	DataMatrixIJV(esint nrows, esint ncols, esint nnz);
	DataMatrixIJV(const DataMatrixIJV &other);
	DataMatrixIJV(DataMatrixIJV &&other);
	DataMatrixIJV& operator=(const DataMatrixIJV &other);
	~DataMatrixIJV();

	void allowUpdating();
	void resize(esint nrows, esint ncols, esint nnz);
	void swap(DataMatrixIJV *other);
	void shallowCopy(const DataMatrixIJV *other);
	void shallowCopyStructure(const DataMatrixIJV *other);
	void deepCopy(const DataMatrixIJV *other);
	void deepCopyStructure(const DataMatrixIJV *other);
	void uniformCombination(const DataMatrixIJV *first, const DataMatrixIJV *second, int nfirst, int nsecond);

	void fillPattern(esint nnz, esint *rows, esint *cols);
	void fill(double value);
	void fillValues(esint nnz, double *vals);
	void fillCombinedValues(const DataMatrixIJV *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void addToCombination(double scale, const DataMatrixIJV *in, esint roffset, esint coffset, esint nsize, esint sumsize);
	void fillDiagonal(DataVectorDense *diagonal) const;
protected:
	_DataMatrixIJV _allocated;
};

}

#endif /* SRC_WRAPPERS_MATH_DATAMATRIXIJV_H_ */
