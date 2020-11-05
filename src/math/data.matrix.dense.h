
#ifndef SRC_WRAPPERS_MATH_DATAMATRIXDENSE_H_
#define SRC_WRAPPERS_MATH_DATAMATRIXDENSE_H_

namespace espreso {

class DataVectorDense;

struct _DataMatrixDense {
	esint nrows, ncols;
	double *vals;

	_DataMatrixDense();
	void realloc(esint nrows, esint ncols);
	void clear();

protected:
	esint _maxvals;
};

class DataMatrixDense: public _DataMatrixDense
{
public:
	DataMatrixDense();
	DataMatrixDense(esint nrows, esint ncols);
	DataMatrixDense(esint nrows, esint ncols, double *vals);
	DataMatrixDense(const DataMatrixDense &other);
	~DataMatrixDense();

	void set(esint nrows, esint ncols, double *vals);
	void resize(esint nrows, esint ncols);
	void swap(DataMatrixDense *other);
	void shallowCopy(const DataMatrixDense *other);
	void deepCopy(const DataMatrixDense *other);
	void uniformCombination(const DataMatrixDense *first, const DataMatrixDense *second, int nfirst, int nsecond);

	void fill(double value);
	void fillValues(double *vals);
	void fillCombinedValues(const DataMatrixDense *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void addToCombination(double alpha, const DataMatrixDense *in, esint roffset, esint coffset, esint nsize, esint sumsize);
	void fillDiagonal(DataVectorDense *diagonal) const;

protected:
	_DataMatrixDense _allocated;
};

}



#endif /* SRC_WRAPPERS_MATH_DATAMATRIXDENSE_H_ */
