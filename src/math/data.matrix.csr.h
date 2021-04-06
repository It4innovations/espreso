
#ifndef SRC_WRAPPERS_MATH_DATAMATRIXCSR_H_
#define SRC_WRAPPERS_MATH_DATAMATRIXCSR_H_

namespace espreso {

class DataVectorDense;
namespace MATH { struct CSRHandler; }

struct _DataMatrixCSR {
	esint nrows, ncols, nnz;
	esint *rows, *cols;
	double *vals;

	_DataMatrixCSR();
	void alloc(esint nrows, esint ncols, esint nnz);
	void clear();
};

class DataMatrixCSR: public _DataMatrixCSR
{
	friend class DataMV;
public:
	static const esint indexing;
	static void combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond);

	DataMatrixCSR();
	DataMatrixCSR(esint nrows, esint ncols, esint nnz);
	DataMatrixCSR(const DataMatrixCSR &other);
	DataMatrixCSR(const MATH::CSRHandler *handler);
	DataMatrixCSR(DataMatrixCSR &&other);
	DataMatrixCSR& operator=(const DataMatrixCSR &other);
	DataMatrixCSR& operator=(const MATH::CSRHandler *handler);
	~DataMatrixCSR();

	void allowUpdating();
	void resize(esint nrows, esint ncols, esint nnz);
	void swap(DataMatrixCSR *other);
	void shallowCopy(const DataMatrixCSR *other);
	void shallowCopyStructure(const DataMatrixCSR *other);
	void deepCopy(const DataMatrixCSR *other);
	void deepCopyStructure(const DataMatrixCSR *other);
	void uniformCombination(const DataMatrixCSR *first, const DataMatrixCSR *second, int nfirst, int nsecond);

	void fillPattern(esint nrows, esint *rows, esint *cols);
	void fill(double value);
	void fillValues(esint nnz, double *vals);
	void fillCombinedValues(const DataMatrixCSR *in, esint roffset, esint coffset, esint nsize, esint sumsize);

	void addToCombination(double scale, const DataMatrixCSR *in, esint roffset, esint coffset, esint nsize, esint sumsize);
	void fillDiagonal(DataVectorDense *diagonal) const;
protected:
	_DataMatrixCSR _allocated;
};

}

#endif /* SRC_WRAPPERS_MATH_DATAMATRIXCSR_H_ */
