
#ifndef SRC_WRAPPERS_MATH_DATAVECTORDENSE_H_
#define SRC_WRAPPERS_MATH_DATAVECTORDENSE_H_

namespace espreso {

class DataVectorSparse;

struct _DataVectorDense {
	esint size;
	double *vals;

	_DataVectorDense();
	void realloc(esint size);
	void clear();

protected:
	esint _maxsize;
};

class DataVectorDense: public _DataVectorDense
{
public:
	DataVectorDense();
	DataVectorDense(esint size);
	DataVectorDense(esint size, double *vals);
	~DataVectorDense();

	void shiftData(esint offset);
	void resize(esint size);
	void swap(DataVectorDense *other);
	void shallowCopy(const DataVectorDense *other);
	void shallowCopyStructure(const DataVectorDense *other);
	void shallowCopyFromHolder(const DataVectorDense *other, esint offset, esint nvectors);
	void deepCopy(const DataVectorDense *other);
	void deepCopyStructure(const DataVectorDense *other);
	void uniformCombination(const DataVectorDense *first, const DataVectorDense *second);

	void fill(double value);
	void fillValues(double *vals);
	void fillSparseValues(esint nnz, esint *indices, double *vals);
	void fillCombinedValues(const DataVectorDense *in, esint offset, esint nsize, esint sumsize);
	void fillValuesFromCombination(const DataVectorDense *in, esint offset, esint nsize, esint sumsize);

	void addToCombination(double alpha, const DataVectorDense *in, esint offset, esint nsize, esint sumsize);

protected:
	_DataVectorDense _allocated;
};

}


#endif /* SRC_WRAPPERS_MATH_DATAVECTORDENSE_H_ */
