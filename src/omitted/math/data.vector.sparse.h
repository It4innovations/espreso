
#ifndef SRC_WRAPPERS_MATH_DATAVECTORSPARSE_H_
#define SRC_WRAPPERS_MATH_DATAVECTORSPARSE_H_

namespace espreso {

class DataVectorDense;

struct _DataVectorSparse {
	esint size, nnz;
	esint *indices;
	double *vals;

	_DataVectorSparse();
	void alloc(esint size, esint nnz);
	void clear();
};

class DataVectorSparse: public _DataVectorSparse
{
public:
	static void combineIndices(esint *result, esint *first, esint *second, esint *firstend, esint *secondend, esint nfirst, esint nsecond);

	DataVectorSparse();
	DataVectorSparse(esint size, esint nnz);

	~DataVectorSparse();

	void shiftData(esint offset);
	void resize(esint size, esint nnz);
	void swap(DataVectorSparse *other);
	void shallowCopy(const DataVectorSparse *other);
	void shallowCopyStructure(const DataVectorSparse *other);
	void shallowCopyFromHolder(const DataVectorSparse *other, esint offset, esint nvectors);
	void deepCopy(const DataVectorSparse *other);
	void deepCopyStructure(const DataVectorSparse *other);
	void uniformCombination(const DataVectorSparse *first, const DataVectorSparse *second, int nfirst, int nsecond);

	void fill(double value);
	void fillPattern(esint *indices);
	void fillValues(double *vals);
	void fillDenseValues(double *vals);
	void fillCombinedValues(const DataVectorSparse *in, esint offset, esint nsize, esint sumsize);
	void fillValuesFromCombination(const DataVectorSparse *in, esint offset, esint nsize, esint sumsize);

	void addToCombination(double alpha, const DataVectorSparse *in, esint offset, esint nsize, esint sumsize);

protected:
	_DataVectorSparse _allocated;
};

}

#endif /* SRC_WRAPPERS_MATH_DATAVECTORSPARSE_H_ */
