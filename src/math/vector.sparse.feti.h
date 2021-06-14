
#ifndef SRC_WRAPPERS_MATH_VECTORSPARSEFETI_H_
#define SRC_WRAPPERS_MATH_VECTORSPARSEFETI_H_

#include "vector.feti.h"
#include "vector.sparse.h"

namespace espreso {

class VectorSparseFETI: public VectorFETI
{
public:
	VectorSparseFETI();
	~VectorSparseFETI();

	VectorSparse& operator[](esint domain) { return *reinterpret_cast<VectorSparse*>(vectors[domain]); }
	const VectorSparse& operator[](esint domain) const { return *reinterpret_cast<VectorSparse*>(vectors[domain]); }

	VectorSparse* at(esint domain) { return reinterpret_cast<VectorSparse*>(vectors[domain]); }
	const VectorSparse* at(esint domain) const { return reinterpret_cast<VectorSparse*>(vectors[domain]); }

	VectorSparseFETI* copy();

	void fromFETI(VectorFETI *other) const;

	void filter(const DataDecomposition *other, esint nindices, esint *indices);

protected:
	Vector* create();
	const char* name() const { return "VectorSparseFETI"; }
};

class VectorsSparseFETI: public VectorsFETI
{
public:
	// used by print
	struct Domain {
		Domain(const VectorsSparseFETI* data, esint domain): data(data), domain(domain) {}
		const VectorsSparseFETI* data;
		esint domain;
	};

	VectorSparseFETI* holder() { return reinterpret_cast<VectorSparseFETI*>(_holder); }
	const VectorSparseFETI* holder() const { return reinterpret_cast<VectorSparseFETI*>(_holder); }

	VectorSparseFETI& operator[](esint index) { return *reinterpret_cast<VectorSparseFETI*>(_vectors[index]); }
	const VectorSparseFETI& operator[](esint index) const { return *reinterpret_cast<VectorSparseFETI*>(_vectors[index]); }

	VectorSparseFETI* at(esint index) { return reinterpret_cast<VectorSparseFETI*>(_vectors[index]); }
	const VectorSparseFETI* at(esint index) const { return reinterpret_cast<VectorSparseFETI*>(_vectors[index]); }

	VectorsSparseFETI();
	~VectorsSparseFETI();

	VectorsSparseFETI* copy();

	void initDomains(DataDecomposition::DUPLICATION duplications, esint domains);
	void resizeDomain(esint domain, esint size, esint nnz);
	void filter(const DataDecomposition *other, esint nindices, esint *indices);

	const char* name() const { return "VectorsSparseFETI"; }

protected:
	Vector* create();
};

}

#endif /* SRC_WRAPPERS_MATH_VECTORSPARSEFETI_H_ */
