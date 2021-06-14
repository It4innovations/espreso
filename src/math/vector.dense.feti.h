
#ifndef SRC_WRAPPERS_MATH_VECTORDENSEFETI_H_
#define SRC_WRAPPERS_MATH_VECTORDENSEFETI_H_

#include "vector.feti.h"
#include "vector.dense.h"

namespace espreso {

class VectorDenseFETI: public VectorFETI
{
public:
	VectorDenseFETI();
	~VectorDenseFETI();

	VectorDense& operator[](esint domain) { return *reinterpret_cast<VectorDense*>(vectors[domain]); }
	const VectorDense& operator[](esint domain) const { return *reinterpret_cast<VectorDense*>(vectors[domain]); }

	VectorDense* at(esint domain) { return reinterpret_cast<VectorDense*>(vectors[domain]); }
	const VectorDense* at(esint domain) const { return reinterpret_cast<VectorDense*>(vectors[domain]); }

	VectorDenseFETI* copy();

	double norm();
	double max();
	double absmax();
	double dot(const Vector *other);

	void allGather() const;
	void averageDuplications();
	void sumDuplications();

	void toFETI(VectorFETI *other) const;

	const char* name() const { return "VectorDenseFETI"; }

protected:
	Vector* create();
};

class VectorsDenseFETI: public VectorsFETI
{
public:
	// used by print
	struct Domain {
		Domain(const VectorsDenseFETI* data, esint domain): data(data), domain(domain) {}
		const VectorsDenseFETI* data;
		esint domain;
	};

	VectorDenseFETI* holder() { return reinterpret_cast<VectorDenseFETI*>(_holder); }
	const VectorDenseFETI* holder() const { return reinterpret_cast<VectorDenseFETI*>(_holder); }

	VectorDenseFETI& operator[](esint index) { return *reinterpret_cast<VectorDenseFETI*>(_vectors[index]); }
	const VectorDenseFETI& operator[](esint index) const { return *reinterpret_cast<VectorDenseFETI*>(_vectors[index]); }

	VectorDenseFETI* at(esint index) { return reinterpret_cast<VectorDenseFETI*>(_vectors[index]); }
	const VectorDenseFETI* at(esint index) const { return reinterpret_cast<VectorDenseFETI*>(_vectors[index]); }

	VectorsDenseFETI();
	~VectorsDenseFETI();

	VectorsDenseFETI* copy();

	void initDomains(DataDecomposition::DUPLICATION duplications, esint domains);
	void resizeDomain(esint domain, esint size);
	void fillDecomposition(esint rank, esint ranks, esint nneighbors, esint *distribution, int *neighbors, const serializededata<esint, DI> *dmap);

	void averageDuplications();
	void sumDuplications();

	const char* name() const { return "VectorsDenseFETI"; }

protected:
	Vector* create();
};

}

#endif /* SRC_WRAPPERS_MATH_VECTORDENSEFETI_H_ */
