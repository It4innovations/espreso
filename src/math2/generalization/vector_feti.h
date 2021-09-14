
#ifndef SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_
#define SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_

#include "vector_base.h"
#include "math2/math2.h"
#include "math2/utils/utils_feti.h"
#include "wrappers/mpi/communication.h"

#include <vector>

namespace espreso {

template <template<typename> typename Vector, typename T>
class Vector_FETI: public Vector_Base<T>
{
public:
	~Vector_FETI()
	{

	}

//	Vector_FETI* copy()
//	{
//		Vector_FETI<Vector, T> *m = new Vector_FETI<Vector, T>();
//		m->domains.reserve(domains.size());
//		for (size_t d = 0; d < domains.size(); ++d) {
//			m->domains.push_back(domains[d]);
//		}
//		return m;
//	}

	Vector_FETI* copyPattern()
	{
		Vector_FETI<Vector, T> *m = new Vector_FETI<Vector, T>();
		m->domains.resize(domains.size());
		for (size_t d = 0; d < domains.size(); ++d) {
			m->domains[d].pattern(domains[d]);
		}
		return m;
	}

	void store(const char *file)
	{
		math::store(*this, file);
	}

	void store(std::vector<T> &output)
	{
		//
	}

	void fill(const T &value)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::fill(domains[d], value);
		}
	}

	void fillData(const Vector_Base<T> *in)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(in)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::copy(domains[d], static_cast<const Vector_FETI<Vector, T>*>(in)->domains[d]);
			}
		}
	}

	void fillData(const Vector_Base<T> *in, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(in)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::copy(domains[d], static_cast<const Vector_FETI<Vector, T>*>(in)->domains[d], offset, size, step);
			}
		}
	}

	void scale(const T &alpha)
	{
		#pragma omp parallel for
		for (size_t d = 0; d < domains.size(); ++d) {
			math::scale(alpha, domains[d]);
		}
	}

	void add(const T &alpha, const Vector_Base<T> *a)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(a)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::add<T>(domains[d], alpha, static_cast<const Vector_FETI<Vector, T>*>(a)->domains[d]);
			}
		}
	}

	void add(const T &alpha, const Vector_Base<T> *a, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(a)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::add<T>(domains[d], alpha, static_cast<const Vector_FETI<Vector, T>*>(a)->domains[d], offset, size, step);
			}
		}
	}

	void sum(const T &alpha, const Vector_Base<T> *a, const T &beta, const Vector_Base<T> *b)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(a) && dynamic_cast<const Vector_FETI<Vector, T>*>(b)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::sum(domains[d], alpha, static_cast<const Vector_FETI<Vector, T>*>(a)->domains[d], beta, static_cast<const Vector_FETI<Vector, T>*>(b)->domains[d]);
			}
		}
	}

	void sum(const T &alpha, const Vector_Base<T> *a, const T &beta, const Vector_Base<T> *b, int offset, int size, int step)
	{
		if (dynamic_cast<const Vector_FETI<Vector, T>*>(a) && dynamic_cast<const Vector_FETI<Vector, T>*>(b)) {
			#pragma omp parallel for
			for (size_t d = 0; d < domains.size(); ++d) {
				math::sum(domains[d], alpha, static_cast<const Vector_FETI<Vector, T>*>(a)->domains[d], beta, static_cast<const Vector_FETI<Vector, T>*>(b)->domains[d], offset, size, step);
			}
		}
	}

	void addTo(const T &alpha, Vector_Sparse<T> *a) const
	{

	}

	T norm()
	{
		return 0;
	}

	T max()
	{
		return 0;
	}

	T absmax()
	{
		return 0;
	}

	T dot(const Vector_Base<T> *other)
	{
		return 0;
	}

	std::vector<Vector<T> > domains;
	DomainDecomposition *decomposition;
};

}

#endif /* SRC_MATH2_GENERALIZATION_VECTOR_FETI_H_ */
