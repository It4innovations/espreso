
#ifndef SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_
#define SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_

#include "elementmapping.h"

#include "math2/primitives/vector_dense.h"
#include "math2/primitives/matrix_info.h"
#include "math2/primitives/matrix_csr.h"
#include "math2/generalization/vector_feti.h"
#include "math2/generalization/matrix_feti.h"
#include "math2/utils/utils_feti.h"

#include <vector>

namespace espreso {

struct UniformNodesFETIPattern {

	struct RegionInfo {
		esint size;
		std::vector<esint> row, column; // row, column indices
		std::vector<esint> K, f; // permutations
	};

	UniformNodesFETIPattern();
	~UniformNodesFETIPattern();

	void set(int dofs, Matrix_Type type);

	template<typename T>
	void fill(Vector_FETI<Vector_Dense, T> &v)
	{
		v.domains.resize(elements.size());
		#pragma omp parallel for
		for (size_t d = 0; d < elements.size(); ++d) {
			v.domains[d].resize(elements[d].size);
		}
	}

	template<typename T>
	void fill(Matrix_FETI<Matrix_CSR, T> &m)
	{
		m.domains.resize(elements.size());
		#pragma omp parallel for
		for (size_t d = 0; d < elements.size(); ++d) {
			m.domains[d].type = m.type;
			m.domains[d].shape = m.shape;
			m.domains[d].resize(elements[d].size, elements[d].size, elements[d].K.size());
			fillCSR(d, m.domains[d].rows, m.domains[d].cols);
		}
	}

	template<typename T>
	void setMap(Matrix_FETI<Matrix_CSR, T> *m) const
	{

	}

	template<typename T>
	void setMap(Vector_FETI<Vector_Dense, T> *m) const
	{

	}

	void fillCSR(int domain, esint *rows, esint *cols);

	serializededata<esint, DI> *dmap;

	std::vector<RegionInfo> elements; // RegionInfo per domain
	std::vector<std::vector<RegionInfo> > bregion; // RegionInfo per domain per boundary region
};

}

#endif /* SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_ */
