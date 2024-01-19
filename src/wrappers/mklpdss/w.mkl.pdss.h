
#ifndef SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_H_
#define SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_H_

#include "config/ecf/linearsolver/mklpdss.h"
#include "analysis/linearsystem/matrices/matrix_distributed.h"
#include "analysis/linearsystem/matrices/vector_distributed.h"

namespace espreso {

template<typename T>
struct MKLPDSSDataHolder;

template<typename T>
class MKLPDSS {
public:
	MKLPDSS(MKLPDSSConfiguration &configuration)
	: configuration(configuration), external(nullptr)
	{
		check();
	}

	~MKLPDSS()
	{
		call(-1);
		clear();
	}

	bool set(const Matrix_Distributed<Matrix_CSR, T> &A);
	bool update(const Matrix_Distributed<Matrix_CSR, T> &A);
	bool solve(const Vector_Distributed<Vector_Dense, T> &b, Vector_Distributed<Vector_Dense, T> &x);

	MKLPDSSConfiguration &configuration;
	MKLPDSSDataHolder<T> *external;

protected:
	bool call(int phase);
	void check();
	void clear();
};

}



#endif /* SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_H_ */
