
#ifndef SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_H_
#define SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_H_

#include "analysis/math/matrix_distributed.h"
#include "analysis/math/vector_distributed.h"
#include "config/ecf/linearsolver/mklpdss.h"

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

    bool set(const Matrix_Distributed<T> &A);
    bool update(const Matrix_Distributed<T> &A);
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
