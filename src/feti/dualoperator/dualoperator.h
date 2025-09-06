
#ifndef SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_
#define SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"
#include "feti/common/matrix_dual.h"
#include "math/wrappers/math.spsolver.h"

#include <vector>

namespace espreso {

struct DualOperatorInfo {
    size_t rows, nnzA, nnzL;
//     size_t memoryL;
    size_t dualA, surfaceA;
};

enum struct DualOperatorStrategy {
    IMPLICIT,
    EXPLICIT
};

template <typename T>
class DualOperator {
public:
    static DualOperator<T>* create(FETI<T> &feti, const step::Step &step);

    DualOperator(FETI<T> &feti): feti(feti), infoPrinted(false) {}
    virtual ~DualOperator() {}

    virtual void setup() {}
    virtual size_t get_wss_gpu_persistent() { return 0; }
    virtual size_t get_wss_gpu_internal() { return 0; }
    virtual void set_ws_gpu_persistent(void * ws_gpu_persistent_) {}

    virtual void info() =0;
    virtual void set(const step::Step &step) =0;
    virtual void update(const step::Step &step) =0;

    // y = F * x
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;
    virtual void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y) =0;
    virtual void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<int> &filter) {}
    virtual void apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y, const std::vector<std::vector<int> > &filter) {}

    // y = K+(f - Bt * x)
    virtual void toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y) =0;
    virtual void BtL(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
    {
        eslog::warning("NOT-IMPLEMENTED BT * LAMBDA\n");
    }

    void estimateMaxEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations);
    void estimateMaxProjectedEigenValue(double &lambda, int &iterations, double epsilon, int maxIterations, double rho = 1, double normPFP = 1);

    FETI<T> &feti;
    Vector_Dual<T> d;

protected:
    bool infoPrinted;

    void getInitVector(Vector_Dual<T> &v);

    void reduceInfo(std::vector<DirectSparseSolver<T> > &KSolver, DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max);
    void printInfo(std::vector<DirectSparseSolver<T> > &KSolver, DualOperatorInfo &sum, DualOperatorInfo &min, DualOperatorInfo &max);
};

}

#endif /* SRC_FETI_DUALOPERATOR_DUALOPERATOR_H_ */
