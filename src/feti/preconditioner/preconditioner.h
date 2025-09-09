
#ifndef SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_
#define SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"

namespace espreso {

template <typename T>
struct Preconditioner {
    static Preconditioner<T>* create(FETI<T> &feti, const step::Step &step);

    Preconditioner(FETI<T> &feti): feti(feti) {}
    virtual ~Preconditioner() {}

    virtual void setup() {}
    virtual size_t get_wss_gpu_persistent() { return 0; }
    virtual size_t get_wss_gpu_internal() { return 0; }
    virtual void set_ws_gpu_persistent(void * ws_gpu_persistent_) {}

    virtual void info() =0;
    virtual void set(const step::Step &step) { }
    virtual void update(const step::Step &step) =0;

    // y = S * x
    virtual void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y) =0;

    virtual bool isset() const { return true; }

protected:
    FETI<T> &feti;
};

}

#endif /* SRC_FETI_PRECONDITIONER_PRECONDITIONER_H_ */
