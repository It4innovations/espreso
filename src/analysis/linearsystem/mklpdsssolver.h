
#ifndef SRC_ANALYSIS_LINEARSYSTEM_MKLPDSSSOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_MKLPDSSSOLVER_H_

#include "directsolver.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "wrappers/mklpdss/w.mkl.pdss.h"

namespace espreso {

template <typename T>
struct MKLPDSSLinearSystemSolver: DirectLinearSystemSolver<T> {

    MKLPDSSLinearSystemSolver(MKLPDSSConfiguration &configuration): mklpdss(configuration) {}
    ~MKLPDSSLinearSystemSolver() {}

    void set(step::Step &step)
    {
        mklpdss.set(this->A);
    }

    void update(step::Step &step)
    {
        if (this->A.updated || this->b.updated || this->dirichlet.updated) {
            this->setDirichlet();
            mklpdss.update(this->A);
            this->A.updated = this->b.updated = true; // set dirichlet update both
        }
        this->dirichlet.updated = false;

        if (info::ecf->output.print_matrices) {
            eslog::storedata(" STORE: system/{A, b, dirichlet}\n");
            math::store(this->A, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
            math::store(this->b, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
            math::store(this->dirichlet, utils::filename(utils::debugDirectory(step) + "/system", "BC").c_str());
        }
    }

    bool solve(step::Step &step)
    {
        if (mklpdss.solve(this->b, this->x)) {
            this->x.scatter();
            if (info::ecf->output.print_matrices) {
                eslog::storedata(" STORE: system/{x}\n");
                math::store(this->x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
            }
            return true;
        }
        return false;
    }

private:
    MKLPDSS<T> mklpdss;
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_MKLPDSSSOLVER_H_ */
