#ifndef SRC_MORPHING_SYSTEM_MATRIX__PRECOND_H_
#define SRC_MORPHING_SYSTEM_MATRIX__PRECOND_H_

#include <vector>
#include <iostream>

#include "math/math.h"

#include "basis/evaluator/evaluator.h"

#include "basis/containers/point.h"

#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"

#include "basis/structures/kdtree.h"


namespace espreso {

class MorphingMatrixPreconditioner{
public:

    MorphingMatrixPreconditioner();
    
    ~MorphingMatrixPreconditioner();
    
    MorphingMatrixPreconditioner(
        esint npoints_total,
        esint nregularization,
        esint dof_shift,
        const std::vector<Point> &points_local,
        const RBFTargetConfiguration &configuration
    );

    void printData() const;
    
    //performs y = this * x * alpha + y * beta
    void apply(
        const double * x_global, 
        double *y_global, 
        double alpha, 
        double beta, 
        bool transpose
    ) const ;
    
    esint getNRows() const;

    esint getNCols() const;
    
private:

    esint nrows = 0;
    esint ncols = 0;
    esint nreg_dim = 0;
    esint dim_id = 0;

    std::vector<double> basis_vector;
    std::vector<double> sparse_system;
    std::vector<double> sparse_system_mem;
    std::vector<double> inverse_row;
    
    std::vector<esint> row_indices;
    std::vector<esint> col_indices;
    std::vector<double> values;
        

    void prepareMatrixSparse(
        const std::vector<Point> &points,
        const RBFTargetConfiguration &configuration
    );
    
    void updateInverseRowSparseOnly(esint idx);
    
    void updateInverseRow(
        const Point & p,
        const std::vector<Point> &sparse_points,
        const RBFTargetConfiguration &configuration
    );
    
};


};//end of namespace espreso

#endif /* SRC_MORPHING_SYSTEM_MATRIX__PRECOND_H_ */
