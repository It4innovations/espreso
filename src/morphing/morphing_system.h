#ifndef SRC_MORPHING_SYSTEM_MATRIX_H_
#define SRC_MORPHING_SYSTEM_MATRIX_H_

#include <limits>

#include "aca.h"
#include "CyclicGraphDecomposition.h"

#include "esinfo/mpiinfo.h"

#include "morphing_system_preconditioning.h"

namespace espreso {

/*
    represents a system matrix consisting of 2 blocks: A, B
    -----
    AB'
    B
    -----
*/
class MorphingMatrix{
public:

    MorphingMatrix();
    
    ~MorphingMatrix();
    
    MorphingMatrix(
        int dim,
        bool tr_translate,
        bool tr_scale,
        bool tr_rotate,
        int regularization_poly_degree,
        const std::vector<Point> &rPoints,
        const std::vector<double> &rDisplacement,
        const RBFTargetConfiguration &configuration,
        double aca_eps,
        double aca_eta,
        esint base_cluster_size,
        
        std::vector<double> &morphing_rhs
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
    
    //performs y = this->P * x * alpha + y * beta
    void applyPreconditioner(
        const double * x_global, 
        double *y_global, 
        double alpha, 
        double beta, 
        bool transpose
    ) const ;
    
    void calculateMorphingError(
        const std::vector<double> &sol,
        const std::vector<double> &rhs,
        double &error_morph_x,
        double &error_morph_y,
        double &error_morph_z,
        double &orthogonality_error_x,
        double &orthogonality_error_y,
        double &orthogonality_error_z
    ) const ;
    
    void transformPoint(
        Point &point_to_morph,
        const Point &point_origin,
        const std::vector<double> &coefficients,
        const std::vector<Point> &rPoints,
        const RBFTargetConfiguration &configuration
    ) const;
    
    esint getNRows() const;

    esint getNCols() const;
    
    double getCompressionRatio() const;
    
    void print(const std::vector<double> &rhs) const;
    
    double getRelativeError(const std::vector<double> &ref_data) const;
    
    double getRelativeErrorMax(const std::vector<double> &ref_data) const;
    
private:

    MorphingMatrixPreconditioner *P = nullptr;

    std::vector<matrix_ACA*> A;
    std::vector<esint> A_inner_block_idx;
    std::vector<esint> A_outer_block_idx;

    std::vector<double> B;
    
    
    std::vector<esint> global_DOFs_start;
    std::vector<esint> global_DOFs_end;
    
    // std::vector<double> apply_data_mem;
    
    esint nrows = 0;
    esint ncols = 0;
    esint dimension = 0;
    esint n_regularization = 0;
    esint regularization_polynomial_degree = 0;
    
    bool use_transformation_translation;
    bool use_transformation_scaling;
    bool use_transformation_rotation;
    
    std::vector<double> scaling = {1.0f, 1.0f, 1.0f};
    std::vector<double> system_shift = {0.0f, 0.0f, 0.0f};
    std::vector<std::vector<double>> rotation_matrix = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};
    
    double calculateConstantShift(
        const Point &Pr
    ) const;
    
    /*
    computes B in row-major format
    */
    void recalculateB(
        const std::vector<Point> &rPoints
    );

    void define_geometry_transformations(
        const std::vector<Point> &sPoints
    );

    void transform_points(
        std::vector<Point> &sPoints
    ) const;

    void transform_points_inverse(
        std::vector<Point> &sPoints
    ) const;

    void transform_point(
        Point &p
    ) const;

    void transform_point_inverse(
        Point &p
    ) const;
    
    void transform_displacements(
        std::vector<double> &displacements
    ) const;

    void transform_displacement(
        double *displacement
    ) const;

};


};//end of namespace espreso

#endif /* SRC_MORPHING_SYSTEM_MATRIX_H_ */
