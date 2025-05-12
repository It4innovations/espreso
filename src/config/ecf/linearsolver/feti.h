
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_

#include "config/description.h"
#include "autoopt.h"
#include "config/ecf/linearsolver/dual_operator_gpu_config.h"

#include <cstddef>

namespace espreso {


struct FETIConfiguration: public ECFDescription {

    enum class METHOD {
        /// Total FETI
        TOTAL_FETI = 0,
        /// Hybrid Total FETI
        HYBRID_FETI = 1,
    };

    enum class DUAL_OPERATOR {
        IMPLICIT = 0,
        EXPLICIT,
        EXPLICIT_GPU,
        IMPLICIT_GPU,
        EXPLICIT_SC,
        EXPLICIT_SC_GPUAPPLY,
        EXPLICIT_SCTRIA,
        EXPLICIT_SCTRIA_GPU,
        EXPLICIT_GENERALSC_CPU,
        EXPLICIT_GENERALSC_GPU
    };

    enum class PROJECTOR {
        ORTHOGONAL = 0,
        CONJUGATE,
    };

    enum PROJECTOR_OPT {
        DEFAULT      = 1 << 0,
        WITH_FACTORS = 1 << 1,
        FULL         = 1 << 2
    };

    enum class ITERATIVE_SOLVER {
        PCG = 0,
        orthogonalPCG,
        GMRES,
        BICGSTAB,
        SMALBE,
        MPRGP
    };

    enum class PRECONDITIONER {
        NONE,
        LUMPED,
        DIRICHLET
    };

    enum class STOPPING_CRITERION {
        RELATIVE,
        ABSOLUTE,
        ARIOLI // https://www.researchgate.net/publication/30411264_A_stopping_criterion_for_the_Conjugate_Gradient_algorithm_in_a_finite_element_method_framework
    };

    enum class REGULARIZATION {
        ANALYTIC = 0,
        ALGEBRAIC = 1,
        SVD
    };

    enum class FIX_POINTS {
        METIS_CENTERS,
        RANDOM,
        SPHERICAL
    };

    double precision;
    size_t print_iteration;
    size_t max_iterations;
    size_t max_stagnation;
    bool exit_on_nonconvergence;
    bool check_input_matrices;

    METHOD method;
    DUAL_OPERATOR dual_operator;
    PROJECTOR projector;
    PROJECTOR_OPT projector_opt;
    ITERATIVE_SOLVER iterative_solver;
    PRECONDITIONER preconditioner;
    REGULARIZATION regularization;
    FIX_POINTS fix_points;
    STOPPING_CRITERION stopping_criterion;

    int exhaustive_info;

    double precision_in, precision_set;
    size_t max_iterations_in;

    int power_maxit;
    double power_precision;

    // SMALSE
    double gamma, M, rho, eta, beta, alpham, delta, rtol;
    bool halfstep, exp_projgrad, prop_projgrad, proj_grad, gradproj, optimset;
    int th, no_enlarg_exp, no_enlarg_prop;

    size_t sc_size;

    AutoOptimizationConfiguration auto_optimization;

    DualOperatorGpuConfig dual_operator_gpu_config;

    FETIConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_LINEARSOLVER_FETI_H_ */
