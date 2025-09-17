
#include "feti.h"
#include "config/configuration.hpp"
#include "esinfo/envinfo.h"

espreso::FETIConfiguration::FETIConfiguration()
{
    method = METHOD::TOTAL_FETI;
    REGISTER(method, ECFMetaData()
            .setdescription({ "Type" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("TOTAL_FETI").setdescription("FETI with Dirichlet in B1."))
            .addoption(ECFOption().setname("HYBRID_FETI").setdescription("HYBRID FETI with Dirichlet in B1.")));

    dual_operator = DUAL_OPERATOR::IMPLICIT;
    REGISTER(dual_operator, ECFMetaData()
            .setdescription({ "Type" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("IMPLICIT").setdescription("Implicit F with sparse direct solver."))
            .addoption(ECFOption().setname("EXPLICIT").setdescription("Explicit F with BLAS."))
            .addoption(ECFOption().setname("EXPLICIT_GPU").setdescription("Explicit F on GPU."))
            .addoption(ECFOption().setname("IMPLICIT_GPU").setdescription("Implicit F on GPU."))
            .addoption(ECFOption().setname("EXPLICIT_GENERALSCHUR_CPU").setdescription("Excplicit F using general SC operation on CPU."))
            .addoption(ECFOption().setname("EXPLICIT_GENERALSCHUR_GPU").setdescription("Excplicit F using general SC operation on GPU."))
            .addoption(ECFOption().setname("IMPLICIT_GENERALSPARSESOLVER_CPU").setdescription("Excplicit F using general sparse solver operation on CPU.")));

    projector = PROJECTOR::ORTHOGONAL;
    REGISTER(projector, ECFMetaData()
            .setdescription({ "Type" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("ORTHOGONAL").setdescription("Orthogonal projector Gt * inv(G * Gt) * G."))
            .addoption(ECFOption().setname("CONJUGATE").setdescription("Conjugate projector Gt * inv(G * F * Gt) * G * F.")));

    projector_opt = PROJECTOR_OPT::DEFAULT;
    REGISTER(projector_opt, ECFMetaData()
            .setdescription({ "Type" })
            .setdatatype({ ECFDataType::ENUM_FLAGS })
            .addoption(ECFOption().setname("DEFAULT").setdescription("Default configuration according to projector type."))
            .addoption(ECFOption().setname("WITH_FACTORS").setdescription("Orthogonal projector Gt * inv(U) * inv(L) * G."))
            .addoption(ECFOption().setname("FULL").setdescription("Total FETI projector.")));

    preconditioner = PRECONDITIONER::DIRICHLET;
    REGISTER(preconditioner, ECFMetaData()
            .setdescription({ "Preconditioner" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("NONE").setdescription("Without precodition."))
            .addoption(ECFOption().setname("LUMPED").setdescription("Lumped precodition."))
            .addoption(ECFOption().setname("DIRICHLET").setdescription("Dirichlet precodition."))
            .addoption(ECFOption().setname("DIRICHLET_GENERALSCHUR_CPU").setdescription("Dirichlet precodition, generalized SC on CPU."))
            .addoption(ECFOption().setname("DIRICHLET_GENERALSCHUR_GPU").setdescription("Dirichlet precodition, generalized SC on GPU.")));

    precision = 1e-5;
    REGISTER(precision, ECFMetaData()
            .setdescription({ "Precision" })
            .setdatatype({ ECFDataType::FLOAT }));

    print_iteration = 10;
    REGISTER(print_iteration, ECFMetaData()
            .setdescription({ "Print only iterations divided by this number." })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    max_iterations = 0;
    REGISTER(max_iterations, ECFMetaData()
            .setdescription({ "Max iterations" })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    max_stagnation = 50;
    REGISTER(max_stagnation, ECFMetaData()
            .setdescription({ "Max stagnation steps." })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    exit_on_nonconvergence = true;
    REGISTER(exit_on_nonconvergence, ECFMetaData()
            .setdescription({ "Finish when FETI solver does not converge." })
            .setdatatype({ ECFDataType::BOOL }));

    check_input_matrices = false;
    REGISTER(check_input_matrices, ECFMetaData()
            .setdescription({ "Check input matrices with eigenvalues and SVD decomposition." })
            .setdatatype({ ECFDataType::BOOL }));

    iterative_solver = ITERATIVE_SOLVER::PCG;
    REGISTER(iterative_solver, ECFMetaData()
            .setdescription({ "Iterative solver" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("PCG").setdescription("Projected conjugate gradients."))
            .addoption(ECFOption().setname("orthogonalPCG").setdescription("Orthogonal PCG."))
            .addoption(ECFOption().setname("GMRES").setdescription("GMRES solver."))
            .addoption(ECFOption().setname("BiCGStab").setdescription("BiCGStab solver."))
            .addoption(ECFOption().setname("SMALBE").setdescription("SMALBE with MPRGP."))
            .addoption(ECFOption().setname("MPRGP").setdescription("MPRGP."))
            );

    regularization = REGULARIZATION::ANALYTIC;
    REGISTER(regularization, ECFMetaData()
            .setdescription({ "Regularization" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("ANALYTIC").setdescription("Analytic regularization provided by a particular physics."))
            .addoption(ECFOption().setname("ALGEBRAIC").setdescription("Regularization based on NULL PIVOTS."))
            .addoption(ECFOption().setname("SVD").setdescription("SVD decomposition.")));

    fix_points = FIX_POINTS::METIS_CENTERS;
    REGISTER(fix_points, ECFMetaData()
            .setdescription({ "Fix point computation" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("METIS_CENTERS").setdescription("Centers of domains computed by METIS."))
            .addoption(ECFOption().setname("RANDOM").setdescription("Random points."))
            .addoption(ECFOption().setname("SPHERICAL").setdescription("Random points in the sphere surface.")));

    stopping_criterion = STOPPING_CRITERION::RELATIVE;
    REGISTER(stopping_criterion, ECFMetaData()
            .setdescription({ "Stopping criterion" })
            .setdatatype({ ECFDataType::OPTION })
            .addoption(ECFOption().setname("RELATIVE").setdescription("Relative stopping."))
            .addoption(ECFOption().setname("ABSOLUTE").setdescription("Absolute stopping.")));

    exhaustive_info = 0;
    REGISTER(exhaustive_info, ECFMetaData()
            .setdescription({ "Print all FETI solver input properties." })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    precision_in = 1e-6;
    REGISTER(precision_in, ECFMetaData()
            .setdescription({ "SMALSE inner precision" })
            .setdatatype({ ECFDataType::FLOAT }));

    precision_set = 1e-9;
    REGISTER(precision_set, ECFMetaData()
            .setdescription({ "Precision for active and free sets." })
            .setdatatype({ ECFDataType::FLOAT }));

    max_iterations_in = 500;
    REGISTER(max_iterations_in, ECFMetaData()
            .setdescription({ "SMALSE inner max iterations" })
            .setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

    gamma = 1;
    REGISTER(gamma, ECFMetaData()
            .setdescription({ "SMALSE gamma" })
            .setdatatype({ ECFDataType::FLOAT }));
    M = 1;
    REGISTER(M, ECFMetaData()
            .setdescription({ "SMALSE M" })
            .setdatatype({ ECFDataType::FLOAT }));
    rho = 1;
    REGISTER(rho, ECFMetaData()
            .setdescription({ "SMALSE rho" })
            .setdatatype({ ECFDataType::FLOAT }));
    eta = 1;
    REGISTER(eta, ECFMetaData()
            .setdescription({ "SMALSE eta" })
            .setdatatype({ ECFDataType::FLOAT }));
    beta = 0.8;
    REGISTER(beta, ECFMetaData()
            .setdescription({ "SMALSE beta" })
            .setdatatype({ ECFDataType::FLOAT }));
    alpham = 2;
    REGISTER(alpham, ECFMetaData()
            .setdescription({ "SMALSE alpham" })
            .setdatatype({ ECFDataType::FLOAT }));
    power_precision = 1e-9;
    REGISTER(power_precision, ECFMetaData()
            .setdescription({ "SMALSE precision for power method" })
            .setdatatype({ ECFDataType::FLOAT }));
    power_maxit = 100;
    REGISTER(power_maxit, ECFMetaData()
            .setdescription({ "SMALSE max iteration for power method" })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));
    delta = 0.25;
    REGISTER(delta, ECFMetaData()
            .setdescription({ "SMALSE delta" })
            .setdatatype({ ECFDataType::FLOAT }));

    rtol = 1e-2;
    REGISTER(rtol, ECFMetaData()
            .setdescription({ "SMALSE rtol" })
            .setdatatype({ ECFDataType::FLOAT }));

    halfstep = true;
    REGISTER(halfstep, ECFMetaData()
            .setdescription({ "SMALSE halfstep" })
            .setdatatype({ ECFDataType::BOOL }));
    exp_projgrad = true;
    REGISTER(exp_projgrad, ECFMetaData()
            .setdescription({ "SMALSE exp_projgrad" })
            .setdatatype({ ECFDataType::BOOL }));
    prop_projgrad = true;
    REGISTER(prop_projgrad, ECFMetaData()
            .setdescription({ "SMALSE prop_projgrad" })
            .setdatatype({ ECFDataType::BOOL }));
    proj_grad = false;
    REGISTER(proj_grad, ECFMetaData()
            .setdescription({ "SMALSE proj_grad" })
            .setdatatype({ ECFDataType::BOOL }));
    gradproj = true;
    REGISTER(gradproj, ECFMetaData()
            .setdescription({ "SMALSE gradproj" })
            .setdatatype({ ECFDataType::BOOL }));
    optimset = true;
    REGISTER(optimset, ECFMetaData()
            .setdescription({ "SMALSE optimset" })
            .setdatatype({ ECFDataType::BOOL }));

    th = 0;
    REGISTER(th, ECFMetaData()
            .setdescription({ "SMALSE th" })
            .setdatatype({ ECFDataType::INTEGER }));
    no_enlarg_exp = 0;
    REGISTER(no_enlarg_exp, ECFMetaData()
            .setdescription({ "SMALSE no_enlarg_exp" })
            .setdatatype({ ECFDataType::INTEGER }));
    no_enlarg_prop = 0;
    REGISTER(th, ECFMetaData()
            .setdescription({ "SMALSE no_enlarg_prop" })
            .setdatatype({ ECFDataType::INTEGER }));

    sc_size = 200;
    REGISTER(sc_size, ECFMetaData()
            .setdescription({ "The size of null pivots for analytics regularization" })
            .setdatatype({ ECFDataType::POSITIVE_INTEGER }));

    REGISTER(auto_optimization, ECFMetaData().setdescription({ "Automatic optimization of FETI solver parameters." }));

    REGISTER(dual_operator_gpu_config, ECFMetaData()
        .setdescription({ "Dual operator on GPU configuration." })
        .setcollapsed());

    REGISTER(dualop_totalfeti_explicit_generalschur_cpu_config, ECFMetaData()
        .setdescription({ "Dual operator total feti explicit general schur on CPU configuration." })
        .setcollapsed());

    REGISTER(dualop_totalfeti_explicit_generalschur_gpu_config, ECFMetaData()
        .setdescription({ "Dual operator total feti explicit general schur on GPU configuration." })
        .setcollapsed());

    REGISTER(dualop_totalfeti_implicit_generalsparsesolver_cpu_config, ECFMetaData()
        .setdescription({ "Dual operator total feti implicit general sparse solver on CPU configuration." })
        .setcollapsed());

    REGISTER(dualop_hybridfeti_explicit_generalschur_cpu_config, ECFMetaData()
        .setdescription({ "Dual operator hybrid feti explicit general schur on CPU configuration." })
        .setcollapsed());

    REGISTER(dirichlet_generalschur_config, ECFMetaData()
        .setdescription({ "Dirichlet preconnditioner using general schur configuration." })
        .setcollapsed());
}



