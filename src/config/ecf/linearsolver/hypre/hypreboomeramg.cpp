
#include "hypreboomeramg.h"

#include "config/configuration.hpp"

using namespace espreso;

HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration::HYPREBoomerAMGCycleRelaxTypeConfiguration()
{
	relax_type = RELAX_TYPE::HGSF;
	REGISTER(relax_type, ECFMetaData()
			.setdescription({ "Defines the smoother to be used. It uses the given smoother on the fine grid, the up and the down cycle and sets the solver on the coarsest level to Gaussian elimination" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("Jacobi").setdescription("Jacobi"))
			.addoption(ECFOption().setname("GSS").setdescription("Gauss-Seidel, sequential (very slow!)"))
			.addoption(ECFOption().setname("GSIP").setdescription("Gauss-Seidel, interior points in parallel, boundary sequential (slow!)"))
			.addoption(ECFOption().setname("HGSF").setdescription("hybrid Gauss-Seidel or SOR, forward solve"))
			.addoption(ECFOption().setname("HGSB").setdescription("hybrid Gauss-Seidel or SOR, backward solve"))
			.addoption(ECFOption().setname("HCHGS").setdescription("hybrid chaotic Gauss-Seidel (works only with OpenMP)"))
			.addoption(ECFOption().setname("HSGS").setdescription("hybrid symmetric Gauss-Seidel or SSOR"))
			.addoption(ECFOption().setname("LSHSGS").setdescription("l1-scaled hybrid symmetric Gauss-Seidel"))
			.addoption(ECFOption().setname("GE").setdescription("Gaussian elimination (only on coarsest level)"))
			.addoption(ECFOption().setname("CG").setdescription("CG (warning - not a fixed smoother - may require FGMRES)"))
			.addoption(ECFOption().setname("Chebyshev").setdescription("Chebyshev"))
			.addoption(ECFOption().setname("FCF").setdescription("FCF-Jacobi"))
			.addoption(ECFOption().setname("LSJ").setdescription("l1-scaled jacobi")));

	relax_type_cycle = RELAX_TYPE_CYCLE::COAREST;
	REGISTER(relax_type_cycle, ECFMetaData()
			.setdescription({ "Defines the smoother at a given cycle" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DOWN").setdescription("the down cycle"))
			.addoption(ECFOption().setname("UP").setdescription("the up cycle"))
			.addoption(ECFOption().setname("COAREST").setdescription("the coarsest level")));

}

HYPREBoomerAMGConfiguration::HYPREBoomerAMGCycleSweepsConfiguration::HYPREBoomerAMGCycleSweepsConfiguration()
{

	sweeps_num_specific = 1;
	REGISTER(sweeps_num_specific, ECFMetaData()
			.setdescription({ "Sets the number of sweeps at a specified cycle" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	relax_type_cycle = RELAX_TYPE_CYCLE::COAREST;
	REGISTER(relax_type_cycle, ECFMetaData()
			.setdescription({ "Defines the smoother at a given cycle" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("DOWN").setdescription("the down cycle"))
			.addoption(ECFOption().setname("UP").setdescription("the up cycle"))
			.addoption(ECFOption().setname("COAREST").setdescription("the coarsest level")));

}

HYPREBoomerAMGConfiguration::HYPREBoomerAMGConfiguration()
{
	convergence_tolerance = 1e-8;
	REGISTER(convergence_tolerance, ECFMetaData()
			.setdescription({ "Set the convergence tolerance" })
			.setdatatype({ ECFDataType::FLOAT }));

	min_iterations = 1;
	REGISTER(min_iterations, ECFMetaData()
			.setdescription({ "Sets minimum number of iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	max_iterations = 20;
	REGISTER(max_iterations, ECFMetaData()
			.setdescription({ "Sets maximum number of iterations" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	max_coarest_grid_size = 9;
	REGISTER(max_coarest_grid_size, ECFMetaData()
			.setdescription({ "Sets maximum size of coarsest grid" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	min_coarest_grid_size = 1;
	REGISTER(min_coarest_grid_size, ECFMetaData()
			.setdescription({ "Sets minimum size of coarsest grid" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	max_multigrid_levels = 25;
	REGISTER(max_multigrid_levels, ECFMetaData()
			.setdescription({ "Sets maximum number of multigrid levels" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	amg_strength_treshold = 0.25;
	REGISTER(amg_strength_treshold, ECFMetaData()
			.setdescription({ "Sets AMG strength threshold. For 2d Laplace operators, 0.25 is a good value, for 3d Laplace operators, 0.5 or 0.6 is a better value. For elasticity problems, a large strength threshold, such as 0.9, is often better" })
			.setdatatype({ ECFDataType::FLOAT }));

	coarsening_type = COARSENING_TYPE::Falgout;
	REGISTER(coarsening_type, ECFMetaData()
			.setdescription({ "Defines which parallel coarsening algorithm is used" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CLJP").setdescription("CLJP-coarsening (a parallel coarsening algorithm using independent sets"))
			.addoption(ECFOption().setname("Ruge_Stuben").setdescription("classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse points on the boundaries"))
			.addoption(ECFOption().setname("Falgout").setdescription("Falgout coarsening (uses 1 first, followed by CLJP using the interior coarse points generated by 1 as its first independent set)"))
			.addoption(ECFOption().setname("PMIS").setdescription("PMIS-coarsening (a parallel coarsening algorithm using independent sets, generating lower complexities than CLJP, might also lead to slower convergence)"))
			.addoption(ECFOption().setname("HMIS").setdescription("HMIS-coarsening (uses one pass Ruge-Stueben on each processor independently, followed by PMIS using the interior C-points generated as its first independent set)"))
			.addoption(ECFOption().setname("CGC").setdescription("CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer"))
			.addoption(ECFOption().setname("CGC_E").setdescription("CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer")));

	interpolation_type = INTERPOLATION_TYPE::CLASSIC_MODIFF;
	REGISTER(interpolation_type, ECFMetaData()
			.setdescription({ "Defines which parallel interpolation operator is used" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("CLASSIC_MODIFF").setdescription("classical modified interpolation"))
			.addoption(ECFOption().setname("LS").setdescription("LS interpolation (for use with GSMG)"))
			.addoption(ECFOption().setname("CLASSIC_HYPERBOLIC").setdescription("classical modified interpolation for hyperbolic PDEs"))
			.addoption(ECFOption().setname("DIRECT").setdescription("direct interpolation (with separation of weights)"))
			.addoption(ECFOption().setname("MULTLIPASS").setdescription("multipass interpolation"))
			.addoption(ECFOption().setname("MULTIPASSS_SEPARATION").setdescription("multipass interpolation (with separation of weights)"))
			.addoption(ECFOption().setname("EXTENDED_I").setdescription("extended+i interpolation"))
			.addoption(ECFOption().setname("EXTENDED_I_NO_NEIGHBOR").setdescription("extended+i (if no common C neighbor) interpolation"))
			.addoption(ECFOption().setname("STANDARD").setdescription("standard interpolation"))
			.addoption(ECFOption().setname("STANDARD_SEPPARATION").setdescription("standard interpolation (with separation of weights)"))
			.addoption(ECFOption().setname("CLASSIC_BLOCK").setdescription("classical block interpolation (for use with nodal systems version only)"))
			.addoption(ECFOption().setname("CLASSIC_BLOCK_DIAG").setdescription("classical block interpolation (for use with nodal systems version only) with diagonalized diagonal blocks"))
			.addoption(ECFOption().setname("FF").setdescription("FF interpolation"))
			.addoption(ECFOption().setname("FF1").setdescription("FF1 interpolation"))
			.addoption(ECFOption().setname("EXTENDED").setdescription("extended interpolation")));

	interp_trunc_factor = 0.0;
	REGISTER(interp_trunc_factor, ECFMetaData()
			.setdescription({ "Defines a truncation factor for the interpolation" })
			.setdatatype({ ECFDataType::FLOAT }));

	max_element_per_row = 0;
	REGISTER(max_element_per_row, ECFMetaData()
			.setdescription({ "Defines the maximal number of elements per row for the interpolation" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	weight_separation = 0;
	REGISTER(weight_separation, ECFMetaData()
			.setdescription({ "Defines whether separation of weights is used when defining strength for standard interpolation or multipass interpolation" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	cycle_type = CYCLE_TYPE::V_CYCLE;
	REGISTER(cycle_type, ECFMetaData()
			.setdescription({ "Defines the type of cycle" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("V_CYCLE").setdescription("V-cycle"))
			.addoption(ECFOption().setname("W_CYCLE").setdescription("W-cycle")));

	sweeps_num = 1;
	REGISTER(sweeps_num, ECFMetaData()
			.setdescription({ "Sets the number of sweeps" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	relax_type = RELAX_TYPE::HGSF;
	REGISTER(relax_type, ECFMetaData()
			.setdescription({ "Defines the smoother to be used. It uses the given smoother on the fine grid, the up and the down cycle and sets the solver on the coarsest level to Gaussian elimination" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("Jacobi").setdescription("Jacobi"))
			.addoption(ECFOption().setname("GSS").setdescription("Gauss-Seidel, sequential (very slow!)"))
			.addoption(ECFOption().setname("GSIP").setdescription("Gauss-Seidel, interior points in parallel, boundary sequential (slow!)"))
			.addoption(ECFOption().setname("HGSF").setdescription("hybrid Gauss-Seidel or SOR, forward solve"))
			.addoption(ECFOption().setname("HGSB").setdescription("hybrid Gauss-Seidel or SOR, backward solve"))
			.addoption(ECFOption().setname("HCHGS").setdescription("hybrid chaotic Gauss-Seidel (works only with OpenMP)"))
			.addoption(ECFOption().setname("HSGS").setdescription("hybrid symmetric Gauss-Seidel or SSOR"))
			.addoption(ECFOption().setname("LSHSGS").setdescription("l1-scaled hybrid symmetric Gauss-Seidel"))
			.addoption(ECFOption().setname("GE").setdescription("Gaussian elimination (only on coarsest level)"))
			.addoption(ECFOption().setname("CG").setdescription("CG (warning - not a fixed smoother - may require FGMRES)"))
			.addoption(ECFOption().setname("Chebyshev").setdescription("Chebyshev"))
			.addoption(ECFOption().setname("FCF").setdescription("FCF-Jacobi"))
			.addoption(ECFOption().setname("LSJ").setdescription("l1-scaled jacobi")));


	allow_cycle_relax_type = false;
	REGISTER(allow_cycle_relax_type, ECFMetaData()
			.setdescription({ "Allow smoother at a given cycle definition" })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(cycle_relax_type, ECFMetaData()
			.setdescription({ "Defines the smoother at a given cycle"}));

	relax_order = 1;
	REGISTER(relax_order, ECFMetaData()
			.setdescription({ "Defines in which order the points are relaxed" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	relax_weight = 1.0;
	REGISTER(relax_weight, ECFMetaData()
			.setdescription({ "Defines the relaxation weight for smoothed Jacobi and hybrid SOR on all levels" })
			.setdatatype({ ECFDataType::FLOAT }));

	allow_cycle_num_sweeps = false;
	REGISTER(allow_cycle_num_sweeps, ECFMetaData()
			.setdescription({ "Allow smoother at a given cycle definition" })
			.setdatatype({ ECFDataType::BOOL }));	

	REGISTER(cycle_sweep_spec, ECFMetaData()
			.setdescription({ "Sets the number of sweeps at a specified cycle"}));


	level_relax_weight = 1.0;
	REGISTER(level_relax_weight, ECFMetaData()
			.setdescription({ "Defines the relaxation weight for smoothed Jacobi and hybrid SOR on the user defined level" })
			.setdatatype({ ECFDataType::FLOAT }));

	out_relax_weight = 1.0;
	REGISTER(out_relax_weight, ECFMetaData()
			.setdescription({ "Defines the outer relaxation weight for hybrid SOR and SSOR on all levels" })
			.setdatatype({ ECFDataType::FLOAT }));

	outer_weight_level = 1.0;
	REGISTER(outer_weight_level, ECFMetaData()
			.setdescription({ "Defines the outer relaxation weight for hybrid SOR or SSOR on the user defined level" })
			.setdatatype({ ECFDataType::FLOAT }));

	cheby_order = 2;
	REGISTER(cheby_order, ECFMetaData()
			.setdescription({ "Defines the Order for Chebyshev smoother.(valid options are 1-4)" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	cheby_fraction = 0.3;
	REGISTER(cheby_fraction, ECFMetaData()
			.setdescription({ "Fraction of the spectrum to use for the Chebyshev smoother" })
			.setdatatype({ ECFDataType::FLOAT }));

	smooth_type = SMOOTH_TYPE::SCHWARZ;
	REGISTER(smooth_type, ECFMetaData()
			.setdescription({ "Enables the use of more complex smoothers" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("SCHWARZ").setdescription("Schwarz smoothers"))
			.addoption(ECFOption().setname("PILUT").setdescription("Pilut"))
			.addoption(ECFOption().setname("PARASAILS").setdescription("ParaSails"))
			.addoption(ECFOption().setname("EUCLID").setdescription("Euclid")));

	smooth_level_num = 0;
	REGISTER(smooth_level_num, ECFMetaData()
			.setdescription({ "Sets the number of levels for more complex smoothers" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	smooth_sweep_num = 1;
	REGISTER(smooth_sweep_num, ECFMetaData()
			.setdescription({ "Sets the number of sweeps for more complex smoothers" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));			

	schwarz_variant = 0;
	REGISTER(schwarz_variant, ECFMetaData()
			.setdescription({ "Defines which variant of the Schwarz method is used" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	schwarz_overlap = 1;
	REGISTER(schwarz_overlap, ECFMetaData()
			.setdescription({ "Defines the overlap for the Schwarz method" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	schwarz_doamin = 2;
	REGISTER(schwarz_doamin, ECFMetaData()
			.setdescription({ "Defines the type of domain used for the Schwarz method" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	schwarz_use_nonsym = 0;
	REGISTER(schwarz_use_nonsym, ECFMetaData()
			.setdescription({ "Indicates that the aggregates may not be SPD for the Schwarz method" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	interp_vec_variant = INTERP_VEC_VARIANT::NONE;
	REGISTER(interp_vec_variant, ECFMetaData()
			.setdescription({ "Defines the interpolation variant" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Omit interpolation"))
			.addoption(ECFOption().setname("GM1").setdescription("GM approach 1"))
			.addoption(ECFOption().setname("GM2").setdescription("GM approach 2 (to be preferred over 1)"))
			.addoption(ECFOption().setname("LM").setdescription("LN approach")));

	measure_type = MEASURE_TYPE::LOCAL;
	REGISTER(measure_type, ECFMetaData()
			.setdescription({ "Defines whether local or global measures are used" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("LOCAL").setdescription("Local measures"))
			.addoption(ECFOption().setname("GLOBAL").setdescription("Global measures")));


	aggressive_coarsening_interp = AGGRESSIVE_COARSENING_INTERP::MULTIPASS;
	REGISTER(aggressive_coarsening_interp, ECFMetaData()
			.setdescription({ "Defines the interpolation used on levels of aggressive coarsening" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("TWO_STAGE_IEXTENDED").setdescription("2-stage extended+i interpolation"))
			.addoption(ECFOption().setname("TWO_STAGE_STANDARD").setdescription("2-stage standard interpolation"))
			.addoption(ECFOption().setname("TWO_STAGE_EXTENDED").setdescription("2-stage extended interpolation"))
			.addoption(ECFOption().setname("MULTIPASS").setdescription("multipass interpolation")));

	aggressive_coarsening_trunc_factor = 0.0;
	REGISTER(aggressive_coarsening_trunc_factor, ECFMetaData()
			.setdescription({ "Defines the truncation factor for the interpolation used for aggressive coarsening" })
			.setdatatype({ ECFDataType::FLOAT }));

	aggressive_coarsening_p12_trunc_fact = 0.0;
	REGISTER(aggressive_coarsening_p12_trunc_fact, ECFMetaData()
			.setdescription({ "Defines the truncation factor for the matrices P1 and P2 which are used to build 2-stage interpolation" })
			.setdatatype({ ECFDataType::FLOAT }));

	aggressive_coarsening_p_max_elements = 0;
	REGISTER(aggressive_coarsening_p_max_elements, ECFMetaData()
			.setdescription({ "Defines the maximal number of elements per row for the interpolation used for aggressive coarsening" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));		

	aggressive_coarsening_p12_max_elements = 0;
	REGISTER(aggressive_coarsening_p12_max_elements, ECFMetaData()
			.setdescription({ "Defines the maximal number of elements per row for the matrices P1 and P2 which are used to build 2-stage interpolation" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));						

	comm_strength_treshold = 1.0;
	REGISTER(comm_strength_treshold, ECFMetaData()
			.setdescription({ "Defines the largest strength threshold for which the strength matrix S uses the communication package of the operator A" })
			.setdatatype({ ECFDataType::FLOAT }));

	diag_dominant_strength = 0.9;
	REGISTER(diag_dominant_strength, ECFMetaData()
			.setdescription({ "Sets a parameter to modify the definition of strength for diagonal dominant portions of the matrix" })
			.setdatatype({ ECFDataType::FLOAT }));

	non_galerkin_drop_tol = 0.0;
	REGISTER(non_galerkin_drop_tol, ECFMetaData()
			.setdescription({ "Defines the non-Galerkin drop-tolerance for sparsifying coarse grid operators and thus re- ducing communication. Value specified here is set on all levels" })
			.setdatatype({ ECFDataType::FLOAT }));

	level_non_galerkin_drop_tol = 0;
	REGISTER(level_non_galerkin_drop_tol, ECFMetaData()
			.setdescription({ "Defines the level specific non-Galerkin drop-tolerances for sparsifying coarse grid operators and thus reducing communication. A drop-tolerance of 0.0 means to skip doing non-Galerkin on that level" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));	

	aggressive_coarsening_levels = 0;
	REGISTER(aggressive_coarsening_levels, ECFMetaData()
			.setdescription({ "Defines the number of levels of aggressive coarsening" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));	

	aggressive_coarsening_degree = 1;
	REGISTER(aggressive_coarsening_degree, ECFMetaData()
			.setdescription({ "Defines the degree of aggressive coarsening" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));	

	cgc_pathes = 0;
	REGISTER(cgc_pathes, ECFMetaData()
			.setdescription({ "Defines the number of pathes for CGC-coarsening" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));	

	nodal_sys_coarsening = NODAL_SYS_COARSENING::UNKNOWN;
	REGISTER(nodal_sys_coarsening, ECFMetaData()
			.setdescription({ "Print solver info" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("UNKNOWN").setdescription("unknown-based coarsening, only coarsens within same function"))
			.addoption(ECFOption().setname("FROBENIUSNORM").setdescription("Frobenius norm"))
			.addoption(ECFOption().setname("SUM_ABSOLUTE_BLOCK").setdescription("sum of absolute values of elements in each block"))
			.addoption(ECFOption().setname("LARGEST_ELEM_BLOCK").setdescription("largest element in each block (not absolute value)"))
			.addoption(ECFOption().setname("SUM_ROW").setdescription("row-sum norm"))
			.addoption(ECFOption().setname("ALL_BLOC_SUM").setdescription("sum of all values in each block")));

	nodal_diag_treatment = 0;
	REGISTER(nodal_diag_treatment, ECFMetaData()
			.setdescription({ "Sets whether to give special treatment to diagonal elements in the nodal systems version" })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));	

	old_default=false;
	REGISTER(old_default, ECFMetaData()
			.setdescription({ "Recovers old default for coarsening and interpolation, ie Falgout coarsening and untruncated modified classical interpolation. This option might be preferred for 2 dimensional problems" })
			.setdatatype({ ECFDataType::BOOL }));



	solver_info = SOLVER_INFO::NO_INFO;
	REGISTER(solver_info, ECFMetaData()
			.setdescription({ "Print solver info" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NO_INFO").setdescription("no printout"))
			.addoption(ECFOption().setname("SETUP_INFO").setdescription("print setup information"))
			.addoption(ECFOption().setname("SOLVE_INFO").setdescription("print solve information"))
			.addoption(ECFOption().setname("SETUP_SOLVE_INFO").setdescription("print both setup and solve information")));

}
