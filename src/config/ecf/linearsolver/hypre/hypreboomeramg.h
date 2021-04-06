
#ifndef SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_
#define SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_

#include "config/description.h"

namespace espreso {

struct HYPREBoomerAMGConfiguration: public ECFDescription {

	double convergence_tolerance;
	int min_iterations, max_iterations, max_coarest_grid_size, min_coarest_grid_size;
	int max_multigrid_levels;

	enum class COARSENING_TYPE {
		CLJP,
		Ruge_Stuben,
		Falgout,
		PMIS,
		HMIS,
		CGC,
		CGC_E
	};
	COARSENING_TYPE coarsening_type;

	enum class INTERPOLATION_TYPE {
		CLASSIC_MODIFF,
		LS,
		CLASSIC_HYPERBOLIC,
		DIRECT,
		MULTLIPASS,
		MULTIPASSS_SEPARATION,
		EXTENDED_I,
		EXTENDED_I_NO_NEIGHBOR,
		STANDARD,
		STANDARD_SEPPARATION,
		CLASSIC_BLOCK,
		CLASSIC_BLOCK_DIAG,
		FF,
		FF1,
		EXTENDED
	};
	INTERPOLATION_TYPE interpolation_type;

	double interp_trunc_factor;
	int max_element_per_row, weight_separation;

	enum class CYCLE_TYPE {
		V_CYCLE,
		W_CYCLE
	};
	CYCLE_TYPE cycle_type;

	enum class RELAX_TYPE {
		Jacobi,
		GSS,
		GSIP,
		HGSF,
		HGSB,
		HCHGS,
		HSGS,
		LSHSGS,
		GE,
		CG,
		Chebyshev,
		FCF,
		LSJ
	};
	RELAX_TYPE relax_type;


	bool allow_cycle_relax_type;

	int sweeps_num;
	bool allow_cycle_num_sweeps;

	struct HYPREBoomerAMGCycleSweepsConfiguration: public ECFDescription {
		
		int sweeps_num_specific;

		enum class RELAX_TYPE_CYCLE {
			DOWN,
			UP,
			COAREST
		};

		RELAX_TYPE_CYCLE relax_type_cycle;
	
		HYPREBoomerAMGCycleSweepsConfiguration();
	};
	HYPREBoomerAMGCycleSweepsConfiguration cycle_sweep_spec;


	struct HYPREBoomerAMGCycleRelaxTypeConfiguration: public ECFDescription {
		
		RELAX_TYPE relax_type;

		enum class RELAX_TYPE_CYCLE {
			DOWN,
			UP,
			COAREST
		};
		RELAX_TYPE_CYCLE relax_type_cycle;
	
		HYPREBoomerAMGCycleRelaxTypeConfiguration();
	};
	HYPREBoomerAMGCycleRelaxTypeConfiguration cycle_relax_type;
	
	int relax_order;
	double relax_weight;
	
	double level_relax_weight;
	
	double out_relax_weight;
	double outer_weight_level;
	int cheby_order;
	double cheby_fraction;

	enum class SMOOTH_TYPE {
		SCHWARZ,
		PILUT,
		PARASAILS,
		EUCLID
	};
	SMOOTH_TYPE smooth_type;

	int smooth_level_num, smooth_sweep_num, schwarz_variant, schwarz_overlap, schwarz_doamin, schwarz_use_nonsym;

	enum class SOLVER_INFO {
		NO_INFO,
		SETUP_INFO,
		SOLVE_INFO,
		SETUP_SOLVE_INFO
	};
	SOLVER_INFO solver_info;

	enum class INTERP_VEC_VARIANT {
		NONE,
		GM1,
		GM2,
		LM
	};
	INTERP_VEC_VARIANT interp_vec_variant;

	bool old_default;

	double amg_strength_treshold, diag_dominant_strength, comm_strength_treshold, non_galerkin_drop_tol;
	int level_non_galerkin_drop_tol;

	enum class MEASURE_TYPE {
		LOCAL,
		GLOBAL
	};
	MEASURE_TYPE measure_type;

	enum class NODAL_SYS_COARSENING {
		UNKNOWN,
		FROBENIUSNORM,
		SUM_ABSOLUTE_BLOCK,
		LARGEST_ELEM_BLOCK,
		SUM_ROW,
		ALL_BLOC_SUM
	};
	NODAL_SYS_COARSENING nodal_sys_coarsening;

	int cgc_pathes, aggressive_coarsening_levels, aggressive_coarsening_degree;

	enum class AGGRESSIVE_COARSENING_INTERP {
		TWO_STAGE_IEXTENDED,
		TWO_STAGE_STANDARD,
		TWO_STAGE_EXTENDED,
		MULTIPASS
	};
	AGGRESSIVE_COARSENING_INTERP aggressive_coarsening_interp;

	int nodal_diag_treatment;
	double aggressive_coarsening_trunc_factor, aggressive_coarsening_p12_trunc_fact;
	int aggressive_coarsening_p_max_elements, aggressive_coarsening_p12_max_elements;

	HYPREBoomerAMGConfiguration();
};

}

#endif /* SRC_CONFIG_ECF_LINEARSOLVER_HYPRE_HYPREBOOMERAMG_H_ */
