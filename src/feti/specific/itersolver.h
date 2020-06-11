
#ifndef SOLVER_SPECIFIC_ITERSOLVER_H_
#define SOLVER_SPECIFIC_ITERSOLVER_H_

//#include <omp.h>
//#include "mpi.h"
//// #include "mkl.h"
//
//#include <string>
//#include <sstream>
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <algorithm>
//#include <math.h>
//#include <iomanip>
//#include <map>
//
//using std::vector;
//using std::map;
//using std::make_pair;

//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>

#include "feti/generic/SparseMatrix.h"
#include "sparsesolvers.h"
#include "clusters.h"
#include "superclusters.h"
#include "feti/generic/utils.h"

namespace espreso {

//class SuperCluster;

class IterSolverBase
{
public:

	// *** Variables *************

	FETIConfiguration &configuration;

	// MPI variables
	int  mpi_rank;
	int  mpi_root;
	int  mpi_size;

	//TODO delete
	int  numClusters;
	SEQ_VECTOR <Cluster*> *clusters;

	// *** solver variables
	SEQ_VECTOR <double> dual_soultion_decompressed_parallel;
	SEQ_VECTOR <double> dual_soultion_compressed_parallel;
	SEQ_VECTOR <double> dual_residuum_compressed_parallel;

	SEQ_VECTOR < SEQ_VECTOR <double> > primal_solution_parallel;
	SEQ_VECTOR <double> amplitudes;

	// Coarse problem variables
	SparseMatrix	GGt_Mat;
	SparseSolverCPU	GGt;
	esint 		GGtsize;

	// *** Setup variables
	esint  USE_KINV;
	esint  USE_GGtINV;
	esint  USE_HFETI;

	FETIConfiguration::PRECONDITIONER  USE_PREC;

	esint  CG_max_iter;

	esint PAR_NUM_THREADS;
	esint SOLVER_NUM_THREADS;

	double precision; // stop condition


	// Timing objects

	// Main timing object for main CG loop
	TimeEval timing; //("Main CG loop timing ");
	TimeEval preproc_timing; // ("Preprocessing timing ");
	TimeEval postproc_timing;

	TimeEval timeEvalAppa; // (string("Apply Kplus timing "));
	TimeEvent apa_B1t; //	  (string("x = B1t * lambda "));
	TimeEvent apa_kplus; //	  (string("multKplus(local or global) "));
	TimeEvent apa_B1; //	  (string("lambda = B1 * x "));
	TimeEvent apa_allred; //  (string("All_Reduce_lambdas "));
	//timeEvalAppa.AddEvent(apa_decomp);
	//timeEvalAppa.AddEvent(apa_comp);
	//timeEvalAppa.AddEvent(apa_B1t);
	//timeEvalAppa.AddEvent(apa_kplus);
	//timeEvalAppa.AddEvent(apa_B1);
	//timeEvalAppa.AddEvent(apa_allred);

	TimeEval  timeEvalPrec;		// (string("Apply Precond. timing "));
	TimeEvent prec_kplus;		//  (string("B1 * P * B1t "));
	TimeEvent prec_allred;		// (string("All_Reduce_lambdas "));


	TimeEval timeEvalProj; // (string("Projector timing "));
	TimeEvent proj_G1t; //	  (string("x = G1 * lambda "));
	TimeEvent proj_Gthr; //	  (string("MPI_gather - collective "));
	TimeEvent proj_GGt; //	  (string("GGt Solve on master node "));
	TimeEvent proj_Sctr; //	  (string("MPI_Scatter - collective "));
	TimeEvent proj_Gx; //	  (string("lambda = G1t * x "));
	TimeEvent proj_allred; // (string("All_Reduce_lambdas "));
	//timeEvalProj.AddEvent(proj_decomp);
	//timeEvalProj.AddEvent(proj_comp);
	//timeEvalProj.AddEvent(proj_G1t);
	//timeEvalProj.AddEvent(proj_Gthr);
	//timeEvalProj.AddEvent(proj_GGt);
	//timeEvalProj.AddEvent(proj_Sctr);
	//timeEvalProj.AddEvent(proj_Gx);
	//timeEvalProj.AddEvent(proj_allred);

	TimeEvent ddot_time; // (string("Parallel DDOT - alpha and gamma"));
	TimeEvent proj_time; // (string("Projector_l "));
	TimeEvent appA_time; // (string("ApplyA_l "));
	TimeEvent vec_time; //  (string("vector processing in CG "));
	TimeEvent norm_time; // (string("parallel DDOT - norm "));

	TimeEvent proj1_time; // (string("Projector_l - before PREC "));
	TimeEvent proj2_time; // (string("Projector_l - after PREC "));
	TimeEvent prec_time; //  (string("Preconditioner "));
	TimeEvent ddot_alpha; // (string("2x ddot for Alpha "));
	TimeEvent ddot_beta; //  (string("2x ddot for Beta "));

	//preproc_timing.totalTime.AddEnd(omp_get_wtime());



	// *** Members ***************

	//Constructor
	IterSolverBase(FETIConfiguration &configuration);

	//Destructor
	virtual ~IterSolverBase() {};

	// *** Coarse problem related members
	void CreateGGt    ( SuperCluster & cluster ); //, int mpi_rank, int mpi_root, int mpi_size, SparseSolverCPU & GGt );
	void CreateGGt_Inv  ( SuperCluster & cluster );
	void CreateConjGGt_Inv  ( SuperCluster & cluster );
	void CreateGGt_Inv_old( SuperCluster & cluster );

	// *** Projectors
	void Projector    ( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 ); // int mpi_rank, SparseSolverCPU & GGt,
	void Projector_Inv( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );
	void Projector_Inv_old( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );

	void ConjProjector_Inv( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );
	void ConjProjector_Inv2( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );
	void ConjProjector_Inv3( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, esint  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );


	void CreateConjProjector(Cluster & cluster);
	void ConjProj(  Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);
	void ConjProj_t(Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

	void ConjProj_lambda0(  Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);



	// *** Apply A embers - moved to children
	virtual void apply_A_l_comp_dom_B        ( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) =0;
	virtual void apply_A_l_comp_dom_B_P      ( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) =0;
	virtual void apply_A_l_comp_dom_B_P_local( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) =0;
	virtual void apply_A_l_comp_dom_B_P_local_sparse( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<esint> & tmp_in_indices, SEQ_VECTOR<double> & tmp_in_values, SEQ_VECTOR<esint> & tmp_out_indices, SEQ_VECTOR<double> & tmp_out_values) =0;

	void apply_A_l_Mat		 ( TimeEval & time_eval, SuperCluster & cluster, SparseMatrix       & X_in, SparseMatrix       & Y_out) ;
	void apply_A_l_Mat_local ( TimeEval & time_eval, SuperCluster & cluster, SparseMatrix       & X_in, SparseMatrix       & Y_out) ;
	void apply_A_l_Mat_local_sparse ( TimeEval & time_eval, SuperCluster & cluster, SparseMatrix       & X_in, SparseMatrix       & Y_out) ;


	// *** Apply preconditioner
	virtual void Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );

	// *** Public functions
	void Preprocessing  ( SuperCluster & cluster );

	void Solve     ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel, SEQ_VECTOR < SEQ_VECTOR <double> > & out_dual_solution_parallel );


	// *** Power Method - Estimation of maximum eigenvalue of matrix
	double Solve_power_method ( SuperCluster & cluster, double tol, esint maxit, esint method);

	//  *** Projected gradient and its components
	void proj_gradient ( SEQ_VECTOR <double> & x,SEQ_VECTOR <double> & g, SEQ_VECTOR <double> & lb,
			double alpha, double prec, SEQ_VECTOR <double> & g_til, SEQ_VECTOR <double> & fi_til, SEQ_VECTOR <double> & beta_til,
			SEQ_VECTOR <bool> & free );



	int Solve_QPCE_singular_dom  ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );

	// *** CG solvers
	int Solve_RegCG  ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_RegCG_ConjProj  ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_full_ortho_CG_singular_dom ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_GMRES_singular_dom ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_GMRES_ConjProj ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_BICGSTAB_singular_dom ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_new_CG_singular_dom ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_PipeCG_singular_dom ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	int Solve_full_ortho_CG_singular_dom_geneo ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal);

//	// *** Dynamic solvers
//	void Solve_RegCG_nonsingular  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel);
//	void Solve_PipeCG_nonsingular ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel);


	// *** Functions related to getting solution from the solver
	void GetSolution_Dual_singular_parallel ( SuperCluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out );
	void GetResiduum_Dual_singular_parallel ( SuperCluster & cluster, SEQ_VECTOR <double> & dual_residuum_out );

	void MakeSolution_Primal_singular_parallel ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out );
	void GetSolution_Primal_singular_parallel  ( SuperCluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out, SEQ_VECTOR < SEQ_VECTOR <double> > & dual_solution_out, int iters);

};



//Utilities

void SendMatrix  (esint  rank, esint  source_rank, SparseMatrix & A_in, esint  dest_rank, SparseMatrix & B_out);

void SendMatrix2 (esint  rank, esint  source_rank, SparseMatrix & A_in, esint  dest_rank, SparseMatrix & B_out);

void RecvMatrix   ( SparseMatrix & B_out, esint  source_rank);
void SendMatrix   ( SparseMatrix & A_in, esint  dest_rank );

void ExchangeMatrices (SparseMatrix & A_in, SEQ_VECTOR <SparseMatrix> & B_out, SEQ_VECTOR <esint> neighbor_ranks );
void ExchangeMatrices2 (SEQ_VECTOR <SparseMatrix> & A_in, SEQ_VECTOR <SparseMatrix> & B_out, SEQ_VECTOR <esint> neighbor_ranks );

void ExchangeVector (SEQ_VECTOR <esint> & vec_in, SEQ_VECTOR <SEQ_VECTOR<esint>> & vec_out, SEQ_VECTOR <esint> neighbor_ranks );

void BcastMatrix(esint  rank, esint  mpi_root, esint  source_rank, SparseMatrix & A);

void All_Reduce_lambdas_compB( SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );
//void All_Reduce_lambdas_compB2( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );

void All_Reduce_lambdas      ( SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ); // POZOR - musi jit pryc

void compress_lambda_vector(SuperCluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda);

void decompress_lambda_vector(SuperCluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda);

double parallel_norm_compressed( SuperCluster & cluster, SEQ_VECTOR<double> & input_vector );

double parallel_ddot_compressed_double( SuperCluster & cluster, double *input_vector1, double *input_vector2 );
double parallel_ddot_compressed( SuperCluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 );

void parallel_ddot_compressed_non_blocking( SuperCluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	SEQ_VECTOR<double> & input_norm_vec,
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf) ;

void parallel_ddot_compressed_non_blocking( SuperCluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf) ;

}

#endif /* SOLVER_SPECIFIC_ITERSOLVER_H_ */
