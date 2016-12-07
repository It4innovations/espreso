
#ifndef SOLVER_SPECIFIC_ITERSOLVER_H_
#define SOLVER_SPECIFIC_ITERSOLVER_H_

#include <omp.h>
#include "mpi.h"
// #include "mkl.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>

using std::vector;
using std::map;
using std::make_pair;

//#include <cilk/cilk.h>
//#include <cilk/cilk_api.h>

#include "../generic/SparseMatrix.h"
#include "sparsesolvers.h"
#include "clusters.h"

#include "../generic/utils.h"

#include "esbasis.h"

namespace espreso {

class IterSolverBase
{
public:

	// *** Variables *************

	ESPRESOSolver configuration;

	// MPI variables
	int  mpi_rank;
	int  mpi_root;
	int  mpi_size;

	// *** Main cluster object associated with iteration solver
	// Cluster & cluster;

	// *** solver variables
	SEQ_VECTOR <double> dual_soultion_decompressed_parallel;
	SEQ_VECTOR <double> dual_soultion_compressed_parallel;
	SEQ_VECTOR <double> dual_residuum_compressed_parallel;

	SEQ_VECTOR < SEQ_VECTOR <double> > primal_solution_parallel;
	SEQ_VECTOR <double> amplitudes;

	// Coarse problem variables
	SparseMatrix	GGt_Mat;
	SparseSolverCPU	GGt;
	eslocal 				GGtsize;

	// *** Setup variables
	eslocal  USE_DYNAMIC;
	eslocal  USE_KINV;
	eslocal  USE_GGtINV;
	eslocal  USE_HFETI;

	ESPRESO_PRECONDITIONER  USE_PREC;

	eslocal  CG_max_iter;

	eslocal PAR_NUM_THREADS;
	eslocal SOLVER_NUM_THREADS;

	double epsilon; // stop condition


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
	IterSolverBase(const ESPRESOSolver &configuration);

	//Destructor
	virtual ~IterSolverBase() {};


	// *** Coarse problem related members
	void CreateGGt    ( Cluster & cluster ); //, int mpi_rank, int mpi_root, int mpi_size, SparseSolverCPU & GGt );
	void CreateGGt_inv_dist  ( Cluster & cluster );
	void CreateGGt_inv_dist_d( Cluster & cluster );

	// *** Projectors
	void Projector_l_compG    ( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 ); // int mpi_rank, SparseSolverCPU & GGt,
	void Projector_l_inv_compG( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );
	void Projector_l_inv_compG_d( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal  output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );

	// *** Apply A embers - moved to children
  virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) =0;

  // *** Apply preconditioner 
  virtual void apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );

	// *** Public functions
	void Setup          ( SEQ_VECTOR <double> & parameters , Cluster & cluster_in );
	void Preprocessing  ( Cluster & cluster );

	void Solve_singular     ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel );
	void Solve_non_singular ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel );



	// *** Power Method - Estimation of maximum eigenvalue of matrix
	double Solve_power_method ( Cluster & cluster, double tol, eslocal maxit, eslocal method);

	//  *** Projected gradient and its components
	void proj_gradient ( SEQ_VECTOR <double> & x,SEQ_VECTOR <double> & g, SEQ_VECTOR <double> & lb,
			double alpha, double prec, SEQ_VECTOR <double> & g_til, SEQ_VECTOR <double> & fi_til, SEQ_VECTOR <double> & beta_til,
			SEQ_VECTOR <bool> & free );



	void Solve_QPCE_singular_dom  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );

	// *** CG solvers
	void Solve_RegCG_singular_dom  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	void Solve_full_ortho_CG_singular_dom ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	void Solve_GMRES_singular_dom ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	void Solve_BICGSTAB_singular_dom ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal ); 
	void Solve_new_CG_singular_dom ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );
	void Solve_PipeCG_singular_dom ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal );

	// *** Dynamic solvers
	void Solve_RegCG_nonsingular  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel);
	void Solve_PipeCG_nonsingular ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel);


	// *** Functions related to getting solution from the solver
	void GetSolution_Dual_singular_parallel ( Cluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out );
	void GetResiduum_Dual_singular_parallel ( Cluster & cluster, SEQ_VECTOR <double> & dual_residuum_out );

	void MakeSolution_Primal_singular_parallel ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out );
	void GetSolution_Primal_singular_parallel  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out );

};

//Utilities

void SendMatrix  (eslocal  rank, eslocal  source_rank, SparseMatrix & A_in, eslocal  dest_rank, SparseMatrix & B_out);

void SendMatrix2 (eslocal  rank, eslocal  source_rank, SparseMatrix & A_in, eslocal  dest_rank, SparseMatrix & B_out);

void RecvMatrix   ( SparseMatrix & B_out, eslocal  source_rank);
void SendMatrix   ( SparseMatrix & A_in, eslocal  dest_rank );

void ExchangeMatrices (SparseMatrix & A_in, SEQ_VECTOR <SparseMatrix> & B_out, SEQ_VECTOR <eslocal> neighbor_ranks );

void BcastMatrix(eslocal  rank, eslocal  mpi_root, eslocal  source_rank, SparseMatrix & A);

void All_Reduce_lambdas_compB( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );
void All_Reduce_lambdas_compB2( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );

void All_Reduce_lambdas      ( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ); // POZOR - musi jit pryc

void compress_lambda_vector(Cluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda);

void decompress_lambda_vector(Cluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda);

double parallel_norm_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector );

double parallel_ddot_compressed_double( Cluster & cluster, double *input_vector1, double *input_vector2 );
double parallel_ddot_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 );

void parallel_ddot_compressed_non_blocking( Cluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	SEQ_VECTOR<double> & input_norm_vec,
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf) ;

void parallel_ddot_compressed_non_blocking( Cluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf) ;

}

#endif /* SOLVER_SPECIFIC_ITERSOLVER_H_ */
