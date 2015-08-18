#ifdef WIN32	 
	#include "stdafx.h"
#endif

#include <omp.h>
#include "mpi.h"
#include "mkl.h"

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
using std::cout;
using std::map; 
using std::make_pair; 

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "SparseMatrix.h"
#include "FEM_Assembler.h"
#include "SparseSolver.h"
#include "TimeEval.h"
#include "Cluster.h"

#include "utils.h"

#pragma once


class IterSolver
{
public:
	
	// *** Variables *************

	// MPI variables 
	int mpi_rank;
	int mpi_root;
	int mpi_size; 

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
	SparseSolver	GGt; 
	int				GGtsize; 
	
	// *** Setup variables 
	int USE_DYNAMIC; 
	int USE_KINV;
	int USE_GGtINV; 
	int USE_HFETI; 

	int USE_PREC;
	int USE_PIPECG; 
	
	int FIND_SOLUTION; 
	
	int CG_max_iter; 
	int NumberOfTimeIterations; 

	double epsilon; // stop condition  
	
	
	// DYNAMIC variables 
	double const_beta;
	double const_gama;
	double const_deltat; 
	

	// Timing objects 
	
	// Main timing object for main CG loop
	TimeEval timing; //("Main CG loop timing "); 
	TimeEval preproc_timing; // ("Preprocessing timing ");


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
	IterSolver(void);

	//Destructor
	~IterSolver(void);


	// *** Coarse problem related members  
	void CreateGGt    ( Cluster & cluster ); //, int mpi_rank, int mpi_root, int mpi_size, SparseSolver & GGt ); 
	void CreateGGt_inv( Cluster & cluster ); //, int mpi_rank, int mpi_root, int mpi_size, SparseSolver & GGt ) 
	void CreateGGt_inv_dist( Cluster & cluster ); 

	// *** Projectors 
	void Projector_l_compG    ( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, int output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 ); // int mpi_rank, SparseSolver & GGt, 
	void Projector_l_inv_compG( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, int output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 );

	// *** Apply A embers 
	void apply_A_l_compB     ( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out); 
	void apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

	// *** Preconditioner members
	void apply_prec_compB( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );
	void apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );

	// *** Public functions
	void Setup          ( SEQ_VECTOR <double> & parameters , Cluster & cluster_in ); 
	void Preprocessing  ( Cluster & cluster ); 
	void Solve_singular ( Cluster & cluster, string & result_file ); 

	
	// *** CG solvers 
	void Solve_RegCG_singular      ( Cluster & cluster ); //, vector <double> & x_l); // dual_soultion_in = x_l
	void Solve_RegCG_singular_dom  ( Cluster & cluster ); //, vector <double> & x_l); // dual_soultion_in = x_l


	void Solve_PipeCG_singular     ( Cluster & cluster ); //, vector <double> & x_l); // dual_soultion_in = x_l
	void Solve_PipeCG_singular_dom ( Cluster & cluster ); //, vector <double> & x_l); // dual_soultion_in = x_l


	// *** Functions related to getting solution from the solver 
	void GetSolution_Dual_singular_parallel ( Cluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out );
	void GetResiduum_Dual_singular_parallel ( Cluster & cluster, SEQ_VECTOR <double> & dual_residuum_out );
	
	void MakeSolution_Primal_singular_parallel ( Cluster & cluster);
	void GetSolution_Primal_singular_parallel ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out );
	
	void Save_to_Ensight_file (Cluster & cluster, string & result_file); 
	void Save_to_Ensight_file (Cluster & cluster, string & result_file, SEQ_VECTOR < SEQ_VECTOR < double> > & in_primal_solution_parallel);

	// *** Dynamic solvers 
	void Solve_RegCG_nonsingular  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel); 
	void Solve_PipeCG_nonsingular ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal, SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel); 
	
	// *** 
	//void Solve_Dynamic ( Cluster & cluster, string result_file );
	void Solve_Dynamic ( Cluster & cluster, string result_file, SEQ_VECTOR < SEQ_VECTOR < SEQ_VECTOR <double> > > & prim_solution);


};

//Utilities 

void SendMatrix(int rank, int source_rank, SparseMatrix & A_in, int dest_rank, SparseMatrix & B_out);

void SendMatrix_2(int source_rank, SparseMatrix & A_in, int dest_rank, SparseMatrix & B_out);
void RecvMatrix   ( SparseMatrix & B_out, int source_rank); 
void SendMatrix   ( SparseMatrix & A_in, int dest_rank );

void BcastMatrix(int rank, int mpi_root, int source_rank, SparseMatrix & A);

void All_Reduce_lambdas_compB( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ); 

void All_Reduce_lambdas      ( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ); // POZOR - musi jit pryc 

void compress_lambda_vector(Cluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda);

void decompress_lambda_vector(Cluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda); 

double parallel_norm_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector ); 

double parallel_ddot_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 );

void parallel_ddot_compressed_non_blocking( Cluster & cluster, 
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b, 
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b, 
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf) ;
