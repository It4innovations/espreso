// Start - 32 bit short and 64 bit long 
//typedef long long		longInt;  //POZOR - XC
typedef int				ShortInt; 
#define MPI_ShortInt MPI_INT  
#define MPI_LongInt  MPI_LONG_LONG
// END - 16bit short and 32 but long 


// Start - 16bit short and 32 but long 
typedef int longInt;
//typedef short ShortInt; 
//#define MPI_ShortInt MPI_SHORT  
//#define MPI_LongInt  MPI_INT
// END - 16bit short and 32 but long 


// start - 16bit short and 64 bit long 
//typedef long long		longInt;
//typedef short			ShortInt; 
//#define MPI_ShortInt	MPI_SHORT  
//#define MPI_LongInt		MPI_LONG_LONG
// end - 16bit short and 64 bit long 

//// 32 bit 
//typedef int		longInt;
//typedef int 	ShortInt; 
//#define MPI_ShortInt	MPI_INT  
//#define MPI_LongInt		MPI_INT



#define MKL_CBWR
#define TIME_MEAS		1
#define DEBUG			0
#define OMP_PARALLEL	1

#ifdef WIN32	 
	#include "stdafx.h"
#endif

#ifndef WIN32	 
 #include "sys/types.h"
 #include "sys/sysinfo.h"
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

#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include <ctime>
#include <stack>
#include <time.h>
std::stack<clock_t> tictoc_stack;

using std::vector;
using std::cout;
using std::map; 
using std::make_pair; 

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "SparseMatrix.h"
#include "FEM_Assembler.h"
#include "SparseDSS.h"
#include "SparseSolver.h"
#include "TimeEval.h"
#include "Domain.h"
#include "Cluster.h"
#include "IterSolver.h"
#include "ConfigFile.h"

#include "utils.h"


// Test routines for debugging 
double vsumhex(SEQ_VECTOR <double> & vec) {
	double sum = 0; 

	for (int i = 0; i < vec.size(); i++) {
		sum = sum + vec[i];
		//printf(" %d - %llx \n", i, vec[i]); 
	}
	return sum;
}

double vsum(SEQ_VECTOR <double> & vec) {
	double sum = 0; 
	
	for (int i = 0; i < vec.size(); i++)
		sum = sum + vec[i];
	
	return sum;
}

double vsumv(double * vec, int size) {
	double sum = 0; 

	for (int i = 0; i < size; i++)
		sum = sum + vec[i];

	return sum;
}
// END - Test routines for debugging 

// *** Memory functions ************************************************************
void GetMemoryStat( ) 
{
#ifndef WIN32

	struct sysinfo memInfo;
	sysinfo (&memInfo);

	long long totalPhysMem;
	long long physMemUsed; 

	totalPhysMem = memInfo.totalram;
	//Multiply in next statement to avoid int overflow on right hand side...
	totalPhysMem *= memInfo.mem_unit;

	physMemUsed	= memInfo.totalram - memInfo.freeram;
	//Multiply in next statement to avoid int overflow on right hand side...
	physMemUsed *= memInfo.mem_unit;

//	cout << endl; 
//	cout << " ******************************************************************************************************************************* " << endl; 
//	cout << " *** Memory Info ... " << endl; 
//	cout << "  - Total RAM memory : " << totalPhysMem << endl; 
//	cout << "  - Used RAM  memory : " << physMemUsed << endl;
//	cout << "  - Usage            : " << 100.0 * (double)physMemUsed/(double)totalPhysMem<< " % " << endl ; 
//	cout << " ******************************************************************************************************************************* " << endl; 
//	cout << endl; 
	cout << " - Total used RAM : " << 100.0 * (double)physMemUsed/(double)totalPhysMem<< " %  - " << physMemUsed/1024/1024 << " MB of " << totalPhysMem/1024/1024 << " MB" << endl;
	 
#endif
}

int parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}
void GetProcessMemoryStat ( ) {

#ifndef WIN32

	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);

	int MPIrank; 

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	//cout << endl; 
	//cout << " ******************************************************************************************************************************* " << endl; 
	cout << " - Memory used by process " << MPIrank << " : " << result / 1024.0 << " MB"<< endl; 
	//cout << " ******************************************************************************************************************************* " << endl; 
	//cout << endl; 

#endif

}

double GetProcessMemory ( ) {

#ifndef WIN32

	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);

	return result / 1024.0;
#else 
	return 0.0;
#endif

}

// *** END - Memory functions ******************************************************
 

void SetMatrixFromCSR   ( SparseMatrix & Mat, ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type ) {

	int nnz = rows[n_rows];  
	int offset = (rows[0]) ? 0 : 1;
	nnz -= rows[0];

	Mat.CSR_I_row_indices.resize(n_rows+1);
	Mat.CSR_J_col_indices.resize(nnz);
	Mat.CSR_V_values	 .resize(nnz); 

	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin()); 
	for (int i = 0; i < Mat.CSR_I_row_indices.size(); i++)
		Mat.CSR_I_row_indices[i] = rows[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin()); 
	for (int i = 0; i < Mat.CSR_J_col_indices.size(); i++)
		Mat.CSR_J_col_indices[i] = cols[i] + offset;

	copy(vals, vals + nnz, Mat.CSR_V_values.begin()); 

	Mat.cols = n_cols;
	Mat.rows = n_rows;
	Mat.nnz  = nnz; 
	Mat.type = type; 

	//SpyText(Mat); 

	char f = 'u'; 
}

void SetMatrixFromCOO   ( SparseMatrix & Mat, ShortInt n_rows, ShortInt n_cols, ShortInt nnz, ShortInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing ) {

	Mat.I_row_indices.resize(nnz);
	Mat.J_col_indices.resize(nnz);
	Mat.V_values	 .resize(nnz); 
	int offset = indexing ? 0 : 1;

	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin()); 
	for (int i = 0; i < Mat.I_row_indices.size(); i++)
		Mat.I_row_indices[i] = I_rows[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin()); 
	for (int i = 0; i < Mat.J_col_indices.size(); i++)
		Mat.J_col_indices[i] = J_cols[i] + offset;

	copy(V_vals, V_vals + nnz, Mat.V_values.begin()); 

	Mat.cols = n_cols;
	Mat.rows = n_rows;
	Mat.nnz  = nnz; 
	Mat.type = type; 

	//Mat.ConvertToCSR( 0 );

	//SpyText(Mat); 

	char f = 'u';  
}

//POZOR - XC
//void SetMatrixFromCOO   ( SparseMatrix & Mat, longInt n_rows, ShortInt n_cols, ShortInt nnz, longInt * I_rows, ShortInt * J_cols, double * V_vals, char type ) {
//
//	Mat.I_row_indices.resize(nnz);
//	Mat.J_col_indices.resize(nnz);
//	Mat.V_values	 .resize(nnz); 
//
//	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin()); 
//	for (int i = 0; i < Mat.I_row_indices.size(); i++)
//		Mat.I_row_indices[i] = I_rows[i] + 1; 
//
//	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin()); 
//	for (int i = 0; i < Mat.J_col_indices.size(); i++)
//		Mat.J_col_indices[i] = J_cols[i] + 1; 
//
//	copy(V_vals, V_vals + nnz, Mat.V_values.begin()); 
//
//	Mat.cols = n_cols;
//	Mat.rows = n_rows;
//	Mat.nnz  = nnz; 
//	Mat.type = type; 
//
//	//Mat.ConvertToCSR( 0 );
//
//	//SpyText(Mat); 
//
//	char f = 'u';  
//}


void SetMatrixFromDense ( SparseMatrix & Mat, ShortInt n_cols, ShortInt n_rows, double * vals, char type ) {

	int dense_size = n_cols * n_rows;  

	Mat.dense_values.resize( dense_size ); 
	copy(vals, vals + dense_size, Mat.dense_values.begin()); 

	Mat.rows = n_rows;
	Mat.cols = n_cols; 
	Mat.type = type; 

	Mat.ConvertDenseToCSR(1); 

	//SpyText(Mat); 

	char f = 'u'; 

}

void SetVecInt ( SEQ_VECTOR <int> & vec, ShortInt incerement_by, ShortInt nnz, ShortInt * vals) {
	vec.resize(nnz);
	//copy( vals, vals + nnz, vec.begin() ); 
	for (int i = 0; i < nnz; i++)
		vec[i] = vals[i] + incerement_by; 

	char f = 'u'; 
}

void SetVecDbl ( SEQ_VECTOR <double> & vec, ShortInt nnz, double * vals) {
	vec.resize(nnz);
	copy( vals, vals + nnz, vec.begin() ); 
	char f = 'u'; 
}


void SetCluster ( Cluster & cluster, ShortInt * subdomains_global_indices_o, ShortInt number_of_subdomains, ShortInt MPI_rank) {
	
	SEQ_VECTOR             <int>   subdomains_global_indices  ( number_of_subdomains );
	
	for (int i = 0; i < number_of_subdomains; i++)
		subdomains_global_indices[i] = (int)subdomains_global_indices_o[i];
	
	cluster.cluster_global_index = MPI_rank + 1; 
	cluster.InitClusterPC(&subdomains_global_indices[0], number_of_subdomains);  

}

void SetMatrixB1_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster, 
						   longInt   n_rows, ShortInt n_cols, ShortInt nnz, 
						   longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing ) {

	int MPIrank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

#if DEBUG == 1
	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << " *** SetMatrixB1_fromCOO ******************************************************************************************************* " << endl; 
	}
	
	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Loading B1 ... " << endl; 
	}
#endif
							   
	// *** Setup B1 matrix ********************************************************************************************
	SetMatrixFromCOO( cluster.domains[domain_index_in_cluster].B1,
		n_rows, n_cols, nnz, 
		I_rows, J_cols, V_vals, type, indexing );
	
	// min version - should work 
	//TODO:Sort nejak blbne - nekdy sekne program - pocat na ten od Ondry
	//cluster.domains[domain_index_in_cluster].B1.sortInCOO();
	cluster.domains[domain_index_in_cluster].B1t = cluster.domains[domain_index_in_cluster].B1; 
	cluster.domains[domain_index_in_cluster].B1t.MatTransposeCOO();
	cluster.domains[domain_index_in_cluster].B1t.ConvertToCSRwithSort(1); 
	// end - min version - should work 
	

	////// Toto tak nejak funguje - ale je to divne 
	////cluster.domains[domain_index_in_cluster].B1t = cluster.domains[domain_index_in_cluster].B1;

	////SparseMatrix B1_tmp_bs;
	////B1_tmp_bs = cluster.domains[domain_index_in_cluster].B1;
	////
	////sortMatrixInCOO(cluster.domains[domain_index_in_cluster].B1); 

	////int test1 = B1_tmp_bs.MatCompareCOO(cluster.domains[domain_index_in_cluster].B1); 
	//////printf("B1 convert test 1 = %d \n",test1);
	////	
	////cluster.domains[domain_index_in_cluster].B1t.I_row_indices.swap(cluster.domains[domain_index_in_cluster].B1t.J_col_indices);
	////int tmp = cluster.domains[domain_index_in_cluster].B1t.rows;
	////cluster.domains[domain_index_in_cluster].B1t.rows = cluster.domains[domain_index_in_cluster].B1t.cols; 
	////cluster.domains[domain_index_in_cluster].B1t.cols = tmp;  
	////cluster.domains[domain_index_in_cluster].B1t.ConvertToCSRwithSort(1);

	////// POZOR !! - funkcni ale narocne na pamet pro velke ulohy
	////cluster.domains[domain_index_in_cluster].B1 = cluster.domains[domain_index_in_cluster].B1t; 
	////cluster.domains[domain_index_in_cluster].B1.ConvertToCOO(1);
	////cluster.domains[domain_index_in_cluster].B1.I_row_indices.swap(cluster.domains[domain_index_in_cluster].B1.J_col_indices);
	////tmp = cluster.domains[domain_index_in_cluster].B1.rows;
	////cluster.domains[domain_index_in_cluster].B1.rows = cluster.domains[domain_index_in_cluster].B1.cols; 
	////cluster.domains[domain_index_in_cluster].B1.cols = tmp;  

	////sortMatrixInCOO(cluster.domains[domain_index_in_cluster].B1); 

	////SparseMatrix B1_tmp; 
	////B1_tmp = cluster.domains[domain_index_in_cluster].B1; 

	////cluster.domains[domain_index_in_cluster].B1.ConvertToCSRwithSort(1); 
	////cluster.domains[domain_index_in_cluster].B1.ConvertToCOO(1);

	////int test2 = B1_tmp.MatCompareCOO(cluster.domains[domain_index_in_cluster].B1); 
	////if (test2 != 0) {
	////	//printf("Iter MPI 1 - B1 convert test 2 = %d \n",test2);
	////	//exit(2);  
	////}
	//// 
	////cluster.domains[domain_index_in_cluster].B1 = B1_tmp;  

	//////printf("B1 convert test 2 = %d \n",test2);
	////
	////// END - Toto tak nejak funguje - ale je to divne 

	// *** END - Setup B1 matrix ****************************************************************************************

#if DEBUG == 1
	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
	}
	
	if (MPIrank == 0) {
		cout << " *** END - SetMatrixB1_fromCOO ************************************************************************************************* " << endl; 
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << endl; 
	}
#endif

} 

void SetMatrixB0_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster, 
	longInt   n_rows, ShortInt   n_cols, ShortInt nnz, 
	longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing ) {

		int MPIrank; 
		MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

#if DEBUG == 1
		if (MPIrank == 0) {
			cout << " ******************************************************************************************************************************* " << endl; 
			cout << " *** SetMatrixB0_fromCOO ******************************************************************************************************* " << endl; 
		}

		if (MPIrank == 0) {
			GetProcessMemoryStat ( );
			GetMemoryStat( );
			cout << " Solver - Loading B0 ... " << endl; 
		}
#endif

		// *** Setup B0 matrix ********************************************************************************************
		SetMatrixFromCOO( cluster.domains[domain_index_in_cluster].B0,
			n_rows, n_cols, nnz, 
			I_rows, J_cols, V_vals, type, indexing );

		cluster.domains[domain_index_in_cluster].B0.ConvertToCSRwithSort(1);

		//cluster.domains[domain_index_in_cluster].B1t = cluster.domains[domain_index_in_cluster].B1; 
		//cluster.domains[domain_index_in_cluster].B1t.I_row_indices.swap(cluster.domains[domain_index_in_cluster].B1t.J_col_indices);
		//int tmp = cluster.domains[domain_index_in_cluster].B1t.rows;
		//cluster.domains[domain_index_in_cluster].B1t.rows = cluster.domains[domain_index_in_cluster].B1t.cols; 
		//cluster.domains[domain_index_in_cluster].B1t.cols = tmp; 
		////cluster.domains[domain_index_in_cluster].B1t.sortCOObyROW();
		//cluster.domains[domain_index_in_cluster].B1t.ConvertToCSR(0);
		// *** END - Setup B1 matrix ****************************************************************************************

#if DEBUG == 1
		if (MPIrank == 0) {
			GetProcessMemoryStat ( );
			GetMemoryStat( );
		}

		if (MPIrank == 0) {
			cout << " *** END - SetMatrixB1_fromCOO ************************************************************************************************* " << endl; 
			cout << " ******************************************************************************************************************************* " << endl; 
			cout << endl; 
		}
#endif

} 


void SetMatrixR_fromDense( Cluster & cluster, ShortInt domain_index_in_cluster,
						   ShortInt n_cols, ShortInt n_rows, double * vals, char type ) {

	int MPIrank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

#if DEBUG == 1	
	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << " *** SetMatrixR_fromDense ****************************************************************************************************** " << endl; 
	}

	if (MPIrank == 0 ) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Loading R ... " << endl; 
	}
#endif

	SetMatrixFromDense( cluster.domains[domain_index_in_cluster].Kplus_R, 
						n_cols, n_rows, vals, type ); 
#if DEBUG == 1
	if (MPIrank == 0 ) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
	}


	if (MPIrank == 0) {
		cout << " *** END - SetMatrixR_fromDense ************************************************************************************************ " << endl; 
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << endl; 
	}
#endif

}

void SetSolverPreprocessing ( Cluster & cluster, IterSolver & solver, 
							  SEQ_VECTOR <SEQ_VECTOR <longInt> > & lambda_map_sub_o, SEQ_VECTOR < ShortInt > & neigh_domains_o ) 
{

	SEQ_VECTOR <SEQ_VECTOR <int> > lambda_map_sub ( lambda_map_sub_o.size() );
	for (int d = 0; d < lambda_map_sub_o.size(); d++) {
		lambda_map_sub[d].resize( lambda_map_sub_o[d].size() );
		for (int i = 0; i < lambda_map_sub_o[d].size(); i++) {
			lambda_map_sub[d][i] = lambda_map_sub_o[d][i]; 
		}
	}

	SEQ_VECTOR             <int>   neigh_domains  ( neigh_domains_o.size()  );
	for (int i = 0; i < neigh_domains_o.size(); i++)
		neigh_domains[i] = (int)neigh_domains_o[i];





	int MPIrank; 
	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << " *** SetSolverPreprocessing **************************************************************************************************** " << endl; 
	}

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Creating G1 per cluster ... " << endl; 
	}


	TimeEvent G1_perCluster_time ("Setup G1 per Cluster time - preprocessing");
	G1_perCluster_time.AddStart(omp_get_wtime());

	TimeEvent G1_perCluster_mem ("Setup G1 per Cluster mem - preprocessing");
	G1_perCluster_mem.AddStartWOBarrier(GetProcessMemory());
	
	cluster.Create_G1_perCluster   ();

	G1_perCluster_time.AddEnd(omp_get_wtime());
	G1_perCluster_time.PrintStatMPI(0.0); 
	
	G1_perCluster_mem.AddEndWOBarrier(GetProcessMemory());
	G1_perCluster_mem.PrintStatMPI(0.0); 

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - CreateVec_d_perCluster ... " << endl; 
	}

	TimeEvent Vec_d_perCluster_time ("Setup Vec d per Cluster - preprocessing");
	Vec_d_perCluster_time.AddStart(omp_get_wtime());

	cluster.CreateVec_d_perCluster ();

	Vec_d_perCluster_time.AddEnd(omp_get_wtime());
	Vec_d_perCluster_time.PrintStatMPI(0.0); 

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Creating Global G1 and Running preprocessing (create GGt) ... " << endl; 
	}

	TimeEvent solver_Preprocessing_time ("Setup solver.Preprocessing() - preprocessing");
	solver_Preprocessing_time.AddStart(omp_get_wtime());

	cluster.my_neighs = neigh_domains; 

	solver.Preprocessing ( cluster );

	solver_Preprocessing_time.AddEnd(omp_get_wtime());
	solver_Preprocessing_time.PrintStatMPI(0.0); 

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - SetClusterPC - lambda map sub a neigh domains ... " << endl; 
	}

	TimeEvent cluster_SetClusterPC_time ("Setup cluster.SetClusterPC() - preprocessing");
	cluster_SetClusterPC_time.AddStart(omp_get_wtime());

	cluster.SetClusterPC( lambda_map_sub, neigh_domains ); // USE_DYNAMIC, USE_KINV 

	cluster_SetClusterPC_time.AddEnd(omp_get_wtime());
	cluster_SetClusterPC_time.PrintStatMPI(0.0); 

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
	}

	if (MPIrank == 0) {
		cout << " *** END - SetSolverPreprocessing ********************************************************************************************** " << endl; 
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << endl; 
	}


}

void SetMatrixK_fromCSR ( Cluster & cluster, ShortInt domain_index_in_cluster, 
						  ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type ) {

    
    int MPIrank = 1;
//	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

#if DEBUG == 1
	if (MPIrank == 0) {
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << " *** SetMatrixK_fromCSR ******************************************************************************************************** " << endl; 
	}

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Loading K make a copy of K ... " << endl;  
	}
#endif

	SetMatrixFromCSR( cluster.domains[domain_index_in_cluster].K, 
 					    n_rows, n_cols, 
					    rows, cols, vals, type); 
	
	if (cluster.USE_DYNAMIC == 1) {
		// suppose M in domain is already set
		double time_const = 1.0 / ( cluster.dynamic_beta * cluster.dynamic_timestep * cluster.dynamic_timestep);
		cluster.domains[domain_index_in_cluster].K.MatAddInPlace(cluster.domains[domain_index_in_cluster].M,'N', time_const);
	}

	if ( cluster.domains[domain_index_in_cluster].K.type == 'G' )
		cluster.domains[domain_index_in_cluster].K.RemoveLower();

#if DEBUG == 1
	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Loading K, regularization and factorization ... " << endl;  
		cluster.domains[domain_index_in_cluster].Kplus.msglvl = 0; 
	}
#endif

	// POZOR - zbytecne kopiruju - pokud se nepouziva LUMPED
	cluster.domains[domain_index_in_cluster].Prec = cluster.domains[domain_index_in_cluster].K;
	//cluster.domains[domain_index_in_cluster].Prec.ConvertCSRToDense( 0 ); 

	cluster.domains[domain_index_in_cluster].K_regularizationFromR( );
	cluster.domains[domain_index_in_cluster].domain_prim_size = cluster.domains[domain_index_in_cluster].Kplus.cols; 

	cluster.domains[domain_index_in_cluster].K.Clear();

	if ( cluster.cluster_global_index == 1 ) { GetMemoryStat( ); GetProcessMemoryStat ( ); }

#if DEBUG == 1
	if (MPIrank == 0) {
		cluster.domains[domain_index_in_cluster].Kplus.msglvl = 0; 
		GetProcessMemoryStat ( );
		GetMemoryStat( );
		cout << " Solver - Running SetClusterPC_AfterKplus ...  " << endl; 
	}

	if (MPIrank == 0) {
		GetProcessMemoryStat ( );
		GetMemoryStat( );
	}
	
	if (MPIrank == 0) {
		cout << " *** END - SetMatrixK_fromCSR ************************************************************************************************** " << endl; 
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << endl; 
	}
#endif
}



void SetMatrixK_fromBEM ( Cluster & cluster, ShortInt domain_index_in_cluster,
                         ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type ) {
    

    int MPIrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
    
#if DEBUG == 1
    if (MPIrank == 0) {
        cout << " ******************************************************************************************************************************* " << endl;
        cout << " *** SetMatrixK_fromCSR ******************************************************************************************************** " << endl;
    }
    
    if (MPIrank == 0) {
        GetProcessMemoryStat ( );
        GetMemoryStat( );
        cout << " Solver - Loading K make a copy of K ... " << endl;
    }
#endif
    
    SetMatrixFromCSR( cluster.domains[domain_index_in_cluster].K,
                     n_rows, n_cols,
                     rows, cols, vals, type);
    
    if ( cluster.domains[domain_index_in_cluster].K.type == 'G' )
        cluster.domains[domain_index_in_cluster].K.RemoveLower();
    
#if DEBUG == 1
    if (MPIrank == 0) {
        GetProcessMemoryStat ( );
        GetMemoryStat( );
        cout << " Solver - Loading K, regularization and factorization ... " << endl;
        cluster.domains[domain_index_in_cluster].Kplus.msglvl = 0;
    }
#endif
    
    // POZOR - zbytecne kopiruju - pokud se nepouziva LUMPED
    cluster.domains[domain_index_in_cluster].Prec = cluster.domains[domain_index_in_cluster].K;
    //cluster.domains[domain_index_in_cluster].Prec.ConvertCSRToDense( 0 );
    
    cluster.domains[domain_index_in_cluster].K_regularizationFromR( );

    //cluster.domains[domain_index_in_cluster].Kplus.ImportMatrix  (cluster.domains[domain_index_in_cluster].K);
    //cluster.domains[domain_index_in_cluster].Kplus.Factorization ();
    
    // POZOR - jen pro Kinv
    if (cluster.domains[domain_index_in_cluster].USE_KINV == 1 || cluster.domains[domain_index_in_cluster].USE_HFETI == 1)
        cluster.domains[domain_index_in_cluster].KplusF.ImportMatrix (cluster.domains[domain_index_in_cluster].K);

    cluster.domains[domain_index_in_cluster].K.Clear();
    
    
    
    
    cluster.domains[domain_index_in_cluster].domain_prim_size = cluster.domains[domain_index_in_cluster].Kplus.cols;
    
#if DEBUG == 1
    if (MPIrank == 0) {
        cluster.domains[domain_index_in_cluster].Kplus.msglvl = 0;
        GetProcessMemoryStat ( );
        GetMemoryStat( );
        cout << " Solver - Running SetClusterPC_AfterKplus ...  " << endl;
    }
    
    if (MPIrank == 0) {
        GetProcessMemoryStat ( );
        GetMemoryStat( );
    }
    
    if (MPIrank == 0) {
        cout << " *** END - SetMatrixK_fromCSR ************************************************************************************************** " << endl;
        cout << " ******************************************************************************************************************************* " << endl; 
        cout << endl; 
    }
#endif
 
 
}




#ifdef LIB
	int Xmain(int argc, char* argv[])
#else 
	int main(int argc, char* argv[])
#endif
{

    int mpi_root = 0;
	int mpi_rank, mpi_size; 

	MPI_Init (&argc, &argv);					// starts MPI 
	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	// get current process id 
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	// get number of processes 

#ifdef _OPENMP
	if (mpi_rank == 0)
		printf("Compiled with OpenMP Support \n");
#else
	if (mpi_rank == 0)
		printf("Compiled without OpenMP Support \n");
#endif
	 
	if (argc < 10) {
		if (mpi_rank == mpi_root)
			cout << "11 input parameters are required - all parameters must be entered in this exact order: " << endl << 
				"  - [0]  directory path                                       " << endl << 
				"  - [1]  CG type:                           PIPECG   or REGCG " << endl << 
				"         - compile with -DUSE_MPI_3 to enable non blocking communication in PIPECG - not set by default" << endl <<
				"  - [2]  Prec type:                         LUMPED   or NONE  " << endl <<
				"  - [3]  Stop condition:                    i.e. 0.0001       " << endl <<
				"  - [4]  Use GGt inverse:                   GGTINV   or NONE  " << endl <<
				"  - [5]  Find solution:                     FINDSOL  or NONE  " << endl <<
				"  - [6]  Turn on MKL CNR (reproducibility): USECNR   or NOCNR " << endl <<
				"  - [7]  Dynamic						   : DYNAMIC  or NONE  " << endl <<
				"  - [8]  Use Hybrid FETI instead of FETI  : HFETI    or FETI  " << endl << 
				"  - [9]  - for FETI only - number of subdomains per MPI process " << endl << 
				"  - [10] - for FETI only - total number of subdomains " << endl; 
				; 
		return -1; 
	}


	int USE_PIPECG = 0;
	if (string(argv[2]) == "PIPECG" || string(argv[2]) == "pipecg") {
		USE_PIPECG = 1;
	} else {
		USE_PIPECG = 0; 
	}


	int USE_PREC   = 0; 
	if (string(argv[3]) == "LUMPED" || string(argv[3]) == "lumped") {
		USE_PREC = 1;
	} else {
		USE_PREC = 0; 
	}


	// Sets the stop criteria from command line 
	double epsilon = 1.0e-4; 
	epsilon = strtod(argv[4], NULL);


	// Disable/Enable use of inverse of GGt matrix 
	int USE_GGtINV = 0; 
	if (string(argv[5]) == "GGTINV" || string(argv[5]) == "ggtinv") {
		USE_GGtINV = 1;
	} else {
		USE_GGtINV = 0; 
	}


	// Disable/Enable find solution at the end 
	int FIND_SOLUTION = 0; 
	if (string(argv[6]) == "FINDSOL" || string(argv[6]) == "findsol") {
		FIND_SOLUTION = 1;
	} else {
		FIND_SOLUTION = 0; 
	}


	int USE_CNR = 0; 
	if (string(argv[7]) == "USECNR" || string(argv[7]) == "usecnr") {
		USE_CNR = 1;
		mkl_cbwr_set(MKL_CBWR_COMPATIBLE); 
		//mkl_cbwr_set(MKL_CBWR_AUTO); 
	} else {
		USE_CNR = 0; 
	}



	int USE_DYNAMIC     = 0; 
	double const_beta   = 0; 
	double const_gama   = 0;
	double const_deltat = 0; 
	
	if (string(argv[8]) == "DYNAMIC" || string(argv[7]) == "dynamic") {
		USE_DYNAMIC = 1;

		const_beta   = 0.25; 
		const_gama   = 0.5; 
		const_deltat = strtod(argv[9], NULL); //0.000000001; 
	} else {
		USE_DYNAMIC = 0; 
	}



	int USE_KINV = 0;
	if (string(argv[10]) == "KINV" || string(argv[10]) == "kinv") {
		USE_KINV = 1;
	} else {
		USE_KINV = 0; 
	}
	


	int USE_HFETI = 0; 
	if (string(argv[11]) == "HFETI" || string(argv[11]) == "hfeti") {
		USE_HFETI = 1;
	} else {
		USE_HFETI = 0; 
	}

	int SUBDOMS_PER_CLUSTER = strtol(argv[12], NULL, 0);
	int SUBDOMS_TOTAL       = strtol(argv[13], NULL, 0);



	// *** Reading the config file **************************************************************
	ConfigFile cf("config.txt");

	// *** Read all sections of the config file *************************************************
	std::SEQ_VECTOR<std::string> sections = cf.GetSections();

	// *** List all the sections in the config file *********************************************
	if (mpi_rank == mpi_root) {
		for (std::SEQ_VECTOR<std::string>::iterator i = sections.begin(); i < sections.end(); i++) {
			std::cout<<"Section: " << *i << std::endl;
		}
	}

	// *** Read particula items from the config file *******************************************
	std::string foo;
	foo = cf.Value("section_1","foo"     , "?");
	
	if (mpi_rank == mpi_root) {
		//std::cout << foo     << std::endl;
	}
	
	
	
	
	// Print configuration as translated from command line arguments 
	if (mpi_rank == mpi_root) {
		cout << "Using directory .................................. : " << argv[1]		 << endl;
		cout << "Using Pipeline CG (1 - PIPECG, 0 - Regular CG)     : " << USE_PIPECG	 << endl;
		cout << "Using Preconditioner (1 - LUMPED, 0 - NONE) ...... : " << USE_PREC		 << endl;
		cout << "Using distributed inverse matrix of GGt  (1 - yes) : " << USE_GGtINV	 << endl; 
		cout << "Find solution  (1 - yes, 0 - no) ................. : " << FIND_SOLUTION << endl;
		cout << "MKL CNR (reproducibility) (1 - on, 0 - off)        : " << USE_CNR       << endl;
		cout << "Stop condition ................................... : " << epsilon		 << endl; 

#ifdef USE_MPI_3
		cout << "USE_MPI_3 is enabled for PIPECG                    : " << endl;
#else
		cout << "USE_MPI_3 is disabled for PIPECG                   : " << endl;
#endif
		cout << "USE_DYNAMIC ...................................... : " << USE_DYNAMIC	 << endl;
		cout << "USE_KINV                                           : " << USE_KINV		 << endl;
		cout << "USE_HFETI ........................................ : " << USE_HFETI	 << endl;
		cout << "FOR FETI ONLY - SUBDOMAINS PER MPI PROCESS         : " << SUBDOMS_PER_CLUSTER << endl;
		cout << "FOR FETI ONLY - TOTAL NUM OF SUBDOMAINS .......... : " << SUBDOMS_TOTAL << endl;
		
		cout << endl << endl; 
	}
	// END - Print configuration as translated from command line arguments 

	TimeEval preproc_timing ("Preprocessing timing ");

	// Define where the input data and log file and result file are stored
	string result_file;
	string log_file; 
#ifdef WIN32
	string directory = argv[1] + string("/");
	result_file = string(directory) + string("/PrePost_data/MATSOL_SVN_Displacement.Nvec"); 
	log_file    = string(directory) + string("log_file.txt"); 
#else
	string directory = argv[1] + string("/");
	result_file = string(directory) + string("/PrePost_data/MATSOL_SVN_Displacement.Nvec"); 
	log_file    = string(directory) + string("log_file.txt"); 
#endif

	//// Open LOG file 
	//FILE *stream;
	//int log_active = 0; 

	//if (rank == mpi_root) {
	//		
	//	// Open for write 
	//	if( (stream = fopen( log_file.c_str(), "a" )) == NULL ) {
	//	
	//		cout << "Log file cannot be created " << log_file << endl;
	//		log_active = 0; 

	//	} else {
	//		
	//		cout << "Log file was created : " << log_file << endl;
	//		log_active = 1; 

	//		char   s0[] = "\n\n******************************************************************** \n";
	//		fprintf( stream, "%s", s0 );

	//		time_t now = time(0);
	//		char* dt = ctime(&now);
	//		cout << "The local date and time is: " << dt << endl;
	//		fprintf( stream, "The local date and time is: %s \n\n", dt );

	//	}
	//}

	preproc_timing.totalTime.AddStart(omp_get_wtime());
	
	// ***************************************************************************
	// *** Load the cluster - one for each MPI rank ******************************
	TimeEvent loadcluster_time("Time to load cluster"); 
	loadcluster_time.AddStart(omp_get_wtime());
	
	Cluster cluster(mpi_rank + 1); // Clusters are numbered from 1 not from 0
	cluster.USE_HFETI		   = USE_HFETI;
	cluster.SUBDOM_PER_CLUSTER = SUBDOMS_PER_CLUSTER; 

	cluster.SetDynamicParameters(const_deltat, const_beta, const_gama); 
	cluster.LoadCluster(directory, USE_DYNAMIC, USE_KINV);	
	 

	loadcluster_time.AddEnd(omp_get_wtime());
	loadcluster_time.PrintStatMPI(0.0); 

	preproc_timing.AddEvent(loadcluster_time); 
	// *** END - Load the cluster - one for each MPI rank ************************
	// ***************************************************************************

	IterSolver solver;  
	solver.epsilon		 = epsilon; 
	solver.USE_HFETI     = USE_HFETI;
	solver.USE_KINV		 = USE_KINV;
	solver.USE_DYNAMIC	 = USE_DYNAMIC;
	solver.USE_GGtINV	 = USE_GGtINV; 
	solver.USE_PIPECG	 = USE_PIPECG;
	solver.USE_PREC		 = USE_PREC; 
	solver.FIND_SOLUTION = FIND_SOLUTION;					

	
	SEQ_VECTOR<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );
	solver.Preprocessing ( cluster );
	

	if (USE_DYNAMIC == 1)
		; //solver.Solve_Dynamic  ( cluster, result_file );
	else 
		solver.Solve_singular ( cluster, result_file );
   	
	
 	MPI_Finalize();
	return 0;
}






