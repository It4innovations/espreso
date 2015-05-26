
#include "utility.h"
#include "Fem.h"
#include "Cluster_g.h"
#include "SaveFEM.h"
#include "Solver.h"
#include "Data.h"


#ifdef HFETI_SOLVER
#include "Cluster.h"
#include "IterSolver.h"
#include "SparseMatrix.h"
#include "TimeEval.h"
#endif



/* -----------------------------------------------------------------------*/
/* FOTRAN lib.
 * 0  - - - library - off - - -
 * 1  MPI_Fint                            main
 * 2  init_general_                       main
 * 3  fe2feti_init_symbolic_              fe_symbolic_Form_matrix
 * 4  fe2feti_symbolic_map_element_       fe_symbolic_form_matrix (in loop)
 * 5  fe2feti_symbolic_map_global_bc_     fe_symbolic_form_matrix
 * 6  fe2feti_symbolic_finalize_          fe_symbolic_form_matrix
 * 7  fe2feti_symbolic_factorize_         fe_symbolic_form_matrix
 * 8  fe2feti_init_numeric_               fe_symbolic_form_matrix ???
 * 9  fe2feti_numeric_map_element_           stima_subdomain         (in loop)
 * 10 fe2feti_numeric_factorize_          stima_subdomain
 * 11 fe2feti_numeric_map_global_bc_      definitionAndSolving
 * 12 fe2feti_map_global_rhs_             definitionAndSolving
 * 13 fe2feti_solve_                      definitionAndSolving
 * 14 fe2feti_finalize_                   main
 * 15 fe2feti_free_data_                  main
 */
/* -----------------------------------------------------------------------*/
// 

int main(int argc, char *argv[]) {
//
#ifdef HFETI_SOLVER
	extern void GetMemoryStat	   ( ); 
	extern void GetProcessMemoryStat ( ); 
	extern double GetProcessMemory ( ); 
#endif
//
#ifdef FLLOP_ENABLED
  FllopAIFInitialize(&argc, &argv, "flloprc");
  PetscLogStageRegister("PermonCube Ass.", &assemble_stage);
  PetscLogStageRegister("PermonCube Sol.", &solve_stage);
  PetscLogStagePush(assemble_stage);
  MPI_Comm comm=PETSC_COMM_WORLD;
#else
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
#endif
//
	int MPIrank;
  MPI_Comm_rank(comm, &MPIrank);

#ifndef WIN32
  if (MPIrank == 0) { cout << " Generator started - System MPI initialization finished                   "; system("date +%T.%6N"); }
#endif

//
	CClust_g clust(comm,argc, argv);
  clust.getGlobalIndexOfSubdomMeshGen();
	//CFem fem(comm);
	//CData data(comm);
	//
  //GetMemoryStat(MPIrank);
 // 
#ifdef HFETI_SOLVER
	if (MPIrank == 0) {
		cout << endl; 
		cout << " ******************************************************************************************************************************* " << endl; 
		cout << " start permoncube ... " << endl; 
		GetProcessMemoryStat ( ); GetMemoryStat( );
		cout << endl; 
	}
#endif
	//
#ifdef flag_Fortran
	MPI_Fint c_com = MPI_Comm_c2f(MPI_COMM_WORLD);
	init_general_(c_com);
#endif
	//
  //set DIR indexes - - s
  int i_faceDir[1] ={0};// 
  int dirDirBC_xyz[3] ={1,1,1};// 1 - fixed displ., 0 - free displ.
  //set DIR indexes - - s
  int i_faceCon[1] ={1};// 
  int dirConBC_xyz[3] ={0,0,1};// 1 - fixed displ., 0 - free displ.

// -------------------------------------------------------------------

#ifndef WIN32
  if (MPIrank == 0) { cout << " mesh initialize start                                                    "; system("date +%T.%6N"); }
#endif

  cilk_for (int i = 0;i<clust.domainG->n_subdomOnClust;i++){
	  clust.fem[i]->mesh.initialize(clust.domainG,
                                  clust.domainG->vec_localSubNumbering[i],
                                  clust.domainG->vec_globalSubNumbering[i],
                                  clust.domainG->vec_globalSubNumberingMeshGen[i]);


    //printf("MPI=%d, localSub,globalSub,globalSubMeshGen %d,%d,%d\n",
    //MPIrank,clust.domainG->vec_localSubNumbering[i],
    //clust.domainG->vec_globalSubNumbering[i],
    //clust.domainG->vec_globalSubNumberingMeshGen[i]
    //);

  }

//#ifndef WIN32
//  if (MPIrank == 0) { cout << " Loop start                                                               "; system("date +%T.%6N"); }
//#endif
//  for (int i = 0;i<clust.domainG->n_subdomOnClust;i++){

#ifndef WIN32
	  if (MPIrank == 0) { cout << " Mesh Generator start                                                     "; system("date +%T.%6N"); }
#endif
	cilk_for (int i = 0;i<clust.domainG->n_subdomOnClust;i++)
		clust.fem[i]->mesh_generator3d(clust.domainG);

	  // creating element, coordinate matrices ...

#ifndef WIN32
	  if (MPIrank == 0) { cout << " DataDirBCSub start                                                       "; system("date +%T.%6N"); }
#endif

    cilk_for (int i = 0;i<clust.domainG->n_subdomOnClust;i++)
		clust.fem[i]->dataDirBCSub(&i_faceDir[0],1, clust.domainG->n_facesSub[clust.fem[i]->i_domOnClust],dirDirBC_xyz);

#ifndef WIN32
	if (MPIrank == 0) { cout << " dataConBCSub start                                                       "; system("date +%T.%6N"); }
#endif


    cilk_for (int i = 0;i<clust.domainG->n_subdomOnClust;i++)
		clust.fem[i]->dataConBCSub(clust.domainG,&i_faceCon[0],1, clust.domainG->n_facesSub[clust.fem[i]->i_domOnClust],dirConBC_xyz);

#ifndef WIN32
	if (MPIrank == 0) { cout << " initialize start                                                         "; system("date +%T.%6N"); }
#endif

	for (int i = 0;i<clust.domainG->n_subdomOnClust;i++) {
		clust.data[i]->initialize(comm, 
						clust.domainG->n_nodsSub[clust.fem[i]->i_domOnClust]*3,
						clust.domainG->neqAll,
						clust.domainG->max_nDOFs_u_is_stored,
						clust.domainG->n_elementsSub[clust.fem[i]->i_domOnClust],
						clust.domainG->n_exterDOFs[clust.fem[i]->i_domOnClust]);
	}

#ifndef WIN32
	if (MPIrank == 0) { cout << " extDofs start                                                            "; system("date +%T.%6N"); }
#endif

	cilk_for (int i = 0;i<clust.domainG->n_subdomOnClust;i++)
		clust.data[i]->extDofs(clust.fem[i], clust.domainG);

#ifndef WIN32
	if (MPIrank == 0) { cout << " extdofs - end                                                            "; system("date +%T.%6N"); }
#endif



//  }
//#ifndef WIN32
//  if (MPIrank == 0) { cout << " Loop - end                                                               "; system("date +%T.%6N"); }
//#endif

// -------------------------------------------------------------------
	CSolver::solving(clust.fem[0],clust.domainG, clust.data[0], clust);
#ifdef HFETI_SOLVER
	if (MPIrank == 0) {
		cout << endl; 
		cout << " after solver... " << endl; 
		GetProcessMemoryStat ( ); GetMemoryStat( );
		cout << endl; 
	}
#endif
  
//  CSaveFEM saveFEM(comm,clust.fem[0],clust.data[0],clust.domainG);

	if (clust.domainG->flag_store_VTK==1){
		clust.GatherDataFemToMaster();
  }
	if (clust.domainG->flag_store_VTK==2){
    clust.createVTK_per_cluster();
  }
	if (clust.domainG->flag_store_VTK==3){
		clust.GatherDataFemToMaster();
        clust.createVTK_per_cluster();
  }

#ifdef FLLOP_ENABLED
	TRY( PetscLogStagePop() );
	FllopAIFFinalize();
#else
#ifdef flag_Fortran
	fe2feti_finalize_(); 
	fe2feti_free_data_();
#endif
	MPI_Finalize();
#endif
	return 0;
	//
}
