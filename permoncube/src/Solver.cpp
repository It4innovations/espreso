//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#ifdef HFETI_SOLVER
#include "Cluster.h"
#include "IterSolver.h"
#include "SparseMatrix.h"
#include "TimeEval.h"
#endif

#include "Solver.h"
#include "LinearAlgebra.h"
#include "Assembler.h"
#include "SaveFEM.h"

CSolver::CSolver() {
}
CSolver::~CSolver() {
}

#ifdef HFETI_SOLVER
extern void GetMemoryStat	   ( ); 
extern void GetProcessMemoryStat ( ); 
extern double GetProcessMemory ( ); 
#endif

void CSolver::definitionAndSolving(CData *data, CFem *fem, CDomain *domainG,CClust_g & clust) {

  double fz_total = domainG->fz_total; 
//  CDomain *domainG = domainG;

  int MPIrank;
  int MPIsize;
  MPI_Comm_rank(fem->comm, &MPIrank);
  MPI_Comm_size(fem->comm, &MPIsize);

#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'definitionAndSolving', before 'Assembler.'.. " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif

#ifndef WIN32
  if (MPIrank == 0) { cout << " stima subdomain start                                                    "; system("date +%T.%6N"); }
#endif

//  CSolid45EqNumbGlobal *stif_glob_number = data->stif_glob_number;
  cilk_for (int i = 0; i < clust.domainG->n_subdomOnClust;i++) {
  CAssembler::stima_subdomain(clust.data[i]->f_subdom,   
                              clust.data[i]->f_subdom_internal,
                              clust.data[i]->u, 
                              clust.data[i]->du,
                              clust.data[i]->stif_glob_number,
                              clust.fem [i], domainG, 
                              clust.data[i]->KSparse);
//  CAssembler::stima_subdomain(data->f_subdom, data->f_subdom_internal,
//      data->u, data->du,
//      stif_glob_number,
//      fem, domainG, data->KSparse);
      }


#ifndef WIN32
  if (MPIrank == 0) { cout << " stima subdomain - end                                                    "; system("date +%T.%6N"); }
#endif


#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'definitionAndSolving', after 'Assembler.'.. " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif

  bool first_step = (domainG->i_load_step == 0 && domainG->i_sub_step == -1);
  // Equality conditions are created once ()
  if (first_step) {
#ifndef FLLOP_NEW_API
    //data->extDofs_l2g(domainG->n_exterDOFs[fem->i_domOnClust],fem->mesh.l2g);    

#ifndef WIN32
	  if (MPIrank == 0) { cout << " ExtDofs_l2g - start                                                      "; system("date +%T.%6N"); }
#endif

	  for (int i=0;i<clust.domainG->n_subdomOnClust;i++){
      clust.data[i]->extDofs_l2g(
                  domainG->n_exterDOFs[fem->i_domOnClust],
                  clust.fem[i]->mesh.l2g);    
	  }

#ifndef WIN32
	  if (MPIrank == 0) { cout << " get_DP_DOFs - start                                                      "; system("date +%T.%6N"); }
#endif

	  for (int i=0;i<clust.domainG->n_subdomOnClust;i++){
		clust.fem[i]->get_DP_DOFs(clust.domainG);
	  }

#ifndef WIN32
	  if (MPIrank == 0) { cout << " collecting external dofs - start                                         "; system("date +%T.%6N"); }
#endif

    
//#ifndef HFETI_SOLVER
    clust.collectingExternalDOFs();
//    data->B->createB(fem,domainG,data->KSparse->n_row,data->indExterDOFs); //!! POZOR

#ifndef WIN32
	if (MPIrank == 0) { cout << " collecting external DOFs - end                                           "; system("date +%T.%6N"); }
#endif

	MPI_Barrier(MPI_COMM_WORLD); 

#ifndef WIN32
	if (MPIrank == 0) { cout << " After MPI  Barrier                                                       "; system("date +%T.%6N"); }
#endif

	//MPI_Finalize();
	//exit(0); 


    clust.B1->createB1(clust.domainG,
                       clust.fem,clust.data,
                       clust.indExterDOFsOnClust,
                       clust.domainG->n_exterDOFsClust,
                       clust.neighbClst,
                       clust.domainG->n_neighbClst,
                       clust.domainG->flag_redund_lagr_mult);



//#endif
#endif
  }
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " 'definitionAndSolving' after 'createB'" << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif

  for (int i = 0;i<clust.domainG->n_subdomOnClust;i++){
  //TODO having function in 'clust.fem[i]->dataNeuBC' is not the best 
    clust.fem[i]->dataNeuBC(clust.data[i]->fE,domainG->neqSub[fem->i_domOnClust], 
              fz_total,5,domainG->n_facesSub[clust.fem[i]->i_domOnClust]);
    CLinearAlgebra::add_vec_a2b(clust.data[i]->f_subdom,clust.data[i]->fE, 1.0, 1.0,
        domainG->neqSub[clust.fem[i]->i_domOnClust]);
    if(domainG->i_sub_step==-1) {
      CLinearAlgebra::add_vec_a2b(clust.data[i]->f_subdom_internal,clust.data[i]->fE, 0.0,
          domainG->del_deltaVec,domainG->neqSub[clust.fem[i]->i_domOnClust]);
    } else {
      CLinearAlgebra::add_vec_a2b(clust.data[i]->f_subdom_internal,clust.data[i]->fE, -1.0,
          domainG->deltaVec[domainG->i_load_step],
          domainG->neqSub[clust.fem[i]->i_domOnClust]);
    }
  }
  //
#ifdef flag_Fortran
  if (!MPIrank && first_step) printf(" * * * fortran solver DDsolv used * * *\n");
  // 
  data->fGlobal = new double[domainG->neqAll];
  memset(&data->fGlobal[0],0,domainG->neqAll*sizeof(double));
  //
  double *tmp_fGlobal   = new double[domainG->neqAll];
  memset(&tmp_fGlobal[0],0,domainG->neqAll*sizeof(double));
  for (int i = 0;i<domainG->neqSub;i++){
    tmp_fGlobal[fem->mesh.l2g[i]] = data->fE[i]; 
  }
  MPI_Allreduce(tmp_fGlobal,data->fGlobal, domainG->neqAll, MPI_DOUBLE, MPI_SUM, 
      MPI_COMM_WORLD);
  delete [] tmp_fGlobal; tmp_fGlobal = NULL;
  //
  fe2feti_numeric_map_global_bc_(fem->bound_cond->dirBC_global->n,
      fem->bound_cond->dirBC_global->ind,  
      fem->bound_cond->dirBC_global->val);
  fe2feti_map_global_rhs_(domainG->neqAll,  data->fGlobal);
  fe2feti_solve_();
  fe2feti_gather_solution_(data->u_restrict,NULL);
  MPI_Bcast(data->u_restrict,domainG->neqAll,MPI_DOUBLE,0,MPI_COMM_WORLD);
  //
  for (int i=0;i<domainG->neqSub;i++){
    data->ddu[i] = data->u_restrict[fem->mesh.l2g[i]];
  }
  //
#else //#ifdef flag_Fortran
  //  solving problem
  // in fE Neumann BC, in f_subdom volume force
  //		CKSparse::saveStiffnessMatrix(fem);
//  if (domainG->flag_storeCompleteMatrices == 1) {
//    CSaveFEM saveFEM(fem->comm, fem, data);
//    saveFEM.save_data();
//  }

#ifdef PERMONCUBE_DUMP    
  // ------------------------------------------------------------------------
  int n_subdom =  clust.domainG->number_of_subdomains;//domainG->Nx*domainG->Ny*domainG->Nz;
  int maxDirichletIndex = ((domainG->nx + 1) * (domainG->ny + 1)) * 3;
  ofstream f_infoFile("data/infoFile.dat");
  // 
  f_infoFile << maxDirichletIndex << endl << domainG->neqAll << endl
    << domainG->n_elementsAll << endl << n_subdom <<endl 
    << domainG->Nx<<endl<<domainG->Ny<<endl<<domainG->Nz<<endl
    << domainG->nx<<endl<<domainG->ny<<endl<<domainG->nz<<endl;
  f_infoFile.close();
  // ------------------------------------------------------------------------
#endif

#ifdef FLLOP_ENABLED
#ifndef FLLOP_SKIP
  if (!MPIrank && first_step) printf(" * * * PETSc solver FLLOP used * * *\n");
  solver_fllop(data,fem,domainG);
#else
#warning "FLLOP solver [solver_fllop(data,fem)] will be skipped"
#endif
#else //FLLOP_ENABLED

#ifdef HFETI_SOLVER
  // *** Interface to HFETI solver ************************************************************
  if ( MPIrank == 0 )
	  cout << "Running HFETI Solver" << endl; 
  //solver_hfeti(data,fem,domainG);
  solver_hfeti(clust);
  // ******************************************************************************************
#else
  if (MPIsize==1){
	  double epsCG=domainG->eps0;
	  int maxIter=domainG->maxIter;
	  //    fem->applicationDirBCtoAxb(data->KSparse, data->fE);
	  if (!MPIrank && first_step) printf(" * * * built-in solver used * * *\n");
	  solver_simple_cg(domainG->neqSub[fem->i_domOnClust],
		  data->ddu,data->fE,epsCG,maxIter,data->KSparse,fem,domainG->verbose);
  } else {
	  cerr << "built-in CG is implemented only for sequential use" << endl;
	  throw -2;
  }
#endif // HFETI_SOLVER

#endif //FLLOP_ENABLED
 
#endif //#ifdef flag_Fortran

  if (domainG->neqAll>domainG->max_nDOFs_u_is_stored) {
    if (!MPIrank && first_step){
      printf(" ... u_restricted is not created  ... \n");
      printf(" ... currently nDOFs set to %d, don't forget pre-set to (at least ) 1e7 ...\n",
          domainG->max_nDOFs_u_is_stored);
    }
  }

  ++domainG->flag_i_callingSolver;  
}	// definitionAndSolving

void CSolver::solver_simple_cg(int neqSub, double *x, double *b, const double &tol,
    const int MaxIter, CKSparse *KSparse, CFem * fem,int verbose) {
  int MPIrank;
  MPI_Comm_rank(fem->comm, &MPIrank);
  clock_t begin = clock();
  int ITER = 0;
  int k = 0;
  double alpha = 0.;
  double beta = 0.;
  double tmp_p;
  double *r = new double[neqSub];
  double *rb = new double[neqSub];
  double *p = new double[neqSub];
  double *Ap = new double[neqSub];
  double dot_r_r = 0.0;
  double dot_r_r0 = 0.0;
  double dot_p_Ap = 0.0;
  double dot_rb_rb = 0.0;

  // zeroing x, r, rb
  memset(x, 0,neqSub * sizeof(double));
  memset(r, 0,neqSub * sizeof(double));
  memset(rb,0,neqSub * sizeof(double));

  // r0
  // via function multAx, r = A*x, r - output vector, x - input vector
  KSparse->multAx(r, x, neqSub, 1);

  dot_r_r = 0.0;
  for (int i = 0; i < neqSub; i++)
  {
    tmp_p = r[i];
    r[i] = (b[i] - tmp_p);
    dot_r_r += r[i] * r[i];
  }
  dot_r_r0 = dot_r_r;
  // copy r0 to p0
  memcpy(p,r,neqSub*sizeof(double));
  // iteration
  while ((ITER < MaxIter) && (sqrt(dot_r_r) > tol*sqrt(dot_r_r0))) {
    ITER++;
    if (verbose>3) {
      printf("iter; %d ||res||:  \t%1.2e\n", ITER,
          sqrt(dot_r_r / dot_r_r0));
    }
    // Ap
    KSparse->multAx(Ap, p, neqSub, 1);
    // dot products
    dot_r_r = 0.;
    dot_p_Ap = 0.;
    for (int i = 0; i < neqSub; i++) {
      dot_r_r += r[i] * r[i];
      dot_p_Ap += Ap[i] * p[i];
    }
    // alpha coefficient
    alpha = dot_r_r / dot_p_Ap;
    // x updating, rb (previous residual) storing
    for (int i = 0; i < neqSub; i++) {
      x[i] += alpha * p[i];
    }
    memcpy(rb,r,neqSub*sizeof(double));
    // iteration
    // dot product storing
    dot_rb_rb = dot_r_r;
    // residual updating
    for (int i = 0; i < neqSub; i++) {
      r[i] -= alpha * Ap[i];
    }
    // residual dot product
    dot_r_r = 0;
    for (int i = 0; i < neqSub; i++) {
      dot_r_r += r[i] * r[i];
    }
    // beta coefficient
    beta = dot_r_r / dot_rb_rb;
    // conjugate vector updating
    for (int i = 0; i < neqSub; i++) {
      tmp_p = p[i];
      p[i] = r[i] + beta * tmp_p;
    }
  } // end of loop

  // last residual
  if (verbose>2 ) {
    printf("iter; %d ||res||:  \t%1.2e\n", ITER + 1,
        sqrt(dot_r_r / dot_r_r0));
  }
  // fields zeroing
  delete[] r;     r   = NULL;
  delete[] rb;    rb  = NULL;
  delete[] p;     p   = NULL;
  delete[] Ap;    Ap  = NULL;

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  if (verbose>1 && !MPIrank) {
    printf("%d iter.\tSolver time ... %3.1e s\n", ITER+1, elapsed_secs);
  }
}

void CSolver::solver_hfeti(CClust_g & clust_g) { //(CData *data,CFem *fem, CDomain *domainG) {

#ifdef HFETI_SOLVER

  CData *data; 
  CFem *fem; 
  CDomain *domainG; 
  int number_of_subdomains_per_cluster = clust_g.data.size(); 

  data    = clust_g.data[0]; 
  fem     = clust_g.fem[0]; 
  domainG = clust_g.domainG; 

  int MPIrank; MPI_Comm_rank(fem->comm, &MPIrank);
  int MPIsize; MPI_Comm_size(fem->comm, &MPIsize);

  extern void SetCluster		  ( Cluster & cluster, int * subdomains_global_indices, int number_of_subdomains, int MPI_rank); 


  extern void SetMatrixB1_fromCOO ( Cluster & cluster, int domain_index_in_cluster, 
      int n_rows, int n_cols, int nnz, 
      int * I_rows, int * J_cols, double * V_vals, char type ); 

  extern void SetMatrixB0_fromCOO ( Cluster & cluster, int domain_index_in_cluster, 
	  int n_rows, int n_cols, int nnz, 
	  int * I_rows, int * J_cols, double * V_vals, char type );

  extern void SetMatrixR_fromDense( Cluster & cluster, int domain_index_in_cluster,
      int n_cols, int n_rows, double * vals, char type ); 

  extern void SetMatrixK_fromCSR ( Cluster & cluster, int domain_index_in_cluster, 
      int n_rows, int n_cols, int * rows, int * cols, double * vals, char type ); 

  extern void SetSolverPreprocessing ( Cluster & cluster, IterSolver & solver, 
      vector <vector <int> > & lambda_map_sub, vector < int > & neigh_domains ); 

  extern void SetMatrixFromCSR   ( SparseMatrix    & Mat, int n_rows, int n_cols, int * rows, int * cols, double * vals, char type );
  extern void SetMatrixFromDense ( SparseMatrix    & Mat, int n_cols, int n_rows, double * vals, char type ); 
  extern void SetMatrixFromCOO   ( SparseMatrix    & Mat, int n_rows, int n_cols, int nnz, int * I_rows, int * J_cols, double * V_vals, char type );
  extern void SetVecInt          ( vector <int>    & vec, int incerement_by, int nnz, int * vals);
  extern void SetVecDbl          ( vector <double> & vec, int nnz,	double * vals); 

  if (MPIrank == 0) {
    cout << endl; 
    cout << " ******************************************************************************************************************************* " << endl; 
    cout << " Solver initialization ... " << endl; 
    GetProcessMemoryStat ( ); GetMemoryStat( );
    cout << endl; 
  }



  Cluster cluster(MPIrank + 1); 
  cluster.USE_DYNAMIC		 = 0;
  cluster.USE_HFETI			 = 0;  
  cluster.USE_KINV			 = 0;   
  cluster.SUBDOM_PER_CLUSTER = number_of_subdomains_per_cluster;
  cluster.NUMBER_OF_CLUSTERS = MPIsize; 

  IterSolver solver;
  solver.CG_max_iter    = 50; 
  solver.USE_GGtINV		= 1; 
  solver.epsilon		= 0.00001; 
  solver.USE_HFETI		= cluster.USE_HFETI;
  solver.USE_KINV		= cluster.USE_KINV;
  solver.USE_DYNAMIC	= 0;  
  solver.USE_PIPECG		= 0; 
  solver.USE_PREC		= 0; 
  solver.FIND_SOLUTION	= 0;					

  mkl_cbwr_set(MKL_CBWR_COMPATIBLE); 


  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }


  // *** Get options from command line *****************************************
#ifndef WIN32
  
  
  optind = 0; 

  int opt= 0;

  static struct option long_options[] = {
	  {"Nx",          required_argument, 0,  'a' },
	  {"Ny",          required_argument, 0,  'b' },
	  {"Nz",          required_argument, 0,  'c' },
	  {"nx",          required_argument, 0,  'd' },
	  {"ny",          required_argument, 0,  'e' },
	  {"nz",          required_argument, 0,  'f' },
	  {"Clx",         required_argument, 0,  'g' },
	  {"Cly",         required_argument, 0,  'h' },
	  {"Clz",         required_argument, 0,  'i' },
	  {"Cox",         required_argument, 0,  'j' },
	  {"Coy",         required_argument, 0,  'k' },
	  {"Coz",         required_argument, 0,  'l' },
	  {"EDGE_COR",    no_argument,       0,  'm' },
	  {"FACE_COR",    no_argument,       0,  'n' },
	  {"VTK",         no_argument,       0,  'o' },
	  {"vtk_max_rank",required_argument, 0,  'p' },
	  {"vtk_min_rank",required_argument, 0,  'q' },
	  {"HFETI",       required_argument, 0,  'A' },
	  {"KINV",        no_argument,       0,  'B' },
	  {"PIPECG",      required_argument, 0,  'C' },
	  {"GGTINV",      required_argument, 0,  'D' },
	  {"ITERS",       required_argument, 0,  'E' },
	  {"TOL",         required_argument, 0,  'F' },
	  {0,             0,                 0,   0  }
  };			    


  int long_index =0;
  while ((opt = getopt_long(clust_g.argc_l, clust_g.argv_l,"a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:A:B:C:D:E:F", 
	  long_options, &long_index )) != -1) {
		  switch (opt) {
		  case 'a' : ;//Nxyz[0] = atoi(optarg);
			  break;
		  case 'b' : ;//Nxyz[1] = atoi(optarg);
			  break;
		  case 'c' : ;//Nxyz[2] = atoi(optarg);
			  break;
		  case 'd' : ;//nxyz[0] = atoi(optarg);
			  break;
		  case 'e' : ;//nxyz[1] = atoi(optarg); 
			  break;
		  case 'f' : ;//nxyz[2] = atoi(optarg);
			  break;
		  case 'g' : ;//Cxyz[0] = atoi(optarg);
			  break;
		  case 'h' : ;//Cxyz[1] = atoi(optarg);
			  break;
		  case 'i' : ;//Cxyz[2] = atoi(optarg);
			  break;
		  case 'j' : ;//nCorners_X = atoi(optarg);
			  break;
		  case 'k' : ;//nCorners_Y = atoi(optarg);
			  break;
		  case 'l' : ;//nCorners_Z = atoi(optarg);
			  break;
		  case 'm' : ;//flag_DP_eges  = true;
			  break;
		  case 'n' : ;//flag_DP_inner = true;
			  break;
		  case 'o' : ;//flag_store_VTK = true;
			  break;
		  case 'p' : ;//vtk_max_rank=atoi(optarg);
			  break;
		  case 'q' : ;//vtk_min_rank = atoi(optarg);
			  break;
		  case 'A' : cluster.USE_HFETI = atoi(optarg);
			  break;
		  case 'B' : cluster.USE_KINV  = 1;
			  break;
		  case 'C' : solver.USE_PIPECG  = atoi(optarg);
			  break;
		  case 'D' : solver.USE_GGtINV  = atoi(optarg);
			  break;
		  case 'E' : solver.CG_max_iter = atoi(optarg);
			  break;
		  case 'F' : solver.epsilon     = atof(optarg);
			  break;
		  default: 
			  goto loopend;
		  }
  }
  loopend:

  solver.USE_KINV		= cluster.USE_KINV;
  solver.USE_HFETI		= cluster.USE_HFETI;

#endif
  // *** END - Get options from command line ***********************************

  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... SetCluster() " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }

  //SetCluster( cluster, MPIrank + 1, 1, MPIrank); 
  SetCluster( cluster, clust_g.domainG->vec_globalSubNumbering, number_of_subdomains_per_cluster, MPIrank); 

  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... solver.Setup() " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }
  
  vector<double> solver_parameters ( 10 );
  solver.Setup ( solver_parameters, cluster );

  TimeEval solver_memory ("Solver runtime memory usage "); 
  solver_memory.totalTime.AddStartWOBarrier(0.0001); //GetProcessMemory()); 

  TimeEvent before_solver_mem ("Memory used before solver "); 
  before_solver_mem.AddStartWOBarrier(0.0001); 
  before_solver_mem.AddEndWOBarrier(GetProcessMemory());
  solver_memory.AddEvent(before_solver_mem); 

  TimeEval solver_timing ("FETI Solver overall timing ");
  solver_timing.totalTime.AddStart(omp_get_wtime());

  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... B0 setup " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }

  // *** Setup B0 matrix *******************************************************************************************
  if (cluster.USE_HFETI == 1 ) {
	  TimeEvent B0_solver_mem ("Memory used by matrix B0 "); 
	  B0_solver_mem.AddStartWOBarrier(GetProcessMemory()); 

	  TimeEvent setupB0_time ("Setup B0 matrix");
	  setupB0_time.AddStart(omp_get_wtime()); 
	  // POZOR
	  for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
		  int domain_index_in_cluster = i; 
		  SetMatrixB0_fromCOO( cluster, domain_index_in_cluster,
			  clust_g.data[i]->B->B0_rows,    // B_full_rows, //n_row_eq, 
			  clust_g.data[i]->B->B0_cols,    // B_full_cols, //n_col, 
			  clust_g.data[i]->B->B0_nnz,     // B_full_nnz,  //nnz_eq, 
			  &clust_g.data[i]->B->B0_I[0],   // BI_full[0], //Bi_coo, 
			  &clust_g.data[i]->B->B0_J[0],   // BJ_full[0], //Bj_coo, 
			  &clust_g.data[i]->B->B0_V[0],   // BV_full[0], //Bv_coo, 
			  'G' ); 
	  }
	  setupB0_time.AddEnd(omp_get_wtime());
	  solver_timing.AddEvent(setupB0_time);

	  B0_solver_mem.AddEndWOBarrier(GetProcessMemory());
	  solver_memory.AddEvent(B0_solver_mem); 
  }
  // *** END - Setup B0 matrix *************************************************************************************

  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... B1 setup " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }

  // *** Setup B1 matrix *******************************************************************************************
  TimeEvent B1_solver_mem ("Memory used by matrix B1 "); 
  B1_solver_mem.AddStartWOBarrier(GetProcessMemory()); 

  TimeEvent setupB1_time ("Setup B1 matrix");
  setupB1_time.AddStart(omp_get_wtime()); 
  cilk_for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
	int domain_index_in_cluster = i; 
	SetMatrixB1_fromCOO( cluster, domain_index_in_cluster,
      clust_g.data[i]->B->B_full_rows, //n_row_eq, 
	  clust_g.data[i]->B->B_full_cols, //n_col, 
	  clust_g.data[i]->B->B_full_nnz,  //nnz_eq, 
      &clust_g.data[i]->B->BI_full[0], //Bi_coo, 
	  &clust_g.data[i]->B->BJ_full[0], //Bj_coo, 
	  &clust_g.data[i]->B->BV_full[0], //Bv_coo, 
	  'G' );
  }
  setupB1_time.AddEnd(omp_get_wtime());
  solver_timing.AddEvent(setupB1_time);

  B1_solver_mem.AddEndWOBarrier(GetProcessMemory());
  solver_memory.AddEvent(B1_solver_mem); 
  // *** END - Setup B1 matrix *************************************************************************************

  if (MPIrank == 0) {
	  cout << endl; 
	  cout << " ******************************************************************************************************************************* " << endl; 
	  cout << " Solver initialization ... " << endl; 
	  GetProcessMemoryStat ( ); GetMemoryStat( );
	  cout << endl; 
  }

  // *** Setup R matrix ********************************************************************************************
  TimeEvent R_solver_mem ("Memory used by matrix R "); 
  R_solver_mem.AddStartWOBarrier(GetProcessMemory()); 

  TimeEvent setupR_time ("Setup R matrix");
  setupR_time.AddStart(omp_get_wtime());
  cilk_for(int i = 0; i < number_of_subdomains_per_cluster; i++) {
	CRSparse R(clust_g.fem[i],domainG->neqSub[clust_g.fem[i]->i_domOnClust]);
	int domain_index_in_cluster = i;
	SetMatrixR_fromDense( cluster, domain_index_in_cluster, R.n_col, R.n_row, R.Rfull, 'G' );
  }
  setupR_time.AddEnd(omp_get_wtime());
  solver_timing.AddEvent(setupR_time);

  R_solver_mem.AddEndWOBarrier(GetProcessMemory());
  solver_memory.AddEvent(R_solver_mem); 
  // *** END - Setup R matrix **************************************************************************************

  // *** Load RHS and fix points for K regularization **************************************************************
  if (MPIrank == 0) {
    cout << " ******************************************************************************************************************************* " << endl; 
    cout << " *** Load RHS and fix points for K regularization ******************************************************************************* " << endl; 
  }
  if (MPIrank == 0) {
    cout << " Solver - Loading f and fix DOFs ... " << endl; 
    GetProcessMemoryStat ( );
    GetMemoryStat( );
  }
    
  cilk_for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
	  //fem  = clust_g.fem[i];
	  //data = clust_g.data[i];
	  //domain_index_in_cluster = i;
	  SetVecDbl( cluster.domains[i].f,        clust_g.data[i]->KSparse->n_row, clust_g.data[i]->fE ); 
      SetVecInt( cluster.domains[i].fix_dofs, 1,                           24, clust_g.fem[i]->mesh.fixingDOFs ); 

	  //for (int j = 0; j < cluster.domains[i].f.size(); j++)
		 // cluster.domains[i].f[j] += rand()

  }

  if (MPIrank == 0) {
    cout << " Solver - before 'stif_glob_number' deleting ... " << endl; 
    GetProcessMemoryStat ( ); GetMemoryStat( );
  }

  if (MPIrank == 0) {
    cout << " *** END - Load RHS and fix points for K regularization ************************************************************************ " << endl; 
    cout << " ******************************************************************************************************************************* " << endl; 
    cout << endl; 
  }
  // *** END - Load RHS and fix points for K regularization ********************************************************

  bool eraseData = !domainG->flag_store_VTK && gv_flag_linear_system;

  if (MPIrank == 0) {
    cout << " Solver - before 'stif_glob_number' deleting ... " << endl; 
    GetProcessMemoryStat ( ); GetMemoryStat( );
  }

  if (eraseData){
    if (data->stif_glob_number!=NULL) { delete [] data->stif_glob_number; data->stif_glob_number=NULL; }
  }
  if (MPIrank == 0) {
    cout << " Solver - after 'stif_glob_number' deleting ... " << endl; 
    GetProcessMemoryStat ( ); GetMemoryStat( ); 
  }

  // *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 ******************* 
  TimeEvent prec_solver_mem ("Memory used by solver - preprocessing "); 
  prec_solver_mem.AddStartWOBarrier(GetProcessMemory()); 

  TimeEvent setupSol_time ("Setup FETI solver - preprocessing");
  setupSol_time.AddStart(omp_get_wtime());
   
  cilk_for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
	  cluster.domains[i].lambda_map_sub = clust_g.data[i]->B->lambda_map_sub;
  }
  
  SetSolverPreprocessing ( cluster, solver, 
						   //clust_g.data[0]->B->lambda_map_sub, 
						   clust_g.domainG->lambda_map_sub,
						   clust_g.data[0]->B->neigh_clusters ); 
  setupSol_time.AddEnd(omp_get_wtime());
  solver_timing.AddEvent(setupSol_time);

  prec_solver_mem.AddEndWOBarrier(GetProcessMemory());
  solver_memory.AddEvent(prec_solver_mem); 
  // *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************

  // *** Load Matrix K and regularization ******************************************************************************
  TimeEvent K_solver_mem ("Memory used by Mat K factorization "); 
  K_solver_mem.AddStartWOBarrier(GetProcessMemory()); 

  // POZOR !! - zakomentoval jsem 
  
  //if (eraseData ){
  //  if (fem->bound_cond != NULL)          { delete    fem->bound_cond;          fem->bound_cond = NULL;}
  //  if (data->fE!=NULL)                 { delete [] data->fE;                 data->fE=NULL;}
  //  if (data->fGlobal!=NULL)            { delete [] data->fGlobal;            data->fGlobal=NULL;}
  //  if (data->Ku!=NULL)                 { delete [] data->Ku;                 data->Ku=NULL;}
  //  if (data->f_subdom!=NULL)           { delete [] data->f_subdom;           data->f_subdom=NULL;}
  //  if (data->f_subdom_internal!=NULL)  { delete [] data->f_subdom_internal;  data->f_subdom_internal=NULL;}
  //  if (data->indExterDOFs!=NULL)       { delete [] data->indExterDOFs;       data->indExterDOFs=NULL;}
  //  if (fem->mesh.element!=NULL)       { delete [] fem->mesh.element;       fem->mesh.element=NULL;}
  //  if (fem->mesh.coordinateSub!=NULL) { delete [] fem->mesh.coordinateSub; fem->mesh.coordinateSub=NULL;}
  //  if (fem->mesh.faceSub!=NULL)       { delete [] fem->mesh.faceSub;       fem->mesh.faceSub=NULL;}
  //  if (fem->mesh.nodOnEdg!=NULL)      { delete [] fem->mesh.nodOnEdg;      fem->mesh.nodOnEdg=NULL;}
  //  if (fem->mesh.neighbSub!=NULL)     { delete [] fem->mesh.neighbSub;     fem->mesh.neighbSub=NULL;}
  //  if (fem->mesh.fixingDOFs!=NULL)    { delete [] fem->mesh.fixingDOFs;    fem->mesh.fixingDOFs=NULL;}
  //}

  // *** Load Matrix K and regularization ******************************************************************************
  TimeEvent setupK_time ("K regularization and factorization");
  setupK_time.AddStart(omp_get_wtime());
  if (MPIrank == 0)
	  cout << "K import and factorization: " << endl; 
  
  MKL_Set_Num_Threads(1);
  cilk_for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
	  //data = clust_g.data[i]; 
	  //domain_index_in_cluster = i;
	  SetMatrixK_fromCSR ( cluster, i, 
		clust_g.data[i]->KSparse->n_row,   clust_g.data[i]->KSparse->n_row, 
		clust_g.data[i]->KSparse->row_ptr, clust_g.data[i]->KSparse->col_ind, clust_g.data[i]->KSparse->val, 'S'); 
	   
	  if (MPIrank == 0 ) { //|| MPIrank == 1) {
		  cout << i << " ";
		  //cout << MPIrank << " : " << i << " "; 
		  //cout << endl;
	  }
  }
  if (MPIrank == 0)
	  cout << endl; 

  if (cluster.USE_KINV == 1)
	  cluster.Create_Kinv_perDomain(); 
  
  setupK_time.AddEnd(omp_get_wtime());
  solver_timing.AddEvent(setupK_time);

  K_solver_mem.AddEndWOBarrier(GetProcessMemory());
  solver_memory.AddEvent(K_solver_mem); 

  setupK_time.PrintStatMPI(0.0);

  if (cluster.USE_HFETI == 1) {
	  TimeEvent HFETI_setup_time ("Setup Hybrid FETI solver - preprocessing"); HFETI_setup_time.AddStart(omp_get_wtime());
	  TimeEvent HFETI_setup_mem ("Setup Hybrid FETI solver - preprocessing");  HFETI_setup_mem.AddStartWOBarrier(GetProcessMemory());

	  cluster.SetClusterHFETI(); 

	  HFETI_setup_time.AddEnd(omp_get_wtime());            solver_timing.AddEvent(HFETI_setup_time);
	  HFETI_setup_mem.AddEndWOBarrier(GetProcessMemory()); solver_memory.AddEvent(HFETI_setup_mem);

	  HFETI_setup_mem.PrintStatMPI(0.0); 
  }


  cluster.SetClusterPC_AfterKplus(); 


  // *** END - Load Matrix K and regularization  ***********************************************************************


  if (eraseData){
    if (data->KSparse!=NULL)            { delete    data->KSparse; data->KSparse=NULL; }
//    if (data->B!=NULL)                  { delete    data->B;                  data->B=NULL;}
  }



  // *** Running Solver ************************************************************************************************
  TimeEvent solverRun_time ("FETI solver - runtime");
  solverRun_time.AddStart(omp_get_wtime());
  if (MPIrank == 0) {
    cout << " ******************************************************************************************************************************* " << endl; 
    cout << " *** Running Solver ************************************************************************************************************ " << endl; 
  }

  string result_file("MATSOL_SVN_Displacement.Nvec");
  if (MPIrank == 0) {
    cout << " Solver - Running Solve ... " << endl; 
    GetProcessMemoryStat ( );
    GetMemoryStat( ); }

  solver.Solve_singular ( cluster, result_file );

  if (MPIrank == 0) {
    cout << " Solver - Solution copy ... " << endl; 
    GetProcessMemoryStat ( ); GetMemoryStat( ); 
  }

  if (MPIrank == 0) {
    GetProcessMemoryStat ( ); GetMemoryStat( ); 
  }

  if (MPIrank == 0) {
    cout << " *** END - Running Solver ****************************************************************************************************** " << endl; 
    cout << " ******************************************************************************************************************************* " << endl; 
    cout << endl; 
  }
  solverRun_time.AddEnd(omp_get_wtime());
  solver_timing.AddEvent(solverRun_time);

  solver_timing.totalTime.AddEnd(omp_get_wtime());
  solver_timing.PrintStatsMPI(); 

  solver_memory.totalTime.AddEndWOBarrier(GetProcessMemory()); 
  solver_memory.PrintStatsMPI(); 




  vector < vector < double > > prim_solution; 
  solver.GetSolution_Primal_singular_parallel(cluster, prim_solution);
  double max_v = 0.0; 

  for (int i = 0; i < number_of_subdomains_per_cluster; i++) 
	  for (int j = 0; j < prim_solution[i].size(); j++) 
		  if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] ); 
	  
  TimeEvent max_sol_ev ("Max solution value"); max_sol_ev.AddStartWOBarrier(0.0); max_sol_ev.AddEndWOBarrier(max_v);            
  
  //printf("Rank %d - max solution value = %f \n", MPIrank, max_v); 
  double max_vg;  
  MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD ); 
  if (MPIrank == 0)
	  cout << " Max value in_solution = " << max_vg << endl; 
  
  max_sol_ev.PrintLastStatMPI_PerNode(max_vg);


   

  if (clust_g.domainG->flag_store_VTK) 
  {
	   for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
		  for (int j = 0; j < prim_solution[i].size(); j++) {
			  if (prim_solution[i][j] > max_v) max_v = prim_solution[i][j]; 
		  }
		  copy(prim_solution[i].begin(), prim_solution[i].end(), clust_g.data[i]->ddu); 
	  }

  }


  // *** END - Running Solver ************************************************************************************************

#endif
}



#ifdef FLLOP_ENABLED
PetscErrorCode ierr;
#define FLLTRY(f) ierr = f; if (ierr) throw ierr;

void CSolver::solver_fllop(CData *data,CFem *fem,CDomain *domainG) {



//printf("------------------------------------B->n_row_eq=%d\n",data->B->n_row_eq);
//for (int i = 0;i<data->B->i_eq[data->B->n_row_eq];i++){
//  printf("[j,v]=[%d,%3.1e]\n",data->B->j_eq[i],data->B->v_eq[i]);
//}


  PetscLogStagePush(solve_stage);
  int MPIrank;
  MPI_Comm_rank(fem->comm, &MPIrank);
  //CRSparse R(fem);
  CRSparse R(fem,domainG->neqSub[fem->i_domOnClust]);
  CBSparse * B = data->B;
  int n = data->KSparse->n_row;
  int d = R.n_col;

  FLLTRY( FllopAIFSetFETIOperator(n, data->KSparse->row_ptr, data->KSparse->col_ind, data->KSparse->val, QP_SYM_UPPER_TRIANGULAR) );
  FLLTRY( FllopAIFSetFETIOperatorNullspaceDense(n, d, R.Rfull) );
  FLLTRY( FllopAIFSetRhs(n, data->fE) );
  FLLTRY( FllopAIFSetSolutionVector(n, data->ddu) );
  //FLLTRY( FllopAIFSetRestrictedSolutionVector(n, u_restrict) );
  int *l2g = &fem->mesh.l2g[0];

#ifdef FLLOP_NEW_API
  /* test new FLLOP API for FETI gluing and Dirichlets */
  if (domainG->verbose>0 && !MPIrank && domainG->i_load_step == 0 && domainG->i_sub_step == -1) {
    cout << "FLLOP_NEW_API  activated" << endl;
  }
  int n_db = fem->bound_cond->dirBCSub->n; //number of dir bc
  QP aif_qp;
  IS l2g_dof_map, neighbor_subdomains, dbcis;
  FLLTRY( FllopAIFGetQP(&aif_qp) );
  FLLTRY( ISCreateGeneral(fem->comm, n, l2g, PETSC_USE_POINTER, &l2g_dof_map) );
  FLLTRY( ISCreateGeneral(fem->comm, domainG->n_neighbSub, fem->mesh.neighbSub, PETSC_USE_POINTER, &neighbor_subdomains) );
  FLLTRY( ISCreateGeneral(fem->comm, n_db, n_db?&fem->bound_cond->dirBCSub->ind[0]:NULL, PETSC_OWN_POINTER, &dbcis) ); 
  FLLTRY( QPFetiSetLocalToUndecomposedGlobalMapping(aif_qp, l2g_dof_map) );
  FLLTRY( QPFetiSetNeighborRanks(aif_qp, neighbor_subdomains) );
  FLLTRY( QPFetiAddDirichlet(aif_qp, dbcis, FETI_LOCAL, PETSC_FALSE) );
  FLLTRY( QPFetiSetUp(aif_qp) );
  FLLTRY( ISDestroy(&l2g_dof_map) );
  FLLTRY( ISDestroy(&neighbor_subdomains) );
  FLLTRY( ISDestroy(&dbcis) );
#else /* FLLOP_NEW_API */

#ifdef FLLOP_NEW_SEPARATED_B
  if (!MPIrank) cout << "FLLOP_NEW_SEPARATED_B  activated" << endl; 
	int *ibd = new int[B->n_row_bd+1];
	int *ibg = new int[B->n_row_bg+1];
	CBSparse::coo2crs(B->n_row_bd, B->nnz_bd, B->i_bd, ibd);
	CBSparse::coo2crs(B->n_row_bg, B->nnz_bg, B->i_bg, ibg);
  FLLTRY( FllopAIFAddEq(B->n_row_bg, B->n_col, PETSC_FALSE, PETSC_TRUE, ibg, B->j_bg, B->v_bg, NULL) );
  FLLTRY( FllopAIFAddEq(B->n_row_bd, B->n_col, PETSC_FALSE, PETSC_TRUE, ibd, B->j_bd, B->v_bd, NULL) );
#else  /* FLLOP_NEW_SEPARATED_B */
  FLLTRY( FllopAIFSetEq(B->n_row_eq, B->n_col, PETSC_FALSE, PETSC_TRUE, B->i_eq, B->j_eq, B->v_eq, NULL) );
#endif /* FLLOP_NEW_SEPARATED_B */

#endif /* FLLOP_NEW_API */

  if (domainG->flag_contact){
    FLLTRY( FllopAIFSetIneq(B->n_row_bc, B->n_col, PETSC_FALSE, PETSC_TRUE, B->i_bc_crs, B->j_bc, B->v_bc, NULL) );
  }

#ifdef FLLOP_NEW_PRECONDITIONING
  FLLTRY( FllopAIFSetPCScalingEdgeMultiplicity(data->B->nnz_eq, data->B->i_multiplicityDD, data->B->multiplicityDD) );
#endif
  if (domainG->verbose>2 && !MPIrank) cout << "      Calling FllopAIFFetiPrepare()" << endl;
  FLLTRY( FllopAIFFetiPrepare(PETSC_TRUE) );

  if (domainG->verbose>2 && !MPIrank) cout << "      Calling FllopAIFSolve()" << endl;
  FLLTRY( FllopAIFSolve() );

  QPS qps;
  int iter;
  KSPConvergedReason reason;
  FLLTRY( FllopAIFGetQPS(&qps) );
  FLLTRY( QPSGetIterationNumber(qps,&iter) );
  FLLTRY( QPSGetConvergedReason(qps,&reason) );
  if (domainG->verbose>2 && !MPIrank) cout << "    FLLOP " 
    << ((reason>0)?"converged":"diverged") << " in " << iter << " iterations" << endl;

  FLLTRY( FllopAIFReset() );

#ifdef FLLOP_NEW_SEPARATED_B
	delete [] ibd;
	delete [] ibg;
#endif
  PetscLogStagePop();
}
#undef FLLTRY
#endif

void CSolver::solving(CFem *fem, CDomain *domainG, CData *data, CClust_g & clust) {
  int MPIrank;
  MPI_Comm_rank(fem->comm, &MPIrank);

#ifndef WIN32
  if (MPIrank == 0) { cout << " Solving start                                                            "; system("date +%T.%6N"); }
#endif

#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << endl; 
  //  cout << " in 'solving', beginning ... " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //  cout << endl; 
  //}
#endif

//  CDomain *domainG = domainG;

  //Prepare deltaVec
  double *deltaVec = domainG->deltaVec;
  double del_deltaVec;
  int flagSubstepConverget=1, numbInIt=0;
  double norm_ddu=0.0, norm_dduGlob=0.0;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  double stressStrein[domainG->n_load_steps*2];
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  for (int i = 0; i < domainG->n_load_steps; i++) {
    deltaVec[i] = ((double) (i + 1)) / domainG->n_load_steps;
  }

  double it1;
  /*	for (int i = 0; i < domainG->n_load_steps; i++) {
      it1 = 2.0*asin(sin(1.5*CONST_PI*(deltaVec[i])))/CONST_PI;
      deltaVec[i] = it1*(pow(1.0+0.2*deltaVec[i],2.0)) + deltaVec[i]*0.1;
      printf("%3.3e\n",deltaVec[i] );
      }*/
  // 
//  CSolid45EqNumbGlobal *stif_glob_number = data->stif_glob_number;

#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << endl; 
  //  cout << " in 'solving', after 'stif_glob_number' ... " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //  cout << endl; 
  //}
#endif
  // *********************************************************
  //							   quasi time iteration 
  // *********************************************************
  for (int i_load_step = 0; i_load_step < domainG->n_load_steps;
      i_load_step++) {
    //
    domainG->i_load_step = i_load_step;
    domainG->i_sub_step = -1;
    /*                                                                        */
    del_deltaVec = deltaVec[i_load_step];    //
    if (i_load_step > 0) {
      del_deltaVec = del_deltaVec - deltaVec[i_load_step - 1];
    } //(i_load_step>0)
    //
    domainG->del_deltaVec = del_deltaVec;
    //
    // *** solving A*ddx = b, first increment of increment ***
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << endl; 
    //  cout << " in 'solving', before 'definitionAndSolving .'.. " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //  cout << endl; 
    //}
#endif

#ifndef WIN32
	if (MPIrank == 0) { cout << " Definition and solving start                                             "; system("date +%T.%6N"); }
#endif
	
	definitionAndSolving(data, fem, domainG, clust); //->duu
    for (int v=0;v<domainG->n_subdomOnClust;v++){
      memcpy(clust.data[v]->du, clust.data[v]->ddu, clust.domainG->neqSub[v]*sizeof(double));
    }

#ifndef WIN32
	if (MPIrank == 0) { cout << " Definition and solving - end                                             "; system("date +%T.%6N"); }
#endif

    //
    // *********************************************************
    //							    nonlinear iteration - START
    // *********************************************************
    for (int i_sub_step = 0; i_sub_step < domainG->n_sub_steps;
        i_sub_step++) {
      // 
      domainG->i_sub_step = i_sub_step;
      // *** solving A*ddx = b, increment of displacements increment  updating 
      definitionAndSolving(data, fem, domainG, clust);
      //
      norm_ddu = CLinearAlgebra::norm_v(data->ddu, domainG->neqSub[fem->i_domOnClust]);
      norm_ddu *= norm_ddu;
      MPI_Allreduce(&norm_ddu,&norm_dduGlob, 1,
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      norm_dduGlob = sqrt(norm_dduGlob);

      if (!MPIrank && domainG->verbose>1) {
        printf("  substep %d abs err %3.3e\n",
            i_sub_step+1,norm_dduGlob);
      }
      flagSubstepConverget=0;
      numbInIt = i_sub_step;
      if (norm_dduGlob < domainG->eps1) {
        flagSubstepConverget=1;
        break;
      } //
      // du = du + ddu;  disilacements increment updating
      for (int v=0;v<clust.domainG->n_subdomOnClust;v++){
        CLinearAlgebra::add_vec_a2b(clust.data[v]->ddu,
                                    clust.data[v]->du,1.0, 1.0,
                                    clust.domainG->neqSub[v]);
      }
////      CLinearAlgebra::add_vec_a2b(data->ddu,
////                                  data->du,1.0, 1.0,
////                                  domainG->neqSub[fem->i_domOnClust]);
    } // loop sub step
    // *********************************************************
    //							    nonlinear iteration - END
    // *********************************************************
    // copy temporery data on elements calculated during subs iterations 
    // (tmp_stress->stress, tmp_epel->epel, etc ...)
    //
    data->copyDataToNextIter(domainG->n_elementsSub[fem->i_domOnClust]);
    //
    if (!MPIrank && domainG->verbose>0) {
      if (flagSubstepConverget==1){
        printf("load step [%d/%d] converged after %d substeps with abs err %3.3e\n",
            i_load_step+1,domainG->n_load_steps,numbInIt+1,norm_dduGlob);
      }
      else{
        printf("load step [%d/%d] did NOT converge after %d substeps with abs err %3.3e\n",
            i_load_step+1,domainG->n_load_steps,numbInIt+1,norm_dduGlob);
      }
    }       
    // u = u + du; displacement updating
////    CLinearAlgebra::add_vec_a2b(data->du,
////                                data->u,1.0,1.0,
////                                domainG->neqSub[fem->i_domOnClust]);
    for (int v=0;v<clust.domainG->n_subdomOnClust;v++){
      CLinearAlgebra::add_vec_a2b(clust.data[v]->du,
                                  clust.data[v]->u,1.0, 1.0,
                                  clust.domainG->neqSub[v]);
    }

    //
    if (!gv_flag_linear_system){
      int iel = domainG->n_elementsSub[fem->i_domOnClust]-1;//124
      int ipp = 2+3*6;
      stressStrein[2*i_load_step] = data->stif_glob_number[iel].stif_loc->epel[ipp] +
        data->stif_glob_number[iel].stif_loc->eppl[ipp]; 
      //    stressStrein[2*i_load_step] = data->u[23]; 
      stressStrein[2*i_load_step+1] = data->stif_glob_number[iel].stif_loc->stress[ipp];
    }

  } // loop load step

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (!gv_flag_linear_system){
    char filename[128];
    sprintf(filename, "data/strainStress_%d.txt",MPIrank);
    ofstream f_001(filename);
    for (int i = 0; i < domainG->n_load_steps; i++) {
      f_001 <<stressStrein[2*i] <<" "<<stressStrein[2*i+1] << endl;
    }
    f_001.close();
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef WIN32
  if (MPIrank == 0) { cout << " Solving - end                                                            "; system("date +%T.%6N"); }
#endif

}
