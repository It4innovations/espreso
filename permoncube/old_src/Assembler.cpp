
#include "Assembler.h"
#include "LinearAlgebra.h"
#ifdef HFETI_SOLVER
#include "TimeEval.h"
#endif

static int event1=0;
#ifdef HFETI_SOLVER
	extern void GetMemoryStat	   ( ); 
	extern void GetProcessMemoryStat ( ); 
	extern double GetProcessMemory ( ); 
#endif

void CAssembler::stima_subdomain(double *f_subdom, double *f_subdom_internal,
        double * u, double *du, CSolid45EqNumbGlobal *stif_glob_number,
        CFem * fem, CDomain *domainG, CKSparse * Ksparse) {
  if (!event1) LogEventRegister("stima_subdomain",&event1);
  LogEventBegin(event1);
  
	int MPIrank;
  MPI_Comm_rank(fem->comm, &MPIrank);
  
	CCoordinate *coordinatesLocal = new CCoordinate[8];
	double *C = new double[36];
	double *gaussPoints = new double[24];
	double e = domainG->youngsModulus;
	double mu = domainG->poissonsRatio;
	double const1 = e / ((1. + mu) * (1. - 2. * mu));
	double mu2 = 1. - mu;
	double mu3 = 0.5 - mu;
	double pt = 0.577350269189626;
  int    i_sub_step = domainG->i_sub_step;
	//
  memset(&(gaussPoints[0]),0, 24 * sizeof(double));
  memset(&(C[0]),0, 36* sizeof(double));
	//
	C[0] = mu2 * const1;
	C[6] = mu * const1;
	C[12] = mu * const1;
	C[1] = mu * const1;
	C[7] = mu2 * const1;
	C[13] = mu * const1;
	C[2] = mu * const1;
	C[8] = mu * const1;
	C[14] = mu2 * const1;
	C[21] = mu3 * const1;
	C[28] = mu3 * const1;
	C[21] = mu3 * const1;
	C[35] = mu3 * const1;
	//
	gaussPoints[0] = -pt; // 0
	gaussPoints[1] = -pt;
	gaussPoints[2] = -pt;
	gaussPoints[3] = -pt;
	gaussPoints[4] = pt;
	gaussPoints[5] = pt;
	gaussPoints[6] = pt;
	gaussPoints[7] = pt;
	gaussPoints[8] = -pt; // 1
	gaussPoints[9] = -pt;
	gaussPoints[10] = pt;
	gaussPoints[11] = pt;
	gaussPoints[12] = -pt;
	gaussPoints[13] = -pt;
	gaussPoints[14] = pt;
	gaussPoints[15] = pt;
	gaussPoints[16] = -pt; // 2
	gaussPoints[17] = pt;
	gaussPoints[18] = -pt;
	gaussPoints[19] = pt;
	gaussPoints[20] = -pt;
	gaussPoints[21] = pt;
	gaussPoints[22] = -pt;
	gaussPoints[23] = pt;
	//
	//
	//
	double howMuchPerc = 0.0;
	int* ieq_glob_fortran = new int[24];
	//
	clock_t begin = clock();
	//
#ifdef HFETI_SOLVER
	//if (MPIrank == 0) {
	//	cout << " in 'Assembler', before 'fe_symbo.'.. " << endl; 
	//	GetProcessMemoryStat ( );
	//	GetMemoryStat( );
	//}
#endif

  if (domainG->i_load_step == 0 && domainG->i_sub_step == -1) {
	  Ksparse->fe_symbolic_form_matrix(fem,domainG, stif_glob_number);
  }
  else {
 #ifndef flag_Fortran
    int nOfEl = sizeof(double) * Ksparse->row_ptr[Ksparse->n_row];
	  memset(Ksparse->val, 0, nOfEl);
 #endif
  }
#ifdef HFETI_SOLVER
	//if (MPIrank == 0) {
	//	cout << " in 'Assembler', after 'fe_symbo.'.. " << endl; 
	//	GetProcessMemoryStat ( );
	//	GetMemoryStat( );
	//}
#endif
  // right hand side vector
  memset(f_subdom,0, domainG->neqSub[fem->i_domOnClust] * sizeof(double));
  memset(f_subdom_internal,0, domainG->neqSub[fem->i_domOnClust]  * sizeof(double));
	
  int	  iLocalNumb;
  longInt iGlobalNumb;
	for (int i = 0, cnt=0; i < domainG->n_elementsSub[fem->i_domOnClust] ; cnt++,i++) {
    CStiffnessLocal *stif_loc = stif_glob_number[i].stif_loc;
    for (int j = 0; j < 8; j++) {
      iGlobalNumb						= fem->mesh.element[i].inod_glob[j];
      iLocalNumb						= fem->mesh.g2l_nodes[iGlobalNumb];
      coordinatesLocal[j]		= fem->mesh.coordinateSub[iLocalNumb];
    }
    for (int j=0;j<24;j++){
      stif_loc->value_u[j] = u[stif_glob_number[cnt].ieq[j]];
      stif_loc->value_du[j] = du[stif_glob_number[cnt].ieq[j]];
    }
//
    CElementLibrary::stima_solid45(i_sub_step,e,mu,
																stif_loc, coordinatesLocal, gaussPoints, C,
  	                            domainG->acceleration,domainG->density);
#ifdef flag_Fortran
		for (int j = 0; j < 24; j++) {
			ieq_glob_fortran[j] = fem->mesh.l2g[stif_glob_number[cnt].ieq[j]] + 1;
		}
		fe2feti_numeric_map_element_(24, &(ieq_glob_fortran[0]), &(stif_loc->value_K[0]));
#else
    Ksparse->fe_numeric_form_matrix(&stif_glob_number[cnt], stif_loc,fem);
#endif
    for (int j = 0; j < 24; j++) {
      f_subdom[stif_glob_number[cnt].ieq[j]] 	
                       += stif_loc->value_f[j];
      f_subdom_internal[stif_glob_number[cnt].ieq[j]] 	
                       += stif_loc->value_f_internal[j];
    }
  
  }	 //
#ifdef flag_Fortran
		fe2feti_numeric_factorize_();
#endif
	//
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	if (!MPIrank & false) {
		printf("K, f .... %3.1e s\n", elapsed_secs);
	}
  

	//
	delete[] gaussPoints; 				gaussPoints = NULL; 
	delete[] C; 									C = NULL;
	delete[] ieq_glob_fortran;    ieq_glob_fortran = NULL;
	delete[] coordinatesLocal;    coordinatesLocal = NULL;
	//
	if (!MPIrank & false) {
		printf("stiffness subdomain: \tdone\n");
	}
  LogEventEnd(event1);
}
