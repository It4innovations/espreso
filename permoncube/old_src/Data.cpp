#include "Data.h"
#ifdef HFETI_SOLVER
#include "TimeEval.h"
#endif

CData::CData(MPI_Comm comm, int _i_domOnClust) :
  comm(comm),
  KSparse(NULL),
  B(NULL),
  stif_glob_number(NULL),
//  u_restrict(NULL),
  u(NULL),
  du(NULL),
  ddu(NULL),
  f(NULL),
  fE(NULL),
  fGlobal(NULL),
  Ku(NULL),
  f_subdom(NULL),
  f_subdom_internal(NULL)
{
}

CData::~CData() {

  if (KSparse)            {delete KSparse;              KSparse=NULL;}
  if (B)                  {delete B;                    B=NULL;}
  if (stif_glob_number)   {delete[] stif_glob_number;   stif_glob_number=NULL; }
  if (u)                  {delete[] u;                  u=NULL;}
  if (du)                 {delete[] du;                 du=NULL;}
  if (ddu)                {delete[] ddu;                ddu=NULL;}
  if (fE)                 {delete[] fE;                 fE=NULL;}
  if (fGlobal)            {delete[] fGlobal;            fGlobal=NULL;}
  if (Ku)                 {delete[] Ku;                 Ku=NULL;}
  if (f_subdom)           {delete[] f_subdom;           f_subdom=NULL;}
//  if (u_restrict)         {delete[] u_restrict;         u_restrict=NULL;}
  if (f_subdom_internal)  {delete[] f_subdom_internal;  f_subdom_internal=NULL;}
  if (indExterDOFs)       {delete[] indExterDOFs;       indExterDOFs=NULL;}
}



#ifdef HFETI_SOLVER
extern void GetMemoryStat	   ( ); 
extern void GetProcessMemoryStat ( ); 
extern double GetProcessMemory ( ); 
#endif


void CData::initialize(MPI_Comm comm, int neqSub, int neqAll, int max_nDOFs_u_is_stored, 
    int n_elementsSub, int n_indExterDOFs)
{
  int MPIrank;
  int MPIsize;
  MPI_Comm_rank(comm, &MPIrank);
  MPI_Comm_size(comm, &MPIsize);

#ifdef HFETI_SOLVER
//  if (MPIrank == 0) {
//    cout << " in 'Data' before KSparse init " << endl; 
//    GetProcessMemoryStat ( ); GetMemoryStat( );
//  }
#endif

  KSparse           = new CKSparse(comm);
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'Data' after KSparse init " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif
  B                 = new CBSparse(comm);
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'Data' before stif_glob_number " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif
  stif_glob_number  = new CSolid45EqNumbGlobal[n_elementsSub];
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'Data' after stif_glob_number " << endl; 
  //  GetProcessMemoryStat ( );
  //  GetMemoryStat( );
  //}
#endif
  u                 = new double[neqSub];
  du                = new double[neqSub];
  ddu               = new double[neqSub];
  f_subdom          = new double[neqSub];
  f_subdom_internal = new double[neqSub];
  fE                = new double[neqSub];
  Ku                = new double[neqSub];
#if defined(FLLOP_ENABLED) || defined(flag_Fortran)
  //u_restrict      = new double[neqAll];
#else
  if (neqAll<max_nDOFs_u_is_stored){
    //u_restrict      = new double[neqAll];
  }
#endif
  indExterDOFs      = new longInt[n_indExterDOFs];

  memset(u,0,neqSub * sizeof(double));
  memset(du,0,neqSub * sizeof(double));
  memset(ddu,0,neqSub * sizeof(double));
  memset(fE,0,neqSub * sizeof(double));
  memset(Ku,0,neqSub * sizeof(double));
  memset(f_subdom,0,neqSub * sizeof(double)); 
  memset(f_subdom_internal,0,neqSub * sizeof(double));
#if defined(FLLOP_ENABLED) || defined(flag_Fortran)
//  memset(u_restrict,0,neqAll * sizeof(double));
#endif
}

void CData::extDofs(CFem *fem,CDomain *domainG)
{
  for (int i = 0;i<domainG->n_nodOnEdg[fem->i_domOnClust];i++){
    this->indExterDOFs[3*i]   =	(longInt) 3*fem->mesh.nodOnEdg[i];
    this->indExterDOFs[3*i+1] = (longInt) 3*fem->mesh.nodOnEdg[i]+1;
    this->indExterDOFs[3*i+2] = (longInt) 3*fem->mesh.nodOnEdg[i]+2;
  }
  
//  extDofs_l2g(domainG->n_exterDOFs[fem->i_domOnClust],fem->mesh.l2g);
}

void CData::extDofs_l2g(int nExtDofs,std::vector<longInt> l2g){
  for (int i = 0;i<nExtDofs;i++){
    this->indExterDOFs[i] = l2g[this->indExterDOFs[i]];	
  }
}

void CData::copyDataToNextIter(int n_elementsSub) {
  {
    if (!gv_flag_linear_system)
      for (int kk=0;kk<n_elementsSub;kk++){
        memcpy( this->stif_glob_number[kk].stif_loc->epel,
            this->stif_glob_number[kk].stif_loc->tmp_epel,
            8*6*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->eppl,
            this->stif_glob_number[kk].stif_loc->tmp_eppl,
            8*6*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->stress,
            this->stif_glob_number[kk].stif_loc->tmp_stress,
            8*6*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->currentIsotropicHardening,
            this->stif_glob_number[kk].stif_loc->tmp_currentIsotropicHardening,
            8*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->statev,
            this->stif_glob_number[kk].stif_loc->tmp_statev,
            8*6*6*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->plwork,
            this->stif_glob_number[kk].stif_loc->tmp_plwork,
            8*sizeof(double));
        memcpy( this->stif_glob_number[kk].stif_loc->epeq,
            this->stif_glob_number[kk].stif_loc->tmp_epeq,
            8*sizeof(double));
      }
  }
}
