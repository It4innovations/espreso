
#include "KSparse.h"
#ifdef HFETI_SOLVER
#include "TimeEval.h"
#endif

static int event1=0, event2=0;




#ifdef HFETI_SOLVER
extern void GetMemoryStat	   ( ); 
extern void GetProcessMemoryStat ( ); 
extern double GetProcessMemory ( ); 
#endif




CKSparse::CKSparse(MPI_Comm comm)
{
  this->comm = comm;
  val  		= NULL;
  col_ind = NULL;
  row_ptr = NULL;
  n_row   = 0;
}

CKSparse::~CKSparse() {
  if (val) 	  	 {delete [] val;     val=NULL;        }
  if (col_ind)   {delete [] col_ind; col_ind=NULL;    }
  if (row_ptr)   {delete [] row_ptr; row_ptr=NULL;    }
}

void CKSparse::multAx(double *Ax, double *x,int neqSub,int regularized) {
  // via function multAx, Ax = A*x, Ax - output vector, x - input vector
  // multAx(Ax, x, fem, 1);
  //int k = 0;
  for (int i = 0; i < neqSub; i++){
    Ax[i] = 0.0;
  }
  //
  for (int i = 0; i < neqSub; i++) {
    for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      Ax[i] += val[j] * x[col_ind[j]];
      if (j > row_ptr[i]) {
        Ax[col_ind[j]] +=	val[j] * x[col_ind[row_ptr[i]]];
      }
    }
  }
}

void CKSparse::fe_symbolic_form_matrix(CFem *fem, CDomain *domainG, CSolid45EqNumbGlobal* stif_glob_number) {
  if (!event1) LogEventRegister("fe_symbolic_form",&event1);
  LogEventBegin(event1);

  int MPIrank;
  MPI_Comm_rank(comm, &MPIrank);

  //XXXXXXXXXXXXXXXXXXXXXXXXX maping vectors g2l, l2g, start XXXXXXXXXXXXXXXXXXXXXXXXX
  const int n_elmt_subdom = domainG->n_elementsSub[fem->i_domOnClust];
  clock_t begin = clock();
  int neqSub = 3 * domainG->n_nodsSub[fem->i_domOnClust];
  domainG->neqSub[fem->i_domOnClust] = neqSub;

  int m = 24;
  int cnt = 0;
  //
  this->n_row = neqSub;
  //
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'KSparse', before 'g2l' and 'l2g' ... " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif
  const int ndofsOnSubdom = n_elmt_subdom*24;
  std::map<longInt,int> & g2l = fem->mesh.g2l;
  std::vector<longInt> & l2g = fem->mesh.l2g;
  int *ieq_glob= new int[24];
  int cntijv=0;
  l2g.resize(ndofsOnSubdom);
  //creation of DOFs indexes in global numbering
  for (int i = 0; i < domainG->n_elementsSub[fem->i_domOnClust]; i++) {
    CSolid45 &element_i = fem->mesh.element[i];
    for (int j = 0; j < 8; j++) {
      for (int k = 0; k < 3; k++) {
        stif_glob_number[cnt].ieq[j + k * 8] = 3 * ((element_i.inod_glob[j])) + k;
        ieq_glob[j + k * 8] = stif_glob_number[cnt].ieq[j + k * 8] + 0;
        l2g[cntijv] = ieq_glob[j + k * 8];
        cntijv++;
      } // loop k
    } // loop j
    cnt++;
  } // loop i
  delete [] ieq_glob;ieq_glob = NULL;

  //
  qsort(&(l2g[0]),ndofsOnSubdom, sizeof (longInt), CLinearAlgebra::compareLongInt);
  std::vector<longInt>::iterator itv;
  itv = std::unique (l2g.begin(), l2g.end());
  l2g.resize( std::distance(l2g.begin(),itv) );
  //
  for (unsigned int ii = 0; ii < l2g.size(); ii++){
    g2l.insert ( std::pair<longInt,int>(l2g[ii],ii) );
  }
  //XXXXXXXXXXXXXXXXXXXXXXXXX maping vectors g2l, l2g, start XXXXXXXXXXXXXXXXXXXXXXXXX
  // renumbering 'ieq' indexes (in global) to local 
  for (int k = 0; k < domainG->n_elementsSub[fem->i_domOnClust]; k++) {
    CSolid45 &element_i = fem->mesh.element[k];
    for (int j = 0; j<24; j++){
      stif_glob_number[k].ieq[j] = g2l[stif_glob_number[k].ieq[j]];
    } 
  }
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'KSparse', after 'g2l' and 'l2g' ... " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif

#ifdef flag_Fortran
  //
  int *ieq_glob_fortran = new int[24];
  fe2feti_init_symbolic_(fem->domainG->neqAll, n_elmt_subdom, 24 * n_elmt_subdom, 1);

  cnt=0;
  for (int i = 0; i < fem->domainG->n_elementsSub; i++) {
    for (int j = 0;j<24;j++){
      ieq_glob_fortran[j] = fem->mesh.l2g[stif_glob_number[i].ieq[j]]+1;
    }
    fe2feti_symbolic_map_element_(24, ieq_glob_fortran, 1);
  } 
  //
  delete [] ieq_glob_fortran; ieq_glob_fortran=NULL;
  //
  int *tmp_ind    = new int[fem->domainG->neqAll];
  double *tmp_val = new double[fem->domainG->neqAll];
  memset(tmp_ind, 0,fem->domainG->neqAll * sizeof(int));
  memset(tmp_val, 0,fem->domainG->neqAll * sizeof(double));
  if (fem->bound_cond->dirBCSub->n>0){
    for (int i = 0;i<fem->bound_cond->dirBCSub->n;i++){
      tmp_ind[fem->mesh.l2g[fem->bound_cond->dirBCSub->ind[i]]] = 1;
      tmp_val[fem->mesh.l2g[fem->bound_cond->dirBCSub->ind[i]]] = 
        fem->bound_cond->dirBCSub->val[i];
    }
  }
  int *tmp_ind2 = new int[fem->domainG->neqAll];
  double *tmp_val2 = new double[fem->domainG->neqAll];
  memset(tmp_ind2, 0,fem->domainG->neqAll * sizeof(int));
  memset(tmp_val2, 0,fem->domainG->neqAll * sizeof(double));

  MPI_Allreduce(tmp_ind,tmp_ind2,fem->domainG->neqAll,MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(tmp_val,tmp_val2,fem->domainG->neqAll,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  int n_DirBC_global = 0;
  for (int i = 0;i<fem->domainG->neqAll;i++){
    if (tmp_ind2[i]!=0) n_DirBC_global++;
  }
  fem->bound_cond->dirBC_global->n    = n_DirBC_global;
  fem->bound_cond->dirBC_global->ind  = new int [n_DirBC_global];
  fem->bound_cond->dirBC_global->val  = new double[n_DirBC_global];
  cnt=0;
  for (int i = 0;i<fem->domainG->neqAll;i++){
    if (tmp_ind2[i]!=0){
      fem->bound_cond->dirBC_global->ind[cnt]=i+1;
      fem->bound_cond->dirBC_global->val[cnt]=tmp_val2[i];    
      cnt++;
    }
  }
  delete [] tmp_ind2;  tmp_ind2=NULL;
  delete [] tmp_val2;  tmp_val2=NULL;
  delete [] tmp_ind;   tmp_ind=NULL;
  delete [] tmp_val;   tmp_val=NULL;

  fe2feti_symbolic_map_global_bc_(n_DirBC_global,fem->bound_cond->dirBC_global->ind);
  fe2feti_symbolic_finalize_();
  fe2feti_symbolic_factorize_();
  fe2feti_init_numeric_(1);
  return;
  //
#else
  
  //int cnt0 = 0;
  int nnz_K= 0;
  if (0) {
    if (MPIrank == 0) printf(" - - - - - - - - - - - - - - Ats. CSR ver.\n");
    // +++ Atsushi version +++
    vector< list < int > > _ind_cols_tmp(neqSub);// = new list<int>[neqSub];
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', before _ind_cols_tmp fields ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
    longInt *elm;
    for (int k = 0; k < domainG->n_elementsSub[fem->i_domOnClust]; k++) {
      elm = stif_glob_number[k].ieq;
      for (int i = 0; i < m; i++) {
        int ii = elm[i];// elm[i];
        for (int j = 0; j < m; j++) {
          int jj = elm[j];//elm[j]];
          if (jj >= ii) {
            if (_ind_cols_tmp[ii].empty()) {
              _ind_cols_tmp[ii].push_back(jj); 
            }	else {
              if (_ind_cols_tmp[ii].back() < jj) {
                _ind_cols_tmp[ii].push_back(jj);
              }	else {
                for (list<int>::iterator it = _ind_cols_tmp[ii].begin();
                    it != _ind_cols_tmp[ii].end(); ++it) {
                  if (*it == jj) {
                    break; // non-zero entry is already registered
                  }
                  if (*it > jj) {
                    _ind_cols_tmp[ii].insert(it, jj);
                    break;
                  }
                }
              } // if (_ind_cols_tmp[ii].back() < jj)
            } // if (_ind_cols_tmp[ii].empty())
          } // if (jj >= ii)
        } // loop : j
      } // loop : i
    } // loop : k
    //
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', after _ind_cols_tmp fields ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
    //
    for (int i = 0; i < neqSub; i++) {
      nnz_K += _ind_cols_tmp[i].size();
    }
    //
    //cout << "nnz = " << nnz_K << endl;
    col_ind = new int[nnz_K];
    val 		= new double[nnz_K];
    row_ptr = new int[neqSub + 1];

#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', col_ind, val and row_ptr created ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
    //
    nnz_K = 0;
    for (int i = 0; i < neqSub; i++) {
      row_ptr[i] = nnz_K;
      for (list<int>::iterator it = _ind_cols_tmp[i].begin();
          it != _ind_cols_tmp[i].end(); ++it) {
        col_ind[nnz_K] = *it;
        val[nnz_K] = 0.0;
        ++nnz_K;
      }
    } // loop i
    //
    //cout << "nnz = " << nnz_K << " /control count/" << endl;
    row_ptr[neqSub] = nnz_K;
    //
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //
    //if (!MPIrank) {
    //	printf("CSR preprocessing ...%3.1es\n", elapsed_secs);
    //}
    list <int> tmp_list(0);
    for (int i = 0; i < neqSub;i++){// = new list<int>[neqSub];
      _ind_cols_tmp[i].clear();
      _ind_cols_tmp[i] = tmp_list;
    }
    //
  }
  else
  { // +++ my version +++
    //if (MPIrank == 0) printf(" - - - - - - - - - - - - - - my CSR ver.\n");
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', before 'forCSRformat' ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
    // variable 'bool symMatrix' prepared for nonsymm. cases
    bool symMatrix=true;
    longInt *indK;
    vector < vector < int > > forCSRformat;
    forCSRformat.resize(domainG->neqSub[fem->i_domOnClust],vector<int>(0));
//
    int elemDOFs = 24;
    for (int i=0;i<domainG->n_elementsSub[fem->i_domOnClust];i++){
      indK = stif_glob_number[i].ieq;
      for (int j=0;j<elemDOFs;j++){
        for (int k=0;k<elemDOFs;k++){
          if (symMatrix && (indK[k]>=indK[j]) || !symMatrix){
            forCSRformat[indK[j]].push_back((int)indK[k]);
          }
        }
      }
    }
    nnz_K = 0;
    std::vector<int>::iterator itv_int;
    for (int i = 0;i<forCSRformat.size();i++){
      sort(forCSRformat[i].begin(),forCSRformat[i].end());
      itv_int = std::unique (forCSRformat[i].begin(), forCSRformat[i].end());
      forCSRformat[i].resize( std::distance(forCSRformat[i].begin(),itv_int) );
      nnz_K+=forCSRformat[i].size();
    }
//
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', after 'forCSRformat' before CSR fields ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
    col_ind    = new int[nnz_K];
    val 		   = new double[nnz_K];
    row_ptr    = new int[neqSub + 1];
    memset(val,0,sizeof(double)*nnz_K);
    int cnt_it = 0;
    row_ptr[0] = 0;
    for (int i=0;i<forCSRformat.size();i++){
      for (vector<int>::iterator it1 = forCSRformat[i].begin();it1!=forCSRformat[i].end();it1++){
        col_ind[cnt_it]	= *it1;
        cnt_it++; 
        row_ptr[i+1]=cnt_it;
      }
    }
#ifdef HFETI_SOLVER
    //if (MPIrank == 0) {
    //  cout << " in 'KSparse', after CSR fields ... " << endl; 
    //  GetProcessMemoryStat ( ); GetMemoryStat( );
    //}
#endif
  }
  //if (MPIrank == 0) printf("\t\t\t neqSub=%d, nnz_K=%d \n",neqSub,nnz_K);
#endif
#ifdef HFETI_SOLVER
  //if (MPIrank == 0) {
  //  cout << " in 'KSparse', leaving of function ... " << endl; 
  //  GetProcessMemoryStat ( ); GetMemoryStat( );
  //}
#endif


  LogEventEnd(event1);
}

void CKSparse::fe_numeric_form_matrix(CSolid45EqNumbGlobal * stif_loc_numb, CStiffnessLocal * stif_loc,CFem * fem) {
  if (!event2) LogEventRegister("fe_numeric_form",&event2);
  LogEventBegin(event2);

  // routine fills structures
  std::map<longInt,int> & g2l = fem->mesh.g2l;
  double K_ij;
  int i_ind, j_ind;
  //
  for (int i = 0; i < 24; i++) {
    for (int j = i; j < 24; j++) {
      j_ind = stif_loc_numb->ieq[j];
      K_ij = stif_loc->value_K[i * 24 + j];
      if (j_ind < stif_loc_numb->ieq[i]) {
        i_ind = stif_loc_numb->ieq[j];
        j_ind = stif_loc_numb->ieq[i];
      } else {
        i_ind = stif_loc_numb->ieq[i];
      }
      for (int k = row_ptr[i_ind]; k < row_ptr[i_ind + 1]; k++) {
        if ((col_ind[k]) == j_ind) {
          val[k] += K_ij;
          break;
        } // if
      } // loop: k
    } // loop: j
  } // loop: i

  LogEventEnd(event2);
} //fe_numeric_form_matrix
//
//

