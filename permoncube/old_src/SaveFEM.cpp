#include "SaveFEM.h"
// not necessary to change 
CSaveFEM::CSaveFEM(MPI_Comm comm, CFem* fem, CData* data, CDomain *domainG) {
  this->comm = comm;
  this->fem = fem;
  this->data = data;
  this->domainG= domainG;
  MPI_Comm_rank(comm,&this->MPIrank);
  MPI_Comm_size(comm,&this->MPIsize);
}
CSaveFEM::~CSaveFEM() {
}

void CSaveFEM::GatherDataFemToMaster(){
//
//  int MPIrank, MPIsize;
//  MPI_Comm comm=MPI_COMM_WORLD;
//  MPI_Comm_rank(comm, &MPIrank);
//  MPI_Comm_size(comm, &MPIsize);
//  MPI_Status mpi_stat; 
//  MPI_Request mpi_req; 
//
//
//  bool b1 =( !domainG->flag_store_VTK ||
//      (domainG->neqAll > domainG->max_nDOFs_u_is_stored));
//  if (b1) {
//  	if (!MPIrank){
//		  printf("VTK is not creating ... \n");
//    }
//	  return;
//  }
//
//
//  int size_vector_DirBC;
//  int nDOFs,neqSub,n_elementsSub,n_nodsSub;
////
//  double  *coordinatesOnMaster;
//  int     *elements_all;
//  int     *ptr_el;
//  int     *ind_subdom_ddm;
//  
//  coordinatesOnMaster = NULL;
//  elements_all        = NULL;
//  ptr_el              = NULL;
//  ind_subdom_ddm      = NULL;       
//
//  if (!MPIrank) printf("VTK is creating ... \n");
//  // number of subdom DOFs sent to rank=0 ||| START ||| - - - - - - - -
//  if (!MPIrank) printf("VTK is creating: start ... \n");
//  int n_buff          = 4;
//  int *mpi_buff= new int[n_buff];
//  if (MPIrank>0){                                          
//    mpi_buff[0]  = domainG->neqSub[fem->i_domOnClust];
//    mpi_buff[1]  = domainG->n_elementsSub[fem->i_domOnClust];
//    mpi_buff[2]  = domainG->n_nodsSub[fem->i_domOnClust];
//    mpi_buff[3]  = domainG->n_elementsSub[fem->i_domOnClust]*8;
//    MPI_Send(mpi_buff, n_buff, MPI_INT, 0 , 1, MPI_COMM_WORLD);
//  }                                                        
//  else {                                                   
//    int n_subdomains = domainG->n_subdomains;
//    //
//    domainG->vec_neqSub           = new int[n_subdomains];
//    domainG->vec_neqSub[0]        = domainG->neqSub; 
//    
//    domainG->vec_n_elementsSub    = new int[n_subdomains];
//    domainG->vec_n_elementsSub[0] = domainG->n_elementsSub;
//    
//    domainG->vec_n_nodsSub        = new int[n_subdomains];
//    domainG->vec_n_nodsSub[0]     = domainG->n_nodsSub;
//    
//    domainG->vec_numelSub         = new int[n_subdomains];
//    domainG->vec_numelSub[0]      = domainG->n_elementsSub*8;
//    //
//    for (int i = 1;i<n_subdomains;i++){
//      MPI_Recv(mpi_buff, n_buff,MPI_INT, i, 1, MPI_COMM_WORLD, &mpi_stat);
//      domainG->vec_neqSub[i]        = mpi_buff[0];
//      domainG->vec_n_elementsSub[i] = mpi_buff[1];
//      domainG->vec_n_nodsSub[i]     = mpi_buff[2];
//      domainG->vec_numelSub[i]      = mpi_buff[3];
//    }
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    if (false) {
//						for (int i = 0;i<MPIsize;i++){
//							printf("neqSub = %d ",       domainG->vec_neqSub[i]);
//							printf("n_elementsSub = %d ",domainG->vec_n_elementsSub[i]);
//							printf("n_nodsSub = %d",     domainG->vec_n_nodsSub[i]);
//							printf("vec_numelSub = %d\n",domainG->vec_numelSub[i]);
//						}
//    }
//  }
//  delete [] mpi_buff; mpi_buff = NULL;
//  // number of subdom DOFs sent to rank=0 ||| END   ||| - - - - - - - -
//  if (!MPIrank) printf("VTK is creating: sending of maping vector to master \n");
//  // 
//  // gathering of l2g map vecs to  rank=0 ||| START ||| - - - - - - - -
//  vector < vector < int > > l2g_vector;  
//  if (MPIrank>0){                                          
//    neqSub = domainG->neqSub;
//    MPI_Send(&fem->mesh.l2g[0], neqSub, MPI_INT, 0 , 1, MPI_COMM_WORLD);
//  }                                                        
//  else {                                                   
//    l2g_vector.resize( domainG->n_subdomains, vector<int>( 0 , 0 ) );
//   	l2g_vector[0].insert(l2g_vector[0].begin(),
//        &fem->mesh.l2g[0] , &fem->mesh.l2g[0] + domainG->neqSub);
//    for (int i = 1;i < domainG->n_subdomains;i++){
//      neqSub= domainG->vec_neqSub[i];
//      int *mpi_buff = new int [neqSub];
//      MPI_Recv(mpi_buff,neqSub, MPI_INT,i, 1, MPI_COMM_WORLD,&mpi_stat);
//    	l2g_vector[i].insert(l2g_vector[i].begin(), mpi_buff, mpi_buff + neqSub);
//      delete [] mpi_buff; mpi_buff = NULL;
//    }
//  }
//  // gathering of l2g map vecs to  rank=0 ||| END   ||| - - - - - - - -
//  if (!MPIrank) printf("VTK is creating: primal. sol. contribution sent to master \n");
//  //
//  // gathering of u_restrict to    rank=0 ||| START ||| - - - - - - - -
//  if (MPIrank>0){
//    neqSub= domainG->neqSub;
//    double *data_u_sub = new double[neqSub];
//    for (int i=0;i<domainG->neqSub;i++){
//      data_u_sub[i] = data->u[i]/data->B->multiplicityPD[i];
//    } 
//    MPI_Send(&data_u_sub[0], neqSub, MPI_DOUBLE, 0 , 1, MPI_COMM_WORLD);
//    delete [] data_u_sub;  data_u_sub=NULL;
//  }
//  else {
////    if (data->u_restrict==NULL){
//      memset(data->u_restrict,0,domainG->neqAll*sizeof(double));
// //   }
//    for (int j=0;j<domainG->neqSub;j++){
//      data->u_restrict[l2g_vector[0][j]] += 
//          data->u[j]/data->B->multiplicityPD[j];
//    }
//    double *mpi_buff;
//    for (int i = 1;i < domainG->n_subdomains;i++){
//      neqSub= domainG->vec_neqSub[i];
//      mpi_buff = new double[neqSub];
//      MPI_Recv(mpi_buff,neqSub, MPI_DOUBLE,i, 1, MPI_COMM_WORLD,&mpi_stat);
//      for (int j = 0; j<neqSub;j++){
//        data->u_restrict[l2g_vector[i][j]]+=mpi_buff[j];
//      }
//      delete [] mpi_buff;mpi_buff=NULL;
//    }
//  }
//  // gathering of u_restrict to    rank=0 ||| END   ||| - - - - - - - -
//  if (!MPIrank) printf("VTK is creating: sending of elements to master \n");
//  //
//  // gathering of mesh to rank=0          ||| START ||| - - - - - - - -
//  if (MPIrank==0){
//    int numel_glob_el_table=0;
//    for (int i = 0;i<domainG->n_subdomains;i++){
//      numel_glob_el_table += domainG->vec_numelSub[i];
//    }
//    elements_all    = new int [numel_glob_el_table];
//    ptr_el          = new int [domainG->n_elementsAll+1];
//    ind_subdom_ddm  = new int [domainG->n_elementsAll];
//    ptr_el[0] = 0;
//    int i_ptr_el=0;
//    int n_nodesOnEl,cnt_ij=0;
//    for (int i = 0;i<domainG->n_elementsSub;i++){
//      n_nodesOnEl = 8;
//      for (int j = 0;j<n_nodesOnEl;j++){
//        elements_all[cnt_ij]=fem->mesh.element[i].inod_glob[j]; 
//        cnt_ij++;
//      }
//      ptr_el[i+1] = cnt_ij; 
//      ind_subdom_ddm[i_ptr_el] = 0;
//      i_ptr_el++;
//    }
//    int *mpi_buff;
//    int *mpi_buff2;
//    int cnt_ptr = 0,tmp_int=cnt_ij;
//    for (int i = 1;i < domainG->n_subdomains;i++){
//      int numelSub = domainG->vec_numelSub[i];
//      mpi_buff = new int [numelSub];
//      MPI_Recv(mpi_buff,numelSub, MPI_INT,i, 1, MPI_COMM_WORLD,&mpi_stat);
//      for (int j = 0; j<numelSub;j++){
//        elements_all[cnt_ij] = mpi_buff[j];
//        cnt_ij++;
//      }
//      int n_elementsSub = domainG->vec_n_elementsSub[i];
//      mpi_buff2 = new int [numelSub];
//      int cnt_j=0;
//      MPI_Recv(mpi_buff2,n_elementsSub+1, MPI_INT,i, 2, MPI_COMM_WORLD,&mpi_stat);
//      cnt_ptr+=domainG->vec_n_elementsSub[i];
//      for (int j = 0; j<n_elementsSub;j++){
//        ptr_el[i_ptr_el] = mpi_buff2[j] + tmp_int;
//        ind_subdom_ddm[i_ptr_el] = i;
//        i_ptr_el++;
//      }
//      tmp_int += mpi_buff2[n_elementsSub];
//      delete [] mpi_buff;    mpi_buff  = NULL;
//      delete [] mpi_buff2;   mpi_buff2 = NULL;
//    }
//    ptr_el[i_ptr_el] = tmp_int;
//  }
//  else {
//    int numelSub        = domainG->n_elementsSub*8;
//    int *elements_sub   = new int [numelSub];
//    int *ptr_el_sub     = new int [domainG->n_elementsSub+1];
//    ptr_el_sub[0] = 0;
//    int n_nodesOnEl;
//    int cnt_ij = 0;;
//    for (int i = 0;i<domainG->n_elementsSub;i++){
//      n_nodesOnEl = 8;
//      for (int j = 0;j<n_nodesOnEl;j++){
//        elements_sub[cnt_ij]=fem->mesh.element[i].inod_glob[j]; 
//        cnt_ij++;
//      }
//      ptr_el_sub[i+1] = cnt_ij; 
//    }
//    MPI_Send(elements_sub, numelSub, MPI_INT, 0 , 1, MPI_COMM_WORLD);
//    MPI_Send(ptr_el_sub, domainG->n_elementsSub+1, 
//                                    MPI_INT, 0 , 2, MPI_COMM_WORLD);
//    delete [] elements_sub;  elements_sub=NULL;
//    delete [] ptr_el_sub;    ptr_el_sub  = NULL;
//  }
//  // gathering of mesh to rank=0          ||| END   ||| - - - - - - - -
//  if (!MPIrank) printf("VTK is creating: sending of coordinates to master \n");
//  //
//  // gathering of coordinates to   rank=0 ||| START ||| - - - - - - - -
//  if (MPIrank>0){
//    neqSub= domainG->neqSub;
//    double *data_coord= new double[neqSub];
//    for (int i=0;i<domainG->n_nodsSub;i++){
//      data_coord[3*i  ] = fem->mesh.coordinateSub[i].x;
//      data_coord[3*i+1] = fem->mesh.coordinateSub[i].y;
//      data_coord[3*i+2] = fem->mesh.coordinateSub[i].z;
//    } 
//    MPI_Send(&data_coord[0], neqSub, MPI_DOUBLE, 0 , 1, MPI_COMM_WORLD);
//    delete [] data_coord;data_coord = NULL;
//  }
//  else {
//    coordinatesOnMaster = new double [domainG->neqAll]; 
//    for (int j=0;j<domainG->n_nodsSub;j++){
//      coordinatesOnMaster[l2g_vector[0][3*j  ]] = fem->mesh.coordinateSub[j].x;
//      coordinatesOnMaster[l2g_vector[0][3*j+1]] = fem->mesh.coordinateSub[j].y;
//      coordinatesOnMaster[l2g_vector[0][3*j+2]] = fem->mesh.coordinateSub[j].z;
//    }
//    double *mpi_buff;
//    for (int i = 1;i < domainG->n_subdomains;i++){
//      neqSub= domainG->vec_neqSub[i];
//      mpi_buff = new double[neqSub];
//      MPI_Recv(mpi_buff,neqSub, MPI_DOUBLE,i, 1, MPI_COMM_WORLD,&mpi_stat);
//      for (int j = 0; j<neqSub;j++){
//        coordinatesOnMaster[l2g_vector[i][j]]=mpi_buff[j];
//      }
//      delete [] mpi_buff; mpi_buff = NULL;
//    }
//  }
//  if (!MPIrank) printf("VTK is creating: all data are on the master \n");
//  // gathering of coordinates to   rank=0 ||| END   ||| - - - - - - - -
//  //
////########$########$###########################$########$########
////########$########$#### #### #      # ### ####$########$########
////########$########$#### ### #### #### ## #####$########$########
////########$########$#### ## #### ####   #######$########$########
////########$########$#### # ##### #### ## ######$########$########
////########$########$####  ##### #### #### #####$########$########
////########$########$###########################$########$########
//  if (MPIrank==0){
//#ifndef WIN32
//    mode_t mode = 0777;
//    mkdir("data", mode);
//#endif
//    if (!MPIrank) printf("VTK is creating: printing into file\n");
//    FILE *fVTK = fopen("data/box_new.vtk", "w");
//    fprintf(fVTK, "# vtk DataFile Version 3.0\n");
//    fprintf(fVTK, "vtk output\n");
//    fprintf(fVTK, "ASCII\n\n");
//    fprintf(fVTK, "DATASET UNSTRUCTURED_GRID\n");
//    fprintf(fVTK, "POINTS %d float\n",domainG->n_nodsAll);
//    for (int i = 0; i <domainG->n_nodsAll; i++) {
//      fprintf(fVTK, "%f %f %f\n", coordinatesOnMaster[3*i+0],
//                                  coordinatesOnMaster[3*i+1],
//                                  coordinatesOnMaster[3*i+2]);
//    }
//    int k1 = domainG->n_elementsAll * 9;
//    fprintf(fVTK, "CELLS %d %d\n",domainG->n_elementsAll, k1);
//    for (int i = 0; i < domainG->n_elementsAll; i++) {
//      fprintf(fVTK, "%d ", ptr_el[i+1]-ptr_el[i]);//fprintf(fVTK, "%d ", 8);
//      for (int j = ptr_el[i]; j < ptr_el[i+1]; j++) {
//        fprintf(fVTK, "% d", elements_all[j]);
//      }
//      fprintf(fVTK, "\n");
//    }
//    fprintf(fVTK, "CELL_TYPES %d\n", domainG->n_elementsAll);
//    for (int i = 0; i < domainG->n_elementsAll; i++) {
//      fprintf(fVTK, "%d\n", 12);
//    }
//    fprintf(fVTK, "POINT_DATA %d\n",domainG->n_nodsAll);
//    fprintf(fVTK, "SCALARS displacements float 3\n");
//    fprintf(fVTK, "LOOKUP_TABLE my_table\n");
//    for (int i = 0; i < domainG->n_nodsAll; i++) {
//      fprintf(fVTK, "%f %f %f\n", data->u_restrict[3 * i + 0], 
//                                  data->u_restrict[3 * i + 1],  
//                                  data->u_restrict[3 * i + 2]);
//    }
//    fprintf(fVTK, "\nCELL_DATA %d\n",domainG->n_elementsAll);
//    fprintf(fVTK, "SCALARS decomposition int 1\n");
//    fprintf(fVTK, "LOOKUP_TABLE decomposition\n");
//    for (int i = 0; i < domainG->n_elementsAll; i++) {
//      fprintf(fVTK, "%d\n", ind_subdom_ddm[i]+1);
//    }
//    fclose(fVTK);
//  }
//  if (!elements_all)  { delete [] elements_all;    elements_all=NULL;  }   
//  if (!ptr_el)        { delete [] ptr_el;          ptr_el=NULL;        }
//  if (!ind_subdom_ddm){ delete [] ind_subdom_ddm;  ind_subdom_ddm=NULL;}
//  if (!MPIrank) printf("VTK is creating: finish\n");
}

void CSaveFEM::save_data() {
  if (!MPIrank) {
    printf(" +++ data saving +++\n");
  }
  saveElements();
  saveCoordinates();
  saveLocalStifnessMat(data->KSparse);
}

void CSaveFEM::saveElements() {
//  // saving text "elements0.dat" --------------------------------------
//  clock_t begin = clock();
//  ofstream fel("data/elements0.dat");
//  for (int i = 0; i < domainG->n_elementsSub; i++) {
//    for (int j = 0; j < 8; j++) {
//      fel << fem->mesh.element[i].inod_glob[j] << "\t";
//    }
//    fel << fem->mesh.element[i].indSubdomain << "\t";
//    fel << endl;
//  }
//  fel.close(); //---------------------------------------------------------
//  clock_t end = clock();
//  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//  if (!MPIrank) {
//    printf("elem0  ... %3.1e s\n", elapsed_secs);
//  }
}

void CSaveFEM::saveCoordinates() {
//  // saving text "coordinates0.dat" ------------------------------------
//  clock_t begin = clock();
//  ofstream fcoord("data/coordinates0.dat");
//  for (int m = 0; m < fem->domainG->n_nodsAll; m++)
//    fcoord << fem->mesh.coordinate[m].x << "\t"
//      << fem->mesh.coordinate[m].y << "\t"
//      << fem->mesh.coordinate[m].z << endl;
//  fcoord.close(); //----------------------------------------------------------
//  clock_t end = clock();
//  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//  if (!MPIrank) {
//    printf("coord0 ... %3.1e s\n", elapsed_secs);
//  }
}




void CSaveFEM::saveLocalStifnessMat(CKSparse * Ksparse) {
//  //++++++++storing of local stiffness matrices++++++++++++++++
//  clock_t begin = clock();
//  char filenameKe[128];
//  sprintf(filenameKe, "data/K_%d.dat", 0);
//  ofstream f_loc_stif_mat(filenameKe);
//  for (int i = 0; i < domainG->neqSub; i++) {
//    for (int j = Ksparse->row_ptr[i]; j < Ksparse->row_ptr[i + 1]; j++) {
//      f_loc_stif_mat << i << "\t" << Ksparse->col_ind[j] << "\t"
//        << setprecision(12) << Ksparse->val[j] << endl;
//    }
//  }
//  f_loc_stif_mat.close(); //                                000
//  clock_t end = clock();
//  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//  if (!MPIrank) {
//    printf("K      ... %3.1e s\n", elapsed_secs);
//  }
//  //++++++++storing of local stiffness matrices++++++++++++++++
}

void CSaveFEM::saveVTK() {
//  bool b1 = MPIrank || 
//       !domainG->flag_store_VTK ||
//      (domainG->neqAll > domainG->max_nDOFs_u_is_stored);
//  if (b1) return;
//
//  // saving VTK --------------------------------------
//  //clock_t begin = clock();
//#ifndef WIN32
//  mode_t mode = 0777;
//  mkdir("data", mode);
//#endif
//
//  FILE *fVTK = fopen("data/box.vtk", "w");
//  fprintf(fVTK, "# vtk DataFile Version 3.0\n");
//  fprintf(fVTK, "vtk output\n");
//  fprintf(fVTK, "ASCII\n\n");
//  fprintf(fVTK, "DATASET UNSTRUCTURED_GRID\n");
//  fprintf(fVTK, "POINTS %d float\n", domainG->n_nodsAll);
//  for (int i = 0; i < domainG->n_nodsSub; i++) {
//    fprintf(fVTK, "%f %f %f\n", fem->mesh.coordinateSub[i].x,
//        fem->mesh.coordinateSub[i].y, fem->mesh.coordinateSub[i].z);
//  }
//  int k1 = domainG->n_elementsAll * 9;
//  fprintf(fVTK, "CELLS %d %d\n", domainG->n_elementsAll, k1);
//  for (int i = 0; i < domainG->n_elementsAll; i++) {
//    fprintf(fVTK, "%d ", 8);
//    for (int j = 0; j < 8; j++) {
//      fprintf(fVTK, "% d", fem->mesh.element[i].inod_glob[j]);
//    }
//    fprintf(fVTK, "\n");
//  }
//  fprintf(fVTK, "CELL_TYPES %d\n", domainG->n_elementsAll);
//  for (int i = 0; i < domainG->n_elementsAll; i++) {
//    fprintf(fVTK, "%d\n", 12);
//  }
//  fprintf(fVTK, "POINT_DATA %d\n", domainG->n_nodsAll);
//  fprintf(fVTK, "SCALARS displacements float 3\n");
//  fprintf(fVTK, "LOOKUP_TABLE my_table\n");
////#ifdef flag_Fortran
//  for (int i = 0; i < domainG->n_nodsAll; i++) {
//    fprintf(fVTK, "%f %f %f\n", data->u_restrict[3 * i + 0], data->u_restrict[3 * i + 1], data->u_restrict[3 * i + 2]);
//  }
////#endif 
//  fprintf(fVTK, "\nCELL_DATA %d\n", domainG->n_elementsAll);
//  fprintf(fVTK, "SCALARS decomposition int 1\n");
//  fprintf(fVTK, "LOOKUP_TABLE decomposition\n");
//  for (int i = 0; i < domainG->n_elementsAll; i++) {
//    fprintf(fVTK, "%d\n", fem->mesh.element[i].indSubdomain + 1);
//  }
//  fclose(fVTK);

}

void CSaveFEM::saveStiffnessMatrix(){
  //
//  char filenameKe[128]; 
//  sprintf(filenameKe, "data/K_%d.dat", MPIrank);
//  ofstream f_loc_stif_mat(filenameKe);
//  int neqSub = domainG->neqSub;
//  for (int i = 0; i < neqSub; i++) {
//    for (int j = data->KSparse->row_ptr[i]; j < data->KSparse->row_ptr[i + 1]; j++) {
//      f_loc_stif_mat << i
//        << "\t" << data->KSparse->col_ind[j] 
//        << "\t" << setprecision(12) << data->KSparse->val[j] << endl;
//    }
//  }
//  f_loc_stif_mat.close(); //
}

void CSaveFEM::saveVector_int(int *x,int n,char * name){
  //    ofstream frhs_f("data/f_RHS.dat");
  ofstream frhs_f(name);

  for (int i = 0; i < n; i++) {

    frhs_f << x[i] << endl;
  }
  frhs_f.close(); //--
}

void CSaveFEM::saveVector_double(double *x,int n,char * name){
  //    ofstream frhs_f("data/f_RHS.dat");
  ofstream frhs_f(name);

  for (int i = 0; i < n; i++) {

    frhs_f << x[i] << endl;
  }
  frhs_f.close(); //--
}
