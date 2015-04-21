#include "Cluster_g.h"

CClust_g::CClust_g(MPI_Comm comm,int argc, char** argv) {
  this->comm  = comm;
	domainG      = new CDomain(comm);
  B1           = new CBSparseC(comm);
  domainG->setFromOptions(&argc,&argv);
  domainG->print();
  fem.reserve(domainG->n_subdomOnClust);
  data.reserve(domainG->n_subdomOnClust);
//--
  for (int i =0;i<domainG->n_subdomOnClust;i++){
    fem.push_back(new CFem(comm,i));
    data.push_back(new CData(comm,i));
  }
  argc_l = argc;
  argv_l = argv; 
}

CClust_g::~CClust_g() {
	if (domainG) 		{	delete domainG;        domainG      = NULL;}
  delete [] neighbClst;
  delete [] indExterDOFsOnClust;
  for (int i =0;i<fem.size();i++){
    delete fem[i];
    delete data[i];
  }
}


void CClust_g::collectingExternalDOFs(){
// function is designed to collect external of each subdomains beloning to particular cluster


  int cnt=0;
  for (int i = 0;i<this->domainG->n_subdomOnClust;i++){
    //printf("n_exterDOFs              %d\n",this->domainG->n_exterDOFs[i]);
    cnt+=this->domainG->n_exterDOFs[i];
    //this->fem[i]->mesh.externalDOFs.reserve(this->domainG->n_exterDOFs[i]);
    this->fem[i]->mesh.externalDOFs.resize(this->domainG->n_exterDOFs[i]);
  }
  this->domainG->n_exterDOFsClust = cnt;

  this->indExterDOFsOnClust = new int [this->domainG->n_exterDOFsClust];

  cnt = 0;
  for (int i = 0;i<this->domainG->n_subdomOnClust;i++){
    for (int j = 0;j<this->domainG->n_exterDOFs[i];j++){
      this->indExterDOFsOnClust[j+cnt] = this->data[i]->indExterDOFs[j];
      this->fem[i]->mesh.externalDOFs[j] = (this->data[i]->indExterDOFs[j]);
    }
  cnt+=this->domainG->n_exterDOFs[i];
  }
// A2,start: Finding neigbours clusters of each cluster
  // + + + only for the cube - start
  int Cx = this->domainG->Cx;
  int Cy = this->domainG->Cy;
  int Cz = this->domainG->Cz;
//  int NxClst = this->domainG->NxClst;
//  int NyClst = this->domainG->NyClst;
//  int NzClst = this->domainG->NzClst;
  int Iclst,Jclst,Kclst;
  int globoalIndClst=this->domainG->MPIrank;


  Kclst = ceil(double (globoalIndClst+1)/(Cx*Cy));
  int inXYplane = (globoalIndClst+1)-(Kclst-1)*Cx*Cy;
  Jclst = ceil(double( inXYplane)/Cx);
  Iclst = inXYplane - Cx*(Jclst-1);
  // + + + only for the cube - end

  int *neighbClst = new int[26];
  cnt = 0;
  int tmpRankGuessed;
//
  for (int k = Kclst-1;k<=Kclst+1;k++){
		for (int j = Jclst-1;j<=Jclst+1;j++){
		  for (int i = Iclst-1;i<=Iclst+1;i++){
				tmpRankGuessed = Cx*Cy*(k-1)+Cx*(j-1)+i-1;
        if ((i>0 && j>0 && k>0) && 
						(i<Cx+1 && j<Cy+1 && k<Cz+1) &&
						(tmpRankGuessed != globoalIndClst)) {
			     neighbClst[cnt] =	tmpRankGuessed;
           cnt++;
				}
		  }
		}
  }

  this->domainG->n_neighbClst = cnt;
  this->neighbClst = neighbClst;


}

void CClust_g::getGlobalIndexOfSubdomMeshGen(){
//  
  int Iclst,Jclst,Kclst;
  int globoalIndClst=this->domainG->MPIrank;
  int Cx = this->domainG->Cx;
  int Cy = this->domainG->Cy;
  int Cz = this->domainG->Cz;
//
  Kclst = ceil(double (globoalIndClst+1)/(Cx*Cy));
  int inXYplane = (globoalIndClst+1)-(Kclst-1)*Cx*Cy;
  Jclst = ceil(double( inXYplane)/Cx);
  Iclst = inXYplane - Cx*(Jclst-1);
//
  int NxClst = this->domainG->NxClst;
  int NyClst = this->domainG->NyClst;
  int NzClst = this->domainG->NzClst;
//
  int Iclst0 = Iclst-1;
  int Jclst0 = Jclst-1;
  int Kclst0 = Kclst-1;
  int ddd;
  int cnt = 0;
  for (int k = Kclst0*NzClst;k<(Kclst0+1)*NzClst;k++){
    for (int j = Jclst0*NyClst;j<(Jclst0+1)*NyClst;j++){
      for (int i = Iclst0*NxClst;i<(Iclst0+1)*NxClst;i++){
        ddd = NxClst*Cx*NyClst*Cy*k + NxClst*Cx*j + i;
        this->domainG->vec_globalSubNumberingMeshGen[cnt]=ddd;
        cnt++;
		  }
		}
  }

}


void CClust_g::GatherDataFemToMaster(){
//
  int MPIrank, MPIsize,cnt;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &MPIrank);
  MPI_Comm_size(comm, &MPIsize);
  MPI_Status mpi_stat; 
  MPI_Request mpi_req; 
//
  bool b1 =( !this->domainG->flag_store_VTK ||
      (this->domainG->neqAll > this->domainG->max_nDOFs_u_is_stored));
  if (b1) {
  	if (!MPIrank){
		  printf("VTK is not creating ... \n"); 
    }
	  return;
  }
  if (!MPIrank){
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  }
//
//  int size_vector_DirBC;
  int nDOFs,neqSub,n_elementsSub,n_nodsSub,neqClst_i;
//
  double  *coordinatesOnMaster;
  int     *elements_all;
  int     *ptr_el;
  int     *sbd_ptr_el;
  int     *ind_subdom_ddm;
  int     FOUR=4;
// 
  coordinatesOnMaster = NULL;
  elements_all        = NULL;
  ptr_el              = NULL;
  ind_subdom_ddm      = NULL;       
//
  if (!MPIrank) printf("VTK is creating ... \n");
// number of subdom DOFs sent to rank=0 ||| START ||| - - - - - - - -
  if (!MPIrank) printf("VTK is creating: start ... \n");
  int *mpi_buff= new int[FOUR];
  if (MPIrank>0){                                          
    for (int i=0;i<FOUR;i++){
      mpi_buff[i]=0;
    }
//
    for (int i = 0;i<this->domainG->n_subdomOnClust;i++){
      mpi_buff[0] += this->domainG->neqSub[i];
      mpi_buff[1] += this->domainG->n_elementsSub[i];
      mpi_buff[2] += this->domainG->n_nodsSub[i];
      mpi_buff[3] += this->domainG->n_elementsSub[i]*8;
    }
    MPI_Send(mpi_buff, FOUR, MPI_INT, 0 , 1, MPI_COMM_WORLD);
  }                                                        
  else {                                                   
    int n_clusters = this->domainG->n_clusters;
//
    this->domainG->vec_neqClst          = new int[n_clusters];
    this->domainG->vec_n_elementsClst   = new int[n_clusters];
    this->domainG->vec_n_nodsClst       = new int[n_clusters];
    this->domainG->vec_numelClst        = new int[n_clusters];
//
    this->domainG->vec_neqClst[0]       = 0; 
    this->domainG->vec_n_elementsClst[0]= 0;
    this->domainG->vec_n_nodsClst[0]    = 0;
    this->domainG->vec_numelClst[0]     = 0;
//
    for (int i=0;i<this->domainG->n_subdomOnClust;i++){
      this->domainG->vec_neqClst[0]        += this->domainG->neqSub[i];
      this->domainG->vec_n_elementsClst[0] += this->domainG->n_elementsSub[i];
      this->domainG->vec_n_nodsClst[0]     += this->domainG->n_nodsSub[i];
      this->domainG->vec_numelClst[0]      += this->domainG->n_elementsSub[i]*8;
    }
    //
    printf("n_clusters = %d\n",n_clusters);
    for (int i = 1;i<n_clusters;i++){
      MPI_Recv(mpi_buff, FOUR,MPI_INT, i, 1, MPI_COMM_WORLD, &mpi_stat);
      this->domainG->vec_neqClst[i]         = mpi_buff[0];
      this->domainG->vec_n_elementsClst[i]  = mpi_buff[1];
      this->domainG->vec_n_nodsClst[i]      = mpi_buff[2];
      this->domainG->vec_numelClst[i]       = mpi_buff[3];
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (true) {
			//for (int i = 0;i<MPIsize;i++){
			//	printf("MPIrank= %d ",        i); 
			//	printf("neqClst = %d ",       this->domainG->vec_neqClst[i]);
			//	printf("n_elementsClst = %d ",this->domainG->vec_n_elementsClst[i]);
			//	printf("n_nodsClst = %d ",    this->domainG->vec_n_nodsClst[i]);
			//	printf("vec_numelClst = %d\n",this->domainG->vec_numelClst[i]);
			//}
    }
  }
  delete [] mpi_buff; mpi_buff = NULL;
// number of subdom DOFs sent to rank=0 ||| END   ||| - - - - - - - -
// checked until this line  
// ----------------------------------------------------------------------------------------
  
  
  
  
  if (!MPIrank) printf("VTK is creating: sending of maping vector to master \n");
  // 
  // gathering of l2g map vecs to  rank=0 ||| START ||| - - - - - - - -
  vector < vector < int > > l2g_vector;  
  int neqClst=0;
  for (int i=0;i<this->domainG->n_subdomOnClust;i++){
    neqClst += this->domainG->neqSub[i];
  }
  if (MPIrank>0){                                          
    vector < int >  l2g_clust;
    l2g_clust.reserve(neqClst);
    cnt=0;
    //printf(" aaa: neqClst=%d\n",neqClst);
    for (int i=0;i<this->domainG->n_subdomOnClust;i++){
      //printf("cnt = %d, size = %d\n",cnt,this->fem[i]->mesh.l2g.size());



    	l2g_clust.insert( l2g_clust.end(),
                        this->fem[i]->mesh.l2g.begin(),
                        this->fem[i]->mesh.l2g.end()    );
    }
    MPI_Send(&l2g_clust[0], neqClst, MPI_INT, 0 , 1, MPI_COMM_WORLD);
//    printf(" bbb: neqClst=%d\n",l2g_clust.size());
//    for (int i=0;i<l2g_clust.size();i++){
//      priniiitf("l2g_clust = %d\n",l2g_clust[i]);
//    }
  }                                                        
  else {                                                   
    l2g_vector.resize( this->domainG->n_clusters, vector<int>( 0 , 0 ) );
    cnt=0;
    // probably should be used resize for l2g_vector[0].resize()
  for (int i=0;i<this->domainG->n_subdomOnClust;i++){
    l2g_vector[0].insert( l2g_vector[0].end(),
                          this->fem[i]->mesh.l2g.begin(),  
                          this->fem[i]->mesh.l2g.end());
  }
    for (int i = 1;i < this->domainG->n_clusters;i++){
      neqClst_i=this->domainG->vec_neqClst[i];
      int *mpi_buff = new int [neqClst_i];
      MPI_Recv(mpi_buff,neqClst_i, MPI_INT,i, 1, MPI_COMM_WORLD,&mpi_stat);
//      l2g_vector[i].reserve(neqClst_i);
    	l2g_vector[i].insert(l2g_vector[i].end(), mpi_buff, mpi_buff + neqClst_i);
      delete [] mpi_buff; mpi_buff = NULL;
    }
  }


//  for (int i=0;i<this->domainG->n_subdomOnClust;i++){
//    for (int k=0;k<this->fem[i]->mesh.l2g.size();k++){
//      printf("\nk=%d, l2g = %d",k,this->fem[i]->mesh.l2g[k]);
//    }
//    printf("\n");
//  }


//  if (MPIrank==0){
//    for (int k = 0;k < this->domainG->n_clusters;k++){
//      printf("size(l2g_vector)=%d\n",l2g_vector[k].size());
//      for (int i=0;i<l2g_vector[k].size();i++){
//        printf("n_clust=%d, l2g_vector = %d\n",k,l2g_vector[k][i]);
//      }
//    }
//  }
  // gathering of l2g map vecs to  rank=0 ||| END   ||| - - - - - - - -
  // ---------------------------------------------------------------------
  // CODE is OK until this line.
  if (!MPIrank) printf("VTK is creating: primal. sol. contribution sent to master \n");
  //
  // gathering of u_restrict to    rank=0 ||| START ||| - - - - - - - -
  if (MPIrank>0){
    double *data_u_clust = new double[neqClst];
    cnt=0;
    for (int i=0;i<this->domainG->n_subdomOnClust;i++){
      for (int j=0;j<this->domainG->neqSub[i];j++){
        data_u_clust[j+cnt] = this->data[i]->u[j];///data->B->multiplicityPD[i];
      } 
     cnt+=this->domainG->neqSub[i];
    }
    MPI_Send(&data_u_clust[0], neqClst, MPI_DOUBLE, 0 , 1, MPI_COMM_WORLD);
    delete [] data_u_clust;  data_u_clust=NULL;
  }
  else {
    this->u_restrict = new double[this->domainG->neqAll];
    memset(this->u_restrict,0,this->domainG->neqAll*sizeof(double));
    cnt=0;
    for (int i=0;i<this->domainG->n_subdomOnClust;i++){
      for (int j=0;j<this->domainG->neqSub[i];j++){
//TODO  multiplicityPD does not work properly, that's why
//      next line is commented out (see the difference '=' instead of '+=')
//        this->u_restrict[l2g_vector[0][j+cnt]] += 
        this->u_restrict[l2g_vector[0][j+cnt]] = 
            this->data[i]->u[j];///data->B->multiplicityPD[j];
////        printf("dom=%d, u=%f\n",i,this->data[i]->u[j]);
      }
      cnt+=this->domainG->neqSub[i];
    }
//    double *mpi_buff;
    for (int i = 1;i < this->domainG->n_clusters;i++){
      neqClst_i=this->domainG->vec_neqClst[i];
      //printf("neqClst_i=%d\n",neqClst_i);
      double *mpi_buff =new double[neqClst_i];
      MPI_Recv(mpi_buff,neqClst_i, MPI_DOUBLE,i, 1, MPI_COMM_WORLD,&mpi_stat);
      for (int j = 0; j<neqClst_i;j++){
        this->u_restrict[l2g_vector[i][j]]=mpi_buff[j];
      }
     delete [] mpi_buff;mpi_buff=NULL;
    }
  }



  // gathering of u_restrict to    rank=0 ||| END   ||| - - - - - - - -
  if (!MPIrank) printf("VTK is creating: sending of elements to master \n");
  //
  // gathering of mesh to rank=0          ||| START ||| - - - - - - - -
  if (MPIrank==0){
    int numel_glob_el_table=0;
    for (int i = 0;i<domainG->n_clusters;i++){
      numel_glob_el_table+=domainG->vec_numelClst[i];
    }
    elements_all    = new int [numel_glob_el_table];
    ptr_el          = new int [domainG->n_elementsAll+1];
    sbd_ptr_el      = new int [domainG->n_elementsAll+1];
    ind_subdom_ddm  = new int [domainG->n_elementsAll];
    ptr_el[0] = 0;
    sbd_ptr_el[0] = 0;
    int i_ptr_el=0;
    int n_nodesOnEl,cnt_ij=0,cnt_i=1;
    for (int k=0;k<this->domainG->n_subdomOnClust;k++){
      for (int i = 0;i<domainG->n_elementsSub[k];i++){
        n_nodesOnEl = 8;
        for (int j = 0;j<n_nodesOnEl;j++){
          //TODO - longInt
			elements_all[cnt_ij]=this->fem[k]->mesh.element[i].inod_glob[j]; 
 ////         printf("%d,",elements_all[cnt_ij]);
          cnt_ij++;
        }
////        printf("\n");
        ptr_el[cnt_i] = cnt_ij; 
        sbd_ptr_el[cnt_i] = 0;
////        if (!MPIrank) printf("ptr_el[cnt_i]=%d\n",ptr_el[cnt_i]);
        cnt_i++;
        ind_subdom_ddm[i_ptr_el] = 0;
        i_ptr_el++;
      }
    }
//    printf("i_ptr_el=%d\n",i_ptr_el);
    int *mpi_buff;
    int *mpi_buff2;
    int cnt_ptr = 0,tmp_int=cnt_ij;
    int n_elementsClst=0,ptr_el_last=0;
    for (int i = 1;i < domainG->n_clusters;i++){
      int numelClst= domainG->vec_numelClst[i];
      mpi_buff = new int [numelClst];
      MPI_Recv(mpi_buff,numelClst, MPI_INT,i, 1, MPI_COMM_WORLD,&mpi_stat);
      for (int j = 0; j<numelClst;j++){
        elements_all[cnt_ij] = mpi_buff[j];
        cnt_ij++;
      }
      n_elementsClst=domainG->vec_n_elementsClst[i];
      mpi_buff2 = new int [numelClst];
      int cnt_j=0;
      MPI_Recv(mpi_buff2,n_elementsClst+1, MPI_INT,i, 2, MPI_COMM_WORLD,&mpi_stat);
      cnt_ptr+=domainG->vec_n_elementsClst[i];

      for (int j = 0; j<n_elementsClst;j++){
        ptr_el[i_ptr_el] = mpi_buff2[j] + tmp_int;
        sbd_ptr_el[i_ptr_el] = i;
        ind_subdom_ddm[i_ptr_el] = i;
        i_ptr_el++;
      }
//      printf("domainG->n_elementsAll,i_ptr_el =[%d,%d]\n",domainG->n_elementsAll,i_ptr_el);
//      tmp_int += mpi_buff2[n_elementsClst];

      ptr_el_last=mpi_buff2[n_elementsClst-1] + tmp_int;
      tmp_int += ptr_el[i_ptr_el-1]; // mpi_buff2[n_elementsClst];
      //printf("mpi_buff2[n_elementsClst] = %d\n",mpi_buff2[n_elementsClst]);


      delete [] mpi_buff;    mpi_buff  = NULL;
      delete [] mpi_buff2;   mpi_buff2 = NULL;
    }
//    ptr_el[i_ptr_el] = ptr_el_last;


   ptr_el[0]=0;
    for (int k=1;k<domainG->n_elementsAll+1;k++){
      ptr_el[k]=ptr_el[k-1]+8;
//      printf("clst=%d, ptr_el[%d]=%d\n",sbd_ptr_el[k-1],k,ptr_el[k]);
    }

  }
  else {             //n_elementsClst

    int numelClst=0,n_elementsClst=0;       
    for (int k=0;k<this->domainG->n_subdomOnClust;k++){
      numelClst+= domainG->n_elementsSub[k]*8;
      n_elementsClst+=domainG->n_elementsSub[k];
    }
    //printf("numelClst=%d, n_elementsClst=%d\n",numelClst,n_elementsClst);
    int *elements_clst  = new int [numelClst];
    int *ptr_el_clst    = new int [n_elementsClst+1];
    ptr_el_clst[0] = 0;
    int n_nodesOnEl;
    int cnt_ij=0,cnt_i=1;
    for (int k=0;k<this->domainG->n_subdomOnClust;k++){
      for (int i = 0;i<domainG->n_elementsSub[k];i++){
        n_nodesOnEl = 8;
        for (int j = 0;j<n_nodesOnEl;j++){
          //TODO longInt
			elements_clst[cnt_ij]=this->fem[k]->mesh.element[i].inod_glob[j]; 
          cnt_ij++;
        }
        ptr_el_clst[cnt_i] = cnt_ij; 
        cnt_i++;
      }
    }
////    for (int k=0;k<numelClst;k++){printf("%d,",elements_clst[k]); if (!(k%8)) printf("\n");}
    MPI_Send(elements_clst, numelClst, MPI_INT, 0 , 1, MPI_COMM_WORLD);
    MPI_Send(ptr_el_clst,n_elementsClst+1, 
                                    MPI_INT, 0 , 2, MPI_COMM_WORLD);
    delete [] elements_clst;  elements_clst=NULL;
    delete [] ptr_el_clst;    ptr_el_clst  = NULL;
  }




////  if (!MPIrank) {
////    for (int i=0;i<domainG->n_elementsAll;i++){
////      printf("ptr_el[cnt_i]=%d\n",ptr_el[i]);
////    }
////  }





  // gathering of mesh to rank=0          ||| END   ||| - - - - - - - -
  if (!MPIrank) printf("VTK is creating: sending of coordinates to master \n");
  //
  // gathering of coordinates to   rank=0 ||| START ||| - - - - - - - -
  if (MPIrank>0){
//    neqClst = domainG->neqClst;
    double *data_coord= new double[neqClst];
    cnt=0;
    for (int k=0;k<domainG->n_subdomOnClust;k++){
      for (int i=0;i<domainG->n_nodsSub[k];i++){
        data_coord[3*cnt  ] = this->fem[k]->mesh.coordinateSub[i].x;
        data_coord[3*cnt+1] = this->fem[k]->mesh.coordinateSub[i].y;
        data_coord[3*cnt+2] = this->fem[k]->mesh.coordinateSub[i].z;
        cnt++;
      } 
    }
    MPI_Send(&data_coord[0], 3*cnt, MPI_DOUBLE, 0 , 1, MPI_COMM_WORLD);
    delete [] data_coord;data_coord = NULL;
  }
  else {
    coordinatesOnMaster = new double [domainG->neqAll]; 
    cnt=0;
    for (int k=0;k<domainG->n_subdomOnClust;k++){
      for (int j=0;j<domainG->n_nodsSub[k];j++){
        coordinatesOnMaster[l2g_vector[0][3*cnt  ]] = this->fem[k]->mesh.coordinateSub[j].x;
        coordinatesOnMaster[l2g_vector[0][3*cnt+1]] = this->fem[k]->mesh.coordinateSub[j].y;
        coordinatesOnMaster[l2g_vector[0][3*cnt+2]] = this->fem[k]->mesh.coordinateSub[j].z;
        cnt++;
      }
    }
    double *mpi_buff;
    for (int i = 1;i < domainG->n_clusters;i++){
      neqClst_i= domainG->vec_neqClst[i];
      mpi_buff = new double[neqClst_i];
      MPI_Recv(mpi_buff,neqClst_i, MPI_DOUBLE,i, 1, MPI_COMM_WORLD,&mpi_stat);
      for (int j = 0; j<neqClst_i;j++){
        coordinatesOnMaster[l2g_vector[i][j]]=mpi_buff[j];
      }
      delete [] mpi_buff; mpi_buff = NULL;
    }
  }
  if (!MPIrank) printf("VTK is creating: all data are on the master \n");
  // gathering of coordinates to   rank=0 ||| END   ||| - - - - - - - -
  //
//########$########$###########################$########$########
//########$########$#### #### #      # ### ####$########$########
//########$########$#### ### #### #### ## #####$########$########
//########$########$#### ## #### ####   #######$########$########
//########$########$#### # ##### #### ## ######$########$########
//########$########$####  ##### #### #### #####$########$########
//########$########$###########################$########$########
  if (MPIrank==0){
#ifndef WIN32
    mode_t mode = 0777;
    mkdir("data", mode);
#endif
    if (!MPIrank) printf("VTK is creating: printing into file\n");
    FILE *fVTK = fopen("data/box_new.vtk", "w");
    fprintf(fVTK, "# vtk DataFile Version 3.0\n");
    fprintf(fVTK, "vtk output\n");
    fprintf(fVTK, "ASCII\n\n");
    fprintf(fVTK, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(fVTK, "POINTS %d float\n",domainG->n_nodsAll);
    for (int i = 0; i <domainG->n_nodsAll; i++) {
      fprintf(fVTK, "%f %f %f\n", coordinatesOnMaster[3*i+0],
                                  coordinatesOnMaster[3*i+1],
                                  coordinatesOnMaster[3*i+2]);
    }
    int k1 = domainG->n_elementsAll * 9;
    fprintf(fVTK, "CELLS %d %d\n",domainG->n_elementsAll, k1);
    for (int i = 0; i < domainG->n_elementsAll; i++) {
      fprintf(fVTK, "%d ", ptr_el[i+1]-ptr_el[i]);//fprintf(fVTK, "%d ", 8);
      for (int j = ptr_el[i]; j < ptr_el[i+1]; j++) {
        fprintf(fVTK, "% d", elements_all[j]);
      }
      fprintf(fVTK, "\n");
    }
    fprintf(fVTK, "CELL_TYPES %d\n", domainG->n_elementsAll);
    for (int i = 0; i < domainG->n_elementsAll; i++) {
      fprintf(fVTK, "%d\n", 12);
    }
    fprintf(fVTK, "POINT_DATA %d\n",domainG->n_nodsAll);
    fprintf(fVTK, "SCALARS displacements float 3\n");
    fprintf(fVTK, "LOOKUP_TABLE my_table\n");
    for (int i = 0; i < domainG->n_nodsAll; i++) {
      fprintf(fVTK, "%f %f %f\n", this->u_restrict[3 * i + 0], 
                                  this->u_restrict[3 * i + 1],  
                                  this->u_restrict[3 * i + 2]);
    }
    fprintf(fVTK, "\nCELL_DATA %d\n",domainG->n_elementsAll);
    fprintf(fVTK, "SCALARS decomposition int 1\n");
    fprintf(fVTK, "LOOKUP_TABLE decomposition\n");
    for (int i = 0; i < domainG->n_elementsAll; i++) {
      fprintf(fVTK, "%d\n", ind_subdom_ddm[i]+1);
    }
    fclose(fVTK);
  }
  if (!elements_all)  { delete [] elements_all;    elements_all=NULL;  }   
  if (!ptr_el)        { delete [] ptr_el;          ptr_el=NULL;        }
  if (!sbd_ptr_el)    { delete [] sbd_ptr_el;      sbd_ptr_el=NULL;    }
  if (!ind_subdom_ddm){ delete [] ind_subdom_ddm;  ind_subdom_ddm=NULL;}
  if (!MPIrank) printf("VTK is creating: finish\n");
}



void CClust_g::createVTK_per_cluster(){

  int MPIrank, MPIsize,cnt,cnt1,cnt2;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &MPIrank);
  MPI_Comm_size(comm, &MPIsize);

  if (this->domainG->vtk_min_rank <= MPIrank && this->domainG->vtk_max_rank >= MPIrank){
#ifndef WIN32
    mode_t mode = 0777;
    mkdir("data", mode);
#endif
    if (!MPIrank){
       printf("---------------------------------------------------------------------\n");
       printf("---------------------------------------------------------------------\n");
       printf("---------------------------------------------------------------------\n");
       printf("VTK is creating: printing into file\n");
    } 

    char filenameKe[128]; 
    sprintf(filenameKe, "data/VTK_%d.vtk", MPIrank);

    FILE *fVTK = fopen(filenameKe, "w");
    fprintf(fVTK, "# vtk DataFile Version 3.0\n");
    fprintf(fVTK, "vtk output\n");
    fprintf(fVTK, "ASCII\n\n");
    fprintf(fVTK, "DATASET UNSTRUCTURED_GRID\n");

//--------------------------------------------------------------------------    
//                        +++ MY_COMM_LAMBDAS_INDICES +++           start
//--------------------------------------------------------------------------    

    vector < int > my_neighs; // my neighboring domains 
    my_neighs.insert(my_neighs.begin(),	this->neighbClst, this->neighbClst + this->domainG->n_neighbClst); 
    std::sort (my_neighs.begin(), my_neighs.end());

    vector < vector <int> >		my_comm_lambdas_indices_comp;
    my_comm_lambdas_indices_comp.resize(this->domainG->n_neighbClst);
    int max_dom_index=MPIsize;
    vector < vector <int> > lambdas_per_subdomain (max_dom_index);
    vector <int> my_lamdas_indices;

    for (int i = 0; i < domainG->lambda_map_sub.size(); i++) 
      my_lamdas_indices.push_back(domainG->lambda_map_sub[i][0]);

    for (int i = 0; i < domainG->lambda_map_sub.size(); i++) {
      if ( domainG->lambda_map_sub[i].size() > 2 ) {
        if ( domainG->lambda_map_sub[i][1] < domainG->lambda_map_sub[i][2] )
         lambdas_per_subdomain[domainG->lambda_map_sub[i][1]].push_back(domainG->lambda_map_sub[i][0]);
         lambdas_per_subdomain[domainG->lambda_map_sub[i][2]].push_back(domainG->lambda_map_sub[i][0]);
      }
    }

    vector < vector <int> >		my_comm_lambdas_indices;
    my_comm_lambdas_indices .resize(my_neighs.size());

    for (int i = 0; i < my_neighs.size(); i++) {
      my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]]; 
    }
    // mapping/compression vector for cluster 
    map <int,int> my_lamdas_map_indices; 
    for (int i = 0; i <my_lamdas_indices.size(); i++)
      my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i)); 
    //// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
    my_comm_lambdas_indices_comp.resize(my_neighs.size());
    for (int i = 0; i < my_neighs.size(); i++) {
      my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
      for (int j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ ) 
        my_comm_lambdas_indices_comp[i][j] = my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]]; 
    }
    //// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *
  //
  //
  //
  //
  //uk
    int  n_elementsClst=0,n_nodsClst=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      n_nodsClst+=domainG->n_nodsSub[i];
      n_elementsClst+=domainG->n_elementsSub[i];
    }

    int n_nodsForVisualofB=0;
    for (int i=0;i<domainG->lambda_map_sub.size();i++){
      if (domainG->lambda_map_sub[i].size()==3){
        n_nodsForVisualofB++;
      }
    }

    int *arrows_B = new int[2*n_nodsForVisualofB];
    double *u_ForVisualofB=new double[3*n_nodsForVisualofB];
    cnt=0;cnt1=0;cnt2=0;
    printf("n_nodsClst = %d\n",n_nodsClst);


    vector <int> tmp_vec;
    int tmp_int,int3;

    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->data[i]->B->BI_full.size();j++){
        tmp_int = my_lamdas_map_indices[this->data[i]->B->BI_full[j]];
        tmp_vec=this->domainG->lambda_map_sub[tmp_int];
        if (tmp_vec.size()==3){
          int3=(int)floorf(this->data[i]->B->BJ_full[j]/3);
          arrows_B[cnt]   = int3 + cnt2;
          arrows_B[cnt+1] = n_nodsClst+cnt1;
          u_ForVisualofB[3*cnt1+0] = this->data[i]->u[3 * int3 + 0];
          u_ForVisualofB[3*cnt1+1] = this->data[i]->u[3 * int3 + 1];
          u_ForVisualofB[3*cnt1+2] = this->data[i]->u[3 * int3 + 2];
          cnt1++;
          cnt+=2;
        }
      }
      cnt2+=domainG->n_nodsSub[i];
    }

    int Isub,Jsub,Ksub;
    int Iclst,Jclst,Kclst;
    int globoalIndSub,globoalIndClst=MPIrank;
    //
    fprintf(fVTK, "POINTS %d float\n",n_nodsClst + n_nodsForVisualofB);
    double xT,yT,zT;


// TODO ad SCALe parameter into 'domainG'
//    double scale=0.80;
    double scale=0.80;



    getMeshCoordsOfClust(globoalIndClst,&Iclst,&Jclst,&Kclst);


    for (int i = 0; i <domainG->n_subdomOnClust; i++) {

////      globoalIndSub=this->fem[i]->mesh.i_domGlobalMeshGen;
////      getMeshCoordsOfSubdom(globoalIndSub,&Isub,&Jsub,&Ksub);
////
////      xT = (0.5+Isub)*this->domainG->Lx/this->domainG->Nx;
////      yT = (0.5+Jsub)*this->domainG->Ly/this->domainG->Ny;
////      zT = (0.5+Ksub)*this->domainG->Lz/this->domainG->Nz;


      xT = (0.5+Iclst)*this->domainG->Lx/this->domainG->Cx;
      yT = (0.5+Jclst)*this->domainG->Ly/this->domainG->Cy;
      zT = (0.5+Kclst)*this->domainG->Lz/this->domainG->Cz;


      for (int j=0;j<this->domainG->n_nodsSub[i];j++){
        fprintf(fVTK, "%f %f %f\n",scale*(this->fem[i]->mesh.coordinateSub[j].x-xT)+xT,
                                   scale*(this->fem[i]->mesh.coordinateSub[j].y-yT)+yT,
                                   scale*(this->fem[i]->mesh.coordinateSub[j].z-zT)+zT);
      }
    }

    int Isub0,Jsub0,Ksub0;
    int Iclst0,Jclst0,Kclst0;
    int globoalIndSub0,indNod;
    int globoalIndClst0;
    double xx,yy,zz,xT0,yT0,zT0;
    int _n_nodsForVisualofB=0;


    globoalIndClst=MPIrank;
    getMeshCoordsOfClust(globoalIndClst,&Iclst,&Jclst,&Kclst);

    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
////      globoalIndSub=this->fem[i]->mesh.i_domGlobalMeshGen;
////      getMeshCoordsOfSubdom(globoalIndSub,&Isub,&Jsub,&Ksub);
      for (int j=0;j<this->data[i]->B->BI_full.size();j++){
        tmp_int = my_lamdas_map_indices[this->data[i]->B->BI_full[j]];
        tmp_vec=this->domainG->lambda_map_sub[tmp_int];
//        for (int k=0;k<tmp_vec.size();k++){
//          printf("tmp[%d]=%d",k,tmp_vec[k]);
//        }
//        printf("\n");
        if (tmp_vec.size()==3){
          _n_nodsForVisualofB++;
          indNod=(int)floorf(this->data[i]->B->BJ_full[j]/3);
          globoalIndClst0=tmp_vec[2];
////          globoalIndSub0=tmp_vec[2];
          getMeshCoordsOfClust(globoalIndClst0,&Iclst0,&Jclst0,&Kclst0);
          xT0 = (0.5+Iclst0)*this->domainG->Lx/this->domainG->Cx;
          yT0 = (0.5+Jclst0)*this->domainG->Ly/this->domainG->Cy;
          zT0 = (0.5+Kclst0)*this->domainG->Lz/this->domainG->Cz;
          xx=this->fem[i]->mesh.coordinateSub[indNod].x;
          yy=this->fem[i]->mesh.coordinateSub[indNod].y;
          zz=this->fem[i]->mesh.coordinateSub[indNod].z;

 //         fprintf(fVTK, "%f %f %f\n",((scale*(xx-xT0)+xT0)),
 //                                    ((scale*(yy-yT0)+yT0)),
 //                                    ((scale*(zz-zT0)+zT0)));

          fprintf(fVTK, "%f %f %f\n",0.5*((scale*(xx-xT)+xT)+(scale*(xx-xT0)+xT0)),
                                     0.5*((scale*(yy-yT)+yT)+(scale*(yy-yT0)+yT0)),
                                     0.5*((scale*(zz-zT)+zT)+(scale*(zz-zT0)+zT0)));
        }
      }
    }

    int n_DP_pointsOnClst=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      n_DP_pointsOnClst+=int(this->fem[i]->mesh.DP_DOFsAll.size()/3);
    }
    printf("n_DP_pointsOnClst=%d\n",n_DP_pointsOnClst);
    if (!MPIrank) printf("++++++++++++++++++++ %d %d \n",n_nodsForVisualofB,_n_nodsForVisualofB);
    int k1 = n_elementsClst * 9 + n_DP_pointsOnClst*2 + n_nodsForVisualofB*3;
    fprintf(fVTK, "CELLS %d %d\n",n_elementsClst+n_DP_pointsOnClst+n_nodsForVisualofB, k1);
    int cnt_el=0;
    tmp_int=0;
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
        fprintf(fVTK, "%d ", 8);
        for (int k=0;k<8;k++){
          tmp_int=this->fem[i]->mesh.g2l_nodes[this->fem[i]->mesh.element[j].inod_glob[k]];
          tmp_int+=cnt_el;
          fprintf(fVTK, "% d",tmp_int);
        }
        fprintf(fVTK, "\n");
      }
      cnt_el+=domainG->n_nodsSub[i];
    }


    cnt_el=0;
    int *DP_nodes = new int[n_DP_pointsOnClst];
    cnt1=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      for (int j=0;j<int(fem[i]->mesh.DP_DOFsAll.size()/3);j++){
        tmp_int=int(fem[i]->mesh.DP_DOFsAll[3*j]/3);
        tmp_int=fem[i]->mesh.g2l_nodes[tmp_int] +
                cnt_el;
        fprintf(fVTK, "%d %d \n",1,tmp_int);
        DP_nodes[cnt1]=tmp_int;
        cnt1++;
      }
      cnt_el+=domainG->n_nodsSub[i];
    }

    for (int i=0;i<n_nodsForVisualofB;i++){
        fprintf(fVTK, "%d %d %d \n",2,arrows_B[2*i],arrows_B[2*i+1]);
    }
    delete [] arrows_B;



    fprintf(fVTK, "CELL_TYPES %d\n", n_elementsClst+n_DP_pointsOnClst+n_nodsForVisualofB);
    for (int i = 0; i < n_elementsClst; i++) {
      fprintf(fVTK, "%d\n", 12);
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", 1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", 3);
    }
    fprintf(fVTK, "POINT_DATA %d\n",n_nodsClst + n_nodsForVisualofB);
    fprintf(fVTK, "SCALARS displacements float 3\n");
    fprintf(fVTK, "LOOKUP_TABLE my_table\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_nodsSub[i];j++){
        fprintf(fVTK, "%f %f %f\n", this->data[i]->u[3 * j + 0], 
                                    this->data[i]->u[3 * j + 1],  
                                    this->data[i]->u[3 * j + 2]);
      }
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%f %f %f\n", u_ForVisualofB[3*i+0],
                                  u_ForVisualofB[3*i+1],
                                  u_ForVisualofB[3*i+2]);
    }

    delete [] u_ForVisualofB;


    fprintf(fVTK, "SCALARS radius float 1\n");
    fprintf(fVTK, "LOOKUP_TABLE default\n");
    cnt1=0;
    for (int i=0;i<n_nodsClst;i++){
      if (i==DP_nodes[cnt1]){
        fprintf(fVTK, "%f \n",1.0); 
        cnt1++;
      }
      else
      {
        fprintf(fVTK, "%f \n",0.0); 
      }
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%f \n",0.0); 
    }

    delete [] DP_nodes;
// displacements
    fprintf(fVTK, "\nCELL_DATA %d\n",n_elementsClst+n_nodsForVisualofB+n_DP_pointsOnClst);
    fprintf(fVTK, "SCALARS decomposition int 1\n");
    fprintf(fVTK, "LOOKUP_TABLE decomposition\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
      fprintf(fVTK, "%d\n",i+1);
      }
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", -1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", -2);
    }

    fprintf(fVTK, "SCALARS decomposition_clust int 1\n");
    fprintf(fVTK, "LOOKUP_TABLE default\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
      fprintf(fVTK, "%d\n",MPIrank);
      }
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", -1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", -2);
    }







    fclose(fVTK);
  }

}



void CClust_g::createVTK_per_cluster_new(){

  int MPIrank, MPIsize,cnt,cnt1,cnt2;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &MPIrank);
  MPI_Comm_size(comm, &MPIsize);

//  if (this->domainG->vtk_min_rank<MPIrank && this->domainG->vtk_max_rank > MPIrank){
#ifndef WIN32
    mode_t mode = 0777;
    mkdir("data", mode);
#endif
    if (!MPIrank){
       printf("---------------------------------------------------------------------\n");
       printf("---------------------------------------------------------------------\n");
       printf("---------------------------------------------------------------------\n");
       printf("VTK is creating: printing into file\n");
    } 

    char filenameKe[128]; 
    sprintf(filenameKe, "data/VTK_%d.vtk", MPIrank);

    FILE *fVTK = fopen(filenameKe, "w");
    fprintf(fVTK, "# vtk DataFile Version 3.0\n");
    fprintf(fVTK, "vtk output\n");
    fprintf(fVTK, "ASCII\n\n");
    fprintf(fVTK, "DATASET UNSTRUCTURED_GRID\n");

//--------------------------------------------------------------------------    
//                        +++ MY_COMM_LAMBDAS_INDICES +++           start
//--------------------------------------------------------------------------    

    vector < int > my_neighs; // my neighboring domains 
    my_neighs.insert(my_neighs.begin(),	this->neighbClst, this->neighbClst + this->domainG->n_neighbClst); 
    std::sort (my_neighs.begin(), my_neighs.end());

    vector < vector <int> >		my_comm_lambdas_indices_comp;
    my_comm_lambdas_indices_comp.resize(this->domainG->n_neighbClst);
    int max_dom_index=MPIsize;
    vector < vector <int> > lambdas_per_subdomain (max_dom_index);
    vector <int> my_lamdas_indices;

    for (int i = 0; i < domainG->lambda_map_sub.size(); i++) 
      my_lamdas_indices.push_back(domainG->lambda_map_sub[i][0]);

    for (int i = 0; i < domainG->lambda_map_sub.size(); i++) {
      if ( domainG->lambda_map_sub[i].size() > 2 ) {
        if ( domainG->lambda_map_sub[i][1] < domainG->lambda_map_sub[i][2] )
         lambdas_per_subdomain[domainG->lambda_map_sub[i][1]].push_back(domainG->lambda_map_sub[i][0]);
         lambdas_per_subdomain[domainG->lambda_map_sub[i][2]].push_back(domainG->lambda_map_sub[i][0]);
      }
    }

    vector < vector <int> >		my_comm_lambdas_indices;
    my_comm_lambdas_indices .resize(my_neighs.size());

    for (int i = 0; i < my_neighs.size(); i++) {
      my_comm_lambdas_indices[i] = lambdas_per_subdomain[my_neighs[i]]; 
    }
    // mapping/compression vector for cluster 
    map <int,int> my_lamdas_map_indices; 
    for (int i = 0; i <my_lamdas_indices.size(); i++)
      my_lamdas_map_indices.insert(make_pair(my_lamdas_indices[i],i)); 
    //// *** Create a vector of communication pattern needed for AllReduceLambdas function *******
    my_comm_lambdas_indices_comp.resize(my_neighs.size());
    for (int i = 0; i < my_neighs.size(); i++) {
      my_comm_lambdas_indices_comp[i].resize( lambdas_per_subdomain[my_neighs[i]].size() );
      for (int j = 0; j < lambdas_per_subdomain[my_neighs[i]].size(); j++ ) 
        my_comm_lambdas_indices_comp[i][j] = my_lamdas_map_indices[lambdas_per_subdomain[my_neighs[i]][j]]; 
    }
    //// *** END - Create a vector of communication pattern needed for AllReduceLambdas function *
  //
  //
  //
  //
  //uk
    int  n_elementsClst=0,n_nodsClst=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      n_nodsClst+=domainG->n_nodsSub[i];
      n_elementsClst+=domainG->n_elementsSub[i];
    }

    int n_nodsForVisualofB=0;
    for (int i=0;i<domainG->lambda_map_sub.size();i++){
      if (domainG->lambda_map_sub[i].size()==3){
        n_nodsForVisualofB++;
      }
    }

    int *arrows_B = new int[2*n_nodsForVisualofB];
    double *u_ForVisualofB=new double[3*n_nodsForVisualofB];
    cnt=0;cnt1=0;cnt2=0;
    printf("n_nodsClst = %d\n",n_nodsClst);


    vector <int> tmp_vec;
    int tmp_int,int3;

    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->data[i]->B->BI_full.size();j++){
        tmp_int = my_lamdas_map_indices[this->data[i]->B->BI_full[j]];
        tmp_vec=this->domainG->lambda_map_sub[tmp_int];
        if (tmp_vec.size()==3){
          int3=(int)floorf(this->data[i]->B->BJ_full[j]/3);
          arrows_B[cnt]   = int3 + cnt2;
          arrows_B[cnt+1] = n_nodsClst+cnt1;
          u_ForVisualofB[3*cnt1+0] = this->data[i]->u[3 * int3 + 0];
          u_ForVisualofB[3*cnt1+1] = this->data[i]->u[3 * int3 + 1];
          u_ForVisualofB[3*cnt1+2] = this->data[i]->u[3 * int3 + 2];
          cnt1++;
          cnt+=2;
        }
      }
      cnt2+=domainG->n_nodsSub[i];
    }

    int Isub,Jsub,Ksub;
    int Iclst,Jclst,Kclst;
    int globoalIndSub,globoalIndClst=MPIrank;
    //
    fprintf(fVTK, "POINTS %d float\n",n_nodsClst + n_nodsForVisualofB);
    double xT,yT,zT;


// TODO ad SCALe parameter into 'domainG'
//    double scale=0.80;
    double scale=0.80;



    getMeshCoordsOfClust(globoalIndClst,&Iclst,&Jclst,&Kclst);


    for (int i = 0; i <domainG->n_subdomOnClust; i++) {

////      globoalIndSub=this->fem[i]->mesh.i_domGlobalMeshGen;
////      getMeshCoordsOfSubdom(globoalIndSub,&Isub,&Jsub,&Ksub);
////
////      xT = (0.5+Isub)*this->domainG->Lx/this->domainG->Nx;
////      yT = (0.5+Jsub)*this->domainG->Ly/this->domainG->Ny;
////      zT = (0.5+Ksub)*this->domainG->Lz/this->domainG->Nz;


      xT = (0.5+Iclst)*this->domainG->Lx/this->domainG->Cx;
      yT = (0.5+Jclst)*this->domainG->Ly/this->domainG->Cy;
      zT = (0.5+Kclst)*this->domainG->Lz/this->domainG->Cz;


      for (int j=0;j<this->domainG->n_nodsSub[i];j++){
        fprintf(fVTK, "%f %f %f\n",scale*(this->fem[i]->mesh.coordinateSub[j].x-xT)+xT,
                                   scale*(this->fem[i]->mesh.coordinateSub[j].y-yT)+yT,
                                   scale*(this->fem[i]->mesh.coordinateSub[j].z-zT)+zT);
      }
    }


    int Isub0,Jsub0,Ksub0;
    int Iclst0,Jclst0,Kclst0;
    int globoalIndSub0,indNod;
    int globoalIndClst0;
    double xx,yy,zz,xT0,yT0,zT0;
    int _n_nodsForVisualofB=0;


    globoalIndClst=MPIrank;
    getMeshCoordsOfClust(globoalIndClst,&Iclst,&Jclst,&Kclst);

    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
////      globoalIndSub=this->fem[i]->mesh.i_domGlobalMeshGen;
////      getMeshCoordsOfSubdom(globoalIndSub,&Isub,&Jsub,&Ksub);
      for (int j=0;j<this->data[i]->B->BI_full.size();j++){
        tmp_int = my_lamdas_map_indices[this->data[i]->B->BI_full[j]];
        tmp_vec=this->domainG->lambda_map_sub[tmp_int];
//        for (int k=0;k<tmp_vec.size();k++){
//          printf("tmp[%d]=%d",k,tmp_vec[k]);
//        }
//        printf("\n");
        if (tmp_vec.size()==3){
          _n_nodsForVisualofB++;
          indNod=(int)floorf(this->data[i]->B->BJ_full[j]/3);
          globoalIndClst0=tmp_vec[2];
////          globoalIndSub0=tmp_vec[2];
          getMeshCoordsOfClust(globoalIndClst0,&Iclst0,&Jclst0,&Kclst0);
          xT0 = (0.5+Iclst0)*this->domainG->Lx/this->domainG->Cx;
          yT0 = (0.5+Jclst0)*this->domainG->Ly/this->domainG->Cy;
          zT0 = (0.5+Kclst0)*this->domainG->Lz/this->domainG->Cz;
          xx=this->fem[i]->mesh.coordinateSub[indNod].x;
          yy=this->fem[i]->mesh.coordinateSub[indNod].y;
          zz=this->fem[i]->mesh.coordinateSub[indNod].z;

 //         fprintf(fVTK, "%f %f %f\n",((scale*(xx-xT0)+xT0)),
 //                                    ((scale*(yy-yT0)+yT0)),
 //                                    ((scale*(zz-zT0)+zT0)));

          fprintf(fVTK, "%f %f %f\n",0.5*((scale*(xx-xT)+xT)+(scale*(xx-xT0)+xT0)),
                                     0.5*((scale*(yy-yT)+yT)+(scale*(yy-yT0)+yT0)),
                                     0.5*((scale*(zz-zT)+zT)+(scale*(zz-zT0)+zT0)));
        }
      }
    }

    int n_DP_pointsOnClst=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      n_DP_pointsOnClst+=int(this->fem[i]->mesh.DP_DOFsAll.size()/3);
    }
    printf("n_DP_pointsOnClst=%d\n",n_DP_pointsOnClst);
    if (!MPIrank) printf("++++++++++++++++++++ %d %d \n",n_nodsForVisualofB,_n_nodsForVisualofB);
    int k1 = n_elementsClst * 9 + n_DP_pointsOnClst*2 + n_nodsForVisualofB*3;
    fprintf(fVTK, "CELLS %d %d\n",n_elementsClst+n_DP_pointsOnClst+n_nodsForVisualofB, k1);
    int cnt_el=0;
    tmp_int=0;
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
        fprintf(fVTK, "%d ", 8);
        for (int k=0;k<8;k++){
          tmp_int=this->fem[i]->mesh.g2l_nodes[this->fem[i]->mesh.element[j].inod_glob[k]];
          tmp_int+=cnt_el;
          fprintf(fVTK, "% d",tmp_int);
        }
        fprintf(fVTK, "\n");
      }
      cnt_el+=domainG->n_nodsSub[i];
    }


    cnt_el=0;
    int *DP_nodes = new int[n_DP_pointsOnClst];
    cnt1=0;
    for (int i=0;i<domainG->n_subdomOnClust;i++){
      for (int j=0;j<int(fem[i]->mesh.DP_DOFsAll.size()/3);j++){
        tmp_int=int(fem[i]->mesh.DP_DOFsAll[3*j]/3);
        tmp_int=fem[i]->mesh.g2l_nodes[tmp_int] +
                cnt_el;
        fprintf(fVTK, "%d %d \n",1,tmp_int);
        DP_nodes[cnt1]=tmp_int;
        cnt1++;
      }
      cnt_el+=domainG->n_nodsSub[i];
    }

    for (int i=0;i<n_nodsForVisualofB;i++){
        fprintf(fVTK, "%d %d %d \n",2,arrows_B[2*i],arrows_B[2*i+1]);
    }
    delete [] arrows_B;



    fprintf(fVTK, "CELL_TYPES %d\n", n_elementsClst+n_DP_pointsOnClst+n_nodsForVisualofB);
    for (int i = 0; i < n_elementsClst; i++) {
      fprintf(fVTK, "%d\n", 12);
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", 1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", 3);
    }
    fprintf(fVTK, "POINT_DATA %d\n",n_nodsClst + n_nodsForVisualofB);
    fprintf(fVTK, "SCALARS displacements float 3\n");
    fprintf(fVTK, "LOOKUP_TABLE my_table\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_nodsSub[i];j++){
        fprintf(fVTK, "%f %f %f\n", this->data[i]->u[3 * j + 0], 
                                    this->data[i]->u[3 * j + 1],  
                                    this->data[i]->u[3 * j + 2]);
      }
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%f %f %f\n", u_ForVisualofB[3*i+0],
                                  u_ForVisualofB[3*i+1],
                                  u_ForVisualofB[3*i+2]);
    }

    delete [] u_ForVisualofB;


    fprintf(fVTK, "SCALARS radius float 1\n");
    fprintf(fVTK, "LOOKUP_TABLE default\n");
    cnt1=0;
    for (int i=0;i<n_nodsClst;i++){
      if (i==DP_nodes[cnt1]){
        fprintf(fVTK, "%f \n",1.0); 
        cnt1++;
      }
      else
      {
        fprintf(fVTK, "%f \n",0.0); 
      }
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%f \n",0.0); 
    }

    delete [] DP_nodes;
// displacements
    fprintf(fVTK, "\nCELL_DATA %d\n",n_elementsClst+n_nodsForVisualofB+n_DP_pointsOnClst);
    fprintf(fVTK, "SCALARS decomposition int 1\n");
    fprintf(fVTK, "LOOKUP_TABLE decomposition\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
      fprintf(fVTK, "%d\n",i+1);
      }
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", -1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", -2);
    }

    fprintf(fVTK, "SCALARS decomposition_clust int 1\n");
    fprintf(fVTK, "LOOKUP_TABLE default\n");
    for (int i = 0; i <domainG->n_subdomOnClust; i++) {
      for (int j=0;j<this->domainG->n_elementsSub[i];j++){
      fprintf(fVTK, "%d\n",MPIrank);
      }
    }
    for (int i=0;i<n_DP_pointsOnClst;i++){
      fprintf(fVTK, "%d\n", -1);
    }
    for (int i=0;i<n_nodsForVisualofB;i++){
      fprintf(fVTK, "%d\n", -2);
    }







    fclose(fVTK);
  //}

}



void CClust_g::getMeshCoordsOfSubdom(int globoalIndSub,
                                     int *Isub, int *Jsub, int *Ksub){
  int isub,jsub,ksub;
  ksub = ceil(double (globoalIndSub+1)/(this->domainG->Nx*this->domainG->Ny));
  int inXYplane = (globoalIndSub+1)-(ksub-1)*this->domainG->Nx*this->domainG->Ny;
  jsub = ceil(double( inXYplane)/this->domainG->Nx);
  isub = inXYplane - this->domainG->Nx*(jsub-1);
  *Isub=isub-1;
  *Jsub=jsub-1;
  *Ksub=ksub-1;

}

void CClust_g::getMeshCoordsOfClust(int globoalIndClst,
                                    int *I, int *J, int *K){
  int i,j,k;
  k= ceil(double (globoalIndClst+1)/(this->domainG->Cx*this->domainG->Cy));
  int inXYplane = (globoalIndClst+1)-(k-1)*this->domainG->Cx*this->domainG->Cy;
  j= ceil(double( inXYplane)/this->domainG->Cx);
  i= inXYplane - this->domainG->Cx*(j-1);
  *I=i-1;
  *J=j-1;
  *K=k-1;

}
