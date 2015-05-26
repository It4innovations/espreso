
#include "BSparseC.h"

static int event1=0;

CBSparseC::CBSparseC(MPI_Comm comm) 
{
	this->comm = comm;
	i_eq = NULL;
	j_eq = NULL;
	v_eq = NULL;

	i_bg = NULL;
	j_bg = NULL;
	v_bg = NULL; 
	i_bd = NULL;
	j_bd = NULL;
	v_bd = NULL;
	i_bc = NULL;
	i_bc_crs = NULL;
	j_bc = NULL;
	v_bc = NULL; 
	multiplicityPD = NULL;
	multiplicityDD = NULL;
	n_row_bg = 0;
	nnz_bg   = 0;
	n_row_bd = 0;
	nnz_bd   = 0;
	n_row_bc = 0;
	nnz_bc   = 0;
	n_row_eq = 0;
	n_col    = 0;

	Bi_coo = NULL; 
	Bj_coo = NULL;
	Bv_coo = NULL;
}

CBSparseC::~CBSparseC() {
  cout << "BSparseC destructor " <<endl;
  /* meaningless ??? */
	BI.clear();  BI.push_back(0);
	BJ.clear();  BJ.push_back(0);
	BV.clear();  BV.push_back(0);

	BI_l.clear();  BI_l.push_back(0);
	BJ_l.clear();  BJ_l.push_back(0);
	BV_l.clear();  BV_l.push_back(0);

    lambda_map_sub.clear();
	
	if (i_eq) {delete [] i_eq; i_eq=NULL;}  
	if (j_eq) {delete [] j_eq; j_eq=NULL;}
	if (v_eq) {delete [] v_eq; v_eq=NULL;}
	if (i_bc) {delete [] i_bc; i_bc=NULL;}
	if (j_bc) {delete [] j_bc; j_bc=NULL;}
	if (v_bc) {delete [] v_bc; v_bc=NULL;}
	if (i_bd) {delete [] i_bd; i_bd=NULL;}
	if (j_bd) {delete [] j_bd; j_bd=NULL;}
	if (v_bd) {delete [] v_bd; v_bd=NULL;}

	if (multiplicityPD)   {delete [] multiplicityPD;    multiplicityPD=NULL;   }
	if (multiplicityDD)   {delete [] multiplicityDD;    multiplicityDD=NULL;   }
	if (Bi_coo) 	        {delete [] Bi_coo;            Bi_coo=NULL;           }

}

int CBSparseC::coo2crs(int m, int nnz, int *row_ind, int *row_ptr) {
	int i,row;

	for (i=0,row=0; row<=m; row++) {
		row_ptr[row]=i;
		while (i < nnz && row == row_ind[i]) i++;
	}
	return 0;
}

//void CBSparseC::createB(CFem *fem,CDomain *domainG, int nPrimalDOFs,int *borderDofs) {
//void CBSparseC::createB(CFem *fem,
//                CDomain *domainG, 
//                int nPrimalDOFs,  
//                int *borderDofs) {


bool comp_vf(const vector<int> &a, const vector<int> &b)
{
	return a[0] < b[0]; 
}

	vector<CFem*>  fem;
	vector<CData*>  data;


void CBSparseC::createB1(CDomain * domainG,
                         vector<CFem*> &fem, vector<CData*> &data,
                         int *indExterDOFsOnClust,
                         int n_exterDOFsClust,
                         int *neighbClst,
                         int n_neighbClst,
                         bool flag_redund_lagr_mult){

	int cnt;
	int MPIrank, MPIsize;
	MPI_Comm_rank(comm, &MPIrank);
	MPI_Comm_size(comm, &MPIsize);

	if (!event1) LogEventRegister("createB1",&event1);
	LogEventBegin(event1);

	// Dirichlet BC setup 

	int lo,hi;
	int total_number_of_dirichlet_lambdas = 0; 

	int nLamDir=0;
	for (int i = 0;i<domainG->n_subdomOnClust;i++){
		nLamDir+=fem[i]->bound_cond.dirBCSub->n;
	}

	MPI_Status status;
	lo = 0;
	if (MPIrank > 0){
		MPI_Recv(&lo, 1, MPI_INT, MPIrank-1, 123, MPI_COMM_WORLD, &status);
	}
	
	hi = lo + nLamDir; 
	if (MPIrank < MPIsize-1) {
		MPI_Send(&hi, 1, MPI_INT, MPIrank+1, 123, MPI_COMM_WORLD);
	}

	int strtI = lo;

	int buf;    
	buf = lo + nLamDir;
	MPI_Bcast(&buf, 1, MPI_INT, MPIsize-1, MPI_COMM_WORLD);
	n_row_bd = buf;
	total_number_of_dirichlet_lambdas = n_row_bd; 

	for (int k = 0; k < domainG->n_subdomOnClust; k++) {
		int nLamDirSub=fem[k]->bound_cond.dirBCSub->n;
		data[k]->B->n_bd = nLamDirSub;

		data[k]->B->n_row_bd = n_row_bd;
		data[k]->B->n_col_bd = domainG->neqSub[fem[k]->i_domOnClust];

		data[k]->B->BI_dir.resize(nLamDirSub);
		data[k]->B->BJ_dir.resize(nLamDirSub);
		data[k]->B->BV_dir.resize(nLamDirSub);

	}

	cilk_for (int k = 0; k < domainG->n_subdomOnClust; k++) {
		int nLamDirSub=fem[k]->bound_cond.dirBCSub->n;
		
		for (int i = 0; i < nLamDirSub; i++){
			data[k]->B->BI_dir[i] = (k * nLamDirSub + i) + strtI; //cnt1+strtI; 
			data[k]->B->BJ_dir[i] = fem[k]->bound_cond.dirBCSub->ind[i]; 
			data[k]->B->BV_dir[i] = 1; 
		}    
	}
	// END - Dirichlet BC setup 




	// Global B1 - start 

//	int *Weight = new int [domainG->neqSub[fem->i_domOnClust]];  //X
//	for (int i = 0;i<domainG->neqSub[fem->i_domOnClust];i++){   //X
//		Weight[i] = 1;
//	}


	if (MPIrank == 0) cout << "Creating Global B1 " << endl; 
#ifndef WIN32
	if (MPIrank == 0) { cout << " point 1                                                                  "; system("date +%T.%6N"); }
#endif

	int neighDomNum  = n_neighbClst;      // number of neighboring sub-domains for current sub-domainG 
	int borderDofNum = n_exterDOFsClust; 
	
	vector < vector < longInt > > neighBorderDofs;        // 2D vector for DOFs on borders of neighboring sub-domains 
	neighBorderDofs.resize( neighDomNum, vector<longInt>( 0 , 0 ) );

	vector < longInt > myBorderDOFs;  // my border DOFs are here 
	std::vector< longInt >::iterator it;
	
	myBorderDOFs.insert(myBorderDOFs.begin(), indExterDOFsOnClust, indExterDOFsOnClust + n_exterDOFsClust); //X
	std::sort (myBorderDOFs.begin(), myBorderDOFs.end());
    it = std::unique (myBorderDOFs.begin(), myBorderDOFs.end());
    myBorderDOFs.resize( std::distance(myBorderDOFs.begin(), it) );

	vector < int > myNeighSubDomains; // my neighboring domains 
	myNeighSubDomains.insert(myNeighSubDomains.begin(),	neighbClst, neighbClst+ n_neighbClst); 
	std::sort (myNeighSubDomains.begin(), myNeighSubDomains.end());

	neigh_domains = myNeighSubDomains;

	MPI_Request mpi_req;
	MPI_Status mpi_stat;

	MPI_Request * mpi_send_req  = new MPI_Request [neighDomNum];
	MPI_Request * mpi_recv_req  = new MPI_Request [neighDomNum];
	MPI_Status  * mpi_recv_stat = new MPI_Status  [neighDomNum];

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 2                                                                  "; system("date +%T.%6N"); }
#endif

	// sending my DOFs on border to all neighboring sub-domains 
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		MPI_Isend(&myBorderDOFs[0], myBorderDOFs.size(), MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]); 
	}

	// receiving all border DOFs from all neighboring sub-domains 
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		neighBorderDofs[i].resize(myBorderDOFs.size()); 
		MPI_Irecv(&neighBorderDofs[i][0], myBorderDOFs.size(),MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_recv_req[i]);
	}

	MPI_Waitall(neighDomNum, mpi_recv_req, mpi_recv_stat); 

    delete [] mpi_recv_req;
    delete [] mpi_recv_stat;
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
	//free(mpi_buff);

	// NOW : 
	// my neighboring sub-domains are in : myNeighSubDomains
	// my border DOFs are in : myBorderDOFs
	// neighboring sub-domains border DOFs are in neighBorderDofs[i] 

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 3                                                                  "; system("date +%T.%6N"); }
#endif

 	vector < vector < int > > myNeighsSparse; 
	myNeighsSparse.resize( myBorderDOFs.size() ); //, vector<int>( 2 , 0 ) );

	for (int j = 0; j < myBorderDOFs.size(); j++) { 			
		myNeighsSparse[j].reserve(5);
		myNeighsSparse[j].push_back(myBorderDOFs[j]);
	}
	
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		int d = myNeighSubDomains[i];
		int mydofs_index = 0; 
		int nedofs_index = 0; 
	
		if ( i == 0 && MPIrank < myNeighSubDomains[i] ) {
			for (unsigned int j = 0; j < myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank); 
		} 

		do {	
		
			if ( neighBorderDofs[i][nedofs_index] == myBorderDOFs[mydofs_index] ) {
				myNeighsSparse[mydofs_index].push_back(d);
				mydofs_index++;
				nedofs_index++;
			} else {
				if ( neighBorderDofs[i][nedofs_index] > myBorderDOFs[mydofs_index] ) {
					mydofs_index++;
				} else {
					nedofs_index++;
				}
			}
			
		} while (mydofs_index < myBorderDOFs.size() && nedofs_index < neighBorderDofs[i].size() );
		

		if ( i < myNeighSubDomains.size() - 1) {
			if ( MPIrank > myNeighSubDomains[i] && MPIrank < myNeighSubDomains[i+1] )
				for (unsigned int j = 0; j < myBorderDOFs.size(); j++)
					myNeighsSparse[j].push_back(MPIrank); 
		}


		if ( i == myNeighSubDomains.size()-1 && MPIrank > myNeighSubDomains[i] ) {
			for (unsigned int j = 0; j < myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank); 
		} 


	}

	// temp variable(s) for sort and unique functions 

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 4                                                                  "; system("date +%T.%6N"); }
#endif

	for (unsigned int i = 0; i < myNeighsSparse.size(); i++)
		if (myNeighsSparse[i].size() < 3)
			myNeighsSparse[i].clear(); 
	
	vector < vector < int > > myLambdas;
	//myLambdas.resize( myNeighsSparse.size(), vector<int>( 5 , 0 ) );

	int lambdaCount = 0; 
	for (unsigned int i = 0; i < myNeighsSparse.size(); i++) {
		int add = 0;
		if (myNeighsSparse[i].size() > 0) {
			int dof = myNeighsSparse[i][0];

			//int cnt_tmp = 0; 
			int cnt_tmp = myNeighsSparse[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb
			//for (int t = myNeighsSparse[i].size() - 2; t > 0; t--) 
			//	cnt_tmp += t; 

			for (unsigned int j = 1; j < myNeighsSparse[i].size(); j++) {
				int neighSD = myNeighsSparse[i][j];
		        //domainG.flag_redund_lagr_mult

				if (add == 1) {

					myLambdas.push_back ( vector <int> () );
					myLambdas[lambdaCount].reserve(6);
					myLambdas[lambdaCount].resize(5);
					myLambdas[lambdaCount][0] = lambdaCount;	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter  
					myLambdas[lambdaCount][1] = dof;			// dof in global numbering 
					myLambdas[lambdaCount][2] = MPIrank;		// my sub-domainG 
					myLambdas[lambdaCount][3] = neighSD;		// neigh. sub-domainG 
					myLambdas[lambdaCount][4] = cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
					lambdaCount++; 
					
					if (!flag_redund_lagr_mult) 
						break;
				}

				if (neighSD == MPIrank)
					add = 1; 
			}
		}
	}

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 5                                                                  "; system("date +%T.%6N"); }
#endif

	// Create lambda global numbering 
	int lambdaGlobalCount = 0; 
	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
	if (MPIrank == 0) lambdaGlobalCount = 0; 
	int lambdaNum = lambdaCount + lambdaGlobalCount;  
	MPI_Bcast(&lambdaNum,1,MPI_INT,MPIsize-1,MPI_COMM_WORLD);
	int total_number_of_global_B1_lambdas = lambdaNum; 

	if (myLambdas.size() > 0) 
		//cilk_
		for (unsigned int i = 0; i < myLambdas.size(); i++) 
			myLambdas[i][0] = myLambdas[i][0] + total_number_of_dirichlet_lambdas + lambdaGlobalCount; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index



#ifndef WIN32
	if (MPIrank == 0) { cout << " point 6                                                                  "; system("date +%T.%6N"); }
#endif

	vector < vector < int > > mpi_send_buff;
	mpi_send_buff.resize( myNeighSubDomains.size(), vector<int>( 0 , 0 ) );
	cilk_for (int i = 0; i < myNeighSubDomains.size(); i++) {
		int index = 0;  
		if (myLambdas.size() > 0 )
		{
			for (unsigned int j = 0; j < myLambdas.size(); j++) {
				if( myLambdas[j][3] == myNeighSubDomains[i] ) { 
					mpi_send_buff[i].push_back(myLambdas[j][0]);
					mpi_send_buff[i].push_back(myLambdas[j][1]);
					mpi_send_buff[i].push_back(myLambdas[j][2]);
					mpi_send_buff[i].push_back(myLambdas[j][3]);  
					mpi_send_buff[i].push_back(myLambdas[j][4]);  
					index++; 
				}
			}
			if (index == 0)
				mpi_send_buff[i].push_back(0);
		} 
		else 
		{
			mpi_send_buff[i].push_back(0);
		}
	}

	
#ifndef WIN32
	if (MPIrank == 0) { cout << " point 7                                                                  "; system("date +%T.%6N"); }
#endif	
	
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) 
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
	
	
#ifndef WIN32
	if (MPIrank == 0) { cout << " point 8                                                                  "; system("date +%T.%6N"); }
#endif

	vector < vector < int > > mpi_recv_buff;
	mpi_recv_buff.resize( myNeighSubDomains.size(), vector<int>( 0 , 0 ) );
	delete [] mpi_send_req;

	// receiving all border DOFs from all neighboring sub-domains 
	
	int messages_received = 0; 
	while ( messages_received < myNeighSubDomains.size() ) {
		for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		
			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
			int flag; 		
			MPI_Iprobe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat ); 

			if (flag) {
				int recv_msg_size = 0; 

				MPI_Get_count(&mpi_stat, MPI_INT, &recv_msg_size);
				int* mpi_tmp_recv_buff = (int*)malloc(sizeof(int) * recv_msg_size);

				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);
				//neighBorderDofs[i].insert(neighBorderDofs[i].begin(), mpi_buff, mpi_buff + myBorderDOFs.size());

				free(mpi_tmp_recv_buff);
				messages_received++;
			}
		}
	}
	// decode received lambdas 

	
#ifndef WIN32
	if (MPIrank == 0) { cout << " point 9                                                                  "; system("date +%T.%6N"); }
#endif

	int recv_lamba_count = 0; 
	for (int i = 0; i < myNeighSubDomains.size(); i++) 
		if (mpi_recv_buff[i].size() > 1)
			recv_lamba_count += mpi_recv_buff[i].size(); 
	
	recv_lamba_count = recv_lamba_count / 5; 
	
	int l_i = myLambdas.size(); 
	myLambdas.resize( myLambdas.size() + recv_lamba_count, vector<int>( 5 , 0 ) );
		
	for (int i = 0; i < myNeighSubDomains.size(); i++) {
		if (mpi_recv_buff[i].size() > 1) {
			for (int j = 0; j < mpi_recv_buff[i].size() / 5; j++) {
				myLambdas[l_i][0] = mpi_recv_buff[i][5*j + 0];
				myLambdas[l_i][1] = mpi_recv_buff[i][5*j + 1];
				myLambdas[l_i][2] = mpi_recv_buff[i][5*j + 3];
				myLambdas[l_i][3] = mpi_recv_buff[i][5*j + 2];
				myLambdas[l_i][4] = mpi_recv_buff[i][5*j + 4];
				l_i++;
			}
		}
	}

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 10                                                                 "; system("date +%T.%6N"); }
#endif

#ifdef HFETI_SOLVER
	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
#else
	std::sort         (myLambdas.begin(), myLambdas.end(), comp_vf); 
#endif

#ifndef WIN32
	if (MPIrank == 0) { cout << " point 11                                                                 "; system("date +%T.%6N"); }
#endif

	int lambda; 
	int DOFNumber; 
	int n_myLambdas = myLambdas.size();
	
#ifdef HFETI_SOLVER
	tbb::mutex m;
	cilk_for (int j = 0; j < n_myLambdas; j++) 
#else
	for (int j = 0; j < n_myLambdas; j++)
#endif
	{
		std::map<longInt,int> :: iterator itm;

		int lambda    = myLambdas[j][0];
		int DOFNumber = myLambdas[j][1];
		int cnt       = myLambdas[j][4];
		double B_value;
		if (myLambdas[j][2] < myLambdas[j][3])
			B_value = 1.0; 
		else 
			B_value = -1.0; 

		for (int i = 0; i < data.size(); i++) {
			itm = fem[i]->mesh.g2l.find(DOFNumber);		//( BJ[j] );
			if ( itm != fem[i]->mesh.g2l.end() ) {
			   #ifdef HFETI_SOLVER
				m.lock();
			   #endif
				 data[i]->B->BI.push_back(lambda);		//( BI[j] );
				 data[i]->B->BJ.push_back(itm->second);
				 data[i]->B->BV.push_back(B_value);		//( BV[j] );
				 myLambdas[j].push_back(i);
			   #ifdef HFETI_SOLVER
				m.unlock();
			   #endif
				break; 
			}
		}
	}


	for (int i = 0; i < data.size(); i++){
		data[i]->B->n_row_bg = lambdaNum;
		data[i]->B->n_col_bg = domainG->neqSub[fem[i]->i_domOnClust];
	}


#ifndef WIN32
	if (MPIrank == 0) { cout << " point 12                                                                 "; system("date +%T.%6N"); }
#endif

	n_row_bg = lambdaNum;
	nnz_bg   = BI.size();

	// -------------------------------------------------------------------------


#ifndef WIN32
	if (MPIrank == 0) { cout << "END - Creating Global B1                                                  "; system("date +%T.%6N"); }
#endif

	MPI_Barrier(MPI_COMM_WORLD); 

	if (MPIrank == 0) cout << endl << "Creating Local  B1" << endl; 
#ifndef WIN32
	if (MPIrank == 0) { cout << "BEG - Creating Local  B1                                                  "; system("date +%T.%6N"); }
#endif

// *** Local B1 - begin ******************************************************************  

  vector < vector < longInt > > neighBorderDofs_l;	// 2D vector for DOFs on borders of neighboring sub-domains 
  neighBorderDofs_l.resize( data.size(), vector<longInt>( 0 , 0 ) );

  cilk_for (int ii = 0; ii < data.size(); ii++) {
	  neighBorderDofs_l[ii].insert(neighBorderDofs_l[ii].begin(), // copy neigh border DOFs from original array to local vector 
		  data[ii]->indExterDOFs, 
		  data[ii]->indExterDOFs + domainG->n_exterDOFs[ii]);
	  
	  std::sort (neighBorderDofs_l[ii].begin(), neighBorderDofs_l[ii].end()); // sorting neigh DOFs 

	  std::vector<longInt>::iterator it_l;
	  it_l = std::unique (neighBorderDofs_l[ii].begin(), neighBorderDofs_l[ii].end()); // this does not have to be here most probably - there should be no duplicate DOFs 
	  neighBorderDofs_l[ii].resize( std::distance(neighBorderDofs_l[ii].begin(), it_l) );
  }

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 13                                                                 "; system("date +%T.%6N"); }
#endif



  vector < vector < vector < int > > > myLambdas_B0; 
  for      (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  myLambdas_B0.push_back(vector < vector < int > > () );
  }

#ifdef HFETI_SOLVER
  tbb::concurrent_vector < tbb::concurrent_vector < tbb::concurrent_vector < int > > > myLambdas_l;
  for      (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  myLambdas_l.push_back(tbb::concurrent_vector < tbb::concurrent_vector < int > > () );
  }
#else
  vector < vector < vector < int > > > myLambdas_l;
  for      (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  myLambdas_l.push_back(vector < vector < int > > () );
  }
#endif

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 14                                                                 "; system("date +%T.%6N"); }
#endif
  cilk_for (int domain_index = 0; domain_index < data.size(); domain_index++) {

	  int neighDomNum_l = data.size() - 1;  // number of subdomains in the cluster = number of neighboring domains 
	  // -1 - means = minus the current domain 
	  vector < int > myNeighSubDomains_l; // my neighboring domains 
	  for (int ii = 0; ii < data.size(); ii++) {
		  if (ii != domain_index) {
			  myNeighSubDomains_l.push_back(ii); // !! POZOR - toto je fiktivni lokalni cislovani 
		  }
	  }
	  std::sort (myNeighSubDomains_l.begin(), myNeighSubDomains_l.end());

	  vector < longInt > myBorderDOFs_l;                // my border DOFs are here 
	  myBorderDOFs_l = neighBorderDofs_l[domain_index]; 



	  vector < vector < int > > myNeighsSparse_l; 
	  myNeighsSparse_l.resize( myBorderDOFs_l.size() );

	  for (unsigned int j = 0; j < myBorderDOFs_l.size(); j++) { 			
		  int t = myBorderDOFs_l[j];  
		  myNeighsSparse_l[j].reserve(5);
		  myNeighsSparse_l[j].resize(1); 
		  myNeighsSparse_l[j][0] = t;
	  }

	  for (unsigned int i = 0; i < myNeighSubDomains_l.size(); i++) {
		  int d = myNeighSubDomains_l[i];
		  int mydofs_index = 0; 
		  int nedofs_index = 0; 

		  //if ( i == 0 && MPIrank < myNeighSubDomains[i] ) {
		  if ( i == 0 && domain_index < myNeighSubDomains_l[i] ) {
			  for (unsigned int j = 0; j < myBorderDOFs_l.size(); j++)
				  myNeighsSparse_l[j].push_back(domain_index); //myNeighsSparse_l[j].push_back(MPIrank); 
		  } 

		  do {	

			  if ( neighBorderDofs_l[myNeighSubDomains_l[i]][nedofs_index] == myBorderDOFs_l[mydofs_index] ) { 
				  myNeighsSparse_l[mydofs_index].push_back(d);
				  mydofs_index++;
				  nedofs_index++;
			  } else {
				  if ( neighBorderDofs_l[myNeighSubDomains_l[i]][nedofs_index] > myBorderDOFs_l[mydofs_index] ) {
					  mydofs_index++;
				  } else {
					  nedofs_index++;
				  }
			  }

		  } while (mydofs_index < myBorderDOFs_l.size() && nedofs_index < neighBorderDofs_l[myNeighSubDomains_l[i]].size() );
		  //} while (mydofs_index < myBorderDOFs.size() && nedofs_index < neighBorderDofs[i].size() );


		  if ( i < myNeighSubDomains_l.size() - 1) {
			  if ( domain_index > myNeighSubDomains_l[i] && domain_index < myNeighSubDomains_l[i+1] ) // if ( MPIrank > myNeighSubDomains[i] && MPIrank < myNeighSubDomains[i+1] )
				  for (unsigned int j = 0; j < myBorderDOFs_l.size(); j++)
					  myNeighsSparse_l[j].push_back(domain_index); // myNeighsSparse_l[j].push_back(MPIrank); 
		  }


		  if ( i == myNeighSubDomains_l.size()-1 && domain_index > myNeighSubDomains_l[i] ) { // if ( i == myNeighSubDomains.size()-1 && MPIrank > myNeighSubDomains[i] ) {
			  for (unsigned int j = 0; j < myBorderDOFs_l.size(); j++)
				  myNeighsSparse_l[j].push_back(domain_index); 
		  } 


	  }

	  // temp variable(s) for sort and unique functions 

	  for (unsigned int i = 0; i < myNeighsSparse_l.size(); i++) {
		  //it = std::unique (myNeighsSparse[i].begin() + 1, myNeighsSparse[i].end());
		  //myNeighsSparse[i].resize( std::distance(myNeighsSparse[i].begin() + 1, it) );

		  if (myNeighsSparse_l[i].size() < 3) {
			  myNeighsSparse_l[i].clear(); 
		  } /*else {
			std::sort (myNeighsSparse[i].begin() + 1, myNeighsSparse[i].end());
			}*/
		  //myNeighsSparse.erase(myNeighsSparse.begin() + i); //  clear();
	  }

	  //for (int i = 0; i < myNeighs.size(); i++) {
	  //	std::sort (myNeighs[i].begin(), myNeighs[i].end());
	  //	it = std::unique (myNeighs[i].begin(), myNeighs[i].end());
	  //	myNeighs[i].resize( std::distance(myNeighs[i].begin(), it) );
	  //	if (myNeighs[i].size() == 1)
	  //		myNeighs[i].clear();
	  //}	


	  //vector < vector < int > > myLambdas_l;
	  //myLambdas.resize( myNeighsSparse.size(), vector<int>( 5 , 0 ) );

	  //myLambdas_l.push_back(vector < vector < int > > () );

	  int lambdaCount_l   = 0; 
	  
	  int lambda_B0_count = 0; 
	  int DP_index        = 0; 
	  vector < longInt > DP_DOFs; 
	  DP_DOFs = fem[domain_index]->mesh.DP_DOFsAll; 
	  std::sort (DP_DOFs.begin(), DP_DOFs.end()); 
	  	  
	  for (int i = 0; i < myNeighsSparse_l.size(); i++) {
		  int add    = 0;
		  
		  if (myNeighsSparse_l[i].size() > 0) {
			  int dof = myNeighsSparse_l[i][0];

			  int cnt_tmp = myNeighsSparse_l[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb
			  //for (int t = myNeighsSparse[i].size() - 2; t > 0; t--) 
			  //	cnt_tmp += t; 

			  for (unsigned int j = 1; j < myNeighsSparse_l[i].size(); j++) {
				  int neighSD = myNeighsSparse_l[i][j];
				  //domainG.flag_redund_lagr_mult

				  // BEG - For B0 local 
				  if (DP_DOFs[DP_index] == dof)  
				  {
				  	  
					  if (add == 1) {
						  vector <int> tmp; 
						  tmp.reserve  (6); 
						  tmp.push_back(lambda_B0_count); // at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter  
						  tmp.push_back(dof);			  // dof in global numbering 
						  tmp.push_back(domain_index);    // MPIrank  // my sub-domainG 
						  tmp.push_back(neighSD);		  // neigh. sub-domainG 
						  tmp.push_back(cnt_tmp);		  // delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
						  myLambdas_B0[domain_index].push_back( tmp );
						  lambda_B0_count++; 

						  //if (flag_redund_lagr_mult) break;
					  } 
				  }
				  // END - For B0 local 
				  //else 
				  //{ // co jde do B0 nejde do B1
				  // BEG - For B1 local 
					  if (add == 1) {
#ifdef HFETI_SOLVER
						  tbb::concurrent_vector <int> tmp;
#else
						  vector <int> tmp; 
#endif
						  tmp.reserve  (6); 
						  tmp.push_back(lambdaCount_l);   // at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter  
						  tmp.push_back(dof);			  // dof in global numbering 
						  tmp.push_back(domain_index);    // MPIrank  // my sub-domainG 
						  tmp.push_back(neighSD);		  // neigh. sub-domainG 
						  tmp.push_back(cnt_tmp);		  // delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
						  myLambdas_l[domain_index].push_back( tmp );
						  lambdaCount_l++; 

						  if (flag_redund_lagr_mult) break;
					  }
				  //}
				  
				  if (neighSD == domain_index) // if (neighSD == MPIrank)
					  add = 1; 
				  // END - For B1 local 
			  }

			  if (DP_DOFs[DP_index] == dof) 
				  DP_index++; 

		  }
	  }


  }

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 15                                                                 "; system("date +%T.%6N"); }
#endif
  

  //if (MPIrank == 0) cout << "BEG - Create lambda global numbering " << endl; 
  // Create lambda global numbering
  int lambdaLocalCount_l = 0;
  for (int domain_index=0; domain_index < data.size(); domain_index++) {
	for (int i = 0; i < myLambdas_l[domain_index].size(); i++) {
		myLambdas_l[domain_index][i][0] = lambdaLocalCount_l;															
		lambdaLocalCount_l++;
	}
  }

  //B0
  int lambdaLocalCount_B0 = 0;
  for (int domain_index=0; domain_index < data.size(); domain_index++) {
	  for (int i = 0; i < myLambdas_B0[domain_index].size(); i++) {
		  myLambdas_B0[domain_index][i][0] = lambdaLocalCount_B0;															
		  lambdaLocalCount_B0++;
	  }
  }

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 16                                                                 "; system("date +%T.%6N"); }
#endif

  // Renumbering of the local B1 matrix - it is placed after dirichlet and B1 from cluster 
  // Create lambda global counting 
  int lambdaGlobalCount_l;   
  MPI_Exscan(&lambdaLocalCount_l, &lambdaGlobalCount_l, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
  if (MPIrank == 0) lambdaGlobalCount_l = 0; 
  int lambdaNum_l = lambdaGlobalCount_l + lambdaLocalCount_l;  

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 17                                                                 "; system("date +%T.%6N"); }
#endif

  MPI_Bcast(&lambdaNum_l,1,MPI_INT,MPIsize-1,MPI_COMM_WORLD);
  int total_number_of_local_B1_lambdas = lambdaNum_l;
  
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 18                                                                 "; system("date +%T.%6N"); }
#endif

  cilk_for (int domain_index=0; domain_index < data.size(); domain_index++) {
	  for (int i = 0; i < myLambdas_l[domain_index].size(); i++) {
		  myLambdas_l[domain_index][i][0] += total_number_of_dirichlet_lambdas + 
											 total_number_of_global_B1_lambdas + 
											 lambdaGlobalCount_l;															
		  
	  }
  }
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 19                                                                 "; system("date +%T.%6N"); }
#endif
  // decode neighbor's lambdas 
  vector < int > lam_sizes ( data.size() );
  int lam_max_size = 0; 
  for (int i = 0; i < data.size(); i++) {
	  lam_sizes[i] = myLambdas_l[i].size();
	  if (lam_sizes[i] > lam_max_size) lam_max_size = lam_sizes[i] ; 
  }

  lam_max_size = lam_max_size * 2; 
  for (int i = 0; i < data.size(); i++) {
	  myLambdas_l[i].reserve(lam_max_size); 
  }
   
  cilk_for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  for (int neigh_domain_index = 0; neigh_domain_index < data.size(); neigh_domain_index++) {
		  if (domain_index != neigh_domain_index) {
			  //for (int neighLambda_index = 0; neighLambda_index < myLambdas_l[neigh_domain_index].size(); neighLambda_index++) {
			  for (  int neighLambda_index = 0; neighLambda_index < lam_sizes[neigh_domain_index];          neighLambda_index++) {
				  if (myLambdas_l[neigh_domain_index][neighLambda_index][3] == domain_index) {
#ifdef HFETI_SOLVER
					  tbb::concurrent_vector <int> tmp;
#else
					  vector <int> tmp; 
#endif					  
					  tmp.reserve(6); 
					  tmp.push_back(myLambdas_l[neigh_domain_index][neighLambda_index][0]); // lambda global number  
					  tmp.push_back(myLambdas_l[neigh_domain_index][neighLambda_index][1]);	// dof in global numbering 
					  tmp.push_back(myLambdas_l[neigh_domain_index][neighLambda_index][3]); // my sub-domainG 
					  tmp.push_back(myLambdas_l[neigh_domain_index][neighLambda_index][2]);	// negh. sub-domainG 
					  tmp.push_back(myLambdas_l[neigh_domain_index][neighLambda_index][4]);	// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
					  
					  myLambdas_l[domain_index].push_back( tmp );
					   
				  }
			  }
		  }
	  }
  }
   
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 20                                                                 "; system("date +%T.%6N"); }
#endif

  //B0
  for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  for (int neigh_domain_index = 0; neigh_domain_index < data.size(); neigh_domain_index++) {
		  if (domain_index != neigh_domain_index) {
			  for (int neighLambda_index = 0; neighLambda_index < myLambdas_B0[neigh_domain_index].size(); neighLambda_index++) {
				  if (myLambdas_B0[neigh_domain_index][neighLambda_index][3] == domain_index) {
					  vector <int> tmp; 
					  tmp.reserve(6); 
					  tmp.push_back(myLambdas_B0[neigh_domain_index][neighLambda_index][0]);    // lambda global number  
					  tmp.push_back(myLambdas_B0[neigh_domain_index][neighLambda_index][1]);	// dof in global numbering 
					  tmp.push_back(myLambdas_B0[neigh_domain_index][neighLambda_index][3]);    // my sub-domainG 
					  tmp.push_back(myLambdas_B0[neigh_domain_index][neighLambda_index][2]);	// negh. sub-domainG 
					  tmp.push_back(myLambdas_B0[neigh_domain_index][neighLambda_index][4]);	// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
					  myLambdas_B0[domain_index].push_back( tmp );
				  }
			  }
		  }
	  }
  }

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 21                                                                 "; system("date +%T.%6N"); }
#endif
  cilk_for (int domain_index = 0; domain_index < data.size(); domain_index++) {

	  int lambda_l; 
	  int DOFNumber_l; 
	  int n_myLambdas_l = myLambdas_l[domain_index].size();
	  int cnt_l; 

	  for (int i = 0; i < n_myLambdas_l; i++) 
	  {
		  lambda_l	  = myLambdas_l[domain_index][i][0];
		  //DOFNumber = fem->mesh.g2l[myLambdas[i][1]]; // maping g2l (global to local)
		  DOFNumber_l = myLambdas_l[domain_index][i][1];
		  cnt_l		  = myLambdas_l[domain_index][i][4];

		  data[domain_index]->B->BI_l.push_back(lambda_l); 
		  
		  //data[domain_index]->B->BJ_l.push_back(DOFNumber_l);
		  data[domain_index]->B->BJ_l.push_back( fem[domain_index]->mesh.g2l[DOFNumber_l] ); 

		  //Weight[DOFNumber] = cnt_l;

		  if (myLambdas_l[domain_index][i][2] < myLambdas_l[domain_index][i][3])
			  data[domain_index]->B->BV_l.push_back( 1.0); 
		  else 
			  data[domain_index]->B->BV_l.push_back(-1.0); 
	  }

	  if (n_myLambdas_l == 0) {
		  data[domain_index]->B->BI_l.push_back(0);
		  data[domain_index]->B->BJ_l.push_back(0);
		  data[domain_index]->B->BV_l.push_back(0);
	  }

  }

#ifndef WIN32
  if (MPIrank == 0) { cout << " point 22                                                                 "; system("date +%T.%6N"); }
#endif  
  //BO

  cilk_for (int domain_index = 0; domain_index < data.size(); domain_index++) {

	  int lambda_B0; 
	  int DOFNumber_B0; 
	  int n_myLambdas_B0 = myLambdas_B0[domain_index].size();
	  int cnt_B0; 

	  for (int i = 0; i < n_myLambdas_B0; i++) 
	  {
		  lambda_B0	   = myLambdas_B0[domain_index][i][0];
		  DOFNumber_B0 = myLambdas_B0[domain_index][i][1];
		  cnt_B0	   = myLambdas_B0[domain_index][i][4];

		  data[domain_index]->B->B0_I.push_back( lambda_B0 ); 
		  data[domain_index]->B->B0_J.push_back( fem[domain_index]->mesh.g2l[DOFNumber_B0] );  // maping g2l (global to local)

		  //Weight[DOFNumber] = cnt_l;
		  if (myLambdas_B0[domain_index][i][2] < myLambdas_B0[domain_index][i][3])
			  data[domain_index]->B->B0_V.push_back( 1.0); 
		  else 
			  data[domain_index]->B->B0_V.push_back(-1.0); 
	  }

	  data[domain_index]->B->B0_cols = data[domain_index]->KSparse->n_row;  
	  data[domain_index]->B->B0_rows = lambdaLocalCount_B0;
	  data[domain_index]->B->B0_nnz  = n_myLambdas_B0;


	  if (n_myLambdas_B0 == 0) {
		  data[domain_index]->B->B0_I.push_back(0);
		  data[domain_index]->B->B0_J.push_back(0);
		  data[domain_index]->B->B0_V.push_back(0);
	  }
  }
	

#ifndef WIN32
  if (MPIrank == 0) { cout << "END - Creating Local  B1                                                  "; system("date +%T.%6N"); }
#endif
  // output variables for sparse matrix B 
  // data[domain_index]->B->BI_l
  // data[domain_index]->B->BJ_l
  // data[domain_index]->B->BV_l
  

  // *** Local B1 - end ********************************************************************

#ifndef WIN32
  if (MPIrank == 0) { cout << endl  << "BEG - Putting B1 together                                                 "; system("date +%T.%6N"); }
#endif

  for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  data[domain_index]->B->BI_full.reserve(data[domain_index]->B->BI_dir.size() +  data[domain_index]->B->BI.size() +     data[domain_index]->B->BI_l.size() );
	  data[domain_index]->B->BI_full.insert( data[domain_index]->B->BI_full.begin(), data[domain_index]->B->BI_dir.begin(), data[domain_index]->B->BI_dir.end() ); 
	  data[domain_index]->B->BI_full.insert( data[domain_index]->B->BI_full.end()  , data[domain_index]->B->BI    .begin(), data[domain_index]->B->BI    .end() ); 
	  data[domain_index]->B->BI_full.insert( data[domain_index]->B->BI_full.end()  , data[domain_index]->B->BI_l  .begin(), data[domain_index]->B->BI_l  .end() ); 

	  data[domain_index]->B->BJ_full.reserve(data[domain_index]->B->BJ_dir.size() +  data[domain_index]->B->BJ.size() +     data[domain_index]->B->BJ_l.size() );
	  data[domain_index]->B->BJ_full.insert( data[domain_index]->B->BJ_full.begin(), data[domain_index]->B->BJ_dir.begin(), data[domain_index]->B->BJ_dir.end() ); 
	  data[domain_index]->B->BJ_full.insert( data[domain_index]->B->BJ_full.end()  , data[domain_index]->B->BJ    .begin(), data[domain_index]->B->BJ    .end() ); 
	  data[domain_index]->B->BJ_full.insert( data[domain_index]->B->BJ_full.end()  , data[domain_index]->B->BJ_l  .begin(), data[domain_index]->B->BJ_l  .end() ); 

	  data[domain_index]->B->BV_full.reserve(data[domain_index]->B->BV_dir.size() +  data[domain_index]->B->BV.size() +     data[domain_index]->B->BV_l.size() );
	  data[domain_index]->B->BV_full.insert( data[domain_index]->B->BV_full.begin(), data[domain_index]->B->BV_dir.begin(), data[domain_index]->B->BV_dir.end() ); 
	  data[domain_index]->B->BV_full.insert( data[domain_index]->B->BV_full.end()  , data[domain_index]->B->BV    .begin(), data[domain_index]->B->BV    .end() ); 
	  data[domain_index]->B->BV_full.insert( data[domain_index]->B->BV_full.end()  , data[domain_index]->B->BV_l  .begin(), data[domain_index]->B->BV_l  .end() ); 

	  data[domain_index]->B->B_full_rows = total_number_of_dirichlet_lambdas + 
										   total_number_of_global_B1_lambdas + 
										   total_number_of_local_B1_lambdas;
	  data[domain_index]->B->B_full_cols = data[domain_index]->KSparse->n_row; 

	  data[domain_index]->B->B_full_nnz	 = data[domain_index]->B->BV_full.size(); 

  }

   

  if (MPIrank == 0)
	  cout << "Creating lambda_map_sub vector of vectors - Dirichlet" << endl; 
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 23                                                                 "; system("date +%T.%6N"); }
#endif

  // for Dirichlet 
  for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  for (int i = 0; i < data[domain_index]->B->BI_dir.size(); i++) {
		  vector <int> tmp_vec (2,0); 
		  tmp_vec[0] = data[domain_index]->B->BI_dir[i]; 
		  tmp_vec[1] = MPIrank; 
		  
		  //for (int domain_index2 = 0; domain_index2 < data.size(); domain_index2++) 
		  //	  data[domain_index2]->B->lambda_map_sub.push_back(tmp_vec); 

		  //data[0]->B->lambda_map_sub.push_back(tmp_vec);
		  domainG->lambda_map_sub.push_back(tmp_vec);


		  //tmp_vec[0] = data[domain_index]->B->BI_dir[i]; 
		  //tmp_vec[1] = domain_index; 
		  //tmp_vec[2] = domainG->vec_globalSubNumbering[domain_index];
		  data[domain_index]->B->lambda_map_sub.push_back( data[domain_index]->B->BI_dir[i] ); 



	  }
  }

  if (MPIrank == 0)
	  cout << "Creating lambda_map_sub vector of vectors - Global B1" << endl; 
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 24                                                                 "; system("date +%T.%6N"); }
#endif

  // for global B1 
  for (int i = 0; i < myLambdas.size(); i++) {
	  vector <int> tmp_vec (3,0); 
	  tmp_vec[0] = myLambdas[i][0]; 
	  tmp_vec[1] = myLambdas[i][2]; 
	  tmp_vec[2] = myLambdas[i][3];
	  //data[0]->B->lambda_map_sub.push_back(tmp_vec);
	  domainG->lambda_map_sub.push_back(tmp_vec);

	  //tmp_vec.resize(2);
	  //tmp_vec[1] = myLambdas[i][5]; 
	  data[myLambdas[i][5]]->B->lambda_map_sub.push_back( myLambdas[i][0] );

  }

  if (MPIrank == 0)
	  cout << "Creating lambda_map_sub vector of vectors - Local B1" << endl; 
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 25                                                                 "; system("date +%T.%6N"); }
#endif

  //for local B1 
  for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	  for (int i = 0; i < myLambdas_l[domain_index].size(); i++) {
		  if ( myLambdas_l[domain_index][i][2] < myLambdas_l[domain_index][i][3]) {
			  vector <int> tmp_vec (2,0); 
			  tmp_vec[0] = myLambdas_l[domain_index][i][0]; 
			  tmp_vec[1] = MPIrank; //myLambdas_l[domain_index][i][2]; 
			  //tmp_vec[2] = myLambdas_l[domain_index][i][3];
			  
			  //for (int domain_index2 = 0; domain_index2 < data.size(); domain_index2++) 
			  //	  data[domain_index2]->B->lambda_map_sub.push_back(tmp_vec); 
			  
			  //data[0]->B->lambda_map_sub.push_back(tmp_vec);
			  domainG->lambda_map_sub.push_back(tmp_vec);

			  //data[domain_index]->B->lambda_map_sub.push_back( myLambdas_l[myLambdas_l[domain_index][i][2]][i][0] ); 
			  //data[domain_index]->B->lambda_map_sub.push_back( myLambdas_l[myLambdas_l[domain_index][i][3]][i][0] ); 

		  }

		  data[domain_index]->B->lambda_map_sub.push_back( myLambdas_l[domain_index][i][0] ); 
	  }
	  std::sort( data[domain_index]->B->lambda_map_sub.begin(), data[domain_index]->B->lambda_map_sub.end() );
	
  }
  
  

  if (MPIrank == 0)
	  cout << "Creating lambda_map_sub vector of vectors - END" << endl; 
#ifndef WIN32
  if (MPIrank == 0) { cout << " point 26                                                                 "; system("date +%T.%6N"); }
#endif

  for (int domain_index = 0; domain_index < data.size(); domain_index++) {
	data[domain_index]->B->neigh_clusters = myNeighSubDomains; // in case of Hybrid FETI or more domains per cluster in general, this is neighboring clusters 
  }

	LogEventEnd(event1);
} 
