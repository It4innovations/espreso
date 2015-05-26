
#include "BSparse.h"

static int event1=0;

CBSparse::CBSparse(MPI_Comm comm) 
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
  n_col_bg = 0;
	nnz_bg   = 0;
	n_row_bd = 0;
	n_col_bd = 0;
	nnz_bd   = 0;
	n_row_bc = 0;
	nnz_bc   = 0;
	n_row_eq = 0;
	n_col    = 0;
  n_bd     = 0;

	Bi_coo = NULL; 
	Bj_coo = NULL;
	Bv_coo = NULL;
}

CBSparse::~CBSparse() {
  //cout << "BSparse destructor " <<endl;
  /* meaningless ??? */
	BI.clear();  BI.push_back(0);
	BJ.clear();  BJ.push_back(0);
	BV.clear();  BV.push_back(0);

	BI_l.clear();  BI_l.push_back(0);
	BJ_l.clear();  BJ_l.push_back(0);
	BV_l.clear();  BV_l.push_back(0);

	BI_dir.clear();  BI_dir.push_back(0);
	BJ_dir.clear();  BJ_dir.push_back(0);
	BV_dir.clear();  BV_dir.push_back(0);

	BI_full.clear();  BI_full.push_back(0);
	BJ_full.clear();  BJ_full.push_back(0);
	BV_full.clear();  BV_full.push_back(0);

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

int CBSparse::coo2crs(int m, int nnz, int *row_ind, int *row_ptr) {
	int i,row;

	for (i=0,row=0; row<=m; row++) {
		row_ptr[row]=i;
		while (i < nnz && row == row_ind[i]) i++;
	}
	return 0;
}

void CBSparse::createB(CFem *fem,CDomain *domainG, int nPrimalDOFs,int *borderDofs) {
	if (!event1) LogEventRegister("createB",&event1);
	LogEventBegin(event1);

	int *Weight = new int [domainG->neqSub[fem->i_domOnClust]];

	for (int i = 0;i<domainG->neqSub[fem->i_domOnClust];i++){
		Weight[i] = 1;
	}

	int cnt;
	int MPIrank, MPIsize;
	MPI_Comm_rank(comm, &MPIrank);
	MPI_Comm_size(comm, &MPIsize);

	int neighDomNum = domainG->n_neighbSub;      // number of neighboring sub-domains for current sub-domainG 
	int borderDofNum = domainG->n_exterDOFs[fem->i_domOnClust];       // number of DOFs on the border for current sub-domainG



	vector < vector < int > > neighBorderDofs;        // 2D vector for DOFs on borders of neighboring sub-domains 
	neighBorderDofs.resize( neighDomNum, vector<int>( 0 , 0 ) );

	vector < int > myBorderDOFs;  // my border DOFs are here 
	myBorderDOFs.insert(myBorderDOFs.begin(),
            borderDofs, borderDofs + domainG->n_exterDOFs[fem->i_domOnClust]); 
	std::sort (myBorderDOFs.begin(), myBorderDOFs.end());

	vector < int > myNeighSubDomains; // my neighboring domains 
	myNeighSubDomains.insert(myNeighSubDomains.begin(),	fem->mesh.neighbSub, fem->mesh.neighbSub + domainG->n_neighbSub); 
	std::sort (myNeighSubDomains.begin(), myNeighSubDomains.end());

	neigh_domains = myNeighSubDomains;

	MPI_Request mpi_req;
	MPI_Status mpi_stat;

	MPI_Request * mpi_send_req  = new MPI_Request [neighDomNum];
	MPI_Request * mpi_recv_req  = new MPI_Request [neighDomNum];

	MPI_Status  * mpi_recv_stat = new MPI_Status  [neighDomNum];

	

	// sending my DOFs on border to all neighboring sub-domains 
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		MPI_Isend(&myBorderDOFs[0], myBorderDOFs.size(), MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]); 
	}

	//int * mpi_buff; 
	//mpi_buff = (int*)malloc(borderDofNum * sizeof(int));
	// receiving all border DOFs from all neighboring sub-domains 
	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		//MPI_Recv(mpi_buff,myBorderDOFs.size(),MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD,&mpi_stat);
		//neighBorderDofs[i].insert(neighBorderDofs[i].begin(), mpi_buff, mpi_buff + myBorderDOFs.size());
		neighBorderDofs[i].resize(myBorderDOFs.size()); 
		MPI_Irecv(&neighBorderDofs[i][0], myBorderDOFs.size(),MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_recv_req[i]);
		
	}

// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------
	MPI_Waitall(neighDomNum, mpi_recv_req, mpi_recv_stat); 

  delete [] mpi_recv_req;
  delete [] mpi_recv_stat;
	//free(mpi_buff);

	// NOW : 
	// my neighboring sub-domains are in : myNeighSubDomains
	// my border DOFs are in : myBorderDOFs
	// neighboring sub-domains border DOFs are in neighBorderDofs[i] 

 	vector < vector < int > > myNeighsSparse; 
	myNeighsSparse.resize( myBorderDOFs.size() ); //, vector<int>( 2 , 0 ) );

	//int d = MPIrank;
	for (unsigned int j = 0; j < myBorderDOFs.size(); j++) { 			
		int t = myBorderDOFs[j];  // !!!!!!!!!!!!! - korekce na cislovani subdomen od 0 namisto od 1
		myNeighsSparse[j].reserve(5);
		myNeighsSparse[j].resize(1); 
		myNeighsSparse[j][0] = t;
		//myNeighsSparse[j][1] = d; 
		//vector <int> tmp; //tmp.push_back(t); //tmp.push_back(d); //myNeighsSparse.push_back( tmp );
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
	std::vector<int>::iterator it;

	for (unsigned int i = 0; i < myNeighsSparse.size(); i++) {
		//it = std::unique (myNeighsSparse[i].begin() + 1, myNeighsSparse[i].end());
		//myNeighsSparse[i].resize( std::distance(myNeighsSparse[i].begin() + 1, it) );

		if (myNeighsSparse[i].size() < 3) {
			myNeighsSparse[i].clear(); 
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


	vector < vector < int > > myLambdas;
	//myLambdas.resize( myNeighsSparse.size(), vector<int>( 5 , 0 ) );

	int lambdaCount = 0; 
	for (unsigned int i = 0; i < myNeighsSparse.size(); i++) {
		int add = 0;
		if (myNeighsSparse[i].size() > 0) {
			int dof = myNeighsSparse[i][0];

			int cnt_tmp = myNeighsSparse[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb
			//for (int t = myNeighsSparse[i].size() - 2; t > 0; t--) 
			//	cnt_tmp += t; 

			for (unsigned int j = 1; j < myNeighsSparse[i].size(); j++) {
				int neighSD = myNeighsSparse[i][j];
//        domainG.flag_redund_lagr_mult

				if (add == 1) {
					vector <int> tmp; 
					tmp.reserve(6); 
					tmp.push_back(lambdaCount);	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter  
					tmp.push_back(dof);			// dof in global numbering 
					tmp.push_back(MPIrank);		// my sub-domainG 
					tmp.push_back(neighSD);		// negh. sub-domainG 
					tmp.push_back(cnt_tmp);		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb 
					myLambdas.push_back( tmp );
					lambdaCount++; 
          if (!domainG->flag_redund_lagr_mult) break;
				}

				if (neighSD == MPIrank)
					add = 1; 
			}
		}
	}

	// Create lambda global counting 
	int lambdaGlobalCount = 0; 
	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
	if (MPIrank == 0) lambdaGlobalCount = 0; 
	int lambdaNum = lambdaCount + lambdaGlobalCount;  
	MPI_Bcast(&lambdaNum,1,MPI_INT,MPIsize-1,MPI_COMM_WORLD);

	if (myLambdas.size() > 0) 
		for (unsigned int i = 0; i < myLambdas.size(); i++) 
			myLambdas[i][0] = myLambdas[i][0] + lambdaGlobalCount; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index

	vector < vector < int > > mpi_send_buff;
	mpi_send_buff.resize( myNeighSubDomains.size(), vector<int>( 0 , 0 ) );

	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) {
		int index = 0;  

		if (myLambdas.size() > 0 )
		{
			for (unsigned int j = 0; j < myLambdas.size(); j++) {
				if(myLambdas[j][3] == myNeighSubDomains[i] ) { // POZOR - oprava pro pocitani pod oblasti od 0 misto od 1 
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

		//MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_req); // POZOR - oprava pro pocitani pod oblasti od 0 misto od 1
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), MPI_INT, myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
	}

	vector < vector < int > > mpi_recv_buff;
	mpi_recv_buff.resize( myNeighSubDomains.size(), vector<int>( 0 , 0 ) );

  delete [] mpi_send_req;
	//int *mpi_tmp_recv_buff; // = (int*)malloc(borderDofNum * sizeof(int));
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

	for (unsigned int i = 0; i < myNeighSubDomains.size(); i++) 
	{
		if (mpi_recv_buff[i].size() > 1)
		{
			for (unsigned int j = 0; j < mpi_recv_buff[i].size(); j=j+5) 
			{
				vector < int > tmp; 
				int ltmp = mpi_recv_buff[i][j + 0];  
				tmp.reserve(5);
				tmp.push_back(mpi_recv_buff[i][j + 0]);
				tmp.push_back(mpi_recv_buff[i][j + 1]);
				tmp.push_back(mpi_recv_buff[i][j + 3]);	// to change the sign in the lambda 
				tmp.push_back(mpi_recv_buff[i][j + 2]);	// to change the sign in the lambda 
				tmp.push_back(mpi_recv_buff[i][j + 4]);

				if (myLambdas.size() == 0) 
				{
					myLambdas.push_back(tmp);
				} 
				else 
				{			
					for (unsigned int s = 0; s < myLambdas.size(); s++) 
					{
						if (myLambdas[s][0] > ltmp) 
						{
							myLambdas.insert(myLambdas.begin() + s, tmp); 
							break;
						}
						if (s == myLambdas.size() - 1) {
							myLambdas.push_back(tmp);
							break;
						}	
					}
				}
			}
		}
	}



	// output variables for sparse matrix B 
	//vector < int > Bi; 
	//vector < int > Bj;
	//vector < double > Bv;


	int lambda; 
	int DOFNumber; 
	int n_myLambdas=myLambdas.size();
	for (int i = 0; i < n_myLambdas; i++) 
	{
		lambda = myLambdas[i][0];
		DOFNumber = fem->mesh.g2l[myLambdas[i][1]]; // maping g2l (global to local)
		cnt = myLambdas[i][4];

		BI.push_back(lambda); 
		BJ.push_back(DOFNumber);

		Weight[DOFNumber] = cnt;

		if (myLambdas[i][2] < myLambdas[i][3])
			BV.push_back(1.0); 
		else 
			BV.push_back(-1.0); 
	}

	if (n_myLambdas == 0) {
		BI.push_back(0);
		BJ.push_back(0);
		BV.push_back(0);
	}


	/* multiplicityPD function */
	multiplicityPD = &Weight[0];
	i_bg = &BI[0];
	j_bg = &BJ[0];
	v_bg = &BV[0];
	n_row_bg = lambdaNum;
	nnz_bg   = BI.size();

#if 0
	//  +++++ print out to ascii- file
	char filenameKe[128];
	sprintf(filenameKe, "data/Bf_%d.dat", MPIrank);
	FILE *fVTK = fopen(filenameKe, "w");
	for (int i = 0;i<BI.size();i++){
		fprintf(fVTK,"%d  %d  %3.3e\n",BI[i],BJ[i],BV[i]);
	}
	fclose(fVTK);
#endif


#ifdef DEBUG
	if (MPIrank==0) printf("define lambda and create B - end\n");
#endif


	// -------------------------------------------------------------------------
#ifdef DEBUG
	if (MPIrank==0) printf("Dirichlet BC - start\n");
#endif




	int lo,hi;

	int nLamDir = fem->bound_cond.dirBCSub->n;
	MPI_Status status;
	lo = 0;
	if (MPIrank>0){
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


	int *I_bd = new int[nLamDir];  
	int *J_bd = new int[nLamDir];  
	double *V_bd = new double[nLamDir];  

	int cnt1=0;
	for (int i = 0;i<nLamDir;i++){
		I_bd[cnt1]=cnt1+strtI;
		J_bd[i]=fem->bound_cond.dirBCSub->ind[i];
		V_bd[i]=1;
		cnt1++;
	}    
	i_bd = I_bd; 
	j_bd = J_bd; 
	v_bd = V_bd; 
	nnz_bd = nLamDir; 


	for (int i = 0; i < nLamDir; i++) {

		vector <int> tmp_vec; 

		//tmp_vec.push_back(n_row_bg + i_bd[i]);
		tmp_vec.push_back( i_bd[i] );
		tmp_vec.push_back( MPIrank ); 

		// POZOR
		//lambda_map_sub.push_back(tmp_vec); 
	}


	for (int i = 0; i < myLambdas.size(); i++) {

		vector < int > tmp_vec; 
		tmp_vec.push_back(n_row_bd + myLambdas[i][0]);
		tmp_vec.push_back(myLambdas[i][2]);
		tmp_vec.push_back(myLambdas[i][3]);

		//POZOR
		//lambda_map_sub.push_back(tmp_vec); 
	}



#ifdef DEBUG
	if (MPIrank==0) printf("Dirichlet BC - end\n");
#endif


#ifdef DEBUG
	if (MPIrank==0) printf("cont BC - start\n");
#endif

	if (domainG->flag_contact){
		int nLamCon= fem->bound_cond.conBCSub->n;
		MPI_Status status1;
		lo = 0;
		if (MPIrank>0){
			MPI_Recv(&lo, 1, MPI_INT, MPIrank-1, 124, MPI_COMM_WORLD, &status1);
		}
		hi = lo + nLamCon; 
		if (MPIrank < MPIsize-1) {
			MPI_Send(&hi, 1, MPI_INT, MPIrank+1, 124, MPI_COMM_WORLD);
		}
		strtI = lo;

		buf = lo + nLamCon;
		MPI_Bcast(&buf, 1, MPI_INT, MPIsize-1, MPI_COMM_WORLD);
		n_row_bc = buf;


		int *I_bc = new int[nLamCon];  
		int *J_bc = new int[nLamCon];  
		double *V_bc = new double[nLamCon];  

		cnt1=0;
		for (int i = 0;i<nLamCon;i++){
			I_bc[cnt1]=cnt1+strtI;
			J_bc[i]=fem->bound_cond.conBCSub->ind[i];
			V_bc[i]=fem->bound_cond.conBCSub->val[i];
			cnt1++;
		}    
		i_bc = I_bc; 
		j_bc = J_bc; 
		v_bc = V_bc; 
		nnz_bc = nLamCon; 
	}
	else{
		i_bc = NULL; 
		j_bc = NULL; 
		v_bc = NULL; 
		nnz_bc = 0; 
	}




#ifdef DEBUG
	if (MPIrank==0) printf("cont BC - end\n");
#endif
	///////////////////////////
#ifdef DEBUG
	if (MPIrank==0) printf("concatenating glueing and Dirichlet BC - start\n");
#endif
	//equality BC
	int Bd_m    = n_row_bd;
	int Bg_m    = n_row_bg;
	n_row_eq      = Bd_m + Bg_m;
	int Bd_nnz  = nnz_bd;
	int Bg_nnz  = nnz_bg;
	int B_nnz   = Bd_nnz + Bg_nnz;
#ifdef FLLOP_ENABLED
	int *ib     = new int[n_row_eq+1];
#endif
	int *Bi_coo = new int[B_nnz];
	int *jb     = new int[B_nnz];
	double *vb  = new double[B_nnz];

	memcpy(jb, j_bd, Bd_nnz*sizeof(int));
	memcpy(&jb[Bd_nnz], j_bg, Bg_nnz*sizeof(int));
	memcpy(vb, v_bd, Bd_nnz*sizeof(double));
	memcpy(&(vb[Bd_nnz]), v_bg, Bg_nnz*sizeof(double));
	memcpy(Bi_coo, i_bd, Bd_nnz*sizeof(int));
	memcpy(&(Bi_coo[Bd_nnz]), i_bg, Bg_nnz*sizeof(int));

	for (int i=Bd_nnz; i<B_nnz; i++) Bi_coo[i] = Bi_coo[i] + Bd_m;

#ifdef FLLOP_ENABLED
	coo2crs(n_row_eq, B_nnz, Bi_coo, ib);
#endif
	i_multiplicityDD = &Bi_coo[0];
	// inequality BC 
	if (domainG->flag_contact) {
		i_bc_crs = new int[n_row_bc+1];
		coo2crs(n_row_bc, nnz_bc, i_bc, i_bc_crs);
	}


	nnz_eq = n_myLambdas+nLamDir; //should be the same like 'nnz_eq = nnz_bg+nnz_bd'
#ifdef FLLOP_ENABLED
	int *WeightDual = new int [nnz_eq];
	for (int i = 0; i < nLamDir; i++) {
		WeightDual[i] = 1;
	}

	for (int i = 0; i < n_myLambdas; i++) {
		WeightDual[i+nLamDir] = int (myLambdas[i][4]);
	}
	/* multiplicityPD function */
	multiplicityDD = &WeightDual[0];
#endif
	//
#ifdef DEBUG
	if (MPIrank==0) printf("concatenating glueing and Dirichlet BC - end\n");
#endif

#ifdef FLLOP_ENABLED
	this->i_eq = ib;
	this->j_eq = jb;
	this->v_eq = vb;
	this->n_col = nPrimalDOFs;
#endif

	this->Bi_coo = Bi_coo;
	this->Bj_coo = jb;
	this->Bv_coo = vb;
	this->n_col = nPrimalDOFs;

	LogEventEnd(event1);
} 
