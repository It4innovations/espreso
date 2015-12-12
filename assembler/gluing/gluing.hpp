
#include "gluing.h"

namespace assembler {

template <class TInput>
Gluing<TInput>::Gluing(TInput &input): Assembler<TInput>(input) {
	_localB.resize(this->subdomains());
	_globalB.resize(this->subdomains());
	_B0.resize(this->subdomains());
	_B1.resize(this->subdomains());

	_lambda_map_sub_B1.resize(this->subdomains());
	_lambda_map_sub_B0.resize(this->subdomains());
	_B1_duplicity.resize(this->subdomains());
	_vec_c.resize(this->subdomains());
};

template <>
void Gluing<API>::computeSubdomainGluing()
{
	esint index_offset = _input.indexing;

//	std::map <esint, double> dirichlet;
//	for (esint i = 0; i < _input.dirichlet.size; i++)
//		dirichlet.insert(std::make_pair(_input.dirichlet.indices[i], _input.dirichlet.values[i] ));

	std::map <esint, esint> g2l;
	for (esint i = 0; i < _input.l2g.size; i++)
		g2l.insert(std::make_pair(_input.l2g.values[i],i));

	std::map<esint, esint> :: iterator map_it;

	std::vector<SparseDOKMatrix<eslocal> > gB(this->subdomains());
	std::vector<SparseDOKMatrix<eslocal> > lB(this->subdomains());

	eslocal lambda_count_B0 = 0, lambda_count_B1 = 0;

	//std::vector<eslocal> local_prim_numbering_d(this->subdomains(), 0);
	for (size_t i = 0; i < _input.dirichlet.size; i++) {

		esint loc_DOF_num = -1;
		map_it = g2l.find(_input.dirichlet.indices[i]);
		if (  map_it != g2l.end() ) {
			loc_DOF_num = map_it->second;

			gB[0](lambda_count_B1, loc_DOF_num) =  1.0;

			_lambda_map_sub_clst.push_back(std::vector<eslocal>({ lambda_count_B1, 0 }));
			_lambda_map_sub_B1[0].push_back(lambda_count_B1);
			_B1_duplicity[0].push_back(1.0);
			_vec_c[0].push_back(_input.dirichlet.values[i]);

			lambda_count_B1++;
		}
	}

	for (eslocal d = 0; d < this->subdomains(); d++) {
		gB[d].resize(lambda_count_B1, _input.l2g.size); //TODO: _K.rows or _K.cols
		_B1[d] = gB[d];

		lB[d].resize(lambda_count_B0, _input.l2g.size); //TODO: _K.rows or _K.cols
		_B0[d] = lB[d];
	}
}

template <>
void Gluing<API>::computeClusterGluing(std::vector<size_t> &rows)
{

	esint index_offset = _input.indexing;

	local_B1_global_resize();

	for (esint i = 0; i < _input.neighbourRanks.size; i++)
		_neighClusters.push_back(_input.neighbourRanks.values[i] - index_offset);

	_myBorderDOFs.resize(_input.l2g.size);
	for (esint i = 0; i < _input.l2g.size; i++)
		_myBorderDOFs[i] = _input.l2g.values[i];

	std::map <esint, esint> g2l;
	for (esint i = 0; i < _input.l2g.size; i++)
		g2l.insert(std::make_pair(_input.l2g.values[i],i));

	std::map<esint, esint> :: iterator map_it;

	int MPIrank = this->rank();
	int MPIsize = this->size();
	std::vector<SparseIJVMatrix<eslocal> > &B1 = _B1;
	std::vector<SparseIJVMatrix<eslocal> > &B0 = _B0;

	bool flag_redund_lagr_mult = true;

	std::vector < SparseDOKMatrix<eslocal> > B1_DOK_tmp(this->subdomains());
	MPI_Request   mpi_req;
	MPI_Status 	  mpi_stat;
	std::set<eslocal>::const_iterator it_set;
	std::set<eslocal>::const_iterator it_set_l;


	// Total number of neighboring clusters
	eslocal neighClustNum;	// number of neighboring sub-domains for current sub-domainG
	neighClustNum = _neighClusters.size();

    if (MPIrank == 0) { std::cout << " Global B - myNeighDOFs arrays are transfered to neighbors                "; system("date +%T.%6N"); }

    MPI_Request * mpi_send_req  = new MPI_Request [_neighClusters.size()];
	MPI_Request * mpi_recv_req  = new MPI_Request [_neighClusters.size()];
	MPI_Status  * mpi_recv_stat = new MPI_Status  [_neighClusters.size()];

	std::vector < std::vector < esglobal > > neighBorderDofs;   // 2D vector for DOFs on borders of neighboring sub-domains
	neighBorderDofs.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );

	// sending my border DOFs to all neighboring clusters
	for (eslocal i = 0; i < _neighClusters.size(); i++) {
		MPI_Isend(&_myBorderDOFs[0], _myBorderDOFs.size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
	}

	// receiving all border DOFs from all neighboring sub-domains
	eslocal messages_received_bd = 0;
	while ( messages_received_bd < _neighClusters.size() ) {
		for (eslocal i = 0; i < _neighClusters.size(); i++) {

			int flag;
			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

			if (flag) {
				int recv_msg_size = 0;

				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);

				neighBorderDofs[i].resize(recv_msg_size);
				if (recv_msg_size < 0)
					std::cout << "error in msg size !";
				MPI_Recv(&neighBorderDofs[i][0], recv_msg_size, esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);

				messages_received_bd++;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);


	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }


	// NOW :
	// my neighboring clusters are in :     	_neighClusters
	// my border DOFs are in : 					myBorderDOFs
	// neighboring clusters border DOFs are in:	neighBorderDofs[i]

 	std::vector < std::vector < esglobal > > myNeighsSparse;
	myNeighsSparse.resize( _myBorderDOFs.size() );

	for (eslocal j = 0; j < _myBorderDOFs.size(); j++) {
		myNeighsSparse[j].reserve(5);
		myNeighsSparse[j].push_back(_myBorderDOFs[j]);
	}

	for (eslocal i = 0; i < _neighClusters.size(); i++) {
		int d = _neighClusters[i];
		int mydofs_index = 0;
		int nedofs_index = 0;

		if ( i == 0 && MPIrank < _neighClusters[i] ) {
			for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank);
		}

		do {

			if ( neighBorderDofs[i][nedofs_index] == _myBorderDOFs[mydofs_index] ) {
				myNeighsSparse[mydofs_index].push_back(d);
				mydofs_index++;
				nedofs_index++;
			} else {
				if ( neighBorderDofs[i][nedofs_index] > _myBorderDOFs[mydofs_index] ) {
					mydofs_index++;
				} else {
					nedofs_index++;
				}
			}

		} while (mydofs_index < _myBorderDOFs.size() && nedofs_index < neighBorderDofs[i].size() );


		if ( i < _neighClusters.size() - 1)
			if ( MPIrank > _neighClusters[i] && MPIrank < _neighClusters[i+1] )
				for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
					myNeighsSparse[j].push_back(MPIrank);



		if ( i == _neighClusters.size() -1 && MPIrank > _neighClusters[i] )
			for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank);

	}

	for (eslocal i = 0; i < myNeighsSparse.size(); i++)
		if (myNeighsSparse[i].size() < 3)
			myNeighsSparse[i].clear();

	if (MPIrank == 0) { std::cout << " Global B - myNeighSparse assembled                                       "; system("date +%T.%6N"); }

	std::vector < std::vector < esglobal > > myLambdas;

	esglobal lambdaCount = 0;

	for (eslocal i = 0; i < myNeighsSparse.size(); i++) {
		bool add = false;
		if (myNeighsSparse[i].size() > 0) {
			esglobal dof = myNeighsSparse[i][0];

			eslocal cnt_tmp = myNeighsSparse[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb

			for (eslocal j = 1; j < myNeighsSparse[i].size(); j++) {
				esglobal neighSD = myNeighsSparse[i][j];

				if ( add ) {

					myLambdas.push_back ( std::vector < esglobal > () );
					myLambdas[lambdaCount].reserve(6);
					myLambdas[lambdaCount].resize(5);
					myLambdas[lambdaCount][0] = lambdaCount;			// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
					myLambdas[lambdaCount][1] = dof;					// dof in global numbering
					myLambdas[lambdaCount][2] = (esglobal)MPIrank;		// my sub-domainG
					myLambdas[lambdaCount][3] = neighSD;				// neigh. sub-domainG
					myLambdas[lambdaCount][4] = (esglobal)cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb
					lambdaCount++;

					if (!flag_redund_lagr_mult)
						break;
				}

				if (neighSD == (esglobal)MPIrank)
					add = true;
			}
		}
	}

	if (MPIrank == 0) { std::cout << " Global B - Create global lambda numbering                                "; system("date +%T.%6N"); }

	esglobal lambdaGlobalCount = 0;

	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);

	if (MPIrank == 0) lambdaGlobalCount = 0;
	esglobal lambdaNum = lambdaCount + lambdaGlobalCount;
	MPI_Bcast(&lambdaNum, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
	esglobal total_number_of_global_B1_lambdas = lambdaNum;

	if (myLambdas.size() > 0)
		//cilk_
		for (eslocal i = 0; i < myLambdas.size(); i++)
			myLambdas[i][0] = myLambdas[i][0]  + lambdaGlobalCount + total_number_of_B1_l_rows; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index


	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }


	std::vector < std::vector < esglobal > > mpi_send_buff;
	mpi_send_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
//#ifdef DEBUG
//	for (int i = 0; i < _neighClusters.size(); i++) {
//#else
	cilk_for (int i = 0; i < _neighClusters.size(); i++) {
//#endif
		int index = 0;
		if ( myLambdas.size() > 0 )
		{
			for (int j = 0; j < myLambdas.size(); j++) {
				if( myLambdas[j][3] == _neighClusters[i] ) {
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

	if (MPIrank == 0) { std::cout << " Global B - Isend                                                         "; system("date +%T.%6N"); }

	for (int i = 0; i < _neighClusters.size(); i++)
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);


	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }

	std::vector < std::vector < esglobal > > mpi_recv_buff;
	mpi_recv_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
	//delete [] mpi_send_req;


	// receiving all border DOFs from all neighboring sub-domains
	eslocal messages_received = 0;
	while ( messages_received < _neighClusters.size() ) {
		for (eslocal i = 0; i < _neighClusters.size(); i++) {

			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
			int flag;
			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

			if (flag) {
				int recv_msg_size = 0;

				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);

				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);

				free(mpi_tmp_recv_buff);
				messages_received++;
			}
		}
	}

	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }

	// decode received lambdas
	eslocal recv_lamba_count = 0;
	for (eslocal i = 0; i < _neighClusters.size(); i++)
		if (mpi_recv_buff[i].size() > 1)
			recv_lamba_count += mpi_recv_buff[i].size();

	recv_lamba_count = recv_lamba_count / 5;

	eslocal l_i = myLambdas.size();
	myLambdas.resize( myLambdas.size() + recv_lamba_count, std::vector< esglobal >( 5 , 0 ) );

	for (eslocal i = 0; i < _neighClusters.size(); i++) {
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

	if (MPIrank == 0) { std::cout << " Global B - myLambdas - sort or tbb:sort                                  "; system("date +%T.%6N"); }

	auto comp_vf = [](const std::vector<esglobal> &a, const std::vector<esglobal> &b) {return a[0] < b[0];};

//#ifdef USE_TBB
//	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
//#else
	std::sort         (myLambdas.begin(), myLambdas.end(), comp_vf);
//#endif

	if (MPIrank == 0) { std::cout << " Global B - Final B assembling with g2l mapping using std::map            "; system("date +%T.%6N"); }

	esglobal lambda;
	esglobal DOFNumber;
	eslocal n_myLambdas = myLambdas.size();

//#ifdef USE_TBB
//	tbb::mutex m;
//	cilk_for (int j = 0; j < n_myLambdas; j++)
//#else
	for (eslocal j = 0; j < n_myLambdas; j++)
//#endif
	{

		esglobal lambda         = myLambdas[j][0];
		esglobal DOFNumber      = myLambdas[j][1];

		//esglobal dofNODEnumber = DOFNumber / this->DOFs();
		//eslocal  dofNODEoffset = DOFNumber % this->DOFs();
		//eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );

		//const std::set < eslocal > & subs_with_element = this->_input.mesh.subdomainBoundaries()[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel

		eslocal  cnt            = myLambdas[j][4];
		double B_value;
		if (myLambdas[j][2] < myLambdas[j][3])
			B_value = 1.0;
		else
			B_value = -1.0;

		//TODO:Tuto smycku prepsat, aby fungovala jen podoblasti, ve ktery je uzel
		// Pomala a otestovana verze
		for (eslocal d = 0; d < this->subdomains(); d++) {

        //TODO: rychle verze nutno otestovat na realne uloze
		//for (it_set = localBoundaries[clustDofNODENumber].begin(); it_set != localBoundaries[clustDofNODENumber].end(); ++it_set) {
		//	eslocal d = *it_set;

			//eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
			//if ( domDofNODENumber != -1 ) {
				//vychazi z mapovani g2l pro uzly a ne DOFy

			map_it = g2l.find(DOFNumber);

				eslocal domDOFNumber = map_it->second; //this->DOFs() * domDofNODENumber + dofNODEoffset;
//			   #ifdef USE_TBB
//				m.lock();
//			   #endif
				B1_DOK_tmp[d](lambda, domDOFNumber) = B_value;
				myLambdas[j].push_back(d);
				_B1_duplicity[d].push_back(1.0 / cnt);
				_vec_c[d].push_back(0.0);
//			   #ifdef USE_TBB
//				m.unlock();
//			   #endif
			//	break;

			//}
		}
	}

	for (eslocal d = 0; d < this->subdomains(); d++){
		//TODO: lambdaNum muze byt 64bit integer
		//TODO: matice B - pocet radku muze byt 64bit int
		B1_DOK_tmp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas , rows[d]);
		SparseIJVMatrix<eslocal> ijv = B1_DOK_tmp[d];
		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
	}


	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }


	// for global B1
	for (int i = 0; i < myLambdas.size(); i++) {
		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
		tmp_vec[0] = myLambdas[i][0];
		tmp_vec[1] = myLambdas[i][2];
		tmp_vec[2] = myLambdas[i][3];

		_lambda_map_sub_clst.push_back(tmp_vec);

		eslocal lam_tmp = myLambdas [i][0]; //TODO: esglobal
		_lambda_map_sub_B1[ myLambdas[i][5] ].push_back( lam_tmp );
	}

	if (MPIrank == 0) { std::cout << " Dual size: " <<  total_number_of_B1_l_rows + total_number_of_global_B1_lambdas  << std::endl; }

        MPI_Barrier(MPI_COMM_WORLD);

   	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }



	for (int i = 0; i < _vec_c.size(); i++) {
		 _vec_c[i].resize(_B1_duplicity[i].size(), 0.0);
	}

	for (size_t p = 0; p < this->subdomains(); p++) {
		_localB[p] = _B0[p];
		_globalB[p] = _B1[p];
	}

}

template <class TInput>
void Gluing<TInput>::computeSubdomainGluing()
{

//#define PARALLEL
#define SEQUENTIAL


	const mesh::Boundaries &localBoundaries = this->_input.mesh.subdomainBoundaries();
	// TODO: cluster boundaries on surface
	const mesh::Boundaries &globalBoundaries = this->_input.mesh.clusterBoundaries();



	if (_B1.size() != this->subdomains())
		_B1.resize(this->subdomains());

	if (_B0.size() != this->subdomains())
		_B0.resize(this->subdomains());

	eslocal lambda_count_B0 = 0, lambda_count_B1 = 0;

	// TODO: make it more general
	std::vector<size_t> properties;
	properties.resize(this->DOFs());
	if (this->DOFs() == 1) {
		properties[0] = mesh::DIRICHLET_X;
	}
	if (this->DOFs() == 3) {
		properties[0] = mesh::DIRICHLET_X;
		properties[1] = mesh::DIRICHLET_Y;
		properties[2] = mesh::DIRICHLET_Z;
	}

#ifdef PARALLEL

	/* Numbers of processors, value of PAR_NUM_THREADS */
	int threads;
	char * var = getenv("PAR_NUM_THREADS");
    if(var != NULL)
    	sscanf( var, "%d", &threads );
	else {
    	printf("Parallel generator of B matrix - Please set environment PAR_NUM_THREADS to 1");
        exit(1);
	}

	std::vector < std::vector < std::vector < eslocal > > > mat_B_dir 				 (threads);
	std::vector < std::vector<eslocal> > 					local_prim_numbering_d_l (threads);
	std::vector<eslocal> 									lambda_count_B1_l 		 (threads, 0);
	std::vector<eslocal> 									lambda_count_B1_sum 	 (threads, 0);




	std::vector < std::vector < std::vector < double > > > 	vec_c_l ( threads );

#pragma omp parallel num_threads(threads)
	{
		int this_thread = omp_get_thread_num();
		int num_threads = omp_get_num_threads();


		std::vector<size_t>::const_iterator vi_l;
		std::set<eslocal>::const_iterator si1_l, si2_l;

		mat_B_dir[this_thread].			  	  resize( this->subdomains()    );
		local_prim_numbering_d_l[this_thread].resize( this->subdomains(), 0 );
		vec_c_l[this_thread].                 resize( this->subdomains()    );

		int chunk_size; 
                if (num_threads > 1) 
                    chunk_size = localBoundaries.size() / (num_threads-1);
                else 
                    chunk_size = localBoundaries.size(); 
                    
                for (int d = 0; d < this->subdomains(); d++) {
			mat_B_dir[this_thread][d].reserve(chunk_size);
			vec_c_l[this_thread][d].reserve(chunk_size);
		}

//#pragma omp for schedule (static, chunk_size)
//		for (size_t i = 0; i < localBoundaries.size(); i++) {

		size_t lo =  this_thread      * chunk_size;
		size_t hi = (this_thread + 1) * chunk_size;

		if (hi > localBoundaries.size())    hi = localBoundaries.size();
		if (this_thread == num_threads - 1) hi = localBoundaries.size();

		for (size_t i = lo; i < hi; i++) {

			std::vector<bool> is_dirichlet(this->DOFs(), false);
			for (si1_l = localBoundaries[i].begin(); si1_l != localBoundaries[i].end(); ++si1_l) {
				for (vi_l = properties.begin(); vi_l != properties.end(); ++vi_l) {
					const std::map<eslocal, double> &property = this->_input.mesh.coordinates().property(*vi_l).values();
					if (property.find(i) != property.end()) {

						double dirichlet_value = property.at(i);

						//is_dirichlet[0] = true;

						mat_B_dir[this_thread][*si1_l].push_back( lambda_count_B1_l[this_thread] );
						mat_B_dir[this_thread][*si1_l].push_back( local_prim_numbering_d_l[this_thread][*si1_l] + (vi_l - properties.begin()) );
						mat_B_dir[this_thread][*si1_l].push_back( 1 );

						vec_c_l   [this_thread][*si1_l].push_back( dirichlet_value );

						//gB[*si1](lambda_count_B1, local_prim_numbering_d[*si1] + (vi - properties.begin())) =  1.0;

						//_lambda_map_sub_clst.		push_back(std::vector<eslocal>({ lambda_count_B1, 0 }));
						//_lambda_map_sub_B1[*si1].	push_back(lambda_count_B1);
						//_B1_duplicity[*si1].		push_back(1.0);
						//_vec_c[*si1].				push_back(dirichlet_value);

						lambda_count_B1_l[this_thread]++;
					}
				}
			}

			for (si1_l = localBoundaries[i].begin(); si1_l != localBoundaries[i].end(); ++si1_l) {
				local_prim_numbering_d_l[this_thread][*si1_l] += this->DOFs();
			}


		}

#pragma omp barrier

#pragma omp master
		{
			lambda_count_B1_sum[0] = lambda_count_B1_l[0];
			for (size_t th = 1; th < threads; th++) {
				lambda_count_B1_sum[th] = lambda_count_B1_sum[th-1] + lambda_count_B1_l[th];
				for (size_t d = 0; d < this->subdomains(); d++) {
					local_prim_numbering_d_l[th][d] += local_prim_numbering_d_l[th - 1][d];
				}
			}
		}


#pragma omp for
		for (size_t d = 0; d < this->subdomains(); d++) {

			for (size_t th = 0; th < threads; th++) {
				eslocal lambda_off = 0;
				eslocal locDOF_off = 0;

				if (th > 0) lambda_off = lambda_count_B1_sum[th-1];

				for (size_t i = 0; i < mat_B_dir[th][d].size(); i = i + 3 ) {
					if (th > 0) locDOF_off = local_prim_numbering_d_l[th-1][d];

					eslocal lambda = mat_B_dir[th][d][i + 0] + lambda_off;
				    eslocal locDOF = mat_B_dir[th][d][i + 1] + locDOF_off;

				    //gB[d]( lambda , locDOF  ) = (double)mat_B_dir[th][d][i + 2];
				    _B1[d].set ( lambda , locDOF, (double)mat_B_dir[th][d][i + 2]);

					_lambda_map_sub_B1[d].	push_back(lambda);
					_B1_duplicity[d]. 		push_back(1.0);

				}

				_vec_c[d].insert( _vec_c[d].end(), vec_c_l[th][d].begin(), vec_c_l[th][d].end() );

			}
		}

#pragma omp master
		{
//			for (eslocal i = 0; i < lambda_count_B1_sum[threads - 1]; i++)
//				_lambda_map_sub_clst.		push_back(std::vector<eslocal>({ i, 0 }));
			lambda_count_B1 = lambda_count_B1_sum[threads - 1];
		}
	}
#endif //PARALLEL

#ifdef SEQUENTIAL

	std::vector<SparseDOKMatrix<eslocal> > gB(this->subdomains());
	std::vector<SparseDOKMatrix<eslocal> > lB(this->subdomains());

	//std::vector<SparseIJVMatrix<eslocal> > gBij(this->subdomains());

 //Dirichelt - seq version
	std::vector<size_t>::const_iterator vi;
	std::set<eslocal>::const_iterator si1, si2;

	std::vector<std::vector<eslocal> > t_lambda_map_sub_clst;
	std::vector<std::vector<eslocal> > t_lambda_map_sub_B1 ( this->subdomains() );
	std::vector<std::vector<double> >  t_B1_duplicity ( this->subdomains() );
	std::vector<std::vector<double> >  t_vec_c ( this->subdomains() );
	std::vector<SparseDOKMatrix<eslocal> > t_gB(this->subdomains());
	int 								t_lambda_count_B1 = 0;

	std::vector<eslocal> t_local_prim_numbering_d(this->subdomains(), 0);
	for (size_t i = 0; i < localBoundaries.size(); i++) {

		std::vector<bool> is_dirichlet(this->DOFs(), false);
		for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
			for (vi = properties.begin(); vi != properties.end(); ++vi) {
				const std::map<eslocal, double> &property = this->_input.mesh.coordinates().property(*vi).values();
				if (property.find(i) != property.end()) {


					double dirichlet_value = property.at(i);

					is_dirichlet[0] = true;

					t_gB[*si1](t_lambda_count_B1, t_local_prim_numbering_d[*si1] + (vi - properties.begin())) =  1.0;

					t_lambda_map_sub_clst.		push_back(std::vector<eslocal>({ t_lambda_count_B1, 0 }));
					t_lambda_map_sub_B1[*si1].	push_back(t_lambda_count_B1);
					t_B1_duplicity[*si1].		push_back(1.0);
					t_vec_c[*si1].				push_back(dirichlet_value);

					t_lambda_count_B1++;
				}
			}
		}

		for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
			t_local_prim_numbering_d[*si1] += this->DOFs();
		}

	}
#endif //SEQUENTIAL

#ifdef PARALLEL

	std::vector < std::vector < std::vector < eslocal > > > mat_B0 				 	 (threads);
	std::vector<eslocal> 									lambda_count_B0_l 		 (threads, 0);

#pragma omp parallel num_threads(threads)
	{
		int this_thread = omp_get_thread_num();
		int num_threads = omp_get_num_threads();

		std::vector<size_t>::const_iterator vi_l;
		std::set<eslocal>::const_iterator si1_l, si2_l;

		mat_B_dir[this_thread].			  	  clear();
		local_prim_numbering_d_l[this_thread].clear();

		mat_B_dir[this_thread].			  resize( this->subdomains() );
                mat_B0   [this_thread].                   resize( this->subdomains()  );
		local_prim_numbering_d_l[this_thread].resize( this->subdomains(),0 );

		int chunk_size; 
                if (num_threads > 1) 
                    chunk_size = localBoundaries.size() / (num_threads-1);
                else 
                    chunk_size = localBoundaries.size(); 
                    
		for (int d = 0; d < this->subdomains(); d++) {
			mat_B_dir[this_thread][d].reserve(chunk_size);
		}

//#pragma omp for schedule (static, chunk_size)
//		for (size_t i = 0; i < localBoundaries.size(); i++) {

		size_t lo =  this_thread      * chunk_size;
		size_t hi = (this_thread + 1) * chunk_size;

		if (hi > localBoundaries.size())    hi = localBoundaries.size();
		if (this_thread == num_threads - 1) hi = localBoundaries.size();

		for (size_t i = lo; i < hi; i++) {

			std::vector<bool> is_dirichlet(this->DOFs(), false);
			for (si1_l = localBoundaries[i].begin(); si1_l != localBoundaries[i].end(); ++si1_l) {
				for (vi_l = properties.begin(); vi_l != properties.end(); ++vi_l) {
					const std::map<eslocal, double> &property = this->_input.mesh.coordinates().property(*vi_l).values();
					if (property.find(i) != property.end()) {
						is_dirichlet[0] = true;
					}
				}
			}

			if (localBoundaries[i].size() > 1) {

				if ( globalBoundaries[i].size() == 1
					 && ( !localBoundaries.isCorner(i) || esconfig::solver::FETI_METHOD != 1 )
					 ) {	//TODO: Pozor Cornery nejdou do B1
					for (si1_l = localBoundaries[i].begin(); si1_l != localBoundaries[i].end(); ++si1_l) {
						for (si2_l = si1_l,++si2_l; si2_l != localBoundaries[i].end(); ++si2_l) {
							for (eslocal d = 0; d < this->DOFs(); d++) {
								if (is_dirichlet[d]) {
									continue;
								}

								mat_B_dir[this_thread][*si1_l].push_back( lambda_count_B1_l[this_thread] );
								mat_B_dir[this_thread][*si1_l].push_back( local_prim_numbering_d_l[this_thread][*si1_l] + d );
								mat_B_dir[this_thread][*si1_l].push_back( 1 * localBoundaries[i].size() );

								mat_B_dir[this_thread][*si2_l].push_back( lambda_count_B1_l[this_thread] );
								mat_B_dir[this_thread][*si2_l].push_back( local_prim_numbering_d_l[this_thread][*si2_l] + d );
								mat_B_dir[this_thread][*si2_l].push_back( -1 * localBoundaries[i].size() );

								lambda_count_B1_l[this_thread]++;
							}
						}
					}
				}

				//if ( 1 == 1 ) {
				if (localBoundaries.isCorner(i)) {
					for (si1_l = localBoundaries[i].begin(), si2_l = si1_l, ++si1_l; si1_l != localBoundaries[i].end(); ++si1_l) {
						for (eslocal d = 0; d < this->DOFs(); d++) {

							//lB[*si2](lambda_count_B0, local_prim_numbering[*si2] + d) =  1.0;
							//lB[*si1](lambda_count_B0, local_prim_numbering[*si1] + d) = -1.0;
							//_lambda_map_sub_B0[*si2].push_back(lambda_count_B0);
							//_lambda_map_sub_B0[*si1].push_back(lambda_count_B0);
							//lambda_count_B0++;

							mat_B0[this_thread][*si1_l].push_back( lambda_count_B0_l[this_thread] );
							mat_B0[this_thread][*si1_l].push_back( local_prim_numbering_d_l[this_thread][*si1_l] + d );
							mat_B0[this_thread][*si1_l].push_back( 1 );

							mat_B0[this_thread][*si2_l].push_back( lambda_count_B0_l[this_thread] );
							mat_B0[this_thread][*si2_l].push_back( local_prim_numbering_d_l[this_thread][*si2_l] + d );
							mat_B0[this_thread][*si2_l].push_back( -1 );

							lambda_count_B0_l[this_thread]++;

						}
						si2_l = si1_l;
					}
				}



			}

			for (si1_l = localBoundaries[i].begin(); si1_l != localBoundaries[i].end(); ++si1_l) {
				local_prim_numbering_d_l[this_thread][*si1_l] += this->DOFs();
			}

		}

#pragma omp barrier
#pragma omp single
		{
			for (size_t th = 1; th < threads; th++) {
				lambda_count_B1_l[th] += lambda_count_B1_l[th - 1];
				lambda_count_B0_l[th] += lambda_count_B0_l[th - 1];
				//std::cout << lambda_count_B1_l[th] << " ";
				for (size_t d = 0; d < this->subdomains(); d++) {
					local_prim_numbering_d_l[th][d] += local_prim_numbering_d_l[th - 1][d];
				}
			}
		}

#pragma omp for
		for (size_t d = 0; d < this->subdomains(); d++) {
			for (size_t th = 0; th < threads; th++) {
				eslocal lambda_off = 0;
				eslocal locDOF_off = 0;

				if (th > 0) lambda_off = lambda_count_B1_l[th-1];

				for (size_t i = 0; i < mat_B_dir[th][d].size(); i = i + 3 ) {
					if (th > 0)
						locDOF_off = local_prim_numbering_d_l[th-1][d];

					eslocal lambda = mat_B_dir[th][d][i + 0] + lambda_off;
				    eslocal locDOF = mat_B_dir[th][d][i + 1] + locDOF_off;

				    double duplicity;
				    double value;
				    eslocal tmp = mat_B_dir[th][d][i + 2];

				    if (tmp < 0) {
				    	duplicity = -1.0 / (double)tmp;
				    	value 	   = -1.0;
				    } else {
				    	duplicity = 1.0 / (double)tmp;
				    	value 	   = 1.0;
				    }

				    //gB[d]( lambda , locDOF  ) = value;
				    _B1[d].set ( lambda , locDOF, value);

					_lambda_map_sub_B1[d].	push_back(lambda);
					_B1_duplicity[d]. 		push_back(duplicity);
					_vec_c[d]. 				push_back(0.0);

				}

			}
		}


#pragma omp for
		for (size_t d = 0; d < this->subdomains(); d++) {
			for (size_t th = 0; th < threads; th++) {
				eslocal lambda_off = 0;
				eslocal locDOF_off = 0;

				if (th > 0) lambda_off = lambda_count_B0_l[th-1];

				for (size_t i = 0; i < mat_B0[th][d].size(); i = i + 3 ) {

					if (th > 0)
						locDOF_off = local_prim_numbering_d_l[th-1][d];

					eslocal lambda = mat_B0[th][d][i + 0] + lambda_off;
				    eslocal locDOF = mat_B0[th][d][i + 1] + locDOF_off;

				    double duplicity;
				    double value = (double) mat_B0[th][d][i + 2];

				    _B0[d].set ( lambda , locDOF, value);
					_lambda_map_sub_B0[d].	push_back(lambda);

				}
			}
		}

#pragma omp master
		{
			for (eslocal i = 0; i < lambda_count_B1_l[threads - 1]; i++) {
				_lambda_map_sub_clst.		push_back(std::vector<eslocal>({ i, 0 }));
			}
			lambda_count_B1 = lambda_count_B1_l[threads - 1];
			lambda_count_B0 = lambda_count_B0_l[threads - 1];
		}

	}


	for (eslocal d = 0; d < this->subdomains(); d++) {
//		//gB[d].resize(lambda_count_B1, local_prim_numbering_d_l[threads - 1][d]);
//		gBij[d].resize(lambda_count_B1, local_prim_numbering_d_l[threads - 1][d]);
//		//_B1[d] = gB[d];
//		_B1[d] = gBij[d];
////		lB[d].resize(lambda_count_B0, local_prim_numbering[threads - 1][d]);
////		_B0[d] = lB[d];
		_B0[d].resize(lambda_count_B0, local_prim_numbering_d_l[threads - 1][d]);
		_B1[d].resize(lambda_count_B1, local_prim_numbering_d_l[threads - 1][d]);

	}


#endif // PARALLEL


#ifdef SEQUENTIAL
	std::vector<eslocal> local_prim_numbering(this->subdomains(), 0);
	for (size_t i = 0; i < localBoundaries.size(); i++) {

		std::vector<bool> is_dirichlet(this->DOFs(), false);
		for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
			for (vi = properties.begin(); vi != properties.end(); ++vi) {
				const std::map<eslocal, double> &property = this->_input.mesh.coordinates().property(*vi).values();
				if (property.find(i) != property.end()) {
					//TODO: Toto asi neni dobre - jak resit DOFs ?
					is_dirichlet[0] = true;
				}
			}
		}

		if (localBoundaries[i].size() > 1) {

			if ( globalBoundaries[i].size() == 1
				 && ( !localBoundaries.isCorner(i) || esconfig::solver::FETI_METHOD != 1 )
				 ) {	//TODO: Pozor Cornery nejdou do B1
				for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
					for (si2 = si1,++si2; si2 != localBoundaries[i].end(); ++si2) {
						for (eslocal d = 0; d < this->DOFs(); d++) {
							if (is_dirichlet[d]) {
								continue;
							}
							t_gB[*si1](t_lambda_count_B1, local_prim_numbering[*si1] + d) =  1.0;
							t_gB[*si2](t_lambda_count_B1, local_prim_numbering[*si2] + d) = -1.0;

							t_lambda_map_sub_B1[*si1].push_back(t_lambda_count_B1);
							t_lambda_map_sub_B1[*si2].push_back(t_lambda_count_B1);
							t_lambda_map_sub_clst.	 push_back(std::vector<eslocal>({ t_lambda_count_B1, 0 }));
							t_B1_duplicity[*si1].	 push_back( 1.0 / (double) localBoundaries[i].size() );
							t_B1_duplicity[*si2].	 push_back( 1.0 / (double) localBoundaries[i].size() );
							t_vec_c[*si1].			 push_back(0.0);
							t_vec_c[*si2].			 push_back(0.0);

							t_lambda_count_B1++;
						}
					}
				}
			}

			//if ( 1 == 1 ) {
			if (localBoundaries.isCorner(i)) {
				for (si1 = localBoundaries[i].begin(), si2 = si1, ++si1; si1 != localBoundaries[i].end(); ++si1) {
					for (eslocal d = 0; d < this->DOFs(); d++) {
						lB[*si2](lambda_count_B0, local_prim_numbering[*si2] + d) =  1.0;
						lB[*si1](lambda_count_B0, local_prim_numbering[*si1] + d) = -1.0;
						_lambda_map_sub_B0[*si2].push_back(lambda_count_B0);
						_lambda_map_sub_B0[*si1].push_back(lambda_count_B0);
						lambda_count_B0++;
					}
					si2 = si1;
				}
			}

		}

		for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
			local_prim_numbering[*si1] += this->DOFs();
		}

	}

	for (eslocal d = 0; d < this->subdomains(); d++) {
		t_gB[d].resize(t_lambda_count_B1, local_prim_numbering[d]);
		_B1[d] = t_gB[d];

		lB[d].resize(lambda_count_B0, local_prim_numbering[d]);
		_B0[d] = lB[d];
	}

//	for (int i = 0; i < t_lambda_map_sub_B1.size(); i++ )
//	{
//		for (int j = 0; j < t_lambda_map_sub_B1[i].size(); j++)
//		{
//			if ( _lambda_map_sub_B1[i][j] != t_lambda_map_sub_B1[i][j] )
//				cout << i << ":" << j << " - " <<_lambda_map_sub_B1[i][j] << " - " << t_lambda_map_sub_B1[i][j] << std::endl;
//		}
//		//if ( _lambda_map_sub_B1[i].size() == t_lambda_map_sub_B1[i].size() )
//		//	cout << i << ":" <<_lambda_map_sub_B1[i].size() << " - " << t_lambda_map_sub_B1[i].size() << std::endl;
//	}

	_lambda_map_sub_B1 = t_lambda_map_sub_B1;
	_lambda_map_sub_clst = t_lambda_map_sub_clst;
	_B1_duplicity = t_B1_duplicity;
	_vec_c = t_vec_c;
#endif //SEQUENTIAL

	if ( this->rank() == 0) { std::cout << " B0 size: " <<  lambda_count_B0  << std::endl; }


}

template <class TInput>
void Gluing<TInput>::local_B1_global_resize() {

	int MPIrank = this->rank();
	int MPIsize = this->size();

	std::vector<SparseIJVMatrix<eslocal> > &B1 = _B1;

	// Local B1 - further processing - update row numbering based on all clusters
    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - start                                   "; system("date +%T.%6N"); }

	// Create lambda global numbering
	esglobal localB1_l_rows = B1[0].rows(); // number of rows B1 for all domains must be the same

	// Renumbering of the local B1 matrix including Dirichlet BC
	// Create lambda global counting for local B1
	esglobal global_B1_l_rows;
	MPI_Exscan(&localB1_l_rows, &global_B1_l_rows, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);
	if (MPIrank == 0)
		global_B1_l_rows = 0; //localB1_l_rows;

	global_B1_l_rows = localB1_l_rows + global_B1_l_rows;
	total_number_of_B1_l_rows = global_B1_l_rows;
	MPI_Bcast(&total_number_of_B1_l_rows, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);

    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - EXscan and Bcast done                  "; system("date +%T.%6N"); }

	cilk_for (eslocal domain_index=0; domain_index < this->subdomains(); domain_index++) {
		//TODO: lambda muze byt esglobal ale IJV matice je jen eslocal
		esglobal row_offset = global_B1_l_rows - localB1_l_rows;
		B1[domain_index].ShiftRowIndex(row_offset);

		eslocal  cols_tmp = B1[domain_index].columns();
		B1[domain_index].resize(total_number_of_B1_l_rows, cols_tmp); // TODO: prvni je esglobal

		for (eslocal i = 0; i < _lambda_map_sub_B1[domain_index].size(); i++) {
			_lambda_map_sub_B1[domain_index][i] += row_offset;
		}

	}

    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - end of renumbering of rows of local B1   "; system("date +%T.%6N"); }


	esglobal row_offset = global_B1_l_rows - localB1_l_rows;
	for (eslocal i = 0; i < _lambda_map_sub_clst.size(); i++) {
		_lambda_map_sub_clst[i][0] += row_offset;
		_lambda_map_sub_clst[i][1] = MPIrank;
	}

	// END - Local B1 - further processing - update row numbering based on all clusters

}

template <class TInput>
void Gluing<TInput>::get_myBorderDOFs_from_mesh() {

	int MPIrank = this->rank();
	int MPIsize = this->size();
	std::vector<SparseIJVMatrix<eslocal> > &B1 = _B1;
	std::vector<SparseIJVMatrix<eslocal> > &B0 = _B0;
	const mesh::Coordinates &coordinates = this->_input.mesh.coordinates();
	// TODO: cluster boundaries in surface
	const mesh::Boundaries &globalBoundaries = this->_input.mesh.clusterBoundaries();
	const mesh::Boundaries &localBoundaries = this->_input.mesh.subdomainBoundaries();


	eslocal borderDofNum; 	// number of DOFs on surface of my cluster
	std::vector< esglobal >::iterator it_vec;

	//std::vector < std::vector < esglobal > > neighBorderDofs;   // 2D vector for DOFs on borders of neighboring sub-domains
	//std::vector < esglobal > myBorderDOFs;  					// my border DOFs are here
	//std::vector < eslocal  > _neighClusters; 					// my neighboring clusters




	// Find all MPIranks of all neighboring clusters in globalBoundaries
	std::vector < eslocal > neigh_tmp  (MPIsize, 0);
	std::set<eslocal>::const_iterator it_set;
	std::set<eslocal>::const_iterator it_set_l;

    if (MPIrank == 0) { std::cout << " Global B - Blobal B1 neighdofs and neigh dofs array building             "; system("date +%T.%6N"); }

    // Now this loop is used to get my Border DOFs and my neighboring subdomains
    // information about neighnoring DOFS is here, but it is not used
	// neighBorderDofs.resize( MPIsize, std::vector< esglobal >( 0 , 0 ) );
	for (size_t i = 0; i < globalBoundaries.size(); i++) {

		//
		if ( globalBoundaries[i].size() > 1 && localBoundaries[i].size() == 1 ) {
			for (it_set = globalBoundaries[i].begin(); it_set != globalBoundaries[i].end(); ++it_set) {
				if (*it_set == MPIrank) { // if it points to local cluster
					for (int d_i = 0; d_i < this->DOFs(); d_i++ ) {
						_myBorderDOFs.push_back( this->DOFs() * coordinates.globalIndex(i) + d_i ); // mapping local local to global
					}
				} else {
					neigh_tmp[*it_set] = 1;
				}
	        }
		}

		//
		if ( globalBoundaries[i].size() > 1 && localBoundaries[i].size() > 1 ) {
			for (it_set = globalBoundaries[i].begin(); it_set != globalBoundaries[i].end(); ++it_set) {
				if (*it_set == MPIrank) {
					// this loop goes over the nodes that are in multiple domains on this one cluster
					for (it_set_l = localBoundaries[i].begin(); it_set_l != localBoundaries[i].end(); ++it_set_l) {
						for (int d_i = 0; d_i < this->DOFs(); d_i++ ) {
							_myBorderDOFs_sp      .push_back( this->DOFs() * coordinates.globalIndex(i) + d_i ); // mapping local local to global
							//myBorderDOFs_sp_loc_n.push_back( this->DOFs() *                          i + d_i ); // in local numbering
						}
					}
				}
	        }
		}

	}

	// Detecting neighboring clusters
	for (eslocal i = 0; i < neigh_tmp.size(); i++)
		if (neigh_tmp[i] == 1)
			_neighClusters.push_back(i);

	//std::sort (myBorderDOFs.begin(), myBorderDOFs.end());
	// removes the duplicit DOFs from myBorderDofs
	// POZOR
	//std::vector< esglobal >::iterator itx;
    //itx = std::unique (myBorderDOFs.begin(), myBorderDOFs.end());
    //myBorderDOFs.resize( std::distance(myBorderDOFs.begin(), itx) );

	// sort my neigh DOFs sp
	std::sort (_myBorderDOFs_sp.begin(), _myBorderDOFs_sp.end());

}

template <class TInput>
void Gluing<TInput>::computeClusterGluing(std::vector<size_t> &rows)
{
	int MPIrank = this->rank();
	int MPIsize = this->size();
	std::vector<SparseIJVMatrix<eslocal> > &B1 = _B1;
	std::vector<SparseIJVMatrix<eslocal> > &B0 = _B0;
	const mesh::Coordinates &coordinates = this->_input.mesh.coordinates();
	// TODO: cluster boundaries in surface
	const mesh::Boundaries &globalBoundaries = this->_input.mesh.clusterBoundaries();
	const mesh::Boundaries &localBoundaries = this->_input.mesh.subdomainBoundaries();

	bool flag_redund_lagr_mult = true;

	std::vector < SparseDOKMatrix<eslocal> > B1_DOK_tmp(this->subdomains());
	MPI_Request   mpi_req;
	MPI_Status 	  mpi_stat;
	std::set<eslocal>::const_iterator it_set;
	std::set<eslocal>::const_iterator it_set_l;


	// Update of local B1 matrices into global numbering
	this->local_B1_global_resize();

	// Get myBorderDOFs from mesh
	this->get_myBorderDOFs_from_mesh();



	// Total number of neighboring clusters
	eslocal neighClustNum;	// number of neighboring sub-domains for current sub-domainG
	neighClustNum = _neighClusters.size();

    if (MPIrank == 0) { std::cout << " Global B - myNeighDOFs arrays are transfered to neighbors                "; system("date +%T.%6N"); }

    MPI_Request * mpi_send_req  = new MPI_Request [_neighClusters.size()];
	MPI_Request * mpi_recv_req  = new MPI_Request [_neighClusters.size()];
	MPI_Status  * mpi_recv_stat = new MPI_Status  [_neighClusters.size()];

	std::vector < std::vector < esglobal > > neighBorderDofs;   // 2D vector for DOFs on borders of neighboring sub-domains
	neighBorderDofs.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );

	// sending my border DOFs to all neighboring clusters
	for (eslocal i = 0; i < _neighClusters.size(); i++) {
		MPI_Isend(&_myBorderDOFs[0], _myBorderDOFs.size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
	}

	// receiving all border DOFs from all neighboring sub-domains
	eslocal messages_received_bd = 0;
	while ( messages_received_bd < _neighClusters.size() ) {
		for (eslocal i = 0; i < _neighClusters.size(); i++) {

			int flag;
			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

			if (flag) {
				int recv_msg_size = 0;

				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);

				neighBorderDofs[i].resize(recv_msg_size);
				if (recv_msg_size < 0)
					std::cout << "error in msg size !";
				MPI_Recv(&neighBorderDofs[i][0], recv_msg_size, esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);

				messages_received_bd++;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);


	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }


	// NOW :
	// my neighboring clusters are in :     	_neighClusters
	// my border DOFs are in : 					myBorderDOFs
	// neighboring clusters border DOFs are in:	neighBorderDofs[i]

 	std::vector < std::vector < esglobal > > myNeighsSparse;
	myNeighsSparse.resize( _myBorderDOFs.size() );

	for (eslocal j = 0; j < _myBorderDOFs.size(); j++) {
		myNeighsSparse[j].reserve(5);
		myNeighsSparse[j].push_back(_myBorderDOFs[j]);
	}

	for (eslocal i = 0; i < _neighClusters.size(); i++) {
		int d = _neighClusters[i];
		int mydofs_index = 0;
		int nedofs_index = 0;

		if ( i == 0 && MPIrank < _neighClusters[i] ) {
			for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank);
		}

		do {

			if ( neighBorderDofs[i][nedofs_index] == _myBorderDOFs[mydofs_index] ) {
				myNeighsSparse[mydofs_index].push_back(d);
				mydofs_index++;
				nedofs_index++;
			} else {
				if ( neighBorderDofs[i][nedofs_index] > _myBorderDOFs[mydofs_index] ) {
					mydofs_index++;
				} else {
					nedofs_index++;
				}
			}

		} while (mydofs_index < _myBorderDOFs.size() && nedofs_index < neighBorderDofs[i].size() );


		if ( i < _neighClusters.size() - 1)
			if ( MPIrank > _neighClusters[i] && MPIrank < _neighClusters[i+1] )
				for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
					myNeighsSparse[j].push_back(MPIrank);



		if ( i == _neighClusters.size() -1 && MPIrank > _neighClusters[i] )
			for (eslocal j = 0; j < _myBorderDOFs.size(); j++)
				myNeighsSparse[j].push_back(MPIrank);

	}

	for (eslocal i = 0; i < myNeighsSparse.size(); i++)
		if (myNeighsSparse[i].size() < 3)
			myNeighsSparse[i].clear();

	if (MPIrank == 0) { std::cout << " Global B - myNeighSparse assembled                                       "; system("date +%T.%6N"); }

	std::vector < std::vector < esglobal > > myLambdas;

	esglobal lambdaCount = 0;

	for (eslocal i = 0; i < myNeighsSparse.size(); i++) {
		bool add = false;
		if (myNeighsSparse[i].size() > 0) {
			esglobal dof = myNeighsSparse[i][0];

			eslocal cnt_tmp = myNeighsSparse[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb

			for (eslocal j = 1; j < myNeighsSparse[i].size(); j++) {
				esglobal neighSD = myNeighsSparse[i][j];

				if ( add ) {

					myLambdas.push_back ( std::vector < esglobal > () );
					myLambdas[lambdaCount].reserve(6);
					myLambdas[lambdaCount].resize(5);
					myLambdas[lambdaCount][0] = lambdaCount;			// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
					myLambdas[lambdaCount][1] = dof;					// dof in global numbering
					myLambdas[lambdaCount][2] = (esglobal)MPIrank;		// my sub-domainG
					myLambdas[lambdaCount][3] = neighSD;				// neigh. sub-domainG
					myLambdas[lambdaCount][4] = (esglobal)cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb
					lambdaCount++;

					if (!flag_redund_lagr_mult)
						break;
				}

				if (neighSD == (esglobal)MPIrank)
					add = true;
			}
		}
	}

	if (MPIrank == 0) { std::cout << " Global B - Create global lambda numbering                                "; system("date +%T.%6N"); }

	esglobal lambdaGlobalCount = 0;

	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);

	if (MPIrank == 0) lambdaGlobalCount = 0;
	esglobal lambdaNum = lambdaCount + lambdaGlobalCount;
	MPI_Bcast(&lambdaNum, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
	esglobal total_number_of_global_B1_lambdas = lambdaNum;

	if (myLambdas.size() > 0)
		//cilk_
		for (eslocal i = 0; i < myLambdas.size(); i++)
			myLambdas[i][0] = myLambdas[i][0]  + lambdaGlobalCount + total_number_of_B1_l_rows; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index


	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }


	std::vector < std::vector < esglobal > > mpi_send_buff;
	mpi_send_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
//#ifdef DEBUG
//	for (int i = 0; i < _neighClusters.size(); i++) {
//#else
	cilk_for (int i = 0; i < _neighClusters.size(); i++) {
//#endif
		int index = 0;
		if ( myLambdas.size() > 0 )
		{
			for (int j = 0; j < myLambdas.size(); j++) {
				if( myLambdas[j][3] == _neighClusters[i] ) {
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

	if (MPIrank == 0) { std::cout << " Global B - Isend                                                         "; system("date +%T.%6N"); }

	for (int i = 0; i < _neighClusters.size(); i++)
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);


	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }

	std::vector < std::vector < esglobal > > mpi_recv_buff;
	mpi_recv_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
	//delete [] mpi_send_req;


	// receiving all border DOFs from all neighboring sub-domains
	eslocal messages_received = 0;
	while ( messages_received < _neighClusters.size() ) {
		for (eslocal i = 0; i < _neighClusters.size(); i++) {

			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
			int flag;
			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

			if (flag) {
				int recv_msg_size = 0;

				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);

				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);

				free(mpi_tmp_recv_buff);
				messages_received++;
			}
		}
	}

	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }

	// decode received lambdas
	eslocal recv_lamba_count = 0;
	for (eslocal i = 0; i < _neighClusters.size(); i++)
		if (mpi_recv_buff[i].size() > 1)
			recv_lamba_count += mpi_recv_buff[i].size();

	recv_lamba_count = recv_lamba_count / 5;

	eslocal l_i = myLambdas.size();
	myLambdas.resize( myLambdas.size() + recv_lamba_count, std::vector< esglobal >( 5 , 0 ) );

	for (eslocal i = 0; i < _neighClusters.size(); i++) {
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

	if (MPIrank == 0) { std::cout << " Global B - myLambdas - sort or tbb:sort                                  "; system("date +%T.%6N"); }

	auto comp_vf = [](const std::vector<esglobal> &a, const std::vector<esglobal> &b) {return a[0] < b[0];};
//#define USE_TBB
#ifdef USE_TBB
	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
#else
	std::sort         (myLambdas.begin(), myLambdas.end(), comp_vf);
#endif

	if (MPIrank == 0) { std::cout << " Global B - Final B assembling with g2l mapping using std::map            "; system("date +%T.%6N"); }

	esglobal lambda;
	esglobal DOFNumber;
	eslocal n_myLambdas = myLambdas.size();

#ifdef USE_TBB
	tbb::mutex m;
	cilk_for (int j = 0; j < n_myLambdas; j++)
#else
	for (eslocal j = 0; j < n_myLambdas; j++)
#endif
	{

		esglobal lambda         = myLambdas[j][0];
		esglobal DOFNumber      = myLambdas[j][1];

		esglobal dofNODEnumber = DOFNumber / this->DOFs();
		eslocal  dofNODEoffset = DOFNumber % this->DOFs();
		eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );

		const std::set < eslocal > & subs_with_element = this->_input.mesh.subdomainBoundaries()[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel

		eslocal  cnt            = myLambdas[j][4];
		double B_value;
		if (myLambdas[j][2] < myLambdas[j][3])
			B_value = 1.0;
		else
			B_value = -1.0;

		//TODO:Tuto smycku prepsat, aby fungovala jen podoblasti, ve ktery je uzel
		// Pomala a otestovana verze
		//for (eslocal d = 0; d < this->subdomains(); d++) {
 
        //TODO: rychle verze nutno otestovat na realne uloze
		for (it_set = localBoundaries[clustDofNODENumber].begin(); it_set != localBoundaries[clustDofNODENumber].end(); ++it_set) {
			eslocal d = *it_set;

			eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
			if ( domDofNODENumber != -1 ) {
				//vychazi z mapovani g2l pro uzly a ne DOFy
				eslocal domDOFNumber = this->DOFs() * domDofNODENumber + dofNODEoffset;
			   #ifdef USE_TBB
				m.lock();
			   #endif
				B1_DOK_tmp[d](lambda, domDOFNumber) = B_value;
				myLambdas[j].push_back(d);
				_B1_duplicity[d].push_back(1.0 / cnt);
				_vec_c[d].push_back(0.0);
			   #ifdef USE_TBB
				m.unlock();
			   #endif
				break;

			}
		}
	}

	cilk_for (eslocal d = 0; d < this->subdomains(); d++){
		//TODO: lambdaNum muze byt 64bit integer
		//TODO: matice B - pocet radku muze byt 64bit int
		B1_DOK_tmp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas , rows[d]);
		SparseIJVMatrix<eslocal> ijv = B1_DOK_tmp[d];
		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
	}


	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }


	// for global B1
	for (int i = 0; i < myLambdas.size(); i++) {
		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
		tmp_vec[0] = myLambdas[i][0];
		tmp_vec[1] = myLambdas[i][2];
		tmp_vec[2] = myLambdas[i][3];

		_lambda_map_sub_clst.push_back(tmp_vec);

		eslocal lam_tmp = myLambdas [i][0]; //TODO: esglobal
		_lambda_map_sub_B1[ myLambdas[i][5] ].push_back( lam_tmp );
	}

	if (MPIrank == 0) { std::cout << " Dual size: " <<  total_number_of_B1_l_rows + total_number_of_global_B1_lambdas  << std::endl; }

        MPI_Barrier(MPI_COMM_WORLD);

   	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }


// *****************************************************************************************************************************************************
// Creating B on corners of clusters and domains


    	flag_redund_lagr_mult = true;
    	eslocal neighClustNum_sp = 0;	// number of neighboring sub-domains for current sub-domainG
    	eslocal borderDofNum_sp  = 0; 	// number of DOFs on surface on my cluster

    	std::vector < std::vector < esglobal > > neighBorderDofs_sp;    // 2D vector for DOFs on borders of neighboring sub-domains
    	//std::vector < esglobal > 				 myBorderDOFs_sp;  	    // my border DOFs are here - in global numbering
        std::vector < esglobal > 				 myBorderDOFs_sp_nr;  	// my border DOFs are here - in global numbering
        std::vector < esglobal > 				 myBorderDOFs_sp_loc_n; // my border DOFs are here - in local numbering


        myBorderDOFs_sp_nr = _myBorderDOFs_sp;
    	neighClustNum 	   = _neighClusters.size();

        // removes the duplicit DOFs from myBorderDofs
    	std::vector< esglobal >::iterator itx;
        itx = std::unique (myBorderDOFs_sp_nr.begin(), myBorderDOFs_sp_nr.end());
        myBorderDOFs_sp_nr.resize( std::distance(myBorderDOFs_sp_nr.begin(), itx) );

        if (MPIrank == 0) { std::cout << " Global B - myNeighDOFs arrays are transfered to neighbors                "; system("date +%T.%6N"); }

       	neighBorderDofs_sp.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );

    	// sending my DOFs on border to all neighboring sub-domains
    	for (eslocal i = 0; i < _neighClusters.size(); i++) {
    		MPI_Isend(&_myBorderDOFs_sp[0], _myBorderDOFs_sp.size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
    	}

    	// receiving all border DOFs from all neighboring sub-domains
    	eslocal messages_received_spp = 0;
    	while ( messages_received_spp < _neighClusters.size() ) {
    		for (eslocal i = 0; i < _neighClusters.size(); i++) {

    			int flag;
    			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

    			if (flag) {
    				int recv_msg_size = 0;

    				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);

    				neighBorderDofs_sp[i].resize(recv_msg_size);

    				MPI_Recv(&neighBorderDofs_sp[i][0], recv_msg_size,esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);

    				messages_received_spp++;
    			}
    		}
    	}

    	MPI_Barrier(MPI_COMM_WORLD);

    	// END - Find all MPIranks of all neighboring clusters in globalBoundaries

    	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }

    	// NOW :
    	// my neighboring sub-domains are in : 	_neighClusters
    	// my border DOFs are in : 				myBorderDOFs
    	// neighboring sub-domains border DOFs are in neighBorderDofs[i]

    	// removes the duplicit DOFs from myBorderDofs
    	//std::vector< esglobal >::iterator it;
        //it = std::unique (myBorderDOFs.begin(), myBorderDOFs.end());
        //myBorderDOFs.resize( std::distance(myBorderDOFs.begin(), it) );

     	std::vector < std::vector < esglobal > > myNeighsSparse_sp;
    	myNeighsSparse_sp.resize( myBorderDOFs_sp_nr.size() );

    	for (eslocal j = 0; j < myBorderDOFs_sp_nr.size(); j++) {
   			myNeighsSparse_sp[j].reserve(5);
   			myNeighsSparse_sp[j].push_back(myBorderDOFs_sp_nr[j]);
    	}


    	//////////////////

    	std::vector < std::vector < eslocal  > > neighBorderDofs_sp_cnt;    	// number of duplicity in each element
    	std::vector < std::vector < esglobal > > neighBorderDofs_sp_red;    	// neigh DOFs wo duplicity

    	neighBorderDofs_sp_cnt.resize( neighClustNum, std::vector< eslocal >( 0 , 0 ) );
    	neighBorderDofs_sp_red.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );

    	std::vector               < eslocal  > 	myBorderDOFs_sp_cnt;    	// the same but for neigh DOFs
    	std::vector               < esglobal > 	myBorderDOFs_sp_red;    	// my border dofs w/o duplicity

       	eslocal dup_cnt_sp = 1;
       	for (eslocal i = 1; i < _myBorderDOFs_sp.size(); i++) {

       		if (_myBorderDOFs_sp[i-1] != _myBorderDOFs_sp[i]) {
       			myBorderDOFs_sp_cnt.push_back(dup_cnt_sp);
       			myBorderDOFs_sp_red.push_back(_myBorderDOFs_sp[i-1]);
       			dup_cnt_sp = 1;
       		} else {
       			dup_cnt_sp++;
       		}

       	}

       	if (dup_cnt_sp > 1) {
   			myBorderDOFs_sp_cnt.push_back(dup_cnt_sp);
   			myBorderDOFs_sp_red.push_back(_myBorderDOFs_sp[_myBorderDOFs_sp.size() - 1]);
       	}

    	for (eslocal j = 0; j < _neighClusters.size(); j++) {
			dup_cnt_sp = 1;
			for (eslocal i = 1; i < neighBorderDofs_sp[j].size(); i++) {
				if (neighBorderDofs_sp[j][i-1] != neighBorderDofs_sp[j][i]) {
					neighBorderDofs_sp_cnt[j].push_back(dup_cnt_sp);
					neighBorderDofs_sp_red[j].push_back(neighBorderDofs_sp[j][i-1]);
					dup_cnt_sp = 1;
				} else {
					dup_cnt_sp++;
				}
			}

			if (dup_cnt_sp > 1 ){
				neighBorderDofs_sp_cnt[j].push_back(dup_cnt_sp);
                neighBorderDofs_sp_red[j].push_back(neighBorderDofs_sp[j][neighBorderDofs_sp[j].size() - 1]);
			}
    	}


    	for (eslocal i = 0; i < _neighClusters.size(); i++) {

    		int d = _neighClusters[i];
    		int mydofs_index = 0;
    		int nedofs_index = 0;

    		if ( i == 0 && MPIrank < _neighClusters[i] ) {
    			for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
    				for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
    					myNeighsSparse_sp[j].push_back(MPIrank);
    		}

    		do {
    			if ( neighBorderDofs_sp_red[i].size() > 0 && myBorderDOFs_sp_red.size() > 0 )
    				if ( neighBorderDofs_sp_red[i][nedofs_index] == myBorderDOFs_sp_red[mydofs_index] ) {
    					for (eslocal j = 0; j < neighBorderDofs_sp_cnt[i][nedofs_index]; j++)
    						myNeighsSparse_sp[mydofs_index].push_back(d);
    					mydofs_index++;
    					nedofs_index++;
    				} else {
    					if ( neighBorderDofs_sp_red[i][nedofs_index] > myBorderDOFs_sp_red[mydofs_index] ) {
    						mydofs_index++;
    					} else {
    						nedofs_index++;
    					}
    				}

    		} while (mydofs_index < myBorderDOFs_sp_red.size() && nedofs_index < neighBorderDofs_sp_red[i].size() );

    		if ( i < _neighClusters.size() - 1)
    			if ( MPIrank > _neighClusters[i] && MPIrank < _neighClusters[i+1] )
    				for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
    					for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
    						myNeighsSparse_sp[j].push_back(MPIrank);


    		if ( i == _neighClusters.size() -1 && MPIrank > _neighClusters[i] )
    			for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
    				for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
    					myNeighsSparse_sp[j].push_back(MPIrank);
    	}


    	for (eslocal i = 0; i < myNeighsSparse_sp.size(); i++)
    		if (myNeighsSparse_sp[i].size() < 3)
    			myNeighsSparse_sp[i].clear();


    	if (MPIrank == 0) { std::cout << " Global B - myNeighSparse assembled                                       "; system("date +%T.%6N"); }

    	std::vector < std::vector < esglobal > > myLambdas_sp;

    	esglobal lambdaCount_sp = 0;

    	for (eslocal i = 0; i < myNeighsSparse_sp.size(); i++) {
    		bool add = false;
    		if (myNeighsSparse_sp[i].size() > 0) {

    			esglobal dof = myNeighsSparse_sp[i][0];
    			eslocal  cnt_tmp = myNeighsSparse_sp[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb

    			int min_index = 1;
    			int max_index = myNeighsSparse_sp[i].size();

    			for (eslocal i_t = 1; i_t < myNeighsSparse_sp[i].size(); i_t++) {
    				if ( myNeighsSparse_sp[i][i_t] == esglobal(MPIrank) ) {
    					min_index = i_t;
    					break;
    				}
    			}

    			for (eslocal i_t = 1; i_t < myNeighsSparse_sp[i].size(); i_t++) {
    				if ( myNeighsSparse_sp[i][i_t] > esglobal(MPIrank) ) {
    					max_index = i_t;
    					break;
    				}
    			}


    			for (eslocal k = min_index; k < max_index; k++) {
    				int tmp_cl = k-min_index;
					for (eslocal j = k + 1; j < myNeighsSparse_sp[i].size(); j++) {

						tmp_cl++;
						if (myNeighsSparse_sp[i][j] != myNeighsSparse_sp[i][j-1])
							tmp_cl = 0;

						esglobal neighSD = myNeighsSparse_sp[i][j];

						//if ( add ) {

							myLambdas_sp.push_back ( std::vector < esglobal > () );
							myLambdas_sp[lambdaCount_sp].reserve(7);
							myLambdas_sp[lambdaCount_sp].resize(7);
							myLambdas_sp[lambdaCount_sp][0] = lambdaCount_sp;	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
							myLambdas_sp[lambdaCount_sp][1] = dof;			// dof in global numbering
							myLambdas_sp[lambdaCount_sp][2] = (esglobal)MPIrank;		// my sub-domainG
							myLambdas_sp[lambdaCount_sp][3] = neighSD;		// neigh. sub-domainG
							myLambdas_sp[lambdaCount_sp][4] = (esglobal)cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb
							myLambdas_sp[lambdaCount_sp][5] = (esglobal)(k - min_index);
							myLambdas_sp[lambdaCount_sp][6] = (esglobal)tmp_cl;
							lambdaCount_sp++;

							//if (!flag_redund_lagr_mult)
							//	break;
						//}

						//if (neighSD == (esglobal)MPIrank)
						//	add = true;
					}
					//add = false;
    			}
    		}
    	}

    	if (MPIrank == 0) { std::cout << " Global B - Create global lambda numbering                                "; system("date +%T.%6N"); }

    	esglobal lambdaGlobalCount_sp = 0;

    	MPI_Exscan(&lambdaCount_sp, &lambdaGlobalCount_sp, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);

    	if (MPIrank == 0) lambdaGlobalCount_sp = 0;
    	esglobal lambdaNum_sp = lambdaCount_sp + lambdaGlobalCount_sp;
    	MPI_Bcast(&lambdaNum_sp, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
    	esglobal total_number_of_global_B1_lambdas_sp = lambdaNum_sp;

    	if (myLambdas_sp.size() > 0)
    		//cilk_
    		for (eslocal i = 0; i < myLambdas_sp.size(); i++)
    			myLambdas_sp[i][0] = myLambdas_sp[i][0]  + lambdaGlobalCount_sp + total_number_of_global_B1_lambdas + total_number_of_B1_l_rows; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index

    	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }

    	mpi_send_buff.clear();
    	mpi_send_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );

    	cilk_for (int i = 0; i < _neighClusters.size(); i++) {
    		int index = 0;
    		if ( myLambdas_sp.size() > 0 )
    		{
    			for (int j = 0; j < myLambdas_sp.size(); j++) {
    				if( myLambdas_sp[j][3] == _neighClusters[i] ) {
    					mpi_send_buff[i].push_back(myLambdas_sp[j][0]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][1]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][2]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][3]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][4]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][5]);
    					mpi_send_buff[i].push_back(myLambdas_sp[j][6]);
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

    	if (MPIrank == 0) { std::cout << " Global B - Isend                                                         "; system("date +%T.%6N"); }

    	for (int i = 0; i < _neighClusters.size(); i++)
    		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);


    	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }

    	mpi_recv_buff.clear();
        mpi_recv_buff.resize( _neighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
    	delete [] mpi_send_req;


    	// receiving all border DOFs from all neighboring sub-domains
    	eslocal messages_received_sp = 0;
    	while ( messages_received_sp < _neighClusters.size() ) {
    		for (eslocal i = 0; i < _neighClusters.size(); i++) {

    			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
    			int flag;
    			MPI_Iprobe( _neighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

    			if (flag) {
    				int recv_msg_size = 0;

    				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
    				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);

    				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, _neighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
    				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);

    				free(mpi_tmp_recv_buff);
    				messages_received_sp++;
    			}
    		}
    	}

    	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }

    	// decode received lambdas
    	eslocal recv_lamba_count_sp = 0;
    	for (eslocal i = 0; i < _neighClusters.size(); i++)
    		if (mpi_recv_buff[i].size() > 1)
    			recv_lamba_count_sp += mpi_recv_buff[i].size();

    	recv_lamba_count_sp = recv_lamba_count_sp / 7;

    	eslocal l_i_sp = myLambdas_sp.size();
    	myLambdas_sp.resize( myLambdas_sp.size() + recv_lamba_count_sp, std::vector< esglobal >( 7 , 0 ) );

    	for (eslocal i = 0; i < _neighClusters.size(); i++) {
    		if (mpi_recv_buff[i].size() > 1) {
    			for (int j = 0; j < mpi_recv_buff[i].size() / 7; j++) {
    				myLambdas_sp[l_i_sp][0] = mpi_recv_buff[i][7*j + 0];
    				myLambdas_sp[l_i_sp][1] = mpi_recv_buff[i][7*j + 1];
    				myLambdas_sp[l_i_sp][2] = mpi_recv_buff[i][7*j + 3];
    				myLambdas_sp[l_i_sp][3] = mpi_recv_buff[i][7*j + 2];
    				myLambdas_sp[l_i_sp][4] = mpi_recv_buff[i][7*j + 4];
    				myLambdas_sp[l_i_sp][5] = mpi_recv_buff[i][7*j + 6];
    				myLambdas_sp[l_i_sp][6] = mpi_recv_buff[i][7*j + 5];
    				l_i_sp++;
    			}
    		}
    	}

    	if (MPIrank == 0) { std::cout << " Global B - myLambdas - sort or tbb:sort                                  "; system("date +%T.%6N"); }

    //	auto comp_vf = [](const std::vector<esglobal> &a, const std::vector<esglobal> &b) {return a[0] < b[0];};

    //#ifdef USE_TBB
    //	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
    //#else
    	std::sort         (myLambdas_sp.begin(), myLambdas_sp.end(), comp_vf);
    //#endif

    	if (MPIrank == 0) { std::cout << " Global B - Final B assembling with g2l mapping using std::map            "; system("date +%T.%6N"); }

    	esglobal lambda_sp;
    	esglobal DOFNumber_sp;
    	eslocal  n_myLambdas_sp = myLambdas_sp.size();

    	std::vector < SparseDOKMatrix<eslocal> > B1_DOK_sp(this->subdomains());

    //#ifdef USE_TBB
    //	tbb::mutex m;
    //	cilk_for (int j = 0; j < n_myLambdas; j++)
    //#else
    	for (eslocal j = 0; j < n_myLambdas_sp; j++)
    //#endif
    	{

    		esglobal lambda_sp         = myLambdas_sp[j][0];
    		esglobal DOFNumber_sp      = myLambdas_sp[j][1];

    		//vychazi z mapovani g2l pro uzly a ne DOFy
    		esglobal dofNODEnumber      = DOFNumber_sp / this->DOFs();
    		eslocal  dofNODEoffset      = DOFNumber_sp % this->DOFs();

    		if ( myLambdas_sp[j][2] == myLambdas_sp[j][3] ) { // resim vazby mezi domenama uvnitr clusteru

    			eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );
				const std::set    < eslocal >  & subs_with_element = localBoundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel
				std::vector < eslocal >    subs_with_elem;
				for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set)
					subs_with_elem.push_back( *it_set );

				eslocal d1                = subs_with_elem[ myLambdas_sp[j][5] ];
				eslocal d2                = subs_with_elem[ myLambdas_sp[j][6] ];

				eslocal domDofNODENumber1 =coordinates.localIndex(clustDofNODENumber, d1);
				eslocal domDofNODENumber2 =coordinates.localIndex(clustDofNODENumber, d2);

				eslocal domDOFNumber1     = this->DOFs() * domDofNODENumber1 + dofNODEoffset;
				eslocal domDOFNumber2     = this->DOFs() * domDofNODENumber2 + dofNODEoffset;

				B1_DOK_sp[d1](lambda_sp, domDOFNumber1) =  1.0;
				B1_DOK_sp[d2](lambda_sp, domDOFNumber2) = -1.0;

				myLambdas_sp[j].push_back(d1);
				myLambdas_sp[j].push_back(d2);

				eslocal  cnt            = myLambdas_sp[j][4];
				_B1_duplicity[d1].push_back(1.0 / cnt);
				_B1_duplicity[d2].push_back(1.0 / cnt);

				_vec_c[d1].push_back(0.0);
				_vec_c[d2].push_back(0.0);

	    		std::vector < eslocal > tmp_vec1 (2,0);		//TODO: must be esglobal
	    		tmp_vec1[0] = myLambdas_sp[j][0];
	    		tmp_vec1[1] = myLambdas_sp[j][2];

	    		_lambda_map_sub_clst.push_back(tmp_vec1);

	    		std::vector < eslocal > tmp_vec2 (2,0);		//TODO: must be esglobal
	    		tmp_vec2[0] = myLambdas_sp[j][0];
	    		tmp_vec2[1] = myLambdas_sp[j][3];

	    		_lambda_map_sub_clst.push_back(tmp_vec2);

	    		eslocal lam_tmp = myLambdas_sp [j][0]; 		//TODO: esglobal
	    		_lambda_map_sub_B1[ d1 ].push_back( lam_tmp );
	        	_lambda_map_sub_B1[ d2 ].push_back( lam_tmp );


    		} else { // resim vazby mezi clustery

    			eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );
				const std::set    < eslocal >  & subs_with_element = localBoundaries[clustDofNODENumber]; //globalBoundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel
				std::vector < eslocal >    subs_with_elem;
				for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set)
					subs_with_elem.push_back( *it_set );


				eslocal  cnt            = myLambdas_sp[j][4];
				double B_value;
				if (myLambdas_sp[j][2] < myLambdas_sp[j][3])
					B_value = 1.0;
				else
					B_value = -1.0;

				eslocal d                = subs_with_elem[ myLambdas_sp[j][5] ];
				//eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
				eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
				eslocal domDOFNumber     = this->DOFs() * domDofNODENumber + dofNODEoffset;
				// #ifdef USE_TBB
				// m.lock();
				// #endif
				B1_DOK_sp[d](lambda_sp, domDOFNumber) = B_value;
				myLambdas_sp[j].push_back(d);
				_B1_duplicity[d].push_back(1.0 / cnt);
				_vec_c[d].push_back(0.0);
				// #ifdef USE_TBB
				// m.unlock();
				// #endif

	    		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
	    		tmp_vec[0] = myLambdas_sp[j][0];
	    		tmp_vec[1] = myLambdas_sp[j][2];
	    		tmp_vec[2] = myLambdas_sp[j][3];

	    		_lambda_map_sub_clst.push_back(tmp_vec);

	    		eslocal lam_tmp = myLambdas_sp [j][0]; 		//TODO: esglobal
	    		_lambda_map_sub_B1[ d ].push_back( lam_tmp );

    		}


    	}


    	for (eslocal d = 0; d < this->subdomains(); d++){
    		//TODO: lambdaNum muze byt 64bit integer
    		//TODO: matice B - pocet radku muze byt 64bit int
    		B1_DOK_sp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas + total_number_of_global_B1_lambdas_sp , rows[d]);
    		SparseIJVMatrix<eslocal> ijv = B1_DOK_sp[d];
    		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
    	}

    	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }

    	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }
    	if (MPIrank == 0) { std::cout << " END - Creating lambda_map_sub vector of vectors - Global B1              "; system("date +%T.%6N"); }

    	MPI_Barrier(MPI_COMM_WORLD);

        if (MPIrank == 0) { std::cout << " Dual size: " <<  total_number_of_B1_l_rows + total_number_of_global_B1_lambdas + total_number_of_global_B1_lambdas_sp  << std::endl; }

    	cilk_for (int i = 0; i < _vec_c.size(); i++) {
    		 _vec_c[i].resize(_B1_duplicity[i].size(), 0.0);
    	}

    	cilk_for (size_t p = 0; p < this->subdomains(); p++) {
    		_localB[p]  = _B0[p];
    		_globalB[p] = _B1[p];
    	}
}

}
