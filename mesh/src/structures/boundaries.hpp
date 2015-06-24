
#include "boundaries.h"
#include <mpi.h>

namespace mesh {

template<typename T>
void Boundaries::create_B1_l(	std::vector < SparseIJVMatrix >         & B1_local,
								std::vector < SparseIJVMatrix >         & B0_local,
								std::vector < std::vector <T> >    		& l2g_vec,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B0,
								std::vector < std::vector <double> >    & B1_l_duplicity,
								const eslocal domains_num)
{

	const std::map<eslocal, double> &dirichlet_x = _mesh.coordinates().property(CP::DIRICHLET_X).values();
	const std::map<eslocal, double> &dirichlet_y = _mesh.coordinates().property(CP::DIRICHLET_Y).values();
	const std::map<eslocal, double> &dirichlet_z = _mesh.coordinates().property(CP::DIRICHLET_Z).values();

	l2g_vec.resize(domains_num);

	std::vector < SparseDOKMatrix > B1_loc(domains_num);
	std::vector < SparseDOKMatrix > B0_loc(domains_num);

	lambda_map_sub_B1.resize(domains_num);
	lambda_map_sub_B0.resize(domains_num);
	B1_l_duplicity.resize(domains_num);
	std::vector<eslocal> local_prim_numbering(domains_num, 0);
	B1_local.resize(domains_num);
	B0_local.resize(domains_num);

	std::set<eslocal>::const_iterator it;
	std::set<eslocal>::const_iterator it1;
	std::set<eslocal>::const_iterator it2;

	eslocal lambda_count_B1 = 0;
	eslocal lambda_count_B0 = 0;

	for (T i = 0; i < _boundaries.size(); i++) {
		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			if ( dirichlet_x.find(index(i)) != dirichlet_x.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 0) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector <eslocal> (1, lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( dirichlet_y.find(index(i)) != dirichlet_y.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 1) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
			if ( dirichlet_z.find(index(i)) != dirichlet_z.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 2) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				lambda_map_sub_clst.push_back( std::vector < eslocal > (1, lambda_count_B1) );
				B1_l_duplicity[*it].push_back( 1.0 / (double)_boundaries[i].size() );
				lambda_count_B1++;
			}
		}

		if ( _boundaries[i].size() > 1 ) {

			// with duplicity
			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
					for (eslocal d_i = 0; d_i < 3; d_i++) {
						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0;
						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0;

						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
						lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );

						B1_l_duplicity[*it1].push_back( 1.0 / (double)_boundaries[i].size() );
						B1_l_duplicity[*it2].push_back( 1.0 / (double)_boundaries[i].size() );

						lambda_count_B1++;
					}
				}
			}


			// no duplicity
			bool is_corner = true;
			if ( is_corner ) {
				for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
					if (it1 != _boundaries[i].begin()) {
						for (eslocal d_i = 0; d_i < 3; d_i++) {
							B0_loc[*it2](lambda_count_B0, local_prim_numbering[*it2] + d_i) =  1.0;
							B0_loc[*it1](lambda_count_B0, local_prim_numbering[*it1] + d_i) = -1.0;
							lambda_map_sub_B0[*it2].push_back(lambda_count_B0);
							lambda_map_sub_B0[*it1].push_back(lambda_count_B0);
							lambda_count_B0++;
						}
					}
					it2 = it1;
				}
			}
		}
		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			local_prim_numbering[*it] += 3;
			l2g_vec[*it].push_back(index(i));
		}
	}

	for (eslocal d = 0; d < domains_num; d++) {
		B1_loc[d].resize(lambda_count_B1, local_prim_numbering[d]);
		B1_local[d] = B1_loc[d];

		B0_loc[d].resize(lambda_count_B0, local_prim_numbering[d]);
		B0_local[d] = B0_loc[d];
	}
}



template<typename T>
void Boundaries::create_B1_g(	std::vector < SparseIJVMatrix >         & B1_local,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
								std::vector < std::vector <double> >    & B1_l_duplicity,
								const eslocal MPIrank,
								const eslocal MPIsize)
{

	eslocal dofs_per_node = 3;
	bool flag_redund_lagr_mult = false;

	eslocal neighClustNum;	// number of neighboring sub-domains for current sub-domainG
	eslocal borderDofNum; 	// number of DOFs on surface on my cluster
	std::vector< esglobal >::iterator it_vec;

	std::vector < std::vector < esglobal > > neighBorderDofs;        			// 2D vector for DOFs on borders of neighboring sub-domains
	std::vector < esglobal > myBorderDOFs;  	// my border DOFs are here
	std::vector < eslocal  > myNeighClusters; 	// my neighboring clusters

	neighBorderDofs.resize( MPIsize, std::vector< esglobal >( 0 , 0 ) );

	MPI_Request   mpi_req;
	MPI_Status 	  mpi_stat;

	// Find all MPIranks of all neighboring clusters in _boundaries
	std::vector < eslocal > neigh_tmp  (MPIsize, 0);
	std::set<eslocal>::const_iterator it_set;

	for (T i = 0; i < _boundaries.size(); i++) {
		if ( _boundaries[i].size() > 1 ) {
			for (it_set = _boundaries[i].begin(); it_set != _boundaries[i].end(); ++it_set) {
				if (*it_set != MPIrank) { // if it does point non local cluster = points to one of neighboring clusters
					neigh_tmp[*it_set] = 1;
					for (int d_i = 0; d_i < dofs_per_node; d_i++ ) {
						neighBorderDofs[*it_set].push_back( dofs_per_node * _mesh.coordinates().globalIndex(i) + d_i ); // mapping local local to global
					}
				} else {
					for (int d_i = 0; d_i < dofs_per_node; d_i++ ) {
						myBorderDOFs.push_back( dofs_per_node * _mesh.coordinates().globalIndex(i) + d_i ); // mapping local local to global
					}
				}
	        }
		}
	}

	for (eslocal i = 0; i < neigh_tmp.size(); i++)
		if (neigh_tmp[i] == 1)
			myNeighClusters.push_back(i);

	neighClustNum = myNeighClusters.size();

	for (int i = 0; i < myNeighClusters.size(); i++) {
		neighBorderDofs[myNeighClusters[i]].swap(neighBorderDofs[i]);
	}
	neighBorderDofs.resize(neighClustNum);

	// END - Find all MPIranks of all neighboring clusters in _boundaries

	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }

	//neighBorderDofs.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );

	MPI_Request * mpi_send_req  = new MPI_Request [neighClustNum];
	MPI_Request * mpi_recv_req  = new MPI_Request [neighClustNum];
	MPI_Status  * mpi_recv_stat = new MPI_Status  [neighClustNum];


	// NOW :
	// my neighboring sub-domains are in : 	myNeighClusters
	// my border DOFs are in : 				myBorderDOFs
	// neighboring sub-domains border DOFs are in neighBorderDofs[i]

 	std::vector < std::vector < esglobal > > myNeighsSparse;
	myNeighsSparse.resize( myBorderDOFs.size() );

	for (eslocal j = 0; j < myBorderDOFs.size(); j++) {
		myNeighsSparse[j].reserve(5);
		myNeighsSparse[j].push_back(myBorderDOFs[j]);
	}

	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
		int d = myNeighClusters[i];
		int mydofs_index = 0;
		int nedofs_index = 0;

		if ( i == 0 && MPIrank < myNeighClusters[i] ) {
			for (eslocal j = 0; j < myBorderDOFs.size(); j++)
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


		if ( i < myNeighClusters.size() - 1)
			if ( MPIrank > myNeighClusters[i] && MPIrank < myNeighClusters[i+1] )
				for (eslocal j = 0; j < myBorderDOFs.size(); j++)
					myNeighsSparse[j].push_back(MPIrank);



		if ( i == myNeighClusters.size() -1 && MPIrank > myNeighClusters[i] )
			for (eslocal j = 0; j < myBorderDOFs.size(); j++)
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
					myLambdas[lambdaCount][0] = lambdaCount;	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
					myLambdas[lambdaCount][1] = dof;			// dof in global numbering
					myLambdas[lambdaCount][2] = (esglobal)MPIrank;		// my sub-domainG
					myLambdas[lambdaCount][3] = neighSD;		// neigh. sub-domainG
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

	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (MPIrank == 0) lambdaGlobalCount = 0;
	esglobal lambdaNum = lambdaCount + lambdaGlobalCount;
	MPI_Bcast(&lambdaNum, 1, MPI_LONG, MPIsize-1, MPI_COMM_WORLD);
	esglobal total_number_of_global_B1_lambdas = lambdaNum;

	if (myLambdas.size() > 0)
		//cilk_
		for (eslocal i = 0; i < myLambdas.size(); i++)
			myLambdas[i][0] = myLambdas[i][0]  + lambdaGlobalCount; //+ total_number_of_dirichlet_lambdas// create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index


	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }


	std::vector < std::vector < esglobal > > mpi_send_buff;
	mpi_send_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
	cilk_for (int i = 0; i < myNeighClusters.size(); i++) {
		int index = 0;
		if ( myLambdas.size() > 0 )
		{
			for (int j = 0; j < myLambdas.size(); j++) {
				if( myLambdas[j][3] == myNeighClusters[i] ) {
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

	for (int i = 0; i < myNeighClusters.size(); i++)
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), MPI_LONG, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);


	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }

	std::vector < std::vector < esglobal > > mpi_recv_buff;
	mpi_recv_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
	delete [] mpi_send_req;


	// receiving all border DOFs from all neighboring sub-domains
	eslocal messages_received = 0;
	while ( messages_received < myNeighClusters.size() ) {
		for (eslocal i = 0; i < myNeighClusters.size(); i++) {

			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
			int flag;
			MPI_Iprobe( myNeighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );

			if (flag) {
				int recv_msg_size = 0;

				MPI_Get_count(&mpi_stat, MPI_LONG, &recv_msg_size);
				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);

				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, MPI_LONG, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);

				free(mpi_tmp_recv_buff);
				messages_received++;
			}
		}
	}

	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }

	// decode received lambdas
	eslocal recv_lamba_count = 0;
	for (eslocal i = 0; i < myNeighClusters.size(); i++)
		if (mpi_recv_buff[i].size() > 1)
			recv_lamba_count += mpi_recv_buff[i].size();

	recv_lamba_count = recv_lamba_count / 5;

	eslocal l_i = myLambdas.size();
	myLambdas.resize( myLambdas.size() + recv_lamba_count, std::vector< esglobal >( 5 , 0 ) );

	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
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

	//bool comp_vf(const std::vector<esglobal> &a, const std::vector<esglobal> &b)
	//{
	//	return a[0] < b[0];
	//}

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

////#ifdef USE_TBB
////	tbb::mutex m;
////	cilk_for (int j = 0; j < n_myLambdas; j++)
////#else
//	for (int j = 0; j < n_myLambdas; j++)
////#endif
//	{
////		std::map< esglobal,int > :: iterator itm;
//
//		esglobal lambda    = myLambdas[j][0];
//		esglobal DOFNumber = myLambdas[j][1];
//		eslocal  cnt       = myLambdas[j][4];
//		double B_value;
//		if (myLambdas[j][2] < myLambdas[j][3])
//			B_value = 1.0;
//		else
//			B_value = -1.0;
//
//		for (int i = 0; i < data.size(); i++) {
//			itm = fem[i]->mesh.g2l.find(DOFNumber);		//( BJ[j] );
//			if ( itm != fem[i]->mesh.g2l.end() ) {
////			   #ifdef USE_TBB
////				m.lock();
////			   #endif
//				 data[i]->B->BI.push_back(lambda);		//( BI[j] );
//				 data[i]->B->BJ.push_back(itm->second);
//				 data[i]->B->BV.push_back(B_value);		//( BV[j] );
//				 myLambdas[j].push_back(i);
////			   #ifdef USE_TBB
////				m.unlock();
////			   #endif
//				break;
//			}
//		}
//	}
//
//
//	for (int i = 0; i < data.size(); i++){
//		data[i]->B->n_row_bg = lambdaNum;
//		data[i]->B->n_col_bg = domainG->neqSub[fem[i]->i_domOnClust];
//	}

	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }







//
//
//	std::vector < SparseDOKMatrix > B1_loc(MPIsize);
//
//	lambda_map_sub_B1.resize(MPIsize);
//
//	B1_l_duplicity.resize(MPIsize);
//	std::vector<eslocal> local_prim_numbering(MPIsize, 0);
//	B1_local.resize(MPIsize);
//
//	std::set<eslocal>::const_iterator it1;
//	std::set<eslocal>::const_iterator it2;
//
//	eslocal lambda_count_B1 = 0;
//
//


//	for (T i = 0; i < _boundaries.size(); i++) {
//
//		if ( _boundaries[i].size() > 1 ) {
//
////			// with duplicity
////			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
////				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
////					for (eslocal d_i = 0; d_i < 3; d_i++) {
////						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0;
////						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0;
////
////						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
////						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
////						lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );
////
////						B1_l_duplicity[*it1].push_back( 1.0 / (double)_boundaries[i].size() );
////						B1_l_duplicity[*it2].push_back( 1.0 / (double)_boundaries[i].size() );
////
////						lambda_count_B1++;
////					}
////				}
////			}
//
//
//			// no duplicity
//			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
//				if (it1 != _boundaries[i].begin()) {
//					for (eslocal d_i = 0; d_i < 3; d_i++) {
//						B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) =  1.0;
//						B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) = -1.0;
//
//						lambda_map_sub_B1[*it2].push_back(lambda_count_B1);
//						lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
//						lambda_map_sub_clst.push_back( std::vector < eslocal > (1,lambda_count_B1) );
//
//						B1_l_duplicity[*it1].push_back( 1.0 / (double)_boundaries[i].size() );
//						B1_l_duplicity[*it2].push_back( 1.0 / (double)_boundaries[i].size() );
//
//						lambda_count_B1++;
//
//						_mesh._coordinates.globalIndex()
//
//					}
//				}
//				it2 = it1;
//			}
//
//
//		}
//
//		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
//			local_prim_numbering[*it] += 3;
//        }
//	}
//
//	for (eslocal d = 0; d < MPIsize; d++) {
//		B1_loc[d].resize(lambda_count_B1, local_prim_numbering[d]);
//		B1_local[d] = B1_loc[d];
//	}
//
//        int a = 10;
//

}



}



