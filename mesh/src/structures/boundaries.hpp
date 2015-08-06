
#include "boundaries.h"
#include <mpi.h>

namespace mesh {


template<typename T>
void Boundaries::create_B1_l(	std::vector < SparseIJVMatrix<T> >         & B1_local,
								std::vector < SparseIJVMatrix<T> >         & B0_local,
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

	std::vector < SparseDOKMatrix<T> > B1_loc(domains_num);
	std::vector < SparseDOKMatrix<T> > B0_loc(domains_num);

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

		std::vector < bool > is_dirichlet (3, false); // TODO: 3 is number of DOFs per node

		for (it = _boundaries[i].begin(); it != _boundaries[i].end(); ++it) {
			if ( dirichlet_x.find(index(i)) != dirichlet_x.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 0) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				std::vector < eslocal > tmp_vec (2);
				tmp_vec[0] = lambda_count_B1;
				tmp_vec[1] = 0;
				lambda_map_sub_clst.push_back( tmp_vec );
				is_dirichlet[0] = true;
				B1_l_duplicity[*it].push_back( 1.0 );

//				if ( _boundaries[i].size() > 1 ) {
//					B1_l_duplicity[*it].push_back( 1.0 / ((double)_boundaries[i].size() + 1.0) );
//
// 				}
//				else {
//					B1_l_duplicity[*it].push_back( 1.0 );
//				}
				lambda_count_B1++;
			}
			if ( dirichlet_y.find(index(i)) != dirichlet_y.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 1) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				std::vector < eslocal > tmp_vec (2);
				tmp_vec[0] = lambda_count_B1;
				tmp_vec[1] = 0;
				lambda_map_sub_clst.push_back( tmp_vec );
				is_dirichlet[1] = true;
				B1_l_duplicity[*it].push_back( 1.0 );

//				if ( _boundaries[i].size() > 1 ) {
//					B1_l_duplicity[*it].push_back( 1.0 / ((double)_boundaries[i].size() + 1.0) );
//
// 				}
//				else {
//					B1_l_duplicity[*it].push_back( 1.0 );
//				}
				lambda_count_B1++;
			}
			if ( dirichlet_z.find(index(i)) != dirichlet_z.end() ) {
				B1_loc[*it](lambda_count_B1, local_prim_numbering[*it] + 2) =  1.0;  // 3*i + d_i
				lambda_map_sub_B1[*it].push_back(lambda_count_B1);
				std::vector < eslocal > tmp_vec (2);
				tmp_vec[0] = lambda_count_B1;
				tmp_vec[1] = 0;
				lambda_map_sub_clst.push_back( tmp_vec );
				is_dirichlet[2] = true;
				B1_l_duplicity[*it].push_back( 1.0 );

//				if ( _boundaries[i].size() > 1 ) {
//					B1_l_duplicity[*it].push_back( 1.0 / ((double)_boundaries[i].size() + 1.0) );
//
// 				}
//				else {
//					B1_l_duplicity[*it].push_back( 1.0 );
//				}
				lambda_count_B1++;
			}
		}

		if ( _boundaries[i].size() > 1 ) {

			// with duplicity
			for (it1 = _boundaries[i].begin(); it1 != _boundaries[i].end(); ++it1) {
				for (it2 = it1,++it2; it2 != _boundaries[i].end(); ++it2) {
					for (eslocal d_i = 0; d_i < 3; d_i++) { //TODO: 3 DOFS per ndoe
						if (!is_dirichlet[d_i]) {
							B1_loc[*it1](lambda_count_B1, local_prim_numbering[*it1] + d_i) =  1.0;
							B1_loc[*it2](lambda_count_B1, local_prim_numbering[*it2] + d_i) = -1.0;

							lambda_map_sub_B1[*it1].push_back(lambda_count_B1);
							lambda_map_sub_B1[*it2].push_back(lambda_count_B1);

							std::vector < eslocal > tmp_vec (2);
							tmp_vec[0] = lambda_count_B1;
							tmp_vec[1] = 0;
							lambda_map_sub_clst.push_back( tmp_vec );

							B1_l_duplicity[*it1].push_back( 1.0 / (double) _boundaries[i].size() );
							B1_l_duplicity[*it2].push_back( 1.0 / (double) _boundaries[i].size() );

							lambda_count_B1++;
						}
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



struct Comp_vf
{
	bool operator() (const std::vector<esglobal> &a, const std::vector<esglobal> &b)
	{
		return a[0] < b[0];
	}
};

//TODO: potrebujeme mapovani pro DOFs global(ne cluster) to local(subdomains)
template<typename T>
void Boundaries::create_B1_g(	std::vector < SparseIJVMatrix<T> >         & B1,
								const std::vector < SparseCSRMatrix<T> >   & K_mat,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
								std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
								std::vector < std::vector <double> >    & B1_duplicity,
								const eslocal MPIrank,
								const eslocal MPIsize,
								const eslocal subDomPerCluster,
								std::vector < eslocal  > & myNeighClusters,
								mesh::Boundaries & local_boundaries)
{


	// Local B1 - further processing - update row numbering based on all clusters
    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - start                                   "; system("date +%T.%6N"); }


	// Create lambda global numbering
	esglobal localB1_l_rows = B1[0].rows();
	//for (eslocal domain_index = 0; domain_index < subDomPerCluster; domain_index++)
	//	localB1_l_rows += B1[domain_index].rows();

	// Renumbering of the local B1 matrix including Dirichlet BC
	// Create lambda global counting for local B1
	esglobal global_B1_l_rows;
	MPI_Exscan(&localB1_l_rows, &global_B1_l_rows, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);
	if (MPIrank == 0)
		global_B1_l_rows = 0; //localB1_l_rows;

	global_B1_l_rows = localB1_l_rows + global_B1_l_rows;
	esglobal total_number_of_B1_l_rows = global_B1_l_rows;
	MPI_Bcast(&total_number_of_B1_l_rows, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);

    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - EXscan and Bcast done                  "; system("date +%T.%6N"); }


	cilk_for (eslocal domain_index=0; domain_index < subDomPerCluster; domain_index++) {
		//TODO: lambda muze byt esglobal ale IJV matice je jen eslocal
		esglobal row_offset = global_B1_l_rows - localB1_l_rows;
		B1[domain_index].ShiftRowIndex(row_offset);

		eslocal  cols_tmp = B1[domain_index].columns();
		B1[domain_index].resize(total_number_of_B1_l_rows, cols_tmp); // TODO: prvni je esglobal

		for (eslocal i = 0; i < lambda_map_sub_B1[domain_index].size(); i++)
			lambda_map_sub_B1[domain_index][i] += row_offset;

//		myLambdas_l[domain_index][i][0] += 	total_number_of_dirichlet_lambdas +
//											total_number_of_global_B1_lambdas +
//											lambdaGlobalCount_l;

	}

        if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - end of renumbering of rows of local B1   "; system("date +%T.%6N"); }


	esglobal row_offset = global_B1_l_rows - localB1_l_rows;
	for (eslocal i = 0; i < lambda_map_sub_clst.size(); i++) {
		lambda_map_sub_clst[i][0] += row_offset;
		lambda_map_sub_clst[i][1] = MPIrank;
	}

	// END - Local B1 - further processing - update row numbering based on all clusters


	std::vector < SparseDOKMatrix<T> > B1_DOK_tmp(subDomPerCluster);

	eslocal dofs_per_node = 3;
	bool flag_redund_lagr_mult = true;

	eslocal neighClustNum;	// number of neighboring sub-domains for current sub-domainG
	eslocal borderDofNum; 	// number of DOFs on surface on my cluster
	std::vector< esglobal >::iterator it_vec;

	std::vector < std::vector < esglobal > > neighBorderDofs;        			// 2D vector for DOFs on borders of neighboring sub-domains
	std::vector < esglobal > myBorderDOFs;  	// my border DOFs are here
	//std::vector < eslocal  > myNeighClusters; 	// my neighboring clusters

	neighBorderDofs.resize( MPIsize, std::vector< esglobal >( 0 , 0 ) );

	MPI_Request   mpi_req;
	MPI_Status 	  mpi_stat;

	// Find all MPIranks of all neighboring clusters in _boundaries
	std::vector < eslocal > neigh_tmp  (MPIsize, 0);
	std::set<eslocal>::const_iterator it_set;

    if (MPIrank == 0) { std::cout << " Global B - Blobal B1 neighdofs and neigh dofs array building             "; system("date +%T.%6N"); }


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

        if (MPIrank == 0) { std::cout << " Global B - neighDOFs array swapping based neigh cluster indexes          "; system("date +%T.%6N"); }


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
	mpi_send_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
#ifdef DEBUG
	for (int i = 0; i < myNeighClusters.size(); i++) {
#else
	cilk_for (int i = 0; i < myNeighClusters.size(); i++) {
#endif
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
		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);


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

				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);

				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
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

//#ifdef USE_TBB
//	tbb::mutex m;
//	cilk_for (int j = 0; j < n_myLambdas; j++)
//#else
	for (eslocal j = 0; j < n_myLambdas; j++)
//#endif
	{
		esglobal lambda         = myLambdas[j][0];
		esglobal DOFNumber      = myLambdas[j][1];

		//TODO: predpoklada 3 stupne volnosti na uzel - vychazi z mapovani g2l pro uzly a ne DOFy
		esglobal dofNODEnumber = DOFNumber / 3;
		eslocal  dofNODEoffset = DOFNumber % 3;
		eslocal  clustDofNODENumber = _mesh.coordinates().clusterIndex( dofNODEnumber );

		std::set < eslocal > & subs_with_element = _boundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel

		eslocal  cnt            = myLambdas[j][4];
		double B_value;
		if (myLambdas[j][2] < myLambdas[j][3])
			B_value = 1.0;
		else
			B_value = -1.0;

		//TODO:Tuto smycku prepsat, aby fungovala jen podoblasti, ve ktery je uzel
		for (eslocal d = 0; d < subDomPerCluster; d++) {
		//for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set) {
		//	eslocal d = *it_set;
			eslocal domDofNODENumber = _mesh.coordinates().localIndex(clustDofNODENumber, d);
			if ( domDofNODENumber != -1 ) {
				//TODO: predpoklada 3 stupne volnosti na uzel - vychazi z mapovani g2l pro uzly a ne DOFy
				eslocal domDOFNumber = 3 * domDofNODENumber + dofNODEoffset;
//			   #ifdef USE_TBB
//				m.lock();
//			   #endif
				//data[d]->B->BI.push_back(lambda);		//( BI[j] );
				//data[d]->B->BJ.push_back(domDOFNumber);
				//data[d]->B->BV.push_back(B_value);		//( BV[j] );
				B1_DOK_tmp[d](lambda, domDOFNumber) = B_value;
				myLambdas[j].push_back(d);
				B1_duplicity[d].push_back(1.0 / cnt);
//			   #ifdef USE_TBB
//				m.unlock();
//			   #endif
				break;
			}
		}
	}

	for (eslocal d = 0; d < subDomPerCluster; d++){
		//TODO: lambdaNum muze byt 64bit integer
		//TODO: matice B - pocet radku muze byt 64bit int
		B1_DOK_tmp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas , K_mat[d].rows());
		SparseIJVMatrix<T> ijv = B1_DOK_tmp[d];
		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
	}

	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }

	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }


	// for global B1
	for (int i = 0; i < myLambdas.size(); i++) {
		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
		tmp_vec[0] = myLambdas[i][0];
		tmp_vec[1] = myLambdas[i][2];
		tmp_vec[2] = myLambdas[i][3];

		lambda_map_sub_clst.push_back(tmp_vec);

		eslocal lam_tmp = myLambdas [i][0]; //TODO: esglobal
		lambda_map_sub_B1[ myLambdas[i][5] ].push_back( lam_tmp );
	}


	if (MPIrank == 0) { std::cout << " END - Creating lambda_map_sub vector of vectors - Global B1              "; system("date +%T.%6N"); }

}


}



