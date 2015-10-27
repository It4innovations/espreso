
#include "linear.h"

namespace physics {

template <>
size_t Linear<FEM>::subdomains()
{
	return this->_mesh.parts();
}

double determinant3x3(DenseMatrix &m)
{
	const double *values = m.values();
	return fabs(
		values[0] * values[4] * values[8] +
		values[1] * values[5] * values[6] +
		values[2] * values[3] * values[7] -
		values[2] * values[4] * values[6] -
		values[1] * values[3] * values[8] -
		values[0] * values[5] * values[7]
   );
}

void inverse(const DenseMatrix &m, DenseMatrix &inv, double det)
{
	const double *values = m.values();
	inv.resize(m.rows(), m.columns());
	double *invj = inv.values();
	double detJx = 1 / det;
	invj[0] = detJx * (values[8] * values[4] - values[7] * values[5]);
	invj[1] = detJx * (-values[8] * values[1] + values[7] * values[2]);
	invj[2] = detJx * (values[5] * values[1] - values[4] * values[2]);
	invj[3] = detJx * (-values[8] * values[3] + values[6] * values[5]);
	invj[4] = detJx * (values[8] * values[0] - values[6] * values[2]);
	invj[5] = detJx * (-values[5] * values[0] + values[3] * values[2]);
	invj[6] = detJx * (values[7] * values[3] - values[6] * values[4]);
	invj[7] = detJx * (-values[7] * values[0] + values[6] * values[1]);
	invj[8] = detJx * (values[4] * values[0] - values[3] * values[1]);
}

// B =
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
void distribute(DenseMatrix &B, DenseMatrix &dND)
{
	eslocal columns = dND.rows() * dND.columns();
	const double *dNDx = dND.values();
	const double *dNDy = dND.values() + dND.columns();
	const double *dNDz = dND.values() + 2 * dND.columns();

	double *v = B.values();

	memcpy(&v[0], dNDx,                               sizeof(double) * dND.columns());
	memcpy(&v[3 * columns + dND.columns()],     dNDx, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns + 2 * dND.columns()], dNDx, sizeof(double) * dND.columns());

	memcpy(&v[1 * columns + dND.columns()],     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[3 * columns],                     dNDy, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + 2 * dND.columns()], dNDy, sizeof(double) * dND.columns());

	memcpy(&v[2 * columns + 2 * dND.columns()], dNDz, sizeof(double) * dND.columns());
	memcpy(&v[4 * columns + dND.columns()],     dNDz, sizeof(double) * dND.columns());
	memcpy(&v[5 * columns],                     dNDz, sizeof(double) * dND.columns());
}


template <>
void Linear<FEM>::KeMefe(
		DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
		DenseMatrix &Ce, const mesh::Element *e, size_t part, bool dynamics)
{
	const std::vector<DenseMatrix> &dN = e->dN();
	const std::vector<DenseMatrix> &N = e->N();
	const std::vector<double> &weighFactor = e->weighFactor();
	std::vector<double> inertia;
	this->inertia(inertia);

	DenseMatrix coordinates(e->size(), mesh::Point::size());
	for (size_t i = 0; i < e->size(); i++) {
		coordinates.values() + i * mesh::Point::size() << _mesh.coordinates().get(e->node(i), part);
	}

	eslocal Ksize = mesh::Point::size() * this->DOFs();
	double detJ;
	DenseMatrix J, invJ, dND;

	Ke.resize(Ksize, Ksize);
	Ke = 0;
	fe.resize(Ksize);
	fill(fe.begin(), fe.end(), 0);
	if (dynamics) {
		Me.resize(e->size(), e->size());
		Me = 0;
	}

	for (eslocal gp = 0; gp < e->gpSize(); gp++) {
		J.multiply(dN[gp], coordinates);
		detJ = determinant3x3(J);
		inverse(J, invJ, detJ);

		dND.multiply(invJ, dN[gp]);

		// TODO: make it more general
		if (this->DOFs() == 3) {
			DenseMatrix B(Ce.rows(), Ksize);
			distribute(B, dND);
			Ke.multiply(B, Ce * B, detJ * weighFactor[gp], 1, true);
		} else {
			Ke.multiply(dND, Ce * dND, detJ * weighFactor[gp], 1, true);
		}

		for (eslocal i = 0; i < Ksize; i++) {
			fe[i] += detJ * weighFactor[gp] * N[gp](0, i % e->size()) * inertia[i / e->size()];
		}

		if (dynamics) {
			// Me = Me + WF * (DENS * dJ) * (N' * N);
			Me.multiply(N[gp], N[gp], this->rho() * detJ * weighFactor[gp] * this->CP(), 1, true);
		}
	}
}

template <>
void Linear<FEM>::integrate(
			DenseMatrix &Ke, DenseMatrix &Me, std::vector<double> &fe,
			SparseVVPMatrix<eslocal> &K, SparseVVPMatrix<eslocal> &M, std::vector<double> &f,
			const mesh::Element *e, bool dynamics)
{
	// Element ordering: xxxx, yyyy, zzzz,...
	// Global ordering:  xyz, xyz, xyz, xyz, ...
	size_t row, column;
	size_t s = this->DOFs();

	for (size_t i = 0; i < s * e->size(); i++) {
		row = s * (e->node(i % e->size())) + i / e->size();
		for (size_t j = 0; j < s * e->size(); j++) {
			column = s * (e->node(j % e->size())) + j / e->size();
			K(row, column) = Ke(i, j);
		}
		f[row] += fe[i];
	}
	if (!dynamics) {
		return;
	}
	for (size_t i = 0; i < e->size(); i++) {
		row = s * (e->node(i));
		for (size_t j = 0; j < e->size(); j++) {
			column = s * (e->node(j)); //i
			for (size_t k = 0; k < s; k++) {
				M(row + k, column + k) += Me(i, j);
			}
		}
	}
}


template <>
void Linear<FEM>::KMf(size_t part, bool dynamics)
{
	SparseVVPMatrix<eslocal> _K;
	SparseVVPMatrix<eslocal> _M;
	eslocal nK = _mesh.coordinates().localSize(part) * this->DOFs();
	_K.resize(nK, nK);
	if (dynamics) {
		_M.resize(nK, nK);
	}
	_f[part].resize(nK);

	DenseMatrix Ke, Me, Ce;
	std::vector<double> fe;

	this->C(Ce);

	const std::vector<eslocal> &partition = _mesh.getPartition();
	const std::vector<mesh::Element*> &elements = _mesh.getElements();
	for (eslocal i = partition[part]; i < partition[part + 1]; i++) {
		KeMefe(Ke, Me, fe, Ce, elements[i], part, dynamics);
		integrate(Ke, Me, fe, _K, _M, _f[part], elements[i], dynamics);
	}

	// TODO: make it direct
	SparseCSRMatrix<eslocal> csrK = _K;
	SparseCSRMatrix<eslocal> csrM = _M;
	this->_K[part] = csrK;
	this->_M[part] = csrM;
}

template <>
void Linear<FEM>::localB()
{
	const mesh::Boundaries &localBoundaries = this->_mesh.subdomainBoundaries();
	const mesh::Boundaries &globalBoundaries = this->_mesh.clusterBoundaries();

	std::vector<SparseDOKMatrix<eslocal> > gB(subdomains());
	std::vector<SparseDOKMatrix<eslocal> > lB(subdomains());

	std::vector<eslocal> local_prim_numbering(subdomains(), 0);

	std::vector<size_t>::const_iterator vi;
	std::set<eslocal>::const_iterator si1;
	std::set<eslocal>::const_iterator si2;

	eslocal lambda_count_B1 = 0;
	eslocal lambda_count_B0 = 0;

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

	for (size_t i = 0; i < localBoundaries.size(); i++) {

		std::vector<bool> is_dirichlet(this->DOFs(), false);
		for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
			for (vi = properties.begin(); vi != properties.end(); ++vi) {
				const std::map<eslocal, double> &property = this->_mesh.coordinates().property(*vi).values();
				if (property.find(i) != property.end()) {
					is_dirichlet[0] = true;

					gB[*si1](lambda_count_B1, local_prim_numbering[*si1] + (vi - properties.begin())) =  1.0;

					_lambda_map_sub_clst.push_back(std::vector<eslocal>({ lambda_count_B1, 0 }));
					_lambda_map_sub_B1[*si1].push_back(lambda_count_B1);
					_B1_duplicity[*si1].push_back(1.0);
					_vec_c[*si1].push_back(0.0);

					lambda_count_B1++;
				}
			}
		}

		if (localBoundaries[i].size() > 1) {
			if (globalBoundaries[i].size() == 1) {
				for (si1 = localBoundaries[i].begin(); si1 != localBoundaries[i].end(); ++si1) {
					for (si2 = si1,++si2; si2 != localBoundaries[i].end(); ++si2) {
						for (eslocal d = 0; d < this->DOFs(); d++) {
							if (is_dirichlet[d]) {
								continue;
							}
							gB[*si1](lambda_count_B1, local_prim_numbering[*si1] + d) =  1.0;
							gB[*si2](lambda_count_B1, local_prim_numbering[*si2] + d) = -1.0;

							_lambda_map_sub_B1[*si1].push_back(lambda_count_B1);
							_lambda_map_sub_B1[*si2].push_back(lambda_count_B1);
							_lambda_map_sub_clst.push_back(std::vector<eslocal>({ lambda_count_B1, 0 }));
							_B1_duplicity[*si1].push_back( 1.0 / (double) localBoundaries[i].size() );
							_B1_duplicity[*si2].push_back( 1.0 / (double) localBoundaries[i].size() );
							_vec_c[*si1].push_back(0.0);
							_vec_c[*si2].push_back(0.0);

							lambda_count_B1++;
						}
					}
				}
			}
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

	// TODO: make it direct
	for (eslocal d = 0; d < subdomains(); d++) {
		gB[d].resize(lambda_count_B1, local_prim_numbering[d]);
		SparseIJVMatrix<eslocal> ijvGB = gB[d];
		_globalB[d] = ijvGB;

		lB[d].resize(lambda_count_B0, local_prim_numbering[d]);
		SparseIJVMatrix<eslocal> ijvLB = lB[d];
		_localB[d] = ijvLB;
	}
}

template <>
void Linear<FEM>::globalB()
{
//	// Local B1 - further processing - update row numbering based on all clusters
//    if (this->_verbose && this->_mesh.rank() == 0) {
//    	std::cout << " Global B - Local preprocessing - start                                   "; system("date +%T.%6N");
//    }
//
//
//	// Create lambda global numbering
//	esglobal localB1_l_rows = B1[0].rows();
//	//for (eslocal domain_index = 0; domain_index < subDomPerCluster; domain_index++)
//	//	localB1_l_rows += B1[domain_index].rows();
//
//	// Renumbering of the local B1 matrix including Dirichlet BC
//	// Create lambda global counting for local B1
//	esglobal global_B1_l_rows;
//	MPI_Exscan(&localB1_l_rows, &global_B1_l_rows, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);
//	if (MPIrank == 0)
//		global_B1_l_rows = 0; //localB1_l_rows;
//
//	global_B1_l_rows = localB1_l_rows + global_B1_l_rows;
//	esglobal total_number_of_B1_l_rows = global_B1_l_rows;
//	MPI_Bcast(&total_number_of_B1_l_rows, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
//
//    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - EXscan and Bcast done                  "; system("date +%T.%6N"); }
//
//
//	cilk_for (eslocal domain_index=0; domain_index < subDomPerCluster; domain_index++) {
//		//TODO: lambda muze byt esglobal ale IJV matice je jen eslocal
//		esglobal row_offset = global_B1_l_rows - localB1_l_rows;
//		B1[domain_index].ShiftRowIndex(row_offset);
//
//		eslocal  cols_tmp = B1[domain_index].columns();
//		B1[domain_index].resize(total_number_of_B1_l_rows, cols_tmp); // TODO: prvni je esglobal
//
//		for (eslocal i = 0; i < lambda_map_sub_B1[domain_index].size(); i++)
//			lambda_map_sub_B1[domain_index][i] += row_offset;
//
////		myLambdas_l[domain_index][i][0] += 	total_number_of_dirichlet_lambdas +
////											total_number_of_global_B1_lambdas +
////											lambdaGlobalCount_l;
//
//	}
//
//    if (MPIrank == 0) { std::cout << " Global B - Local preprocessing - end of renumbering of rows of local B1   "; system("date +%T.%6N"); }
//
//
//	esglobal row_offset = global_B1_l_rows - localB1_l_rows;
//	for (eslocal i = 0; i < lambda_map_sub_clst.size(); i++) {
//		lambda_map_sub_clst[i][0] += row_offset;
//		lambda_map_sub_clst[i][1] = MPIrank;
//	}
//
//	// END - Local B1 - further processing - update row numbering based on all clusters
//
//
//	std::vector < SparseDOKMatrix<T> > B1_DOK_tmp(subDomPerCluster);
//
//	eslocal dofs_per_node = DOFS_PER_NODE;
//	bool flag_redund_lagr_mult = true;
//
//	eslocal neighClustNum;	// number of neighboring sub-domains for current sub-domainG
//	eslocal borderDofNum; 	// number of DOFs on surface on my cluster
//	std::vector< esglobal >::iterator it_vec;
//
//	std::vector < std::vector < esglobal > > neighBorderDofs;        			// 2D vector for DOFs on borders of neighboring sub-domains
//	std::vector < esglobal > myBorderDOFs;  	// my border DOFs are here
//	//std::vector < eslocal  > myNeighClusters; 	// my neighboring clusters
//
//
//	MPI_Request   mpi_req;
//	MPI_Status 	  mpi_stat;
//
//	// Find all MPIranks of all neighboring clusters in _boundaries
//	std::vector < eslocal > neigh_tmp  (MPIsize, 0);
//	std::set<eslocal>::const_iterator it_set;
//	std::set<eslocal>::const_iterator it_set_l;
//
//    if (MPIrank == 0) { std::cout << " Global B - Blobal B1 neighdofs and neigh dofs array building             "; system("date +%T.%6N"); }
//
//    // Now this loop is used to get my Border DOFs and my neighboring subdomains
//    // information about neighnoring DOFS is here, but it is not used
//	// neighBorderDofs.resize( MPIsize, std::vector< esglobal >( 0 , 0 ) );
//	for (T i = 0; i < _boundaries.size(); i++) {
//		if ( _boundaries[i].size() > 1 && local_boundaries[i].size() == 1 ) {
//			for (it_set = _boundaries[i].begin(); it_set != _boundaries[i].end(); ++it_set) {
//				if (*it_set != MPIrank) { // if it does point non local cluster = points to one of neighboring clusters
//					//neigh_tmp[*it_set] = 1;
//					// TODO: If the neighborring nodes are assembled localy - info from Mesh - use next 3 lines - but it does not have information on duplicity
//					//for (int d_i = 0; d_i < dofs_per_node; d_i++ ) {
//					//	neighBorderDofs[*it_set].push_back( dofs_per_node * coordinates.globalIndex(i) + d_i ); // mapping local local to global
//					//}
//				} else {
//					// this loop goes over the nodes that are in multiple domains on this one cluster
//					//for (it_set_l = local_boundaries[i].begin(); it_set_l != local_boundaries[i].end(); ++it_set_l) {
//						for (int d_i = 0; d_i < dofs_per_node; d_i++ ) {
//							myBorderDOFs.push_back( dofs_per_node * coordinates.globalIndex(i) + d_i ); // mapping local local to global
//						}
//					//}
//				}
//	        }
//		}
//	}
//
//
//	for (T i = 0; i < _boundaries.size(); i++) {
//		if ( _boundaries[i].size() > 1 && local_boundaries[i].size() == 1 ) {
//			for (it_set = _boundaries[i].begin(); it_set != _boundaries[i].end(); ++it_set) {
//				if (*it_set != MPIrank) { // if it does point non local cluster = points to one of neighboring clusters
//					neigh_tmp[*it_set] = 1;
//				}
//			}
//		}
//	}
//
//	for (eslocal i = 0; i < neigh_tmp.size(); i++)
//		if (neigh_tmp[i] == 1)
//			myNeighClusters.push_back(i);
//
//	neighClustNum = myNeighClusters.size();
//
//
//    if (MPIrank == 0) { std::cout << " Global B - myNeighDOFs arrays are transfered to neighbors                "; system("date +%T.%6N"); }
//
//    MPI_Request * mpi_send_req  = new MPI_Request [neighClustNum];
//	MPI_Request * mpi_recv_req  = new MPI_Request [neighClustNum];
//	MPI_Status  * mpi_recv_stat = new MPI_Status  [neighClustNum];
//
//	neighBorderDofs.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );
//
//
//	neighBorderDofs.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );
//
//	// sending my DOFs on border to all neighboring sub-domains
//	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//		MPI_Isend(&myBorderDOFs[0], myBorderDOFs.size(), esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
//	}
//
//	// receiving all border DOFs from all neighboring sub-domains
//	eslocal messages_received_bd = 0;
//	while ( messages_received_bd < myNeighClusters.size() ) {
//		for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//
//			int flag;
//			MPI_Iprobe( myNeighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );
//
//			if (flag) {
//				int recv_msg_size = 0;
//
//				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
//
//				neighBorderDofs[i].resize(recv_msg_size);
//				if (recv_msg_size < 0)
//					std::cout << "error in msg size !";
//				MPI_Recv(&neighBorderDofs[i][0], recv_msg_size, esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
//
//				messages_received_bd++;
//			}
//		}
//	}
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//
//	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }
//
//
//
// 	std::vector < std::vector < esglobal > > myNeighsSparse;
//	myNeighsSparse.resize( myBorderDOFs.size() );
//
//	for (eslocal j = 0; j < myBorderDOFs.size(); j++) {
//		myNeighsSparse[j].reserve(5);
//		myNeighsSparse[j].push_back(myBorderDOFs[j]);
//
//	}
//
//	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//		int d = myNeighClusters[i];
//		int mydofs_index = 0;
//		int nedofs_index = 0;
//
//		if ( i == 0 && MPIrank < myNeighClusters[i] ) {
//			for (eslocal j = 0; j < myBorderDOFs.size(); j++)
//				myNeighsSparse[j].push_back(MPIrank);
//		}
//
//		do {
//
//			if ( neighBorderDofs[i][nedofs_index] == myBorderDOFs[mydofs_index] ) {
//				myNeighsSparse[mydofs_index].push_back(d);
//				mydofs_index++;
//				nedofs_index++;
//			} else {
//				if ( neighBorderDofs[i][nedofs_index] > myBorderDOFs[mydofs_index] ) {
//					mydofs_index++;
//				} else {
//					nedofs_index++;
//				}
//			}
//
//		} while (mydofs_index < myBorderDOFs.size() && nedofs_index < neighBorderDofs[i].size() );
//
//
//		if ( i < myNeighClusters.size() - 1)
//			if ( MPIrank > myNeighClusters[i] && MPIrank < myNeighClusters[i+1] )
//				for (eslocal j = 0; j < myBorderDOFs.size(); j++)
//					myNeighsSparse[j].push_back(MPIrank);
//
//
//
//		if ( i == myNeighClusters.size() -1 && MPIrank > myNeighClusters[i] )
//			for (eslocal j = 0; j < myBorderDOFs.size(); j++)
//				myNeighsSparse[j].push_back(MPIrank);
//
//	}
//
//	for (eslocal i = 0; i < myNeighsSparse.size(); i++)
//		if (myNeighsSparse[i].size() < 3)
//			myNeighsSparse[i].clear();
//
//	if (MPIrank == 0) { std::cout << " Global B - myNeighSparse assembled                                       "; system("date +%T.%6N"); }
//
//	std::vector < std::vector < esglobal > > myLambdas;
//
//	esglobal lambdaCount = 0;
//
//	for (eslocal i = 0; i < myNeighsSparse.size(); i++) {
//		bool add = false;
//		if (myNeighsSparse[i].size() > 0) {
//			esglobal dof = myNeighsSparse[i][0];
//
//			eslocal cnt_tmp = myNeighsSparse[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb
//
//			for (eslocal j = 1; j < myNeighsSparse[i].size(); j++) {
//				esglobal neighSD = myNeighsSparse[i][j];
//
//				if ( add ) {
//
//					myLambdas.push_back ( std::vector < esglobal > () );
//					myLambdas[lambdaCount].reserve(6);
//					myLambdas[lambdaCount].resize(5);
//					myLambdas[lambdaCount][0] = lambdaCount;	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
//					myLambdas[lambdaCount][1] = dof;			// dof in global numbering
//					myLambdas[lambdaCount][2] = (esglobal)MPIrank;		// my sub-domainG
//					myLambdas[lambdaCount][3] = neighSD;		// neigh. sub-domainG
//					myLambdas[lambdaCount][4] = (esglobal)cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb
//					lambdaCount++;
//
//					if (!flag_redund_lagr_mult)
//						break;
//				}
//
//				if (neighSD == (esglobal)MPIrank)
//					add = true;
//			}
//		}
//	}
//
//	if (MPIrank == 0) { std::cout << " Global B - Create global lambda numbering                                "; system("date +%T.%6N"); }
//
//	esglobal lambdaGlobalCount = 0;
//
//	MPI_Exscan(&lambdaCount, &lambdaGlobalCount, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);
//
//	if (MPIrank == 0) lambdaGlobalCount = 0;
//	esglobal lambdaNum = lambdaCount + lambdaGlobalCount;
//	MPI_Bcast(&lambdaNum, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
//	esglobal total_number_of_global_B1_lambdas = lambdaNum;
//
//	if (myLambdas.size() > 0)
//		//cilk_
//		for (eslocal i = 0; i < myLambdas.size(); i++)
//			myLambdas[i][0] = myLambdas[i][0]  + lambdaGlobalCount + total_number_of_B1_l_rows; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index
//
//
//	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }
//
//
//	std::vector < std::vector < esglobal > > mpi_send_buff;
//	mpi_send_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
////#ifdef DEBUG
////	for (int i = 0; i < myNeighClusters.size(); i++) {
////#else
//	cilk_for (int i = 0; i < myNeighClusters.size(); i++) {
////#endif
//		int index = 0;
//		if ( myLambdas.size() > 0 )
//		{
//			for (int j = 0; j < myLambdas.size(); j++) {
//				if( myLambdas[j][3] == myNeighClusters[i] ) {
//					mpi_send_buff[i].push_back(myLambdas[j][0]);
//					mpi_send_buff[i].push_back(myLambdas[j][1]);
//					mpi_send_buff[i].push_back(myLambdas[j][2]);
//					mpi_send_buff[i].push_back(myLambdas[j][3]);
//					mpi_send_buff[i].push_back(myLambdas[j][4]);
//					index++;
//				}
//			}
//			if (index == 0)
//				mpi_send_buff[i].push_back(0);
//		}
//		else
//		{
//			mpi_send_buff[i].push_back(0);
//		}
//	}
//
//	if (MPIrank == 0) { std::cout << " Global B - Isend                                                         "; system("date +%T.%6N"); }
//
//	for (int i = 0; i < myNeighClusters.size(); i++)
//		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
//
//
//	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }
//
//	std::vector < std::vector < esglobal > > mpi_recv_buff;
//	mpi_recv_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
//	//delete [] mpi_send_req;
//
//
//	// receiving all border DOFs from all neighboring sub-domains
//	eslocal messages_received = 0;
//	while ( messages_received < myNeighClusters.size() ) {
//		for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//
//			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
//			int flag;
//			MPI_Iprobe( myNeighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );
//
//			if (flag) {
//				int recv_msg_size = 0;
//
//				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
//				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);
//
//				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
//				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);
//
//				free(mpi_tmp_recv_buff);
//				messages_received++;
//			}
//		}
//	}
//
//	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }
//
//	// decode received lambdas
//	eslocal recv_lamba_count = 0;
//	for (eslocal i = 0; i < myNeighClusters.size(); i++)
//		if (mpi_recv_buff[i].size() > 1)
//			recv_lamba_count += mpi_recv_buff[i].size();
//
//	recv_lamba_count = recv_lamba_count / 5;
//
//	eslocal l_i = myLambdas.size();
//	myLambdas.resize( myLambdas.size() + recv_lamba_count, std::vector< esglobal >( 5 , 0 ) );
//
//	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//		if (mpi_recv_buff[i].size() > 1) {
//			for (int j = 0; j < mpi_recv_buff[i].size() / 5; j++) {
//				myLambdas[l_i][0] = mpi_recv_buff[i][5*j + 0];
//				myLambdas[l_i][1] = mpi_recv_buff[i][5*j + 1];
//				myLambdas[l_i][2] = mpi_recv_buff[i][5*j + 3];
//				myLambdas[l_i][3] = mpi_recv_buff[i][5*j + 2];
//				myLambdas[l_i][4] = mpi_recv_buff[i][5*j + 4];
//				l_i++;
//			}
//		}
//	}
//
//	if (MPIrank == 0) { std::cout << " Global B - myLambdas - sort or tbb:sort                                  "; system("date +%T.%6N"); }
//
//	//bool comp_vf(const std::vector<esglobal> &a, const std::vector<esglobal> &b)
//	//{
//	//	return a[0] < b[0];
//	//}
//
//	auto comp_vf = [](const std::vector<esglobal> &a, const std::vector<esglobal> &b) {return a[0] < b[0];};
//
////#ifdef USE_TBB
////	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
////#else
//	std::sort         (myLambdas.begin(), myLambdas.end(), comp_vf);
////#endif
//
//	if (MPIrank == 0) { std::cout << " Global B - Final B assembling with g2l mapping using std::map            "; system("date +%T.%6N"); }
//
//	esglobal lambda;
//	esglobal DOFNumber;
//	eslocal n_myLambdas = myLambdas.size();
//
////#ifdef USE_TBB
////	tbb::mutex m;
////	cilk_for (int j = 0; j < n_myLambdas; j++)
////#else
//	for (eslocal j = 0; j < n_myLambdas; j++)
////#endif
//	{
//
//		esglobal lambda         = myLambdas[j][0];
//		esglobal DOFNumber      = myLambdas[j][1];
//
//		esglobal dofNODEnumber = DOFNumber / DOFS_PER_NODE;
//		eslocal  dofNODEoffset = DOFNumber % DOFS_PER_NODE;
//		eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );
//
//		const std::set < eslocal > & subs_with_element = _boundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel
//
//		eslocal  cnt            = myLambdas[j][4];
//		double B_value;
//		if (myLambdas[j][2] < myLambdas[j][3])
//			B_value = 1.0;
//		else
//			B_value = -1.0;
//
//		//TODO:Tuto smycku prepsat, aby fungovala jen podoblasti, ve ktery je uzel
//		for (eslocal d = 0; d < subDomPerCluster; d++) {
//		//for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set) {
//		//	eslocal d = *it_set;
//			eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
//			if ( domDofNODENumber != -1 ) {
//				//vychazi z mapovani g2l pro uzly a ne DOFy
//				eslocal domDOFNumber = DOFS_PER_NODE * domDofNODENumber + dofNODEoffset;
////			   #ifdef USE_TBB
////				m.lock();
////			   #endif
//				B1_DOK_tmp[d](lambda, domDOFNumber) = B_value;
//				myLambdas[j].push_back(d);
//				B1_duplicity[d].push_back(1.0 / cnt);
//				vec_c[d].push_back(0.0);
////			   #ifdef USE_TBB
////				m.unlock();
////			   #endif
//				break;
//
//			}
//		}
//	}
//
//	for (eslocal d = 0; d < subDomPerCluster; d++){
//		//TODO: lambdaNum muze byt 64bit integer
//		//TODO: matice B - pocet radku muze byt 64bit int
//		B1_DOK_tmp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas , K_mat[d].rows());
//		SparseIJVMatrix<T> ijv = B1_DOK_tmp[d];
//		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
//	}
//
//	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }
//
//	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }
//
//
//	// for global B1
//	for (int i = 0; i < myLambdas.size(); i++) {
//		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
//		tmp_vec[0] = myLambdas[i][0];
//		tmp_vec[1] = myLambdas[i][2];
//		tmp_vec[2] = myLambdas[i][3];
//
//		lambda_map_sub_clst.push_back(tmp_vec);
//
//		eslocal lam_tmp = myLambdas [i][0]; //TODO: esglobal
//		lambda_map_sub_B1[ myLambdas[i][5] ].push_back( lam_tmp );
//	}
//
//
//	if (MPIrank == 0) { std::cout << " END - Creating lambda_map_sub vector of vectors - Global B1              "; system("date +%T.%6N"); }
//
//	if (MPIrank == 0) { std::cout << " Dual size: " <<  total_number_of_B1_l_rows + total_number_of_global_B1_lambdas  << std::endl; }
//
//        MPI_Barrier(MPI_COMM_WORLD);
//
//
//    	flag_redund_lagr_mult = true;
//    	eslocal neighClustNum_sp = 0;	// number of neighboring sub-domains for current sub-domainG
//    	eslocal borderDofNum_sp  = 0; 	// number of DOFs on surface on my cluster
//
//    	std::vector < std::vector < esglobal > > neighBorderDofs_sp;    // 2D vector for DOFs on borders of neighboring sub-domains
//    	std::vector < esglobal > 				 myBorderDOFs_sp;  	// my border DOFs are here - in global numbering
//        std::vector < esglobal > 				 myBorderDOFs_sp_nr;  	// my border DOFs are here - in global numbering
//        std::vector < esglobal > 				 myBorderDOFs_sp_loc_n; // my border DOFs are here - in local numbering
//
//    	//std::vector < eslocal  > myNeighClusters; 	// my neighboring clusters
//
//
//    	//MPI_Request   mpi_req;
//    	//MPI_Status 	  mpi_stat;
//
//    	// Find all MPIranks of all neighboring clusters in _boundaries
//    	//std::vector < eslocal > neigh_tmp_sp  (MPIsize, 0);
//    	//std::set<eslocal>::const_iterator it_set;
//    	//std::set<eslocal>::const_iterator it_set_l;
//
//        if (MPIrank == 0) { std::cout << " Global B SP - Blobal B1 neighdofs and neigh dofs array building          "; system("date +%T.%6N"); }
//
//        // Now this loop is used to get my Border DOFs and my neighboring subdomains
//        // information about neighnoring DOFS is here, but it is not used
//    	// neighBorderDofs.resize( MPIsize, std::vector< esglobal >( 0 , 0 ) );
//    	for (T i = 0; i < _boundaries.size(); i++) {
//    		if ( _boundaries[i].size() > 1 && local_boundaries[i].size() > 1 ) {
//    			for (it_set = _boundaries[i].begin(); it_set != _boundaries[i].end(); ++it_set) {
//    				if (*it_set == MPIrank) {
//    					// this loop goes over the nodes that are in multiple domains on this one cluster
//    					for (it_set_l = local_boundaries[i].begin(); it_set_l != local_boundaries[i].end(); ++it_set_l) {
//    						for (int d_i = 0; d_i < dofs_per_node; d_i++ ) {
//    							myBorderDOFs_sp      .push_back( dofs_per_node * coordinates.globalIndex(i) + d_i ); // mapping local local to global
//    							myBorderDOFs_sp_loc_n.push_back( dofs_per_node *                          i + d_i ); // in local numbering
//    						}
//    					}
//    				}
//    	        }
//    		}
//    	}
//
//    	// sort my neigh DOFs sp
//    	std::sort (myBorderDOFs_sp.begin(), myBorderDOFs_sp.end());
//        myBorderDOFs_sp_nr = myBorderDOFs_sp;
//    	// removes the duplicit DOFs from myBorderDofs
//    	std::vector< esglobal >::iterator itx;
//        itx = std::unique (myBorderDOFs_sp_nr.begin(), myBorderDOFs_sp_nr.end());
//        myBorderDOFs_sp_nr.resize( std::distance(myBorderDOFs_sp_nr.begin(), itx) );
//
//    	neighClustNum = myNeighClusters.size();
//
//
//        if (MPIrank == 0) { std::cout << " Global B - myNeighDOFs arrays are transfered to neighbors                "; system("date +%T.%6N"); }
//
//
//       	neighBorderDofs_sp.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );
//
//    	// sending my DOFs on border to all neighboring sub-domains
//    	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//    		MPI_Isend(&myBorderDOFs_sp[0], myBorderDOFs_sp.size(), esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
//    	}
//
//    	// receiving all border DOFs from all neighboring sub-domains
//    	eslocal messages_received_spp = 0;
//    	while ( messages_received_spp < myNeighClusters.size() ) {
//    		for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//
//    			int flag;
//    			MPI_Iprobe( myNeighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );
//
//    			if (flag) {
//    				int recv_msg_size = 0;
//
//    				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
//
//    				neighBorderDofs_sp[i].resize(recv_msg_size);
//
//    				MPI_Recv(&neighBorderDofs_sp[i][0], recv_msg_size,esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
//
//    				messages_received_spp++;
//    			}
//    		}
//    	}
//
//    	MPI_Barrier(MPI_COMM_WORLD);
//
//    	if (MPIrank == 0) { std::cout << " Global B - Local preprocessing done                                      "; system("date +%T.%6N"); }
//
//
//     	std::vector < std::vector < esglobal > > myNeighsSparse_sp;
//    	myNeighsSparse_sp.resize( myBorderDOFs_sp_nr.size() );
//
//    	for (eslocal j = 0; j < myBorderDOFs_sp_nr.size(); j++) {
//   			myNeighsSparse_sp[j].reserve(5);
//   			myNeighsSparse_sp[j].push_back(myBorderDOFs_sp_nr[j]);
//    	}
//
//
//    	//////////////////
//
//    	std::vector < std::vector < eslocal  > > neighBorderDofs_sp_cnt;    	// number of duplicity in each element
//    	std::vector < std::vector < esglobal > > neighBorderDofs_sp_red;    	// neigh DOFs wo duplicity
//
//    	neighBorderDofs_sp_cnt.resize( neighClustNum, std::vector< eslocal >( 0 , 0 ) );
//    	neighBorderDofs_sp_red.resize( neighClustNum, std::vector< esglobal >( 0 , 0 ) );
//
//    	std::vector               < eslocal  > 	myBorderDOFs_sp_cnt;    	// the same but for neigh DOFs
//    	std::vector               < esglobal > 	myBorderDOFs_sp_red;    	// my border dofs w/o duplicity
//
//       	eslocal dup_cnt_sp = 1;
//       	for (eslocal i = 1; i < myBorderDOFs_sp.size(); i++) {
//
//       		if (myBorderDOFs_sp[i-1] != myBorderDOFs_sp[i]) {
//       			myBorderDOFs_sp_cnt.push_back(dup_cnt_sp);
//       			myBorderDOFs_sp_red.push_back(myBorderDOFs_sp[i-1]);
//       			dup_cnt_sp = 1;
//       		} else {
//       			dup_cnt_sp++;
//       		}
//
//       	}
//
//       	if (dup_cnt_sp > 1) {
//   			myBorderDOFs_sp_cnt.push_back(dup_cnt_sp);
//   			myBorderDOFs_sp_red.push_back(myBorderDOFs_sp[myBorderDOFs_sp.size() - 1]);
//       	}
//
//
//    	for (eslocal j = 0; j < myNeighClusters.size(); j++) {
//
//    		//if ( neighBorderDofs_sp[j].size() && myBorderDOFs_sp.size() ) {
//
//				dup_cnt_sp = 1;
//				for (eslocal i = 1; i < neighBorderDofs_sp[j].size(); i++) {
//
//					if (neighBorderDofs_sp[j][i-1] != neighBorderDofs_sp[j][i]) {
//						neighBorderDofs_sp_cnt[j].push_back(dup_cnt_sp);
//						neighBorderDofs_sp_red[j].push_back(neighBorderDofs_sp[j][i-1]);
//						dup_cnt_sp = 1;
//					} else {
//						dup_cnt_sp++;
//					}
//
//				}
//
//				if (dup_cnt_sp > 1 ){
//					neighBorderDofs_sp_cnt[j].push_back(dup_cnt_sp);
//					//neighBorderDofs_sp_red[j].push_back(myBorderDOFs_sp[neighBorderDofs_sp_red[j].size() - 1]);
//                                                  neighBorderDofs_sp_red[j].push_back(neighBorderDofs_sp[j][neighBorderDofs_sp[j].size() - 1]);
//				}
//
//    		//}
//    	}
//
//
//    	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//
//    		//if ( neighBorderDofs_sp[i].size() && myBorderDOFs_sp.size() ) {
//
//				int d = myNeighClusters[i];
//				int mydofs_index = 0;
//				int nedofs_index = 0;
//
//				if ( i == 0 && MPIrank < myNeighClusters[i] ) {
//					for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
//						for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
//							myNeighsSparse_sp[j].push_back(MPIrank);
//				}
//
//				do {
//				  if ( neighBorderDofs_sp_red[i].size() > 0 && myBorderDOFs_sp_red.size() > 0 )
//					if ( neighBorderDofs_sp_red[i][nedofs_index] == myBorderDOFs_sp_red[mydofs_index] ) {
//						for (eslocal j = 0; j < neighBorderDofs_sp_cnt[i][nedofs_index]; j++)
//							myNeighsSparse_sp[mydofs_index].push_back(d);
//						mydofs_index++;
//						nedofs_index++;
//					} else {
//						if ( neighBorderDofs_sp_red[i][nedofs_index] > myBorderDOFs_sp_red[mydofs_index] ) {
//							mydofs_index++;
//						} else {
//							nedofs_index++;
//						}
//					}
//
//				} while (mydofs_index < myBorderDOFs_sp_red.size() && nedofs_index < neighBorderDofs_sp_red[i].size() );
//
//
//				if ( i < myNeighClusters.size() - 1)
//					if ( MPIrank > myNeighClusters[i] && MPIrank < myNeighClusters[i+1] )
//						for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
//							for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
//								myNeighsSparse_sp[j].push_back(MPIrank);
//
//
//
//				if ( i == myNeighClusters.size() -1 && MPIrank > myNeighClusters[i] )
//					for (eslocal j = 0; j < myBorderDOFs_sp_red.size(); j++)
//						for (eslocal k = 0; k < myBorderDOFs_sp_cnt[j]; k++)
//							myNeighsSparse_sp[j].push_back(MPIrank);
//    		//}
//    	}
//
//
//
//    	for (eslocal i = 0; i < myNeighsSparse_sp.size(); i++)
//    		if (myNeighsSparse_sp[i].size() < 3)
//    			myNeighsSparse_sp[i].clear();
//
//
//    	if (MPIrank == 0) { std::cout << " Global B - myNeighSparse assembled                                       "; system("date +%T.%6N"); }
//
//    	std::vector < std::vector < esglobal > > myLambdas_sp;
//
//    	esglobal lambdaCount_sp = 0;
//
//    	for (eslocal i = 0; i < myNeighsSparse_sp.size(); i++) {
//    		bool add = false;
//    		if (myNeighsSparse_sp[i].size() > 0) {
//
//    			esglobal dof = myNeighsSparse_sp[i][0];
//    			eslocal  cnt_tmp = myNeighsSparse_sp[i].size() - 1; // pocet uzlu ucastnicich se duplicitnich vazeb
//
//    			int min_index = 1;
//    			int max_index = myNeighsSparse_sp[i].size();
//
//    			for (eslocal i_t = 1; i_t < myNeighsSparse_sp[i].size(); i_t++) {
//    				if ( myNeighsSparse_sp[i][i_t] == esglobal(MPIrank) ) {
//    					min_index = i_t;
//    					break;
//    				}
//    			}
//
//    			for (eslocal i_t = 1; i_t < myNeighsSparse_sp[i].size(); i_t++) {
//    				if ( myNeighsSparse_sp[i][i_t] > esglobal(MPIrank) ) {
//    					max_index = i_t;
//    					break;
//    				}
//    			}
//
//
//    			for (eslocal k = min_index; k < max_index; k++) {
//    				int tmp_cl = k-min_index;
//					for (eslocal j = k + 1; j < myNeighsSparse_sp[i].size(); j++) {
//
//						tmp_cl++;
//						if (myNeighsSparse_sp[i][j] != myNeighsSparse_sp[i][j-1])
//							tmp_cl = 0;
//
//						esglobal neighSD = myNeighsSparse_sp[i][j];
//
//						//if ( add ) {
//
//							myLambdas_sp.push_back ( std::vector < esglobal > () );
//							myLambdas_sp[lambdaCount_sp].reserve(7);
//							myLambdas_sp[lambdaCount_sp].resize(7);
//							myLambdas_sp[lambdaCount_sp][0] = lambdaCount_sp;	// at this point local index of this lambda - needs to be updated after MPIgather and MPIscatter
//							myLambdas_sp[lambdaCount_sp][1] = dof;			// dof in global numbering
//							myLambdas_sp[lambdaCount_sp][2] = (esglobal)MPIrank;		// my sub-domainG
//							myLambdas_sp[lambdaCount_sp][3] = neighSD;		// neigh. sub-domainG
//							myLambdas_sp[lambdaCount_sp][4] = (esglobal)cnt_tmp;		// delitel lambdy pri zapisu do B matice - odpovida poctu duplicitnich vazeb
//							myLambdas_sp[lambdaCount_sp][5] = (esglobal)(k - min_index);
//							myLambdas_sp[lambdaCount_sp][6] = (esglobal)tmp_cl;
//							lambdaCount_sp++;
//
//							//if (!flag_redund_lagr_mult)
//							//	break;
//						//}
//
//						//if (neighSD == (esglobal)MPIrank)
//						//	add = true;
//					}
//					//add = false;
//    			}
//    		}
//    	}
//
//    	if (MPIrank == 0) { std::cout << " Global B - Create global lambda numbering                                "; system("date +%T.%6N"); }
//
//    	esglobal lambdaGlobalCount_sp = 0;
//
//    	MPI_Exscan(&lambdaCount_sp, &lambdaGlobalCount_sp, 1, esglobal_mpi, MPI_SUM, MPI_COMM_WORLD);
//
//    	if (MPIrank == 0) lambdaGlobalCount_sp = 0;
//    	esglobal lambdaNum_sp = lambdaCount_sp + lambdaGlobalCount_sp;
//    	MPI_Bcast(&lambdaNum_sp, 1, esglobal_mpi, MPIsize-1, MPI_COMM_WORLD);
//    	esglobal total_number_of_global_B1_lambdas_sp = lambdaNum_sp;
//
//    	if (myLambdas_sp.size() > 0)
//    		//cilk_
//    		for (eslocal i = 0; i < myLambdas_sp.size(); i++)
//    			myLambdas_sp[i][0] = myLambdas_sp[i][0]  + lambdaGlobalCount_sp + total_number_of_global_B1_lambdas + total_number_of_B1_l_rows; // create global lambda numbering <=> increment lambda numbering by number of lambdas created by all subdomains with smaller index
//
//    	if (MPIrank == 0) { std::cout << " Global B - Assembling messages with lambdas for MPI                      "; system("date +%T.%6N"); }
//
//
//    	//std::vector < std::vector < esglobal > > mpi_send_buff;
//    	mpi_send_buff.clear();
//    	mpi_send_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
//    //#ifdef DEBUG
//    //	for (int i = 0; i < myNeighClusters.size(); i++) {
//    //#else
//    	cilk_for (int i = 0; i < myNeighClusters.size(); i++) {
//    //#endif
//    		int index = 0;
//    		if ( myLambdas_sp.size() > 0 )
//    		{
//    			for (int j = 0; j < myLambdas_sp.size(); j++) {
//    				if( myLambdas_sp[j][3] == myNeighClusters[i] ) {
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][0]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][1]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][2]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][3]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][4]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][5]);
//    					mpi_send_buff[i].push_back(myLambdas_sp[j][6]);
//    					index++;
//    				}
//    			}
//    			if (index == 0)
//    				mpi_send_buff[i].push_back(0);
//    		}
//    		else
//    		{
//    			mpi_send_buff[i].push_back(0);
//    		}
//    	}
//
//    	if (MPIrank == 0) { std::cout << " Global B - Isend                                                         "; system("date +%T.%6N"); }
//
//    	for (int i = 0; i < myNeighClusters.size(); i++)
//    		MPI_Isend(&mpi_send_buff[i][0], mpi_send_buff[i].size(), esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_send_req[i]);
//
//
//    	if (MPIrank == 0) { std::cout << " Global B - Iprobe and MPIrecv                                            "; system("date +%T.%6N"); }
//
//    	//std::vector < std::vector < esglobal > > mpi_recv_buff;
//    	mpi_recv_buff.clear();
//        mpi_recv_buff.resize( myNeighClusters.size(), std::vector< esglobal >( 0 , 0 ) );
//    	delete [] mpi_send_req;
//
//
//    	// receiving all border DOFs from all neighboring sub-domains
//    	eslocal messages_received_sp = 0;
//    	while ( messages_received_sp < myNeighClusters.size() ) {
//    		for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//
//    			//MPI_Probe( myNeighSubDomains[i], 0, MPI_COMM_WORLD, &mpi_stat);
//    			int flag;
//    			MPI_Iprobe( myNeighClusters[i], 0, MPI_COMM_WORLD, &flag, &mpi_stat );
//
//    			if (flag) {
//    				int recv_msg_size = 0;
//
//    				MPI_Get_count(&mpi_stat, esglobal_mpi, &recv_msg_size);
//    				esglobal* mpi_tmp_recv_buff = (esglobal*)malloc(sizeof(esglobal) * recv_msg_size);
//
//    				MPI_Recv(mpi_tmp_recv_buff, recv_msg_size, esglobal_mpi, myNeighClusters[i], 0, MPI_COMM_WORLD, &mpi_stat);
//    				mpi_recv_buff[i].insert(mpi_recv_buff[i].begin(), mpi_tmp_recv_buff, mpi_tmp_recv_buff + recv_msg_size);
//
//    				free(mpi_tmp_recv_buff);
//    				messages_received_sp++;
//    			}
//    		}
//    	}
//
//    	if (MPIrank == 0) { std::cout << " Global B - Decode received lambdas                                       "; system("date +%T.%6N"); }
//
//    	// decode received lambdas
//    	eslocal recv_lamba_count_sp = 0;
//    	for (eslocal i = 0; i < myNeighClusters.size(); i++)
//    		if (mpi_recv_buff[i].size() > 1)
//    			recv_lamba_count_sp += mpi_recv_buff[i].size();
//
//    	recv_lamba_count_sp = recv_lamba_count_sp / 7;
//
//    	eslocal l_i_sp = myLambdas_sp.size();
//    	myLambdas_sp.resize( myLambdas_sp.size() + recv_lamba_count_sp, std::vector< esglobal >( 7 , 0 ) );
//
//    	for (eslocal i = 0; i < myNeighClusters.size(); i++) {
//    		if (mpi_recv_buff[i].size() > 1) {
//    			for (int j = 0; j < mpi_recv_buff[i].size() / 7; j++) {
//    				myLambdas_sp[l_i_sp][0] = mpi_recv_buff[i][7*j + 0];
//    				myLambdas_sp[l_i_sp][1] = mpi_recv_buff[i][7*j + 1];
//    				myLambdas_sp[l_i_sp][2] = mpi_recv_buff[i][7*j + 3];
//    				myLambdas_sp[l_i_sp][3] = mpi_recv_buff[i][7*j + 2];
//    				myLambdas_sp[l_i_sp][4] = mpi_recv_buff[i][7*j + 4];
//    				myLambdas_sp[l_i_sp][5] = mpi_recv_buff[i][7*j + 6];
//    				myLambdas_sp[l_i_sp][6] = mpi_recv_buff[i][7*j + 5];
//    				l_i_sp++;
//    			}
//    		}
//    	}
//
//    	if (MPIrank == 0) { std::cout << " Global B - myLambdas - sort or tbb:sort                                  "; system("date +%T.%6N"); }
//
//    //	auto comp_vf = [](const std::vector<esglobal> &a, const std::vector<esglobal> &b) {return a[0] < b[0];};
//
//    //#ifdef USE_TBB
//    //	tbb::parallel_sort(myLambdas.begin(), myLambdas.end(), comp_vf);
//    //#else
//    	std::sort         (myLambdas_sp.begin(), myLambdas_sp.end(), comp_vf);
//    //#endif
//
//    	if (MPIrank == 0) { std::cout << " Global B - Final B assembling with g2l mapping using std::map            "; system("date +%T.%6N"); }
//
//    	esglobal lambda_sp;
//    	esglobal DOFNumber_sp;
//    	eslocal  n_myLambdas_sp = myLambdas_sp.size();
//
//    	std::vector < SparseDOKMatrix<T> > B1_DOK_sp(subDomPerCluster);
//
//    //#ifdef USE_TBB
//    //	tbb::mutex m;
//    //	cilk_for (int j = 0; j < n_myLambdas; j++)
//    //#else
//    	for (eslocal j = 0; j < n_myLambdas_sp; j++)
//    //#endif
//    	{
//
//    		esglobal lambda_sp         = myLambdas_sp[j][0];
//    		esglobal DOFNumber_sp      = myLambdas_sp[j][1];
//
//    		//vychazi z mapovani g2l pro uzly a ne DOFy
//    		esglobal dofNODEnumber      = DOFNumber_sp / DOFS_PER_NODE;
//    		eslocal  dofNODEoffset      = DOFNumber_sp % DOFS_PER_NODE;
//
//    		if ( myLambdas_sp[j][2] == myLambdas_sp[j][3] ) { // resim vazby mezi domenama uvnitr clusteru
//
//    			eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );
//				const std::set    < eslocal >  & subs_with_element = local_boundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel
//				std::vector < eslocal >    subs_with_elem;
//				for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set)
//					subs_with_elem.push_back( *it_set );
//
//				eslocal d1                = subs_with_elem[ myLambdas_sp[j][5] ];
//				eslocal d2                = subs_with_elem[ myLambdas_sp[j][6] ];
//
//				eslocal domDofNODENumber1 =coordinates.localIndex(clustDofNODENumber, d1);
//				eslocal domDofNODENumber2 =coordinates.localIndex(clustDofNODENumber, d2);
//
//				eslocal domDOFNumber1     = DOFS_PER_NODE * domDofNODENumber1 + dofNODEoffset;
//				eslocal domDOFNumber2     = DOFS_PER_NODE * domDofNODENumber2 + dofNODEoffset;
//
//				B1_DOK_sp[d1](lambda_sp, domDOFNumber1) =  1.0;
//				B1_DOK_sp[d2](lambda_sp, domDOFNumber2) = -1.0;
//
//				myLambdas_sp[j].push_back(d1);
//				myLambdas_sp[j].push_back(d2);
//
//				eslocal  cnt            = myLambdas_sp[j][4];
//				B1_duplicity[d1].push_back(1.0 / cnt);
//				B1_duplicity[d2].push_back(1.0 / cnt);
//
//				vec_c[d1].push_back(0.0);
//				vec_c[d2].push_back(0.0);
//
//	    		std::vector < eslocal > tmp_vec1 (2,0);		//TODO: must be esglobal
//	    		tmp_vec1[0] = myLambdas_sp[j][0];
//	    		tmp_vec1[1] = myLambdas_sp[j][2];
//
//	    		lambda_map_sub_clst.push_back(tmp_vec1);
//
//	    		std::vector < eslocal > tmp_vec2 (2,0);		//TODO: must be esglobal
//	    		tmp_vec2[0] = myLambdas_sp[j][0];
//	    		tmp_vec2[1] = myLambdas_sp[j][3];
//
//	    		lambda_map_sub_clst.push_back(tmp_vec2);
//
//	    		eslocal lam_tmp = myLambdas_sp [j][0]; 		//TODO: esglobal
//	    		lambda_map_sub_B1[ d1 ].push_back( lam_tmp );
//	        	lambda_map_sub_B1[ d2 ].push_back( lam_tmp );
//
//
//    		} else { // resim vazby mezi clustery
//
//    			eslocal  clustDofNODENumber = coordinates.clusterIndex( dofNODEnumber );
//				const std::set    < eslocal >  & subs_with_element = local_boundaries[clustDofNODENumber]; //_boundaries[clustDofNODENumber]; // mnozina podoblasti na ktery je tento uzel
//				std::vector < eslocal >    subs_with_elem;
//				for (it_set = subs_with_element.begin(); it_set != subs_with_element.end(); ++it_set)
//					subs_with_elem.push_back( *it_set );
//
//
//				eslocal  cnt            = myLambdas_sp[j][4];
//				double B_value;
//				if (myLambdas_sp[j][2] < myLambdas_sp[j][3])
//					B_value = 1.0;
//				else
//					B_value = -1.0;
//
//				eslocal d                = subs_with_elem[ myLambdas_sp[j][5] ];
//				//eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
//				eslocal domDofNODENumber = coordinates.localIndex(clustDofNODENumber, d);
//				eslocal domDOFNumber     = DOFS_PER_NODE * domDofNODENumber + dofNODEoffset;
//				// #ifdef USE_TBB
//				// m.lock();
//				// #endif
//				B1_DOK_sp[d](lambda_sp, domDOFNumber) = B_value;
//				myLambdas_sp[j].push_back(d);
//				B1_duplicity[d].push_back(1.0 / cnt);
//				vec_c[d].push_back(0.0);
//				// #ifdef USE_TBB
//				// m.unlock();
//				// #endif
//
//	    		std::vector < eslocal > tmp_vec (3,0);		//TODO: must be esglobal
//	    		tmp_vec[0] = myLambdas_sp[j][0];
//	    		tmp_vec[1] = myLambdas_sp[j][2];
//	    		tmp_vec[2] = myLambdas_sp[j][3];
//
//	    		lambda_map_sub_clst.push_back(tmp_vec);
//
//	    		eslocal lam_tmp = myLambdas_sp [j][0]; 		//TODO: esglobal
//	    		lambda_map_sub_B1[ d ].push_back( lam_tmp );
//
//    		}
//    	}
//
//
//    	for (eslocal d = 0; d < subDomPerCluster; d++){
//    		//TODO: lambdaNum muze byt 64bit integer
//    		//TODO: matice B - pocet radku muze byt 64bit int
//    		B1_DOK_sp[d].resize( total_number_of_B1_l_rows + total_number_of_global_B1_lambdas + total_number_of_global_B1_lambdas_sp , K_mat[d].rows());
//    		SparseIJVMatrix<T> ijv = B1_DOK_sp[d];
//    		B1[d].AppendMatrix(ijv); //    = B1_DOK_tmp[d];
//    	}
//
//    	if (MPIrank == 0) { std::cout << " Global B - END                                                           "; system("date +%T.%6N"); }
//
//    	if (MPIrank == 0) { std::cout << " Creating lambda_map_sub vector of vectors - Global B1                    "; system("date +%T.%6N"); }
//
//
//
//
//    	if (MPIrank == 0) { std::cout << " END - Creating lambda_map_sub vector of vectors - Global B1              "; system("date +%T.%6N"); }
//
//    	MPI_Barrier(MPI_COMM_WORLD);
//
//        if (MPIrank == 0) { std::cout << " Dual size: " <<  total_number_of_B1_l_rows + total_number_of_global_B1_lambdas + total_number_of_global_B1_lambdas_sp  << std::endl; }
}

template <>
void Linear<FEM>::RHS()
{
	//	const std::map<eslocal, double> &forces_x = this->_mesh.coordinates().property(mesh::FORCES_X).values();
	//	const std::map<eslocal, double> &forces_y = this->_mesh.coordinates().property(mesh::FORCES_Y).values();
	//	const std::map<eslocal, double> &forces_z = this->_mesh.coordinates().property(mesh::FORCES_Z).values();
	//
	//	for (size_t p = 0; p < this->_mesh.parts(); p++) {
	//		for (eslocal i = 0; i < this->_mesh.coordinates().localSize(p); i++) {
	//			if (forces_x.find(l2g_vec[d][iz]) != forces_x.end()) {
	//				f_vec[d][3 * iz + 0] = forces_x.at(l2g_vec[d][iz]);
	//			}
	//			if (forces_y.find(l2g_vec[d][iz]) != forces_y.end()) {
	//				f_vec[d][3 * iz + 1] = forces_y.at(l2g_vec[d][iz]);
	//			}
	//			if (forces_z.find(l2g_vec[d][iz]) != forces_z.end()) {
	//				f_vec[d][3 * iz + 2] = forces_z.at(l2g_vec[d][iz]);
	//			}
	//		}
	//	}
}

}
