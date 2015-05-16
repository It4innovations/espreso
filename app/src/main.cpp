#include "mpi.h"

#include "esmesh.h"
#include "essolver.h"
#include "espmcube.h"

#include <vector>
#include <iostream>

struct FEMInput {

	FEMInput(): mesh(coordinates) { };

	Coordinates coordinates;
	Mesh mesh;
	std::map<int, double> dirichlet_x;
	std::map<int, double> dirichlet_y;
	std::map<int, double> dirichlet_z;
};

struct FEMParams {

	FEMParams(): generateMesh(false), subdomains(3), elementsInSub(3) { };

	bool generateMesh;
	std::vector<int> subdomains;
	std::vector<int> elementsInSub;
};

FEMInput input;
FEMParams params;


void setParams(int argc, char** argv)
{
	if (argc != 7) {
		return;
	}

	params.generateMesh = true;
	int subdomains;
	int elementsInSub;

	for (int i = 0; i < 3; i++) {
		sscanf(argv[2 * i + 1], "%i", &subdomains);
		sscanf(argv[2 * i + 2], "%i", &elementsInSub);
		std::cout << subdomains << " " << elementsInSub << " ";
		params.subdomains[i] = subdomains;
		params.elementsInSub[i] = elementsInSub;
	}
	std::cout << "\n";
}

void testFEM(int argc, char** argv);

void testBEM(int argc, char** argv);

void load_mesh();

void generate_mesh();

int main(int argc, char** argv)
{
	setParams(argc, argv);
	if (params.generateMesh) {
		generate_mesh();
	} else {
		load_mesh();
	}

	testFEM(argc, argv);
}


void load_mesh()
{
	input.coordinates = Coordinates("matrices/TET/10/coord");
	input.mesh = Mesh("matrices/TET/10/elem", input.coordinates, 4, 8);

	// fix down face
	for (int i = 0; i < 11 * 11; i++) {
		input.dirichlet_x[i + input.coordinates.getOffset()] = 0;
		input.dirichlet_y[i + input.coordinates.getOffset()] = 0;
		input.dirichlet_z[i + input.coordinates.getOffset()] = 0;
	}
}

void generate_mesh()
{
	std::cout << "mesh_generator3d" << std::endl;
	CFem::mesh_generator3d(input.mesh, input.coordinates, &(params.subdomains[0]), &(params.elementsInSub[0]));
	std::cout << "dirichlet" << std::endl;
	CFem::dirichlet(input.dirichlet_x, input.dirichlet_y, input.dirichlet_z, &(params.subdomains[0]), &(params.elementsInSub[0]));
	std::cout << "fix points" << std::endl;

	// TODO: set fix points in PERMONCUBE
	input.mesh.computeFixPoints(4);
	std::cout << "fix points - end" << std::endl;
}


template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v)
{
	//out << "[";
	size_t last = v.size() - 1;
	for(size_t i = 0; i < v.size(); ++i) {
		out << v[i];
		if (i != last)
			out << ", ";
	}
	//out << "]";
	return out;
}

void testFEM(int argc, char** argv)
{
	double start;
	start = omp_get_wtime();
	std::cout.precision(15);

	size_t partsCount = input.mesh.getPartsCount();
	size_t fixPointsCount = input.mesh.getFixPointsCount();

	std::cout << "4 : " << omp_get_wtime() - start<< std::endl;

	// TODO: fill boundaries in PERMONCUBE
	Boundaries boundaries(input.mesh, input.coordinates);

	std::cout << "5 : " << omp_get_wtime() - start<< std::endl;

	//Faces faces(mesh, coordinates);

	std::cout << "6 : " << omp_get_wtime() - start<< std::endl;

	//Corners corners(faces.getFaces(), coordinates);

	std::cout << "7 : " << omp_get_wtime() - start<< std::endl;

	std::vector < SparseCSRMatrix >			K_mat;
	std::vector < SparseCSRMatrix >			M_mat;
	std::vector < SparseIJVMatrix >			B1_mat;
	std::vector < SparseIJVMatrix >			B0_mat;

	std::vector < std::vector <int> >		lambda_map_sub_B1;
	std::vector < std::vector <int> >		lambda_map_sub_B0;
	std::vector < std::vector <int> >		lambda_map_sub_clst;
	std::vector < std::vector <double> >	B1_l_duplicity;

	std::vector < std::vector < double > >	f_vec     (partsCount);
	std::vector < std::vector < int > >		fix_nodes (partsCount);
	std::vector < std::vector <int> >		l2g_vec;

	std::cout << "8 : " << omp_get_wtime() - start<< std::endl;

	K_mat.reserve(partsCount);
	M_mat.reserve(partsCount);
	for (int d = 0; d < partsCount; d++) {
		K_mat.push_back( SparseCSRMatrix (0,0) );
		M_mat.push_back( SparseCSRMatrix (0,0) );
	}

	std::cout << "9 : " << omp_get_wtime() - start<< std::endl;

	cilk_for (int d = 0; d < partsCount; d++) {

		int dimension = input.mesh.getPartNodesCount(d) * Point::size();
		std::vector<double> f(dimension);

		input.mesh.elasticity(K_mat[d], M_mat[d], f, d);

		//K_mat[d] = K;
		//M_mat[d] = M;
		f_vec[d].swap(f);

		std::cout << d << " " << std::endl;
	}

	std::cout << "10: " << omp_get_wtime() - start<< std::endl;

	const std::vector<idx_t> fixPoints = input.mesh.getFixPoints();

	cilk_for (int d = 0; d < partsCount; d++) {
		for (int fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
			fix_nodes[d].push_back(fixPoints[d * fixPointsCount + fixPoint]);
		}
		std::sort ( fix_nodes[d].begin(), fix_nodes[d].end() );
	}

	std::cout << "11: " << omp_get_wtime() - start<< std::endl;
	boundaries.create_B1_l(
		B1_mat,
		B0_mat,
		l2g_vec,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		lambda_map_sub_B0,
		B1_l_duplicity,
		input.dirichlet_x,
		input.dirichlet_y,
		input.dirichlet_z,
		partsCount,
		input.mesh
	);

	std::cout << "12: " << omp_get_wtime() - start<< std::endl;

	std::cout.precision(10);

	// Start - Stupid version of ESPRESO interface

	MPI_Init (&argc, &argv);					// starts MPI

	typedef int       ShortInt ;
	typedef int       longInt  ;


	int MPIrank = 0; //MPI_Comm_rank(fem->comm, &MPIrank);
	int MPIsize = 1; //MPI_Comm_size(fem->comm, &MPIsize);
	int number_of_subdomains_per_cluster = partsCount;


	extern void SetCluster		  ( Cluster & cluster, ShortInt * subdomains_global_indices, ShortInt number_of_subdomains, ShortInt MPI_rank);

	extern void SetMatrixB1_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type );

	extern void SetMatrixB0_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type );

	extern void SetMatrixR_fromDense( Cluster & cluster, ShortInt domain_index_in_cluster,
		ShortInt n_cols, ShortInt n_rows, double * vals, char type );

	extern void SetMatrixK_fromCSR ( Cluster & cluster, ShortInt domain_index_in_cluster,
		ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type );

	extern void SetSolverPreprocessing ( Cluster & cluster, IterSolver & solver,
		vector <vector <longInt> > & lambda_map_sub, vector < ShortInt > & neigh_domains );

	extern void SetMatrixFromCSR   ( SparseMatrix    & Mat, ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type );
	extern void SetMatrixFromDense ( SparseMatrix    & Mat, ShortInt n_cols, ShortInt n_rows, double * vals, char type );
	extern void SetMatrixFromCOO   ( SparseMatrix    & Mat, ShortInt n_rows, ShortInt n_cols, ShortInt nnz, ShortInt * I_rows, ShortInt * J_cols, double * V_vals, char type );
	extern void SetVecInt          ( vector <int>    & vec, ShortInt incerement_by, ShortInt nnz, ShortInt * vals);
	extern void SetVecDbl          ( vector <double> & vec, ShortInt nnz,	double * vals);



	Cluster cluster(MPIrank + 1);
	cluster.USE_DYNAMIC			= 0;
	cluster.USE_HFETI			= 0;
	cluster.USE_KINV			= 0;
	cluster.SUBDOM_PER_CLUSTER	= number_of_subdomains_per_cluster;
	cluster.NUMBER_OF_CLUSTERS	= MPIsize;

	IterSolver solver;
	solver.CG_max_iter	 = 100;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = 0.0001;
	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = 0;
	solver.USE_PIPECG	 = 1;
	solver.USE_PREC		 = 0;
	solver.FIND_SOLUTION = 0;


	std::vector <int> domain_list (number_of_subdomains_per_cluster,0);
	for (int i = 0; i<number_of_subdomains_per_cluster; i++)
		domain_list[i] = i;

	SetCluster( cluster, &domain_list[0], number_of_subdomains_per_cluster, MPIrank);

	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );

	// *** Setup B0 matrix *******************************************************************************************
	if (cluster.USE_HFETI == 1 ) {

		cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
			ShortInt domain_index_in_cluster = i;

			SetMatrixB0_fromCOO( cluster, domain_index_in_cluster,
				B0_mat[i].rows(),			//clust_g.data[i]->B->B0_rows,		// B_full_rows, //n_row_eq,
				B0_mat[i].columns(),		//.data[i]->B->B0_cols,				// B_full_cols, //n_col,
				B0_mat[i].nonZeroValues(),	//.data[i]->B->B0_nnz,				// B_full_nnz,  //nnz_eq,
				B0_mat[i].rowIndices(),		//&clust_g.data[i]->B->B0_I[0],		// BI_full[0], //Bi_coo,
				B0_mat[i].columnIndices(),	//&clust_g.data[i]->B->B0_J[0],		// BJ_full[0], //Bj_coo,
				B0_mat[i].values(),			//&clust_g.data[i]->B->B0_V[0],		// BV_full[0], //Bv_coo,
				'G' );
		}
	}
	// *** END - Setup B0 matrix *************************************************************************************

	// *** Setup B1 matrix *******************************************************************************************
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
		ShortInt domain_index_in_cluster = i;
		SetMatrixB1_fromCOO( cluster, domain_index_in_cluster,
			B1_mat[i].rows(),			//clust_g.data[i]->B->B_full_rows, //n_row_eq,
			B1_mat[i].columns(),		//clust_g.data[i]->B->B_full_cols, //n_col,
			B1_mat[i].nonZeroValues(),	//clust_g.data[i]->B->B_full_nnz,  //nnz_eq,
			B1_mat[i].rowIndices(),		//&clust_g.data[i]->B->BI_full[0], //Bi_coo,
			B1_mat[i].columnIndices(),	//&clust_g.data[i]->B->BJ_full[0], //Bj_coo,
			B1_mat[i].values(),			//&clust_g.data[i]->B->BV_full[0], //Bv_coo,
			'G' );
	}

	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
		cluster.domains[i].B1_scale_vec = B1_l_duplicity[i];
	}
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	cilk_for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
		for (int i = 0; i < l2g_vec[d].size(); i++) {
			std::vector <double> tmp_vec (3,0);
			tmp_vec[0] = input.coordinates[l2g_vec[d][i]].x;
			tmp_vec[1] = input.coordinates[l2g_vec[d][i]].y;
			tmp_vec[2] = input.coordinates[l2g_vec[d][i]].z;
			cluster.domains[d].coordinates.push_back(tmp_vec);
		}
		cluster.domains[d].CreateKplus_R();
		//cluster.domains[d].Kplus_R.ConvertCSRToDense(0);
	}
	// *** END - Setup R matrix **************************************************************************************

	// *** Load RHS and fix points for K regularization **************************************************************
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
		//SetVecDbl( cluster.domains[i].f,        clust_g.data[i]->KSparse->n_row, clust_g.data[i]->fE );
		cluster.domains[d].f = f_vec[d];

		//SetVecInt( cluster.domains[i].fix_dofs, 1,                           24, clust_g.fem[i]->mesh.fixingDOFs );
		for (int i = 0; i < fix_nodes[d].size(); i++) {
 			for (int d_i = 0; d_i < 3; d_i++) {
				cluster.domains[d].fix_dofs.push_back( 3 * fix_nodes[d][i] + d_i);
			}
		}
	}
	// *** END - Load RHS and fix points for K regularization ********************************************************

	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
		cluster.domains[i].lambda_map_sub = lambda_map_sub_B1[i];
	}

	std::vector < int > neigh_clusters;
	//neigh_clusters.push_back(0);

	SetSolverPreprocessing ( cluster, solver, lambda_map_sub_clst, neigh_clusters );
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// *** Load Matrix K and regularization ******************************************************************************
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
		SetMatrixK_fromCSR ( cluster, d,
			K_mat[d].rows(), K_mat[d].columns(), //  .data[i]->KSparse->n_row,   clust_g.data[i]->KSparse->n_row,
			K_mat[d].rowPtrs(), K_mat[d].columnIndices(), K_mat[d].values(), //clust_g.data[i]->KSparse->row_ptr, clust_g.data[i]->KSparse->col_ind, clust_g.data[i]->KSparse->val,
			'G');
	}

	if (cluster.USE_HFETI == 1)
		cluster.SetClusterHFETI();

	cluster.SetClusterPC_AfterKplus();
	// *** END - Load Matrix K and regularization  ***********************************************************************


	// *** Running Solver ************************************************************************************************
	string result_file("MATSOL_SVN_Displacement.Nvec");
	solver.Solve_singular ( cluster, result_file );

	vector < vector < double > > prim_solution;
	solver.GetSolution_Primal_singular_parallel(cluster, prim_solution);
	double max_v = 0.0;

	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++)
		for (ShortInt j = 0; j < prim_solution[i].size(); j++)
			if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.AddStartWOBarrier(0.0); max_sol_ev.AddEndWOBarrier(max_v);

	std::cout.precision(15);

	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	if (MPIrank == 0)
		std::cout << " Max value in_solution = " << max_vg << std::endl;

	max_sol_ev.PrintLastStatMPI_PerNode(max_vg);

	input.mesh.saveVTK(prim_solution, l2g_vec);




	//if (clust_g.domainG->flag_store_VTK)
	//{
	//	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
	//		for (ShortInt j = 0; j < prim_solution[i].size(); j++) {
	//			if (prim_solution[i][j] > max_v) max_v = prim_solution[i][j];
	//		}
	//		copy(prim_solution[i].begin(), prim_solution[i].end(), clust_g.data[i]->ddu);
	//	}

	//}


	// *** END - Running Solver ************************************************************************************************






	// END - Stupid version of ESPRESO interface




	//// Export matrices to FILE

	//for (int d = 0; d < partsCount; d++) {

	//	std::stringstream ss;
	//	ss << "K_" << d;
	//
	//	std::string name = ss.str();
	//
	//	std::ofstream myfile;
	//	myfile.open ( name );
	//	myfile << K_mat[d];
	//	myfile.close();

	//}

	//for (int d = 0; d < partsCount; d++) {

	//	std::stringstream ss;
	//	ss << "B_" << d;

	//	std::string name = ss.str();

	//	std::ofstream myfile;
	//	myfile.open ( name );
	//	myfile << B1_mat[d];
	//	myfile.close();

	//}

	//for (int d = 0; d < partsCount; d++) {

	//	std::stringstream ss;
	//	ss << "f_" << d;

	//	std::string name = ss.str();

	//	std::ofstream myfile;
	//	myfile.open ( name );
	//	myfile << f_vec[d];
	//	myfile.close();

	//}

	//for (int d = 0; d < partsCount; d++) {

	//	std::stringstream ss;
	//	ss << "l2g_" << d;

	//	std::string name = ss.str();

	//	std::ofstream myfile;
	//	myfile.open ( name );
	//	myfile << l2g_vec[d];
	//	myfile.close();

	//}

	//SparseDOKMatrix dok(4, 4);
	//dok(1, 1) = 5;
	//dok(1, 2) += 5;
	//std::std::cout << dok;

	//std::std::cout << "Fix points:\n";
	//for (int part = 0; part < partsCount; part++) {
	//	std::std::cout << "Part: " << part << ": ";
	//	for (int fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
	//		std::std::cout << fixPoints[part * fixPointsCount + fixPoint] << " ";
	//	}
	//	std::std::cout << "\n";
	//}

	//std::vector <int> fix_points;

	// pocet podoblasti * pocet fix bodu / 8 * 4
	//mesh.getFixPoints();

	//std::ofstream myfile;
	//myfile.open ("K.txt");
	//myfile << K;
	//myfile.close();

	//std::cout << K;
}




