#include "mpi.h"

#include "esmesh.h"
#include "essolver.h"
#include "espmcube.h"

#include <vector>
#include <iostream>

enum {
	HEXA8,		// 0 OK
	HEXA20,		// 1 OK
	TETRA4,		// 2 OK
	TETRA10,	// 3 OK
	PRISMA6,	// 4 OK
	PRISMA15,	// 5 OK
	PYRAMID5,	// 6 OK
	PYRAMID13	// 7 incorrect indexing and global B
};

struct FEMInput {

	FEMInput(): mesh(NULL), localBoundaries(NULL), globalBoundaries(NULL) {};
	~FEMInput() {
		if (mesh != NULL) { delete mesh; }
		if (localBoundaries != NULL) { delete localBoundaries; }
		if (globalBoundaries != NULL) { delete globalBoundaries; }
	}

	mesh::Mesh *mesh;
	mesh::Boundaries *localBoundaries;
	mesh::Boundaries *globalBoundaries;
};

struct FEMParams {

	FEMParams(): type(HEXA8), generateMesh(false) {
	};

	int type;
	bool generateMesh;
	permoncube::Settings settings;
};

permoncube::Generator *generator;
FEMInput input;
FEMParams params;

void setParams(int argc, char** argv)
{
	if (argc != 11) {
		return;
	}

	params.generateMesh = true;

	int type;
	sscanf(argv[1], "%i", &type);
	params.type = type;

	int clusters;
	int subdomains;
	int elementsInSub;

	for (int i = 0; i < 3; i++) {
		sscanf(argv[i + 2], "%i", &clusters);
		sscanf(argv[i + 5], "%i", &subdomains);
		sscanf(argv[i + 8], "%i", &elementsInSub);
		params.settings.clusters[i] = clusters;
		params.settings.subdomainsInCluster[i] = subdomains;
		params.settings.elementsInSubdomain[i] = elementsInSub;
	}
}

void testFEM(int argc, char** argv);

void testBEM(int argc, char** argv);

void testMPI(int argc, char** argv, int MPIrank, int MPIsize);

void load_mesh();

void generate_mesh( int MPIrank );

int main(int argc, char** argv)
{
	setParams(argc, argv);

	MPI_Init(&argc, &argv);					// starts MPI

	int MPIrank;
	int MPIsize;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

	if (params.settings.clusters[0] * params.settings.clusters[1] * params.settings.clusters[2] != MPIsize) {
		std::cerr << "Invalid number of processes.\n";
		exit(EXIT_FAILURE);
	}

	if (params.generateMesh) {
		generate_mesh( MPIrank );
	} else {
		load_mesh();
	}

	if (params.settings.clusters[0] * params.settings.clusters[1] * params.settings.clusters[2] == 1) {
		//testBEM(argc, argv);
		//testFEM(argc, argv);
		testMPI(argc, argv, MPIrank, MPIsize);
	} else {
		testMPI(argc, argv, MPIrank, MPIsize);
	}

	if (params.generateMesh) {
		delete generator;
	}

	MPI_Finalize();
}


void load_mesh()
{
	mesh::Ansys ansys("matrices/spanner/Model");
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_X) = "BC/Elasticity/NUX.dat";
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_Y) = "BC/Elasticity/NUY.dat";
	ansys.coordinatesProperty(mesh::CP::DIRICHLET_Z) = "BC/Elasticity/NUZ.dat";
	ansys.coordinatesProperty(mesh::CP::FORCES_X) = "BC/Elasticity/NFX.dat";
	ansys.coordinatesProperty(mesh::CP::FORCES_Y) = "BC/Elasticity/NFY.dat";
	ansys.coordinatesProperty(mesh::CP::FORCES_Z) = "BC/Elasticity/NFZ.dat";

	input.mesh = new mesh::Mesh(ansys, 2, 8);
	input.localBoundaries = new mesh::Boundaries(*input.mesh);
}

void generate_mesh( int MPIrank )
{
	if (MPIrank == 0) { std::cout << "Permoncube - start                                                       "; system("date +%T.%6N"); }

	switch (params.type) {
	case HEXA8: {
		generator = new permoncube::ElementGenerator<permoncube::Hexahedron8>(params.settings);
		break;
	}
	case TETRA10: {
		generator = new permoncube::ElementGenerator<permoncube::Tetrahedron10>(params.settings);
		break;
	}
	case TETRA4: {
		generator = new permoncube::ElementGenerator<permoncube::Tetrahedron4>(params.settings);
		break;
	}
	case HEXA20: {
		generator = new permoncube::ElementGenerator<permoncube::Hexahedron20>(params.settings);
		break;
	}
	case PRISMA6: {
		generator = new permoncube::ElementGenerator<permoncube::Prisma6>(params.settings);
		break;
	}
	case PRISMA15: {
		generator = new permoncube::ElementGenerator<permoncube::Prisma15>(params.settings);
		break;
	}
	case PYRAMID5: {
		generator = new permoncube::ElementGenerator<permoncube::Pyramid5>(params.settings);
		break;
	}
	case PYRAMID13: {
		generator = new permoncube::ElementGenerator<permoncube::Pyramid13>(params.settings);
		break;
	}
	}

	size_t index = MPIrank;

	int globoalIndClst = MPIrank;
	int Cx = params.settings.clusters[0];
	int Cy = params.settings.clusters[1];
	int Cz = params.settings.clusters[2];

	int Iclst,Jclst,Kclst, inXYplane;

	Kclst = ceil(double (globoalIndClst+1)/(Cx*Cy));
	inXYplane = (globoalIndClst+1)-(Kclst-1)*Cx*Cy;
	Jclst = ceil(double( inXYplane)/Cx);
	Iclst = inXYplane - Cx*(Jclst-1);

	size_t cluster[3];
	cluster[0] = Iclst - 1;
	cluster[1] = Jclst - 1;
	cluster[2] = Kclst - 1;

	input.mesh = new mesh::Mesh();
	input.globalBoundaries = new mesh::Boundaries(*input.mesh);
	generator->mesh(*input.mesh, cluster);
	input.localBoundaries = new mesh::Boundaries(*input.mesh);
	//generator->fixZeroPlanes(input.mesh, cluster);
	generator->fixBottom(*input.mesh, cluster);

	generator->fillGlobalBoundaries(*input.globalBoundaries, cluster);

	generator->setFixPoints(*input.mesh);
	//input.mesh->computeFixPoints(8);

	size_t number[3] = { 3, 3, 3 };
	generator->setCorners(*input.localBoundaries, number, true, true, true);

	if (MPIrank == 0) { std::cout << "Permoncube - end                                                         "; system("date +%T.%6N"); }
}


void testMPI(int argc, char** argv, int MPIrank, int MPIsize)
{


	TimeEval timeEvalMain(string("ESPRESO Solver Overal Timing"));
	timeEvalMain.totalTime.AddStartWithBarrier();



	eslocal clust_index = MPIrank; // clust_index = cislo clusteru

	double start;
	start = omp_get_wtime();
	std::cout.precision(15);

	size_t partsCount 	  = input.mesh->parts();
	size_t fixPointsCount = input.mesh->getFixPointsCount();

	//if (MPIrank == 0) std::cout << "4 : " << omp_get_wtime() - start<< std::endl;

	 TimeEvent timeBoundaries(string("Create Boundaries from Mesh"));
	 timeBoundaries.AddStart();

	// TODO: fill boundaries in PERMONCUBE
	mesh::Boundaries boundaries(*input.mesh);

	 timeBoundaries.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeBoundaries);

	//if (MPIrank == 0) std::cout << "5 : " << omp_get_wtime() - start<< std::endl;

	 TimeEvent timeFaces(string("Create Faces"));
	 timeFaces.AddStart();
	//Faces faces(mesh, coordinates);
	 timeFaces.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeFaces);

	//if (MPIrank == 0) std::cout << "6 : " << omp_get_wtime() - start<< std::endl;

	 TimeEvent timeCorners(string("Create Corners"));
	 timeCorners.AddStart();
	//Corners corners(faces.getFaces(), coordinates);
	 timeCorners.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeCorners);

	//if (MPIrank == 0) std::cout << "7 : " << omp_get_wtime() - start<< std::endl;

	 TimeEvent timeKasm(string("Create K matrices"));
	 timeKasm.AddStart();

	std::vector < SparseCSRMatrix<eslocal> >			K_mat;
	std::vector < SparseCSRMatrix<eslocal> >			M_mat;
	std::vector < SparseIJVMatrix<eslocal> >			B1_mat;
	std::vector < SparseIJVMatrix<eslocal> >			B0_mat;

	std::vector < std::vector <eslocal> >	lambda_map_sub_B1;
	std::vector < std::vector <eslocal> >	lambda_map_sub_B0;
	std::vector < std::vector <eslocal> >	lambda_map_sub_clst;
	std::vector < std::vector <double> >	B1_duplicity;

	std::vector < std::vector <double > >	f_vec     (partsCount);
	std::vector < std::vector <eslocal > >	fix_nodes (partsCount);
	std::vector < std::vector <eslocal> >	l2g_vec;

	//if (MPIrank == 0) std::cout << "8 : " << omp_get_wtime() - start<< std::endl;

	K_mat.reserve(partsCount);
	M_mat.reserve(partsCount);
	for (eslocal d = 0; d < partsCount; d++) {
		K_mat.push_back( SparseCSRMatrix<eslocal>(0,0) );
		M_mat.push_back( SparseCSRMatrix<eslocal>(0,0) );
	}

	//if (MPIrank == 0) std::cout << "9 : " << omp_get_wtime() - start<< std::endl;

#ifndef DEBUG
	cilk_for (eslocal d = 0; d < partsCount; d++) {
#else
	for (eslocal d = 0; d < partsCount; d++) {
#endif
		eslocal dimension = input.mesh->getPartNodesCount(d) * mesh::Point::size();
		std::vector<double> f(dimension);

		//input.mesh.elasticity(K_mat[d], M_mat[d], f, d);
		input.mesh->elasticity(K_mat[d], f, d);

        f_vec[d].swap(f);

        if (MPIrank == 0) std::cout << d << " " ; //<< std::endl;
	}

	//if (MPIrank == 0) std::cout << std::endl << "10: " << omp_get_wtime() - start<< std::endl;

	 timeKasm.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeKasm);



	TimeEvent timeFnodes(string("Create Fix nodes"));
	timeFnodes.AddStart();
	const std::vector<eslocal> fixPoints = input.mesh->getFixPoints();

#ifndef DEBUG
	cilk_for (eslocal d = 0; d < partsCount; d++) {
#else
	for (eslocal d = 0; d < partsCount; d++) {
#endif
		for (eslocal fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
			fix_nodes[d].push_back(fixPoints[d * fixPointsCount + fixPoint]);
		}
		std::sort ( fix_nodes[d].begin(), fix_nodes[d].end() );
	}

	 timeFnodes.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeFnodes);

	 TimeEvent timeB1loc(string("Create B1 local"));
	 timeB1loc.AddStart();

	//if (MPIrank == 0) std::cout << "11: " << omp_get_wtime() - start<< std::endl;
	input.localBoundaries->create_B1_l<eslocal>(
		B1_mat,
		B0_mat,
		l2g_vec,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		lambda_map_sub_B0,
		B1_duplicity,
		partsCount
	);

	 timeB1loc.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1loc);

	std::vector < eslocal > neigh_clusters;

	//if (MPIrank == 0) std::cout << "11.1: " << omp_get_wtime() - start<< std::endl;

	 TimeEvent timeB1glob(string("Create B1 global"));
	 timeB1glob.AddStart();

	input.globalBoundaries->create_B1_g<eslocal>(
		B1_mat,
		K_mat,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		B1_duplicity,
		MPIrank,
		MPIsize,
		partsCount,
		neigh_clusters,
		boundaries
	);

	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);

	 TimeEvent timeBforces(string("Create boundary forces ??"));
	 timeBforces.AddStart();

	const std::map<eslocal, double> &forces_x = input.mesh->coordinates().property(mesh::CP::FORCES_X).values();
	const std::map<eslocal, double> &forces_y = input.mesh->coordinates().property(mesh::CP::FORCES_Y).values();
	const std::map<eslocal, double> &forces_z = input.mesh->coordinates().property(mesh::CP::FORCES_Z).values();

	for (eslocal d = 0; d < partsCount; d++) {
		for (eslocal iz = 0; iz < l2g_vec[d].size(); iz++) {
			if (forces_x.find(l2g_vec[d][iz]) != forces_x.end()) {
				f_vec[d][3 * iz + 0] = forces_x.at(l2g_vec[d][iz]);
			}
			if (forces_y.find(l2g_vec[d][iz]) != forces_y.end()) {
				f_vec[d][3 * iz + 1] = forces_y.at(l2g_vec[d][iz]);
			}
			if (forces_z.find(l2g_vec[d][iz]) != forces_z.end()) {
				f_vec[d][3 * iz + 2] = forces_z.at(l2g_vec[d][iz]);
			}
		}
	}

	 timeBforces.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeBforces);

	//if (MPIrank == 0) std::cout << "12: " << omp_get_wtime() - start<< std::endl;

	std::cout.precision(10);


	// Start - Stupid version of ESPRESO interface

	typedef int       ShortInt ;
	typedef int       longInt  ;

	int number_of_subdomains_per_cluster = partsCount;

	extern void SetCluster		  ( Cluster & cluster, ShortInt * subdomains_global_indices, ShortInt number_of_subdomains, ShortInt MPI_rank);

	extern void SetMatrixB1_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );

	extern void SetMatrixB0_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );

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
	cluster.USE_HFETI			= 1;
	cluster.USE_KINV			= 0;
	cluster.SUBDOM_PER_CLUSTER	= number_of_subdomains_per_cluster;
	cluster.NUMBER_OF_CLUSTERS	= MPIsize;

	IterSolver solver;
	solver.CG_max_iter	 = 1000;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = 0.0001;
	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = 0;
	solver.USE_PIPECG	 = 0;
	solver.USE_PREC		 = 1;
	solver.FIND_SOLUTION = 0;

	 TimeEvent timeSetClust(string("Solver - Set cluster"));
	 timeSetClust.AddStart();

	std::vector <int> domain_list (number_of_subdomains_per_cluster,0);
	for (int i = 0; i<number_of_subdomains_per_cluster; i++)
		domain_list[i] = i;

	SetCluster( cluster, &domain_list[0], number_of_subdomains_per_cluster, MPIrank);

	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );

	 timeSetClust.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetClust);

	// *** Setup B0 matrix *******************************************************************************************
	 if (cluster.USE_HFETI == 1 ) {

		 TimeEvent timeSetB0(string("Solver - Set B0"));
		 timeSetB0.AddStart();

#ifndef DEBUG
		cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
		for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
			ShortInt domain_index_in_cluster = i;

			SetMatrixB0_fromCOO( cluster, domain_index_in_cluster,
				B0_mat[i].rows(),			//clust_g.data[i]->B->B0_rows,		// B_full_rows, //n_row_eq,
				B0_mat[i].columns(),		//.data[i]->B->B0_cols,				// B_full_cols, //n_col,
				B0_mat[i].nonZeroValues(),	//.data[i]->B->B0_nnz,				// B_full_nnz,  //nnz_eq,
				B0_mat[i].rowIndices(),		//&clust_g.data[i]->B->B0_I[0],		// BI_full[0], //Bi_coo,
				B0_mat[i].columnIndices(),	//&clust_g.data[i]->B->B0_J[0],		// BJ_full[0], //Bj_coo,
				B0_mat[i].values(),			//&clust_g.data[i]->B->B0_V[0],		// BV_full[0], //Bv_coo,
				'G',
				B0_mat[i].indexing() );
		}

		 timeSetB0.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeSetB0);
	}
	// *** END - Setup B0 matrix *************************************************************************************

	// *** Setup B1 matrix *******************************************************************************************
	 TimeEvent timeSetB1(string("Solver - Set B1"));
	 timeSetB1.AddStart();

#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		ShortInt domain_index_in_cluster = i;
		SetMatrixB1_fromCOO( cluster, domain_index_in_cluster,
			B1_mat[i].rows(),			//clust_g.data[i]->B->B_full_rows, //n_row_eq,
			B1_mat[i].columns(),		//clust_g.data[i]->B->B_full_cols, //n_col,
			B1_mat[i].nonZeroValues(),	//clust_g.data[i]->B->B_full_nnz,  //nnz_eq,
			B1_mat[i].rowIndices(),		//&clust_g.data[i]->B->BI_full[0], //Bi_coo,
			B1_mat[i].columnIndices(),	//&clust_g.data[i]->B->BJ_full[0], //Bj_coo,
			B1_mat[i].values(),			//&clust_g.data[i]->B->BV_full[0], //Bv_coo,
			'G',
			B1_mat[i].indexing() );
	}

#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		cluster.domains[i].B1_scale_vec = B1_duplicity[i];
	}

	 timeSetB1.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetB1);
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	 TimeEvent timeSetR(string("Solver - Set R"));
	 timeSetR.AddStart();
#ifndef DEBUG
	cilk_for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
		for (int i = 0; i < l2g_vec[d].size(); i++) {
			std::vector <double> tmp_vec (3,0);
			tmp_vec[0] = input.mesh->coordinates()[l2g_vec[d][i]].x;
			tmp_vec[1] = input.mesh->coordinates()[l2g_vec[d][i]].y;
			tmp_vec[2] = input.mesh->coordinates()[l2g_vec[d][i]].z;
			cluster.domains[d].coordinates.push_back(tmp_vec);
		}
		cluster.domains[d].CreateKplus_R();
		//cluster.domains[d].Kplus_R.ConvertCSRToDense(0);
	}
	 timeSetR.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetR);
	// *** END - Setup R matrix **************************************************************************************

	// *** Load RHS and fix points for K regularization **************************************************************
	 TimeEvent timeSetRHS(string("Solver - Set RHS and Fix points"));
	 timeSetRHS.AddStart();
#ifndef DEBUG
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
		//SetVecDbl( cluster.domains[i].f,        clust_g.data[i]->KSparse->n_row, clust_g.data[i]->fE );
		cluster.domains[d].f = f_vec[d];

		//SetVecInt( cluster.domains[i].fix_dofs, 1,                           24, clust_g.fem[i]->mesh.fixingDOFs );
		for (int i = 0; i < fix_nodes[d].size(); i++) {
 			for (int d_i = 0; d_i < 3; d_i++) {
				cluster.domains[d].fix_dofs.push_back( 3 * fix_nodes[d][i] + d_i);
			}
		}
	}
	 timeSetRHS.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetRHS);
	// *** END - Load RHS and fix points for K regularization ********************************************************

	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	 TimeEvent timeSolPrec(string("Solver - FETI Preprocessing"));
	 timeSolPrec.AddStart();
#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		cluster.domains[i].lambda_map_sub = lambda_map_sub_B1[i];
	}

	SetSolverPreprocessing ( cluster, solver, lambda_map_sub_clst, neigh_clusters );

	 timeSolPrec.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolPrec);
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// *** Load Matrix K and regularization ******************************************************************************
	 TimeEvent timeSolKproc(string("Solver - K regularization and factorization"));
	 timeSolKproc.AddStart();

	if (MPIrank == 0) std::cout << "K regularization and factorization ... " << std::endl ;
	#ifndef DEBUG
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif

		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=1;
		if (MPIrank == 0) std::cout << d << " " ;

		SetMatrixK_fromCSR ( cluster, d,
			K_mat[d].rows(), K_mat[d].columns(), //  .data[i]->KSparse->n_row,   clust_g.data[i]->KSparse->n_row,
			K_mat[d].rowPtrs(), K_mat[d].columnIndices(), K_mat[d].values(), //clust_g.data[i]->KSparse->row_ptr, clust_g.data[i]->KSparse->col_ind, clust_g.data[i]->KSparse->val,
			'G');
	}

	 timeSolKproc.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolKproc);

	if ( cluster.USE_KINV == 1 ) {
		 TimeEvent timeSolSC1(string("Solver - Schur Complement asm. - using solver"));
		 timeSolSC1.AddStart();
		cluster.Create_Kinv_perDomain();
		 timeSolSC1.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeSolSC1);

		 TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC"));
		 timeSolSC2.AddStart();
		cluster.Create_SC_perDomain();
		 timeSolSC2.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeSolSC2);
	}

	if (MPIrank == 0) std::cout << std::endl ;

	if (cluster.USE_HFETI == 1) {
		 TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		 timeHFETIprec.AddStart();
		cluster.SetClusterHFETI();
		 timeHFETIprec.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeHFETIprec);
	}

	 TimeEvent timeSolAkpl(string("Solver - Set Solver After Kplus"));
	 timeSolAkpl.AddStart();
	cluster.SetClusterPC_AfterKplus();
	 timeSolAkpl.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolAkpl);

	// *** END - Load Matrix K and regularization  ***********************************************************************


	// *** Running Solver ************************************************************************************************
	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 timeSolCG.AddStart();
    string result_file("MATSOL_SVN_Displacement.Nvec");
	solver.Solve_singular ( cluster, result_file );
	 timeSolCG.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolCG);

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.AddStart();
	vector < vector < double > > prim_solution;
	solver.GetSolution_Primal_singular_parallel(cluster, prim_solution);
	 timeGetSol.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeGetSol);

	double max_v = 0.0;
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++)
		for (ShortInt j = 0; j < prim_solution[i].size(); j++)
			if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.AddStartWOBarrier(0.0); max_sol_ev.AddEndWOBarrier(max_v);

	std::cout.precision(12);

	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	if (MPIrank == 0)
		std::cout << " Max value in_solution = " << max_vg << std::endl;

	max_sol_ev.PrintLastStatMPI_PerNode(max_vg);

	std::stringstream ss;
	ss << "mesh_" << MPIrank << ".vtk";

	 TimeEvent timeSaveVTK(string("Solver - Save VTK"));
	 timeSaveVTK.AddStart();
	input.mesh->saveVTK(ss.str().c_str(), prim_solution, l2g_vec, *input.localBoundaries, *input.globalBoundaries, 0.95, 0.9);
	 timeSaveVTK.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSaveVTK);

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


	timeEvalMain.totalTime.AddEndWithBarrier();
	timeEvalMain.PrintStatsMPI();



















}


void testBEM(int argc, char** argv)
{
	double start = omp_get_wtime();
	size_t partsCount = input.mesh->parts();
	size_t fixPointsCount = 4;

	std::cout << "1 : " << omp_get_wtime() - start << std::endl;

	mesh::SurfaceMesh sMesh(*input.mesh);

	std::cout << "2 : " << omp_get_wtime() - start << std::endl;

	sMesh.computeFixPoints(fixPointsCount);


	std::cout << "4 : " << omp_get_wtime() - start << std::endl;

	std::vector<DenseMatrix> K_mat_dense;

	K_mat_dense.reserve(partsCount);
	for (int d = 0; d < partsCount; d++) {
		K_mat_dense.push_back( DenseMatrix (0, 0) );
	}

    
#ifndef DEBUG
    cilk_for (size_t d = 0; d < partsCount; d++) {
#else
    for (size_t d = 0; d < partsCount; d++) {
#endif

/*
        std::ofstream Kmat_file_o;
        std::ofstream Kmat_file_p;
        
        Kmat_file_p.precision(15);
        Kmat_file_o.precision(15);
       4
        Kmat_file_o << std::scientific;
        Kmat_file_p << std::scientific;
        
        
        std::stringstream Kmat_file_name_o;
        Kmat_file_name_o << "Kmat_o_" << d;
        
        Kmat_file_o.open ( Kmat_file_name_o.str().c_str() );
        Kmat_file_p.open ("Kmat_p.txt");
        
 */
        DenseMatrix K_tmp;
        
        sMesh.elasticity(K_mat_dense[d], d);
        std::cout << d << " " << std::endl;
        
//        Kmat_file_o << K_mat_dense[d];
        
        int n = K_mat_dense[d].rows();
        
        K_tmp = K_mat_dense[d];
        
        for (int i = 0; i < n/3; i++) {
            for (int j = 0; j < n; j++) {
                K_tmp( 3*i+0,j) = K_mat_dense[d](0*(n/3) + i ,j);
                K_tmp( 3*i+1,j) = K_mat_dense[d](1*(n/3) + i ,j);
                K_tmp( 3*i+2,j) = K_mat_dense[d](2*(n/3) + i ,j);
            }
        }
 
        for (int i = 0; i < n/3; i++) {
            for (int j = 0; j < n; j++) {
                K_mat_dense[d]( j, 3*i+0) = K_tmp(j, 0*(n/3) + i );
                K_mat_dense[d]( j, 3*i+1) = K_tmp(j, 1*(n/3) + i );
                K_mat_dense[d]( j, 3*i+2) = K_tmp(j, 2*(n/3) + i );
            }
        }

/*        
        Kmat_file_p << K_mat_dense[d];
        
        Kmat_file_o.close();
        Kmat_file_p.close();
*/
        
    }

	std::cout << "5 : " << omp_get_wtime() - start << std::endl;


	// TODO:

    std::vector < SparseCSRMatrix<eslocal> >			K_mat;
    std::vector < SparseIJVMatrix<eslocal> >			B1_mat;
    std::vector < SparseIJVMatrix<eslocal> >			B0_mat;
    
    std::vector < std::vector <eslocal> >		lambda_map_sub_B1;
    std::vector < std::vector <eslocal> >		lambda_map_sub_B0;
    std::vector < std::vector <eslocal> >		lambda_map_sub_clst;
    std::vector < std::vector <double> >	B1_l_duplicity;
    
    std::vector < std::vector < double > >	f_vec     (partsCount);
    std::vector < std::vector < eslocal > >		fix_nodes (partsCount);
    std::vector < std::vector <eslocal> >		l2g_vec;
    
    std::cout << "BEM 8 : " << omp_get_wtime() - start<< std::endl;

    K_mat.reserve(partsCount);
    for (size_t d = 0; d < partsCount; d++) {
        K_mat.push_back( SparseCSRMatrix<eslocal> (0,0) );
    }

#ifndef DEBUG
    cilk_for (size_t d = 0; d < partsCount; d++) {
#else
    for (size_t d = 0; d < partsCount; d++) {
#endif
        K_mat[d] = K_mat_dense[d];
        f_vec[d].resize(K_mat_dense[d].rows() , 0.0);
    }
    K_mat_dense.clear();
    
    for (int d = 0; d < partsCount; d++) {
//      std::cout<< "d: "<< d <<std::endl;
      sMesh.integrateUpperFaces(f_vec[d],d);
    }
        
    
    std::cout << "9 : " << omp_get_wtime() - start<< std::endl;

    const std::vector<eslocal> fixPoints = sMesh.getFixPoints(); // input.mesh->getFixPoints();

#ifndef DEBUG
    cilk_for (eslocal d = 0; d < partsCount; d++) {
#else
    for (eslocal d = 0; d < partsCount; d++) {
#endif
        for (eslocal fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
            fix_nodes[d].push_back(fixPoints[d * fixPointsCount + fixPoint]);
        }
        std::sort ( fix_nodes[d].begin(), fix_nodes[d].end() );
    }

    std::cout << "11: " << omp_get_wtime() - start<< std::endl;
    input.localBoundaries->create_B1_l<eslocal>(
                            B1_mat,
                            B0_mat,
                            l2g_vec,
                            lambda_map_sub_clst,
                            lambda_map_sub_B1,
                            lambda_map_sub_B0,
                            B1_l_duplicity,
                            partsCount
                        );
//    for (int d = 0; d < partsCount; d++) {
//        for (int iz = 0; iz < l2g_vec[d].size(); iz++) {
//            if ( fabs( 30.0 - sMesh.coordinates()[l2g_vec[d][iz]].z ) < 0.00001 )
//                f_vec[d][3 * iz + 2] = 1.0;
//        }
//    }

    
    
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
                                     longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );
    
    extern void SetMatrixB0_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
                                     longInt n_rows, ShortInt n_cols, ShortInt nnz,
                                     longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );
    
    extern void SetMatrixR_fromDense( Cluster & cluster, ShortInt domain_index_in_cluster,
                                     ShortInt n_cols, ShortInt n_rows, double * vals, char type );
    
    extern void SetMatrixK_fromCSR ( Cluster & cluster, ShortInt domain_index_in_cluster,
                                    ShortInt n_rows, ShortInt n_cols, ShortInt * rows, ShortInt * cols, double * vals, char type );
    
    extern void SetMatrixK_fromBEM ( Cluster & cluster, ShortInt domain_index_in_cluster,
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
    solver.CG_max_iter	 = 1000;
    solver.USE_GGtINV	 = 1;
    solver.epsilon		 = 0.0000001;
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
        
#ifndef DEBUG
    cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
    for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
        ShortInt domain_index_in_cluster = i;
                
        SetMatrixB0_fromCOO( cluster, domain_index_in_cluster,
                B0_mat[i].rows(),			//clust_g.data[i]->B->B0_rows,		// B_full_rows, //n_row_eq,
                B0_mat[i].columns(),		//.data[i]->B->B0_cols,				// B_full_cols, //n_col,
                B0_mat[i].nonZeroValues(),	//.data[i]->B->B0_nnz,				// B_full_nnz,  //nnz_eq,
                B0_mat[i].rowIndices(),		//&clust_g.data[i]->B->B0_I[0],		// BI_full[0], //Bi_coo,
                B0_mat[i].columnIndices(),	//&clust_g.data[i]->B->B0_J[0],		// BJ_full[0], //Bj_coo,
                B0_mat[i].values(),			//&clust_g.data[i]->B->B0_V[0],		// BV_full[0], //Bv_coo,
                'G', B0_mat[i].indexing() );
        }
    }
    // *** END - Setup B0 matrix *************************************************************************************
        
    // *** Setup B1 matrix *******************************************************************************************
#ifndef DEBUG
    cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
    for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
        ShortInt domain_index_in_cluster = i;
        SetMatrixB1_fromCOO( cluster, domain_index_in_cluster,
                B1_mat[i].rows(),			//clust_g.data[i]->B->B_full_rows, //n_row_eq,
                B1_mat[i].columns(),		//clust_g.data[i]->B->B_full_cols, //n_col,
                B1_mat[i].nonZeroValues(),	//clust_g.data[i]->B->B_full_nnz,  //nnz_eq,
                B1_mat[i].rowIndices(),		//&clust_g.data[i]->B->BI_full[0], //Bi_coo,
                B1_mat[i].columnIndices(),	//&clust_g.data[i]->B->BJ_full[0], //Bj_coo,
                B1_mat[i].values(),			//&clust_g.data[i]->B->BV_full[0], //Bv_coo,
                'G', B1_mat[i].indexing() );
    }
            
#ifndef DEBUG
    cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
    for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
            cluster.domains[i].B1_scale_vec = B1_l_duplicity[i];
    }
    // *** END - Setup B1 matrix *************************************************************************************

                
                
    // *** Setup R matrix ********************************************************************************************
#ifndef DEBUG
    cilk_for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
    for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
        for (int i = 0; i < l2g_vec[d].size(); i++) {
            std::vector <double> tmp_vec (3,0);
            
            tmp_vec[0] = sMesh.coordinates()[l2g_vec[d][i]].x;
            tmp_vec[1] = sMesh.coordinates()[l2g_vec[d][i]].y;
            tmp_vec[2] = sMesh.coordinates()[l2g_vec[d][i]].z;
            cluster.domains[d].coordinates.push_back(tmp_vec);
        }
        cluster.domains[d].CreateKplus_R();
        //cluster.domains[d].Kplus_R.ConvertCSRToDense(0);
    }
    // *** END - Setup R matrix **************************************************************************************
                    
    // *** Load RHS and fix points for K regularization **************************************************************

    for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
        //SetVecDbl( cluster.domains[i].f,        clust_g.data[i]->KSparse->n_row, clust_g.data[i]->fE );
        cluster.domains[d].f = f_vec[d];
                            
        //SetVecInt( cluster.domains[i].fix_dofs, 1,                           24, clust_g.fem[i]->mesh.fixingDOFs );
        for (size_t i = 0; i < fix_nodes[d].size(); i++) {
            for (size_t d_i = 0; d_i < 3; d_i++) {
                cluster.domains[d].fix_dofs.push_back( 3 * fix_nodes[d][i] + d_i);
            }
        }
    }
    // *** END - Load RHS and fix points for K regularization ********************************************************
                        
    // *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
    for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
        cluster.domains[i].lambda_map_sub = lambda_map_sub_B1[i];
    }
                            
    std::vector < int > neigh_clusters;
    //neigh_clusters.push_back(0);
                            
    SetSolverPreprocessing ( cluster, solver, lambda_map_sub_clst, neigh_clusters );
    // *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************
                            
     
        
        
    // *** Load Matrix K and regularization ******************************************************************************
#ifndef DEBUG
    cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
    for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
        SetMatrixK_fromCSR ( cluster, d,
        K_mat[d].rows(), K_mat[d].columns(), //  .data[i]->KSparse->n_row,   clust_g.data[i]->KSparse->n_row,
        K_mat[d].rowPtrs(), K_mat[d].columnIndices(), K_mat[d].values(), //clust_g.data[i]->KSparse->row_ptr, clust_g.data[i]->KSparse->col_ind, clust_g.data[i]->KSparse->val,
                                                        'G');
    }
    K_mat.clear();

        
        
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
                                
    //input.mesh->saveVTK(prim_solution, l2g_vec);
                                
    sMesh.saveVTK("mesh.vtk", prim_solution, l2g_vec, *input.localBoundaries, *input.globalBoundaries, 0.95, 0.9);
    
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
        
                
    

}

void testFEM(int argc, char** argv)
{
	double start;
	start = omp_get_wtime();
	std::cout.precision(15);

	size_t partsCount = input.mesh->parts();
	size_t fixPointsCount = input.mesh->getFixPointsCount();


	std::cout << "5 : " << omp_get_wtime() - start<< std::endl;

	//Faces faces(mesh, coordinates);

	std::cout << "6 : " << omp_get_wtime() - start<< std::endl;

	//Corners corners(faces.getFaces(), coordinates);

	std::cout << "7 : " << omp_get_wtime() - start<< std::endl;

	std::vector < SparseCSRMatrix<eslocal> >			K_mat;
	std::vector < SparseCSRMatrix<eslocal> >			M_mat;
	std::vector < SparseIJVMatrix<eslocal> >			B1_mat;
	std::vector < SparseIJVMatrix<eslocal> >			B0_mat;

	std::vector < std::vector <eslocal> >		lambda_map_sub_B1;
	std::vector < std::vector <eslocal> >		lambda_map_sub_B0;
	std::vector < std::vector <eslocal> >		lambda_map_sub_clst;
	std::vector < std::vector <double> >	B1_l_duplicity;

	std::vector < std::vector < double > >	f_vec     (partsCount);
	std::vector < std::vector < eslocal > >		fix_nodes (partsCount);
	std::vector < std::vector <eslocal> >		l2g_vec;

	std::cout << "8 : " << omp_get_wtime() - start<< std::endl;

	K_mat.reserve(partsCount);
	M_mat.reserve(partsCount);
	for (eslocal d = 0; d < partsCount; d++) {
		K_mat.push_back( SparseCSRMatrix<eslocal> (0,0) );
		M_mat.push_back( SparseCSRMatrix<eslocal> (0,0) );
	}

	std::cout << "9 : " << omp_get_wtime() - start<< std::endl;

#ifndef DEBUG
	cilk_for (eslocal d = 0; d < partsCount; d++) {
#else
	for (eslocal d = 0; d < partsCount; d++) {
#endif
		eslocal dimension = input.mesh->getPartNodesCount(d) * mesh::Point::size();
		std::vector<double> f(dimension);

		input.mesh->elasticity(K_mat[d], M_mat[d], f, d);

		//K_mat[d] = K;
		//M_mat[d] = M;

        f_vec[d].swap(f);
        //f_vec[d].resize(K_mat[d].rows() , 0.0);

		std::cout << d << " " << std::endl;
	}
    //f_vec[partsCount-1][f_vec[partsCount-1].size() - 1] = 1.0;

        
        
        
        
	std::cout << "10: " << omp_get_wtime() - start<< std::endl;

	const std::vector<eslocal> fixPoints = input.mesh->getFixPoints();

#ifndef DEBUG
	cilk_for (eslocal d = 0; d < partsCount; d++) {
#else
	for (eslocal d = 0; d < partsCount; d++) {
#endif
		for (eslocal fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
			fix_nodes[d].push_back(fixPoints[d * fixPointsCount + fixPoint]);
		}
		std::sort ( fix_nodes[d].begin(), fix_nodes[d].end() );
	}

	std::cout << "11: " << omp_get_wtime() - start<< std::endl;
	input.localBoundaries->create_B1_l<eslocal>(
		B1_mat,
		B0_mat,
		l2g_vec,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		lambda_map_sub_B0,
		B1_l_duplicity,
		partsCount
	);

	//std::cout << B1_mat[0];

	const std::map<eslocal, double> &forces_x = input.mesh->coordinates().property(mesh::CP::FORCES_X).values();
	const std::map<eslocal, double> &forces_y = input.mesh->coordinates().property(mesh::CP::FORCES_Y).values();
	const std::map<eslocal, double> &forces_z = input.mesh->coordinates().property(mesh::CP::FORCES_Z).values();

	for (eslocal d = 0; d < partsCount; d++) {
		for (eslocal iz = 0; iz < l2g_vec[d].size(); iz++) {
			if (forces_x.find(l2g_vec[d][iz]) != forces_x.end()) {
				f_vec[d][3 * iz + 0] = forces_x.at(l2g_vec[d][iz]);
			}
			if (forces_y.find(l2g_vec[d][iz]) != forces_y.end()) {
				f_vec[d][3 * iz + 1] = forces_y.at(l2g_vec[d][iz]);
			}
			if (forces_z.find(l2g_vec[d][iz]) != forces_z.end()) {
				f_vec[d][3 * iz + 2] = forces_z.at(l2g_vec[d][iz]);
			}
		}
	}
//    for (eslocal d = 0; d < partsCount; d++) {
//        for (eslocal iz = 0; iz < l2g_vec[d].size(); iz++) {
//            if ( fabs( 30.0 - input.mesh->coordinates()[l2g_vec[d][iz]].z ) < 0.00001 ) {
//                //f_vec[d][3 * iz + 2] = 1.0;
//            }
//        }
//    }

/*
	for (eslocal d = 0; d < partsCount; d++) {
    // K
    SparseIJVMatrix tmpK = K_mat[d];
    std::ofstream (K_mat_file);
    K_mat_file.precision(15);
    K_mat_file << std::scientific;
    std::stringstream K_mat_file_name;
    K_mat_file_name << "dumped_files/K_mat_" << d;
    K_mat_file.open ( K_mat_file_name.str().c_str() );
    K_mat_file << tmpK;
    K_mat_file.close();
    // f
    std::ofstream (f_vec_file);
    f_vec_file.precision(15);
    f_vec_file << std::scientific;
    std::stringstream f_vec_file_name;
    f_vec_file_name << "dumped_files/f_vec_" << d;
    f_vec_file.open ( f_vec_file_name.str().c_str() );
    f_vec_file << f_vec[d];
    f_vec_file.close();
    // B
    std::ofstream (B0_mat_file);
    B0_mat_file.precision(15);
    B0_mat_file << std::scientific;
    std::stringstream B0_mat_file_name;
    B0_mat_file_name << "dumped_files/B0_mat_" << d;
    B0_mat_file.open ( B0_mat_file_name.str().c_str() );
    B0_mat_file << B0_mat[d];
    B0_mat_file.close();
  }
*/



        
        
	std::cout << "12: " << omp_get_wtime() - start<< std::endl;

	std::cout.precision(10);

	// Start - Stupid version of ESPRESO interface

	//MPI_Init (&argc, &argv);					// starts MPI

	typedef int       ShortInt ;
	typedef int       longInt  ;


	int MPIrank = 0; //MPI_Comm_rank(fem->comm, &MPIrank);
	int MPIsize = 1; //MPI_Comm_size(fem->comm, &MPIsize);
	int number_of_subdomains_per_cluster = partsCount;


	extern void SetCluster		  ( Cluster & cluster, ShortInt * subdomains_global_indices, ShortInt number_of_subdomains, ShortInt MPI_rank);

	extern void SetMatrixB1_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );

	extern void SetMatrixB0_fromCOO ( Cluster & cluster, ShortInt domain_index_in_cluster,
		longInt n_rows, ShortInt n_cols, ShortInt nnz,
		longInt * I_rows, ShortInt * J_cols, double * V_vals, char type, int indexing );

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
	solver.CG_max_iter	 = 1000;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = 0.00001;
	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = 0;
	solver.USE_PIPECG	 = 0;
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

#ifndef DEBUG
		cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
		for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
			ShortInt domain_index_in_cluster = i;

			SetMatrixB0_fromCOO( cluster, domain_index_in_cluster,
				B0_mat[i].rows(),			//clust_g.data[i]->B->B0_rows,		// B_full_rows, //n_row_eq,
				B0_mat[i].columns(),		//.data[i]->B->B0_cols,				// B_full_cols, //n_col,
				B0_mat[i].nonZeroValues(),	//.data[i]->B->B0_nnz,				// B_full_nnz,  //nnz_eq,
				B0_mat[i].rowIndices(),		//&clust_g.data[i]->B->B0_I[0],		// BI_full[0], //Bi_coo,
				B0_mat[i].columnIndices(),	//&clust_g.data[i]->B->B0_J[0],		// BJ_full[0], //Bj_coo,
				B0_mat[i].values(),			//&clust_g.data[i]->B->B0_V[0],		// BV_full[0], //Bv_coo,
				'G', B0_mat[i].indexing() );
		}
	}
	// *** END - Setup B0 matrix *************************************************************************************

	// *** Setup B1 matrix *******************************************************************************************
#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		ShortInt domain_index_in_cluster = i;
		SetMatrixB1_fromCOO( cluster, domain_index_in_cluster,
			B1_mat[i].rows(),			//clust_g.data[i]->B->B_full_rows, //n_row_eq,
			B1_mat[i].columns(),		//clust_g.data[i]->B->B_full_cols, //n_col,
			B1_mat[i].nonZeroValues(),	//clust_g.data[i]->B->B_full_nnz,  //nnz_eq,
			B1_mat[i].rowIndices(),		//&clust_g.data[i]->B->BI_full[0], //Bi_coo,
			B1_mat[i].columnIndices(),	//&clust_g.data[i]->B->BJ_full[0], //Bj_coo,
			B1_mat[i].values(),			//&clust_g.data[i]->B->BV_full[0], //Bv_coo,
			'G', B1_mat[i].indexing() );
	}

#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		cluster.domains[i].B1_scale_vec = B1_l_duplicity[i];
	}
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
#ifndef DEBUG
	cilk_for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for(ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
		for (int i = 0; i < l2g_vec[d].size(); i++) {
			std::vector <double> tmp_vec (3,0);
			tmp_vec[0] = input.mesh->coordinates()[l2g_vec[d][i]].x;
			tmp_vec[1] = input.mesh->coordinates()[l2g_vec[d][i]].y;
			tmp_vec[2] = input.mesh->coordinates()[l2g_vec[d][i]].z;
			cluster.domains[d].coordinates.push_back(tmp_vec);
		}
		cluster.domains[d].CreateKplus_R();
		//cluster.domains[d].Kplus_R.ConvertCSRToDense(0);
	}
	// *** END - Setup R matrix **************************************************************************************

	// *** Load RHS and fix points for K regularization **************************************************************
#ifndef DEBUG
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif
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
#ifndef DEBUG
	cilk_for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#else
	for (ShortInt i = 0; i < number_of_subdomains_per_cluster; i++) {
#endif
		cluster.domains[i].lambda_map_sub = lambda_map_sub_B1[i];
	}

	std::vector < int > neigh_clusters;
	//neigh_clusters.push_back(0);

	SetSolverPreprocessing ( cluster, solver, lambda_map_sub_clst, neigh_clusters );
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// *** Load Matrix K and regularization ******************************************************************************
#ifndef DEBUG
	cilk_for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#else
	for (ShortInt d = 0; d < number_of_subdomains_per_cluster; d++) {
#endif

		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=1;

		SetMatrixK_fromCSR ( cluster, d,
			K_mat[d].rows(), K_mat[d].columns(), //  .data[i]->KSparse->n_row,   clust_g.data[i]->KSparse->n_row,
			K_mat[d].rowPtrs(), K_mat[d].columnIndices(), K_mat[d].values(), //clust_g.data[i]->KSparse->row_ptr, clust_g.data[i]->KSparse->col_ind, clust_g.data[i]->KSparse->val,
			'G');
	}

	//std::cout << std::endl;

	//cluster.Create_Kinv_perDomain();

	//cluster.Create_SC_perDomain();


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

	input.mesh->saveVTK("mesh.vtk", prim_solution, l2g_vec, *input.localBoundaries, *input.globalBoundaries, 1.0, 0.95);

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
}




