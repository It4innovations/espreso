/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#include "LinearSolver.h"

//#include <Eigen/Dense>
//using Eigen::MatrixXd;

using namespace espreso;

LinearSolver::LinearSolver(Instance *instance, const ESPRESOSolver &configuration)
: instance(instance),
  configuration(configuration),
  physics(NULL),
  constraints(NULL),
  timeEvalMain("ESPRESO Solver Overal Timing"),
  cluster(NULL),
  solver(NULL)
{
}


LinearSolver::LinearSolver(const ESPRESOSolver &configuration, OldPhysics &physics, Constraints &constraints)
: instance(NULL),
  configuration(configuration),
  physics(&physics),
  constraints(&constraints),
  timeEvalMain("ESPRESO Solver Overal Timing")
{
	cluster  = new Cluster(configuration, instance);
	solver = new IterSolver(configuration);

//	numClusters = 1 + *std::max_element(instance->clustersMap.begin(), instance->clustersMap.end());
//
//	clusters.resize(numClusters);
//	for (eslocal c = 0; c < clusters.size(); c++) {
//		clusters[c] = new Cluster(configuration, instance);
//	}

}

LinearSolver::~LinearSolver() {
	delete cluster;

	for (eslocal c = 0; c < clusters.size(); c++) {
		delete clusters[c];
	}

	delete solver;
}


// make full initialization of solver
void LinearSolver::init()
{
	if (cluster != NULL) {
		delete cluster;
		for (eslocal c = 0; c < clusters.size(); c++) {
			delete clusters[c];
		}
	}

	if (solver != NULL) {
		delete solver;
	}

	cluster = new Cluster(configuration, instance);
	solver  = new IterSolver(configuration);

	numClusters = 1 + *std::max_element(instance->clustersMap.begin(), instance->clustersMap.end());

	clusters.resize(numClusters);
	for (eslocal c = 0; c < clusters.size(); c++) {
		clusters[c] = new Cluster(configuration, instance);
	}

	solver->numClusters = numClusters;
	solver->clusters = &clusters;

	init(instance->neighbours);


}

// make partial initialization according to updated matrices
void LinearSolver::update(Matrices matrices)
{
	// TODO update appropriate solver objects and stop steeling matrices! :)

	if (matrices & Matrices::K) {
		// factorization and preconditioners and HFETI preprocessing

		delete cluster;
		for (eslocal c = 0; c < clusters.size(); c++) {
			delete clusters[c];
		}

		delete solver;

		cluster = new Cluster(configuration, instance);
		solver  = new IterSolver(configuration);

		for (eslocal c = 0; c < clusters.size(); c++) {
			clusters[c] = new Cluster(configuration, instance);
		}

		solver->numClusters = numClusters;
		solver->clusters = &clusters;

		init(instance->neighbours);

	}

	if (matrices & (Matrices::N | Matrices::B1)) { // N is kernel of matrix K
		// updateGGt();
		setup_CreateG_GGt_CompressG();
	}

	if (matrices & (Matrices::B1c | Matrices::f)) {
		// update dual RHS
		//for f:
		// - updated by solve

		//for B1c:
		setup_SetDirichletBoundaryConditions();
	}

	if (matrices & Matrices::B0) {
		// HFETI preprocessing
		setup_HTFETI();
	}

	if (matrices & Matrices::B1duplicity) {
		// update duplicity vector
		#pragma omp parallel for
		for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
			cluster->domains[d].B1_scale_vec = instance->B1duplicity[d];
		}
	}

	// TODO: remove full re-initialization
//	delete cluster;
//	delete solver;
//	cluster = new Cluster(configuration);
//	solver = new IterSolver(configuration);
//	setup();
//	init(instance->neighbours);
}

// run solver and store primal and dual solution
void LinearSolver::run()
{
	Solve(instance->f, instance->primalSolution, instance->dualSolution);
}


void LinearSolver::setup_HTFETI() {
// Setup Hybrid FETI part of the solver

	if (cluster->USE_HFETI == 1) {
			TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
			timeHFETIprec.start();

			if (numClusters == 1) {
				cluster->SetClusterHFETI();
			} else {
				for (eslocal c = 0; c < clusters.size(); c++) {
					clusters[c]->SetClusterHFETI();
				}
			}

			timeHFETIprec.endWithBarrier();
			timeEvalMain.addEvent(timeHFETIprec);

			ESLOG(MEMORY) << "After HFETI preprocessing process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
			ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	}
}

void LinearSolver::setup_LocalSchurComplement() {
// Computation of the Local Schur Complements

	if (cluster->USE_KINV == 1) {
			TimeEvent KSCMem(string("Solver - SC asm. w PARDISO-SC mem [MB]")); KSCMem.startWithoutBarrier(GetProcessMemory_u());
			TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC"));timeSolSC2.start();
		bool USE_FLOAT = false;
		if (configuration.schur_precision == FLOAT_PRECISION::SINGLE) {
			USE_FLOAT = true;
		}
		cluster->Create_SC_perDomain(USE_FLOAT);
			timeSolSC2.endWithBarrier(); timeEvalMain.addEvent(timeSolSC2);
			KSCMem.endWithoutBarrier(GetProcessMemory_u()); //KSCMem.printLastStatMPIPerNode();
			ESLOG(MEMORY) << "After K inv. process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
			ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	} else {
		for (size_t d = 0; d < cluster->domains.size(); d++) {
			cluster->domains[d].isOnACC = 0;
		}
	}

}

void LinearSolver::setup_Preconditioner() {
// Load Matrix K, Regularization
		TimeEvent timeRegKproc(string("Solver - Setup preconditioners")); timeRegKproc.start();
		ESLOG(MEMORY) << "Before - Setup preconditioners - process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent KregMem(string("Solver - Setup preconditioners mem. [MB]")); KregMem.startWithoutBarrier(GetProcessMemory_u());
		ESINFO(PROGRESS3) << "Setup preconditioners";

	if (numClusters == 1) {
		cluster->SetupPreconditioner();
	} else {
		for (eslocal c = 0; c < clusters.size(); c++) {
			clusters[c]->SetupPreconditioner();
		}
	}

		KregMem.endWithoutBarrier(GetProcessMemory_u()); //KregMem.printLastStatMPIPerNode();

		ESLOG(MEMORY) << "After - Setup preconditioners " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		timeRegKproc.endWithBarrier();
		timeEvalMain.addEvent(timeRegKproc);
}

void LinearSolver::setup_FactorizationOfStiffnessMatrices() {
// K Factorization
		TimeEvent timeSolKproc(string("Solver - K factorization")); timeSolKproc.start();
		TimeEvent KFactMem(string("Solver - K factorization mem. [MB]")); KFactMem.startWithoutBarrier(GetProcessMemory_u());
		ESINFO(PROGRESS3) << "Factorize K";

	if (numClusters == 1) {
		cluster->SetupKsolvers();
	} else {
		for (eslocal c = 0; c < clusters.size(); c++) {
			clusters[c]->SetupKsolvers();
		}
	}

		KFactMem.endWithoutBarrier(GetProcessMemory_u()); //KFactMem.printLastStatMPIPerNode();
		ESLOG(MEMORY) << "After K solver setup process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		timeSolKproc.endWithBarrier();
		timeEvalMain.addEvent(timeSolKproc);
}

void LinearSolver::setup_KernelMatrices() {
// Setup R matrix
//	if (SINGULAR) {
//		TimeEvent timeSetR(string("Solver - Set R"));
//		timeSetR.start();
//		for (int d = 0; d < number_of_subdomains_per_cluster; d++) {
//			// physics.R1[d].GramSchmidtOrtho();
//			cluster->domains[d].Kplus_R = instance->N1[d];
//			cluster->domains[d].Kplus_R2 = instance->N2[d];
//			cluster->domains[d].Kplus_Rb = instance->N1[d];
//			cluster->domains[d].Kplus_Rb2 = instance->N2[d];
//		}
//		timeSetR.endWithBarrier();
//		timeEvalMain.addEvent(timeSetR);
//	}


}

void LinearSolver::setup_B1Matrices() {
// Setup B1 matrix
//	TimeEvent timeSetB1(string("Solver - Set B1"));
//	timeSetB1.start();
//
//	#pragma omp parallel for
//	for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//
//		cluster->domains[d].B1 = instance->B1[d];
//		cluster->domains[d].B1.type = 'G';
//
//		cluster->domains[d].B1t = cluster->domains[d].B1;
//		cluster->domains[d].B1t.MatTransposeCOO();
//		cluster->domains[d].B1t.ConvertToCSRwithSort(1);
//
//	}
//
//	#pragma omp parallel for
//	for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		cluster->domains[d].B1_scale_vec = instance->B1duplicity[d];
//	}
//
//
//	timeSetB1.endWithBarrier();
//	timeEvalMain.addEvent(timeSetB1);
}

void LinearSolver::setup_B0Matrices() {
// Setup B0 matrix

//	if (cluster->USE_HFETI == 1) {
//		TimeEvent timeSetB0(string("Solver - Set B0"));
//		timeSetB0.start();
//
//		#pragma omp parallel for
//		for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//			cluster->domains[d].B0 = instance->B0[d];
//			cluster->domains[d].B0.type = 'G';
//			cluster->domains[d].B0.ConvertToCSRwithSort(1);
//		}
//
//		timeSetB0.endWithBarrier();
//		timeEvalMain.addEvent(timeSetB0);
//	}
}

void LinearSolver::setup_SetDirichletBoundaryConditions() {
// Set Dirichlet Boundary Condition

		TimeEvent timeSetInitialCondition(string("Solver - Set Dirichlet Boundary Condition"));
		timeSetInitialCondition.start();

	#pragma omp parallel for
	for (int d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster->domains[d].vec_c = instance->B1c[d];
		cluster->domains[d].vec_lb = instance->LB[d];
	}

		timeSetInitialCondition.endWithBarrier();
		timeEvalMain.addEvent(timeSetInitialCondition);
}

void LinearSolver::setup_CreateDirichletPreconditioner() {
// Calculate Dirichlet Preconditioner
//	if (    configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET
//		 || configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET) {
//
//			TimeEvent timeDirPrec( string("Solver - Dirichlet Preconditioner calculation")); timeDirPrec.start();
//			ESINFO(PROGRESS3) << "Calculate Dirichlet preconditioner";
//		cluster->CreateDirichletPrec(instance);
//			ESINFO(PROGRESS3);
//			timeDirPrec.endWithBarrier();
//			timeEvalMain.addEvent(timeDirPrec);
//	}
}

void LinearSolver::setup_CreateG_GGt_CompressG() {

		TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.start();

		ESLOG(MEMORY) << "Solver Preprocessing - HFETI with regularization from K matrix";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses "<< Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent G1_perCluster_time("Setup G1 per Cluster time   - preprocessing"); G1_perCluster_time.start();
		TimeEvent G1_perCluster_mem("Setup G1 per Cluster memory - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());


		// *** Multiple clusters support ********************************************
		if (numClusters == 1) {
			cluster->Create_G_perCluster();
		} else {
			for (eslocal c = 0; c < numClusters; c++) {
				clusters[c]->Create_G_perCluster();
				cluster->G1.MatAppend(clusters[c]->G1);
				cluster->G2.MatAppend(clusters[c]->G2);
			}
		}

		G1_perCluster_time.end(); G1_perCluster_time.printStatMPI();
		G1_perCluster_mem.endWithoutBarrier(GetProcessMemory_u()); G1_perCluster_mem.printStatMPI();

		ESLOG(MEMORY) << "Created G1 per cluster";
		ESLOG(MEMORY) << "Before HFETI create GGt process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent solver_Preprocessing_time("Setup GGt time   - preprocessing"); solver_Preprocessing_time.start();
		TimeEvent solver_Preprocessing_mem("Setup GGt memory - preprocessing"); solver_Preprocessing_mem.start();
	solver->Preprocessing(*cluster);
		solver_Preprocessing_time.end(); solver_Preprocessing_time.printStatMPI();
		solver_Preprocessing_mem.end();  solver_Preprocessing_mem.printStatMPI();

		ESLOG(MEMORY) << "Create GGtInv";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent solver_G1comp_time("Setup G1 compression time   - preprocessing"); solver_G1comp_time.start();
		TimeEvent solver_G1comp_mem("Setup G1 compression memory - preprocessing");  solver_G1comp_mem.start();

		if (numClusters == 1) {
			cluster->Compress_G1(); // Compression of Matrix G1 to work with compressed lambda vectors
		} else {
			cluster->Compress_G1(); // Compression of Matrix G1 to work with compressed lambda vectors
			for (eslocal c = 0; c < numClusters; c++) {
				clusters[c]->Compress_G1(); // Compression of Matrix G1 to work with compressed lambda vectors
			}
		}

		solver_G1comp_time.end(); solver_G1comp_time.printStatMPI();
		solver_G1comp_mem.end();  solver_G1comp_mem.printStatMPI();
		ESLOG(MEMORY) << "G1 compression";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		timeSolPrec.endWithBarrier(); timeEvalMain.addEvent(timeSolPrec);


}

void LinearSolver::setup_SetupCommunicationLayer() {
//
//		ESLOG(MEMORY) << "Preprocessing setup comm. layer - start";
//		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
//		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
//
//		TimeEvent cluster_SetClusterPC_time("Setup communication layer time   - pre-processing"); cluster_SetClusterPC_time.start();
//		TimeEvent cluster_SetClusterPC_mem ("Setup communication layer memory - pre-processing"); cluster_SetClusterPC_mem .start();
//
//	#pragma omp parallel for
//	for (int d = 0; d < number_of_subdomains_per_cluster; d++) {
//		cluster->domains[d].lambda_map_sub = instance->B1subdomainsMap[d];
//	}
//
//	cluster->SetClusterPC(); //instance->B1clustersMap); //lambda_map_sub
//
//		cluster_SetClusterPC_time.end(); cluster_SetClusterPC_time.printStatMPI();
//		cluster_SetClusterPC_mem .end(); cluster_SetClusterPC_mem .printStatMPI();
//
//		ESLOG(MEMORY) << "Preprocessing setup comm. layer - stop";
//		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
//		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
}

void LinearSolver::setup_InitClusterAndSolver( )
{
// Setup Cluster and Solver

	TimeEvent timeSetClust(string("Solver - Set cluster")); timeSetClust.start();


	// ***************************************************************************************************************************
	// Cluster structure  setup
	cluster->USE_DYNAMIC = 0;
	SINGULAR = true; // TODO: refactor
	if (instance != NULL) {
		cluster->USE_DYNAMIC = 1;
		SINGULAR = false; // TODO: refactor
		for (size_t d = 0; d < instance->domains; d++) {
			if (instance->N1[d].cols) {
				cluster->USE_DYNAMIC = 0;
				SINGULAR = true; // TODO: refactor
				break;
			}
		}
	}
	switch (configuration.method) {
	case ESPRESO_METHOD::TOTAL_FETI:
		cluster->USE_HFETI = false;
		break;
	case ESPRESO_METHOD::HYBRID_FETI:
		cluster->USE_HFETI = true;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
	}
	cluster->USE_KINV = configuration.use_schur_complement ? 1 : 0;
	cluster->SUBDOM_PER_CLUSTER = number_of_subdomains_per_cluster;
	cluster->NUMBER_OF_CLUSTERS = environment->MPIsize;
	cluster->PAR_NUM_THREADS = environment->PAR_NUM_THREADS;
	cluster->SOLVER_NUM_THREADS = environment->SOLVER_NUM_THREADS;


	eslocal GUSE_DYNAMIC = 0;
	int global_numClusters;

	MPI_Allreduce(&numClusters, &global_numClusters, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// *** Multiple clusters support ********************************************
	for (eslocal c = 0; c < numClusters; c++) {
		clusters[c]->USE_DYNAMIC = 0;
		SINGULAR = true; // TODO: refactor
		if (instance != NULL) {
			clusters[c]->USE_DYNAMIC = 1;
			SINGULAR = false; // TODO: refactor
			for (size_t d = 0; d < instance->domains; d++) {
				if (instance->N1[d].cols) {
					clusters[c]->USE_DYNAMIC = 0;
					SINGULAR = true; // TODO: refactor
					break;
				}
			}
		}

		if (clusters[c]->USE_DYNAMIC == 1)
			GUSE_DYNAMIC = 1;

		switch (configuration.method) {
		case ESPRESO_METHOD::TOTAL_FETI:
			clusters[c]->USE_HFETI = false;
			break;
		case ESPRESO_METHOD::HYBRID_FETI:
			clusters[c]->USE_HFETI = true;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
		}
		clusters[c]->USE_KINV = configuration.use_schur_complement ? 1 : 0;
		clusters[c]->SUBDOM_PER_CLUSTER = number_of_subdomains_per_cluster;

		clusters[c]->NUMBER_OF_CLUSTERS = global_numClusters; //environment->MPIsize;				//TODO: MPC Fix - not true anymore

		clusters[c]->PAR_NUM_THREADS = environment->PAR_NUM_THREADS;
		clusters[c]->SOLVER_NUM_THREADS = environment->SOLVER_NUM_THREADS;
	}
	// *** END - Multiple clusters support ********************************************

	// ***************************************************************************************************************************


	// ***************************************************************************************************************************
	// Iter Solver Set-up
	solver->CG_max_iter = configuration.iterations;
	solver->USE_GGtINV = 1;
	solver->epsilon = configuration.epsilon;
	solver->USE_PREC = configuration.preconditioner;

	//solver->USE_HFETI = cluster->USE_HFETI;
	switch (configuration.method) {
	case ESPRESO_METHOD::TOTAL_FETI:
		solver->USE_HFETI = false;
		break;
	case ESPRESO_METHOD::HYBRID_FETI:
		solver->USE_HFETI = true;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
	}

	solver->USE_KINV = configuration.use_schur_complement ? 1 : 0; // cluster->USE_KINV;
	solver->USE_DYNAMIC = GUSE_DYNAMIC; //TODO: flag per MPI process // cluster->USE_DYNAMIC;

	solver->PAR_NUM_THREADS = environment->PAR_NUM_THREADS;
	solver->SOLVER_NUM_THREADS = environment->SOLVER_NUM_THREADS;
	// ***************************************************************************************************************************



	number_of_subdomains_per_cluster = instance->K.size();
	std::vector<eslocal> domain_list(number_of_subdomains_per_cluster, 0);
	for (int i = 0; i < number_of_subdomains_per_cluster; i++) {
		domain_list[i] = i;
	}
	cluster->cluster_global_index = environment->MPIrank + 1;
	cluster->my_neighs = std::vector<eslocal>(instance->neighbours.begin(), instance->neighbours.end());
	cluster->mtype = instance->K[0].mtype; // TODO: refactor

	switch (cluster->mtype) {
	case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		cluster->SYMMETRIC_SYSTEM = true;
		break;
	case MatrixType::REAL_SYMMETRIC_INDEFINITE:
		cluster->SYMMETRIC_SYSTEM = true;
		break;
	case MatrixType::REAL_UNSYMMETRIC:
		cluster->SYMMETRIC_SYSTEM = false;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown matrix type";
	}

	if (!cluster->SYMMETRIC_SYSTEM
			&& (configuration.solver != ESPRESO_ITERATIVE_SOLVER::GMRES
					&& configuration.solver
							!= ESPRESO_ITERATIVE_SOLVER::BICGSTAB)) {
		ESINFO(GLOBAL_ERROR)
				<< "Only GMRES or BICGSTAB solvers can solve the non-symmetric systems.";
	}

	cluster->InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);

	timeSetClust.endWithBarrier(); timeEvalMain.addEvent(timeSetClust);


		ESLOG(MEMORY) << "Preprocessing setup comm. layer - start";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent cluster_SetClusterPC_time("Setup communication layer time   - pre-processing"); cluster_SetClusterPC_time.start();
		TimeEvent cluster_SetClusterPC_mem ("Setup communication layer memory - pre-processing"); cluster_SetClusterPC_mem .start();

	cluster->SetClusterPC(); //instance->B1clustersMap); //lambda_map_sub

		cluster_SetClusterPC_time.end(); cluster_SetClusterPC_time.printStatMPI();
		cluster_SetClusterPC_mem .end(); cluster_SetClusterPC_mem .printStatMPI();

		ESLOG(MEMORY) << "Preprocessing setup comm. layer - stop";
		ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";



	// *** Multiple clusters support ********************************************
	int glob_clust_index = 0;
	MPI_Exscan(&numClusters, &glob_clust_index, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	for (eslocal c = 0; c < numClusters; c++) {

		eslocal number_of_subdomains_per_cluster = 0;
		std::vector<eslocal> domain_list;
		for (eslocal i = 0; i < instance->clustersMap.size(); i++) {
			if (instance->clustersMap[i] == c) {
				number_of_subdomains_per_cluster++;
				domain_list.push_back(i);
			}
		}
		clusters[c]->cluster_global_index = glob_clust_index + c + 1; // environment->MPIrank + 1;
		clusters[c]->my_neighs = std::vector<eslocal>(instance->neighbours.begin(), instance->neighbours.end()); // TODO: musim kopirovat - nedokaze ondra dodat jinak ?
		clusters[c]->mtype = instance->K[0].mtype; // TODO: opravit

		switch (clusters[c]->mtype) {
		case MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
			clusters[c]->SYMMETRIC_SYSTEM = true;
			break;
		case MatrixType::REAL_SYMMETRIC_INDEFINITE:
			clusters[c]->SYMMETRIC_SYSTEM = true;
			break;
		case MatrixType::REAL_UNSYMMETRIC:
			clusters[c]->SYMMETRIC_SYSTEM = false;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown matrix type";
		}

		if (!clusters[c]->SYMMETRIC_SYSTEM
				&& (configuration.solver != ESPRESO_ITERATIVE_SOLVER::GMRES
						&& configuration.solver
								!= ESPRESO_ITERATIVE_SOLVER::BICGSTAB)) {
			ESINFO(GLOBAL_ERROR)
					<< "Only GMRES or BICGSTAB solvers can solve the non-symmetric systems.";
		}

		clusters[c]->InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);
		clusters[c]->SetClusterPC();

	}
	// *** END - Multiple clusters support ********************************************


}

// TODO: const parameters
void LinearSolver::init(const std::vector<int> &neighbours)
{
	if (physics != NULL) {
		instance = new Instance(physics->K.size(), neighbours);

		physics->K.swap(instance->K);
		for (size_t i = 0; i < physics->K.size(); i++) {
			instance->K[i].mtype = physics->mtype;
		}

		physics->R1.swap(instance->N1);
		physics->R2.swap(instance->N2);
		physics->RegMat.swap(instance->RegMat);

		constraints->B1.swap(instance->B1);
		constraints->B1c.swap(instance->B1c);
		constraints->B1subdomainsMap.swap(instance->B1subdomainsMap);
		constraints->B1clustersMap.swap(instance->B1clustersMap);
		constraints->B1duplicity.swap(instance->B1duplicity);

		constraints->B0.swap(instance->B0);
		constraints->B0subdomainsMap.swap(instance->B0subdomainsMap);

		constraints->LB.swap(instance->LB);
		constraints->inequality.swap(instance->inequality);
		constraints->inequalityC.swap(instance->inequalityC);

		cluster->instance = instance;

	}

	//mkl_cbwr_set(MKL_CBWR_COMPATIBLE);

	// Overall Linear Solver Time measurement structure
	timeEvalMain.totalTime.startWithBarrier();



	// *** Initialize the Cluster and Solver structures
	setup_InitClusterAndSolver();


	// *** Setup B0 matrix
	// setup_B0Matrices();

	// *** Setup B1 matrix
	// setup_B1Matrices();

	// *** Setup R matrix
	// setup_KernelMatrices();

	// *** Set Dirichlet Boundary Condition
	// setup_SetDirichletBoundaryConditions();

	// *** Setup the communication layer using lambda_map_sub vectors
	// ***TODO: Moved to "setup_InitClusterAndSolver"
	// setup_SetupCommunicationLayer();



	// *** if TFETI is used or if HTFETI and analytical kernel are used we can compute GGt here - between solution in terms of peak memory
	if (  !(cluster->USE_HFETI == 1 && configuration.regularization == REGULARIZATION::NULL_PIVOTS)  ) {
		setup_CreateG_GGt_CompressG();
	}

	// *** setup all preconditioners
	setup_Preconditioner();

	// *** Computation of the Schur Complement
	setup_LocalSchurComplement();

	// *** K Factorization
	setup_FactorizationOfStiffnessMatrices();

	// *** Setup Hybrid FETI part of the solver
	setup_HTFETI();

	// *** For algebraic kernels, GGt needs to be computed after HTFETI preprocessing
	if (cluster->USE_HFETI == 1 && configuration.regularization == REGULARIZATION::NULL_PIVOTS) {
		setup_CreateG_GGt_CompressG();
	}


	// Cleanup of unnecessary objects
	// TODO: This can be a problem in some cases - need to be verified
	cluster->_my_lamdas_map_indices.clear();

	#pragma omp parallel for
	for (size_t d = 0; d < cluster->domains.size(); d++) {
		cluster->domains[d].B1.Clear();
	}

	ESLOG(MEMORY) << "End of preprocessing - process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution)
{

	std::vector < std::vector < double > > dual_solution;
	Solve(f_vec, prim_solution, dual_solution);
}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution,
		                  std::vector < std::vector < double > > & dual_solution) {

	 	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 	 timeSolCG.start();

	 solver->Solve_singular    ( *cluster, f_vec, prim_solution, dual_solution );

	 	 timeSolCG.endWithBarrier();
	 	 timeEvalMain.addEvent(timeSolCG);

}

void LinearSolver::Postprocessing( ) {

}

void LinearSolver::finilize() {

	// Show Linear Solver Runtime Evaluation
	solver->preproc_timing.printStatsMPI();
	solver->timing.printStatsMPI();
	solver->postproc_timing.printStatsMPI();
	solver->timeEvalAppa.printStatsMPI();
	solver->timeEvalProj.printStatsMPI();

	if ( solver->USE_PREC != ESPRESO_PRECONDITIONER::NONE ) solver->timeEvalPrec.printStatsMPI();

	if ( cluster->USE_HFETI == 1 ) cluster->ShowTiming();

	 timeEvalMain.totalTime.endWithBarrier();
	 timeEvalMain.printStatsMPI();
}

void LinearSolver::CheckSolution( vector < vector < double > > & prim_solution ) {
    // *** Solutin correctnes test **********************************************************************************************
	double max_v = 0.0;
		for (eslocal i = 0; i < number_of_subdomains_per_cluster; i++)
			for (size_t j = 0; j < prim_solution[i].size(); j++)
				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.startWithoutBarrier(0.0); max_sol_ev.endWithoutBarrier(max_v);

	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	ESINFO(DETAILS) << "Maxvalue in solution = " << std::setprecision(12) << max_vg;

	//max_sol_ev.printLastStatMPIPerNode(max_vg);
	// *** END - Solutin correctnes test ******************************************************************************************

}

void LinearSolver::Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub) {



}
