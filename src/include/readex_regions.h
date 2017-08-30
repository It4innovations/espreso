#ifndef READEX_REGIONS_H_
#define READEX_REGIONS_H_

//	src/app/espreso.cpp
READEX_REGION_DEFINE(REG_Main);

//	src/app/factory/factory.cpp
READEX_REGION_DEFINE(REG_instance_solve);
READEX_REGION_DEFINE(REG_instance_PrepareMesh); // renamed from REG_Assembler_PrepareMesh

//src/assembler/physicssolver/assembler.cpp
//READEX_REGION_DEFINE(REG_Assembler_SaveMeshtoVTK);	// removed
READEX_REGION_DEFINE(REG_Assembler_SolverSolve);
READEX_REGION_DEFINE(REG_Assembler_SaveResults);
READEX_REGION_DEFINE(REG_Assembler_AssembleStiffnesMatrices);
READEX_REGION_DEFINE(REG_Assembler_SaveMatricesToFile);
READEX_REGION_DEFINE(REG_Assembler_SolverInit);
//READEX_REGION_DEFINE(REG_Assembler_K_Regularization);	//solver ?

//	src/assembler/physics/physics.cpp -  originally src/assembler/instance/linear/instance.hpp
READEX_REGION_DEFINE(REG_Assembler_Assemble_B1);
READEX_REGION_DEFINE(REG_Assembler_Assemble_B0);

//	src/solver/generic/FETISolver.cpp - all regions renamed from LinSolver to FETISolver
READEX_REGION_DEFINE(REG_FETISolver_Init);
READEX_REGION_DEFINE(REG_FETISolver_Create_Dirichlet_Prec);	//correct position ?
READEX_REGION_DEFINE(REG_FETISolver_Solve);
//READEX_REGION_DEFINE(REG_FETISolver_SetB1);	// missing
//READEX_REGION_DEFINE(REG_FETISolver_Preprocessing);	// empty

//	src/solver/specific/cluster.cpp
READEX_REGION_DEFINE(REG_Cluster_SetClusterPC);
READEX_REGION_DEFINE(REG_Cluster_HFETIpreprocessing);
READEX_REGION_DEFINE(REG_Cluster_CreateF0_AssembleF0);
READEX_REGION_DEFINE(REG_Cluster_CreateF0_FactF0);
READEX_REGION_DEFINE(REG_Cluster_CreateSa);
READEX_REGION_DEFINE(REG_Cluster_CreateSa_SolveF0vG0);
READEX_REGION_DEFINE(REG_Cluster_CreateSa_SaReg);
READEX_REGION_DEFINE(REG_Cluster_CreateSa_SaFactorization);
READEX_REGION_DEFINE(REG_Cluster_CreateG1_perCluster);

//	src/solver/specific/itersolver.cpp
READEX_REGION_DEFINE(REG_Solve_RegCG); // renamed, originally REG_Solve_RegCG_Singular
READEX_REGION_DEFINE(REG_IterSolver_MakeSolution);
READEX_REGION_DEFINE(REG_RegCG_AllIterations);
READEX_REGION_DEFINE(REG_RegCG_OneIteration);
READEX_REGION_DEFINE(REG_Create_GGT_Inv);
READEX_REGION_DEFINE(REG_Create_GGT_Inv_Exchange_local_GGt_MPI);
READEX_REGION_DEFINE(REG_Create_GGT_Inv_G1_MatMat_MatAdd);
READEX_REGION_DEFINE(REG_Create_GGT_Inv_CollectGGtPiecesToMaster);
READEX_REGION_DEFINE(REG_Create_GGT_Inv_GGt_Factorization);
READEX_REGION_DEFINE(REG_Create_GGT_Inv_GGt_Solve);
READEX_REGION_DEFINE(REG_Projector_Inv);
READEX_REGION_DEFINE(REG_Projector_Inv_AllGather);
READEX_REGION_DEFINE(REG_Projector_Inv_DenseMatVec);

//	src/solver/specific/cpu/clustercpu.cpp
READEX_REGION_DEFINE(REG_Cluster_CreateLSC);
//READEX_REGION_DEFINE(REG_Cluster_Kfactorization);	//deleted

//	src/solver/specific/cpu/itersolvercpu.cpp
READEX_REGION_DEFINE(REG_Apply_A);
READEX_REGION_DEFINE(REG_Apply_Prec);
	
#endif /* READEX_REGIONS_H_ */
