/*
 * linearelasticity.cpp
 *
 *  Created on: Sep 27, 2015
 *      Author: lriha
 */

#include "linearelasticity.h"

void Linear_elasticity::init() {

	DOFS_PER_NODE = 3; //TODO - nacist z config souboru
	timeEvalMain.SetName("Linear Elasticity Solver Overal Timing");

    MPI_rank = _instance.rank();
    MPI_size = _instance.size();

	timeEvalMain.totalTime.AddStartWithBarrier();

	std::cout.precision(15);

	partsCount 	  = _instance.mesh().parts();

	TimeEvent timeKasm(string("Create K matrices ans RHS"));
	timeKasm.AddStart();

    // *************************************************************************************

	K_mat.reserve(partsCount);
	for (eslocal d = 0; d < partsCount; d++) {
		K_mat.push_back( SparseCSRMatrix<eslocal>(0,0) );
		f_vec.push_back( std::vector <double> ());
	}

	cilk_for (eslocal d = 0; d < partsCount; d++) {
		//eslocal dimension = _instance.mesh().getPartNodesCount(d) * mesh::Point::size();

		_instance.mesh().elasticity(K_mat[d], f_vec[d], d);

        if (_instance.rank() == 0) std::cout << d << " " ; //<< std::endl;
	}
	if (_instance.rank() == 0) std::cout << std::endl;
	timeKasm.AddEndWithBarrier();
	timeEvalMain.AddEvent(timeKasm);

    // *************************************************************************************

	TimeEvent timeFnodes(string("Create Fix nodes"));
	timeFnodes.AddStart();

	size_t fixPointsCount = _instance.mesh().getFixPointsCount();
	const std::vector<eslocal> fixPoints = _instance.mesh().getFixPoints();
	fix_nodes.resize(partsCount);

	cilk_for (eslocal d = 0; d < partsCount; d++) {
			for (eslocal fixPoint = 0; fixPoint < fixPointsCount; fixPoint++) {
				fix_nodes[d].push_back(fixPoints[d * fixPointsCount + fixPoint]);
			}
			std::sort ( fix_nodes[d].begin(), fix_nodes[d].end() );
		}

	 timeFnodes.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeFnodes);

	 // *************************************************************************************

	 TimeEvent timeB1loc(string("Create B1 local"));
	 timeB1loc.AddStart();

	 _instance.localBoundaries().create_B1_l<eslocal>(
		B1_mat,
		B0_mat,
		l2g_vec,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		lambda_map_sub_B0,
		B1_duplicity,
		partsCount,
		DOFS_PER_NODE,
		_instance.globalBoundaries()
	);

	 timeB1loc.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1loc);

	 // ************************************************************************************

	 TimeEvent timeB1glob(string("Create B1 global"));
	 timeB1glob.AddStart();

	 _instance.globalBoundaries().create_B1_g<eslocal>(
		B1_mat,
		K_mat,
		lambda_map_sub_clst,
		lambda_map_sub_B1,
		B1_duplicity,
		_instance.rank(),
		_instance.size(),
		partsCount,
		DOFS_PER_NODE,
		neigh_clusters,
        _instance.localBoundaries()
	);

	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);


	 // ************************************************************************************

	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);

	 TimeEvent timeBforces(string("Create boundary forces ??"));
	 timeBforces.AddStart();

	 //TODO: DOFS_PER_NODE
	 const std::map<eslocal, double> &forces_x = _instance.mesh().coordinates().property(mesh::CP::FORCES_X).values();
	 const std::map<eslocal, double> &forces_y = _instance.mesh().coordinates().property(mesh::CP::FORCES_Y).values();
	 const std::map<eslocal, double> &forces_z = _instance.mesh().coordinates().property(mesh::CP::FORCES_Z).values();

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



	 TimeEvent timeMconv(string("Convert M, B1 and B0 to Solver Format"));
	 timeMconv.AddStart();

	K_mat_ls.resize(partsCount);
	B1_mat_ls.resize(partsCount);
	B0_mat_ls.resize(partsCount);

	cilk_for (eslocal d = 0; d < partsCount; d++) {
 		K_mat_ls[d]  = K_mat[d];
 		B1_mat_ls[d] = B1_mat[d];
 		B0_mat_ls[d] = B0_mat[d];
 	}

	K_mat.clear();
	B1_mat.clear();
	B0_mat.clear();

	 timeMconv.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeMconv);

	 lin_solver.setup( _instance.rank(), _instance.size(), true );

	 lin_solver.init(

		_instance.mesh(),

		K_mat_ls,

		B1_mat_ls,
		B0_mat_ls,

		lambda_map_sub_B1,
		lambda_map_sub_B0,
		lambda_map_sub_clst,
		B1_duplicity,

		f_vec,
		fix_nodes,
		l2g_vec,

		neigh_clusters

	);

}


void Linear_elasticity::pre_solve_update() {

}

void Linear_elasticity::post_solve_update() {

	 	 TimeEvent timeSaveVTK(string("Solver - Save VTK"));
	 	 timeSaveVTK.AddStart();
		std::stringstream ss;
		ss << "mesh_" << MPI_rank << ".vtk";
		_instance.mesh().saveVTK(ss.str().c_str(), prim_solution, l2g_vec, _instance.localBoundaries(), _instance.globalBoundaries(), 0.95, 0.9);
		 timeSaveVTK.AddEndWithBarrier();
	  	 timeEvalMain.AddEvent(timeSaveVTK);

}


void Linear_elasticity::solve(){
	lin_solver.Solve(f_vec, prim_solution);
}

void Linear_elasticity::finalize() {
	lin_solver.finilize();
}
