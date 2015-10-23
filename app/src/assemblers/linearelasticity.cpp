/*
 * linearelasticity.cpp
 *
 *  Created on: Sep 27, 2015
 *      Author: lriha
 */

#include "linearelasticity.h"

void Linear_elasticity::init() {

	BEM = true;
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

	if (BEM) {
		cilk_for (eslocal d = 0; d < partsCount; d++) {
			DenseMatrix K_dense, K_tmp;
			_instance.surf_mesh().elasticity(K_dense, d);
			//TODO: not optimal
			K_tmp = K_dense;

			int n = K_dense.rows();

	        for (int i = 0; i < n/3; i++) {
	            for (int j = 0; j < n; j++) {
	                K_tmp( 3*i+0,j) = K_dense(0*(n/3) + i ,j);
	                K_tmp( 3*i+1,j) = K_dense(1*(n/3) + i ,j);
	                K_tmp( 3*i+2,j) = K_dense(2*(n/3) + i ,j);
	            }
	        }

	        for (int i = 0; i < n/3; i++) {
	            for (int j = 0; j < n; j++) {
	                K_dense( j, 3*i+0) = K_tmp(j, 0*(n/3) + i );
	                K_dense( j, 3*i+1) = K_tmp(j, 1*(n/3) + i );
	                K_dense( j, 3*i+2) = K_tmp(j, 2*(n/3) + i );
	            }
	        }

	        K_mat[d] = K_dense;
	        f_vec[d].resize(K_dense.rows() , 0.0);

	        //TODO: Musi byt pro kostku
	        _instance.surf_mesh().integrateUpperFaces(f_vec[d],d);

	        if (_instance.rank() == 0) std::cout << d << " " ; //<< std::endl;

		}

	} else {

		cilk_for (eslocal d = 0; d < partsCount; d++) {
			//eslocal dimension = _instance.mesh().getPartNodesCount(d) * mesh::Point::size();

			_instance.mesh().elasticity(K_mat[d], f_vec[d], d);

			if (_instance.rank() == 0) std::cout << d << " " ; //<< std::endl;
		}

	}

	if (_instance.rank() == 0) std::cout << std::endl;
	timeKasm.AddEndWithBarrier();
	timeEvalMain.AddEvent(timeKasm);

    // *************************************************************************************

	TimeEvent timeFnodes(string("Create Fix nodes"));
	timeFnodes.AddStart();

	size_t fixPointsCount;
	std::vector<eslocal> fixPoints;
	if (BEM) {

		fixPointsCount = _instance.surf_mesh().getFixPointsCount();
		fixPoints      = _instance.surf_mesh().getFixPoints();
	} else {
		fixPoints      = _instance.mesh().getFixPoints();
		fixPointsCount = _instance.mesh().getFixPointsCount();
	}

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

	 //TODO: implement filling of vec_c
	 vec_c.resize(partsCount);

	 if (BEM) {
		 _instance.surf_mesh().subdomainBoundaries().create_B1_l<eslocal>(
			B1_mat,
			B0_mat,
			l2g_vec,
			lambda_map_sub_clst,
			lambda_map_sub_B1,
			lambda_map_sub_B0,
			B1_duplicity,
			vec_c,
			partsCount,
			DOFS_PER_NODE,
			_instance.globalBoundaries(),
			_instance.surf_mesh().coordinates()

		);
	 } else {
		 _instance.localBoundaries().create_B1_l<eslocal>(
			B1_mat,
			B0_mat,
			l2g_vec,
			lambda_map_sub_clst,
			lambda_map_sub_B1,
			lambda_map_sub_B0,
			B1_duplicity,
			vec_c,
			partsCount,
			DOFS_PER_NODE,
			_instance.globalBoundaries(),
			_instance.mesh().coordinates()
		 );
	 }
		 timeB1loc.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeB1loc);

		 // ************************************************************************************

		 TimeEvent timeB1glob(string("Create B1 global"));
		 timeB1glob.AddStart();

	if (BEM) {
		 _instance.globalBoundaries().create_B1_g<eslocal>(
			B1_mat,
			K_mat,
			lambda_map_sub_clst,
			lambda_map_sub_B1,
			B1_duplicity,
			vec_c,
			_instance.rank(),
			_instance.size(),
			partsCount,
			DOFS_PER_NODE,
			neigh_clusters,
			_instance.localBoundaries(),
			_instance.surf_mesh().coordinates()
		);
	} else {
		 _instance.globalBoundaries().create_B1_g<eslocal>(
			B1_mat,
			K_mat,
			lambda_map_sub_clst,
			lambda_map_sub_B1,
			B1_duplicity,
			vec_c,
			_instance.rank(),
			_instance.size(),
			partsCount,
			DOFS_PER_NODE,
			neigh_clusters,
			_instance.localBoundaries(),
			_instance.mesh().coordinates()
		);

	}


	 for (int i = 0; i < vec_c.size(); i++)
		 vec_c[i].resize(B1_duplicity[i].size(), 0.0);

	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);


	 // ************************************************************************************

	 TimeEvent timeBforces(string("Create boundary forces ??"));
	 timeBforces.AddStart();

	 //TODO: DOFS_PER_NODE
	 const std::map<eslocal, double> &forces_x = _instance.mesh().coordinates().property(mesh::FORCES_X).values();
	 const std::map<eslocal, double> &forces_y = _instance.mesh().coordinates().property(mesh::FORCES_Y).values();
	 const std::map<eslocal, double> &forces_z = _instance.mesh().coordinates().property(mesh::FORCES_Z).values();

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

	//std::ofstream file ("K_mat");
	//file << K_mat[0];
	//file.close();

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

	 TimeEvent timeLSconv(string("Linear Solver - preprocessing"));
	 timeLSconv.AddStart();

	 lin_solver.DOFS_PER_NODE = DOFS_PER_NODE;
	 lin_solver.setup( _instance.rank(), _instance.size(), true );

	 if (BEM) {

		 lin_solver.init(

			_instance.surf_mesh(),

			K_mat_ls,

			B1_mat_ls,
			B0_mat_ls,

			lambda_map_sub_B1,
			lambda_map_sub_B0,
			lambda_map_sub_clst,
			B1_duplicity,

			f_vec,
			vec_c,

			fix_nodes,
			l2g_vec,

			neigh_clusters

		);
	 } else {
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
			vec_c,

			fix_nodes,
			l2g_vec,

			neigh_clusters

		);

	 }

	 timeLSconv.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeLSconv);


}


void Linear_elasticity::pre_solve_update() {

}

void Linear_elasticity::post_solve_update() {

	 	 TimeEvent timeSaveVTK(string("Solver - Save VTK"));
	 	 timeSaveVTK.AddStart();

	 	 if (BEM) {
	 		_instance.surf_mesh().store(mesh::VTK, "mesh", prim_solution, 0.95, 0.9);
	 	 } else {
	 	 	_instance.mesh().     store(mesh::VTK, "mesh", prim_solution, 0.95, 0.9);
	 	 }

		 timeSaveVTK.AddEndWithBarrier();
	  	 timeEvalMain.AddEvent(timeSaveVTK);

}


void Linear_elasticity::solve(){

	 TimeEvent timeLSrun(string("Linear Solver - runtime"));
	 timeLSrun.AddStart();
	lin_solver.Solve(f_vec, prim_solution);
	 timeLSrun.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeLSrun);
}

void Linear_elasticity::finalize() {
	lin_solver.finilize();

	 timeEvalMain.totalTime.AddEndWithBarrier();
	 timeEvalMain.PrintStatsMPI();
}
