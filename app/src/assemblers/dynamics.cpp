
#include "dynamics.h"


void Dynamics::init()
{
	DOFS_PER_NODE = 3; //TODO - nacist z config souboru


	dynamic_beta     = 0.25;
	dynamic_gama     = 0.5;
	dynamic_timestep = 0.000001;
	time_const = 1.0 / ( dynamic_beta * dynamic_timestep * dynamic_timestep);

	// END Setup Variables *******************************************************************

	timeEvalMain.SetName("Dynamic Elasticity Solver Overal Timing");
	timeEvalMain.totalTime.AddStartWithBarrier();

	std::cout.precision(15);

	partsCount 	  = _instance.mesh().parts();
    MPI_rank = _instance.rank();
    MPI_size = _instance.size();

	TimeEvent timeKasm(string("Create K matrices ans RHS"));
	timeKasm.AddStart();

	K_mat.reserve(partsCount);
	M_mat.reserve(partsCount);
	for (eslocal d = 0; d < partsCount; d++) {
		K_mat.push_back( SparseCSRMatrix<eslocal>(0,0) );
		M_mat.push_back( SparseCSRMatrix<eslocal>(0,0) );
		f_vec.push_back( std::vector <double> ());
	}

	cilk_for (eslocal d = 0; d < partsCount; d++) {
		_instance.mesh().elasticity(K_mat[d], M_mat[d], f_vec[d], d);
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
		_instance.globalBoundaries(),
		_instance.mesh().coordinates()
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
        _instance.localBoundaries(),
		_instance.mesh().coordinates()
	);


	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);

	 // ************************************************************************************

	 timeB1glob.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeB1glob);

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
	M_mat_ls.resize(partsCount);
	B1_mat_ls.resize(partsCount);
	B0_mat_ls.resize(partsCount);
	domain_prim_size.resize(partsCount);

	cilk_for (eslocal d = 0; d < partsCount; d++) {
 		K_mat_ls[d]  = K_mat[d];
 		M_mat_ls[d]  = M_mat[d];

 		K_mat_ls[d].MatAddInPlace( M_mat_ls[d] ,'N', time_const);
 		domain_prim_size[d] = K_mat[d].rows();
 		B1_mat_ls[d] = B1_mat[d];
 		B0_mat_ls[d] = B0_mat[d];
 	}

	M_mat.clear();
	K_mat.clear();
	B1_mat.clear();
	B0_mat.clear();

	 timeMconv.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeMconv);


	 lin_solver.setup( _instance.rank(), _instance.size(), false );

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



	 vec_u.resize( partsCount );
	 vec_v.resize( partsCount );
	 vec_w.resize( partsCount );

	 vec_u_n.resize( partsCount );
	 vec_v_n.resize( partsCount );
	 vec_w_n.resize( partsCount );

	 vec_b.resize( partsCount );
	 vec_t_tmp.resize( partsCount );

 	for (int d = 0; d < partsCount; d++) {
 		vec_u[d]    .resize(domain_prim_size[d], 0);
 		vec_v[d]    .resize(domain_prim_size[d], 0);
 		vec_w[d]    .resize(domain_prim_size[d], 0);

 		vec_u_n[d]  .resize(domain_prim_size[d], 0);
 		vec_v_n[d]  .resize(domain_prim_size[d], 0);
 		vec_w_n[d]  .resize(domain_prim_size[d], 0);

 		vec_b[d]    .resize(domain_prim_size[d], 0);
 		vec_t_tmp[d].resize(domain_prim_size[d], 0);
 	}

 	// *** Set up the initial acceleration ***********************
 	for (int d = 0; d < partsCount; d++) {
 		for (int i = 2; i < vec_w[d].size(); i=i+3) {
 			vec_w[d][i] = 1.0;
 		}
 	}
 	// *** END - Set up the initial accel. ***********************

	double const_beta   = dynamic_beta;
	double const_deltat = dynamic_timestep;
	double const_gama   = dynamic_gama;


	const_a.resize(8);

	const_a[0] = 1.0 / (const_beta * const_deltat * const_deltat);
	const_a[1] = const_gama / (const_beta * const_deltat);
	const_a[2] = 1.0 / (const_beta * const_deltat);
	const_a[3] = (1.0 / (2 * const_beta)) - 1.0;
	const_a[4] = (const_gama / const_beta) - 1.0;
	const_a[5] = const_deltat * ((const_gama / (2.0 * const_beta)) - 1.0);
	const_a[6] = const_deltat * (1.0 - const_gama);
	const_a[7] = const_deltat * const_gama;

	timeStep = 0;

}

void Dynamics::pre_solve_update()
{
	// TODO
	//_instance

	if (MPI_rank == 0) {
		cout << endl << "Time iter "  << timeStep << "\t" << endl;
	}

	// *** calculate the right hand side in primal ********************************************
	cilk_for (int d = 0; d < partsCount; d++) {
		for(int i = 0; i < vec_u[d].size(); i++) {
			vec_t_tmp[d][i] = const_a[0] * vec_u[d][i] + const_a[2] * vec_v[d][i] + const_a[3] * vec_w[d][i];
		}
		M_mat_ls[d].MatVec(vec_t_tmp[d], vec_b[d],'N');
	}


}

void Dynamics::post_solve_update()
{

	cilk_for (int d = 0; d < partsCount; d++) {
		for(int i = 0; i < vec_u[d].size(); i++) {
			vec_w_n[d][i] = (const_a[0] * (vec_u_n[d][i] - vec_u[d][i])) - (const_a[2] * vec_v[d][i]) - (const_a[3] * vec_w  [d][i]);
			vec_v_n[d][i] = vec_v[d][i]                  + (const_a[6] * vec_w[d][i])                 + (const_a[7] * vec_w_n[d][i]);

			vec_u[d][i] = vec_u_n[d][i];
			vec_v[d][i] = vec_v_n[d][i];
			vec_w[d][i] = vec_w_n[d][i];
		}
	}

#ifdef CATALYST
	unsigned int timeStep = tt;
	double time = timeStep * dynamic_timestep;
	Adaptor::CoProcess(input.mesh,l2g_vec, vec_u_n,  time, timeStep, timeStep == numberOfTimeSteps - 1);
#endif

	std::stringstream ss;
	ss << "mesh_" << _instance.rank() << "_" << timeStep << ".vtk";
	// TODO: return save VTK
	//_instance.mesh().saveVTK(ss.str().c_str(), vec_u_n, l2g_vec, _instance.localBoundaries(), _instance.globalBoundaries(), 0.95, 0.9);

	timeStep++;
}

void Dynamics::solve()
{
	// *** Run the CG solver **************************************************************
	lin_solver.Solve(vec_b, vec_u_n);
	// *** END - Run the CG solver ********************************************************

}

void Dynamics::finalize() {
	lin_solver.finilize();
}


