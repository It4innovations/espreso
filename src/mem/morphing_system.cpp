
#include "morphing_system.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

MorphingMatrix::MorphingMatrix(){}

MorphingMatrix::~MorphingMatrix(){
	
	for(auto &el: this->A){
		if(el){
			delete el;
			el = nullptr;
		}
	}

	if(this->P){
		delete this->P;
		this->P = nullptr;
	}
	
}

MorphingMatrix::MorphingMatrix(
	int dim,
	bool tr_translate,
	bool tr_scale,
	bool tr_rotate,
	int regularization_poly_degree,
	const std::vector<Point> &sPoints,
	const std::vector<double> &sDisplacement,
	const RBFTargetConfiguration &configuration,
	double aca_eps,
	double aca_eta,
	esint base_cluster_size,
		
	std::vector<double> &morphing_rhs
){

	this->dimension = dim;
	this->n_regularization = 0;
	this->regularization_polynomial_degree = regularization_poly_degree;
	
	
	this->use_transformation_translation = tr_translate;
	this->use_transformation_scaling = tr_scale;
	this->use_transformation_rotation = tr_rotate;
	
	auto nranks = info::mpi::size;
	this->A.reserve(nranks);
	this->A_inner_block_idx.reserve(nranks);
	this->A_outer_block_idx.reserve(nranks);

	std::vector<esint> nPoints(nranks);
	std::fill(nPoints.begin(), nPoints.end(), 0);
	nPoints[info::mpi::rank] = sPoints.size();
	Communication::allReduce(nPoints, Communication::OP::SUM);
	
	/*
		block_distribution[c + r*nranks] denotes the index of an MPI process assembling that block
	*/
	std::vector<esint> block_distribution;
	
	/*
		if active_block_indices[i] is true, then this MPI process requires data from the i-th MPI process
	*/
	std::vector<bool> active_block_indices;

	/*
		block_processes[i] denotes a list of MPI processes requiring data from the i-th MPI process
	*/
	std::vector<std::vector<esint>> block_processes;
	
	/*
		process_blocks[i] denotes a list of blocks required by the i-th MPI process
	*/
	std::vector<std::vector<esint>> process_blocks;
	
	CyclicGraphDecomposition::define_workload(
		block_distribution, 
		active_block_indices,
		block_processes,
		process_blocks,
		nranks,
		info::mpi::rank
	);
	
	//transform all points
	this->define_geometry_transformations(
		sPoints
	);

	std::vector<std::vector<Point>> active_point_set(nranks);
	active_point_set[info::mpi::rank] = sPoints;
	this->transform_points(
		active_point_set[info::mpi::rank]
	);
	auto displacementLocal = sDisplacement;
	this->transform_displacements(displacementLocal);

	//synchronizing relevant data among processes
	//TODO move into communication...
	esint nReqs = block_processes[info::mpi::rank].size() + process_blocks[info::mpi::rank].size();
	if(active_block_indices[info::mpi::rank]){
		nReqs -= 2;
	}
	
	std::vector<MPI_Request> requests(nReqs);
	esint req_idx = 0;
	//sending own data
	for(auto el: block_processes[info::mpi::rank]){
		if(el != info::mpi::rank){
			MPI_Isend(
				&active_point_set[info::mpi::rank][0], 
				nPoints[info::mpi::rank] * sizeof(Point), 
				MPI_BYTE, 
				el, 
				0,
				MPI_COMM_WORLD, 
				&requests[req_idx]
			);
			req_idx++;
		}
	}
	
	//recieving foreign data
	// printf("point_data = [\n");
	for(auto el: process_blocks[info::mpi::rank]){
		if(el != info::mpi::rank){
			active_point_set[el].resize(nPoints[el]);
			
			MPI_Irecv(
				&active_point_set[el][0], 
				nPoints[el] * sizeof(Point), 
				MPI_BYTE, 
				el,
				0,
				MPI_COMM_WORLD, 
				&requests[req_idx]
			);
			req_idx++;
		}

		// for(auto &p: active_point_set[el]){
		// 	printf("  %10.8f, %10.8f, %10.8f;\n", p.x, p.y, p.z);
		// }
	}
	MPI_Waitall(req_idx, requests.data(), MPI_STATUSES_IGNORE);
	// printf("];\n");

	
	//construction of relevant trees and global DOFs
	std::vector<ClusterTree*> trees(nranks);
	std::fill(trees.begin(), trees.end(), nullptr);
	for(int i = 0; i < nranks; ++i){
		if(active_block_indices[i]){
			trees[ i ] = new ClusterTree(&active_point_set[ i ]);
			trees[ i ]->createClusterTree(base_cluster_size);
		}
	}
	
	
	this->global_DOFs_start.resize(nranks);
	this->global_DOFs_end.resize(nranks);
	esint npoints_total = 0;
	
	for( int i = 0; i < nranks; ++i ){
		
		global_DOFs_start[i] = npoints_total;
		global_DOFs_end[i] = global_DOFs_start[i] + nPoints[ i ];
		
		npoints_total += nPoints[ i ];
	}

	for(int r = 0; r < nranks; ++r){
		for(int c = 0; c < nranks; ++c){
			if(block_distribution[c + r * nranks] == info::mpi::rank){
				this->A.push_back(
					new matrix_ACA(
						*trees[ r ], 
						*trees[ c ], 
						aca_eps, 
						aca_eta,
						configuration
					)
				);
				this->A_inner_block_idx.push_back(c);
				this->A_outer_block_idx.push_back(r);
			}
		}
	}
	
	this->recalculateB(
		active_point_set[info::mpi::rank]
	);
	
	this->ncols = npoints_total + this->n_regularization;
	this->nrows = this->ncols;


	// this->apply_data_mem.resize(this->nrows);

	// this->P = new MorphingMatrixPreconditioner(
	// 	npoints_total, 
	// 	this->n_regularization,
	// 	this->global_DOFs_start[info::mpi::rank],
	// 	active_point_set[info::mpi::rank], 
	// 	configuration
	// );
	
	morphing_rhs.resize(this->nrows * this->dimension);
	std::fill(morphing_rhs.begin(), morphing_rhs.end(), 0.0f);
	size_t shift_rhs = this->global_DOFs_start[info::mpi::rank];
	for (int d = 0; d < this->dimension; d++) {
		for (size_t r = 0; r < active_point_set[info::mpi::rank].size(); r++) {
			morphing_rhs[this->nrows * d + r + shift_rhs] = displacementLocal[r * this->dimension + d];
		}
	}
	Communication::allReduce(morphing_rhs, Communication::OP::SUM);
	
	for(auto &el:trees){
		if(el){
			delete el;
			el = nullptr;
		}
	}
	
	// double compr = this->getCompressionRatio();
	// if(info::mpi::rank == 0){
		// printf("ACA morphing matrix constructed, compression ratio %f [%%]\n", 100.0f * compr);
	// }

	// this->print(displacementLocal);
	// this->P->printData();
}

void MorphingMatrix::transform_points(
	std::vector<Point> &sPoints
) const {
	for(auto &p: sPoints){
		this->transform_point(p);
	}
}

void MorphingMatrix::transform_points_inverse(
	std::vector<Point> &sPoints
) const {
	for(auto &p: sPoints){
		this->transform_point_inverse(p);
	}
}

void MorphingMatrix::transform_point(
	Point &p
) const {
	
	if(this->use_transformation_translation){
		p.x += this->system_shift[0];
		p.y += this->system_shift[1];
		p.z += this->system_shift[2];
	}
	
	if(this->use_transformation_scaling){
		p.x *= this->scaling[0];
		p.y *= this->scaling[1];
		p.z *= this->scaling[2];
	}
	
	
	if(this->use_transformation_rotation){
		double x = p.x;
		double y = p.y;
		double z = p.z;
		
		p.x = this->rotation_matrix[0][0] * x + this->rotation_matrix[0][1] * y + this->rotation_matrix[0][2] * z;
		p.y = this->rotation_matrix[1][0] * x + this->rotation_matrix[1][1] * y + this->rotation_matrix[1][2] * z;
		p.z = this->rotation_matrix[2][0] * x + this->rotation_matrix[2][1] * y + this->rotation_matrix[2][2] * z;
		
		// printf("(%10.6f, %10.6f, %10.6f) - > (%10.6f, %10.6f, %10.6f)\n", x, y, z, p.x, p.y, p.z);
	}
	
}

void MorphingMatrix::transform_point_inverse(
	Point &p
) const {

	if(this->use_transformation_rotation){
		double x = p.x;
		double y = p.y;
		double z = p.z;
		
		p.x = this->rotation_matrix[0][0] * x + this->rotation_matrix[1][0] * y + this->rotation_matrix[2][0] * z;
		p.y = this->rotation_matrix[0][1] * x + this->rotation_matrix[1][1] * y + this->rotation_matrix[2][1] * z;
		p.z = this->rotation_matrix[0][2] * x + this->rotation_matrix[1][2] * y + this->rotation_matrix[2][2] * z;
	}

	if(this->use_transformation_scaling){
		p.x /= this->scaling[0];
		p.y /= this->scaling[1];
		p.z /= this->scaling[2];
	}

	if(this->use_transformation_translation){
		p.x -= this->system_shift[0];
		p.y -= this->system_shift[1];
		p.z -= this->system_shift[2];
	}
	
}

void MorphingMatrix::transform_displacements(
	std::vector<double> &displacements
) const {
	
	esint npoints = displacements.size() / this->dimension;
	for (size_t r = 0; r < npoints; r++) {
		this->transform_displacement(
			&displacements[r * this->dimension]
		);
	}
}

void MorphingMatrix::transform_displacement(
	double *displacement
) const {
	//shift has no effect on displacement
	
	//scaling
	if(this->use_transformation_scaling){
		for (int d = 0; d < this->dimension; d++) {
			displacement[d] *= this->scaling[d];
		}	
	}
	
	//rotation
	if(this->use_transformation_rotation){
		double x = displacement[0];
		double y = displacement[1];
		double z = 0.0f;
		
		if(this->dimension == 3){
			z = displacement[2];
		}
		
		displacement[0] = this->rotation_matrix[0][0] * x + this->rotation_matrix[0][1] * y + this->rotation_matrix[0][2] * z;
		displacement[1] = this->rotation_matrix[1][0] * x + this->rotation_matrix[1][1] * y + this->rotation_matrix[1][2] * z;
		if(this->dimension == 3){
			displacement[2] = this->rotation_matrix[2][0] * x + this->rotation_matrix[2][1] * y + this->rotation_matrix[2][2] * z;
		}
	}

}


void MorphingMatrix::define_geometry_transformations(
	const std::vector<Point> &sPoints
){
	std::vector<int> npoints_total = {(int)sPoints.size()};
	
	double dmin = std::numeric_limits<double>::min();
	double dmax = std::numeric_limits<double>::max();
	std::vector<double> coordinate_max = {dmin, dmin, dmin};
	std::vector<double> coordinate_min = {dmax, dmax, dmax};
	std::vector<double> transform_centroid = {0.0f, 0.0f, 0.0f};
	
	if(sPoints.size() > 0){
		
		for(auto i = 0; i < sPoints.size(); ++i){
			
			auto &p = sPoints[i];
			
			if(p.x > coordinate_max[0]){
				coordinate_max[0] = p.x;
			}
			if(p.y > coordinate_max[1]){
				coordinate_max[1] = p.y;
			}
			if(p.z > coordinate_max[2]){
				coordinate_max[2] = p.z;
			}
			if(p.x < coordinate_min[0]){
				coordinate_min[0] = p.x;
			}
			if(p.y < coordinate_min[1]){
				coordinate_min[1] = p.y;
			}
			if(p.z < coordinate_min[2]){
				coordinate_min[2] = p.z;
			}
			
			transform_centroid[0] += p.x;
			transform_centroid[1] += p.y;
			transform_centroid[2] += p.z;
		}
	}
	
	Communication::allReduce(coordinate_max, Communication::OP::MAX);
	Communication::allReduce(coordinate_min, Communication::OP::MIN);
	this->scaling[0] = 1.0f/(coordinate_max[0] - coordinate_min[0]);
	this->scaling[1] = 1.0f/(coordinate_max[1] - coordinate_min[1]);
	this->scaling[2] = 1.0f/(coordinate_max[2] - coordinate_min[2]);
	
	

	Communication::allReduce(transform_centroid, Communication::OP::SUM);
	Communication::allReduce(npoints_total, Communication::OP::SUM);
	this->system_shift[0] = -transform_centroid[0] / npoints_total[0];
	this->system_shift[1] = -transform_centroid[1] / npoints_total[0];
	this->system_shift[2] = -transform_centroid[2] / npoints_total[0];
	
	
	//rotation transformation
	std::vector<double> rotation_axis = {1.0f, 0.5f, 0.25f};
	double rotation_angle = 2.0f;//in radians
	
	double axis_norm = std::sqrt(rotation_axis[0] * rotation_axis[0] + rotation_axis[1] * rotation_axis[1] + rotation_axis[2] * rotation_axis[2]);
	
	double ux = rotation_axis[0] / axis_norm;
	double uy = rotation_axis[1] / axis_norm;
	double uz = rotation_axis[2] / axis_norm;
	
	double uxx = ux * ux;
	double uxy = ux * uy;
	double uxz = ux * uz;
	double uyy = uy * uy;
	double uyz = uy * uz;
	double uzz = uz * uz;
	
	double cosa = std::cos(rotation_angle);
	double sina = std::sin(rotation_angle);

	this->rotation_matrix[0][0] = uxx * (1.0f - cosa) + cosa;
	this->rotation_matrix[0][1] = uxy * (1.0f - cosa) - uz * sina;
	this->rotation_matrix[0][2] = uxz * (1.0f - cosa) + uy * sina;
	
	this->rotation_matrix[1][0] = uxy * (1.0f - cosa) + uz * sina;
	this->rotation_matrix[1][1] = uyy * (1.0f - cosa) + cosa;
	this->rotation_matrix[1][2] = uyz * (1.0f - cosa) - ux * sina;
	
	this->rotation_matrix[2][0] = uxz * (1.0f - cosa) - uy * sina;
	this->rotation_matrix[2][1] = uyz * (1.0f - cosa) + ux * sina;
	this->rotation_matrix[2][2] = uzz * (1.0f - cosa) + cosa;
}

double MorphingMatrix::calculateConstantShift(
	const Point &Pr
) const {
	
	double P_ = 0.0f;
	double P_d = 0.0f;
	
	if(this->n_regularization > 1){
		P_ += Pr.x;
		P_d += 1.0f;
	}
	
	if(this->n_regularization > 2){
		P_ += Pr.y;
		P_d += 1.0f;
	}
	
	if(this->n_regularization > 3){
		P_ += Pr.z;
		P_d += 1.0f;
	}
	
	P_ /= P_d;
	
	double out = 1.0f;
	if(this->n_regularization > 1){
		out += (Pr.x + P_)*(Pr.x + P_);
	}
	
	if(this->n_regularization > 2){
		out += (Pr.y + P_)*(Pr.y + P_);
	}
	
	if(this->n_regularization > 3){
		out += (Pr.z + P_)*(Pr.z + P_);
	}
	
	return 1.0f;
	return std::sqrt(out);
}

void MorphingMatrix::recalculateB(
	const std::vector<Point> &sPoints
){
	this->B.resize(0);
	
	if(this->regularization_polynomial_degree > 0){
		this->n_regularization++;//constant term
	}
	if(this->regularization_polynomial_degree > 1){
		this->n_regularization += this->dimension;//linear terms
	}
	if(this->regularization_polynomial_degree > 2){
		//quadratic terms
		if(this->dimension == 2){
			this->n_regularization += 3;
		}
		else if(this->dimension == 3){
			this->n_regularization += 6;
		}
	}
	
	
	this->B.reserve(sPoints.size() * this->n_regularization);
	
	//constant term
	if(this->regularization_polynomial_degree >= 1){
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( 1.0f );
		}
	}
	
	//linear terms
	if(this->regularization_polynomial_degree >= 2){
		//x
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( sPoints[c].x );
		}
		//y
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( sPoints[c].y );
		}
		
		if(this->dimension >= 3){
			//z
			for (size_t c = 0; c < sPoints.size(); c++) {
				this->B.push_back( sPoints[c].z );
			}
		}
	}
	
	//quadratic terms terms
	if(this->regularization_polynomial_degree >= 3){
		//xx
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( sPoints[c].x * sPoints[c].x );
		}
		//xy
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( sPoints[c].y * sPoints[c].x );
		}

		//yy
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back( sPoints[c].y * sPoints[c].y );
		}
		
		
		if(this->dimension >= 3){
			//xz
			for (size_t c = 0; c < sPoints.size(); c++) {
				this->B.push_back( sPoints[c].z * sPoints[c].x );
			}
			//yz
			for (size_t c = 0; c < sPoints.size(); c++) {
				this->B.push_back( sPoints[c].z * sPoints[c].y );
			}
			//zz
			for (size_t c = 0; c < sPoints.size(); c++) {
				this->B.push_back( sPoints[c].z * sPoints[c].z );
			}
		}
	}
}

void MorphingMatrix::calculateMorphingError(
	const std::vector<double> &sol,
	const std::vector<double> &rhs,
	double &error_morph_x,
	double &error_morph_y,
	double &error_morph_z,
	double &orthogonality_error_x,
	double &orthogonality_error_y,
	double &orthogonality_error_z
) const {
	std::vector<double> res_x(this->nrows),  res_y(this->nrows),  res_z(this->nrows);
	
	this->apply(
		&sol[0], 
		&res_x[0], 
		1.0f, 
		0.0f, 
		false
	);
	
	this->apply(
		&sol[this->ncols], 
		&res_y[0], 
		1.0f, 
		0.0f, 
		false
	);
	
	this->apply(
		&sol[this->ncols*2], 
		&res_z[0], 
		1.0f, 
		0.0f, 
		false
	);
	
	error_morph_x = 0.0f;
	error_morph_y = 0.0f;
	error_morph_z = 0.0f;
	double error_morph_x_ref = 0.0f;
	double error_morph_y_ref = 0.0f;
	double error_morph_z_ref = 0.0f;
	
	double displ_local[3];
	
	for(size_t i = 0; i < this->nrows - this->n_regularization; ++i){
		
		displ_local[0] = rhs[i];
		displ_local[1] = rhs[i + this->nrows];
		displ_local[2] = rhs[i + this->nrows*2];
		
		// this->transform_displacement(&displ_local[0]);
		
		double ref_x = displ_local[0];
		double ref_y = displ_local[1];
		double ref_z = displ_local[2];
		
		error_morph_x += (ref_x - res_x[i])*(ref_x - res_x[i]);
		error_morph_y += (ref_y - res_y[i])*(ref_y - res_y[i]);
		error_morph_z += (ref_z - res_z[i])*(ref_z - res_z[i]);
		
		error_morph_x_ref += ref_x * ref_x;
		error_morph_y_ref += ref_y * ref_y;
		error_morph_z_ref += ref_z * ref_z;
		
	}

	error_morph_x = std::sqrt(error_morph_x / (error_morph_x_ref > 0?error_morph_x_ref:1));
	error_morph_y = std::sqrt(error_morph_y / (error_morph_y_ref > 0?error_morph_y_ref:1));
	error_morph_z = std::sqrt(error_morph_z / (error_morph_z_ref > 0?error_morph_z_ref:1));

	orthogonality_error_x = 0.0f;
	orthogonality_error_y = 0.0f;
	orthogonality_error_z = 0.0f;
	for(size_t i = this->nrows - this->n_regularization; i < this->nrows; ++i){
		orthogonality_error_x += res_x[i] * res_x[i];
		orthogonality_error_y += res_y[i] * res_y[i];
		orthogonality_error_z += res_z[i] * res_z[i];
	}
	

	orthogonality_error_x = std::sqrt(orthogonality_error_x);
	orthogonality_error_y = std::sqrt(orthogonality_error_y);
	orthogonality_error_z = std::sqrt(orthogonality_error_z);
	
	// if(info::mpi::rank == 0){
		// printf("Morphing errors:\n");
		// if(error_morph_x_ref > 0.0f){
			// printf("  [x] -> %10.6f [%%]\n", std::sqrt(error_morph_x / error_morph_x_ref) * 100.0f);
		// }
		// else{
			// printf("  [x] -> %20.18f\n", std::sqrt(error_morph_x));
		// }
		// if(error_morph_y_ref > 0.0f){
			// printf("  [y] -> %10.6f [%%]\n", std::sqrt(error_morph_y / error_morph_y_ref) * 100.0f);
		// }
		// else{
			// printf("  [y] -> %20.18f\n", std::sqrt(error_morph_y));
		// }
		// if(error_morph_z_ref > 0.0f){
			// printf("  [z] -> %10.6f [%%]\n", std::sqrt(error_morph_z / error_morph_z_ref) * 100.0f);
		// }
		// else{
			// printf("  [z] -> %20.18f\n", std::sqrt(error_morph_z));
		// }
		
		// printf("Orthogonality errors:\n");
		// printf("  [x] -> %20.18f\n", std::fabs(orthogonality_error_x));
		// printf("  [y] -> %20.18f\n", std::fabs(orthogonality_error_y));
		// printf("  [z] -> %20.18f\n", std::fabs(orthogonality_error_z));
	// }
	
}

void MorphingMatrix::transformPoint(
	Point &point_to_morph,
	const Point &point_origin,
	const std::vector<double> &coefficients,
	const std::vector<Point> &rPoints,
	const RBFTargetConfiguration &configuration
) const {
	esint points_size = rPoints.size();
	
	esint n = this->nrows;
	esint n2 = 2 * n;
	
	Point Po = point_origin;
	this->transform_point(Po);
	
	// printf("(%10.6f, %10.6f, %10.6f)", point_to_morph.x, point_to_morph.y, point_to_morph.z);

	this->transform_point(point_to_morph);
	
	// printf(" -> (%10.6f, %10.6f, %10.6f)", point_to_morph.x, point_to_morph.y, point_to_morph.z);
	
	for (size_t i = 0; i < rPoints.size(); i++) {
		auto Pr = rPoints[i];
		this->transform_point(Pr);
		
		double R = configuration.function.evaluator->evaluate((Pr - Po).length());

		point_to_morph.x += R * coefficients[i ];
		point_to_morph.y += R * coefficients[i + n];
		if (this->dimension == 3) {
			point_to_morph.z += R * coefficients[i + n2];
		}
	}
	// printf(" -> (%10.6f, %10.6f, %10.6f)", point_to_morph.x, point_to_morph.y, point_to_morph.z);
	
	// double cs = this->calculateConstantShift(Po);
	
	//polynomial correction
	
	int shift_b = 0;
	int shift;

	//constant term
	if(this->regularization_polynomial_degree >= 1){
		point_to_morph.x += coefficients[points_size];
		point_to_morph.y += coefficients[points_size + n];
		
		if(this->dimension >= 3){
			point_to_morph.z += coefficients[points_size + n2];
		}
		
		shift_b += 1;
	}

	//linear terms
	if(this->regularization_polynomial_degree >= 2){
		shift = shift_b + points_size;
		point_to_morph.x += Po.x * coefficients[shift++];
		point_to_morph.x += Po.y * coefficients[shift++];
		if(this->dimension >= 3){
			point_to_morph.x += Po.z * coefficients[shift++];
		}

		shift = shift_b + n + points_size;
		point_to_morph.y += Po.x * coefficients[shift++];
		point_to_morph.y += Po.y * coefficients[shift++];
		if(this->dimension >= 3){
			point_to_morph.y += Po.z * coefficients[shift++];
		}

		if(this->dimension >= 3){
			shift = shift_b + n2 + points_size;
			point_to_morph.z += Po.x * coefficients[shift++];
			point_to_morph.z += Po.y * coefficients[shift++];
			point_to_morph.z += Po.z * coefficients[shift++];
		}

		shift_b += this->dimension;
	}

	//quadratic terms
	if(this->regularization_polynomial_degree >= 3){
		shift = shift_b + points_size;
		point_to_morph.x += Po.x * Po.x * coefficients[shift++];
		point_to_morph.x += Po.x * Po.y * coefficients[shift++];
		point_to_morph.x += Po.y * Po.y * coefficients[shift++];
		if(this->dimension >= 3){
			point_to_morph.x += Po.x * Po.z * coefficients[shift++];
			point_to_morph.x += Po.y * Po.z * coefficients[shift++];
			point_to_morph.x += Po.z * Po.z * coefficients[shift++];
		}

		shift = n + shift_b + points_size;
		point_to_morph.y += Po.x * Po.x * coefficients[shift++];
		point_to_morph.y += Po.x * Po.y * coefficients[shift++];
		point_to_morph.y += Po.y * Po.y * coefficients[shift++];
		if(this->dimension >= 3){
			point_to_morph.y += Po.x * Po.z * coefficients[shift++];
			point_to_morph.y += Po.y * Po.z * coefficients[shift++];
			point_to_morph.y += Po.z * Po.z * coefficients[shift++];
		}

		shift = n2 + shift_b + points_size;
		point_to_morph.z += Po.x * Po.x * coefficients[shift++];
		point_to_morph.z += Po.x * Po.y * coefficients[shift++];
		point_to_morph.z += Po.y * Po.y * coefficients[shift++];
		if(this->dimension >= 3){
			point_to_morph.z += Po.x * Po.z * coefficients[shift++];
			point_to_morph.z += Po.y * Po.z * coefficients[shift++];
			point_to_morph.z += Po.z * Po.z * coefficients[shift++];
		}

		if(this->dimension == 2){
			shift_b += 3;
		}
		else if(this->dimension == 3){
			shift_b += 6;
		}
	}

	// printf(" -> (%10.6f, %10.6f, %10.6f)", point_to_morph.x, point_to_morph.y, point_to_morph.z);
	
	this->transform_point_inverse(point_to_morph);
	// printf(" -> (%10.6f, %10.6f, %10.6f)\n", point_to_morph.x, point_to_morph.y, point_to_morph.z);
}

void MorphingMatrix::printData() const {
	//TODO	
}


//performs y = this * x * alpha + y * beta
void MorphingMatrix::apply(
	const double * x_global, 
	double *y_global, 
	double alpha, 
	double beta, 
	bool transpose
) const {
	
	MATH::vecScale((transpose?this->ncols:this->nrows), beta / info::mpi::size, y_global);
	esint bi;
	esint bo;
	
	esint shift_x;
	esint shift_y;
	
	for(int b = 0; b < this->A.size(); ++b){
		
		bi = this->A_inner_block_idx[b];
		bo = this->A_outer_block_idx[b];
		
		shift_x = this->global_DOFs_start[bi];
		shift_y = this->global_DOFs_start[bo];
		
		this->A[b]->apply(x_global + shift_x, y_global + shift_y, alpha, 1.0f, transpose);
		// this->A[b]->apply(x_global + shift_x, &this->apply_data_mem[shift_y], alpha, 0.0f, transpose);
	}
	
	
	#ifdef HAVE_MKL
	// B'
	// cblas_dgemv(
		// CblasRowMajor,
		// CblasNoTrans,
		// this->dimension + 1,
		// this->ncols - this->dimension - 1,
		// alpha,
		// &this->B[0],
		// this->ncols - this->dimension - 1,
		// &x_global[0],
		// 1,
		// 1.0f,
		// &y_global[this->nrows - this->dimension - 1],
		// 1
	// );
	// B'
	// cblas_dgemv(
		// CblasRowMajor,
		// CblasTrans,
		// this->dimension + 1,
		// this->ncols - this->dimension - 1,
		// alpha,
		// &this->B[0],
		// this->dimension + 1,
		// &x_global[this->ncols - this->dimension - 1],
		// 1,
		// 1.0f,
		// &y_global[0],
		// 1
	// );
	#else
	
	esint nc = this->global_DOFs_end[info::mpi::rank] - this->global_DOFs_start[info::mpi::rank];
	esint nr = this->nrows - n_regularization;
	
	shift_x = this->global_DOFs_start[info::mpi::rank];
	shift_y = this->global_DOFs_start[info::mpi::rank];
	for(esint r = 0; r < n_regularization; ++r){
		for(esint c = 0; c < nc; ++c){
			//B
			y_global[nr + r] += alpha * this->B[c + r * nc] * x_global[c + shift_x];
			
			//B'
			y_global[c + shift_y] += alpha * this->B[c + r * nc] * x_global[nr + r];

			// //B
			// this->apply_data_mem[nr + r] += alpha * this->B[c + r * nc] * x_global[c + shift_x];
			
			// //B'
			// this->apply_data_mem[c + shift_y] += alpha * this->B[c + r * nc] * x_global[nr + r];
		}
	}
	
	#endif

	// this->P.apply(&this->apply_data_mem[0], y_global, 1.0f, 1.0f, false);
	
	Communication::allReduce(MPI_IN_PLACE, y_global, this->nrows, MPI_DOUBLE, MPI_SUM);
}

void MorphingMatrix::applyPreconditioner(
	const double * x_global, 
	double *y_global, 
	double alpha, 
	double beta, 
	bool transpose
) const {
	this->P->apply(
		x_global,
		y_global,
		alpha,
		beta,
		transpose
	);
}

esint MorphingMatrix::getNRows() const{
	return this->nrows;
}

esint MorphingMatrix::getNCols() const{
	return this->ncols;
}

double MorphingMatrix::getCompressionRatio() const{
	esint out = this->B.size();
	
	for(auto &B: this->A){
		out += B->getNEntries();
	}
	Communication::allReduce(MPI_IN_PLACE, &out, 1, MPI_INT, MPI_SUM);
	
	return (double)out / (this->getNCols() * this->getNRows());
}

void MorphingMatrix::print(const std::vector<double> &rhs) const {
	std::vector<double> basis_vector( this->ncols );
	std::fill(basis_vector.begin(), basis_vector.end(), 0.0f);
	
	std::vector<std::vector<double>> M(this->ncols);
	for(auto &el: M){
		el.resize(this->nrows);
		std::fill(el.begin(), el.end(), 0.0f);
	}
	
	for(esint i = 0; i < this->ncols; ++i){
		basis_vector[i] = 1.0f;

		this->apply(
			&basis_vector[0],
			&M[i][0],
			1.0f,
			0.0f,
			false
		);

		basis_vector[i] = 0.0f;
	}
	
	if(info::mpi::rank == 0){
		printf("---\n");
		printf("A=[\n");
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c < this->ncols; ++c){
				printf("%10.6f ", M[r][c]);
			}
			printf(";\n");
		}
		printf("];\n");
		printf("---\n");
		printf("b=[\n");
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c < this->dimension; ++c){
				printf("%10.6f ", rhs[this->nrows * c + r]);
			}
			printf(";\n");
		}
		printf("];\n");
		printf("---\n");		
	}
}

double MorphingMatrix::getRelativeError(const std::vector<double> &ref_data) const {
	double diff = 0.0f;
	double n = 0.0f;
	
	std::vector<double> basis_vector( this->ncols );
	std::vector<double> result_vector( this->ncols );
	std::fill(basis_vector.begin(), basis_vector.end(), 0.0f);
	
	std::vector<std::vector<double>> *M = nullptr;
	if(info::mpi::rank == 0){
		M = new std::vector<std::vector<double>>(this->ncols);
	}

	for(auto &el: *M){
		el.resize(this->nrows);
		std::fill(el.begin(), el.end(), 0.0f);
	}
	
	for(esint i = 0; i < this->ncols; ++i){
		basis_vector[i] = 1.0f;
		if(info::mpi::rank == 0){
			this->apply(
				&basis_vector[0],
				&M->at(i)[0],
				1.0f,
				1.0f,
				false
			);
		}
		else{
			this->apply(
				&basis_vector[0],
				&result_vector[0],
				1.0f,
				1.0f,
				false
			);
		}

		basis_vector[i] = 0.0f;
	}
	
	std::vector<double> result(1);
	if(info::mpi::rank == 0){
		esint idx_ = 0;
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c <= r; ++c){
				
				double m = ref_data[idx_++];
				n += m*m;
				
				m = m - M->at(c)[r];
				
				diff += m*m;
				
			}
		}
		result[0] = std::sqrt(diff / n);
	}
	
	
	Communication::broadcastUnknownSize(result);
	
	return result[0];
}

double MorphingMatrix::getRelativeErrorMax(const std::vector<double> &ref_data) const {
	double diff = 0.0f, diff_max = 0.0f;
	double n = 0.0f;
	
	std::vector<double> basis_vector( this->ncols );
	std::vector<double> result_vector( this->ncols );
	std::fill(basis_vector.begin(), basis_vector.end(), 0.0f);
	
	std::vector<std::vector<double>> *M = nullptr;
	if(info::mpi::rank == 0){
		M = new std::vector<std::vector<double>>(this->ncols);
	}

	for(auto &el: *M){
		el.resize(this->nrows);
		std::fill(el.begin(), el.end(), 0.0f);
	}
	
	for(esint i = 0; i < this->ncols; ++i){
		basis_vector[i] = 1.0f;
		if(info::mpi::rank == 0){
			this->apply(
				&basis_vector[0],
				&M->at(i)[0],
				1.0f,
				1.0f,
				false
			);
		}
		else{
			this->apply(
				&basis_vector[0],
				&result_vector[0],
				1.0f,
				1.0f,
				false
			);
		}

		basis_vector[i] = 0.0f;
	}
	
	std::vector<double> result(1);
	if(info::mpi::rank == 0){
		esint idx_ = 0;
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c <= r; ++c){
				
				double m = ref_data[idx_++];
				n += m*m;
				
				m = m - M->at(c)[r];
				
				diff = m*m;
				
				if(diff > diff_max){
					diff_max = diff;
				}
				
			}
		}
		result[0] = diff_max;
	}
	
	
	Communication::broadcastUnknownSize(result);
	
	return result[0];
}

