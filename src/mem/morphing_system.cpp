
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
	
}

MorphingMatrix::MorphingMatrix(
	const std::vector<Point> &sPoints,
	const std::vector<double> &sDisplacement,
	const RBFTargetConfiguration &configuration,
	bool use_x, 
	bool use_y, 
	bool use_z,
	double aca_eps,
	double aca_eta,
	esint base_cluster_size,
		
	std::vector<double> &morphing_rhs
){
	
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
	
	// if(info::mpi::rank == 0){
		// for(int r = 0; r < nranks; ++r){
			// for(int c = 0; c < nranks; ++c){
				// printf("[%2d]", block_distribution[c + r * nranks]);
			// }
			// printf("\n");
		// }
		
		// printf("\n");
		// for(int r = 0; r < nranks; ++r){
			// for(int c = 0; c < block_processes[0].size(); ++c){
				// printf("[%2d]", block_processes[r][c]);
			// }
			// printf("\n");
		// }

		// printf("\n");
		// for(int r = 0; r < nranks; ++r){
			// for(int c = 0; c < process_blocks[0].size(); ++c){
				// printf("[%2d]", process_blocks[r][c]);
			// }
			// printf("\n");
		// }
	// }
	// printf("[%2d] active block indices: ", info::mpi::rank);
	// for(int r = 0; r < nranks; ++r){
		// if(active_block_indices[r]){
			// printf("[%2d]", r);
		// }
	// }
	// printf("\n");

	std::vector<std::vector<Point>> active_point_set(nranks);
	active_point_set[info::mpi::rank] = sPoints;

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
	}
	MPI_Waitall(req_idx, requests.data(), MPI_STATUSES_IGNORE);
	// printf("Number of requests: %3d\n", req_idx);

	
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
	// this->x_local.resize(nranks);
	// this->y_local.resize(nranks);
	esint npoints_total = 0;
	
	for( int i = 0; i < nranks; ++i ){
		
		global_DOFs_start[i] = npoints_total;
		global_DOFs_end[i] = global_DOFs_start[i] + nPoints[ i ];
		
		// if(active_block_indices[i]){
			// this->x_local[i].resize(nPoints[ i ]);
			// this->y_local[i].resize(nPoints[ i ]);
		// }

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
	
	this->dimension = (use_x?1:0) + (use_y?1:0) + (use_z?1:0);
	this->ncols = npoints_total + this->dimension + 1;
	this->nrows = this->ncols;
	
	// printf(" # of points: %d, dim: %d, #of cols: %d\n", npoints_total, this->dimension, this->nrows);
	
	morphing_rhs.resize(this->nrows * this->dimension);
	std::fill(morphing_rhs.begin(), morphing_rhs.end(), 0.0f);
	size_t shift_rhs = this->global_DOFs_start[info::mpi::rank];
	for (int d = 0; d < this->dimension; d++) {
		for (size_t r = 0; r < active_point_set[info::mpi::rank].size(); r++) {
			morphing_rhs[this->ncols * d + r + shift_rhs] = sDisplacement[r * this->dimension + d];
		}
	}	
	Communication::allReduce(morphing_rhs, Communication::OP::SUM);
	
	this->recalculateB(
		active_point_set[info::mpi::rank],
		use_x,
		use_y,
		use_z
	);
	
	for(auto &el:trees){
		if(el){
			delete el;
			el = nullptr;
		}
	}
	
	double compr = this->getCompressionRatio();
	if(info::mpi::rank == 0){
		printf("ACA morphing matrix constructed, compression ratio %f [%%]\n", 100.0f * compr);
	}
}

void MorphingMatrix::recalculateB(
	const std::vector<Point> &sPoints,
	bool use_x, 
	bool use_y, 
	bool use_z
){
	this->B.resize(0);
	this->dimension = (use_x?1:0) + (use_y?1:0) + (use_z?1:0);
	this->B.reserve(sPoints.size() * (this->dimension + 1));
	
	if (use_x) {
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back(sPoints[c].x);
		}
	}

	if (use_y) {
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back(sPoints[c].y);
		}
	}

	if (use_z) {
		for (size_t c = 0; c < sPoints.size(); c++) {
			this->B.push_back(sPoints[c].z);
		}
	}

	for (size_t c = 0; c < sPoints.size(); c++) {
		this->B.push_back(1.0f);
	}
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
	esint nr = this->nrows - this->dimension - 1;
	
	shift_x = this->global_DOFs_start[info::mpi::rank];
	shift_y = this->global_DOFs_start[info::mpi::rank];
	for(esint r = 0; r < this->dimension + 1; ++r){
		for(esint c = 0; c < nc; ++c){
			//B
			y_global[nr + r] += alpha * this->B[c + r * nc] * x_global[c + shift_x];
			
			//B'
			y_global[c + shift_y] += alpha * this->B[c + r * nc] * x_global[nr + r];
		}
	}
	
	#endif
	
	Communication::allReduce(MPI_IN_PLACE, y_global, this->nrows, MPI_DOUBLE, MPI_SUM);
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

void MorphingMatrix::print() const {
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
			1.0f,
			false
		);

		basis_vector[i] = 0.0f;
	}
	
	if(info::mpi::rank == 0){
		printf("---\n");
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c <= r; ++c){
				printf("%10.6f ", M[c][r]);
			}
			printf("\n");
		}
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

