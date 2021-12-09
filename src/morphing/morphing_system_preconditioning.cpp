#include "morphing_system_preconditioning.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

MorphingMatrixPreconditioner::MorphingMatrixPreconditioner(){


}

MorphingMatrixPreconditioner::~MorphingMatrixPreconditioner(){


}

MorphingMatrixPreconditioner::MorphingMatrixPreconditioner(
    esint npoints_total,
    esint nregularization,
    esint dof_shift,
    const std::vector<Point> &points_local,
    const RBFTargetConfiguration &configuration
){

    //find a sparse representation of the local set
    esint nsparse_points = 0;
    if(points_local.size() > 0){
        nsparse_points = std::ceil(std::log2(points_local.size()));
    }
	nsparse_points = 2;

    KDTree<double> T(points_local.data(), points_local.data() + points_local.size(), nsparse_points);
    esint nbuckets = T.getNLeaves();

    std::vector<esint> sparse_point_indices_global(nbuckets);
    std::vector<Point> sparse_points_local(nbuckets);
    std::vector<esint> local_points_in_sparse_set(points_local.size());

    std::fill(local_points_in_sparse_set.begin(), local_points_in_sparse_set.end(), -1);
    for(esint i = 0; i < nbuckets; ++i){
        esint point_index = T.getMidPoint(i);
        local_points_in_sparse_set[point_index] = i;

        sparse_point_indices_global[i] = point_index + dof_shift;
        sparse_points_local[i] = points_local[point_index];
    }
    // printf("%d %d\n", sparse_point_indices_global[0], sparse_point_indices_global[1]);
    Communication::allGatherUnknownSize<esint>(sparse_point_indices_global);
    Communication::allGatherUnknownSize<Point>(sparse_points_local);

	//recieving foreign data
	// printf("sparse_indices = [");
	// for(auto el: sparse_point_indices_global){
	// 	printf("  %d, ", el);
	// }
	// printf("];\n");
	
   
	printf("# of sparse points: %d, total # of local points: %d\n", (int)sparse_points_local.size(), (int)points_local.size());
	// printf("%d %d\n", sparse_point_indices_global[0], sparse_point_indices_global[1]);

    this->nrows = sparse_points_local.size() + 1;

    this->prepareMatrixSparse(sparse_points_local, configuration);

    //for each local point, we assemble local sparse RBF systems and compute the appropriate inverse functions
    for(esint i = 0; i < (esint)points_local.size(); ++i){
        esint global_dof_row = i + dof_shift;

		// printf("Point %3d\n", i);
        if(local_points_in_sparse_set[i] >= 0){
			//no need to modify the sparse system
            this->updateInverseRowSparseOnly(local_points_in_sparse_set[i]);
        }
        else{
			//we update the first column of the sparse system
            this->updateInverseRow(
				points_local[i],
				sparse_points_local,
				configuration
			);
			
			this->row_indices.push_back(global_dof_row);
			this->col_indices.push_back(global_dof_row);
			this->values.push_back(this->inverse_row[0]);
        }

		for(esint j = 0; j < (esint)sparse_point_indices_global.size(); ++j){
			this->row_indices.push_back(global_dof_row);
			this->col_indices.push_back(sparse_point_indices_global[j]);
			this->values.push_back(this->inverse_row[j + 1]);
		}
    }
	
	for(int i = 0; i < nregularization; ++i){
		this->row_indices.push_back(npoints_total - i - 1);
		this->col_indices.push_back(npoints_total - i - 1);
		this->values.push_back(1.0f / info::mpi::size);
	}

    this->ncols = this->nrows;
	this->dim_id = npoints_total;

    this->nreg_dim = nregularization;
}

void MorphingMatrixPreconditioner::updateInverseRowSparseOnly(esint idx){

	std::fill(this->basis_vector.begin(), this->basis_vector.end(), 0.0f);
	this->basis_vector[idx + 1] = 1.0f;

	std::copy(
		this->sparse_system.begin(), 
		this->sparse_system.end(),
        this->sparse_system_mem.begin()
	);

	MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
		this->nrows - 1,  
		&this->sparse_system_mem[0],
		1, 
		&this->basis_vector[1]
	);

	std::copy(
		this->basis_vector.begin(), 
		this->basis_vector.end(),
        this->inverse_row.begin()
	);
}

void MorphingMatrixPreconditioner::updateInverseRow(
	const Point & p,
	const std::vector<Point> &sparse_points,
	const RBFTargetConfiguration &configuration
){
	std::fill(this->basis_vector.begin(), this->basis_vector.end(), 0.0f);
	this->basis_vector[0] = 1.0f;
	
	// std::copy(
		// this->sparse_system.begin() + this->nrows, 
		// this->sparse_system.end(),
        // this->sparse_system_mem.begin() + this->nrows
	// );

	// this->sparse_system_mem[0] = 0.0f;
    // for (size_t r = 0; r < sparse_points.size(); r++) {
		// this->sparse_system_mem[r + 1] = configuration.function.evaluator->evaluate((p - sparse_points[r]).length());
    // }
	
	esint l = 0, k = 0;
	this->sparse_system_mem[0] = 0.0f;
	for(int r = 1; r < this->nrows; ++r){
		this->sparse_system_mem[l++] = configuration.function.evaluator->evaluate((p - sparse_points[r-1]).length());
		
		for(int c = 1; c <= r; ++c){
			this->sparse_system_mem[l++] = this->sparse_system[k++];
		}
	}	

	MATH::SOLVER::directUpperSymetricIndefiniteColumnMajor(
		this->nrows,  
		&this->sparse_system_mem[0],
		1, 
		&this->basis_vector[0]
	);
	
	std::copy(
		this->basis_vector.begin(), 
		this->basis_vector.end(),
        this->inverse_row.begin()
	);
}

	

void MorphingMatrixPreconditioner::prepareMatrixSparse(
    const std::vector<Point> &sparse_points,
    const RBFTargetConfiguration &configuration
) {

    size_t rows = sparse_points.size() + 1;

    this->inverse_row.resize(rows);
	this->basis_vector.resize(rows);

	esint ssize = ((rows + 1)*rows) / 2;
    this->sparse_system.resize(ssize - rows);
    this->sparse_system_mem.resize(ssize);
	
	// printf("Preparation, system size for %d points: %d\n", rows, ssize);

    size_t idx = 0;
    for (size_t r = 0; r < rows - 1; r++) {
        for(size_t c = 0; c <= r; c++) {
            this->sparse_system[idx++] = configuration.function.evaluator->evaluate((sparse_points[r] - sparse_points[c]).length());
        }
    }
	
	// idx = 0;
	// for(int r = 0; r < this->nrows; ++r){
		// for(int c = 0; c <= r; ++c){
			// printf("%10.8f ", this->sparse_system[idx++]);
		// }
		// printf("\n");
	// }	
}

void MorphingMatrixPreconditioner::printData() const {

	printf("prec_row_indices = [\n");
	for(auto &v: this->row_indices){
		printf("%d ", v);
	}
	printf("];\n");
	
	printf("prec_col_indices = [\n");
	for(auto &v: this->col_indices){
		printf("%d ", v);
	}
	printf("];\n");
	
	printf("prec_values = [\n");
	for(auto &v: this->values){
		printf("%f ", v);
	}
	printf("];\n");
	
	std::vector<double> basis_vector( this->dim_id + this->nreg_dim );
	std::fill(basis_vector.begin(), basis_vector.end(), 0.0f);
	
	std::vector<std::vector<double>> M(this->dim_id + this->nreg_dim);
	for(auto &el: M){
		el.resize(this->dim_id + this->nreg_dim);
		std::fill(el.begin(), el.end(), 0.0f);
	}
	
	for(esint i = 0; i < this->dim_id + this->nreg_dim; ++i){
		basis_vector[i] = 1.0f;

		this->apply(
			&basis_vector[0],
			&M[i][0],
			1.0f,
			0.0f,
			true
		);

		basis_vector[i] = 0.0f;
	}
	
	if(info::mpi::rank == 0){
		printf("---\n");
		printf("P=[\n");
		for(esint r = 0; r < this->dim_id + this->nreg_dim; ++r){
			for(esint c = 0; c < this->dim_id + this->nreg_dim; ++c){
				printf("%10.6f ", M[r][c]);
			}
			printf(";\n");
		}
		printf("];\n");
		printf("---\n");		
	}

}
	
//performs y = this * x * alpha + y * beta
void MorphingMatrixPreconditioner::apply(
    const double * x_global, 
    double *y_global, 
    double alpha, 
    double beta, 
    bool transpose
) const {

    // for(esint i = 0; i < this->dim_id + this->nreg_dim; ++i){
        // y_global[i] = x_global[i];
    // }
	
	
	MATH::vecScale(this->dim_id + this->nreg_dim, beta / info::mpi::size, y_global);
	
	size_t ri, ci;
	double v;
	if(transpose){
		for(size_t i = 0; i < this->row_indices.size(); ++i){
			ri = this->col_indices[i];
			ci = this->row_indices[i];
			v = this->values[i];
			
			y_global[ri] += v * alpha * x_global[ci];
		}
	}
	else{
		for(size_t i = 0; i < this->row_indices.size(); ++i){
			ri = this->row_indices[i];
			ci = this->col_indices[i];
			v = this->values[i];
			
			y_global[ri] += v * alpha * x_global[ci];
		}
	}
	
	Communication::allReduce(MPI_IN_PLACE, y_global, this->dim_id + this->nreg_dim, MPI_DOUBLE, MPI_SUM);
}

esint MorphingMatrixPreconditioner::getNRows() const {
    return this->dim_id + this->nreg_dim;
}

esint MorphingMatrixPreconditioner::getNCols() const{
    return this->dim_id + this->nreg_dim;
}
