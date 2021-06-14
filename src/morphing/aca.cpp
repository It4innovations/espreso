#include "aca.h"

using namespace espreso;

FullRankBlock::FullRankBlock(const Cluster *L, const Cluster *R, const RBFTargetConfiguration &configuration){
	this->nrows = L->size();
	this->ncols = R->size();
	
	// printf("Constructing full matrix of size: %d x %d\n", this->nrows, this->ncols);
	
	esint n = std::max(this->nrows, this->ncols);
	
	this->x_local.resize(n);
	this->y_local.resize(n);
	
	this->data.reserve(this->nrows * this->ncols);
	
	double v;
	
	for(esint r = 0; r < this->nrows; ++r){
		for(esint c = 0; c < this->ncols; ++c){
			v = configuration.function.evaluator->evaluate((L->getPoint(r) - R->getPoint(c)).length());
			this->data.push_back(v);
		}
	}
	
	this->DOF_indices_left.reserve( this->nrows );
	this->DOF_indices_right.reserve( this->ncols );
	
	for(esint i = 0; i < this->nrows; ++i){
		this->DOF_indices_left.push_back(L->getPointIndexGlobal(i));
	}
	for(esint i = 0; i < this->ncols; ++i){
		this->DOF_indices_right.push_back(R->getPointIndexGlobal(i));
	}
}

FullRankBlock::~FullRankBlock(){
	this->data.clear();
	this->DOF_indices_left.clear();
	this->DOF_indices_right.clear();
	this->x_local.clear();
	this->y_local.clear();
}

esint FullRankBlock::size() const {
	return this->nrows * this->ncols;
}

//performs y = this * x * alpha + y * beta
void FullRankBlock::apply(const double * x_global, double *y_global, double alpha, double beta, bool transpose){

	// MATH::vecScale((transpose?this->ncols:this->nrows), beta, y_global);
	
	if(!transpose){
		for(esint i = 0; i < this->ncols; ++i){
			this->x_local[i] = x_global[this->DOF_indices_right[i]];
		}
		std::fill(this->y_local.begin(), this->y_local.begin() + this->nrows, 0.0f);
	}
	else{
		for(esint i = 0; i < this->nrows; ++i){
			this->x_local[i] = x_global[this->DOF_indices_left[i]];
		}
		std::fill(this->y_local.begin(), this->y_local.begin() + this->ncols, 0.0f);
	}
	
	#ifdef HAVE_MKL
	// cblas_dgemv(
		// CblasRowMajor,
		// transpose ? CblasTrans : CblasNoTrans,
		// this->nrows,
		// this->ncols,
		// alpha,
		// &this->data[0],
		// this->ncols,
		// &this->x_local[0],
		// 1,
		// 1.0f,
		// &this->y_local[0],
		// 1
	// );
	#else
		
	if(!transpose){
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c < this->ncols; ++c){
				this->y_local[ r ] += alpha * this->data[c + r * this->ncols] * this->x_local[ c ];
			}
		}
	}
	else{
		for(esint r = 0; r < this->nrows; ++r){
			for(esint c = 0; c < this->ncols; ++c){
				this->y_local[ c ] += alpha * this->data[c + r * this->ncols] * this->x_local[ r ];
			}
		}
	}
	
	#endif
	
	esint idx;
	if(!transpose){
		for(esint i = 0; i < this->nrows; ++i){
			idx = this->DOF_indices_left[i];
			y_global[idx] += y_local[i];
		}
	}
	else{
		for(esint i = 0; i < this->ncols; ++i){
			idx = this->DOF_indices_right[i];
			y_global[idx] += y_local[i];
		}
	}
}


LowRankBlock::LowRankBlock(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration){
	this->nrows = L->size();
	this->ncols = R->size();
	
	esint n = std::max(this->nrows, this->ncols);
	
	
	this->generateBlocksImplicit(L, R, eps, configuration);
	// this->generateBlocksExplicit(L, R, eps, configuration);

	this->x_local.resize(n);
	this->y_local.resize(n);
	this->z_local.resize(this->rank);

	this->DOF_indices_left.reserve( this->nrows );
	this->DOF_indices_right.reserve( this->ncols );
	
	for(esint i = 0; i < this->nrows; ++i){
		this->DOF_indices_left.push_back(L->getPointIndexGlobal(i));
	}
	for(esint i = 0; i < this->ncols; ++i){
		this->DOF_indices_right.push_back(R->getPointIndexGlobal(i));
	}	
}

LowRankBlock::~LowRankBlock(){
	this->data_left.clear();
	this->data_right.clear();
	this->DOF_indices_left.clear();
	this->DOF_indices_right.clear();
	this->x_local.clear();
	this->y_local.clear();
	this->z_local.clear();
}

esint LowRankBlock::size() const {
	return this->data_right.size() + this->data_left.size();
}

//performs y = this * x * alpha + y * beta
void LowRankBlock::apply(const double * x_global, double *y_global, double alpha, double beta, bool transpose){
	
	if(this->data_left.size() + this->data_right.size() == 0){
		return;
	}

	// MATH::vecScale((transpose?this->ncols:this->nrows), beta, y_global);
	
	if(!transpose){
		for(esint i = 0; i < this->ncols; ++i){
			this->x_local[i] = x_global[this->DOF_indices_right[i]];
		}
		std::fill(this->y_local.begin(), this->y_local.begin() + this->nrows, 0.0f);
	}
	else{
		for(esint i = 0; i < this->nrows; ++i){
			this->x_local[i] = x_global[this->DOF_indices_left[i]];
		}
		std::fill(this->y_local.begin(), this->y_local.begin() + this->ncols, 0.0f);
	}
	
	if(this->data_right.size() > 0){
		std::fill(this->z_local.begin(), this->z_local.end(), 0.0f);
	}
	
	#ifdef HAVE_MKL
	// cblas_dgemv(
		// CblasRowMajor,
		// transpose ? CblasTrans : CblasNoTrans,
		// this->rank,
		// this->ncols,
		// alpha,
		// &this->data_right[0],
		// this->rank,
		// &this->x_local[0],
		// 1,
		// 1.0f,
		// &this->z_local[0],
		// 1
	// );

	// cblas_dgemv(
		// CblasRowMajor,
		// transpose ? CblasTrans : CblasNoTrans,
		// this->nrows,
		// this->rank,
		// 1.0f,
		// &this->data_left[0],
		// this->nrows,
		// &this->z_local[0],
		// 1,
		// 1.0f,
		// &this->y_local[0],
		// 1
	// );
	#else
		
	if(this->data_right.size() <= 0){
		if(!transpose){
			for(esint r = 0; r < this->nrows; ++r){
				for(esint c = 0; c < this->ncols; ++c){
					this->y_local[ r ] += alpha * this->data_left[c + r * this->ncols] * this->x_local[ c ];
				}
			}
		}
		else{
			for(esint r = 0; r < this->nrows; ++r){
				for(esint c = 0; c < this->ncols; ++c){
					this->y_local[ c ] += alpha * this->data_left[c + r * this->ncols] * this->x_local[ r ];
				}
			}
		}
	}
	else{
		if(!transpose){
			for(esint k = 0; k < this->rank; ++k){
				for(esint c = 0; c < this->ncols; ++c){
					this->z_local[ k ] += alpha * this->data_right[c + k * this->ncols] * this->x_local[ c ];
				}
			}

			for(esint r = 0; r < this->nrows; ++r){
				for(esint k = 0; k < this->rank; ++k){
					this->y_local[ r ] += this->data_left[k + r * this->rank] * this->z_local[ k ];
				}
			}
		}
		else{
			for(esint k = 0; k < this->rank; ++k){
				for(esint c = 0; c < this->nrows; ++c){
					this->z_local[ k ] += alpha * this->data_left[c + k * this->nrows] * this->x_local[ c ];
				}
			}

			for(esint r = 0; r < this->ncols; ++r){
				for(esint k = 0; k < this->rank; ++k){
					this->y_local[ r ] += this->data_right[k + r * this->rank] * this->z_local[ k ];
				}
			}
		}
	}
	
	#endif

	esint idx;
	if(!transpose){
		for(esint i = 0; i < this->nrows; ++i){
			idx = this->DOF_indices_left[i];
			y_global[idx] += y_local[i];
		}
	}
	else{
		for(esint i = 0; i < this->ncols; ++i){
			idx = this->DOF_indices_right[i];
			y_global[idx] += y_local[i];
		}
	}
}


void LowRankBlock::generateBlocksExplicit(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration){
	std::vector<double> data_full;
	data_full.reserve(this->nrows * this->ncols);
	
	double v;
	double ref_norm_init = 0.0f;
	for(esint r = 0; r < this->nrows; ++r){
		for(esint c = 0; c < this->ncols; ++c){
			v = configuration.function.evaluator->evaluate((L->getPoint(r) - R->getPoint(c)).length());
			data_full.push_back(v);
			ref_norm_init += v * v;
		}
	}

	this->rank = 0;
	
	std::vector<std::vector<double>*> data_rows;
	std::vector<std::vector<double>*> data_cols;
	
	std::vector<bool> assembled_rows(this->nrows);
	std::vector<bool> assembled_cols(this->ncols);
	
	std::fill(assembled_rows.begin(), assembled_rows.end(), false);
	std::fill(assembled_cols.begin(), assembled_cols.end(), false);
	
	double ref_norm_curr = 0.0f;
	while(true){
		esint cross_row_idx = -1;
		esint cross_col_idx = -1;
		double val_pivot = 0.0f;
		
		// first we find a pivotal entry of the reference matrix
		for(esint r = 0; r < this->nrows; ++r){
			if(assembled_rows[r]){
				continue;
			}
			
			for(esint c = 0; c < this->ncols; ++c){
				if(assembled_cols[c]){
					continue;
				}
				
				
				if(std::abs(data_full[c + r * this->ncols]) > val_pivot){
					val_pivot = std::abs(data_full[c + r * this->ncols]);
					cross_row_idx = r;
					cross_col_idx = c;
				}
				
			}
		}
		
		if(cross_row_idx < 0){
			// printf("  error: %20.18f\n", std::sqrt(ref_norm_curr / ref_norm_init));
			break;
		}
		
		
		//assemble the new cross
		std::vector<double> *cross_row = new std::vector<double>(this->ncols);
		std::vector<double> *cross_col = new std::vector<double>(this->nrows);

		for(esint c = 0; c < this->ncols; ++c){
			if(assembled_cols[c]){
				cross_row->at(c) = 0.0f;
			}
			else{
				cross_row->at(c) = data_full[c + cross_row_idx * this->ncols];
			}
		}
		for(esint r = 0; r < this->nrows; ++r){
			if(assembled_rows[r]){
				cross_col->at(r) = 0.0f;
			}
			else{
				cross_col->at(r) = data_full[cross_col_idx + r * this->ncols];
			}
		}
		//normalize the product 
		double v = 1.0f/cross_row->at(cross_col_idx);
		for(esint c = 0; c < this->ncols; ++c){
			if(!assembled_cols[c]){
				cross_row->at(c) *= v;
			}
		}
		
		this->rank++;
		assembled_rows[cross_row_idx] = true;
		assembled_cols[cross_col_idx] = true;
		data_rows.push_back(cross_row);
		data_cols.push_back(cross_col);

		//subtract the product of the new cross from the referential matrix
		for(esint r = 0; r < this->nrows; ++r){
			if(assembled_rows[r]){
				continue;
			}
			
			for(esint c = 0; c < this->ncols; ++c){
				if(assembled_cols[c]){
					continue;
				}
				data_full[c + r * this->ncols] -= cross_col->at(r) * cross_row->at(c);
			}
		}
		
		//updated referential matrix norm
		ref_norm_curr = 0.0f;
		for(esint r = 0; r < this->nrows; ++r){
			if(assembled_rows[r]){
				continue;
			}
			
			for(esint c = 0; c < this->ncols; ++c){
				if(assembled_cols[c]){
					continue;
				}
				double m_ = data_full[c + r*this->ncols];
				ref_norm_curr += m_ * m_;
			}
		}
		
		//terminating condition
		if(std::sqrt(ref_norm_curr) < std::sqrt(ref_norm_init) * eps){
			// printf("  ACA finished with error: %20.18f\n", std::sqrt(ref_norm_curr / ref_norm_init));
			// printf("  rank %6d/%6d\n", this->rank, std::min(this->ncols, this->nrows));
			break;
		}
	}
	
	if(this->rank * (this->ncols + this->nrows) > this->ncols * this->nrows ){
		this->data_right.resize(0);
		this->data_left.resize(this->ncols * this->nrows);
		std::fill(this->data_left.begin(), this->data_left.end(), 0.0f);
		
		for(esint k = 0; k < this->rank; ++k){
			for(esint r = 0; r < this->nrows; ++r){
				for(esint c = 0; c < this->ncols; ++c){
					this->data_left[c + r * this->ncols] += data_cols[k]->at(r) * data_rows[k]->at(c);
				}
			}
		}
	}
	else{
		this->data_right.resize(this->rank * this->ncols);
		this->data_left.resize(this->rank * this->nrows);
		
		for(esint k = 0; k < this->rank; ++k){
			for(esint c = 0; c < this->ncols; ++c){
				this->data_right[c + k * this->ncols] = data_rows[k]->at(c);
			}
		}

		for(esint r = 0; r < this->nrows; ++r){
			for(esint k = 0; k < this->rank; ++k){
				this->data_left[k + r * this->rank] = data_cols[k]->at(r);
			}
		}
	}
	
	
	//cleanup
	for(auto &el: data_rows){
		delete el;
	}
	for(auto &el: data_cols){
		delete el;
	}
}

void LowRankBlock::generateBlocksImplicit(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration){

	this->rank = 0;
	
	std::vector<std::vector<double>*> data_rows;
	std::vector<std::vector<double>*> data_cols;
	
	std::vector<bool> assembled_rows(this->nrows);
	std::vector<bool> assembled_cols(this->ncols);
	
	std::fill(assembled_rows.begin(), assembled_rows.end(), false);
	std::fill(assembled_cols.begin(), assembled_cols.end(), false);
	
	double ref_norm_curr = 0.0f;
	double approx_norm = 0.0f;
	esint cross_row_idx = 0;
	esint cross_col_idx = -1;
	double val_pivot = 0.0f;
	double m;
	
	while(true){

		if(cross_row_idx < 0){
			break;
		}

		std::vector<double> *cross_row = new std::vector<double>(this->ncols);
		for(esint c = 0; c < this->ncols; ++c){
				if( assembled_cols[c] ){
					cross_row->at(c) = 0.0f;
				}
				else{
					cross_row->at(c) = configuration.function.evaluator->evaluate((L->getPoint(cross_row_idx) - R->getPoint(c)).length());
				}
		}
		for( esint i = 0; i < this->rank; ++i){
			for(esint c = 0; c < this->ncols; ++c){
				if( !assembled_cols[c] ){
					cross_row->at(c) -= data_rows[i]->at(c) * data_cols[i]->at(cross_row_idx);
				}
			}
		}

		// find the pivotal column in the prescribed row
		cross_col_idx = -1;
		val_pivot = 0.0f;
		for(esint i = 0; i < this->ncols; ++i){
			if(assembled_cols[i]){
				continue;
			}
			
			m = std::abs(cross_row->at(i));
			if(m > val_pivot){
				val_pivot = m;
				cross_col_idx = i;
			}
		}
		
		if(cross_col_idx < 0){
			//no additional column has been found, terminate
			delete cross_row;
			break;
		}
		
		
		//assemble the new column
		std::vector<double> *cross_col = new std::vector<double>(this->nrows);
		for(esint r = 0; r < this->nrows; ++r){
			if( assembled_rows[r] ){
				cross_col->at(r) = 0.0f;
			}
			else{
				cross_col->at(r) = configuration.function.evaluator->evaluate((L->getPoint(r) - R->getPoint(cross_col_idx)).length());
			}
		}
		
		//subtract the product of already computed crosses
		for( esint i = 0; i < this->rank; ++i){
			for(esint r = 0; r < this->nrows; ++r){
				if(!assembled_rows[r]){
					cross_col->at(r) -= data_rows[i]->at(cross_col_idx) * data_cols[i]->at(r);
				}
			}
		}
	
		//normalize the product 
		val_pivot = 1.0f/cross_row->at(cross_col_idx);
		for(esint c = 0; c < this->ncols; ++c){
			cross_row->at(c) *= val_pivot;
		}
		

		//we check for terminating condition satisfaction
		double norm_new_row = espreso::MATH::vecDot(this->ncols, &cross_row->at(0));
		double norm_new_col = espreso::MATH::vecDot(this->nrows, &cross_col->at(0));
		double new_cross_norm = norm_new_row * norm_new_col;
		
		approx_norm += new_cross_norm;
		for ( esint i = 0; i < this->rank; ++i ) {
		  approx_norm += 2.0f * espreso::MATH::vecDot(this->ncols, &cross_row->at(0), &data_rows[i]->at(0)) * espreso::MATH::vecDot(this->nrows, &cross_col->at(0), &data_cols[i]->at(0));
		}
		
		this->rank++;
		data_rows.push_back(cross_row);
		data_cols.push_back(cross_col);
		assembled_rows[cross_row_idx] = true;
		assembled_cols[cross_col_idx] = true;

		//terminating condition
		if(std::sqrt(new_cross_norm) < std::sqrt(approx_norm) * eps){
			break;
		}


		//we check for the new pivotal row
		cross_row_idx = -1;
		val_pivot = 0.0f;
		for(esint i = 0; i < this->nrows; ++i){
			if(assembled_rows[i]){
				continue;
			}
			
			m = std::abs(cross_col->at(i));
			
			if(m > val_pivot){
				val_pivot = m;
				cross_row_idx = i;
			}
		}
	}
	
	if(this->rank * (this->ncols + this->nrows) > this->ncols * this->nrows ){
		this->data_right.resize(0);
		this->data_left.resize(this->ncols * this->nrows);
		std::fill(this->data_left.begin(), this->data_left.end(), 0.0f);
		
		for(esint k = 0; k < this->rank; ++k){
			for(esint r = 0; r < this->nrows; ++r){
				for(esint c = 0; c < this->ncols; ++c){
					this->data_left[c + r * this->ncols] += data_cols[k]->at(r) * data_rows[k]->at(c);
				}
			}
		}
	}
	else{
		this->data_right.resize(this->rank * this->ncols);
		this->data_left.resize(this->rank * this->nrows);
		
		for(esint k = 0; k < this->rank; ++k){
			for(esint c = 0; c < this->ncols; ++c){
				this->data_right[c + k * this->ncols] = data_rows[k]->at(c);
			}
		}

		for(esint r = 0; r < this->nrows; ++r){
			for(esint k = 0; k < this->rank; ++k){
				this->data_left[k + r * this->rank] = data_cols[k]->at(r);
			}
		}
	}
	
	//cleanup
	for(auto &el: data_rows){
		delete el;
	}
	for(auto &el: data_cols){
		delete el;
	}

}



matrix_ACA::matrix_ACA(){
	
}

matrix_ACA::matrix_ACA(
	const ClusterTree &lT, 
	const ClusterTree &rT, 
	double aca_epsilon, 
	double aca_eta,
	const RBFTargetConfiguration &configuration
){
	this->nrows = lT.getRoot()->size();
	this->ncols = rT.getRoot()->size();
	
	if(this->nrows * this->ncols <= 0){
		return;
	}

	esint nthreads = 1;
	#pragma omp parallel
	{
		nthreads = omp_get_num_threads();
	}
	this->y_tmp.resize(nthreads);
	for(auto &y: this->y_tmp){
		y.resize(std::max(this->nrows, this->ncols));
	}
	this->aca_blocks_threaded.resize(nthreads);
	
	BlockClusterTree T(lT, rT, aca_eta);
	this->aca_blocks.resize(T.leaf_size());

	esint size_nonadmissible_ = 0;
	esint size_admissible_ = 0;
	#pragma omp parallel reduction(+:size_nonadmissible_, size_admissible_)
	{
		esint idx;
		esint size_admissible_tmp = 0;
		esint size_nonadmissible_tmp = 0;
		esint T_s = T.leaf_size();

		#pragma omp for schedule(dynamic, 1) nowait
		for(esint i = 0; i < T_s; ++i){
			const BlockCluster* L = T.get_leaf( i );
			
			const Cluster* lC = L->getLeftCluster();
			const Cluster* rC = L->getRightCluster();
			
			if(L->getIsAdmissible()){
				this->aca_blocks[i] = new LowRankBlock(lC, rC, aca_epsilon, configuration);
				size_admissible_tmp += this->aca_blocks[i]->size();
			}
			else{
				this->aca_blocks[i] = new FullRankBlock(lC, rC, configuration);
				size_nonadmissible_tmp += this->aca_blocks[i]->size();
			}
		}
		size_nonadmissible_ = size_nonadmissible_tmp;
		size_admissible_ = size_admissible_tmp;
	}
	this->size_nonadmissible = size_nonadmissible_;
	this->size_admissible = size_admissible_;

	std::vector<esint> aca_blocks_size_threaded(nthreads);
	std::fill(aca_blocks_size_threaded.begin(), aca_blocks_size_threaded.end(), 0);
	for(auto M: this->aca_blocks){
		esint s = M->size();

		esint min_s = aca_blocks_size_threaded[0];
		esint bucket_idx = 0;
		for(esint i = 1; i < aca_blocks_size_threaded.size(); ++i){
			if(min_s > aca_blocks_size_threaded[i]) {
				min_s = aca_blocks_size_threaded[i];
				bucket_idx = i;
			}
		}

		this->aca_blocks_threaded[bucket_idx].push_back(M);
		aca_blocks_size_threaded[bucket_idx] += M->size();
	}
}

matrix_ACA::~matrix_ACA(){
	for(auto el: this->aca_blocks){
		delete el;
	}
}

void matrix_ACA::apply(const double* x, double* y, double alpha, double beta, bool transpose)  {

	if(this->getNEntries() <= 0){
		return;
	}

	esint dim = (transpose?this->ncols:this->nrows);
	MATH::vecScale(dim, beta, y);
	
	#pragma omp parallel
  	{
		esint tid = omp_get_thread_num( );
		MATH::vecScale(dim, 0.0, this->y_tmp[tid].data());

		for(auto &M: this->aca_blocks_threaded[tid]){
			M->apply(x, &this->y_tmp[tid][0], alpha, 1.0f, transpose);
		}
	}
	
	for(esint i = 0; i < this->y_tmp.size(); ++i){
		MATH::vecAdd(dim, y, 1.0, this->y_tmp[i].data());
	}
}

double matrix_ACA::getCompressionRatio() const {
	double size_real = this->size_admissible + this->size_nonadmissible;
	double size_max = this->nrows * this->ncols;
	
	// printf("ACA matrix: real size: %f, max size: %f\n", size_real, size_max);
	return size_real / size_max;
}

esint matrix_ACA::getNEntries() const {
	return this->size_admissible + this->size_nonadmissible;
}