#ifndef SRC_MORPHING_SYSTEM_MATRIX_H_
#define SRC_MORPHING_SYSTEM_MATRIX_H_

#include "aca.h"
#include "CyclicGraphDecomposition.h"

#include "esinfo/mpiinfo.h"
#include "CyclicGraphDecomposition.h"

namespace espreso {

/*
	represents a system matrix consisting of 2 blocks: A, B
	-----
	AB'
	B
	-----
*/
class MorphingMatrix{
public:

	MorphingMatrix();
	
	~MorphingMatrix();
	
	MorphingMatrix(
		const std::vector<Point> &rPoints,
		const std::vector<double> &rDisplacement,
		const RBFTargetConfiguration &configuration,
		bool use_x, 
		bool use_y, 
		bool use_z,
		double aca_eps,
		double aca_eta,
		esint base_cluster_size,
		
		std::vector<double> &morphing_rhs
	);
	
	/*
	computes B in row-major format
	*/
	void recalculateB(
		const std::vector<Point> &rPoints,
		bool use_x, 
		bool use_y, 
		bool use_z
	);

	//performs y = this * x * alpha + y * beta
	void apply(
		const double * x_global, 
		double *y_global, 
		double alpha, 
		double beta, 
		bool transpose
	) const ;
	
	esint getNRows() const;

	esint getNCols() const;
	
	double getCompressionRatio() const;
	
	void print() const;
	
	double getRelativeError(const std::vector<double> &ref_data) const;
	
	double getRelativeErrorMax(const std::vector<double> &ref_data) const;
	
private:

	std::vector<matrix_ACA*> A;
	std::vector<esint> A_inner_block_idx;
	std::vector<esint> A_outer_block_idx;

	std::vector<double> B;
	
	
	std::vector<esint> global_DOFs_start;
	std::vector<esint> global_DOFs_end;
	
	// std::vector<std::vector<double>> x_local;
	// std::vector<std::vector<double>> y_local;
	
	esint nrows = 0;
	esint ncols = 0;
	esint dimension = 0;

};


};//end of namespace espreso

#endif /* SRC_MORPHING_SYSTEM_MATRIX_H_ */
