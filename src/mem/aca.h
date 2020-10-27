#ifndef SRC_ACA_MATRIX_H_
#define SRC_ACA_MATRIX_H_

#include "cluster_tree.h"
#include "block_cluster_tree.h"

#include "esinfo/ecfinfo.h"
#include "basis/evaluator/evaluator.h"


namespace espreso {

struct RBFTargetConfiguration;
struct RBFTargetTransformationConfiguration;
	
class FullRankBlock{
public:

	FullRankBlock(const Cluster *L, const Cluster *R, const RBFTargetConfiguration &configuration);
	
	~FullRankBlock();
	
	esint size() const;
	
	//performs y = this * x * alpha + y * beta
	void apply(const double * x_global, double *y_global, double alpha, double beta, bool transpose);
	
private:
	std::vector<double> data;

	std::vector<esint> DOF_indices_left;
	std::vector<esint> DOF_indices_right;
	
	std::vector<double> x_local;
	std::vector<double> y_local;
	
	esint nrows = 0;
	esint ncols = 0;
};

class LowRankBlock{
public:

	LowRankBlock(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration);
	
	~LowRankBlock();
	
	esint size() const;
	
	//performs y = this * x * alpha + y * beta
	void apply(const double * x_global, double *y_global, double alpha, double beta, bool transpose);
	
private:

	void generateBlocks(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration);

	void generateBlocksExplicit(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration);

	void generateBlocksImplicit(const Cluster *L, const Cluster *R, double eps, const RBFTargetConfiguration &configuration);

	std::vector<double> data_left;
	std::vector<double> data_right;

	std::vector<esint> DOF_indices_left;
	std::vector<esint> DOF_indices_right;
	
	std::vector<double> x_local;
	std::vector<double> y_local;
	std::vector<double> z_local;
	
	esint nrows = 0;
	esint ncols = 0;
	esint rank = 0;
};

class matrix_ACA{
public:
	
	matrix_ACA();

	matrix_ACA(
		const ClusterTree &lT, 
		const ClusterTree &rT, 
		double aca_epsilon, 
		double aca_eta,
		const RBFTargetConfiguration &configuration
	);
	
	~matrix_ACA();
	
	void apply(const double* x, double* y, double alpha, double beta, bool transpose) const;
	
	double getCompressionRatio() const;
	
	esint getNEntries() const;
	
private:

	void createAdmissibleBlock(const Cluster* L, const Cluster* R, double eps, const RBFTargetConfiguration &configuration);
	
	void createNonAdmissibleBlock(const Cluster* L, const Cluster* R, const RBFTargetConfiguration &configuration);
	
	std::vector<LowRankBlock*> admissible_blocks;
	std::vector<FullRankBlock*> nonadmissible_blocks;
	
	esint size_admissible = 0;
	esint size_nonadmissible = 0;
	
	esint nrows = 0;
	esint ncols = 0;
};



};//end of namespace espreso

#endif /* SRC_ACA_MATRIX_H_ */