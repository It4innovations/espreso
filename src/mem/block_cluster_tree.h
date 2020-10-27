#ifndef SRC_BLOCK_CLUSTER_TREE_H_
#define SRC_BLOCK_CLUSTER_TREE_H_

#include "cluster_tree.h"

#include "basis/containers/point.h"
#include "math/math.h"

namespace espreso {

class BlockCluster{
public:
	
	BlockCluster();
	
	BlockCluster(const Cluster* l, const Cluster* r, double eta);
	
	bool getIsAdmissible() const;

	void setIsAdmissible(bool v);

	const Cluster* getLeftCluster() const;

	const Cluster* getRightCluster() const;
	
	esint createChildren();
	
	bool hasChildren();
	
	BlockCluster* getChild(int li, int ri);

	BlockCluster* getChild(int idx);
	

private:

	bool isAdmissible(double eta) const;

	const Cluster *left_cluster = nullptr;
	const Cluster *right_cluster = nullptr;
	
	//this block cluster is subdivided into 4 sub-blocks
	BlockCluster* child_LL = nullptr;
	BlockCluster* child_LR = nullptr;
	BlockCluster* child_RL = nullptr;
	BlockCluster* child_RR = nullptr;

	//this block cluster is subdivided into 2 sub-blocks, LEFt cluster is not subdivided
	BlockCluster* child_lL = nullptr;
	BlockCluster* child_lR = nullptr;
	
	//this block cluster is subdivided into 2 sub-blocks, RIGHT cluster is not subdivided
	BlockCluster* child_Lr = nullptr;
	BlockCluster* child_Rr = nullptr;
	
	bool is_admissible = false;
	
	double eta = 0.0f;
	
};

class BlockClusterTree{
public:
	
	BlockClusterTree(
		const ClusterTree &lT,
		const ClusterTree &rT,
		double eta
	);
	
	esint leaf_size() const;
	
	const BlockCluster* get_leaf(esint idx);
	
private:

	esint createTree( BlockCluster* C, esint depth = 0 );

	BlockCluster *root = nullptr;
	
	std::vector<BlockCluster*> leaf_clusters;

};





};//end of namespace espreso

#endif /* SRC_BLOCK_CLUSTER_TREE_H_ */
