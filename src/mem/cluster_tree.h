#ifndef SRC_CLUSTER_TREE_H_
#define SRC_CLUSTER_TREE_H_

#include "basis/containers/point.h"
#include "math/math.h"

#include <vector>

/*
	premistit do: src/basic/structures
*/

namespace espreso {

class Cluster{
public:
	
	Cluster();
	
	Cluster(const std::vector<Point> *b);
	
	Cluster(const std::vector<Point> *b, const std::vector<esint> *ci);
	
	~Cluster();
	
	void createLeftChild( const std::vector<esint> *ci );
	
	void createRightChild( const std::vector<esint> *ci );
	
	esint size() const;
	
	Point getCentroid() const;
	
	Point getPoint(esint idx) const;
	
	esint getPointIndexGlobal(esint idx_local ) const;
	
	Cluster* getLeftChild() const;
	
	Cluster* getRightChild() const;
	
	double getRadius() const;
	
	bool hasChildren() const;
	
private:

	void recalculateCentroid();

	void recalculateRadius();
	
	//TODO: USE SMART POINTERS INSTEAD
	const std::vector<Point> * base = nullptr;
	std::vector<esint> *cluster_indices = nullptr;
	
	Point centroid;
	
	Cluster *lChild = nullptr;
	Cluster *rChild = nullptr;
	
	double radius = 0.0f;
};

class ClusterTree{
public:
	ClusterTree(Cluster &r);

	ClusterTree(const std::vector<Point> *b);
	
	~ClusterTree();
	
	void createClusterTree(esint base_cluster_size);
	
	Cluster* getRoot() const;
	
private:

	void doCreateClusterTree(Cluster *C, esint s);

	Cluster *root = nullptr;
	
	double varMatrix[9];
	double eigenvalues[3];
	double eigenvectors[9];
	double normal[3];
};




};//end of namespace espreso

#endif /* SRC_CLUSTER_TREE_H_ */
