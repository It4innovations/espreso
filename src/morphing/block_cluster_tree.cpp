
#include "block_cluster_tree.h"

using namespace espreso;

BlockCluster::BlockCluster(){}

BlockCluster::~BlockCluster(){
    if(this->child_LL){
        delete this->child_LL;
        this->child_LL = nullptr;
    }
    if(this->child_LR){
        delete this->child_LR;
        this->child_LR = nullptr;
    }
    if(this->child_RL){
        delete this->child_RL;
        this->child_RL = nullptr;
    }
    if(this->child_RR){
        delete this->child_RR;
        this->child_RR = nullptr;
    }
    if(this->child_lL){
        delete this->child_lL;
        this->child_lL = nullptr;
    }
    if(this->child_lR){
        delete this->child_lR;
        this->child_lR = nullptr;
    }
    if(this->child_Lr){
        delete this->child_Lr;
        this->child_Lr = nullptr;
    }
    if(this->child_Rr){
        delete this->child_Rr;
        this->child_Rr = nullptr;
    }
}

BlockCluster::BlockCluster(const Cluster* l, const Cluster* r, double eta_){
    this->left_cluster = l;
    this->right_cluster = r;
    
    this->eta = eta_;
    this->setIsAdmissible(this->isAdmissible(this->eta));
}

bool BlockCluster::getIsAdmissible() const {
    return this->is_admissible;
}

void BlockCluster::setIsAdmissible(bool v){
    this->is_admissible = v;
}

const Cluster* BlockCluster::getLeftCluster() const {
    return this->left_cluster;
}

const Cluster* BlockCluster::getRightCluster() const {
    return this->right_cluster;
}

esint BlockCluster::createChildren(){
    
    if(this->left_cluster->hasChildren()){
        if(this->right_cluster->hasChildren()){
            this->child_LL = new BlockCluster(this->left_cluster->getLeftChild(), this->right_cluster->getLeftChild(), this->eta);
            this->child_LR = new BlockCluster(this->left_cluster->getLeftChild(), this->right_cluster->getRightChild(), this->eta);
            this->child_RL = new BlockCluster(this->left_cluster->getRightChild(), this->right_cluster->getLeftChild(), this->eta);
            this->child_RR = new BlockCluster(this->left_cluster->getRightChild(), this->right_cluster->getRightChild(), this->eta);
            return 4;
        }
        else{
            this->child_Lr = new BlockCluster(this->left_cluster->getLeftChild(), this->right_cluster, this->eta);
            this->child_Rr = new BlockCluster(this->left_cluster->getRightChild(), this->right_cluster, this->eta);
            return 2;
        }
    }
    else{
        if(this->right_cluster->hasChildren()){
            this->child_lL = new BlockCluster(this->left_cluster, this->right_cluster->getLeftChild(), this->eta);
            this->child_lR = new BlockCluster(this->left_cluster, this->right_cluster->getRightChild(), this->eta);
            return 2;
        }
    }
    
    return 0;
}

bool BlockCluster::hasChildren(){
    return (this->child_LL || this->child_LR || this->child_RL || this->child_RR || this->child_lL || this->child_lR || this->child_Lr || this->child_Rr);
}

BlockCluster* BlockCluster::getChild(int li, int ri){
    if(li == 0 && ri == 0){
        return this->child_LL;
    }
    
    if(li == 0 && ri == 1){
        return this->child_LR;
    }
    
    if(li == 1 && ri == 0){
        return this->child_RL;
    }
    
    if(li == 1 && ri == 1){
        return this->child_RR;
    }
    
    return nullptr;
}


BlockCluster* BlockCluster::getChild(int idx){

    if(this->child_lL){
        
        if(idx == 0){
            return this->child_lL;
        }
        
        if(idx == 1){
            return this->child_lR;
        }
        
    }

    if(this->child_Lr){
        
        if(idx == 0){
            return this->child_Lr;
        }
        if(idx == 1){
            return this->child_Rr;
        }
        
    }
    
    return nullptr;
}

bool BlockCluster::isAdmissible(double eta) const{
    bool is_admissible_ = false;
    
    const Cluster *lC = this->getLeftCluster();
    const Cluster *rC = this->getRightCluster();
    
    const Point lCentroid = lC->getCentroid();
    const Point rCentroid = rC->getCentroid();
    Point D = lCentroid - rCentroid;
    
    double lRadius = lC->getRadius();
    double rRadius = rC->getRadius();
    
    double centroid_distance = D.length();
    
    // we check the admisibility condition, if the clusters are far from each other, we may approximate the relations
    double distance_ = centroid_distance;//std::max( centroid_distance - lRadius - rRadius, (double)0.0f );//0 if close to overlap
    double diam_coef = 2.0f * std::max(lRadius, rRadius);

    if(diam_coef < eta * distance_){
        is_admissible_ = true;
    }
    
    return is_admissible_;
}

BlockClusterTree::BlockClusterTree(
    const ClusterTree &lT,
    const ClusterTree &rT,
    double eta
){
    this->root = new BlockCluster(lT.getRoot(), rT.getRoot(), eta);
    esint e = this->createTree( this->root, 0 );
    
    if(e == 0){
        this->leaf_clusters.push_back( this->root );
    }
}

BlockClusterTree::~BlockClusterTree(){
    if(this->root){
        delete this->root;
        this->root = nullptr;
    }
}

esint BlockClusterTree::leaf_size() const {
    return this->leaf_clusters.size();
}

const BlockCluster* BlockClusterTree::get_leaf(esint idx){
    return this->leaf_clusters.at( idx );
}


esint BlockClusterTree::createTree( BlockCluster* C, esint depth ){
    
    if(C->getIsAdmissible()){
        // admissible leaf
        this->leaf_clusters.push_back( C );

        return 1;
    }

    esint n_children = C->createChildren();
    
    if(n_children == 0){
        return 0;
    }
    
    if(n_children == 2){
        esint res_L = this->createTree( C->getChild(0), depth + 1 );
        esint res_R = this->createTree( C->getChild(1), depth + 1 );
        
        if(res_L + res_R == 0){
            // we let the parent decide what to do with C
            return 0;
        }

        // one child has been added as a leaf, the other was not, we cannot let the parent decide what to do with C
        if(res_L == 0){
            this->leaf_clusters.push_back( C->getChild(0) );
        }
        
        if(res_R == 0){
            this->leaf_clusters.push_back( C->getChild(1) );
        }
        
        return 1;
    }
    else{
        esint res_LL = this->createTree( C->getChild(0, 0), depth + 1 );
        esint res_LR = this->createTree( C->getChild(0, 1), depth + 1 );
        esint res_RL = this->createTree( C->getChild(1, 0), depth + 1 );
        esint res_RR = this->createTree( C->getChild(1, 1), depth + 1 );
        
        if(res_LL + res_LR + res_RL + res_RR == 0){
            // we let the parent decide what to do with C
            return 0;
        }

        // at least one child has been added as a leaf, we cannot let the parent decide what to do with C
        if(res_LL == 0){
            this->leaf_clusters.push_back( C->getChild(0, 0) );
        }
        
        if(res_LR == 0){
            this->leaf_clusters.push_back( C->getChild(0, 1) );
        }
        
        if(res_RL == 0){
            this->leaf_clusters.push_back( C->getChild(1, 0) );
        }
        
        if(res_RR == 0){
            this->leaf_clusters.push_back( C->getChild(1, 1) );
        }
        
        return 1;
    }
}


