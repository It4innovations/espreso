
#include "cluster_tree.h"

using namespace espreso;

Cluster::Cluster(){
    
}

Cluster::Cluster(const std::vector<Point> *b){
    this->base = b;
    this->cluster_indices = new std::vector<esint>();
    this->cluster_indices->reserve( b->size() );

    for(esint i = 0; i < (esint)b->size(); ++i){
        this->cluster_indices->push_back( i );
    }
    
    this->recalculateCentroid();
    this->recalculateRadius();
}

Cluster::Cluster(const std::vector<Point> *b, const std::vector<esint> *ci){
    this->base = b;
    this->cluster_indices = new std::vector<esint>(*ci);

    this->recalculateCentroid();
    this->recalculateRadius();
}

Cluster::~Cluster(){
    
    if( this->base ){
        
    }
    
    if( this->cluster_indices ){
        delete this->cluster_indices;
        this->cluster_indices = nullptr;
    }

    if( this->lChild ){
        delete this->lChild;
        this->lChild = nullptr;
    }

    if( this->rChild ){
        delete this->rChild;
        this->rChild = nullptr;
    }
}

void Cluster::createLeftChild( const std::vector<esint> *ci ) {
    
    if( !ci ){
        return;
    }
    
    if( this->lChild ){
        delete this->lChild;
    }
    
    this->lChild = new Cluster( this->base, ci );
}

void Cluster::createRightChild( const std::vector<esint> *ci ) {
    
    if( !ci ){
        return;
    }
    
    if( this->rChild ){
        delete this->rChild;
    }
    
    this->rChild = new Cluster( this->base, ci );
}

esint Cluster::size() const {
    return this->cluster_indices->size();
}

Point Cluster::getCentroid() const {
    return this->centroid;
}

Point Cluster::getPoint(esint idx) const {
    return this->base->at( this->cluster_indices->at(idx) );
}

esint Cluster::getPointIndexGlobal(esint idx_local ) const {
    return this->cluster_indices->at( idx_local );
}

Cluster* Cluster::getLeftChild() const {
    return this->lChild;
}

Cluster* Cluster::getRightChild() const {
    return this->rChild;
}

double Cluster::getRadius() const {
    return this->radius;
}

bool Cluster::hasChildren() const {
    return (this->lChild && this->rChild);
}



void Cluster::recalculateCentroid(){
    for(auto &el: *this->cluster_indices){
        this->centroid += this->base->at(el);
    }
    
    esint n = this->size();
    if(n > 0){
        this->centroid /= n;
    }
}

void Cluster::recalculateRadius(){
    this->radius = 0.0f;
    
    Point d;
    double m_;
    for(esint i = 0; i < this->size(); ++i){
        d = this->centroid - this->getPoint( i );
        m_ = d.length();
        
        if( m_ > this->radius ){
            this->radius = m_;
        }
    }
    
}


ClusterTree::ClusterTree(const std::vector<Point> *b){
    this->root = new Cluster(b);
}

ClusterTree::ClusterTree(Cluster &r){
    this->root = &r;
}

ClusterTree::~ClusterTree(){
    if( this->root ){
        delete this->root;
        this->root = nullptr;
    }
}

void ClusterTree::createClusterTree(esint base_cluster_size) {
    this->doCreateClusterTree(this->root, base_cluster_size);
}

Cluster* ClusterTree::getRoot() const {
    return this->root;
}


void ClusterTree::doCreateClusterTree(Cluster *C, esint s) {
    
    esint n = C->size();
    
    if(n <= s){
        return;//current cluster C is small enough
    }
    
    //split the current cluster C into two sub-clusters
    Point centroid = C->getCentroid();
    
    // create the covariance matrix of the cluster in the Upper form
    std::fill( &this->varMatrix[0], &this->varMatrix[0] + 9, 0.0f );
    
    esint idx;
    double v1, v2, v3;
    for ( auto i = 0; i < n; i++ ) {
        idx = 0;
        const Point &p = C->getPoint( i );
        
        for ( int c = 0; c < 3; c++ ) {
            v1 = p[c] - centroid[c];
            
            for ( int r = 0; r < 3; r++ ) {
                v2 = p[r] - centroid[r];
                
                this->varMatrix[idx] += v1 * v2;
                
                ++idx;
            }
        }
    }
    eslog::error("call upDense3x3EigenValuesEigenVectors\n");
//    MATH::upDense3x3EigenValuesEigenVectors(&this->varMatrix[0], &this->eigenvalues[0], &this->eigenvectors[0]);
    
    // find the largest eigenvalue
    double maxEig = 0;
    idx = 0;
    for ( int i = 0; i < 3; i++ ) {
        if ( fabs( maxEig ) < fabs( this->eigenvalues[i] ) ) {
            maxEig = this->eigenvalues[i];
            idx = i;
        }
    }

    //normal to the cutting plane
    this->normal[0] = this->eigenvectors[3 * idx];
    this->normal[1] = this->eigenvectors[3 * idx + 1];
    this->normal[2] = this->eigenvectors[3 * idx + 2];

    // split elements into two groups by the plane defined by the normal vector
    std::vector<esint> *leftPointIndices = new std::vector<esint>( );
    std::vector<esint> *rightPointIndices = new std::vector<esint>( );
    
    double w;
    for ( esint i = 0; i < n; i++ ) {
        const Point &p = C->getPoint( i );
        v1 = (p[0] - centroid[0]) * this->normal[0];
        v2 = (p[1] - centroid[1]) * this->normal[1];
        v3 = (p[2] - centroid[2]) * this->normal[2];
        
        w = v1 + v2 + v3;

        if ( w >= 0 ) {
            leftPointIndices->push_back( C->getPointIndexGlobal( i ) );
        } else {
            rightPointIndices->push_back( C->getPointIndexGlobal( i ) );
        }
    }

    // we subdivide this cluster only when it has both non-empty children
    if(leftPointIndices->size() > 0 && rightPointIndices->size() > 0){
        C->createLeftChild( leftPointIndices );
        C->createRightChild( rightPointIndices );

        this->doCreateClusterTree(C->getLeftChild(), s);
        this->doCreateClusterTree(C->getRightChild(), s);
    }
    else{
        delete leftPointIndices;
        delete rightPointIndices;
    }
}

