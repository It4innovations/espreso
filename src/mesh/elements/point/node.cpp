
#include "node.h"

using namespace espreso;

std::vector<Property> Node::_DOFElement;
std::vector<Property> Node::_DOFFace;
std::vector<Property> Node::_DOFEdge;
std::vector<Property> Node::_DOFPoint;
std::vector<Property> Node::_DOFMidPoint;

std::vector<DenseMatrix> Node::_dN;
std::vector<DenseMatrix> Node::_N;
std::vector<double> Node::_weighFactor;


