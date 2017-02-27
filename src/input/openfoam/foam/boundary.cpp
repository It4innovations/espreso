
#include "boundary.h"

#include <algorithm>

#include "../../../basis/utilities/utils.h"

using namespace espreso::input;

Boundary::Boundary(int procNo) {
	this->procNo = procNo;
}

Boundary::~Boundary() {
	// TODO Auto-generated destructor stub
}

void Boundary::prepareNodes() {
	std::sort(nodes.begin(), nodes.end());
	Esutils::removeDuplicity(nodes);
}

