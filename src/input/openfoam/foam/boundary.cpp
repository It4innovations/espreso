/*
 * Boundary.cpp
 *
 *  Created on: Dec 4, 2016
 *      Author: beh01
 */

#include "boundary.h"

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

