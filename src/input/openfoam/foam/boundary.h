/*
 * Boundary.h
 *
 *  Created on: Dec 4, 2016
 *      Author: beh01
 */

#ifndef SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_
#define SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_

#include <vector>
#include <ostream>
#include "../../loader.h"

namespace espreso {
namespace input {

class Boundary {
public:
	Boundary(int procNo);
	virtual ~Boundary();

	int getProcNo() const {
		return this->procNo;
	}

	void add(eslocal node) {
		nodes.push_back(node);
	}

	std::vector< eslocal >& getNodes() {
		return nodes;
	}

	void prepareNodes();

	friend inline std::ostream& operator<<(std::ostream& os, const Boundary& obj)
	{
		os << "Rank: " << obj.procNo <<"{\n";
		for (auto node : obj.nodes) {
			os << node << "\n";
		}
		os << "};\n";
		return os;
	}

private:
	int procNo;
	std::vector< eslocal > nodes;

};
}
}
#endif /* SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_ */
