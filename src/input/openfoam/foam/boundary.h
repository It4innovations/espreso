/*
 * Boundary.h
 *
 *  Created on: Dec 4, 2016
 *      Author: beh01
 */

#ifndef SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_
#define SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_

#include <set>
#include "../../loader.h"

class Boundary {
public:
	Boundary(int procNo);
	virtual ~Boundary();

	int getProcNo() const {
		return this->procNo;
	}

	void add(eslocal node) {
		nodes.insert(node);
	}

	std::set< eslocal >& getNodes() {
		return nodes;
	}


	friend inline std::ostream& operator<<(std::ostream& os,
				const Boundary& obj) {
			os << "Rank: " << obj.procNo <<"{\n";
			for (auto node : obj.nodes) {
				os << node << "\n";
			}
			os << "};\n";
			return os;
	}

private:
	int procNo;
	std::set< eslocal > nodes;

};

#endif /* SRC_INPUT_OPENFOAM_FOAM_BOUNDARY_H_ */
