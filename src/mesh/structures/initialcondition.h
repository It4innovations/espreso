
#ifndef SRC_MESH_STRUCTURES_INITIALCONDITION_H_
#define SRC_MESH_STRUCTURES_INITIALCONDITION_H_

#include "../elements/element.h"

namespace espreso {

class InitialCondition {

public:
	virtual InitialCondition* copy() =0;

	virtual double get(const Point &p) =0;
	virtual ~InitialCondition() {};
};

//class IntervalInitialization: public InitialCondition {
//
//public:
//	std::vector<eslocal> elements;
//
//
//	IntervalInitialization* copy()
//	{
//		return new ElementInitialization(*this);
//	}
//};
}



#endif /* SRC_MESH_STRUCTURES_INITIALCONDITION_H_ */
