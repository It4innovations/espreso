
#ifndef SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_
#define SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_

#include "config/description.h"

namespace espreso {

struct ContactInterfaceConfiguration: public ECFDescription {

	double search_area;
	double max_angle;
	bool self_contact;

	ContactInterfaceConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_ */
