
#ifndef SRC_CONFIG_ECF_INPUT_CONTACT_H_
#define SRC_CONFIG_ECF_INPUT_CONTACT_H_

#include "config/description.h"

namespace espreso {

struct ContactConfiguration: public ECFDescription {

	double search_area;
	double max_angle;
	bool self_contact;

	ContactConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_CONTACT_H_ */
