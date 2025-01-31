
#ifndef SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_
#define SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_

#include "config/description.h"

#include <string>
#include <vector>

namespace espreso {

struct ContactInterfaceConfiguration: public ECFDescription {

	enum class DETECTION {
		ALL_BODIES,
		BODY_LIST,
		CONTACT_PAIR
	};

	enum class CRITERION {
		BOUND,
		GAP,
		SKIP
	};

	DETECTION detection;
	CRITERION criterion;

	std::vector<std::string> body_list;

	bool self_contact;

	float gap;
	float angle;

	std::vector<esint> found_interfaces;

	ContactInterfaceConfiguration();
};

}



#endif /* SRC_CONFIG_ECF_INPUT_CONTACTINTERFACE_H_ */
