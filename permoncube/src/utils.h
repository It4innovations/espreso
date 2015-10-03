
#ifndef PM_UTILS_H_
#define PM_UTILS_H_

#include "settings.h"

namespace esinput {

template <class TElement>
class Utils {

public:
	static void globalNodesCount(const Settings &settings, esglobal nodes[]);
	static void clusterNodesCount(const Settings &settings, eslocal nodes[]);
	static eslocal clusterElementsCount(const Settings &settings);
};

}


#include "utils.hpp"


#endif /* PM_UTILS_H_ */
