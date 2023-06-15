
#ifndef SRC_ANALYSIS_COMPOSER_MAPPING_H_
#define SRC_ANALYSIS_COMPOSER_MAPPING_H_

#include <cstddef>
#include <vector>

namespace espreso {

template <typename T>
struct Mapping {

	struct Map {
		int filter;
		T* data;
		const esint* position;

		Map(): filter(255), data(nullptr), position(nullptr) {}
	};

	std::vector<Map> elements; // per interval map
	std::vector<std::vector<Map> > boundary; // per boundary x interval
};

}

#endif /* SRC_ANALYSIS_COMPOSER_MAPPING_H_ */
