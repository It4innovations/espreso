
#ifndef SRC_PHYSICS_KERNELS_UTILS_MORPHING_H_
#define SRC_PHYSICS_KERNELS_UTILS_MORPHING_H_

#include "basis/containers/point.h"

#include <string>
#include <vector>
#include <map>

#include "math/math.h"

#include "morphing/morphing_system.h"


namespace espreso {

struct RBFTargetConfiguration;
struct RBFTargetTransformationConfiguration;

namespace morphing {
	

void morphRBF(
	const std::string &name, 
	const RBFTargetConfiguration &configuration, 
	int dimension
);

void processMorpher(
	const RBFTargetTransformationConfiguration &target, 
	int dimension, std::vector<Point> &sPoints, 
	size_t startPoint, 
	std::vector<double> &sDisplacement
);

esint prepareMatrixM(
	std::vector<Point> &rPoints, 
	std::vector<double> &rDisplacement, 
	int dimension, 
	const RBFTargetConfiguration &configuration, 
	std::vector<double> &M_values, 
	bool use_x = true,
	bool use_y = true,
	bool use_z = true
);

void readExternalFile(
	const RBFTargetConfiguration &configuration, 
	int dimension, 
	std::map<std::string, 
	std::vector<Point>> &external_data
);

void rbf();

}

}

#endif /* SRC_PHYSICS_KERNELS_UTILS_MORPHING_H_ */
