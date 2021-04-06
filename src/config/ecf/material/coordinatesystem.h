
#ifndef SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_
#define SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_

#include "basis/containers/point.h"

#include "config/holders/expression.h"
#include "config/description.h"

namespace espreso {

struct CoordinateSystemConfiguration: public ECFDescription {

	enum class TYPE {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	TYPE type;
	DIMENSION *dimension;

	ECFExpressionVector rotation;
	ECFExpressionVector center;

	CoordinateSystemConfiguration(DIMENSION *D);

	void createScalingMatrix(std::vector<double> &m, double x, double y, double z=0) const;

	void createTranslationMatrixToCenter(std::vector<double> &m) const;
	void createTranslationMatrixToZero(std::vector<double> &m) const;
	void createTranslationMatrix(std::vector<double> &m, double x, double y, double z=0) const;
	void createRotationMatrix(std::vector<double> &m) const;

	void multiplyTransformationMatrices(std::vector<double> &left, std::vector<double> &result) const;
	Point applyTransformation(std::vector<double> &m, const Point &p) const;

};


}

#endif /* SRC_CONFIG_ECF_MATERIAL_COORDINATESYSTEM_H_ */
