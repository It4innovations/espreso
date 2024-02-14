
#ifndef SRC_CONFIGURATION_MESHMORPHING_H_
#define SRC_CONFIGURATION_MESHMORPHING_H_

#include "config/holders/expression.h"
#include "config/description.h"
#include "config/ecf/material/coordinatesystem.h"

namespace espreso {

struct ECF;

enum class MORPHING_TYPE {
	NONE = 0,
	RBF = 1
};

enum class MORPHING_TRANSFORMATION {
	FIXED,
	OFFSET,
	SCALING,
	TRANSLATION,
	ROTATION
};

struct RBFTargetTransformationConfiguration: public ECFDescription {

	MORPHING_TRANSFORMATION transformation;

	ECFExpression offset;
	ECFExpressionVector scaling, translation;
	CoordinateSystemConfiguration coordinate_system;

	bool override;

	RBFTargetTransformationConfiguration(ECF *ECF);
protected:
	ECF *_ECF;
};

struct ExternalFFDConfiguration: public ECFDescription {

	std::string path;
	std::map<std::string, RBFTargetTransformationConfiguration> morphers;

	ExternalFFDConfiguration(ECF *ECF);
};

struct RBFTargetConfiguration: public ECFDescription {

	ECFExpression function;
	double solver_precision;
	int solver_max_iter;
	int polynomial_regularization_degree;
	bool use_transform_translate;
	bool use_transform_scale;
	bool use_transform_rotate;

	double aca_precision;
	double aca_eta;
	int aca_cluster_tree_leaf_size;

	std::string target;
	std::map<std::string, RBFTargetTransformationConfiguration> morphers;

	ExternalFFDConfiguration external_ffd;

	RBFTargetConfiguration(ECF *ECF);
};

struct MeshMorphing: public ECFDescription {

	MORPHING_TYPE type;
	std::map<std::string, RBFTargetConfiguration> rbf;

	MeshMorphing(ECF *ECF);
};

}

#endif /* SRC_CONFIGURATION_MESHMORPHING_H_ */
