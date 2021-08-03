
#include <config/ecf/ecf.h>
#include "config/configuration.hpp"
#include "esinfo/eslog.hpp"


espreso::RBFTargetTransformationConfiguration::RBFTargetTransformationConfiguration(ECF *ECF)
: dimension(DIMENSION::D3),
  offset(ECFExpression::Scope::NODE),
  scaling(&dimension, "1", ECFExpression::Scope::NODE),
  translation(&dimension, "0", ECFExpression::Scope::NODE),
  coordinate_system(&dimension),
  override(true),
  _ECF(ECF)
{
	transformation = MORPHING_TRANSFORMATION::TRANSLATION;
	REGISTER(transformation, ECFMetaData()
		.setdescription({ "Transformation variant." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("FIXED").setdescription("Fixed position."))
		.addoption(ECFOption().setname("OFFSET").setdescription("OFFSET."))
		.addoption(ECFOption().setname("SCALING").setdescription("Scaling."))
		.addoption(ECFOption().setname("TRANSLATION").setdescription("Translation."))
		.addoption(ECFOption().setname("ROTATION").setdescription("ROTATION."))
		);

	ecfdescription->addSeparator();

	REGISTER(offset, ECFMetaData()
		.setdescription({ "Offset size." })
		.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(scaling, ECFMetaData()
		.setdescription({ "Scale vector." }));

	REGISTER(translation, ECFMetaData()
		.setdescription({ "Translation vector." }));

	ecfdescription->addSeparator();

	REGISTER(coordinate_system, ECFMetaData()
		.setdescription({ "Configuration of coordinate system." }));

	REGISTER(override, ECFMetaData()
		.setdescription({ "Turn morphing target override on/off." })
		.setdatatype({ ECFDataType::BOOL }));

	ECF->ecfdescription->getParameter(&ECF->physics)->addListener(ECFParameter::Event::VALUE_SET, [&] (const std::string &value) {
		switch (ECF->physics) {
		case PhysicsConfiguration::TYPE::HEAT_TRANSFER_2D:
		case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_2D:
			*translation.dimension = DIMENSION::D2;
			break;
		case PhysicsConfiguration::TYPE::HEAT_TRANSFER_3D:
		case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS_3D:
			*translation.dimension = DIMENSION::D3;
			break;
		default:
			eslog::internalFailure("unknown physics while set RBFTargetTransformation.");
		}
	});
}

espreso::RBFTargetConfiguration::RBFTargetConfiguration(ECF *ECF)
: function("R", ECFExpression::Scope::GLOBAL),
  external_ffd(ECF)
{

	solver_precision = 1e-5;
	REGISTER(solver_precision, ECFMetaData()
		.setdescription({ "Solver requested precision." })
		.setdatatype({ ECFDataType::FLOAT }));
		
	solver_max_iter = 600;
		REGISTER(solver_max_iter, ECFMetaData()
			.setdescription({ "Solver requested maximum number of iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	polynomial_regularization_degree = 0;
		REGISTER(polynomial_regularization_degree, ECFMetaData()
			.setdescription({ "Morphing degree of the regularization polynomial." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	use_transform_translate = false;
		REGISTER(use_transform_translate, ECFMetaData()
			.setdescription({ "Morphing auto-translate referential points to the point-cloud centroid." })
			.setdatatype({ ECFDataType::BOOL }));

	use_transform_scale = false;
		REGISTER(use_transform_scale, ECFMetaData()
			.setdescription({ "Morphing auto-rescale referential points and displacements." })
			.setdatatype({ ECFDataType::BOOL }));

	use_transform_rotate = false;
		REGISTER(use_transform_rotate, ECFMetaData()
			.setdescription({ "Morphing auto-rotate referential points and displacements." })
			.setdatatype({ ECFDataType::BOOL }));

	REGISTER(function, ECFMetaData()
		.setdescription({ "Radial basis function." })
		.setdatatype({ ECFDataType::EXPRESSION }));

	aca_precision = 1e-6;
	REGISTER(aca_precision, ECFMetaData()
		.setdescription({ "ACA precision." })
		.setdatatype({ ECFDataType::FLOAT }));

	aca_eta = 2;
	REGISTER(aca_eta, ECFMetaData()
		.setdescription({ "ACA eta." })
		.setdatatype({ ECFDataType::FLOAT }));

	aca_cluster_tree_leaf_size = 10;
		REGISTER(aca_cluster_tree_leaf_size, ECFMetaData()
			.setdescription({ "Base size of cluster tree leaves in ACA." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSeparator();

	REGISTER(target, ECFMetaData()
		.setdescription({ "Set of morphed elements." })
		.setdatatype({ ECFDataType::BOUNDARY_REGION }));

	REGISTER(morphers , ECFMetaData()
		.setdescription({ "Morphed region name.", "Target configuration." })
		.setdatatype({ ECFDataType::BOUNDARY_REGION })
		.setpattern({ "REGION" }),
		ECF);

	REGISTER(external_ffd, ECFMetaData()
		.setdescription({ "Configurations defined by an external FFD file." }));
}


espreso::ExternalFFDConfiguration::ExternalFFDConfiguration(ECF *ECF)
{
	path = "";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Set path to external configuration file with regions." })
			.setdatatype({ ECFDataType::STRING }));

	REGISTER(morphers , ECFMetaData()
			.setdescription({ "Morphed region name.", "Target configuration." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "REGION" }),
			ECF);
}

espreso::MeshMorphing::MeshMorphing(ECF *ECF)
{
	type = MORPHING_TYPE::NONE;
	REGISTER(type, ECFMetaData()
		.setdescription({ "Mesh Morphing type." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("NONE").setdescription("No Mesh Morphing."))
		.addoption(ECFOption().setname("RBF").setdescription("RBF (Radial Base Function) Mesh Morphing.")));

	REGISTER(rbf, ECFMetaData()
		.setdescription({ "Morphing name.", "Named RBF configuration." })
		.setdatatype({ ECFDataType::STRING })
		.setpattern({ "MORPHING_NAME" }),
		ECF);
}


