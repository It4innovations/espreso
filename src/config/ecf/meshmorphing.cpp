
#include <config/ecf/ecf.h>
#include "config/configuration.hpp"
#include "esinfo/eslog.hpp"


espreso::RBFTargetTransformationConfiguration::RBFTargetTransformationConfiguration(ECF *ECF)
: dimension(DIMENSION::D3),
  offset(ECFMetaData::getboundaryconditionvariables()),
  scaling(&dimension, ECFMetaData::getboundaryconditionvariables(), "1"),
  translation(&dimension, ECFMetaData::getboundaryconditionvariables(), "0"),
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
			eslog::globalerror("ESPRESO internal error: unknown physics while set RBFTargetTransformation.");
		}
	});
}

espreso::RBFTargetConfiguration::RBFTargetConfiguration(ECF *ECF)
: function({ "R" }, "R"),
  external_ffd(ECF)
{
	solver = MORPHING_RBF_SOLVER::DIRECT;
	REGISTER(solver, ECFMetaData()
		.setdescription({ "Mesh Morphing solver." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("ITERATIVE").setdescription("Iterative"))
		.addoption(ECFOption().setname("DIRECT").setdescription("Direct")));

	solver_precision = 1e-5;
	REGISTER(solver_precision, ECFMetaData()
		.setdescription({ "Solver requested precision." })
		.setdatatype({ ECFDataType::FLOAT }));

	solver_max_iter = 600;
		REGISTER(solver_max_iter, ECFMetaData()
			.setdescription({ "Solver requested maximum number of iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(function, ECFMetaData()
		.setdescription({ "Radial basis function." })
		.setdatatype({ ECFDataType::EXPRESSION }));

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


