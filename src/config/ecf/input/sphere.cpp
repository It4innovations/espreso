
#include "sphere.h"

#include "config/configuration.hpp"

espreso::SphereGeneratorConfiguration::SphereGeneratorConfiguration()
{
	inner_radius = 5;
	outer_radius = 10;
	clusters = 1;
	layers = 1;

	REGISTER(inner_radius, ECFMetaData()
				.setdescription({ "Inner radius of generater sphere." })
				.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(outer_radius, ECFMetaData()
				.setdescription({ "Outer radius of generater sphere." })
				.setdatatype({ ECFDataType::FLOAT }));

	ecfdescription->addSeparator();

	REGISTER(nodes, ECFMetaData()
			.setdescription({ "The name of generated region.", "A specification of a region." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::INTERVAL })
			.setpattern({ "MY_NODE_REGION", "(0, 1) (0, 1) <0, 0>" }));

	REGISTER(edges, ECFMetaData()
			.setdescription({ "The name of generated region.", "A specification of a region." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::INTERVAL })
			.setpattern({ "MY_EDGE_REGION", "<0, 1> <0, 0> <0, 0>" }));

	REGISTER(faces, ECFMetaData()
			.setdescription({ "The name of generated region.", "A specification of a region." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::INTERVAL })
			.setpattern({ "MY_FACE_REGION", "<0, 1> <0, 1> <0, 0>" }));

	REGISTER(elements, ECFMetaData()
			.setdescription({ "The name of generated region.", "A specification of a region." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::INTERVAL })
			.setpattern({ "MY_ELEMENT_REGION", "CHESSBOARD_WHITE" })
			.addoption(ECFOption().setname("CHESSBOARD_WHITE").setdescription("White parts of chess-board."))
			.addoption(ECFOption().setname("CHESSBOARD_BLACK").setdescription("Black parts of chess-board.")));

	ecfdescription->addSeparator();

	REGISTER(clusters, ECFMetaData()
			.setdescription({ "Number of clusters in x, y-direction in generated SPHERE." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(layers, ECFMetaData()
			.setdescription({ "Number of layers in generated SPHERE." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

}

