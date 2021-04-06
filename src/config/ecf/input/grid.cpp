
#include "grid.h"

#include "config/configuration.hpp"

espreso::GridGeneratorConfiguration::GridGeneratorConfiguration()
{
	chessboard_size = 2;

	blocks_x = blocks_y = blocks_z = 1;
	clusters_x = clusters_y = clusters_z = 1;

	ecfdescription->addSeparator();

	REGISTER(noncontinuous, ECFMetaData()
			.setdescription({ "A cluster ID.", "A number of pseudo non-continuous parts." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::POSITIVE_INTEGER })
			.setpattern({ "0", "2" }));

	nonuniform_nparts = 0;
	REGISTER(nonuniform_nparts, ECFMetaData()
			.setdescription({ "A number of parts when nonuniform decomposition is set (0 for keep the number of parst the same)." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

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

	ecfdescription->addSpace();

	REGISTER(chessboard_size, ECFMetaData()
			.setdescription({ "A number of chess-board blocks." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSpace();

	REGISTER(blocks_x, ECFMetaData()
			.setdescription({ "Number of blocks in x-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(blocks_y, ECFMetaData()
			.setdescription({ "Number of blocks in y-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(blocks_z, ECFMetaData()
			.setdescription({ "Number of blocks in z-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSpace();

	REGISTER(clusters_x, ECFMetaData()
			.setdescription({ "Number of clusters in x-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(clusters_y, ECFMetaData()
			.setdescription({ "Number of clusters in y-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(clusters_z, ECFMetaData()
			.setdescription({ "Number of clusters in z-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSpace();

	REGISTER(blocks, ECFMetaData()
			.setdescription({ "A block ID.", "Turn generation of a block on/off." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::BOOL })
			.setpattern({ "0", "FALSE" }));
}



