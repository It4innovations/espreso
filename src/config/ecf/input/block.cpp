
#include "block.h"

#include "config/configuration.hpp"

espreso::BlockGeneratorConfiguration::BlockGeneratorConfiguration()
{
	element_type = GENERATOR_ELEMENT_TYPE::HEXA8;
	start_x.value = start_y.value = start_z.value = "0";
	length_x.value = length_y.value = length_z.value = "1";
	projection_x.value = "x";
	projection_y.value = "y";
	projection_z.value = "z";

	domains_x = domains_y = domains_z = 2;
	elements_x = elements_y = elements_z = 5;

	ecfdescription->addSeparator();

	REGISTER(start_x, ECFMetaData()
			.setdescription({ "A x-coordinate of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(start_y, ECFMetaData()
			.setdescription({ "A y-coordinate of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(start_z, ECFMetaData()
			.setdescription({ "A z-coordinate of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	ecfdescription->addSpace();

	REGISTER(length_x, ECFMetaData()
			.setdescription({ "A x-length of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(length_y, ECFMetaData()
			.setdescription({ "A y-length of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(length_z, ECFMetaData()
			.setdescription({ "A z-length of generated GRID." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	ecfdescription->addSpace();

	REGISTER(projection_x, ECFMetaData()
			.setdescription({ "Projection of x-coordinate of generated GRID." })
			.setvariables({ "X" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(projection_y, ECFMetaData()
			.setdescription({ "Projection of y-coordinate of generated GRID." })
			.setvariables({ "Y" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(projection_z, ECFMetaData()
			.setdescription({ "Projection of z-coordinate of generated GRID." })
			.setvariables({ "Z" })
			.setdatatype({ ECFDataType::EXPRESSION }));

	ecfdescription->addSeparator();

	REGISTER(element_type, ECFMetaData()
			.setdescription({ "A type of an element generated by GRID GENERATOR." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("HEXA8").setdescription("HEXAHEDRON"))
			.addoption(ECFOption().setname("HEXA20").setdescription("HEXAHEDRON with midpoints"))
			.addoption(ECFOption().setname("TETRA4").setdescription("TETRAHEDRON"))
			.addoption(ECFOption().setname("TETRA10").setdescription("TETRAHEDRON with midpoints"))
			.addoption(ECFOption().setname("PRISMA6").setdescription("PRISMA"))
			.addoption(ECFOption().setname("PRISMA15").setdescription("PRISMA with midpoints"))
			.addoption(ECFOption().setname("PYRAMID5").setdescription("PYRAMID"))
			.addoption(ECFOption().setname("PYRAMID13").setdescription("PYRAMID with midpoints"))
			.addoption(ECFOption().setname("SQUARE4").setdescription("SQUARE"))
			.addoption(ECFOption().setname("SQUARE8").setdescription("SQUARE with midpoints"))
			.addoption(ECFOption().setname("TRIANGLE3").setdescription("TRIANGLE"))
			.addoption(ECFOption().setname("TRIANGLE6").setdescription("TRIANGLE with midpoints")));

	ecfdescription->addSpace();

	REGISTER(domains_x, ECFMetaData()
			.setdescription({ "Number of domains in x-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(domains_y, ECFMetaData()
			.setdescription({ "Number of domains in y-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(domains_z, ECFMetaData()
			.setdescription({ "Number of domains in z-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	ecfdescription->addSpace();

	REGISTER(elements_x, ECFMetaData()
			.setdescription({ "Number of elements in x-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(elements_y, ECFMetaData()
			.setdescription({ "Number of elements in y-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(elements_z, ECFMetaData()
			.setdescription({ "Number of elements in z-direction in generated GRID." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));
}

