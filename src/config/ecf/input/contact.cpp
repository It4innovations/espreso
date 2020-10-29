
#include "contact.h"
#include "config/configuration.hpp"

using namespace espreso;

ContactConfiguration::ContactConfiguration()
{
	search_area = 0.01;
	REGISTER(search_area, ECFMetaData()
			.setdescription({ "Area for contact search algorithm." })
			.setdatatype({ ECFDataType::FLOAT }));

	max_angle = 80;
	REGISTER(max_angle, ECFMetaData()
			.setdescription({ "Maximal angle between faces." })
			.setdatatype({ ECFDataType::FLOAT }));

	self_contact = false;
	REGISTER(self_contact, ECFMetaData()
			.setdescription({ "Include contact within the same bodies." })
			.setdatatype({ ECFDataType::BOOL }));
}
