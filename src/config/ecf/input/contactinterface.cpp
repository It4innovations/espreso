
#include "contactinterface.h"
#include "config/configuration.hpp"

using namespace espreso;

ContactInterfaceConfiguration::ContactInterfaceConfiguration()
{
	detection = DETECTION::ALL_BODIES;
	REGISTER(detection, ECFMetaData()
			.setdescription({ "A set for contact detection." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ALL_BODIES").setdescription("Minimum."))
			.addoption(ECFOption().setname("BODY_LIST").setdescription("Maximum."))
//			.addoption(ECFOption().setname("CONTACT_PAIR").setdescription("Average."))
			);

	criterion = CRITERION::BOUND;
	REGISTER(criterion, ECFMetaData()
			.setdescription({ "Criteria what is considered as contact." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("BOUND").setdescription("Only bound faces are considered as contact."))
			.addoption(ECFOption().setname("GAP").setdescription("All faces within a given gap are considered as contact."))
			.addoption(ECFOption().setname("SKIP").setdescription("Skip this contact interface.")));

	REGISTER(body_list, ECFMetaData()
			.setdescription({ "List of element regions considered for contact." })
			.setdatatype({ ECFDataType::STRING }));

	self_contact = false;
	REGISTER(self_contact, ECFMetaData()
			.setdescription({ "Include contact within the same bodies." })
			.setdatatype({ ECFDataType::BOOL }));

	gap = 0.001;
	REGISTER(gap, ECFMetaData()
			.setdescription({ "Area for contact search algorithm." })
			.setdatatype({ ECFDataType::FLOAT }));

	angle = 30;
	REGISTER(angle, ECFMetaData()
			.setdescription({ "Maximal angle between faces." })
			.setdatatype({ ECFDataType::FLOAT }));
}
