#include "sphereproblem.h"

#include <cmath>

#include "../config/configuration.hpp"

using namespace espreso;

SphereProblem::SphereProblem()
{
    x, y = 0;

    REGISTER(x, ECFMetaData()
            .setname("X")
			.setdescription({ "X" })
			.setdatatype({ ECFDataType::INTEGER }));

	REGISTER(y, ECFMetaData()
			.setname("Y")
			.setdescription({ "Heat capacity" })
			.setdatatype({ ECFDataType::INTEGER }));
	
	// x, y = SphereEnum::FIVE;

	// REGISTER(x, ECFMetaData()
    //         .setname("X")
	// 		.setdescription({ "X" })
	// 		.setdatatype({ ECFDataType::OPTION })
	// 		.addoption(ECFOption().setname("ZERO").setdescription("ZERO"))
	// 		.addoption(ECFOption().setname("ONE").setdescription("ONE"))
	// 		.addoption(ECFOption().setname("TWO").setdescription("TWO"))
	// 		.addoption(ECFOption().setname("THREE").setdescription("THREE"))
	// 		.addoption(ECFOption().setname("FOUR").setdescription("FOUR"))
	// 		.addoption(ECFOption().setname("FIVE").setdescription("FIVE"))
	// );

	// REGISTER(y, ECFMetaData()
	// 		.setname("Y")
	// 		.setdescription({ "Heat capacity" })
	// 		.setdatatype({ ECFDataType::OPTION })
	// 		.addoption(ECFOption().setname("ZERO").setdescription("ZERO"))
	// 		.addoption(ECFOption().setname("ONE").setdescription("ONE"))
	// 		.addoption(ECFOption().setname("TWO").setdescription("TWO"))
	// 		.addoption(ECFOption().setname("THREE").setdescription("THREE"))
	// 		.addoption(ECFOption().setname("FOUR").setdescription("FOUR"))
	// 		.addoption(ECFOption().setname("FIVE").setdescription("FIVE"))
	// );
}

double SphereProblem::evaluate()
{
	// return pow(static_cast<int>(x), 2) + pow(static_cast<int>(y), 2);
	return pow(x, 2) + pow(y, 2);
}