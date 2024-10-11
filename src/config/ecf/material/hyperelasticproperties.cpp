
#include "hyperelasticproperties.h"
#include "config/configuration.hpp"

espreso::HyperElasticPropertiesConfiguration::HyperElasticPropertiesConfiguration()
: model(MODEL::NEO_HOOKEN_CMP)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Material model." })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("KIRCHHOFF").setdescription("Kirchhoff linear"))
			.addoption(ECFOption().setname("NEO_HOOKEN_CMP").setdescription("Neo-Hooken-compressible"))
			.addoption(ECFOption().setname("NEO_HOOKEN_INC").setdescription("Neo-Hooken-incompresible"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_2PARAMS").setdescription("Mooney-Rivlin with 2 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_3PARAMS").setdescription("Mooney-Rivlin with 3 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_5PARAMS").setdescription("Mooney-Rivlin with 5 parameters"))
			.addoption(ECFOption().setname("MOONEY_RIVLIN_9PARAMS").setdescription("Mooney-Rivlin with 9 parameters"))
			.addoption(ECFOption().setname("ARRUDA_BOYCE").setdescription("Arruda-Boyce"))
			.addoption(ECFOption().setname("BLATZ_KO_FOAM").setdescription("Blatz-ko FOAM"))
			.addoption(ECFOption().setname("GENT").setdescription("Gent"))
			.addoption(ECFOption().setname("OGDEN_1").setdescription("Ogden 1"))
			.addoption(ECFOption().setname("OGDEN_2").setdescription("Ogden 2"))
			.addoption(ECFOption().setname("OGDEN_3").setdescription("Ogden 3")));

	ecfdescription->addSeparator();

	REGISTER(C01, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () {
				return
						model == MODEL::MOONEY_RIVLIN_2PARAMS ||
						model == MODEL::MOONEY_RIVLIN_3PARAMS ||
						model == MODEL::MOONEY_RIVLIN_5PARAMS ||
						model == MODEL::MOONEY_RIVLIN_9PARAMS;;
			}));

	REGISTER(C10, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () {
				return
						model == MODEL::MOONEY_RIVLIN_2PARAMS ||
						model == MODEL::MOONEY_RIVLIN_3PARAMS ||
						model == MODEL::MOONEY_RIVLIN_5PARAMS ||
						model == MODEL::MOONEY_RIVLIN_9PARAMS;
			}));

	REGISTER(C11, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () {
				return
						model == MODEL::MOONEY_RIVLIN_3PARAMS ||
						model == MODEL::MOONEY_RIVLIN_5PARAMS ||
						model == MODEL::MOONEY_RIVLIN_9PARAMS;;
			}));

	REGISTER(C02, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () {
				return
					model == MODEL::MOONEY_RIVLIN_5PARAMS ||
					model == MODEL::MOONEY_RIVLIN_9PARAMS;;
			}));

	REGISTER(C20, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () {
				return
					model == MODEL::MOONEY_RIVLIN_5PARAMS ||
					model == MODEL::MOONEY_RIVLIN_9PARAMS;;
			}));

	REGISTER(C30, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return model == MODEL::MOONEY_RIVLIN_9PARAMS; }));

	REGISTER(C21, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return model == MODEL::MOONEY_RIVLIN_9PARAMS; }));

	REGISTER(C12, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return model == MODEL::MOONEY_RIVLIN_9PARAMS; }));

	REGISTER(C03, ECFMetaData()
			.setdescription({ "Material constant." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.allowonly([&] () { return model == MODEL::MOONEY_RIVLIN_9PARAMS; }));

	REGISTER(E, ECFMetaData()
			.setdescription({ "E parameter." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(mu, ECFMetaData()
			.setdescription({ "mu parameter." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(d, ECFMetaData()
			.setdescription({ "Incompressibility parameter." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(G, ECFMetaData()
			.setdescription({ "Initial share modulus." })
			.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(lambda, ECFMetaData()
			.setdescription({ "Limiting network stretch." })
			.setdatatype({ ECFDataType::EXPRESSION }));
}

