
#include "config/configuration.hpp"
#include "config/conditions.h"
#include "thermalconductivity.h"

espreso::ThermalConductivityConfiguration::ThermalConductivityConfiguration(DIMENSION *D)
: model(MODEL::ISOTROPIC),
  dimension(D),
  values(3)
{
	REGISTER(model, ECFMetaData()
			.setdescription({ "Model" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("DIAGONAL").setdescription("Orthotropic model."))
			.addoption(ECFOption().setname("SYMMETRIC").setdescription("Symmetric model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return *dimension == DIMENSION::D3; })));

	ecfdescription->registerParameter("KXX", values.get(0, 0), ECFMetaData()
			.setdescription({ "KXX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 0, 0));
	ecfdescription->registerParameter("KYY", values.get(1, 1), ECFMetaData()
			.setdescription({ "KYY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 1, 1)
			.allowonly([&] () { return model != MODEL::ISOTROPIC; })
			.addconstraint(ECFCondition(model, ECFCondition::NOT_EQUALS, MODEL::ISOTROPIC)));
	ecfdescription->registerParameter("KZZ", values.get(2, 2), ECFMetaData()
			.addconstraint(ECFCondition(model, ECFCondition::NOT_EQUALS, MODEL::ISOTROPIC))
			.setdescription({ "KZZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 2, 2)
			.allowonly([&] () { return model != MODEL::ISOTROPIC && *dimension == DIMENSION::D3; })
			.addconstraint(ECFCondition(model, ECFCondition::NOT_EQUALS, MODEL::ISOTROPIC) & ECFCondition(*dimension, ECFCondition::EQUALS, DIMENSION::D3)));

	ecfdescription->registerParameter("KXY", values.get(0, 1), ECFMetaData()
			.setdescription({ "KXY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 0, 1)
			.allowonly([&] () { return model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC; })
			.addconstraint(ECFCondition(model, ECFCondition::EQUALS, MODEL::SYMMETRIC) | ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC)));
	ecfdescription->registerParameter("KXZ", values.get(0, 2), ECFMetaData()
			.setdescription({ "KXZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 0, 2)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC) && *dimension == DIMENSION::D3; })
			.addconstraint((ECFCondition(model, ECFCondition::EQUALS, MODEL::SYMMETRIC) | ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC)) & ECFCondition(*dimension, ECFCondition::EQUALS, DIMENSION::D3)));
	ecfdescription->registerParameter("KYZ", values.get(1, 2), ECFMetaData()
			.setdescription({ "KYZ" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 1, 2)
			.allowonly([&] () { return (model == MODEL::SYMMETRIC || model == MODEL::ANISOTROPIC)  && *dimension == DIMENSION::D3; })
			.addconstraint((ECFCondition(model, ECFCondition::EQUALS, MODEL::SYMMETRIC) | ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC)) & ECFCondition(*dimension, ECFCondition::EQUALS, DIMENSION::D3)));

	ecfdescription->registerParameter("KYX", values.get(1, 0), ECFMetaData()
			.setdescription({ "KYX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 1, 0)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; })
			.addconstraint(ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC)));
	ecfdescription->registerParameter("KZX", values.get(2, 0), ECFMetaData()
			.setdescription({ "KZX" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 2, 0)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && *dimension == DIMENSION::D3; })
			.addconstraint(ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC) & ECFCondition(*dimension, ECFCondition::EQUALS, DIMENSION::D3)));
	ecfdescription->registerParameter("KZY", values.get(2, 1), ECFMetaData()
			.setdescription({ "KZY" })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(values, 2, 1)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC && *dimension == DIMENSION::D3; })
			.addconstraint(ECFCondition(model, ECFCondition::EQUALS, MODEL::ANISOTROPIC) & ECFCondition(*dimension, ECFCondition::EQUALS, DIMENSION::D3)));
}
