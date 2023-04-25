
#include "linearelasticproperties.h"
#include "config/configuration.hpp"
#include "config/conditions.h"

espreso::LinearElasticPropertiesConfiguration::LinearElasticPropertiesConfiguration(DIMENSION *D)
: orientation(false),
  model(MODEL::ISOTROPIC),
  dimension(D),
  poisson_ratio(3),
  young_modulus(3),
  shear_modulus(3),
  anisotropic(6)
{
	REGISTER(orientation, ECFMetaData()
			.setdescription({ "Orientation from elements' regions." })
			.setdatatype( { ECFDataType::BOOL } ));

	REGISTER(model, ECFMetaData()
			.setdescription({ "Tensors model." })
			.setdatatype({ ECFDataType::OPTION })
			.registertensor(poisson_ratio)
			.registertensor(young_modulus)
			.registertensor(shear_modulus)
			.registertensor(anisotropic)
			.addoption(ECFOption().setname("ISOTROPIC").setdescription("Isotropic model."))
			.addoption(ECFOption().setname("ORTHOTROPIC").setdescription("Orthotropic model."))
			.addoption(ECFOption().setname("ANISOTROPIC").setdescription("Anisotropic model.").allowonly([&] () { return *dimension == DIMENSION::D3; })));

	ecfdescription->addSeparator();

	ecfdescription->registerParameter("MIXY", poisson_ratio.get(0, 0), ECFMetaData()
			.setdescription({ "Poisson ratio XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(poisson_ratio, 0, 0)
			.allowonly([&] () { return model != MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("MIXZ", poisson_ratio.get(1, 1), ECFMetaData()
			.setdescription({ "Poisson ratio XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(poisson_ratio, 1, 1)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));
	ecfdescription->registerParameter("MIYZ", poisson_ratio.get(2, 2), ECFMetaData()
			.setdescription({ "Poisson ratio YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(poisson_ratio, 2, 2)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));

	ecfdescription->addSeparator()->metadata.allowonly([&] () { return model != MODEL::ISOTROPIC; });

	ecfdescription->registerParameter("EX", young_modulus.get(0, 0), ECFMetaData()
			.setdescription({ "Young modulus X." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(young_modulus, 0, 0)
			.allowonly([&] () { return model != MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("EY", young_modulus.get(1, 1), ECFMetaData()
			.setdescription({ "Young modulus Y." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(young_modulus, 1, 1)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC; }));
	ecfdescription->registerParameter("EZ", young_modulus.get(2, 2), ECFMetaData()
			.setdescription({ "Young modulus Z." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(young_modulus, 2, 2)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));

	ecfdescription->addSeparator()->metadata.allowonly([&] () { return model != MODEL::ISOTROPIC; });

	ecfdescription->registerParameter("GXY", shear_modulus.get(0, 0), ECFMetaData()
			.setdescription({ "Shear modulus XY." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(shear_modulus, 0, 0)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));
	ecfdescription->registerParameter("GXZ", shear_modulus.get(1, 1), ECFMetaData()
			.setdescription({ "Shear modulus XZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(shear_modulus, 1, 1)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));
	ecfdescription->registerParameter("GYZ", shear_modulus.get(2, 2), ECFMetaData()
			.setdescription({ "Shear modulus YZ." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(shear_modulus, 2, 2)
			.allowonly([&] () { return model == MODEL::ORTHOTROPIC && *dimension == DIMENSION::D3; }));

	// anisotropic is allowed only in 3D
	ecfdescription->registerParameter("D11", anisotropic.get(0, 0), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 1)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 0)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D12", anisotropic.get(0, 1), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 2)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 1)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D13", anisotropic.get(0, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 2)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D14", anisotropic.get(0, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 3)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D15", anisotropic.get(0, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 4)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D16", anisotropic.get(0, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (1, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 0, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D22", anisotropic.get(1, 1), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 2)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 1, 1)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D23", anisotropic.get(1, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 1, 2)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D24", anisotropic.get(1, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 1, 3)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D25", anisotropic.get(1, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 1, 4)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D26", anisotropic.get(1, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (2, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 1, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D33", anisotropic.get(2, 2), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 3)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 2, 2)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D34", anisotropic.get(2, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 2, 3)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D35", anisotropic.get(2, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 2, 4)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D36", anisotropic.get(2, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (3, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 2, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D44", anisotropic.get(3, 3), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 4)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 3, 3)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D45", anisotropic.get(3, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 3, 4)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D46", anisotropic.get(3, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (4, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 3, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D55", anisotropic.get(4, 4), ECFMetaData()
			.setdescription({ "Anisotropic parameter (5, 5)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 4, 4)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D56", anisotropic.get(4, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (5, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 4, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
	ecfdescription->registerParameter("D66", anisotropic.get(5, 5), ECFMetaData()
			.setdescription({ "Anisotropic parameter (6, 6)." })
			.setdatatype({ ECFDataType::EXPRESSION })
			.settensor(anisotropic, 5, 5)
			.allowonly([&] () { return model == MODEL::ANISOTROPIC; }));
}

bool espreso::LinearElasticPropertiesConfiguration::needCoordinates() const
{
	for (size_t i = 0; i < poisson_ratio.values.size(); ++i) {
		if (poisson_ratio.values[i].evaluator && poisson_ratio.values[i].evaluator->needCoordinates()) {
			return true;
		}
	}
	for (size_t i = 0; i < young_modulus.values.size(); ++i) {
		if (poisson_ratio.values[i].evaluator && poisson_ratio.values[i].evaluator->needCoordinates()) {
			return true;
		}
	}
	for (size_t i = 0; i < shear_modulus.values.size(); ++i) {
		if (poisson_ratio.values[i].evaluator && poisson_ratio.values[i].evaluator->needCoordinates()) {
			return true;
		}
	}
	for (size_t i = 0; i < anisotropic.values.size(); ++i) {
		if (anisotropic.values[i].evaluator && anisotropic.values[i].evaluator->needCoordinates()) {
			return true;
		}
	}

	return false;
}

bool espreso::LinearElasticPropertiesConfiguration::needTemperature() const
{
	for (size_t i = 0; i < poisson_ratio.values.size(); ++i) {
		if (poisson_ratio.values[i].evaluator && poisson_ratio.values[i].evaluator->needTemperature()) {
			return true;
		}
	}
	for (size_t i = 0; i < young_modulus.values.size(); ++i) {
		if (young_modulus.values[i].evaluator && young_modulus.values[i].evaluator->needTemperature()) {
			return true;
		}
	}
	for (size_t i = 0; i < shear_modulus.values.size(); ++i) {
		if (shear_modulus.values[i].evaluator && shear_modulus.values[i].evaluator->needTemperature()) {
			return true;
		}
	}
	for (size_t i = 0; i < anisotropic.values.size(); ++i) {
		if (anisotropic.values[i].evaluator && anisotropic.values[i].evaluator->needTemperature()) {
			return true;
		}
	}
	return false;
}

