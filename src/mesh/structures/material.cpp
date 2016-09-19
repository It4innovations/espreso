
#include "material.h"

using namespace espreso;


bool Material::setParameter(const std::string &parameter, const std::string &value)
{
	std::cout << parameter << ": " << value << "\n";
	bool correctlySet = false;

	auto set = [&] (Evaluator* &evaluator, const std::string &name) {
		if (StringCompare::caseInsensitiveEq(parameter, name)) {
			delete evaluator;
			if (value.find("xyz") == std::string::npos) {
				espreso::Expression expr(value, {});
				evaluator = new ConstEvaluator(expr.evaluate({}));
			} else {
				evaluator = new CoordinatesEvaluator(value, *_coordinates);
			}
			correctlySet = true;
		}
	};

	set(_density            , "DENS");
	set(_poissonRatio[0]    , "MIXY");
	set(_poissonRatio[0]    , "NUXY");
	set(_poissonRatio[1]    , "MIXZ");
	set(_poissonRatio[0]    , "NUXZ");
	set(_poissonRatio[2]    , "MIYZ");
	set(_poissonRatio[0]    , "NUYZ");
	set(_youngModulus[0]    , "EX");
	set(_youngModulus[1]    , "EY");
	set(_youngModulus[2]    , "EZ");
	set(_termalCapacity     , "CP");
	set(_termalCapacity     , "C");
	set(_termalConduction[0], "KXX");
	set(_termalConduction[1], "KXY");
	set(_termalConduction[2], "KXZ");
	set(_termalExpansion[0] , "ALPX");
	set(_termalExpansion[1] , "ALPY");
	set(_termalExpansion[2] , "ALPZ");

	if (StringCompare::caseInsensitiveEq(parameter, "RSVX")) {
		ESINFO(ALWAYS) << TextColor::YELLOW << "Skipped material parameter 'RSVX'";
		return true;
	}
	if (StringCompare::caseInsensitiveEq(parameter, "RSVX")) {
		ESINFO(ALWAYS) << TextColor::YELLOW << "Skipped material parameter 'RSVX'";
		return true;
	}

	if (StringCompare::caseInsensitiveEq(parameter, "MODEL")) {
		Parameter p("MODEL", _model, "model", {
				{"LINEAR_ELASTIC_ISOTROPIC", MODEL::LINEAR_ELASTIC_ISOTROPIC, "Isotropic"},
				{"LINEAR_ELASTIC_ORTHOTROPIC", MODEL::LINEAR_ELASTIC_ORTHOTROPIC, "Orthotropic"},
				{"LINEAR_ELASTIC_ANISOTROPIC", MODEL::LINEAR_ELASTIC_ANISOTROPIC, "Anisotropic"},
		});
		p.set(value);
	}

	return correctlySet;
}




