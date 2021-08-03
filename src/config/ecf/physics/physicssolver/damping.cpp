
#include "damping.h"
#include "config/configuration.hpp"

using namespace espreso;

DirectDampingConfiguration::DirectDampingConfiguration()
: mass(ECFExpression::Scope::EGPS), stiffness(ECFExpression::Scope::EGPS)
{
	stiffness.value = "0";
	REGISTER(stiffness, ECFMetaData()
			.setdescription({ "Stiffnes Damping" })
			.setdatatype({ ECFDataType::FLOAT }));

	mass.value = "0";
	REGISTER(mass, ECFMetaData()
			.setdescription({ "Mass Damping" })
			.setdatatype({ ECFDataType::FLOAT }));
}

RatioDampingConfiguration::RatioDampingConfiguration()
: ratio(ECFExpression::Scope::EGPS), frequency(ECFExpression::Scope::EGPS)
{
	ratio.value = "0";
	REGISTER(ratio, ECFMetaData()
			.setdescription({ "Damping Ratio" })
			.setdatatype({ ECFDataType::FLOAT }));

	frequency.value = "0";
	REGISTER(frequency, ECFMetaData()
			.setdescription({ "Damping frequency" })
			.setdatatype({ ECFDataType::FLOAT }));
}

RayleighDampingConfiguration::RayleighDampingConfiguration()
{
	type = TYPE::NONE;
	REGISTER(type, ECFMetaData()
			.setdescription({ "Damping type" })
			.setdatatype({ ECFDataType::OPTION })
			.addoption(ECFOption().setname("NONE").setdescription("Without damping."))
			.addoption(ECFOption().setname("DIRECT").setdescription("Direct."))
			.addoption(ECFOption().setname("DAMPING_RATIO").setdescription("By ratio.")));


	REGISTER(direct_damping, ECFMetaData()
			.setdescription({ "Direct damping configuration." }));
	REGISTER(ratio_damping, ECFMetaData()
			.setdescription({ "Ratio damping configuration." }));
}

HarmonicRayleighDampingConfiguration::HarmonicRayleighDampingConfiguration()
: structural_damping_coefficient(ECFExpression::Scope::EGPS)
{
	structural_damping_coefficient.value = "0";
	REGISTER(structural_damping_coefficient, ECFMetaData()
			.setdescription({ "Structural Damping Coefficient" })
			.setdatatype({ ECFDataType::FLOAT }));
}

CoriolisEffectAxisConfiguration::CoriolisEffectAxisConfiguration()
: x(0), y(0), z(1)
{
	REGISTER(x, ECFMetaData()
		.setdescription({ "X-axis." })
		.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(y, ECFMetaData()
		.setdescription({ "Y-axis." })
		.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(z, ECFMetaData()
		.setdescription({ "Z-axis." })
		.setdatatype({ ECFDataType::FLOAT }));
}

CoriolisEffectConfiguration::CoriolisEffectConfiguration()
{
	coriolis_damping = spin_softening = false;

	REGISTER(coriolis_damping, ECFMetaData()
		.setdescription({ "Turn on/off coriolis damping." })
		.setdatatype({ ECFDataType::FLOAT }));
	REGISTER(spin_softening, ECFMetaData()
		.setdescription({ "Turn on/off spin softening." })
		.setdatatype({ ECFDataType::FLOAT }));

	REGISTER(rotation_axis, ECFMetaData()
			.setdescription({ "Rotation axis configuration." }));
}

DampingConfiguration::DampingConfiguration()
{
	REGISTER(rayleigh, ECFMetaData()
			.setdescription({ "RAYLEIGH damping configuration." }));
//	REGISTER(coriolis_effect, ECFMetaData()
//			.setdescription({ "Coriolis effect configuration." }));
}

HarmonicDampingConfiguration::HarmonicDampingConfiguration()
{
	REGISTER(rayleigh, ECFMetaData()
			.setdescription({ "RAYLEIGH damping configuration." }));
//	REGISTER(coriolis_effect, ECFMetaData()
//			.setdescription({ "Coriolis effect configuration." }));
}





