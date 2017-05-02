
#ifndef SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_
#define SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_

#include "../configuration.hpp"

namespace espreso {

enum class CONVECTION_TYPE {
	USER,
	EXTERNAL_NATURAL,
	INTERNAL_NATURAL,
	EXTERNAL_FORCED,
	INTERNAL_FORCED
};

enum class CONVECTION_VARIANT {
	VERTICAL_WALL,
	INCLINED_WALL,
	HORIZONTAL_CYLINDER,
	SPHERE,
	HORIZONTAL_PLATE_UP,
	HORIZONTAL_PLATE_DOWN,
	AVERAGE_PLATE,
	PARALLEL_PLATES,
	CIRCULAR_TUBE,
	TUBE
};

enum class CONVECTION_FLUID {
	AIR,
	WATER,
	OIL
};

struct AdvectionDiffusionConvection: public Configuration {

	OPTION(CONVECTION_TYPE, type, "Convection type", CONVECTION_TYPE::USER, OPTIONS({
		{ "USER"            , CONVECTION_TYPE::USER            , "User defined convection." },
		{ "EXTERNAL_NATURAL", CONVECTION_TYPE::EXTERNAL_NATURAL, "..." },
		{ "INTERNAL_NATURAL", CONVECTION_TYPE::INTERNAL_NATURAL, "..." },
		{ "EXTERNAL_FORCED" , CONVECTION_TYPE::EXTERNAL_FORCED , "..." },
		{ "INTERNAL_FORCED" , CONVECTION_TYPE::INTERNAL_FORCED , "..." }
	}));

	OPTION(CONVECTION_VARIANT, variant, "Type variant", CONVECTION_VARIANT::VERTICAL_WALL, OPTIONS({
		{ "VERTICAL_WALL"        , CONVECTION_VARIANT::VERTICAL_WALL        , "User defined convection." },
		{ "INCLINED_WALL"        , CONVECTION_VARIANT::INCLINED_WALL        , "..." },
		{ "HORIZONTAL_CYLINDER"  , CONVECTION_VARIANT::HORIZONTAL_CYLINDER  , "..." },
		{ "SPHERE"               , CONVECTION_VARIANT::SPHERE               , "..." },
		{ "HORIZONTAL_PLATE_UP"  , CONVECTION_VARIANT::HORIZONTAL_PLATE_UP  , "..." },
		{ "HORIZONTAL_PLATE_DOWN", CONVECTION_VARIANT::HORIZONTAL_PLATE_DOWN, "..." },
		{ "AVERAGE_PLATE"        , CONVECTION_VARIANT::AVERAGE_PLATE        , "..." },
		{ "PARALLEL_PLATES"      , CONVECTION_VARIANT::PARALLEL_PLATES      , "..." },
		{ "CIRCULAR_TUBE"        , CONVECTION_VARIANT::CIRCULAR_TUBE        , "..." },
		{ "TUBE"                 , CONVECTION_VARIANT::TUBE                 , "..." }
	}));

	PARAMETER(std::string, heat_transfer_coefficient, "Heat transfer coefficient.", "0");
	PARAMETER(std::string, external_temperature     , "External temperature."     , "0");

	PARAMETER(std::string, wall_height   , "Wall height."   , "0");
	PARAMETER(std::string, tilt_angle    , "Tilt angle."    , "0");
	PARAMETER(std::string, diameter      , "Diameter."      , "0");
	PARAMETER(std::string, plate_length  , "Plate length."  , "0");
	PARAMETER(std::string, fluid_velocity, "Fluid velocity.", "0");
	PARAMETER(std::string, plate_distance, "Plate distance.", "0");
	PARAMETER(std::string, length        , "Length."        , "0");

	OPTION(CONVECTION_FLUID, fluid, "Fluid type", CONVECTION_FLUID::AIR, OPTIONS({
		{ "AIR"  , CONVECTION_FLUID::AIR  , "Air." },
		{ "WATER", CONVECTION_FLUID::WATER, "Water." },
		{ "OIL"  , CONVECTION_FLUID::OIL  , "Oil." }
	}));

	PARAMETER(std::string, absolute_pressure, "Pressure.", "0");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_ */
