
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
};

struct AdvectionDiffusionConvection: public Configuration {

	OPTION(CONVECTION_TYPE, type, "Convection type", CONVECTION_TYPE::USER, OPTIONS({
		{ "USER"            , CONVECTION_TYPE::USER            , "User defined convection." },
		{ "EXTERNAL_NATURAL", CONVECTION_TYPE::EXTERNAL_NATURAL, "External natural convection" },
		{ "INTERNAL_NATURAL", CONVECTION_TYPE::INTERNAL_NATURAL, "Internal natural convection" },
		{ "EXTERNAL_FORCED" , CONVECTION_TYPE::EXTERNAL_FORCED , "External forced convection" },
		{ "INTERNAL_FORCED" , CONVECTION_TYPE::INTERNAL_FORCED , "Internal forced convection" }
	}));

	OPTION(CONVECTION_VARIANT, variant, "Type variant", CONVECTION_VARIANT::VERTICAL_WALL, OPTIONS({
		{ "VERTICAL_WALL"        , CONVECTION_VARIANT::VERTICAL_WALL        , "External natural convection - vertical wall" },
		{ "INCLINED_WALL"        , CONVECTION_VARIANT::INCLINED_WALL        , "External natural convection - inclined wall" },
		{ "HORIZONTAL_CYLINDER"  , CONVECTION_VARIANT::HORIZONTAL_CYLINDER  , "External natural convection - long horizontal cylinder" },
		{ "SPHERE"               , CONVECTION_VARIANT::SPHERE               , "External natural convection - sphere" },
		{ "HORIZONTAL_PLATE_UP"  , CONVECTION_VARIANT::HORIZONTAL_PLATE_UP  , "External natural convection - parallel plate, upside" },
		{ "HORIZONTAL_PLATE_DOWN", CONVECTION_VARIANT::HORIZONTAL_PLATE_DOWN, "External natural convection - parallel plate, downside" },
		{ "AVERAGE_PLATE"        , CONVECTION_VARIANT::AVERAGE_PLATE        , "External forced convection - plate averaged transfer coefficient" },
		{ "PARALLEL_PLATES"      , CONVECTION_VARIANT::PARALLEL_PLATES      , "Internal natural convection - parallel plates" },
		{ "CIRCULAR_TUBE"        , CONVECTION_VARIANT::CIRCULAR_TUBE        , "Internal natural convection - circular tube" },
		{ "TUBE"                 , CONVECTION_VARIANT::TUBE                 , "Internal forced convection - isothermal tube" }
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
	}));

	PARAMETER(std::string, absolute_pressure, "Pressure.", "0");
};

}



#endif /* SRC_CONFIGURATION_PHYSICS_ADVECTIONDIFFUSIONCONVECTION_H_ */
