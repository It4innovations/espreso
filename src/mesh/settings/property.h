
#ifndef SRC_MESH_SETTINGS_PROPERTY_H_
#define SRC_MESH_SETTINGS_PROPERTY_H_

#include <sstream>

namespace espreso {

enum class Property : int {
	UNKNOWN, // API has unknown properties
	DISPLACEMENT_X,
	DISPLACEMENT_Y,
	DISPLACEMENT_Z,
	MOMENTUM_X,
	MOMENTUM_Y,
	MOMENTUM_Z,
	TEMPERATURE,
	PRESSURE,

	THICKNESS,
	INITIAL_TEMPERATURE,
	FORCE_X,
	FORCE_Y,
	FORCE_Z,
	ACCELERATION_X,
	ACCELERATION_Y,
	ACCELERATION_Z,
	HEAT_SOURCE,
	TRANSLATION_MOTION_X,
	TRANSLATION_MOTION_Y,
	TRANSLATION_MOTION_Z,
	HEAT_FLOW,
	HEAT_FLUX,
	HEAT_TRANSFER_COEFFICIENT,
	EXTERNAL_TEMPERATURE,

	WALL_HEIGHT,
	TILT_ANGLE,
	DIAMETER,
	PLATE_LENGTH,
	FLUID_VELOCITY,
	PLATE_DISTANCE,
	LENGTH,
	ABSOLUTE_PRESSURE,

	THERMAL_RESISTANCE,
	SURFACE_TEMPERATURE,
	LAYER_THICKNESS,
	LAYER_CONDUCTIVITY,

	OBSTACLE,
	NORMAL_DIRECTION,
	NONMATCHING_ELEMENT,
	EMPTY,

	SIZE
};

inline std::ostream& operator<<(std::ostream& os, const Property& property)
{
	switch (property) {
		case Property::UNKNOWN: return os << "UNKNOWN";
		case Property::DISPLACEMENT_X: return os << "DISPLACEMENT_X";
		case Property::DISPLACEMENT_Y: return os << "DISPLACEMENT_Y";
		case Property::DISPLACEMENT_Z: return os << "DISPLACEMENT_Z";
		case Property::MOMENTUM_X: return os << "MOMENTUM_X";
		case Property::MOMENTUM_Y: return os << "MOMENTUM_Y";
		case Property::MOMENTUM_Z: return os << "MOMENTUM_Z";
		case Property::TEMPERATURE: return os << "TEMPERATURE";
		case Property::PRESSURE: return os << "PRESSURE";

		case Property::THICKNESS: return os << "THICKNESS";
		case Property::INITIAL_TEMPERATURE: return os << "INITIAL_TEMPERATURE";
		case Property::FORCE_X: return os << "FORCE_X";
		case Property::FORCE_Y: return os << "FORCE_Y";
		case Property::FORCE_Z: return os << "FORCE_Z";
		case Property::ACCELERATION_X: return os << "ACCELERATION_X";
		case Property::ACCELERATION_Y: return os << "ACCELERATION_Y";
		case Property::ACCELERATION_Z: return os << "ACCELERATION_Z";
		case Property::HEAT_SOURCE: return os << "HEAT_SOURCE";
		case Property::HEAT_FLOW: return os << "HEAT_FLOW";
		case Property::TRANSLATION_MOTION_X: return os << "TRANSLATION_MOTION_X";
		case Property::TRANSLATION_MOTION_Y: return os << "TRANSLATION_MOTION_Y";
		case Property::TRANSLATION_MOTION_Z: return os << "TRANSLATION_MOTION_Z";
		case Property::HEAT_FLUX: return os << "HEAT_FLUX";
		case Property::HEAT_TRANSFER_COEFFICIENT: return os << "HEAT_TRANSFER_COEFFICIENT";
		case Property::EXTERNAL_TEMPERATURE: return os << "EXTERNAL_TEMPERATURE";

		case Property::WALL_HEIGHT: return os << "WALL_HEIGHT";
		case Property::TILT_ANGLE: return os << "TILT_ANGLE";
		case Property::DIAMETER: return os << "DIAMETER";
		case Property::PLATE_LENGTH: return os << "PLATE_LENGTH";
		case Property::FLUID_VELOCITY: return os << "FLUID_VELOCITY";
		case Property::PLATE_DISTANCE: return os << "PLATE_DISTANCE";
		case Property::LENGTH: return os << "LENGTH";
		case Property::ABSOLUTE_PRESSURE: return os << "ABSOLUTE_PRESSURE";

		case Property::THERMAL_RESISTANCE: return os << "THERMAL_RESISTANCE";
		case Property::LAYER_THICKNESS: return os << "LAYER_THICKNESS";
		case Property::LAYER_CONDUCTIVITY: return os << "LAYER_CONDUCTIVITY";

		case Property::OBSTACLE: return os << "OBSTACLE";
		case Property::NORMAL_DIRECTION: return os << "NORMAL_DIRECTION";
		case Property::NONMATCHING_ELEMENT: return os << "NONMATCHING_ELEMENT";
		default: return os;
	}
}

}

#endif /* SRC_MESH_SETTINGS_PROPERTY_H_ */
