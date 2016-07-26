
#ifndef SRC_MESH_SETTINGS_SETTING_H_
#define SRC_MESH_SETTINGS_SETTING_H_

#include "evaluator.h"

namespace espreso {

enum class Property : int {
	FIXED_DISPLACEMENT_X,
	FIXED_DISPLACEMENT_Y,
	FIXED_DISPLACEMENT_Z,
	FIXED_TEMPERATURE,
	FIXED_PRESSURE,

	DISPLACEMENT_X,
	DISPLACEMENT_Y,
	DISPLACEMENT_Z,
	TEMPERATURE,
	PRESSURE,

	HEAT_SOURCE,
	TRANSLATION_MOTION_X,
	TRANSLATION_MOTION_Y,
	TRANSLATION_MOTION_Z,
	EMPTY
};

inline std::ostream& espreso::operator<<(std::ostream& os, const Property& property)
{
	switch (property) {
		case Property::FIXED_DISPLACEMENT_X: return os << "FIXED_DISPLACEMENT_X";
		case Property::FIXED_DISPLACEMENT_Y: return os << "FIXED_DISPLACEMENT_Y";
		case Property::FIXED_DISPLACEMENT_Z: return os << "FIXED_DISPLACEMENT_Z";
		case Property::FIXED_TEMPERATURE: return os << "FIXED_TEMPERATURE";
		case Property::FIXED_PRESSURE: return os << "FIXED_PRESSURE";

		case Property::DISPLACEMENT_X: return os << "DISPLACEMENT_X";
		case Property::DISPLACEMENT_Y: return os << "DISPLACEMENT_Y";
		case Property::DISPLACEMENT_Z: return os << "DISPLACEMENT_Z";
		case Property::TEMPERATURE: return os << "TEMPERATURE";
		case Property::PRESSURE: return os << "PRESSURE";

		case Property::HEAT_SOURCE: return os << "HEAT_SOURCE";
		case Property::TRANSLATION_MOTION_X: return os << "TRANSLATION_MOTION_X";
		case Property::TRANSLATION_MOTION_Y: return os << "TRANSLATION_MOTION_Y";
		case Property::TRANSLATION_MOTION_Z: return os << "TRANSLATION_MOTION_Z";
		default: return os;
	}
}


class Settings {

public:
	Settings()
	{
		_settings[Property::EMPTY].push_back(new Evaluator());
	}

	~Settings()
	{
		delete _settings[Property::EMPTY][0];
	}

	std::vector<Evaluator*>& operator[](Property property)
	{
		return _settings[property];
	}

	const std::vector<Evaluator*>& operator[](Property property) const
	{
		if (_settings.find(property) != _settings.end()) {
			return _settings.find(property)->second;
		} else {
			return _settings.find(Property::EMPTY)->second;
		}
	}

	bool isSet(Property property) const
	{
		return _settings.find(property) != _settings.end();
	}

private:
	std::map<Property, std::vector<Evaluator*> > _settings;
};

}



#endif /* SRC_MESH_SETTINGS_SETTING_H_ */
