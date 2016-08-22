
#ifndef SRC_MESH_SETTINGS_SETTING_H_
#define SRC_MESH_SETTINGS_SETTING_H_

#include "evaluator.h"

namespace espreso {

enum class Property : int {
	DISPLACEMENT_X,
	DISPLACEMENT_Y,
	DISPLACEMENT_Z,
	TEMPERATURE,
	PRESSURE,

	HEAT_SOURCE,
	TRANSLATION_MOTION_X,
	TRANSLATION_MOTION_Y,
	TRANSLATION_MOTION_Z,

	NONMATCHING_ELEMENT,
	EMPTY
};

inline std::ostream& espreso::operator<<(std::ostream& os, const Property& property)
{
	switch (property) {
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

	Settings(const Settings &other)
	{
		_settings[Property::EMPTY].push_back(new Evaluator());
	}

	Settings& operator=(const Settings &other)
	{
		if (this != &other) {
			_settings[Property::EMPTY].push_back(new Evaluator());
		}
		return *this;
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

	size_t size() const { return _settings.size() - 1; }

private:
	std::map<Property, std::vector<Evaluator*> > _settings;
};

}



#endif /* SRC_MESH_SETTINGS_SETTING_H_ */
