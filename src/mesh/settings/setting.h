
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

	THICKNESS,
	INITIAL_TEMPERATURE,
	ACCELERATION_X,
	ACCELERATION_Y,
	ACCELERATION_Z,
	HEAT_SOURCE,
	TRANSLATION_MOTION_X,
	TRANSLATION_MOTION_Y,
	TRANSLATION_MOTION_Z,

	NONMATCHING_ELEMENT,
	EMPTY
};

inline std::ostream& operator<<(std::ostream& os, const Property& property)
{
	switch (property) {
		case Property::DISPLACEMENT_X: return os << "DISPLACEMENT_X";
		case Property::DISPLACEMENT_Y: return os << "DISPLACEMENT_Y";
		case Property::DISPLACEMENT_Z: return os << "DISPLACEMENT_Z";
		case Property::TEMPERATURE: return os << "TEMPERATURE";
		case Property::PRESSURE: return os << "PRESSURE";

		case Property::THICKNESS: return os << "THICKNESS";
		case Property::INITIAL_TEMPERATURE: return os << "INITIAL_TEMPERATURE";
		case Property::ACCELERATION_X: return os << "ACCELERATION_X";
		case Property::ACCELERATION_Y: return os << "ACCELERATION_Y";
		case Property::ACCELERATION_Z: return os << "ACCELERATION_Z";
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

	void store(std::ofstream &os, const std::vector<Evaluator*> &evaluators)
	{
		eslocal index, size = _settings.size() - 1;
		os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
		for (auto it = _settings.begin(); it != _settings.end(); ++it) {
			if (it->first == Property::EMPTY) {
				continue;
			}
			os.write(reinterpret_cast<const char*>(&(it->first)), sizeof(Property));
			size = it->second.size();
			os.write(reinterpret_cast<const char*>(&size), sizeof(eslocal));
			for (size_t i = 0; i < it->second.size(); i++) {
				index = std::find(evaluators.begin(), evaluators.end(), it->second[i]) - evaluators.begin();
				os.write(reinterpret_cast<const char*>(&index), sizeof(eslocal));
			}
		}
	}

	void load(std::ifstream &is, const std::vector<Evaluator*> &evaluators)
	{
		eslocal size;
		is.read(reinterpret_cast<char *>(&size), sizeof(eslocal));
		for (eslocal i = 0; i < size; i++) {
			Property property;
			is.read(reinterpret_cast<char*>(&(property)), sizeof(Property));
			eslocal eSize;
			is.read(reinterpret_cast<char *>(&eSize), sizeof(eslocal));
			for (size_t e = 0; e < eSize; e++) {
				eslocal index;
				is.read(reinterpret_cast<char *>(&index), sizeof(eslocal));
				_settings[property].push_back(evaluators[index]);
			}
		}
	}

private:
	std::map<Property, std::vector<Evaluator*> > _settings;
};

}



#endif /* SRC_MESH_SETTINGS_SETTING_H_ */
