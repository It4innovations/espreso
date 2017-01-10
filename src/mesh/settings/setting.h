
#ifndef SRC_MESH_SETTINGS_SETTING_H_
#define SRC_MESH_SETTINGS_SETTING_H_


#include "evaluator.h"

namespace espreso {

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

	std::vector<Evaluator*>& property(Property property)
	{
		return _settings[property];
	}

	const std::vector<Evaluator*>& property(Property property) const
	{
		if (_settings.find(property) != _settings.end()) {
			return _settings.find(property)->second;
		} else {
			return _settings.find(Property::EMPTY)->second;
		}
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
			for (eslocal e = 0; e < eSize; e++) {
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
