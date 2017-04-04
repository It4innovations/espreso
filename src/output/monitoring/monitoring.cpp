
#include "monitoring.h"

#include "../../configuration/output.h"
#include "../../mesh/settings/property.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../assembler/step.h"
#include "../../assembler/solution.h"

#include "../../configuration/environment.h"

using namespace espreso::output;

char Monitoring::delimiter = ';';

static espreso::StatisticalData getStatistics(const std::string &name)
{
	if (espreso::StringCompare::caseInsensitiveEq(name, "AVG")) {
		return espreso::StatisticalData::AVERAGE;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "MIN")) {
		return espreso::StatisticalData::MIN;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "MAX")) {
		return espreso::StatisticalData::MAX;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "NORM")) {
		return espreso::StatisticalData::NORM;
	}

	ESINFO(espreso::GLOBAL_ERROR) << "Unknown operation " << name << "\n";
	return espreso::StatisticalData::AVERAGE;
}

static espreso::Property getProperty(const std::string &name)
{
	if (espreso::StringCompare::caseInsensitiveEq(name, "DISPLACEMENT_X")) {
		return espreso::Property::DISPLACEMENT_X;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "DISPLACEMENT_Y")) {
		return espreso::Property::DISPLACEMENT_Y;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "DISPLACEMENT_Z")) {
		return espreso::Property::DISPLACEMENT_Z;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "TEMPERATURE")) {
		return espreso::Property::TEMPERATURE;
	}

	ESINFO(espreso::GLOBAL_ERROR) << "Cannot monitor " << name << "\n";
	return espreso::Property::UNKNOWN;
}

std::string center(const std::string &value, size_t size)
{
	std::string ret;
	ret.insert(ret.end(), (size - value.size()) / 2 + (size - value.size()) % 2, ' ');
	ret.insert(ret.end(), value.begin(), value.end());
	ret.insert(ret.end(), (size - value.size()) / 2, ' ');
	return ret;
}

std::string right(const std::string &value, size_t size)
{
	std::string ret;
	ret.insert(ret.end(), size - value.size(), ' ');
	ret.insert(ret.end(), value.begin(), value.end());
	return ret;
}

Monitoring::Monitoring(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: Store(output), _mesh(mesh)
{
	size_t length = 0;
	_monitors.reserve(_configuration.monitoring.size());
	for (auto it = _configuration.monitoring.begin(); it != _configuration.monitoring.end(); ++it) {
		_monitors.push_back(Monitor());
		_monitors.back().region = _mesh->region(it->first);
		_mesh->addMonitoredRegion(_monitors.back().region);
		std::vector<std::string> args = Parser::split(it->second, " ");
		if (args.size() != 2) {
			ESINFO(GLOBAL_ERROR) << "Invalid monitoring format: use <REGION> <OPERATION> <VALUE>.";
		}
		_monitors.back().statistics = getStatistics(args[0]);
		_monitors.back().property = getProperty(args[1]);
		_monitors.back().printSize = std::max(std::max((size_t)10, it->first.size()), std::max(args[0].size(), args[1].size()) + 2) + 4;
		length += _monitors.back().printSize;
	}

	size_t rowHeaderSize = 19;

	if (environment->MPIrank) {
		return;
	}

	_os.open(path);
	if (!_os.is_open()) {
		ESINFO(GLOBAL_ERROR) << "Cannot open file " << path << "\n";
	}

	_os << "\n";
	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		_os << center(_monitors[i].region->name, _monitors[i].printSize) << delimiter;
	}
	_os << "\n";

	_os << right("step", 9) << delimiter << right("substep", 9) << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		std::stringstream ss;
		ss << _monitors[i].property;
		_os << center(ss.str(), _monitors[i].printSize) << delimiter;
	}
	_os << "\n";

	_os << std::string(9, ' ') << delimiter << std::string(9, ' ') << delimiter;
	for (size_t i = 0; i < _monitors.size(); i++) {
		switch (_monitors[i].statistics) {
		case StatisticalData::AVERAGE: _os << center("<AVERAGE>", _monitors[i].printSize) << delimiter; break;
		case StatisticalData::MIN:     _os << center("<MIN>"    , _monitors[i].printSize) << delimiter; break;
		case StatisticalData::MAX:     _os << center("<MAX>"    , _monitors[i].printSize) << delimiter; break;
		case StatisticalData::NORM:    _os << center("<NORM>"   , _monitors[i].printSize) << delimiter; break;
		}
	}
	_os << "\n\n";
}

void Monitoring::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	for (size_t i = 0; i < _monitors.size(); i++) {
		for (size_t s = 0; s < solution.size(); s++) {
			if (solution[s]->hasProperty(_monitors[i].property)) {
				solution[s]->computeStatisticalData();
			}
		}
	}

	if (environment->MPIrank) {
		return;
	}
	_os << right(std::to_string(step.step + 1), 8) << " " << delimiter;
	_os << right(std::to_string(step.substep + 1), 8) << " " << delimiter;

	for (size_t i = 0; i < _monitors.size(); i++) {
		double value;
		for (size_t s = 0; s < solution.size(); s++) {
			if (solution[s]->hasProperty(_monitors[i].property)) {
				value = solution[s]->getStatisticalData(_monitors[i].property, _monitors[i].statistics, _monitors[i].region);
			}
		}

		std::stringstream ss;
		ss << std::scientific << value;
		_os << center(ss.str(), _monitors[i].printSize) << delimiter;
	}
	_os << "\n";

}

void Monitoring::finalize()
{
	if (environment->MPIrank) {
		return;
	}
	_os.close();
}



