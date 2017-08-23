
#include "monitoring.h"

#include "../../configuration/output.h"
#include "../../mesh/settings/property.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../assembler/step.h"
#include "../../assembler/solution.h"

#include "../../configuration/environment.h"

using namespace espreso;

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

std::vector<espreso::Property> Monitoring::getProperties(const std::string &name)
{
	for (auto p = _mesh->propertyGroups().begin(); p != _mesh->propertyGroups().end(); ++p) {
		if (p->second.size() > 1) {
			std::stringstream ss; ss << p->first;
			std::string pname = ss.str().substr(0, ss.str().find_last_of("_"));
			if (StringCompare::caseInsensitiveEq(pname, name)) {
				return p->second;
			}
		}
	}

	std::stringstream ss(name);
	espreso::Property property;
	ss >> property;

	return { property };
}

Monitoring::Monitoring(const OutputConfiguration &output, const Mesh *mesh)
: Store(output), _mesh(mesh)
{
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
		_monitors.back().properties = getProperties(args[1]);
		_monitors.back().printSize = std::max(std::max((size_t)10, it->first.size()), std::max(args[0].size(), args[1].size()) + 2) + 4;
	}

	if (environment->MPIrank) {
		return;
	}

	_os.open(Logging::outputRoot() + "/" + Logging::name + ".emr");
	if (!_os.is_open()) {
		ESINFO(GLOBAL_ERROR) << "Cannot open file for storing monitor report\n";
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
		if (_monitors[i].properties.size() > 1) {
			std::stringstream ssp; ssp << _monitors[i].properties[0];
			ss << ssp.str().substr(0, ssp.str().find_last_of("_"));
		} else {
			ss << _monitors[i].properties[0];
		}
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
		default: break;
		}
	}
	_os << "\n\n";
}

void Monitoring::updateMesh()
{
	// probably empty?
}

void Monitoring::storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties)
{
	for (size_t i = 0; i < _monitors.size(); i++) {
		for (size_t s = 0; s < solution.size(); s++) {
			if (solution[s]->hasProperty(_monitors[i].properties[0])) {
				solution[s]->computeStatisticalData(step);
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
		bool found = false;
		for (size_t s = 0; s < solution.size(); s++) {
			if (solution[s]->hasProperty(_monitors[i].properties[0])) {
				value = solution[s]->getStatisticalData(_monitors[i].properties, _monitors[i].statistics, _monitors[i].region);
				found = true;
				break;
			}
		}
		if (!found) {
			ESINFO(GLOBAL_ERROR) << "ESPRESO monitor error: request for unknown property: " << _monitors[i].properties[0];
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



