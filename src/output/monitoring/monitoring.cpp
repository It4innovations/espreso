
#include "monitoring.h"

#include "../../configuration/output.h"
#include "../../mesh/settings/property.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/region.h"
#include "../../assembler/step.h"
#include "../../assembler/solution.h"

using namespace espreso::output;

static espreso::output::Operation getOperation(const std::string &name)
{
	if (espreso::StringCompare::caseInsensitiveEq(name, "AVERAGE")) {
		return espreso::output::Operation::AVERAGE;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "MIN")) {
		return espreso::output::Operation::MIN;
	}
	if (espreso::StringCompare::caseInsensitiveEq(name, "MAX")) {
		return espreso::output::Operation::MAX;
	}

	ESINFO(espreso::GLOBAL_ERROR) << "Unknown operation " << name << "\n";
	return espreso::output::Operation::AVERAGE;
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
	_os.open(path);
	if (!_os.is_open()) {
		ESINFO(GLOBAL_ERROR) << "Cannot open file " << path << "\n";
	}

	size_t length = 0;
	_monitors.reserve(_configuration.monitoring.size());
	for (auto it = _configuration.monitoring.begin(); it != _configuration.monitoring.end(); ++it) {
		_monitors.push_back(Monitor());
		_monitors.back().region = _mesh->region(it->first);
		std::vector<std::string> args = Parser::split(it->second, " ");
		if (args.size() != 2) {
			ESINFO(GLOBAL_ERROR) << "Invalid monitoring format: use <REGION> <OPERATION> <VALUE>.";
		}
		_monitors.back().operation = getOperation(args[0]);
		_monitors.back().property = getProperty(args[1]);
		_monitors.back().printSize = std::max(std::max((size_t)10, it->first.size()), std::max(args[0].size(), args[1].size()) + 2) + 4;
		length += _monitors.back().printSize;
	}

	size_t rowHeaderSize = 19;

	_os << std::string(length + _monitors.size() + 1 + rowHeaderSize, '-') << "\n";
	_os << std::string(rowHeaderSize, ' ') << "|";
	for (size_t i = 0; i < _monitors.size(); i++) {
		_os << center(_monitors[i].region->name, _monitors[i].printSize) << "|";
	}
	_os << "\n";

	_os << right("step", 9) << right("substep", 9) << " |";
	for (size_t i = 0; i < _monitors.size(); i++) {
		std::stringstream ss;
		ss << _monitors[i].property;
		_os << center(ss.str(), _monitors[i].printSize) << "|";
	}
	_os << "\n";

	_os << std::string(rowHeaderSize, ' ') << "|";
	for (size_t i = 0; i < _monitors.size(); i++) {
		switch (_monitors[i].operation) {
		case Operation::AVERAGE: _os << center("<AVERAGE>", _monitors[i].printSize) << "|"; break;
		case Operation::MIN:     _os << center("<MIN>"    , _monitors[i].printSize) << "|"; break;
		case Operation::MAX:     _os << center("<MAX>"    , _monitors[i].printSize) << "|"; break;
		}
	}
	_os << "\n";

	_os << std::string(length + _monitors.size() + 1 + rowHeaderSize, '-') << "\n";
}

void Monitoring::storeSolution(const Step &step, const std::vector<Solution*> &solution)
{
	_os << right(std::to_string(step.step), 9);
	_os << right(std::to_string(step.substep), 9);
	_os << " |";

	for (size_t i = 0; i < _monitors.size(); i++) {
		double value; // TODO: get required value
		std::stringstream ss;
		ss << std::scientific << value;
		_os << center(ss.str(), _monitors[i].printSize) << "|";
	}
	_os << "\n";

}

void Monitoring::finalize()
{
	_os.close();
}



