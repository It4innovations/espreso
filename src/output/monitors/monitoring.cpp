
#include "monitoring.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/parser.h"

#include "config/ecf/output.h"

#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"

#include <iomanip>
#include <sstream>

using namespace espreso;

char Monitoring::delimiter = ';';

int pre(const std::string &value, int size)
{
	return (size - value.size()) / 2 + (size - value.size()) % 2 + value.size();
}

int post(const std::string &value, int size)
{
	return (size - value.size()) / 2 + 1;
}

bool Monitoring::storeStep(const step::Step &step)
{
	if (step.type == step::TYPE::FTT) {
		return true;
	} else {
		switch (info::ecf->output.monitors_store_frequency) {
		case OutputConfiguration::STORE_FREQUENCY::NEVER:
			return false;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_SUBSTEP:
			return true;
		case OutputConfiguration::STORE_FREQUENCY::EVERY_NTH_SUBSTEP:
			return step.substep % info::ecf->output.monitors_nth_stepping == 0;
		case OutputConfiguration::STORE_FREQUENCY::LAST_SUBSTEP:
			return step::isLast(step);
		default:
			return false;
		}
	}
}


Monitoring::Monitoring()
: _runFile(NULL), _fttFile(NULL)
{

}

Monitoring::~Monitoring()
{
	if (_runFile) {
		fclose(_runFile);
	}
}

void Monitoring::updateMonitors(const step::Step &step)
{
	for (auto it = info::ecf->output.monitoring.begin(); it != info::ecf->output.monitoring.end(); ++it) {
		if (it->first <= 0) {
			eslog::globalerror("Invalid column index in monitoring.\n");
		}
		if (it->first > _monitors.size()) {
			_monitors.resize(it->first);
		}

		ElementsRegionStore *estore = NULL;
		BoundaryRegionStore *bstore = NULL;
		bool regionNotFound = true;
		for (size_t r = 0; regionNotFound && r < info::mesh->elementsRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, info::mesh->elementsRegions[r]->name)) {
				estore = info::mesh->elementsRegions[r];
				regionNotFound = false;
			}
		}
		for (size_t r = 0; regionNotFound && r < info::mesh->boundaryRegions.size(); r++) {
			if (StringCompare::caseInsensitiveEq(it->second.region, info::mesh->boundaryRegions[r]->name)) {
				bstore = info::mesh->boundaryRegions[r];
				regionNotFound = false;
			}
		}
		if (regionNotFound) {
			eslog::globalerror("Monitoring contains unknown region '%s'.\n", it->second.region.c_str());
		}

		NodeData *ndata = NULL;
		ElementData *edata = NULL;
		bool propertyNotFound = true;
		for (size_t i = 0; i < info::mesh->nodes->data.size(); i++) {
			if (info::mesh->nodes->data[i]->name.size()) {
				if (!info::mesh->nodes->data[i]->onlySuffixed() && StringCompare::caseInsensitiveEq(it->second.property, info::mesh->nodes->data[i]->name)) {
					ndata = info::mesh->nodes->data[i];
					propertyNotFound = false;
				}
				for (int p = 0; p < info::mesh->nodes->data[i]->dimension; p++) {
					if (StringCompare::caseInsensitiveEq(it->second.property, info::mesh->nodes->data[i]->name + info::mesh->nodes->data[i]->suffix(p))) {
						ndata = info::mesh->nodes->data[i];
						propertyNotFound = false;
					}
				}
			}
		}

		for (size_t i = 0; i < info::mesh->elements->data.size(); i++) {
			if (info::mesh->elements->data[i]->name.size()) {
				if (!info::mesh->elements->data[i]->onlySuffixed() && StringCompare::caseInsensitiveEq(it->second.property, info::mesh->elements->data[i]->name)) {
					edata = info::mesh->elements->data[i];
					propertyNotFound = false;
				}
				for (int p = 0; p < info::mesh->elements->data[i]->dimension; p++) {
					if (StringCompare::caseInsensitiveEq(it->second.property, info::mesh->elements->data[i]->name + info::mesh->elements->data[i]->suffix(p))) {
						edata = info::mesh->elements->data[i];
						propertyNotFound = false;
					}
				}
			}
		}
		if (propertyNotFound) {
			eslog::globalerror("Monitoring contains unknown property '%s'.\n", it->second.property.c_str());
		}

		if (edata != NULL && bstore != NULL) {
			eslog::globalerror("Cannot monitor element property '%s' on element region '%s'.\n", it->second.property.c_str(), it->second.region.c_str());
		}
		if (edata != NULL && estore != NULL) {
			for (size_t i = 0; i < _edata.size(); i++) {
				if (_edata[i].first->name == edata->name && _edata[i].second->name == estore->name) {
					edata = NULL;
					estore = NULL;
					break;
				}
			}
		}
		if (ndata != NULL && bstore != NULL) {
			for (size_t i = 0; i < _nbdata.size(); i++) {
				if (_nbdata[i].first->name == ndata->name && _nbdata[i].second->name == bstore->name) {
					ndata = NULL;
					bstore = NULL;
					break;
				}
			}
		}
		if (ndata != NULL && estore != NULL) {
			for (size_t i = 0; i < _nedata.size(); i++) {
				if (_nedata[i].first->name == ndata->name && _nedata[i].second->name == estore->name) {
					ndata = NULL;
					estore = NULL;
					break;
				}
			}
		}

		if (edata != NULL && estore != NULL) {
			_edata.push_back(std::make_pair(edata, estore));
		}
		if (ndata != NULL && bstore != NULL) {
			_nbdata.push_back(std::make_pair(ndata, bstore));
		}
		if (ndata != NULL && estore != NULL) {
			_nedata.push_back(std::make_pair(ndata, estore));
		}
	}

	for (size_t i = 0; i < _edata.size(); i++) {
		_statistics.resize(_statistics.size() + _edata[i].first->nstatistics());
	}
	for (size_t i = 0; i < _nbdata.size(); i++) {
		_statistics.resize(_statistics.size() + _nbdata[i].first->nstatistics());
	}
	for (size_t i = 0; i < _nedata.size(); i++) {
		_statistics.resize(_statistics.size() + _nedata[i].first->nstatistics());
	}

	for (auto it = info::ecf->output.monitoring.begin(); it != info::ecf->output.monitoring.end(); ++it) {
		_monitors[it->first - 1].name = it->second.region;
		_monitors[it->first - 1].property = it->second.property;
		_monitors[it->first - 1].printSize = std::max(std::max((size_t)11, it->second.region.size()), it->second.property.size()) + 4;

		esint offset = 0;
		for (size_t i = 0; i < _edata.size(); i++) {
			if (!_edata[i].first->onlySuffixed()) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _edata[i].first->name) &&
						StringCompare::caseInsensitiveEq(it->second.region, _edata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
				++offset;
			}
			for (int p = 0; _edata[i].first->withSuffixes() && p < _edata[i].first->dimension; ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _edata[i].first->name + _edata[i].first->suffix(p)) &&
						StringCompare::caseInsensitiveEq(it->second.region, _edata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}
		for (size_t i = 0; i < _nbdata.size(); i++) {
			if (!_nbdata[i].first->onlySuffixed()) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nbdata[i].first->name) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nbdata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
				++offset;
			}
			for (int p = 0; _nbdata[i].first->withSuffixes() && p < _nbdata[i].first->dimension; ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nbdata[i].first->name + _nbdata[i].first->suffix(p)) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nbdata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}
		for (size_t i = 0; i < _nedata.size(); i++) {
			if (!_nedata[i].first->onlySuffixed()) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nedata[i].first->name) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nedata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
				++offset;
			}
			for (int p = 0; _nedata[i].first->withSuffixes() && p < _nedata[i].first->dimension; ++p, ++offset) {
				if (
						StringCompare::caseInsensitiveEq(it->second.property, _nedata[i].first->name + _nedata[i].first->suffix(p)) &&
						StringCompare::caseInsensitiveEq(it->second.region, _nedata[i].second->name)) {

					_monitors[it->first - 1].data = (double*)(_statistics.data() + offset);
				}
			}
		}

		switch (it->second.statistics) {
		case MonitorConfiguration::STATISTICS::MIN:
			_monitors[it->first - 1].stats = "<MIN>"; break;
		case MonitorConfiguration::STATISTICS::MAX:
			_monitors[it->first - 1].data += 1;
			_monitors[it->first - 1].stats = "<MAX>"; break;
		case MonitorConfiguration::STATISTICS::AVG:
			_monitors[it->first - 1].data += 2;
			_monitors[it->first - 1].stats = "<AVERAGE>"; break;
		case MonitorConfiguration::STATISTICS::NORM:
			_monitors[it->first - 1].data += 3;
			_monitors[it->first - 1].stats = "<NORM>"; break;
		case MonitorConfiguration::STATISTICS::ABSMIN:
			_monitors[it->first - 1].data += 4;
			_monitors[it->first - 1].stats = "<ABSMIN>"; break;
		case MonitorConfiguration::STATISTICS::ABSMAX:
			_monitors[it->first - 1].data += 5;
			_monitors[it->first - 1].stats = "<ABSMAX>"; break;
		}
	}

	if (info::mpi::rank == 0) {
		_runFile = fopen((_path + "/" + std::string(info::ecf->name.c_str()) + ".emr").c_str(), "w");

		auto center = [&] (const std::string &value, int size) {
			fprintf(_runFile, "%*s%*c", pre(value, size), value.c_str(), post(value, size), delimiter);
		};

		// 1. line with region names
		fprintf(_runFile, "%10c %10c %14c", delimiter, delimiter, delimiter);
		for (size_t i = 0; i < _monitors.size(); i++) {
			center(_monitors[i].name, _monitors[i].printSize);
		}
		fprintf(_runFile, "\n");

		// 2. line with parameters
		fprintf(_runFile, "%8s %c %8s %c ", "loadstep", delimiter, "substep", delimiter);
		switch (step.type) {
		case step::TYPE::TIME:
			fprintf(_runFile, "%12s %c", "time", delimiter); break;
		case step::TYPE::FREQUENCY:
			fprintf(_runFile, "%12s %c", "frequency", delimiter); break;
		default: break;
		}
		for (size_t i = 0; i < _monitors.size(); i++) {
			center(_monitors[i].property, _monitors[i].printSize);
		}
		fprintf(_runFile, "\n");

		// 3. line with statistics
		fprintf(_runFile, "%10c %10c %14c", delimiter, delimiter, delimiter);
		for (size_t i = 0; i < _monitors.size(); i++) {
			center(_monitors[i].stats, _monitors[i].printSize);
		}
		fprintf(_runFile, "\n\n");

//		for (int i = 0; i < step::outduplicate.offset; i++) {
//			fprintf(_runFile, "%10c %10c %14c", delimiter, delimiter, delimiter);
//			for (size_t i = 0; i < _monitors.size(); i++) {
//				center(_monitors[i].stats, _monitors[i].printSize);
//			}
//			fprintf(_runFile, "\n\n");
//		}

		fflush(_runFile);
	}
}

void Monitoring::updateSolution(const step::Step &step, const step::Time &time)
{
	updateSolution(step);

	if (info::mpi::rank == 0) {
		fprintf(_runFile, "%8d %c %8d %c ", step.loadstep + 1, delimiter, step.substep + 1, delimiter);
		fprintf(_runFile, "%12.6f %c ", time.current, delimiter);
		storeSolution(step);
	}
}

void Monitoring::updateSolution(const step::Step &step, const step::Frequency &frequency)
{
	updateSolution(step);

	if (info::mpi::rank == 0) {
		fprintf(_runFile, "%8d %c %8d %c ", step.loadstep + 1, delimiter, step.substep + 1, delimiter);
		fprintf(_runFile, "%12.4f %c ", frequency.current, delimiter);
		storeSolution(step);
	}
}

void Monitoring::updateSolution(const step::Step &step)
{
	if (!storeStep(step)) {
		return;
	}

	esint offset = 0;
	for (size_t i = 0; i < _edata.size(); offset += _edata[i++].first->nstatistics()) {
		_edata[i].first->statistics(_edata[i].second->elements->datatarray(), _edata[i].second->nodeInfo.totalSize, _statistics.data() + offset);
	}
	for (size_t i = 0; i < _nbdata.size(); offset += _nbdata[i++].first->nstatistics()) {
		_nbdata[i].first->statistics(_nbdata[i].second->nodes->datatarray(), _nbdata[i].second->nodeInfo.totalSize, _statistics.data() + offset);
	}
	for (size_t i = 0; i < _nedata.size(); offset += _nedata[i++].first->nstatistics()) {
		_nedata[i].first->statistics(_nedata[i].second->nodes->datatarray(), _nedata[i].second->nodeInfo.totalSize, _statistics.data() + offset);
	}
}

void Monitoring::storeSolution(const step::Step &step)
{
	for (size_t i = 0; i < _monitors.size(); i++) {
		if (_monitors[i].data != NULL) {
			fprintf(_runFile, "% *e %c ", _monitors[i].printSize - 2, *_monitors[i].data, delimiter);
		} else {
			fprintf(_runFile, "%*c ", _monitors[i].printSize, delimiter);
		}
	}
	fprintf(_runFile, "\n");
	fflush(_runFile);
}







