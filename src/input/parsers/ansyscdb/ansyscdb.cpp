
#include "ansyscdb.h"
#include "parser/blockend.h"
#include "parser/cm.h"
#include "parser/cmblock.h"
#include "parser/eblock.h"
#include "parser/esel.h"
#include "parser/et.h"
#include "parser/nblock.h"
#include "parser/nsel.h"
#include "esinfo/eslog.h"
#include "basis/logging/profiler.h"
#include "config/ecf/input/input.h"
#include "input/parsers/distributedscanner.h"

using namespace espreso;

AnsysCDBLoader::AnsysCDBLoader(const InputConfiguration &configuration)
: _configuration(configuration), _file({ _configuration.path }, 32 * MAX_LINE_STEP * MAX_LINE_SIZE, MAX_LINE_STEP * MAX_LINE_SIZE)
{

}

void AnsysCDBLoader::load()
{
	eslog::startln("ANSYS CDB PARSER: STARTED", "ANSYS CDB PARSER");
	profiler::syncstart("ansys_cdb");

	_file.prepare();
	profiler::synccheckpoint("prepare_reader");
	eslog::checkpointln("ANSYS CDB PARSER: READER PREPARED");

	_file.read();
	profiler::synccheckpoint("read");
	eslog::checkpointln("ANSYS CDB PARSER: DATA READ");

	_file.next();
	DistributedScanner::align(_file, "\n");

	WorkbenchParser::offset = _file.distribution[info::mpi::rank];
	WorkbenchParser::begin = _file.begin;
	WorkbenchParser::end = _file.end;

	scan();
	profiler::synccheckpoint("scan");
	eslog::checkpointln("ANSYS CDB PARSER: DATA SCANNED");

	parse();
	if (!_configuration.keep_material_sets) {
		std::fill(material.begin(), material.end(), 0);
	}
	profiler::synccheckpoint("parse");
	profiler::syncend("ansys_cdb");
	eslog::endln("ANSYS CDB PARSER: DATA PARSED");
}

template<typename T>
void setends(T &t, std::vector<BlockEnd> &ends, std::vector<size_t> &distribution)
{
	for (size_t i = 0; i < t.size(); i++) {
		t[i].fillDistribution(ends, distribution);
	}
}

void AnsysCDBLoader::scan()
{
	DistributedScanner parser;

	parser.add("n,", [&] (const char *c) { _blockEnds.push_back(BlockEnd().parse(c)); });
	parser.add("N,", [&] (const char *c) { _blockEnds.push_back(BlockEnd().parse(c)); });
	parser.add("-1", [&] (const char *c) { _blockEnds.push_back(BlockEnd().parse(c)); });
	parser.add("nblock", [&] (const char *c) { _NBlocks.push_back(NBlock().parse(c)); });
	parser.add("NBLOCK", [&] (const char *c) { _NBlocks.push_back(NBlock().parse(c)); });
	parser.add("eblock", [&] (const char *c) { _EBlocks.push_back(EBlock().parse(c)); });
	parser.add("EBLOCK", [&] (const char *c) { _EBlocks.push_back(EBlock().parse(c)); });
	parser.add("cmblock", [&] (const char *c) { _CMBlocks.push_back(CMBlock().parse(c)); });
	parser.add("CMBLOCK", [&] (const char *c) { _CMBlocks.push_back(CMBlock().parse(c)); });
	parser.add("et,", [&] (const char *c) { _ET.push_back(ET().parse(c)); });
	parser.add("ET,", [&] (const char *c) { _ET.push_back(ET().parse(c)); });
	parser.add("esel,", [&] (const char *c) { _ESel.push_back(ESel().parse(c)); });
	parser.add("ESEL,", [&] (const char *c) { _ESel.push_back(ESel().parse(c)); });
	parser.add("cm,", [&] (const char *c) { _CM.push_back(CM().parse(c)); });
	parser.add("CM,", [&] (const char *c) { _CM.push_back(CM().parse(c)); });

	// old Ansys can have '-1' with leading spaces
	parser.addEnd(" -1", [&] (const char *c) { _blockEnds.push_back(BlockEnd().parse(c)); });

	parser.scanlines(_file);
	parser.synchronize(_blockEnds, _NBlocks, _EBlocks, _CMBlocks, _ET, _ESel, _CM);

	setends(_NBlocks, _blockEnds, _file.distribution);
	setends(_EBlocks, _blockEnds, _file.distribution);
	setends(_CMBlocks, _blockEnds, _file.distribution);
	setends(_ET, _blockEnds, _file.distribution);
	setends(_ESel, _blockEnds, _file.distribution);
	setends(_NSel, _blockEnds, _file.distribution);
	setends(_CM, _blockEnds, _file.distribution);
}

void AnsysCDBLoader::parse()
{
	std::vector<ET> et;
	esint maxet = 0;
	for (size_t i = 0; i < _ET.size(); i++) {
		if (maxet < _ET[i].id) {
			maxet = _ET[i].id;
		}
	}
	et.resize(maxet + 1);
	for (size_t i = 0; i < _ET.size(); i++) {
		if (_ET[i].id >= 0) {
			et[_ET[i].id] = _ET[i];
		}
	}
	_ET.swap(et);

	for (size_t i = 0; i < _NBlocks.size(); i++) {
		if (!_NBlocks[i].readData(*this)) {
			eslog::globalerror("Workbench parser: something wrong happens while read NBLOCK.\n");
		}
	}

	for (size_t i = 0; i < _EBlocks.size(); i++) {
		if (!_EBlocks[i].readData(_ET, *this)) {
			eslog::globalerror("Workbench parser: something wrong happens while read EBLOCK.\n");
		}
	}

	for (size_t i = 0; i < _CMBlocks.size(); i++) {
		switch (_CMBlocks[i].entity) {
		case CMBlock::Entity::NODE:
			if (!_CMBlocks[i].readData(nregions[_CMBlocks[i].name])) {
				eslog::globalerror("Workbench parser: something wrong happens while read CMBLOCK.\n");
			}
			break;
		case CMBlock::Entity::ELEMENT:
			if (!_CMBlocks[i].readData(eregions[_CMBlocks[i].name])) {
				eslog::globalerror("Workbench parser: something wrong happens while read CMBLOCK.\n");
			}
			break;
		default:
			break;
		}
	}

	for (size_t i = 0; i < _CM.size(); i++) {
		_CM[i].addRegion(*this, _ESel, eregions, _NSel, nregions);
	}
}
