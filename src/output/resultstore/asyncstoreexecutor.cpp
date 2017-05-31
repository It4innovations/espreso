
#include <cstring>
#include <cassert>
#include <iostream>
#include "asyncstoreexecutor.h"
#include "vtklegacy.h"
#include "vtkxmlascii.h"
#include "vtkxmlbinary.h"
#include "../../configuration/output.h"

using namespace espreso::output;

void AsyncStoreExecutor::execInit(const async::ExecInfo &info, const OutputConfiguration &config)
{
	assert(config.results || config.settings);

	switch (config.format) {
	case espreso::OUTPUT_FORMAT::VTK_LEGACY:
		_store = new VTKLegacy(config, 0L);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_ASCII:
		_store = new VTKXMLASCII(config, 0L);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_BINARY:
		_store = new VTKXMLBinary(config, 0L);
		break;
	default:
		assert(false);
	}
}

void AsyncStoreExecutor::exec(const async::ExecInfo &info, const Param &param)
{
	const char* regionBuffer = static_cast<const char*>(info.buffer(1));
	for (size_t i = 0; i < info.bufferSize(0) / 1024; i++) {
		std::string names(static_cast<const char*>(info.buffer(0) + i * 1024));
		std::vector<std::string> name = Parser::split(names, ";");

		size_t regions;
		memcpy(&regions, regionBuffer, sizeof(size_t));
		regionBuffer += sizeof(size_t);
		for (size_t r = 0; r < regions; r++) {
			RegionData region;
			region.unpack(regionBuffer);

			_store->store(name[r], region);
		}
	}
}
