
#include <cassert>
#include <iostream>
#include "asyncstoreexecutor.h"
#include "vtklegacy.h"
#include "vtkxmlascii.h"
#include "vtkxmlbinary.h"
#include "../../assembler/step.h"
#include "../../configuration/output.h"

using namespace espreso::output;

void AsyncStoreExecutor::execInit(const async::ExecInfo &info, const OutputConfiguration &config)
{
	assert(config.results || config.settings);

	std::string path(static_cast<const char*>(info.buffer(0)), info.bufferSize(0));

	// Reassemble the mesh
	const Mesh* mesh = *reinterpret_cast<const Mesh* const *>(info.buffer(1));

	std::cout << "execInit format: " << static_cast<int>(config.format) << std::endl;
	std::cout << mesh << ' ' << path << std::endl;
	switch (config.format) {
	case espreso::OUTPUT_FORMAT::VTK_LEGACY:
		_store = new VTKLegacy(config, mesh, path);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_ASCII:
		_store = new VTKXMLASCII(config, mesh, path);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_BINARY:
		_store = new VTKXMLBinary(config, mesh, path);
		break;
	default:
		assert(false);
	}

	std::cout << "execInit. Done." << std::endl;
}

void AsyncStoreExecutor::exec(const async::ExecInfo &info, const Step &step)
{
}
