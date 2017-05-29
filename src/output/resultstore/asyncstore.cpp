
#include <iostream>
#include "logging/logging.h"
#include "vtklegacy.h"
#include "vtkxmlascii.h"
#include "vtkxmlbinary.h"
#include "asyncstore.h"
#include "../../configuration/output.h"

using namespace espreso::output;

void AsyncStore::init(const Mesh *mesh)
{
	async::Module<AsyncStoreExecutor, OutputConfiguration, Param>::init();

	// Setup filename buffer (hopefully 1024 chars are enough)
	addBuffer(0L, 1024);

	callInit(configuration());

	switch (configuration().format) {
	case espreso::OUTPUT_FORMAT::VTK_LEGACY:
		_headerStore = new VTKLegacy(configuration(), 0L, "");
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_ASCII:
		_headerStore = new VTKXMLASCII(configuration(), 0L, "");
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_BINARY:
		_headerStore = new VTKXMLBinary(configuration(), 0L, "");
		break;
	default:
		assert(false);
	}
}

std::string AsyncStore::store(const std::string &name, const RegionData &regionData)
{
	// Create the buffer the first time we send data
	// Note: This has to be done before calling wait() since addBuffer() has an implicit wait()/call()
	if (_bufferSize) {
		if (_bufferSize != regionData.packedSize())
			ESINFO(GLOBAL_ERROR) << "Packed data size changed?";
	} else {
		_bufferSize = regionData.packedSize();
		addBuffer(0L, _bufferSize);
	}

	// Wait for the last I/O call to finish
	wait();

	// Send file name
	if (name.size()+1 > 1024)
		ESINFO(GLOBAL_ERROR) << "File name to long (> 1024)";
	memcpy(managedBuffer<char*>(0), name.c_str(), name.size()+1);
	sendBuffer(0);

	// Pack data
	regionData.pack(managedBuffer<char*>(1));

	// Send data
	sendBuffer(1);

	Param param;
	call(param);

	return name + ".vtu";
}

std::string AsyncStore::linkClusters(const std::string &root, const std::string &name, const RegionData &regionData)
{
	return _headerStore->linkClusters(root, name, regionData);
}

void AsyncStore::linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps)
{
	_headerStore->linkSteps(name, steps);
}

void AsyncStore::finalize()
{
	if (_finalized)
		return;

	// Wait for the last I/O call to finish
	wait();

	// Finalize the I/O module
	async::Module<AsyncStoreExecutor, OutputConfiguration, Param>::finalize();

	// Link files
	ResultStore::finalize();

	_finalized = true;
}
