
#include <iostream>

#include "../../basis/logging/logging.h"
#include "../../configuration/output.h"
#include "../../configuration/environment.h"
#include "../basis/utilities/utils.h"

#include "vtklegacy.h"
#include "vtkxmlascii.h"
#include "vtkxmlbinary.h"
#include "asyncstore.h"


using namespace espreso;

void AsyncStore::init(const Mesh *mesh)
{
	async::Module<AsyncStoreExecutor, OutputConfiguration, Param>::init();

	// Setup filename buffer (hopefully 1024 chars are enough)
	addBuffer(0L, 1024);

	callInit(configuration());

	switch (configuration().format) {
	case espreso::OUTPUT_FORMAT::VTK_LEGACY:
		_headerStore = new VTKLegacy(configuration(), 0L);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_ASCII:
		_headerStore = new VTKXMLASCII(configuration(), 0L);
		break;
	case espreso::OUTPUT_FORMAT::VTK_XML_BINARY:
		_headerStore = new VTKXMLBinary(configuration(), 0L);
		break;
	default:
		assert(false);
	}
}

std::vector<std::string> AsyncStore::store(const std::string &name, const Step &step, const MeshInfo *meshInfo)
{
	std::string root;
	if (_configuration.subsolution) {
		root = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), "iteration" + std::to_string(step.iteration)});
	} else {
		root = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep) });
	}
	std::vector<std::string> files;

	if (meshInfo->distributed()) {
		std::string prefix;
		if (_configuration.subsolution) {
			prefix = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), "iteration" + std::to_string(step.iteration), std::to_string(environment->MPIrank) });
		} else {
			prefix = Esutils::createDirectory({ Logging::outputRoot(), "PRE_POST_DATA", "step" + std::to_string(step.step), "substep" + std::to_string(step.substep), std::to_string(environment->MPIrank) });
		}

		size_t bSize = sizeof(size_t);
		for (size_t r = 0; r < meshInfo->regions(); r++) {
			bSize += meshInfo->region(r).packedSize();
		}

		if (_bufferSize) {
			if (_bufferSize != bSize) {
				ESINFO(ERROR) << "Packed data size changed?";
			}
		} else {
			_bufferSize = bSize;
			addBuffer(0L, _bufferSize);
		}

		wait();

		std::string bName;
		for (size_t r = 0; r < meshInfo->regions(); r++) {
			bName += prefix + name + std::to_string(r) + ";";
		}

		if (bName.size() + 1 > 1024) {
			ESINFO(ERROR) << "File name to long (> 1024)";
		}
		memcpy(managedBuffer<char*>(0), bName.c_str(), bName.size() + 1);
		sendBuffer(0);

		char *p = managedBuffer<char*>(1);
		size_t datasize = meshInfo->regions();
		memcpy(p, &datasize, sizeof(size_t));
		p += sizeof(size_t);

		for (size_t r = 0; r < meshInfo->regions(); r++) {
			meshInfo->region(r).pack(p);
			p += meshInfo->region(r).packedSize();
			if (!environment->MPIrank) {
				files.push_back(linkClusters(root, name + std::to_string(r), meshInfo->region(r)));
			}
		}

		sendBuffer(1);

		Param param;
		call(param);

	} else {
		ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: COLLECTED output is not implemented for ASYNC MODE.";
	}

	return files;
}

std::string AsyncStore::store(const std::string &name, const RegionData &regionData)
{
	// Create the buffer the first time we send data
	// Note: This has to be done before calling wait() since addBuffer() has an implicit wait()/call()
	if (_bufferSize) {
		if (_bufferSize != regionData.packedSize()) {
			ESINFO(ERROR) << "Packed data size changed?";
		}
	} else {
		_bufferSize = regionData.packedSize();
		addBuffer(0L, _bufferSize);
	}

	// Wait for the last I/O call to finish
	wait();

	// Send file name
	if (name.size() + 1 > 1024) {
		ESINFO(ERROR) << "File name to long (> 1024)";
	}
	memcpy(managedBuffer<char*>(0), name.c_str(), name.size() + 1);
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
