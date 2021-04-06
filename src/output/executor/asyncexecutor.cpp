
#include "asyncexecutor.h"

#include "basis/utilities/utils.h"
#include "config/ecf/ecf.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

#ifdef HAVE_ASYNC

std::vector<int> AsyncBufferManager::_indices(Buffer::SIZE, -1);

void AsyncStore::prepareBuffer(AsyncBufferManager::Buffer buffer, size_t size)
{
	int index = AsyncBufferManager::buffer(buffer);
	if (index == -1) {
		AsyncBufferManager::buffer(buffer, index = addBuffer(NULL, size));
	}
	if (bufferSize(index) != size) {
		removeBuffer(index);
		AsyncBufferManager::buffer(buffer, index = addBuffer(NULL, size));
	}
}

void AsyncStore::updateMesh()
{
	wait();

	prepareBuffer(AsyncBufferManager::NODES, _mesh.nodes->packedSize());
	_mesh.nodes->pack(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODES)));

	prepareBuffer(AsyncBufferManager::ELEMENTS, _mesh.elements->packedSize());
	_mesh.elements->pack(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS)));

	{ // ELEMENT REGIONS
		size_t esize = 0;
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			esize += _mesh.elementsRegions[r]->packedSize();
		}
		prepareBuffer(AsyncBufferManager::ELEMENTREGIONS, esize);

		_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTREGIONS));
		for (size_t r = 0; r < _mesh.elementsRegions.size(); r++) {
			_mesh.elementsRegions[r]->pack(_buffer);
		}
	}

	{ // BOUNDARY REGIONS
		size_t bsize = 0;
		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			bsize += _mesh.boundaryRegions[r]->packedSize();
		}
		prepareBuffer(AsyncBufferManager::BOUNDARYREGIONS, bsize);

		_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::BOUNDARYREGIONS));
		for (size_t r = 0; r < _mesh.boundaryRegions.size(); r++) {
			_mesh.boundaryRegions[r]->pack(_buffer);
		}
	}

	call(ExecParameters(AsyncBufferManager::NODES, AsyncBufferManager::ELEMENTS, AsyncBufferManager::ELEMENTREGIONS, AsyncBufferManager::BOUNDARYREGIONS));
}

void AsyncStore::updateMonitors()
{
	wait();

	prepareBuffer(AsyncBufferManager::NODEDATA, _mesh.nodes->packedDataHeaderSize());
	_mesh.nodes->packDataHeader(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA)));

	prepareBuffer(AsyncBufferManager::ELEMENTDATA, _mesh.elements->packedDataHeaderSize());
	_mesh.elements->packDataHeader(_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTDATA)));

	call(ExecParameters(AsyncBufferManager::NODEDATA, AsyncBufferManager::ELEMENTDATA));
}

void AsyncStore::updateSolution()
{
	wait();

	prepareBuffer(AsyncBufferManager::NODESOLUTION, _mesh.nodes->packedDataSize());
	_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::NODESOLUTION));
	_mesh.nodes->packData(_buffer);

	prepareBuffer(AsyncBufferManager::ELEMENTSOLUTION, _mesh.elements->packedDataSize());
	_buffer = managedBuffer<char*>(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTSOLUTION));
	_mesh.elements->packData(_buffer);

	call(ExecParameters(AsyncBufferManager::NODESOLUTION, AsyncBufferManager::ELEMENTSOLUTION));
}

AsyncStore::AsyncStore(const Mesh &mesh)
: ResultStoreExecutor(mesh), _executor(mesh.store), _buffer(NULL)
{
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::init();
	callInit(InitParameters());
}

void AsyncStore::addResultStore(ResultStoreBase *resultStore)
{
	_executor.addResultStore(resultStore);
}

bool AsyncStore::hasStore()
{
	return _executor.hasStore();
}

AsyncStore::~AsyncStore()
{
	wait();
	async::Module<AsyncExecutor, InitParameters, ExecParameters>::finalize();
}

AsyncExecutor::AsyncExecutor(ResultStore *store)
: DirectExecutor(_mesh), _buffer(NULL)
{
	_mesh.store = store;
}

void AsyncExecutor::execInit(const async::ExecInfo &info, const InitParameters &initParameters)
{

}

void AsyncExecutor::exec(const async::ExecInfo &info, const ExecParameters &parameters)
{
	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODES) {
		_mesh.nodes->unpack(_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODES))));
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS) {
		_mesh.elements->unpack(_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTS))));
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTREGIONS) {
		int bufferid = AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTREGIONS);
		_buffer = static_cast<const char*>(info.buffer(bufferid));
		while (_buffer < static_cast<const char*>(info.buffer(bufferid)) + info.bufferSize(bufferid)) {
			_mesh.elementsRegions.push_back(new ElementsRegionStore(""));
			_mesh.elementsRegions.back()->unpack(_buffer);
		}
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::BOUNDARYREGIONS) {
		int bufferid = AsyncBufferManager::buffer(AsyncBufferManager::BOUNDARYREGIONS);
		_buffer = static_cast<const char*>(info.buffer(bufferid));
		while (_buffer < static_cast<const char*>(info.buffer(bufferid)) + info.bufferSize(bufferid)) {
			_mesh.boundaryRegions.push_back(new BoundaryRegionStore(""));
			_mesh.boundaryRegions.back()->unpack(_buffer);
		}
	}

	if (
			(parameters.updatedBuffers & 1 << AsyncBufferManager::NODES) ||
			(parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTS)) {

		updateMesh();
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODEDATA) {
		_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODEDATA)));
		_mesh.nodes->unpackDataHeader(_buffer);
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTDATA) {
		_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTDATA)));
		_mesh.elements->unpackDataHeader(_buffer);
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::NODESOLUTION) {
		_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::NODESOLUTION)));
		_mesh.nodes->unpackData(_buffer);
	}

	if (parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTSOLUTION) {
		_buffer = static_cast<const char*>(info.buffer(AsyncBufferManager::buffer(AsyncBufferManager::ELEMENTSOLUTION)));
		_mesh.elements->unpackData(_buffer);
	}

	if (
			(parameters.updatedBuffers & 1 << AsyncBufferManager::NODESOLUTION) ||
			(parameters.updatedBuffers & 1 << AsyncBufferManager::ELEMENTSOLUTION)) {

		updateSolution();
	}
}

#endif /* HAVE ASYNC */
