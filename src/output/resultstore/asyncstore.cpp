
#include <iostream>
#include "asyncstore.h"

using namespace espreso::output;

void AsyncStore::init(const Mesh *mesh, const std::string &path)
{
	async::Module<AsyncStoreExecutor, OutputConfiguration, Step>::init();

	// Setup buffers
	addSyncBuffer(path.c_str(), path.size(), true);
	addSyncBuffer(&mesh, sizeof(Mesh*), true); // TODO Will not work with MPI
	std::cout << "Mesh: " << mesh << std::endl;

	// Send buffer required for initialization
	sendBuffer(0);
	sendBuffer(1);

	callInit(configuration());

	// Delete buffer no longer in use
	removeBuffer(0);
	removeBuffer(1);

	std::cout << "init. Done." << std::endl;
}

void AsyncStore::finalize()
{
	if (_finalized)
		return;

	std::cout << "AsyncStore::finalize()" << std::endl;
	// Wait for the last I/O call to finish
	wait();

	// Finalize the I/O module
	async::Module<AsyncStoreExecutor, OutputConfiguration, Step>::finalize();

	_finalized = true;
}
