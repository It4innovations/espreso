
#include "mpiinfo.h"
#include "eslog.hpp"
#include "basis/utilities/communication.h"

int espreso::info::mpi::rank = 0;
int espreso::info::mpi::size = 1;
MPI_Comm espreso::info::mpi::comm = MPI_COMM_WORLD;

int espreso::info::mpi::irank = 0;
int espreso::info::mpi::isize = 1;
MPI_Comm espreso::info::mpi::icomm = MPI_COMM_SELF;

int espreso::info::mpi::grank = 0;
int espreso::info::mpi::gsize = 1;
MPI_Comm espreso::info::mpi::gcomm = MPI_COMM_WORLD;

int espreso::info::mpi::threading = 0;

using namespace espreso::info;

void mpi::init(int *argc, char ***argv)
{
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &threading);
	set();
}

void mpi::set()
{
	int initialized;
	MPI_Initialized(&initialized);

	mpi::comm = MPI_COMM_WORLD;
	if (initialized) {
		MPI_Comm_rank(mpi::comm, &mpi::rank);
		MPI_Comm_size(mpi::comm, &mpi::size);
	}

	mpi::grank = mpi::rank;
	mpi::gsize = mpi::size;
	mpi::gcomm = mpi::comm;
}

bool mpi::divide(int meshDuplication)
{
	if (meshDuplication == 1) {
		return true;
	}

	if (espreso::info::mpi::size % meshDuplication != 0) {
		return false;
	}

	int color = mpi::rank / (mpi::size / meshDuplication);

	MPI_Comm_split(mpi::gcomm, color, mpi::grank, &mpi::comm);
	MPI_Comm_rank(mpi::comm, &mpi::rank);
	MPI_Comm_size(mpi::comm, &mpi::size);

	MPI_Comm_split(mpi::gcomm, mpi::rank, mpi::grank, &mpi::icomm);
	MPI_Comm_rank(mpi::icomm, &mpi::irank);
	MPI_Comm_size(mpi::icomm, &mpi::isize);

	MPITools::reinit();
	eslog::reinit();

	return true;
}

void mpi::finish()
{
	if (mpi::isize > 1) {
		MPI_Comm_free(&mpi::comm);
		MPI_Comm_free(&mpi::icomm);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}


