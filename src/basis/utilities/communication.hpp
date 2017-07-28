
#include "mpi.h"

#include "communication.h"
#include "../../configuration/environment.h"

namespace espreso {

template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		MPI_Isend(sBuffer[n].data(), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + n);
	}

	int flag;
	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(Ttype));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] > environment->MPIrank) {
			MPI_Isend(sBuffer[n].data(), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
		if (neighbours[n] < environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(Ttype) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size(); n++) {
		if (neighbours[n] < environment->MPIrank) {
			MPI_Isend(sBuffer[n].data(), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
		if (neighbours[n] > environment->MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(Ttype) * rBuffer[n].size(), MPI_BYTE, neighbours[n], 1, environment->MPICommunicator, req.data() + n);
		}
	}

	MPI_Waitall(neighbours.size(), req.data(), MPI_STATUSES_IGNORE);
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbours)
{
	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbours.size());
	for (size_t n = 0; n < neighbours.size() && neighbours[n] < environment->MPIrank; n++) {
		MPI_Isend(sBuffer[n].data(), sizeof(Ttype) * sBuffer[n].size(), MPI_BYTE, neighbours[n], 0, environment->MPICommunicator, req.data() + rSize++);
	}

	int flag;
	size_t counter = std::lower_bound(neighbours.begin(), neighbours.end(), environment->MPIrank) - neighbours.begin();
	MPI_Status status;
	while (counter < neighbours.size()) {
		MPI_Iprobe(MPI_ANY_SOURCE, 0, environment->MPICommunicator, &flag, &status);
		if (flag) {
			int count;
			MPI_Get_count(&status, MPI_BYTE, &count);
			rBuffer[n2i(status.MPI_SOURCE)].resize(count / sizeof(Ttype));
			MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, MPI_BYTE, status.MPI_SOURCE, 0, environment->MPICommunicator, MPI_STATUS_IGNORE);
			counter++;
		}
	}

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(environment->MPICommunicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	return true;
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer)
{
	std::vector<size_t> offsets;
	return gatherUnknownSize(sBuffer, rBuffer, offsets);
}

template <typename Ttype>
bool Communication::gatherUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, std::vector<size_t> &offsets)
{
	int size = sizeof(Ttype) * sBuffer.size();
	std::vector<int> rSizes(environment->MPIsize), rOffsets(environment->MPIsize);
	MPI_Gather(&size, sizeof(int), MPI_BYTE, rSizes.data(), sizeof(int), MPI_BYTE, 0, environment->MPICommunicator);

	if (!environment->MPIrank) {
		size = 0;
		for (size_t i = 0; i < rSizes.size(); i++) {
			rOffsets[i] = size;
			size += rSizes[i];
		}
		rBuffer.resize(size / sizeof(Ttype));
	}

	MPI_Gatherv(sBuffer.data(), sBuffer.size() * sizeof(Ttype), MPI_BYTE, rBuffer.data(), rSizes.data(), rOffsets.data(), MPI_BYTE, 0, environment->MPICommunicator);

	offsets.resize(environment->MPIsize);
	for (size_t i = 0; i < rOffsets.size(); i++) {
		offsets[i] = rOffsets[i] / sizeof(Ttype);
	}
	return true;
}

template <typename Ttype>
bool Communication::broadcastUnknownSize(std::vector<Ttype> &buffer)
{
	int size = buffer.size();
	MPI_Bcast(&size, sizeof(int), MPI_BYTE, 0, environment->MPICommunicator);
	buffer.resize(size);
	MPI_Bcast(buffer.data(), sizeof(Ttype) * size, MPI_BYTE, 0, environment->MPICommunicator);
	return true;
}

template <typename Ttype>
static void offsetSum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	*(static_cast<Ttype*>(out)) += *(static_cast<Ttype*>(in));
}

template <typename Ttype>
Ttype Communication::exscan(Ttype &value)
{
	size_t size = value;
	if (environment->MPIsize == 1) {
		value = 0;
		return size;
	}

	MPI_Op op;
	MPI_Op_create(offsetSum<Ttype>, 1, &op);
	MPI_Exscan(&size, &value, sizeof(size_t), MPI_BYTE, op, environment->MPICommunicator);

	size = value + size;
	MPI_Bcast(&size, sizeof(size_t), MPI_BYTE, environment->MPIsize - 1, environment->MPICommunicator);
	if (environment->MPIrank == 0) {
		value = 0;
	}

	return size;
}

}



