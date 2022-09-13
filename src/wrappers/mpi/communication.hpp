
#include "communication.h"
#include "basis/logging/profiler.h"

#include <algorithm>
#include <cmath>

namespace espreso {

template <>
inline MPIType MPITools::getType<char>()
{
	return { MPI_CHAR };
}

template <>
inline MPIType MPITools::getType<short>()
{
	return { MPI_SHORT };
}

template <>
inline MPIType MPITools::getType<int>()
{
	return { MPI_INT };
}

template <>
inline MPIType MPITools::getType<uint>()
{
	return { MPI_UNSIGNED };
}

template <>
inline MPIType MPITools::getType<long>()
{
	return { MPI_LONG };
}

template <>
inline MPIType MPITools::getType<size_t>()
{
	return { MPI_UNSIGNED_LONG };
}

template <>
inline MPIType MPITools::getType<float>()
{
	return { MPI_FLOAT };
}

template <>
inline MPIType MPITools::getType<double>()
{
	return { MPI_DOUBLE, };
}

// WARNING: be care of this (MPI can divide long message in arbitrary position)
template <typename Ttype>
inline MPIType MPITools::getType()
{
	return { MPI_BYTE, sizeof(Ttype), false };
}

template <typename Ttype>
bool Communication::exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_exchange_known_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());

	for (size_t n = 0; n < neighbors.size(); n++) {
		if (type.mpisize * sBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_exchange_known_size");
			return false;
		}
		if (type.mpisize * rBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_exchange_known_size");
			return false;
		}
	}
	std::vector<MPI_Request> req(2 * neighbors.size());

	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::EX_KNOWN, group->communicator, req.data() + 2 * n);
		smin.min(type.mpisize * sBuffer[n].size());
		smax.max(type.mpisize * sBuffer[n].size());
	}

	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Irecv(const_cast<Ttype*>(rBuffer[n].data()), type.mpisize * rBuffer[n].size(), type.mpitype, neighbors[n], TAG::EX_KNOWN, group->communicator, req.data() + 2 * n + 1);
		rmin.min(type.mpisize * rBuffer[n].size());
		rmax.max(type.mpisize * rBuffer[n].size());
	}
	profiler::synccheckpoint("msg_prepared");

	MPI_Waitall(2 * neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::EX_KNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_exchange_known_size");
	return true;
}

template <typename Ttype>
bool Communication::exchangeKnownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_exchange_known_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());

	if (type.mpisize * sBuffer.size() > 1 << 30) {
		profiler::syncend("comm_exchange_known_size");
		return false;
	}
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (type.mpisize * rBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_exchange_known_size");
			return false;
		}
	}
	std::vector<MPI_Request> req(2 * neighbors.size());

	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer.data()), type.mpisize * sBuffer.size(), type.mpitype, neighbors[n], TAG::EX_KNOWN, group->communicator, req.data() + 2 * n);
	}

	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Irecv(const_cast<Ttype*>(rBuffer[n].data()), type.mpisize * rBuffer[n].size(), type.mpitype, neighbors[n], TAG::EX_KNOWN, group->communicator, req.data() + 2 * n + 1);
		rmin.min(type.mpisize * rBuffer[n].size());
		rmax.max(type.mpisize * rBuffer[n].size());
	}
	profiler::synccheckpoint("msg_prepared");

	MPI_Waitall(2 * neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::EX_KNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_exchange_known_size");
	return true;
}


template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_exchange_unknown_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	MPIType type(MPITools::getType<Ttype>());

	for (size_t n = 0; n < neighbors.size(); n++) {
		if (type.mpisize * sBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_exchange_unknown_size");
			return false;
		}
	}

	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::EX_UNKNOWN, group->communicator, req.data() + n);
		smin.min(type.mpisize * sBuffer[n].size());
		smax.max(type.mpisize * sBuffer[n].size());
	}
	profiler::synccheckpoint("msg_prepared");

	size_t counter = 0;
	while (counter < neighbors.size()) {
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, TAG::EX_UNKNOWN, group->communicator, &status);
		int count;
		MPI_Get_count(&status, type.mpitype, &count);
		rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.mpisize);
		MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.mpitype, status.MPI_SOURCE, TAG::EX_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
		counter++;
		profiler::checkpoint("receive");
		rmin.min(count * type.mpisize);
		rmax.max(count * type.mpisize);
	}
	profiler::synccheckpoint("receive_finish");

	MPI_Waitall(neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Probe(ANY_SOURCE) can be problem when calling this function more times
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::EX_UNKNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_exchange_unknown_size");
	return true;
}

template <typename Ttype>
bool Communication::exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_exchange_unknown_size");
	profiler::syncparam("neighbors", neighbors.size());
	profiler::syncparam("size", sBuffer.size());
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * sBuffer.size() > 1 << 30) {
		profiler::syncend("comm_exchange_unknown_size");
		return false;
	}

	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer.data()), type.mpisize * sBuffer.size(), type.mpitype, neighbors[n], TAG::EX_UNKNOWN, group->communicator, req.data() + n);
	}
	profiler::synccheckpoint("msg_prepared");

	size_t counter = 0;
	MPI_Status status;
	while (counter < neighbors.size()) {
		MPI_Probe(MPI_ANY_SOURCE, TAG::EX_UNKNOWN, group->communicator, &status);
		int count;
		MPI_Get_count(&status, type.mpitype, &count);
		rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.mpisize);
		MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.mpitype, status.MPI_SOURCE, TAG::EX_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
		counter++;
		profiler::checkpoint("receive");
		rmin.min(count * type.mpisize);
		rmax.max(count * type.mpisize);
	}
	profiler::synccheckpoint("receive_finish");

	MPI_Waitall(neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Probe(ANY_SOURCE) can be problem when calling this function more times
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::EX_UNKNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_exchange_unknown_size");
	return true;
}

template <typename Ttype, typename Talloc>
bool Communication::receiveLower(std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, MPIGroup *group)
{
	profiler::syncstart("comm_receive_lower");
	profiler::syncparam("size", sBuffer.size());
	rBuffer.resize(sBuffer.size());
	if (group->rank % 2 == 1) {
		MPI_Recv(rBuffer.data(), rBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank - 1, 0, group->communicator, MPI_STATUS_IGNORE);
		if (group->rank + 1 != group->size) {
			MPI_Send(sBuffer.data(), sBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank + 1, 0, group->communicator);
		}
	} else {
		if (group->rank + 1 != group->size) {
			MPI_Send(sBuffer.data(), sBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank + 1, 0, group->communicator);
		}
		if (group->rank) {
			MPI_Recv(rBuffer.data(), rBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank - 1, 0, group->communicator, MPI_STATUS_IGNORE);
		}
	}
	profiler::syncend("comm_receive_lower");
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_receive_lower_known_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());

	for (size_t n = 0; n < neighbors.size(); n++) {
		if (neighbors[n] > group->rank) {
			if (type.mpisize * sBuffer[n].size() > 1 << 30) {
				profiler::syncend("comm_receive_lower_known_size");
				return false;
			}
		} else {
			if (type.mpisize * rBuffer[n].size() > 1 << 30) {
				profiler::syncend("comm_receive_lower_known_size");
				return false;
			}
		}
	}

	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (neighbors[n] > group->rank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_LOW_KNOWN, group->communicator, req.data() + n);
			smin.min(type.mpisize * sBuffer[n].size());
			smax.max(type.mpisize * sBuffer[n].size());
		}
		if (neighbors[n] < group->rank) {
			MPI_Irecv(rBuffer[n].data(), type.mpisize * rBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_LOW_KNOWN, group->communicator, req.data() + n);
			rmin.min(type.mpisize * rBuffer[n].size());
			rmax.max(type.mpisize * rBuffer[n].size());
		}
	}
	profiler::synccheckpoint("msg_prepared");

	MPI_Waitall(neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::R_LOW_KNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_receive_lower_known_size");
	return true;
}

template <typename Ttype>
bool Communication::receiveLowerUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_receive_lower_unknown_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (type.mpisize * sBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_receive_lower_unknown_size");
			return false;
		}
	}

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (group->rank < neighbors[n]) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_LOW_UNKNOWN, group->communicator, req.data() + rSize++);
			smin.min(type.mpisize * sBuffer[n].size());
			smax.max(type.mpisize * sBuffer[n].size());
		}
	}
	profiler::synccheckpoint("msg_prepared");

	size_t counter = neighbors.end() - std::lower_bound(neighbors.begin(), neighbors.end(), group->rank);
	MPI_Status status;
	while (counter < neighbors.size()) {
		MPI_Probe(MPI_ANY_SOURCE, TAG::R_LOW_UNKNOWN, group->communicator, &status);
		int count;
		MPI_Get_count(&status, type.mpitype, &count);
		rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.mpisize);
		MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.mpitype, status.MPI_SOURCE, TAG::R_LOW_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
		counter++;
		profiler::checkpoint("receive");
		rmin.min(count * type.mpisize);
		rmax.max(count * type.mpisize);
	}
	profiler::synccheckpoint("receive_finish");

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::R_LOW_UNKNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_receive_lower_unknown_size");
	return true;
}

template <typename Ttype, typename Talloc>
bool Communication::receiveUpper(std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, const MPIGroup *group)
{
	profiler::syncstart("comm_receive_upper");
	profiler::syncparam("size", sBuffer.size());
	rBuffer.resize(sBuffer.size());
	if (group->rank % 2 == 1) {
		MPI_Send(sBuffer.data(), sBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank - 1, 0, group->communicator);
		if (group->rank + 1 != group->size) {
			MPI_Recv(rBuffer.data(), rBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank + 1, 0, group->communicator, MPI_STATUS_IGNORE);
		}
	} else {
		if (group->rank + 1 != group->size) {
			MPI_Recv(rBuffer.data(), rBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank + 1, 0, group->communicator, MPI_STATUS_IGNORE);
		}
		if (group->rank) {
			MPI_Send(sBuffer.data(), sBuffer.size() * sizeof(Ttype), MPI_BYTE, group->rank - 1, 0, group->communicator);
		}
	}
	profiler::syncend("comm_receive_upper");
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_receive_upper_known_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (neighbors[n] < group->rank) {
			if (type.mpisize * sBuffer[n].size() > 1 << 30) {
				return false;
			}
		} else {
			if (type.mpisize * rBuffer[n].size() > 1 << 30) {
				return false;
			}
		}
	}

	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size(); n++) {
		if (neighbors[n] < group->rank) {
			// bullxmpi violate MPI standard (cast away constness)
			MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_UP_KNOWN, group->communicator, req.data() + n);
			smin.min(type.mpisize * sBuffer[n].size());
			smax.max(type.mpisize * sBuffer[n].size());
		}
		if (neighbors[n] > group->rank) {
			MPI_Irecv(rBuffer[n].data(), type.mpisize * rBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_UP_KNOWN, group->communicator, req.data() + n);
			rmin.min(type.mpisize * rBuffer[n].size());
			rmax.max(type.mpisize * rBuffer[n].size());
		}
	}
	profiler::synccheckpoint("msg_prepared");

	MPI_Waitall(neighbors.size(), req.data(), MPI_STATUSES_IGNORE);
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::R_UP_KNOWN;
	}
	profiler::synccheckpoint("msg_waitall");
	profiler::syncend("comm_receive_upper_known_size");
	return true;
}

template <typename Ttype>
bool Communication::receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group)
{
	profiler::syncstart("comm_receive_upper_unknown_size");
	profiler::syncparam("neighbors", neighbors.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());
	for (size_t n = 0; n < neighbors.size() && neighbors[n] < group->rank; n++) {
		if (type.mpisize * sBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_receive_upper_unknown_size");
			return false;
		}
	}

	auto n2i = [ & ] (size_t neighbor) {
		return std::lower_bound(neighbors.begin(), neighbors.end(), neighbor) - neighbors.begin();
	};

	size_t rSize = 0;
	std::vector<MPI_Request> req(neighbors.size());
	for (size_t n = 0; n < neighbors.size() && neighbors[n] < group->rank; n++) {
		// bullxmpi violate MPI standard (cast away constness)
		MPI_Isend(const_cast<Ttype*>(sBuffer[n].data()), type.mpisize * sBuffer[n].size(), type.mpitype, neighbors[n], TAG::R_UP_UNKNOWN, group->communicator, req.data() + rSize++);
		smin.min(type.mpisize * sBuffer[n].size());
		smax.max(type.mpisize * sBuffer[n].size());
	}
	profiler::synccheckpoint("msg_prepared");

	size_t counter = std::lower_bound(neighbors.begin(), neighbors.end(), group->rank) - neighbors.begin();
	MPI_Status status;
	while (counter < neighbors.size()) {
		MPI_Probe(MPI_ANY_SOURCE, TAG::R_UP_UNKNOWN, group->communicator, &status);
		int count;
		MPI_Get_count(&status, type.mpitype, &count);
		rBuffer[n2i(status.MPI_SOURCE)].resize(count / type.mpisize);
		MPI_Recv(rBuffer[n2i(status.MPI_SOURCE)].data(), count, type.mpitype, status.MPI_SOURCE, TAG::R_UP_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
		counter++;
		profiler::checkpoint("receive");
		rmin.min(count * type.mpisize);
		rmax.max(count * type.mpisize);
	}
	profiler::synccheckpoint("receive_finish");

	MPI_Waitall(rSize, req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::R_UP_UNKNOWN;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_receive_upper_unknown_size");
	return true;
}

template <typename Ttype, typename Talloc>
bool Communication::gatherUnknownSize(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, MPIGroup *group)
{
	std::vector<size_t> offsets;
	return gatherUnknownSize(sBuffer, rBuffer, offsets, group);
}

template <typename Ttype, typename Talloc>
bool Communication::gatherUnknownSize(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, std::vector<size_t> &offsets, MPIGroup *group)
{
	profiler::syncstart("comm_gather_unknown_size");
	MPIType type(MPITools::getType<Ttype>());
	profiler::syncparam("size", sBuffer.size() * type.mpisize);

	if (Communication::manual) {
		int left = 0;
		int right = left + group->size;
		int rank = left + group->rank;

		MPIType otype = MPITools::getType<size_t>();
		size_t ssize = sBuffer.size(), rsize = 0;
		offsets.resize(group->size + 1);
		MPI_Allgather(&ssize, 1, otype.mpitype, offsets.data(), 1, otype.mpitype, group->communicator);
		profiler::synccheckpoint("allgather_size");
		for (int i = 0; i <= group->size; ++i) {
			size_t tmp = offsets[i];
			offsets[i] = rsize;
			rsize += tmp;
		}
		rsize -= offsets[group->rank];
		std::vector<Ttype, Talloc> recv;
		recv.reserve(rsize);
		recv.insert(recv.end(), sBuffer.begin(), sBuffer.end());

		RecursiveHalving rh(rank, left, right);
		int levels = std::ceil(std::log2(group->size));

		profiler::synccheckpoint("init_recursion");

		while (levels) {
			if (rh.recurse(levels--)) {
				if (rh.left == rank) {
					int recvsize;
					size_t offset = recv.size();
					MPI_Status status;
					MPI_Probe(rh.mid, TAG::GATHER_UNKNOWN, group->communicator, &status);
					MPI_Get_count(&status, type.mpitype, &recvsize);
					recv.resize(recv.size() + (recvsize / type.mpisize));
					MPI_Recv(recv.data() + offset, recvsize, type.mpitype, rh.mid, TAG::GATHER_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
					profiler::checkpoint("recv");
					profiler::param("size", recvsize * type.mpisize);
				}
				if (rh.mid == rank) {
					MPI_Send(recv.data(), recv.size() * type.mpisize, type.mpitype, rh.left, TAG::GATHER_UNKNOWN, group->communicator);
					profiler::checkpoint("send");
					profiler::param("size", recv.size() * type.mpisize);
					profiler::syncend("comm_gather_unknown_size");
					return true;
				}
			}
		}
		rBuffer.swap(recv);
		profiler::syncend("comm_gather_unknown_size");
		return true;
	}

	int size = type.mpisize * sBuffer.size();
	std::vector<int> rSizes(group->size), rOffsets(group->size);
	MPI_Gather(&size, 1, MPI_INT, rSizes.data(), 1, MPI_INT, 0, group->communicator);

	if (!group->rank) {
		size = 0;
		for (size_t i = 0; i < rSizes.size(); i++) {
			rOffsets[i] = size;
			size += rSizes[i];
		}
		rBuffer.resize(size / type.mpisize);
	}

	// bullxmpi violate MPI standard (cast away constness)
	MPI_Gatherv(const_cast<Ttype*>(sBuffer.data()), type.mpisize * sBuffer.size(), type.mpitype, rBuffer.data(), rSizes.data(), rOffsets.data(), type.mpitype, 0, group->communicator);

	offsets.resize(group->size + 1);
	for (size_t i = 0; i < rOffsets.size(); i++) {
		offsets[i] = rOffsets[i] / type.mpisize;
	}
	offsets.back() = (rOffsets.back() + rSizes.back()) / type.mpisize;
	profiler::syncend("comm_gather_unknown_size");
	return true;
}
template <typename Ttype>
bool Communication::allGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group)
{
	profiler::syncstart("comm_allgather_unknown_size");
	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * data.size() > 1 << 30) {
		profiler::syncend("comm_allgather_unknown_size");
		return false;
	}

	if (Communication::manual) {
		int left = 0;
		int right = left + group->size;
		int rank = left + group->rank;

		std::vector<Ttype> recv;
		RecursiveHalving rh(rank, left, right);
		int levels = std::ceil(std::log2(group->size));
		profiler::synccheckpoint("init_recursion");

		auto receive = [&] (int rank) {
			int recvsize;
			MPI_Status status;
			MPI_Probe(rank, TAG::ALLGATHER_UNKNOWN, group->communicator, &status);
			MPI_Get_count(&status, type.mpitype, &recvsize);
			recv.resize(recvsize / type.mpisize);
			MPI_Recv(recv.data(), recvsize, type.mpitype, rank, TAG::ALLGATHER_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
			profiler::checkpoint("recv");
			profiler::param("size", recvsize * type.mpisize);
		};

		while (levels) {
			if (rh.recurse(levels--)) {
				if (rh.islower()) {
					if (rh.ispaired()) {
						MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.twin, TAG::ALLGATHER_UNKNOWN, group->communicator);
						profiler::checkpoint("send");
						profiler::param("size", data.size() * type.mpisize);
						receive(rh.twin);
					} else {
						receive(rh.mid);
					}
					data.insert(data.end(), recv.begin(), recv.end());
					recv.clear();
				} else {
					receive(rh.twin);
					MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.twin, TAG::ALLGATHER_UNKNOWN, group->communicator);
					profiler::checkpoint("send");
					profiler::param("size", data.size() * type.mpisize);
					if (rh.treatodd()) {
						MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.mid - 1, TAG::ALLGATHER_UNKNOWN, group->communicator);
						profiler::checkpoint("send");
						profiler::param("size", data.size() * type.mpisize);
					}
					data.insert(data.begin(), recv.begin(), recv.end());
					recv.clear();
				}
				profiler::checkpoint("merge");
			}
		}
		profiler::syncend("comm_allgather_unknown_size");
		return true;
	}

	int size = type.mpisize * data.size();
	std::vector<int> rSizes(group->size), rOffsets(group->size);
	MPI_Allgather(&size, 1, MPI_INT, rSizes.data(), 1, MPI_INT, group->communicator);

	std::vector<Ttype> rdata;
	size = 0;
	for (size_t i = 0; i < rSizes.size(); i++) {
		rOffsets[i] = size;
		size += rSizes[i];
	}
	rdata.resize(size / type.mpisize);

	MPI_Allgatherv(data.data(), type.mpisize * data.size(), type.mpitype, rdata.data(), rSizes.data(), rOffsets.data(), type.mpitype, group->communicator);

	rdata.swap(data);
	profiler::syncend("comm_allgather_unknown_size");
	return true;
}

template <typename Ttype>
bool Communication::uniqueAllGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group)
{
	profiler::syncstart("comm_unique_allgather_unknown_size");
	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * data.size() > 1 << 30) {
		profiler::syncend("comm_allgather_unknown_size");
		return false;
	}

	int left = 0;
	int right = left + group->size;
	int rank = left + group->rank;

	std::vector<Ttype> recv, tmp;
	RecursiveHalving rh(rank, left, right);
	int levels = std::ceil(std::log2(group->size));
	profiler::synccheckpoint("init_recursion");

	auto receive = [&] (int rank) {
		int recvsize;
		MPI_Status status;
		MPI_Probe(rank, TAG::ALLGATHER_UNKNOWN, group->communicator, &status);
		MPI_Get_count(&status, type.mpitype, &recvsize);
		recv.resize(recvsize / type.mpisize);
		MPI_Recv(recv.data(), recvsize, type.mpitype, rank, TAG::ALLGATHER_UNKNOWN, group->communicator, MPI_STATUS_IGNORE);
		profiler::checkpoint("recv");
		profiler::param("size", recvsize * type.mpisize);
	};

	auto merge = [&] () {
		tmp.resize(data.size() + recv.size());
		auto end = std::set_union(data.begin(), data.end(), recv.begin(), recv.end(), tmp.begin());
		tmp.resize(end - tmp.begin());
		data.swap(tmp);
	};

	while (levels) {
		if (rh.recurse(levels--)) {
			if (rh.islower()) {
				if (rh.ispaired()) {
					MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.twin, TAG::ALLGATHER_UNKNOWN, group->communicator);
					profiler::checkpoint("send");
					profiler::param("size", data.size() * type.mpisize);
					receive(rh.twin);
				} else {
					receive(rh.mid);
				}
				merge();
				recv.clear();
			} else {
				receive(rh.twin);
				MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.twin, TAG::ALLGATHER_UNKNOWN, group->communicator);
				profiler::checkpoint("send");
				profiler::param("size", data.size() * type.mpisize);
				if (rh.treatodd()) {
					MPI_Send(data.data(), data.size() * type.mpisize, type.mpitype, rh.mid - 1, TAG::ALLGATHER_UNKNOWN, group->communicator);
					profiler::checkpoint("send");
					profiler::param("size", data.size() * type.mpisize);
				}
				merge();
				recv.clear();
			}
			profiler::checkpoint("merge");
		}
	}
	profiler::syncend("comm_allgather_unknown_size");

	profiler::syncend("comm_unique_allgather_unknown_size");
	return true;
}

template <typename Ttype>
bool Communication::broadcastUnknownSize(std::vector<Ttype> &buffer, MPIGroup *group)
{
	profiler::syncstart("comm_broadcast_unknown_size");
	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * buffer.size() > 1 << 30) {
		profiler::syncend("comm_broadcast_unknown_size");
		return false;
	}
	int size = buffer.size();
	MPI_Bcast(&size, 1, MPI_INT, 0, group->communicator);
	buffer.resize(size);
	MPI_Bcast(buffer.data(), type.mpisize * size, type.mpitype, 0, group->communicator);
	profiler::syncend("comm_broadcast_unknown_size");
	return true;
}

template <typename Ttype, typename Tdistribution>
bool Communication::balance(std::vector<Ttype> &buffer, const std::vector<Tdistribution> &currentDistribution, const std::vector<Tdistribution> &targetDistribution, MPIGroup *group)
{
	profiler::syncstart("comm_balance");
	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * buffer.size() > 1 << 30) {
		return false;
	}

	std::vector<Ttype> result(targetDistribution[group->rank + 1] - targetDistribution[group->rank]);
	std::vector<int> ssize(group->size), sdisp(group->size), rsize(group->size), rdisp(group->size);

	auto fill = [&] (
			const std::vector<Tdistribution> &from, const std::vector<Tdistribution> &to,
			std::vector<int> &size, std::vector<int> &disp) {

		Tdistribution offset = 0;
		Tdistribution restSize = from[group->rank + 1] - from[group->rank];
		Tdistribution tIndex = std::lower_bound(to.begin(), to.end(), from[group->rank] + 1) - to.begin() - 1;
		Tdistribution tOffset = from[group->rank] - to[tIndex];
		while (restSize) {
			if (restSize < to[tIndex + 1] - to[tIndex] - tOffset) {
				size[tIndex] = restSize;
				disp[tIndex] = offset;
				restSize = 0;
			} else {
				size[tIndex] = to[tIndex + 1] - to[tIndex] - tOffset;
				disp[tIndex] = offset;
				restSize -= size[tIndex];
				offset += size[tIndex];
				++tIndex;
			}
			tOffset = 0;
		}

//		for (int r = 0; r < group->size; ++r) {
//			size[r] *= sizeof(Ttype);
//			disp[r] *= sizeof(Ttype);
//		}
	};

	fill(currentDistribution, targetDistribution, ssize, sdisp);
	fill(targetDistribution, currentDistribution, rsize, rdisp);

	std::vector<MPI_Request> requests(2 * group->size);
	int nrequests = 0;

	for (int r = 0; r < group->size; ++r) {
		if (rsize[r]) {
			MPI_Irecv(result.data() + rdisp[r], type.mpisize * rsize[r], type.mpitype, r, TAG::BALANCE, group->communicator, requests.data() + nrequests++);
		}
	}

	for (int r = 0; r < group->size; ++r) {
		if (ssize[r]) {
			MPI_Isend(buffer.data() + sdisp[r], type.mpisize * ssize[r], type.mpitype, r, TAG::BALANCE, group->communicator, requests.data() + nrequests++);
		}
	}

	MPI_Waitall(nrequests, requests.data(), MPI_STATUSES_IGNORE);
	buffer.swap(result);
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::BALANCE;
	}
	profiler::syncend("comm_balance");
	return true;
}

template <typename Ttype>
bool Communication::allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize, MPIGroup *group)
{
	profiler::syncstart("comm_alltoallv");
	MPIType type(MPITools::getType<Ttype>());
	if (type.mpisize * sBuffer.size() > 1 << 30) {
		profiler::syncend("comm_alltoallv");
		return false;
	}
	std::vector<int> _ssize = ssize, _rsize = rsize;
	std::vector<int> sdisp(group->size), rdisp(group->size);
	for (int r = 0; r < group->size; r++) {
		_ssize[r] *= type.mpisize;
		_rsize[r] *= type.mpisize;
	}
	for (int r = 1; r < group->size; r++) {
		sdisp[r] = sdisp[r - 1] + _ssize[r - 1];
		rdisp[r] = rdisp[r - 1] + _rsize[r - 1];
	}
	MPI_Alltoallv(sBuffer.data(), _ssize.data(), sdisp.data(), type.mpitype, rBuffer.data(), _rsize.data(), rdisp.data(), type.mpitype, group->communicator);
	profiler::syncend("comm_alltoallv");
	return true;
}

template <typename Ttype>
Ttype Communication::exscan(Ttype &value, MPIGroup *group)
{
	profiler::syncstart("comm_exscan");
	MPIType type(MPITools::getType<Ttype>());
	Ttype size = value;
	if (group->size == 1) {
		value = 0;
		profiler::syncend("comm_exscan");
		return size;
	}

	if (type.isnative) {
		MPI_Exscan(&size, &value, 1, MPITools::getType<Ttype>().mpitype, MPI_SUM, group->communicator);
	} else {
		// TODO:
	}
	profiler::synccheckpoint("mpi_exscan");

	size = value + size;
	MPI_Bcast(&size, type.mpisize, type.mpitype, group->size - 1, group->communicator);
	if (group->rank == 0) {
		value = 0;
	}
	MPI_Barrier(group->communicator);
	profiler::synccheckpoint("get_size");
	profiler::syncend("comm_exscan");
	return size;
}

template <typename Ttype>
bool Communication::exscan(std::vector<Ttype> &sum, std::vector<Ttype> &offset, MPIGroup *group)
{
	profiler::syncstart("comm_exscan");
	profiler::syncparam("size", offset.size());
	MPIType type(MPITools::getType<Ttype>());
	if (!type.isnative) {
		profiler::syncend("comm_exscan");
		return false;
	}
	sum = offset;
	if (group->size == 1) {
		std::fill(offset.begin(), offset.end(), 0);
		profiler::syncend("comm_exscan");
		return true;
	}

	MPI_Exscan(sum.data(), offset.data(), sum.size(), MPITools::getType<Ttype>().mpitype, MPI_SUM, group->communicator);
	profiler::synccheckpoint("mpi_exscan");
	for (size_t i = 0; i < sum.size(); i++) {
		sum[i] += offset[i];
	}
	MPI_Bcast(sum.data(), sum.size(), type.mpitype, group->size - 1, group->communicator);
	if (group->rank == 0) {
		std::fill(offset.begin(), offset.end(), 0);
	}
	MPI_Barrier(group->communicator);
	profiler::synccheckpoint("get_size");
	profiler::syncend("comm_exscan");
	return true;
}

template <typename Ttype>
std::vector<Ttype> Communication::getDistribution(Ttype size, MPIGroup *group)
{
	profiler::syncstart("comm_get_distribution");
	MPIType type(MPITools::getType<Ttype>());
	std::vector<Ttype> result(group->size + 1);
	Ttype esize = size;
	result.back() = Communication::exscan(esize, group);

	profiler::synccheckpoint("exscan");
	MPI_Allgather(&esize, type.mpisize, type.mpitype, result.data(), type.mpisize, type.mpitype, group->communicator);
	profiler::synccheckpoint("allgather");
	profiler::syncend("comm_get_distribution");
	return result;
}

template <typename Ttype>
bool Communication::sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources, MPIGroup *group)
{
	profiler::syncstart("comm_send_various_target");
	profiler::syncparam("targets", targets.size());
	auto &smin = profiler::syncparam<size_t>("smin", 1 << 30);
	auto &smax = profiler::syncparam<size_t>("smax", 0);
	auto &rmin = profiler::syncparam<size_t>("rmin", 1 << 30);
	auto &rmax = profiler::syncparam<size_t>("rmax", 0);
	MPIType type(MPITools::getType<Ttype>());
	for (size_t n = 0; n < targets.size(); n++) {
		if (type.mpisize * sBuffer[n].size() > 1 << 30) {
			profiler::syncend("comm_send_various_target");
			return false;
		}
	}

	std::vector<int> smsgcounter(group->size);
	std::vector<int> rmsgcounter(group->size);
	for (size_t n = 0; n < targets.size(); n++) {
		smsgcounter[targets[n]] = 1;
	}

	MPI_Allreduce(smsgcounter.data(), rmsgcounter.data(), group->size, MPI_INT, MPI_SUM, group->communicator);
	profiler::synccheckpoint("msg_counter");

	std::vector<MPI_Request> req(targets.size());
	for (size_t t = 0; t < targets.size(); t++) {
		MPI_Isend(const_cast<Ttype*>(sBuffer[t].data()), type.mpisize * sBuffer[t].size(), type.mpitype, targets[t], TAG::SEND_VARIOUS, group->communicator, req.data() + t);
		smin.min(type.mpisize * sBuffer[t].size());
		smax.max(type.mpisize * sBuffer[t].size());
	}
	profiler::synccheckpoint("msg_prepared");

	int counter = 0;
	MPI_Status status;
	sources.clear();
	std::vector<std::vector<Ttype> > tmpBuffer;
	tmpBuffer.reserve(rmsgcounter[group->rank]);
	while (counter < rmsgcounter[group->rank]) {
		MPI_Probe(MPI_ANY_SOURCE, TAG::SEND_VARIOUS, group->communicator, &status);
		int count;
		MPI_Get_count(&status, type.mpitype, &count);
		tmpBuffer.push_back(std::vector<Ttype>(count / type.mpisize));
		MPI_Recv(tmpBuffer.back().data(), count, type.mpitype, status.MPI_SOURCE, TAG::SEND_VARIOUS, group->communicator, MPI_STATUS_IGNORE);
		sources.push_back(status.MPI_SOURCE);
		counter++;
		profiler::checkpoint("receive");
		rmin.min(count * type.mpisize);
		rmax.max(count * type.mpisize);
	}
	profiler::synccheckpoint("receive_finish");

	std::vector<int> permutation(sources.size());
	for (size_t i = 0; i < sources.size(); i++) {
		permutation[i] = i;
	}
	std::sort(permutation.begin(), permutation.end(), [&] (int i, int j) { return sources[i] < sources[j]; });
	rBuffer.resize(tmpBuffer.size());
	for (size_t i = 0; i < permutation.size(); i++) {
		rBuffer[i].swap(tmpBuffer[permutation[i]]);
	}

	std::sort(sources.begin(), sources.end());
	profiler::synccheckpoint("permute");

	MPI_Waitall(targets.size(), req.data(), MPI_STATUSES_IGNORE);
	MPI_Barrier(group->communicator); // MPI_Iprobe(ANY_SOURCE) can be problem when calling this function more times
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::SEND_VARIOUS;
	}
	profiler::synccheckpoint("waitall");
	profiler::syncend("comm_send_various_target");
	return true;
}

template <typename Ttype>
void Communication::allReduce(std::vector<Ttype> &data, OP op, int left, int right, MPIGroup *group)
{
	profiler::start("comm_allreduce");
	profiler::param("size", data.size());
	std::vector<Ttype> recv = data;
	MPIType type(MPITools::getType<Ttype>());

	auto merge = [&] () {
		for (size_t i = 0; i < data.size(); ++i) {
			switch (op) {
			case OP::SUM: data[i] += recv[i]; break;
			case OP::MAX: data[i] = std::max(data[i], recv[i]); break;
			case OP::MIN: data[i] = std::min(data[i], recv[i]); break;
			}
		}
	};

	RecursiveHalving rh(group->rank, left, right);
	while (rh.recurse()) {
		if (rh.islower()) {
			// LOWER half to UPPER half
			if (rh.isodd()) {
				MPI_Send(data.data(), type.mpisize * data.size(), type.mpitype, rh.mid, TAG::ALLREDUCE, group->communicator);
			}
			if (rh.ispaired()) {
				MPI_Send(data.data(), type.mpisize * data.size(), type.mpitype, rh.twin, TAG::ALLREDUCE, group->communicator);
				MPI_Recv(recv.data(), type.mpisize * data.size(), type.mpitype, rh.twin, TAG::ALLREDUCE, group->communicator, MPI_STATUS_IGNORE);
				merge();
			}
		} else {
			// UPPER half to LOWER half
			MPI_Recv(recv.data(), type.mpisize * data.size(), type.mpitype, rh.twin, TAG::ALLREDUCE, group->communicator, MPI_STATUS_IGNORE);
			MPI_Send(data.data(), type.mpisize * data.size(), type.mpitype, rh.twin, TAG::ALLREDUCE, group->communicator);
			merge();
			if (rh.treatodd()) {
				MPI_Recv(recv.data(), type.mpisize * data.size(), type.mpitype, rh.mid - 1, TAG::ALLREDUCE, group->communicator, MPI_STATUS_IGNORE);
				merge();
			}
		}
		rh.exchanged();
	}
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::ALLREDUCE;
	}
	profiler::end("comm_allreduce");
}

template <typename Ttype, typename Talloc>
void Communication::scatterv(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, std::vector<size_t> &displacement, MPIGroup *group)
{
	profiler::syncstart("comm_scatterv");
	profiler::syncparam("size", sBuffer.size());
	MPIType type(MPITools::getType<Ttype>());

	if (Communication::manual) {
		int left = 0;
		int right = left + group->size;
		int rank = left + group->rank;

		RecursiveHalving rh(rank, left, right);
		while (rh.recurse()) {
			if (rh.left == rank) {
				if (rank == left) {
					MPI_Send(sBuffer.data() + displacement[rh.mid - left] - displacement[rh.left - left], type.mpisize * (displacement[rh.right - left] - displacement[rh.mid - left]), type.mpitype, rh.mid, TAG::SCATTERV, group->communicator);
				} else {
					MPI_Send(rBuffer.data() + displacement[rh.mid - left] - displacement[rh.left - left], type.mpisize * (displacement[rh.right - left] - displacement[rh.mid - left]), type.mpitype, rh.mid, TAG::SCATTERV, group->communicator);
					rBuffer.resize(displacement[rh.mid - left] - displacement[rh.left - left]);
				}
			}
			if (rh.mid == rank) {
				rBuffer.resize(displacement[rh.right - left] - displacement[rh.mid - left]);
				MPI_Recv(rBuffer.data(), type.mpisize * (displacement[rh.right - left] - displacement[rh.mid - left]), type.mpitype, rh.left, TAG::SCATTERV, group->communicator, MPI_STATUS_IGNORE);
			}
			rh.exchanged();
		}
		if (rank == left) {
			rBuffer.clear();
			rBuffer.insert(rBuffer.end(), sBuffer.begin(), sBuffer.begin() + displacement[1]);
		}
	} else {
		std::vector<int> disp(displacement.begin(), displacement.end()), count(group->size);
		for (size_t i = 1; i <= count.size(); i++) {
			count[i - 1] = disp[i] - disp[i - 1];
		}
		rBuffer.resize(count[group->rank]);
		MPI_Scatterv(sBuffer.data(), count.data(), disp.data(), type.mpitype, rBuffer.data(), count[group->rank], type.mpitype, 0, group->communicator);
	}
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::SCATTERV;
	}
	profiler::syncend("comm_scatterv");
}

template <typename Ttype, typename Talloc>
bool Communication::scatter(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, size_t size, MPIGroup *group)
{
	if (size != (size_t)(int)size) {
		return false;
	}
	profiler::syncstart("comm_scatter");

	MPIType type(MPITools::getType<Ttype>());

	if (Communication::manual) {
		int left = 0;
		int right = left + group->size;
		int rank = left + group->rank;

		RecursiveHalving rh(rank, left, right);
		while (rh.recurse()) {
			if (rh.left == rank) {
				if (rank == left) {
					MPI_Send(sBuffer.data() + (rh.mid - rh.left) * size, type.mpisize * (rh.right - rh.mid) * size, type.mpitype, rh.mid, TAG::SCATTER, group->communicator);
				} else {
					MPI_Send(rBuffer.data() + (rh.mid - rh.left) * size, type.mpisize * (rh.right - rh.mid) * size, type.mpitype, rh.mid, TAG::SCATTER, group->communicator);
					rBuffer.resize((rh.mid - rh.left) * size);
				}
			}
			if (rh.mid == rank) {
				rBuffer.resize((rh.right - rh.mid) * size);
				MPI_Recv(rBuffer.data(), type.mpisize * (rh.right - rh.mid) * size, type.mpitype, rh.left, TAG::SCATTER, group->communicator, MPI_STATUS_IGNORE);
			}
			rh.exchanged();
		}
		if (rank == left) {
			rBuffer.clear();
			rBuffer.insert(rBuffer.end(), sBuffer.begin(), sBuffer.begin() + size);
		}
	} else {
		rBuffer.resize(size);
		MPI_Scatter(sBuffer.data(), type.mpisize * size, type.mpitype, rBuffer.data(), type.mpisize * size, type.mpitype, 0, group->communicator);
	}
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::SCATTER;
	}
	profiler::syncend("comm_scatter");
	return true;
}

template <typename Ttype, typename Talloc>
bool Communication::allToAllWithDataSizeAndTarget(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, int left, int right, MPIGroup *group)
{
	profiler::syncstart("comm_alltoall");
	profiler::syncparam("size", sBuffer.size());
	MPIType type(MPITools::getType<Ttype>());
	std::vector<Ttype, Talloc> prevsend, send, recv;
	recv.reserve(1.2 * sBuffer.size());
	rBuffer.reserve(recv.capacity());
	send.reserve(2.5 * sBuffer.size());
	prevsend.reserve(send.capacity());

	send = sBuffer;

	MPI_Status status;
	int recvsize, recvmidsize;

	auto movebefore = [&] (std::vector<Ttype, Talloc> &data, int rank, size_t begin, size_t end) {
		size_t pos = begin;
		while (pos < end && (int)data[pos + 1] < rank) {
			pos += data[pos];
		}
		return pos;
	};

	size_t mybegin = movebefore(send, group->rank, 0, send.size());
	size_t myend = movebefore(send, group->rank + 1, mybegin, send.size());
	rBuffer.clear();
	rBuffer.insert(rBuffer.end(), send.begin() + mybegin, send.begin() + myend);

	profiler::synccheckpoint("init");

	RecursiveHalving rh(group->rank, left, right);
	while (rh.recurse()) {
		if (rh.islower()) {
			// LOWER half to UPPER half
			if (rh.isodd()) {
				// PRE :
				// send: l1, l2, l3, ME, u1, u2, u3

				// POST:
				// SEND: ME -> u1
				// RECV:

				// send: l1, l2, l3
				size_t my = movebefore(send, group->rank, 0, send.size());
				size_t upper = movebefore(send, rh.mid, my, send.size());

				if ((send.size() - upper) * type.mpisize >= 1UL << 31) {
					profiler::syncend("comm_alltoall");
					return false;
				}
				MPI_Send(send.data() + upper, type.mpisize * (send.size() - upper), type.mpitype, rh.mid, TAG::ALL_TO_ALL_OPT, group->communicator);
				profiler::checkpoint("send");
				profiler::param("size", type.mpisize * (send.size() - upper));
				send.resize(my);
			} else {
				// PRE :
				// send: l1, l2(ME), l3, l4, u1, u2, u3


				// POST:
				// SEND: u1, u2, u3 -> u2;
				// RECV: l1, l2, l3, l4 <- u2

				// sBuffer: l1, l2, l3, l4, u1, u2, u3
				// send: l1, l1, l2, l2, l3, l3, l4, l4

				size_t upper = movebefore(send, rh.mid, 0, send.size());
				if ((send.size() - upper) * type.mpisize >= 1UL << 31) {
					profiler::syncend("comm_alltoall");
					return false;
				}
				MPI_Send(send.data() + upper, type.mpisize * (send.size() - upper), type.mpitype, rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator);
				profiler::checkpoint("send");
				profiler::param("size", type.mpisize * (send.size() - upper));

				MPI_Probe(rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator, &status);
				MPI_Get_count(&status, type.mpitype, &recvsize);
				recv.resize(recvsize / type.mpisize);
				MPI_Recv(recv.data(), recvsize, type.mpitype, rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator, MPI_STATUS_IGNORE);
				profiler::checkpoint("recv");
				profiler::param("size", recvsize * type.mpisize);

				send.swap(prevsend);
				send.clear();

				size_t recvbegin = 0;
				size_t recvend = recvbegin;
				size_t sendbegin = 0;
				size_t sendend = sendbegin;
				for (int r = rh.left; r < rh.mid; r++) {
					recvbegin = recvend;
					recvend = movebefore(recv, r + 1, recvbegin, recv.size());
					sendbegin = sendend;
					sendend = movebefore(prevsend, r + 1, sendbegin, prevsend.size());
					if (r == group->rank) {
						rBuffer.insert(rBuffer.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					} else {
						send.insert(send.end(), recv.begin() + recvbegin, recv.begin() + recvend);
						send.insert(send.end(), prevsend.begin() + sendbegin, prevsend.begin() + sendend);
					}
				}
				profiler::checkpoint("post_process");
			}
		} else {
			// UPPER half to LOWER half

			size_t upper = movebefore(send, rh.mid, 0, send.size());

			MPI_Probe(rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator, &status);
			MPI_Get_count(&status, type.mpitype, &recvsize);
			recv.resize(recvsize / type.mpisize);
			MPI_Recv(recv.data(), recvsize, type.mpitype, rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator, MPI_STATUS_IGNORE);
			profiler::checkpoint("recv");
			profiler::param("size", recvsize * type.mpisize);
			if (upper * type.mpisize >= 1UL << 31) {
				profiler::syncend("comm_alltoall");
				return false;
			}
			MPI_Send(send.data(), type.mpisize * upper, type.mpitype, rh.twin, TAG::ALL_TO_ALL_OPT, group->communicator);
			profiler::checkpoint("send");
			profiler::param("size", type.mpisize * upper);

			recvmidsize = recvsize / type.mpisize;
			if (rh.treatodd()) {
				// l1, l2, l3, l4, u1(ME), u2, u3
				// RECV: l4
				MPI_Probe(rh.mid - 1, TAG::ALL_TO_ALL_OPT, group->communicator, &status);
				MPI_Get_count(&status, type.mpitype, &recvsize);
				recv.resize(recv.size() + recvsize / type.mpisize);
				MPI_Recv(recv.data() + recvmidsize, recvsize, type.mpitype, rh.mid - 1, TAG::ALL_TO_ALL_OPT, group->communicator, MPI_STATUS_IGNORE);
				profiler::checkpoint("recv");
				profiler::param("size", recvsize * type.mpisize);
			}

			// PRE :
			// send: l1, l2, l3, l4, u1(ME), u2, u3

			// POST:
			// SEND: l1, l2, l3, l4 -> l1;
			// RECV: u1, u2, u3 <- l1
			// RECV: u1, u2, u3 <- l4 (recvsize > recvmidsize)

			// sBuffer: l1, l2, l3, l4, u1, u2, u3
			// send: l1, l1, l2, l2, l3, l3, l4, l4

			send.swap(prevsend);
			send.clear();

			size_t recvbegin = 0;
			size_t recvend = recvbegin;
			size_t recvmidbegin = recvmidsize;
			size_t recvmidend = recvmidbegin;
			size_t sendbegin = movebefore(prevsend, rh.mid, 0, prevsend.size());
			size_t sendend = sendbegin;
			for (int r = rh.mid; r < rh.right; r++) {
				recvbegin = recvend;
				recvend = movebefore(recv, r + 1, recvbegin, recvmidsize);
				recvmidbegin = recvmidend;
				recvmidend = movebefore(recv, r + 1, recvmidbegin, recv.size());
				sendbegin = sendend;
				sendend = movebefore(prevsend, r + 1, sendbegin, prevsend.size());
				if (r == group->rank) {
					rBuffer.insert(rBuffer.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					rBuffer.insert(rBuffer.end(), recv.begin() + recvmidbegin, recv.begin() + recvmidend);
				} else {
					send.insert(send.end(), recv.begin() + recvbegin, recv.begin() + recvend);
					send.insert(send.end(), recv.begin() + recvmidbegin, recv.begin() + recvmidend);
					send.insert(send.end(), prevsend.begin() + sendbegin, prevsend.begin() + sendend);
				}
			}
			profiler::checkpoint("post_process");
		}
		rh.exchanged();
	}
	if (group->communicator == MPI_COMM_WORLD) {
		++TAG::ALL_TO_ALL_OPT;
	}
	profiler::syncend("comm_alltoall");
	return true;
}

}



