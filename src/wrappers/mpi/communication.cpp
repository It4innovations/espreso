
#include "communication.h"

#include "basis/containers/tarray.h"
#include "basis/sfc/spacefillingcurve.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/parser.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/statisticsstore.h"

#include <string>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace espreso;

#define __GAP__ 1000

MPIOperations* MPITools::operations = NULL;
MPIGroup* MPITools::procs = NULL;
MPIGroup* MPITools::node = NULL;
MPIGroup* MPITools::instances = NULL;
MPIGroup* MPITools::global = NULL;
MPIGroup* MPITools::asynchronous = NULL;

MPISubset* MPITools::subset = NULL;
MPISubset* MPITools::singleton = NULL;

int Communication::TAG::SFC               =  0 * __GAP__;
int Communication::TAG::EX_KNOWN          =  1 * __GAP__;
int Communication::TAG::EX_UNKNOWN        =  2 * __GAP__;
int Communication::TAG::R_LOW_KNOWN       =  3 * __GAP__;
int Communication::TAG::R_LOW_UNKNOWN     =  4 * __GAP__;
int Communication::TAG::R_UP_KNOWN        =  5 * __GAP__;
int Communication::TAG::R_UP_UNKNOWN      =  6 * __GAP__;
int Communication::TAG::GATHER_UNKNOWN    =  7 * __GAP__;
int Communication::TAG::ALLGATHER_UNKNOWN =  8 * __GAP__;
int Communication::TAG::BCAST_UNKNOWN     =  9 * __GAP__;
int Communication::TAG::BALANCE           = 10 * __GAP__;
int Communication::TAG::ALL_TO_ALLV       = 11 * __GAP__;
int Communication::TAG::EXSCAN            = 12 * __GAP__;
int Communication::TAG::DISTRIBUTION      = 13 * __GAP__;
int Communication::TAG::SEND_VARIOUS      = 14 * __GAP__;
int Communication::TAG::ALL_TO_ALL_OPT    = 15 * __GAP__;
int Communication::TAG::SPLITTERS         = 16 * __GAP__;
int Communication::TAG::ALLREDUCE         = 17 * __GAP__;
int Communication::TAG::SCATTERV          = 18 * __GAP__;
int Communication::TAG::SCATTER           = 19 * __GAP__;

template<typename Ttype>
static void _scan(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) += *(static_cast<Ttype*>(in) + i);
	}
}

static void _mergeStatistics(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		(static_cast<Statistics*>(out) + i)->min = std::min((static_cast<Statistics*>(out) + i)->min, (static_cast<Statistics*>(in) + i)->min);
		(static_cast<Statistics*>(out) + i)->max = std::max((static_cast<Statistics*>(out) + i)->max, (static_cast<Statistics*>(in) + i)->max);
		(static_cast<Statistics*>(out) + i)->avg += (static_cast<Statistics*>(in) + i)->avg;
		(static_cast<Statistics*>(out) + i)->norm += (static_cast<Statistics*>(in) + i)->norm;
		(static_cast<Statistics*>(out) + i)->absmin = std::min(std::fabs((static_cast<Statistics*>(out) + i)->absmin), std::fabs((static_cast<Statistics*>(in) + i)->absmin));
		(static_cast<Statistics*>(out) + i)->absmax = std::max(std::fabs((static_cast<Statistics*>(out) + i)->absmax), std::fabs((static_cast<Statistics*>(in) + i)->absmax));
	}
}

template<typename Ttype>
static void _max(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) = std::max(*(static_cast<Ttype*>(in) + i), *(static_cast<Ttype*>(out) + i));
	}
}

template<typename Ttype>
static void _min(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) = std::min(*(static_cast<Ttype*>(in) + i), *(static_cast<Ttype*>(out) + i));
	}
}

template<typename Ttype>
static void _sum(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	int size = *len / sizeof(Ttype);
	for (int i = 0; i < size; i++) {
		*(static_cast<Ttype*>(out) + i) += *(static_cast<Ttype*>(in) + i);
	}
}

MPIOperations::MPIOperations()
{
	MPI_Op_create(_mergeStatistics, 1, &mergeStatistics);
	MPI_Type_contiguous(sizeof(Statistics), MPI_BYTE, &STATISTICS);
	MPI_Type_commit(&STATISTICS);
}

MPIOperations::~MPIOperations()
{
	MPI_Op_free(&mergeStatistics);
	MPI_Type_free(&STATISTICS);
}

MPIGroup::MPIGroup()
{
	MPI_Comm_dup(info::mpi::comm, &communicator);
	rank = info::mpi::rank;
	size = info::mpi::size;
}

MPIGroup::MPIGroup(MPI_Comm &comm)
{
	MPI_Comm_dup(comm, &communicator);
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);
}

MPIGroup::MPIGroup(MPI_Comm &&comm)
: communicator(std::move(comm))
{
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);
}

void MPIGroup::filter(const std::vector<int> &ranks)
{
	profiler::start("mpi_group_filter");
	MPI_Group commgroup, subgroup;
	MPI_Comm_group(communicator, &commgroup);
	MPI_Comm_free(&communicator);
	communicator = MPI_COMM_NULL;
	MPI_Group_incl(commgroup, ranks.size(), ranks.data(), &subgroup);
	if (std::binary_search(ranks.begin(), ranks.end(), rank)) {
		MPI_Comm_create_group(info::mpi::comm, subgroup, 0, &communicator);
	}
	MPI_Group_free(&commgroup);
	MPI_Group_free(&subgroup);
	if (communicator != MPI_COMM_NULL) {
		MPI_Comm_rank(communicator, &rank);
		MPI_Comm_size(communicator, &size);
	}
	profiler::end("mpi_group_filter");
}

void MPIGroup::split(int color, int key)
{
	profiler::syncstart("mpi_group_split");
	MPI_Comm oldcomm;
	MPI_Comm_dup(communicator, &oldcomm);
	MPI_Comm_free(&communicator);
	MPI_Comm_split(oldcomm, color, key, &communicator);
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &size);
	profiler::syncend("mpi_group_split");
}

MPIGroup::~MPIGroup()
{
	if (communicator != MPI_COMM_NULL) {
		MPI_Comm_free(&communicator);
	}
}

MPISubset::MPISubset(int max_mpi_procs)
{
	profiler::syncstart("mpi_subset");
	int subsize = 1;
	while (max_mpi_procs < (info::mpi::size / subsize) + (info::mpi::size % subsize ? 1 : 0)) {
		subsize += 1;
	}
	int wcolor = info::mpi::rank / subsize;
	int acolor = info::mpi::rank % subsize;

	withinsize = subsize;
	acrosssize = info::mpi::size / subsize + ((info::mpi::size % subsize) ? 1 : 0);
	MPI_Comm_split(info::mpi::comm, wcolor, info::mpi::rank, &within.communicator);
	MPI_Comm_rank(within.communicator, &within.rank);
	MPI_Comm_size(within.communicator, &within.size);

	MPI_Comm_split(info::mpi::comm, acolor, info::mpi::rank, &across.communicator);
	MPI_Comm_rank(across.communicator, &across.rank);
	MPI_Comm_size(across.communicator, &across.size);
	profiler::syncend("mpi_subset");
}

void MPITools::init()
{
	operations = new MPIOperations();
	procs = new MPIGroup();
	instances = new MPIGroup(info::mpi::icomm);
	global = new MPIGroup(info::mpi::gcomm);
	asynchronous = new MPIGroup(info::mpi::comm);

	MPI_Comm nodeComm;
	MPI_Comm_split_type(info::mpi::comm, MPI_COMM_TYPE_SHARED, info::mpi::rank, MPI_INFO_NULL, &nodeComm);
	node = new MPIGroup(nodeComm);

	subset = new MPISubset(info::mpi::size);
	singleton = new MPISubset(1);
}

void MPITools::setSubset(int subsetSize)
{
	if (subset) { delete subset; }
	if (singleton) { delete singleton; }

	subset = new MPISubset(subsetSize);
	singleton = new MPISubset(1);
}

void MPITools::reinit()
{
	if (procs) { delete procs; }
	if (instances) { delete instances; }
	if (global) { delete global; }
	if (asynchronous) { delete asynchronous; }

	procs = new MPIGroup();
	instances = new MPIGroup(info::mpi::icomm);
	global = new MPIGroup(info::mpi::gcomm);
	asynchronous = new MPIGroup(info::mpi::comm);

	MPI_Comm nodeComm;
	MPI_Comm_split_type(info::mpi::comm, MPI_COMM_TYPE_SHARED, info::mpi::rank, MPI_INFO_NULL, &nodeComm);
	node = new MPIGroup(nodeComm);
}

void MPITools::finish()
{
	if (operations) {
		delete operations;
	}
	if (instances) {
		delete instances;
	}
	if (global) {
		delete global;
	}
	if (asynchronous) {
		delete asynchronous;
	}
	if (procs) {
		delete procs;
	}

	if (subset) {
		delete subset;
	}
	if (singleton) {
		delete singleton;
	}
}

bool Communication::RecursiveHalving::recurse()
{
	mid = left + (right - left) / 2 + (right - left) % 2;
	if (islower()) {
		twin = rank + (mid - left);
	} else {
		twin = rank - (mid - left);
	}
	return left + 1 < right;
}

bool Communication::RecursiveHalving::recurse(int level)
{
	int depth = 1;
	left = min, right = max;
	while (recurse() && depth++ < level) {
		exchanged();
	}
	return left + 1 < right;
}

bool Communication::RecursiveHalving::islower()
{
	return rank < mid;
}

bool Communication::RecursiveHalving::isupper()
{
	return !islower();
}

bool Communication::RecursiveHalving::isodd()
{
	return ((right - left) % 2) && rank + 1 == mid;
}

bool Communication::RecursiveHalving::ispaired()
{
	return !isodd();
}

bool Communication::RecursiveHalving::treatodd()
{
	return ((right - left) % 2) && rank == mid;
}

void Communication::RecursiveHalving::exchanged()
{
	if (rank < mid) {
		right = mid;
	} else {
		left = mid;
	}
}

bool Communication::computeSFCBalancedBorders(SpaceFillingCurve &sfc, std::vector<esint> &sfcbuckets, std::vector<esint> &permutation, std::vector<esint> &sfcborders)
{
	eslog::start("SFC BORDERS", "SFC");

	esint esize = sfcbuckets.size();
	esize = Communication::exscan(esize);
	std::vector<esint> targetDistribution = tarray<esint>::distribute(info::mpi::size, esize);

	double PRECISION = 0.001 * std::log2(info::mpi::size);
	if (PRECISION * (targetDistribution.back() / info::mpi::size) < 2) {
		PRECISION = 2.01 / (targetDistribution.back() / info::mpi::size);
	}

	// allowed difference to the perfect distribution
	int ETOLERANCE = PRECISION * targetDistribution.back() / info::mpi::size;

	eslog::param("ETOLERANCE", ETOLERANCE);
	eslog::ln();

	esint DEPTH = 2;
	while (sfc.buckets(DEPTH++) < (size_t)info::mpi::size);

	size_t buckets = sfc.buckets(DEPTH);
	esint bstep = sfc.buckets(sfc.depth()) / buckets;

	std::vector<std::vector<esint> > bucketSum(DEPTH);
	for (esint d = 0; d < DEPTH; d++) {
		bucketSum[d].resize(sfc.buckets(d + 1) + 1);
	}

	for (auto e = sfcbuckets.begin(); e != sfcbuckets.end(); ++e) {
		++bucketSum.back()[*e / bstep];
	}
	utils::sizesToOffsets(bucketSum.back());

	for (esint d = DEPTH - 2; d >= 0; --d) {
		for (size_t b = 0; b < bucketSum[d].size(); ++b) {
			bucketSum[d][b] = bucketSum[d + 1][sfc.bucketSize() * b];
		}
	}

	std::vector<esint> scounts, rcounts;
	sfcborders.resize(info::mpi::size + 1, sfc.buckets(sfc.depth()));
	sfcborders.front() = 0;

	esint LEVEL = 0;
	esint bsize = sfc.bucketSize();
	esint coarsenig = sfc.buckets(sfc.depth());

	eslog::checkpointln("BUCKET PREPROCESSING");
	do {
		coarsenig /= bsize;
		scounts.resize(sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		rcounts.resize(sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		for (size_t b = 0, index = 0; b < sfc.sfcRefined(LEVEL).size(); b++) {
			for (esint i = 0; i <= bsize; i++, index++) {
				scounts[index] = bucketSum[LEVEL][bsize * sfc.sfcRefined(LEVEL)[b] + i];
			}
		}

		eslog::checkpointln("SCOUNT");

		MPI_Allreduce(scounts.data(), rcounts.data(), scounts.size(), MPITools::getType<esint>().mpitype, MPI_SUM, info::mpi::comm);

		eslog::checkpoint("ALLREDUCE");
		eslog::param("COUNT", scounts.size());
		eslog::ln();

		sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < sfc.sfcRefined(LEVEL).size(); b++) {
			esint boffset = b * (bsize + 1);
			esint rbegin = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset]) - targetDistribution.begin();
			esint rend = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset + bsize]) - targetDistribution.begin();

			for (esint r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < targetDistribution[r] && targetDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				sfcborders[r] = coarsenig * (bsize * sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > targetDistribution[r] && rcounts[boffset + i] - targetDistribution[r] > ETOLERANCE) {
					sfc.recurce(bsize * sfc.sfcRefined(LEVEL)[b] + i - 1);
				}
			}
		}
		sfc.finishLevel(LEVEL + 1);

		eslog::checkpointln("REFINE");

	} while (++LEVEL < DEPTH && sfc.hasLevel(LEVEL));

	eslog::checkpointln("COARSE LEVELS");

	// Go deeper if needed
	scounts.resize(sfc.sfcRefined(LEVEL).size() * (bsize + 1));
	std::fill(scounts.begin(), scounts.end(), 0);
	for (size_t b = 0, index = 0; b < sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
		scounts[index] = bucketSum[LEVEL - 1][sfc.sfcRefined(LEVEL)[b]];
	}
	rcounts.resize(scounts.size());
	std::vector<esint> refinedindices;

	eslog::checkpointln("PREPARE DEEPER LEVELS");

	while (LEVEL < SFCDEPTH && sfc.hasLevel(LEVEL)) {
		coarsenig /= bsize;
		bstep /= buckets;
		for (size_t b = 0, index = 0; b < sfc.sfcRefined(LEVEL).size(); b++, index++) {
			auto e = std::lower_bound(permutation.begin(), permutation.end(), coarsenig * bsize * sfc.sfcRefined(LEVEL)[b], [&] (esint i, esint bound) { return sfcbuckets[i] < bound; });
			for (esint i = 0; i < bsize; i++, index++) {
				while (e != permutation.end() && sfcbuckets[*e] < coarsenig * bsize * (esint)sfc.sfcRefined(LEVEL)[b] + (i + 1) * coarsenig) {
					++scounts[index + 1];
					++e;
				}
				scounts[index + 1] += scounts[index];
			}
		}

		eslog::checkpointln("DEEPER SCOUNT");

		MPI_Allreduce(scounts.data(), rcounts.data(), scounts.size(), MPITools::getType<esint>().mpitype, MPI_SUM, info::mpi::comm);

		eslog::checkpoint("DEEPER ALLREDUCE");
		eslog::param("DEEPER COUNT", scounts.size());
		eslog::ln();

		sfc.setLevel(LEVEL + 1);

		for (size_t b = 0; b < sfc.sfcRefined(LEVEL).size(); b++) {
			esint boffset = b * (bsize + 1);
			esint rbegin = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset]) - targetDistribution.begin();
			esint rend = std::lower_bound(targetDistribution.begin(), targetDistribution.end(), rcounts[boffset + bsize]) - targetDistribution.begin();

			for (esint r = rbegin, i = 0; r < rend; r++) {
				while (i <= bsize && rcounts[boffset + i] < targetDistribution[r] && targetDistribution[r] - rcounts[boffset + i] >= ETOLERANCE) {
					++i;
				}
				sfcborders[r] = coarsenig * (bsize * sfc.sfcRefined(LEVEL)[b] + i);
				if (rcounts[boffset + i] > targetDistribution[r] && rcounts[boffset + i] - targetDistribution[r] > ETOLERANCE) {
					sfc.recurce(bsize * sfc.sfcRefined(LEVEL)[b] + i - 1);
					refinedindices.push_back(boffset + i - 1);
				}
			}
		}
		sfc.finishLevel(++LEVEL);
		utils::sortAndRemoveDuplicates(refinedindices);

		rcounts.swap(scounts);
		scounts.resize(sfc.sfcRefined(LEVEL).size() * (bsize + 1));
		std::fill(scounts.begin(), scounts.end(), 0);
		for (size_t b = 0, index = 0; b < sfc.sfcRefined(LEVEL).size(); b++, index += bsize + 1) {
			scounts[index] = rcounts[refinedindices[b]];
		}
		rcounts.resize(scounts.size());
		refinedindices.clear();

		eslog::checkpointln("DEEPER REFINE");
	}
	eslog::endln("FINISH");
	return true;
}

static void getMedian(esint *keys, esint targetsum, esint tolerance, esint begin, esint end, esint lowersum, esint uppersum, esint lowerkey, esint upperkey, esint &localmid, esint &globalmid, esint &keymid, int left, int right)
{
	profiler::start("get_median");
	profiler::param("keys", end - begin);
	const esint maxdepth = 3, maxbuckets = 31;
	esint buckets = maxbuckets;
	std::vector<esint> shistogram(maxbuckets + 1), rhistogram(maxbuckets + 1);
	esint depth = 0, initlowersum = lowersum;

	profiler::checkpoint("init");
	while (lowersum + tolerance < targetsum && targetsum < uppersum - tolerance && depth++ < maxdepth) {
		if (buckets > upperkey - lowerkey) {
			buckets = upperkey - lowerkey;
		}
		esint bsize = (upperkey - lowerkey) / buckets + ((upperkey - lowerkey) % buckets ? 1 : 0);
		std::fill(shistogram.data(), shistogram.data() + buckets + 1, 0);
		for (auto key = keys + begin; key != keys + end; ++key) {
			++shistogram[(*key - lowerkey) / bsize];
		}
		utils::sizesToOffsets(shistogram.data(), shistogram.data() + buckets + 1, begin);

		memcpy(rhistogram.data(), shistogram.data(), sizeof(esint) * (buckets + 1));
		profiler::checkpoint("compute_histogram");
		Communication::allReduce(rhistogram, Communication::OP::SUM, left, right);
		profiler::checkpoint("allreduce");

		auto it  = std::lower_bound(rhistogram.data(), rhistogram.data() + buckets + 1, targetsum - initlowersum);
		begin    = shistogram[it - rhistogram.data() - 1];
		end      = shistogram[it - rhistogram.data()];
		uppersum = initlowersum + *(it);
		lowersum = initlowersum + *(it - 1);
		upperkey = lowerkey + bsize * (it - rhistogram.data());
		lowerkey = lowerkey + bsize * (it - rhistogram.data() - 1);
		profiler::checkpoint("set_splitter");
	}

	if (targetsum - lowersum < uppersum - targetsum) {
		localmid = begin;
		globalmid = lowersum;
		keymid = lowerkey;
	} else {
		localmid = end;
		globalmid = uppersum;
		keymid = upperkey;
	}
	profiler::end("get_median");
}

bool Communication::computeSplitters(std::vector<esint> &keys, std::vector<esint> &permutation, std::vector<esint> &splitters, MPIGroup *group)
{
	profiler::syncstart("compute_splitters");
	profiler::syncparam("keys", keys.size());
	splitters.resize(group->size + 1);
	MPIType type = MPITools::getType<esint>();

	esint size, mysize = keys.size();
	esint max, mymax = 0;
	if (keys.size()) {
		mymax = keys[permutation.back()];
	}
	MPI_Allreduce(&mymax, &max, 1, MPITools::getType<esint>().mpitype, MPI_MAX, group->communicator);
	MPI_Allreduce(&mysize, &size, 1, MPITools::getType<esint>().mpitype, MPI_SUM, group->communicator);
	splitters.back() = max + 1;
	profiler::synccheckpoint("reduce_stats");

	std::vector<esint> targetDistribution = tarray<esint>::distribute(group->size, size), distribution(group->size + 1);
	distribution.back() = size;
	esint tolerance = 0.25 * std::ceil(size / (double)group->size);

	std::vector<esint> prevsend, send, recv, recvmid, recvmerge;
	send.reserve(keys.size());
	recv.reserve(keys.size());

	for (auto i = permutation.begin(); i != permutation.end(); ++i) {
		send.push_back(keys[*i]);
	}

	profiler::synccheckpoint("init");

	auto receive = [&] (int rank) {
		int recvsize;
		MPI_Status status;
		MPI_Probe(rank, TAG::SPLITTERS, group->communicator, &status);
		MPI_Get_count(&status, type.mpitype, &recvsize);
		recv.resize(recvsize / type.mpisize);
		MPI_Recv(recv.data(), recvsize, type.mpitype, rank, TAG::SPLITTERS, group->communicator, MPI_STATUS_IGNORE);
		profiler::checkpoint("recv");
		profiler::param("recvsize", recvsize);
	};

	RecursiveHalving rh(group->rank, 0, group->size);
	while (rh.recurse()) {
		esint localmid, globalmid, keymid;
		getMedian(send.data(), targetDistribution[rh.mid], tolerance, 0, send.size(), distribution[rh.left], distribution[rh.right], splitters[rh.left], splitters[rh.right], localmid, globalmid, keymid, rh.left, rh.right);
		splitters[rh.mid] = keymid;
		distribution[rh.mid] = globalmid;

		profiler::checkpoint("get_splitter");

		if (rh.islower()) {
			// LOWER half to UPPER half
			if (type.mpisize * (send.size() - localmid) > 1 << 30) {
				return false;
			}

			if (rh.isodd()) {
				MPI_Send(send.data() + localmid, type.mpisize * (send.size() - localmid), type.mpitype, rh.mid, TAG::SPLITTERS, group->communicator);
				profiler::checkpoint("send");
				profiler::param("sendsize", type.mpisize * (send.size() - localmid));
				recv.clear();
			}
			if (rh.ispaired()) {
				MPI_Send(send.data() + localmid, type.mpisize * (send.size() - localmid), type.mpitype, rh.twin, TAG::SPLITTERS, group->communicator);
				profiler::checkpoint("send");
				profiler::param("sendsize", type.mpisize * (send.size() - localmid));
				receive(rh.twin);
			}

			prevsend.swap(send);
			send.resize(localmid + recv.size());
			std::merge(prevsend.begin(), prevsend.begin() + localmid, recv.begin(), recv.end(), send.begin());
			profiler::checkpoint("merge");
		} else {
			// UPPER half to LOWER half
			if (type.mpisize * localmid > 1 << 30) {
				return false;
			}

			receive(rh.twin);
			MPI_Send(send.data(), type.mpisize * localmid, type.mpitype, rh.twin, TAG::SPLITTERS, group->communicator);
			profiler::checkpoint("send");
			profiler::param("sendsize", type.mpisize * localmid);

			if (rh.treatodd()) {
				recvmid.swap(recv);
				receive(rh.mid - 1);
				recvmerge.resize(recv.size() + recvmid.size());
				std::merge(recv.begin(), recv.end(), recvmid.begin(), recvmid.end(), recvmerge.begin());
				recv.swap(recvmerge);
				profiler::checkpoint("merge");
			}

			prevsend.swap(send);
			send.resize(prevsend.size() - localmid + recv.size());
			std::merge(prevsend.begin() + localmid, prevsend.end(), recv.begin(), recv.end(), send.begin());
			profiler::checkpoint("merge");
		}
		rh.exchanged();
	}

	profiler::synccheckpoint("rh_finish");

	esint begin = send.size();
	Communication::exscan(begin, group);
	esint end = begin + send.size();

	std::fill(splitters.begin(), splitters.end(), 0);
	for (auto it = std::lower_bound(targetDistribution.begin() + 1, targetDistribution.end(), begin); *it < end; ++it) { // skip the first value that is 0
		splitters[it - targetDistribution.begin()] = send[*it - begin];
	}
	profiler::synccheckpoint("set_splitters");

	Communication::allReduce(splitters, Communication::OP::MAX, 0, group->size);
	splitters.front() = 0;
	for (auto it = splitters.rbegin(); it != splitters.rend(); ++it) {
		if (*it == 0) {
			*it = max + 1;
		} else {
			break;
		}
	}

	profiler::syncend("compute_splitters");
	return true;
}

bool Communication::barrier(MPIGroup *group)
{
	profiler::syncstart("mpi_barrier");
	MPI_Barrier(group->communicator);
	profiler::syncend("mpi_barrier");
	return true;
}

bool Communication::broadcast(void *data, size_t size, MPI_Datatype type, int root, MPIGroup *group)
{
	if (size != (size_t)(int)size) {
		return false;
	}
	profiler::syncstart("mpi_bcast");
	profiler::syncparam("size", size);
	MPI_Bcast(data, size, type, root, group->communicator);
	profiler::syncend("mpi_bcast");
	return true;
}

bool Communication::reduce(void *in, void *out, size_t size, MPI_Datatype type, MPI_Op op, int root, MPIGroup *group)
{
	if (size != (size_t)(int)size) {
		return false;
	}
	profiler::syncstart("mpi_reduce");
	profiler::syncparam("size", size);
	if (out == NULL) {
		if (info::mpi::rank == root) {
			MPI_Reduce(MPI_IN_PLACE, in, size, type, op, root, group->communicator);
		} else {
			MPI_Reduce(in, NULL, size, type, op, root, group->communicator);
		}
	} else {
		MPI_Reduce(in, out, size, type, op, root, group->communicator);
	}
	profiler::syncend("mpi_reduce");
	return true;
}

bool Communication::allReduce(void *in, void *out, size_t size, MPI_Datatype type, MPI_Op op, MPIGroup *group)
{
	if (size != (size_t)(int)size) {
		return false;
	}
	profiler::syncstart("mpi_allreduce");
	profiler::syncparam("size", size);
	if (out == NULL) {
		MPI_Allreduce(MPI_IN_PLACE, in, size, type, op, group->communicator);
	} else {
		MPI_Allreduce(in, out, size, type, op, group->communicator);
	}
	profiler::syncend("mpi_allreduce");
	return true;
}

bool Communication::allGather(void *in, void *out, size_t size, MPI_Datatype type, MPIGroup *group)
{
	if (size != (size_t)(int)size) {
		return false;
	}
	profiler::syncstart("mpi_allgather");
	profiler::syncparam("size", size);
	MPI_Allgather(in, size, type, out, size, type, group->communicator);
	profiler::syncend("mpi_allgather");
	return true;
}

void Communication::serialize(std::function<void(void)> fnc, MPIGroup *group)
{
	for (int r = 0; r < group->size; ++r) {
		if (r == group->rank) {
			fnc();
			utils::callusleep(500);
		}
		MPI_Barrier(group->communicator);
	}
	MPI_Barrier(group->communicator);
}


