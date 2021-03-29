
#ifndef SRC_BASIS_UTILITIES_COMMUNICATION_H_
#define SRC_BASIS_UTILITIES_COMMUNICATION_H_

#include "esinfo/mpiinfo.h"

#include <cstddef>
#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct ProcessesReduction;
class SpaceFillingCurve;

class MPIOperations {
	friend class MPITools;
public:
	MPI_Op mergeStatistics;
	MPI_Datatype STATISTICS;

private:
	MPIOperations();
	~MPIOperations();
public:
	MPIOperations(MPIOperations const&) = delete;
	void operator=(MPIOperations const&) = delete;
};

struct MPIType {
	MPI_Datatype mpitype;
	int mpisize;
	bool isnative;

	MPIType(MPI_Datatype type): mpitype(type), mpisize(1), isnative(true) {}
	MPIType(MPI_Datatype type, int mpisize, bool isnative): mpitype(type), mpisize(mpisize), isnative(isnative) {}
};

struct MPIGroup {
	MPI_Comm communicator;
	int rank, size;

	MPIGroup();
	MPIGroup(MPI_Comm &comm);
	MPIGroup(MPI_Comm &&comm);
	~MPIGroup();

	void filter(const std::vector<int> &ranks);
	void split(int color, int key);
};

class MPISubset {
	friend class MPITools;

public:
	int withinsize, acrosssize;
	MPIGroup within, across;

	MPISubset(int max_mpi_procs);
private:
	MPISubset(MPISubset const&) = delete;
	void operator=(MPISubset const&) = delete;
};

class MPITools
{

public:
	static MPIOperations *operations;
	static MPIGroup *procs;
	static MPIGroup *node;
	static MPIGroup *instances;
	static MPIGroup *global;
	static MPIGroup *asynchronous;

	static MPISubset *subset;
	static MPISubset *singleton;

	template <typename Ttype>
	static MPIType getType();
	template <typename Ttype>
	static MPIType getType(const Ttype &value)
	{
		return getType<Ttype>();
	}

	static void init();
	static void setSubset(int subsetSize);
	static void reinit();
	static void finish();

private:
	MPITools() = delete;
};

struct Communication {

	const static bool manual = true;

	struct RecursiveHalving {
		int rank;
		int left, mid, right;
		int twin; // recv, send rank
		int min, max;

		RecursiveHalving(int rank, int left, int right)
		: rank(rank), left(left), right(right), min(left), max(right) { recurse(); }

		bool recurse();
		bool recurse(int level);
		bool islower();
		bool isupper();
		bool ispaired();
		bool isodd();
		bool treatodd();
		void exchanged();
	};

	enum class OP {
		SUM,
		MAX,
		MIN
	};

	struct TAG {
		static int
		SFC,
		EX_KNOWN, EX_UNKNOWN,
		R_LOW_KNOWN, R_LOW_UNKNOWN,
		R_UP_KNOWN, R_UP_UNKNOWN,
		GATHER_UNKNOWN, ALLGATHER_UNKNOWN, BCAST_UNKNOWN,
		BALANCE, ALL_TO_ALLV,
		EXSCAN, DISTRIBUTION,
		SEND_VARIOUS, ALL_TO_ALL_OPT,
		SPLITTERS, ALLREDUCE, SCATTERV, SCATTER;
	};

	static bool computeSFCBalancedBorders(SpaceFillingCurve &sfc, std::vector<esint> &sfcbuckets, std::vector<esint> &permutation, std::vector<esint> &sfcborders);
	static bool computeSplitters(std::vector<esint> &keys, std::vector<esint> &permutation, std::vector<esint> &splitters, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool allToAllWithDataSizeAndTarget(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, int left = 0, int right = MPITools::procs->size, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exchangeKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exchangeUnknownSize(const std::vector<Ttype> &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool receiveLower(std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveLowerKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveLowerUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool receiveUpper(std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveUpperKnownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool receiveUpperUnknownSize(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &neighbors, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool gatherUnknownSize(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool gatherUnknownSize(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, std::vector<size_t> &offsets, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool allGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool uniqueAllGatherUnknownSize(std::vector<Ttype> &data, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool broadcastUnknownSize(std::vector<Ttype> &buffer, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Tdistribution>
	static bool balance(std::vector<Ttype> &buffer, const std::vector<Tdistribution> &currentDistribution, const std::vector<Tdistribution> &targetDistribution, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool allToAllV(const std::vector<Ttype> &sBuffer, std::vector<Ttype> &rBuffer, const std::vector<int> &ssize, const std::vector<int> &rsize, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static Ttype exscan(Ttype &value, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool exscan(std::vector<Ttype> &sum, std::vector<Ttype> &offset, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static std::vector<Ttype> getDistribution(Ttype size, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, MPIGroup *group = MPITools::procs)
	{
		std::vector<int> sources;
		return sendVariousTargets(sBuffer, rBuffer, targets, sources, group);
	}

	template <typename Ttype>
	static bool sendVariousTargets(const std::vector<std::vector<Ttype> > &sBuffer, std::vector<std::vector<Ttype> > &rBuffer, const std::vector<int> &targets, std::vector<int> &sources, MPIGroup *group = MPITools::procs);

	template <typename Ttype>
	static void allReduce(std::vector<Ttype> &data, OP op, int left = 0, int right = MPITools::procs->size, MPIGroup *group = MPITools::procs);

	static bool barrier(MPIGroup *group = MPITools::procs);
	static bool broadcast(void *data, size_t size, MPI_Datatype type, int root, MPIGroup *group = MPITools::procs);
	static bool reduce(void *in, void *out, size_t size, MPI_Datatype type, MPI_Op op, int root, MPIGroup *group = MPITools::procs);
	static bool allReduce(void *in, void *out, size_t size, MPI_Datatype type, MPI_Op op, MPIGroup *group = MPITools::procs);
	static bool allGather(void *in, void *out, size_t size, MPI_Datatype type, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static void scatterv(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, std::vector<size_t> &displacement, MPIGroup *group = MPITools::procs);

	template <typename Ttype, typename Talloc=std::allocator<Ttype> >
	static bool scatter(const std::vector<Ttype, Talloc> &sBuffer, std::vector<Ttype, Talloc> &rBuffer, size_t size, MPIGroup *group = MPITools::procs);

	static void serialize(std::function<void(void)> fnc, MPIGroup *group = MPITools::procs);
};


}

#include "communication.hpp"




#endif /* SRC_BASIS_UTILITIES_COMMUNICATION_H_ */
