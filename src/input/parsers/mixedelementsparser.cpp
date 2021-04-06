
#include "mixedelementsparser.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"

using namespace espreso;

void MixedElementsParser::add(esint *data, size_t size)
{
	_elements.push_back(ElementData{data, size});
}

void MixedElementsParser::parse(std::function<esint(size_t, esint)> enodes)
{
	profiler::syncstart("parse_mixed_elements");
	std::vector<int> check(_elements.size(), info::mpi::size), me(_elements.size()), last(_elements.size());
	std::vector<esint> esize(_elements.size()), lmissing(_elements.size()), dummy(_elements.size());
	first.resize(_elements.size());
	missing.resize(_elements.size());
	invalid.resize(_elements.size());
	nelements.resize(_elements.size());

	for (size_t i = 0; i < _elements.size(); ++i) {
		recognize(i, esize[i], first[i], missing[i], enodes);
	}

	profiler::synccheckpoint("recognize");

	offset = esize;
	Communication::exscan(dummy, offset);
	Communication::receiveLower(missing, lmissing);
	for (size_t i = 0; i < _elements.size(); ++i) {
		if (first[i] != lmissing[i]) {
			check[i] = info::mpi::rank;
		}
		if (_elements[i].size) {
			me[i] = info::mpi::rank;
		}
	}
	Communication::allReduce(check.data(), invalid.data(), _elements.size(), MPI_INT, MPI_MIN);
	Communication::allReduce(me.data(), last.data(), _elements.size(), MPI_INT, MPI_MAX);
	profiler::synccheckpoint("synchronize");

	fix(lmissing, esize, last, enodes);
	profiler::synccheckpoint("fix");
	for (size_t i = 0; i < _elements.size(); ++i) {
		nelements[i] = offset[i] + esize[i];
	}
	Communication::broadcast(nelements.data(), nelements.size(), MPITools::getType<esint>().mpitype, info::mpi::size - 1);
	profiler::syncend("parse_mixed_elements");
}

void MixedElementsParser::recognize(size_t index, esint &esize, esint &coffset, esint &missing, std::function<esint(size_t, esint)> enodes)
{
	auto tryrecognize = [&] (size_t i) {
		esize = 0;
		esint n = 1;
		while (n && i < _elements[index].size) {
			if ((n = enodes(index, _elements[index].data[i]))) {
				i += 1 + n;
				++esize;
			}
		}
		return i;
	};

	// try to recognize the first cell type description
	for (; coffset < 30 && (size_t)coffset < _elements[index].size; ++coffset) {
		size_t lastrecognized = tryrecognize(coffset);
		if (_elements[index].size <= lastrecognized) {
			missing = lastrecognized - _elements[index].size;
			break;
		}
	}
}

void MixedElementsParser::fix(std::vector<esint> &lmissing, std::vector<esint> &esize, std::vector<int> &last, std::function<esint(size_t, esint)> enodes)
{
	for (size_t i = 0; i < invalid.size(); ++i) {
		if (invalid[i] != info::mpi::size) {
			if (info::mpi::rank <= last[i]) {
				// fix the first cell type description position
				if (invalid[i] < info::mpi::rank) {
					MPI_Recv(&lmissing[i], 1, MPITools::getType(lmissing[i]).mpitype, info::mpi::rank - 1, 0, info::mpi::comm, MPI_STATUS_IGNORE);
					MPI_Recv(&offset[i], 1, MPITools::getType(offset[i]).mpitype, info::mpi::rank - 1, 0, info::mpi::comm, MPI_STATUS_IGNORE);
				}
				if (invalid[i] <= info::mpi::rank) {
					esint _esize = offset[i] + esize[i];
					if (first[i] != lmissing[i]) {
						first[i] = lmissing[i];
						recognize(i, _esize, first[i], missing[i], enodes);
						_esize += offset[i];
					}

					if (info::mpi::rank < last[i]) {
						// in the case of [invalid == info::mpi::rank] it initiates the synchronization
						MPI_Send(&missing[i], 1, MPITools::getType(missing[i]).mpitype, info::mpi::rank + 1, 0, info::mpi::comm);
						MPI_Send(&_esize, 1, MPITools::getType(_esize).mpitype, info::mpi::rank + 1, 0, info::mpi::comm);
					}
				}
			}
		}
	}
}
