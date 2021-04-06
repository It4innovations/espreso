
#include "asciiparser.h"
#include "basis/containers/tarray.h"
#include "basis/io/inputfile.h"
#include "basis/logging/profiler.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

using namespace espreso;

int ASCIIParser::keyend(const char *c)
{
	int i = 0;
	while (isempty(c + i)) { ++i; }
	while (!isempty(c + i)) { ++i; }
	return i;
}

static inline void _push(std::vector<esint> &data, const char* c, size_t &size)
{
	char* next;
	data.push_back(strtol(c, &next, 10));
	size = next - c;
}

static inline void _push(std::vector<double> &data, const char* c, size_t &size)
{
	char* next;
	data.push_back(strtod(c, &next));
	size = next - c;
}

static inline void _next(InputFile &file, size_t &begin)
{
	while (file.begin + begin < file.hardend && ASCIIParser::isempty(file.begin + begin)) {
		++begin;
	}
}

template <typename TType>
void ASCIIParser::_parse(std::vector<TType> &data, InputFile &file, size_t begin, size_t end)
{
	profiler::syncstart("parse_ascii");
	std::vector<size_t> tdistribution;
	size_t mpibegin = std::max(begin, file.distribution[info::mpi::rank]);
	size_t mpiend = std::min(end, file.distribution[info::mpi::rank + 1]);

	if (mpibegin < mpiend) {
		mpibegin -= file.distribution[info::mpi::rank];
		mpiend   -= file.distribution[info::mpi::rank];
		profiler::syncparam("size", mpiend - mpibegin);
		tdistribution = tarray<size_t>::distribute(info::env::OMP_NUM_THREADS, mpiend - mpibegin);
		for (size_t t = 0; t < tdistribution.size(); ++t) {
			tdistribution[t] += mpibegin;
		}
		for (size_t t = 1; t + 1 < tdistribution.size(); ++t) {
			while (!ASCIIParser::isempty(file.begin + tdistribution[t]++));
		}
	} else {
		tdistribution.resize(info::env::OMP_NUM_THREADS + 1, 0);
	}

	std::vector<std::vector<TType> > tdata(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<TType> _data;
		_data.reserve((tdistribution[t + 1] - tdistribution[t]) / 2);
		size_t begin = tdistribution[t], size;
		for (_next(file, begin); begin < tdistribution[t + 1]; _next(file, begin)) {
			_push(_data, file.begin + begin, size);
			begin += size;
		}
		tdata[t].swap(_data);
	}

	size_t dsize = 0;
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		dsize += tdata[t].size();
	}
	data.reserve(dsize + 50);
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		data.insert(data.end(), tdata[t].begin(), tdata[t].end());
	}
	profiler::syncend("parse_ascii");
}

template <typename TType>
void ASCIIParser::_addmore(std::vector<TType> &data, InputFile &file, size_t n, size_t end)
{
	if (n == 0) {
		return;
	}
	size_t size = 0, begin = std::min(end, file.distribution[info::mpi::rank + 1]) - file.distribution[info::mpi::rank];
	for (size_t i = 0; i < n; ++i) {
		_next(file, begin);
		_push(data, file.begin + begin, size);
		begin += size;
	}
}

void ASCIIParser::parse(std::vector<esint> &data, InputFile &file, size_t begin, size_t end)
{
	_parse(data, file, begin, end);
}

void ASCIIParser::parse(std::vector<double> &data, InputFile &file, size_t begin, size_t end)
{
	_parse(data, file, begin, end);
}

void ASCIIParser::addmore(std::vector<esint> &data, InputFile &file, size_t n, size_t end)
{
	_addmore(data, file, n, end);
}

void ASCIIParser::addmore(std::vector<double> &data, InputFile &file, size_t n, size_t end)
{
	_addmore(data, file, n, end);
}
