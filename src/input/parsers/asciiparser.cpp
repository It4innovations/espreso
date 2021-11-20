
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

static inline void _pushWithString(std::vector<esint> &data, const char *begin, const char* c, size_t &size)
{
	char* next;
	data.push_back(strtol(c, &next, 10));
	if (data.back() == 0 && c == next) {
		data.back() = -(c - begin);
		while (!ASCIIParser::isempty(next)) { ++next; }
	}
	size = next - c;
}

static inline void _push(std::vector<float> &data, const char* c, size_t &size)
{
	char* next;
	data.push_back(strtof(c, &next));
	size = next - c;
}

static inline void _push(std::vector<double> &data, const char* c, size_t &size)
{
	char* next;
	data.push_back(strtod(c, &next));
	size = next - c;
}

static inline void _push(std::vector<_Point<esfloat> > &data, const char* c, size_t &size)
{
	char* next;
	const char* _c = c;
	double x = strtod(c, &next); c = next;
	double y = strtod(c, &next); c = next;
	double z = strtod(c, &next); c = next;
	data.push_back(Point(x, y, z));
	size = c - _c;
}

static inline void _next(InputFile &file, size_t &begin)
{
	while (file.begin + begin < file.hardend && ASCIIParser::isempty(file.begin + begin)) {
		++begin;
	}
}

void setDistribution(InputFile &file, size_t begin, size_t end, std::vector<size_t> &tdistribution, bool toLines = false)
{
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
			if (toLines) {
				while (*(file.begin + tdistribution[t]++) != '\n');
			} else {
				while (!ASCIIParser::isempty(file.begin + tdistribution[t]++));
			}
		}
	} else {
		tdistribution.resize(info::env::OMP_NUM_THREADS + 1, 0);
	}
}

template <typename TType>
void ASCIIParser::_parse(std::vector<TType> &data, InputFile &file, size_t begin, size_t end)
{
	profiler::syncstart("parse_ascii");
	std::vector<size_t> tdistribution;
	setDistribution(file, begin, end, tdistribution);

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

void ASCIIParser::parse(std::vector<esint> &ids, std::vector<_Point<esfloat> > &coordinates, InputFile &file, size_t begin, size_t end)
{
	profiler::syncstart("parse_ascii_mixed");
	std::vector<size_t> tdistribution;
	setDistribution(file, begin, end, tdistribution, true);

	std::vector<std::vector<esint> > tids(info::env::OMP_NUM_THREADS);
	std::vector<std::vector<_Point<esfloat> > > tcoords(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<esint> _ids;
		std::vector<_Point<esfloat> > _coords;
		_ids.reserve((tdistribution[t + 1] - tdistribution[t]) / 2);
		_coords.reserve((tdistribution[t + 1] - tdistribution[t]) / 2);
		size_t begin = tdistribution[t], size;
		for (_next(file, begin); begin < tdistribution[t + 1]; _next(file, begin)) {
			_push(_ids, file.begin + begin, size); begin += size;
			_push(_coords, file.begin + begin, size); begin += size;
		}
		tids[t].swap(_ids);
		tcoords[t].swap(_coords);
	}

	size_t isize = 0, dsize = 0;
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		isize += tids[t].size();
	}
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		dsize += tcoords[t].size();
	}
	ids.reserve(isize + 50);
	coordinates.reserve(dsize + 50);
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		ids.insert(ids.end(), tids[t].begin(), tids[t].end());
	}
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		coordinates.insert(coordinates.end(), tcoords[t].begin(), tcoords[t].end());
	}
	profiler::syncend("parse_ascii_mixed");
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

void ASCIIParser::parseWithStrings(std::vector<esint> &data, InputFile &file, size_t begin, size_t end)
{
	profiler::syncstart("parse_with_string");
	std::vector<size_t> tdistribution;
	setDistribution(file, begin, end, tdistribution);

	std::vector<std::vector<esint> > tdata(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		std::vector<esint> _data;
		_data.reserve((tdistribution[t + 1] - tdistribution[t]) / 2);
		size_t begin = tdistribution[t], size;
		for (_next(file, begin); begin < tdistribution[t + 1]; _next(file, begin)) {
			_pushWithString(_data, file.begin, file.begin + begin, size);
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
	profiler::syncend("parse_with_string");
}

void ASCIIParser::parse(std::vector<esint> &data, InputFile &file, size_t begin, size_t end)
{
	_parse(data, file, begin, end);
}

void ASCIIParser::parse(std::vector<float> &data, InputFile &file, size_t begin, size_t end)
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

void ASCIIParser::addmore(std::vector<float> &data, InputFile &file, size_t n, size_t end)
{
	_addmore(data, file, n, end);
}

void ASCIIParser::addmore(std::vector<double> &data, InputFile &file, size_t n, size_t end)
{
	_addmore(data, file, n, end);
}
