
#include "distributedscanner.h"

using namespace espreso;

bool DistributedScanner::check(const char *c, const std::string &key)
{
	return memcmp(c, key.data(), key.size()) == 0;
}

void DistributedScanner::align(InputFile &input, const std::string &chars)
{
	auto isin = [&] (char c) {
		for (size_t i = 0; i < chars.size(); ++i) {
			if (c == chars[i]) {
				return true;
			}
		}
		return false;
	};
	if (input.distribution[info::mpi::rank] != input.distribution[info::mpi::rank + 1]) {
		if (input.distribution[info::mpi::rank] != 0) {
			while (!isin(*input.begin)) {
				++input.begin;
				++input.distribution[info::mpi::rank];
			}
			++input.begin;
			++input.distribution[info::mpi::rank];
		}

		if (input.distribution[info::mpi::rank + 1] != input.distribution.back()) {
			while (input.end < input.hardend && !isin(*input.end)) {
				++input.end;
				++input.distribution[info::mpi::rank + 1];
			}
			if (input.end < input.hardend) {
				++input.end;
				++input.distribution[info::mpi::rank + 1];
			}
		}
	}
}

void DistributedScanner::add(const char* key, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip)
{
	_filter[(unsigned char)key[0]] = 1;
	_keywords.push_back(Keyword{ (int)_keys.size(), (int)strlen(key), constructor, skip });
	_keys += key;
}

void DistributedScanner::add(std::initializer_list<const char*> keys, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip)
{
	for (auto it = keys.begin(); it != keys.end(); ++it) {
		add(*it, constructor, skip);
	}
}

void DistributedScanner::addEnd(const char* key, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip)
{
	_endFilter[(unsigned char)key[(int)strlen(key) - 1]] = 1;
	_endKeywords.push_back(Keyword{ (int)_keys.size(), (int)strlen(key), constructor, skip });
	_keys += key;
}

void DistributedScanner::scan(InputFile &input)
{
	profiler::syncstart("distributed_scan");
	profiler::syncparam("size", input.end - input.begin);
	int threads = info::env::OMP_NUM_THREADS;
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, input.end - input.begin);

	struct position { size_t c; size_t keyindex; };
	std::vector<std::vector<position> > found(threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		std::vector<position> tfound;
		for (size_t c = tdistribution[t]; c < tdistribution[t + 1]; ++c) {
			if (_filter[(unsigned char)*(input.begin + c)]) {
				for (size_t k = 0; k < _keywords.size(); ++k) {
					if (memcmp(input.begin + c, _keys.data() + _keywords[k].offset, _keywords[k].size) == 0) {
						tfound.push_back({ c, k });
						c += _keywords[k].skip(input.begin + c);
						break;
					}
				}
			}
		}
		found[t].swap(tfound);
	}

	for (int t = 0; t < threads; t++) {
		for (size_t i = 0; i < found[t].size(); ++i) {
			_keywords[found[t][i].keyindex].constructor(input.begin + found[t][i].c);
		}
	}
	profiler::syncend("distributed_scan");
}

void DistributedScanner::scanlines(InputFile &input)
{
	profiler::syncstart("distributed_scanlines");
	profiler::syncparam("size", input.end - input.begin);
	int threads = info::env::OMP_NUM_THREADS;
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, input.end - input.begin);

	struct position { size_t c; size_t keyindex; };
	std::vector<std::vector<position> > found(threads);
	size_t cc = 0;
	while (cc < tdistribution.back() && *(input.begin + cc++) != '\n');
	int endsize = 2 + ((cc < tdistribution.back() && *(input.begin + cc - 2) == '\r') ? 1 : 0);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		std::vector<position> tfound;

		size_t c = tdistribution[t];
		if (c && *(input.begin + c - 1) != '\n') {
			while (c < tdistribution[t + 1] && *(input.begin + c++) != '\n'); // start at new line
		}

		while (c < tdistribution[t + 1]) {
			if (_filter[(unsigned char)*(input.begin + c)]) {
				for (size_t i = 0; i < _keywords.size(); ++i) {
					if (memcmp(input.begin + c, _keys.data() + _keywords[i].offset, _keywords[i].size) == 0) {
						tfound.push_back({ c, i });
						c += _keywords[i].skip(input.begin + c);
						break;
					}
				}
			}
			while (input.begin + c < input.hardend && *(input.begin + c++) != '\n'); // to the next line start
			if (input.begin + c - endsize < input.hardend && _endFilter[(unsigned char)*(input.begin + c - endsize)]) {
				for (size_t i = 0; i < _endKeywords.size(); ++i) {
					if (memcmp(input.begin + c - endsize - _endKeywords[i].size + 1, _keys.data() + _endKeywords[i].offset, _endKeywords[i].size) == 0) {
						c -= endsize + _endKeywords[i].size - 1;
						while (c && *(input.begin + c - 1) != '\n') { --c; } // back to line start
						tfound.push_back({ c, i + _keywords.size() });
						c += _endKeywords[i].skip(input.begin + c);
						while (c < tdistribution[t + 1] && *(input.begin + c++) != '\n');
						break;
					}
				}
			}
		}

		found[t].swap(tfound);
	}

	for (int t = 0; t < threads; t++) {
		for (size_t i = 0; i < found[t].size(); ++i) {
			if (found[t][i].keyindex < _keywords.size()) {
				_keywords[found[t][i].keyindex].constructor(input.begin + found[t][i].c);
			} else {
				_endKeywords[found[t][i].keyindex - _keywords.size()].constructor(input.begin + found[t][i].c);
			}
		}
	}
	profiler::syncend("distributed_scanlines");
}
