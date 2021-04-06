
#include "utils.h"

#include <algorithm>

namespace espreso {
namespace utils {

template<typename Ttype>
Ttype sizesToOffsets(Ttype *begin, Ttype *end, Ttype offset)
{
	Ttype sum = offset;
	for (auto it = begin; it != end; ++it) {
		Ttype tmp = *it;
		*it = sum;
		sum += tmp;
	}
	return sum;
}

template<typename Ttype>
std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes, const std::vector<Ttype> &offsets)
{
	std::vector<Ttype> sum(offsets);
	for (size_t i = 0; i < sizes.front().size(); i++) {
		for (size_t t = 0; t < sizes.size(); t++) {
			Ttype tmp = sizes[t][i];
			sizes[t][i] = sum[i];
			sum[i] += tmp;
		}
	}
	return sum;
}

template<typename Ttype, typename Tpermutation>
void permute(std::vector<Ttype> &data, const std::vector<Tpermutation> &permutation, size_t elementsize)
{
	std::vector<Ttype> _data(data.size());
	_data.swap(data);
	for (size_t e = 0; e < elementsize; e++) {
		for (size_t i = 0; i < data.size() / elementsize; i++) {
			data[elementsize * i + e] = _data[elementsize * permutation[i] + e];
		}
	}
}

template<typename Ttype>
void threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution)
{
	std::vector<size_t> offsets;
	for (size_t t = 0; t < distribution.size(); t++) {
		offsets.push_back(distribution[t].size() ? distribution[t].back() : 0);
	}
	sizesToOffsets(offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < distribution.size(); t++) {
		size_t offset = offsets[t];
		for (size_t i = 0; i < distribution[t].size(); i++) {
			distribution[t][i] += offset;
		}
	}
}

template<typename Ttype>
void threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	size_t threads = distribution.size() - 1;
	std::vector<size_t> offsets(distribution.size());
	for (size_t t = 0; t < threads; t++) {
		if (distribution[t] != distribution[t + 1]) {
			offsets[t] = data[distribution[t + 1]];
		}
	}
	sizesToOffsets(offsets);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		size_t offset = offsets[t];
		for (size_t i = distribution[t] + 1; i < distribution[t + 1] + 1; i++) {
			data[i] += offset;
		}
	}
}

template<typename Ttype>
void removeDuplicates(std::vector<Ttype> &data, size_t begin)
{
	if (data.size() == begin) {
		return;
	}
	size_t unique = begin;
	for (size_t d = begin + 1; d < data.size(); d++) {
		if (data[unique] != data[d]) {
			data[++unique] = data[d];
		}
	}

	data.resize(unique + 1);
}

template<typename Ttype>
void sortAndRemoveDuplicates(std::vector<Ttype> &data, size_t begin)
{
	std::sort(data.begin() + begin, data.end());
	removeDuplicates(data, begin);
}

template<typename Ttype>
void sortAndRemoveDuplicates(std::vector<std::vector<Ttype> > &data)
{
	for (size_t n = 0; n < data.size(); n++) {
		sortAndRemoveDuplicates(data[n]);
	}
}

template<typename Ttype>
void mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data)
{
	for (size_t t = 1; t < data.size(); t++) {
		data[0].insert(data[0].end(), data[t].begin(), data[t].end());
	}
	sortAndRemoveDuplicates(data[0]);
}

template<typename Ttype>
void mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data)
{
	#pragma omp parallel for
	for (size_t n = 0; n < data[0].size(); n++) {
		for (size_t t = 1; t < data.size(); t++) {
			data[0][n].insert(data[0][n].end(), data[t][n].begin(), data[t][n].end());
		}
		sortAndRemoveDuplicates(data[0][n]);
	}
}

template<typename Ttype>
void inplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	size_t size = distribution.size() - 1;

	size_t align = 1;
	while (align < distribution.size()) align = align << 1;

	std::vector<size_t> _distribution = distribution;
	_distribution.insert(_distribution.end(), align - size, _distribution.back());

	for (size_t i = 2; i <= align; i *= 2) {
		#pragma omp parallel for
		for (size_t t = 0; t < align / i; t++) {
			std::inplace_merge(
					data.data() + _distribution[i * t],
					data.data() + _distribution[i * t + i / 2],
					data.data() + _distribution[i * t + i]);
		}
	}
}

template<typename Ttype>
void sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	size_t size = distribution.size() - 1;

	#pragma omp parallel for
	for (size_t t = 0; t < size; t++) {
		std::sort(
				data.data() + distribution[t],
				data.data() + distribution[t + 1]);
	}

	inplaceMerge(data, distribution);
}

template<typename Ttype>
void inplaceMerge(std::vector<std::vector<Ttype> > &data)
{
	std::vector<size_t> distribution = { 0, data[0].size() };
	for (size_t t = 1; t < data.size(); t++) {
		data[0].insert(data[0].end(), data[t].begin(), data[t].end());
		distribution.push_back(data[0].size());
	}
	inplaceMerge(data[0], distribution);
}

template<typename Ttype>
void sortWithInplaceMerge(std::vector<std::vector<Ttype> > &data)
{
	std::vector<size_t> distribution = { 0, data[0].size() };
	for (size_t t = 1; t < data.size(); t++) {
		data[0].insert(data[0].end(), data[t].begin(), data[t].end());
		distribution.push_back(data[0].size());
	}
	sortWithInplaceMerge(data[0], distribution);
}

template<typename Ttype>
void sortWithUniqueMerge(std::vector<std::vector<Ttype> > &data)
{
	#pragma omp parallel for
	for (size_t t = 0; t < data.size(); t++) {
		sortAndRemoveDuplicates(data[t]);
	}
	std::vector<size_t> distribution = { 0, data[0].size() };
	for (size_t t = 1; t < data.size(); t++) {
		data[0].insert(data[0].end(), data[t].begin(), data[t].end());
		distribution.push_back(data[0].size());
	}
	inplaceMerge(data[0], distribution);
	sortAndRemoveDuplicates(data[0]);
}

template<typename Ttype>
void mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution)
{
	std::vector<size_t> _distribution(distribution.begin() + 1, distribution.end());

	size_t align = 1;
	while (align < _distribution.size()) align = align << 1;

	_distribution.insert(_distribution.end(), align - distribution.size() + 2, _distribution.back());

	for (size_t i = 2; i <= align; i *= 2) {
		#pragma omp parallel for
		for (size_t t = 0; t < align / i; t++) {
			std::inplace_merge(
					data.data() + _distribution[i * t],
					data.data() + _distribution[i * t + i / 2],
					data.data() + _distribution[i * t + i]);
		}
	}

	std::inplace_merge(
			data.data(),
			data.data() + _distribution.front(),
			data.data() + _distribution.back());
}

}
}
