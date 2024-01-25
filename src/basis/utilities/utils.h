
#ifndef BASIS_UTILITIES_UTILS_H_
#define BASIS_UTILITIES_UTILS_H_

#include <cstddef>
#include <vector>
#include <limits>

namespace espreso {
namespace utils {

	template<typename Ttype>
	Ttype sizesToOffsets(Ttype *begin, Ttype *end, Ttype offset = 0);

	template<typename Ttype>
	Ttype sizesToOffsets(std::vector<Ttype> &sizes, Ttype offset = 0)
	{
		return sizesToOffsets(sizes.data(), sizes.data() + sizes.size(), offset);
	}

	template<typename Ttype>
	std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes, const std::vector<Ttype> &offsets);

	template<typename Ttype>
	std::vector<Ttype> sizesToOffsets(std::vector<std::vector<Ttype> > &sizes)
	{
		return sizesToOffsets(sizes, std::vector<Ttype>(sizes.front().size()));
	}

	template<typename Ttype, typename Tpermutation>
	void permute(std::vector<Ttype> &data, const std::vector<Tpermutation> &permutation, size_t elementsize = 1);

	template<typename Ttype>
	void threadDistributionToFullDistribution(std::vector<std::vector<Ttype> > &distribution);

	template<typename Ttype>
	void threadDistributionToFullDistribution(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void removeDuplicates(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	void sortAndRemoveDuplicates(std::vector<Ttype> &data, size_t begin = 0);

	template<typename Ttype>
	void sortAndRemoveDuplicates(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeThreadedUniqueData(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeThreadedUniqueData(std::vector<std::vector<std::vector<Ttype> > > &data);

	template<typename Ttype>
	void inplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void inplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void sortWithInplaceMerge(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename Ttype>
	void sortWithInplaceMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void sortWithUniqueMerge(std::vector<std::vector<Ttype> > &data);

	template<typename Ttype>
	void mergeAppendedData(std::vector<Ttype> &data, const std::vector<size_t> &distribution);

	template<typename I>
	constexpr I my_int_pow(I base, I exponent)
	{
		if(exponent < 0) return 0;
		if(exponent == 0) return 1;
		if(exponent == 1) return base;
		I tmp = my_int_pow<I>(base, exponent / 2);
		I result = tmp * tmp;
		if(exponent % 2 == 1) result *= base;
		return result;
	}

	template<typename T>
	constexpr size_t get_max_val_no_precision_loss_in_fp()
	{
		constexpr size_t result = my_int_pow<long long>(std::numeric_limits<T>::radix, std::numeric_limits<T>::digits);
		return result;
	}
};

}

#include "utils.hpp"


#endif /* BASIS_UTILITIES_UTILS_H_ */
