
#ifndef BASIS_UTILITIES_MINMAXAVG_H_
#define BASIS_UTILITIES_MINMAXAVG_H_

#include <functional>
#include <type_traits>
#include <algorithm>
#include <limits>

#include "wrappers/mpi/communication.h"



namespace espreso {

    template<typename T>
	struct minmaxavg
	{
		T min;
		T max;
		T avg;

		static minmaxavg<T> mpi_reduce_and_return(T min, T max, T sum, size_t n)
		{
			Communication::allReduce(&min, nullptr, 1, MPITools::getType<T>().mpitype, MPI_MIN);
			Communication::allReduce(&max, nullptr, 1, MPITools::getType<T>().mpitype, MPI_MAX);
			Communication::allReduce(&sum, nullptr, 1, MPITools::getType<T>().mpitype, MPI_SUM);
			Communication::allReduce(&n, nullptr, 1, MPITools::getType<size_t>().mpitype, MPI_SUM);

			minmaxavg<T> mma;
			mma.min = min;
			mma.max = max;
			mma.avg = sum / n;
			return mma;
		}

		static minmaxavg<T> compute_from_allranks(T val)
		{
			return mpi_reduce_and_return(val, val, val, 1);
		}

		template<typename ITER>
		static minmaxavg<T> compute_from_my_rank_only(ITER begin, ITER end, std::function<T(const typename ITER::value_type &)> f = [](typename ITER::value_type val){ return val; })
		{
			T min = std::numeric_limits<T>::max();
			T max = std::numeric_limits<T>::min();
			T sum = 0;
			size_t n = 0;

			for(ITER it = begin; it != end; ++it)
			{
				T val = f(*it);
				min = std::min(min, val);
				max = std::max(max, val);
				sum += val;
				n += 1;
			}

			minmaxavg<T> mma;
			mma.min = min;
			mma.max = max;
			mma.avg = sum / n;
			return mma;
		}

		template<typename ITER>
		static minmaxavg<T> compute_from_allranks(ITER begin, ITER end, std::function<T(const typename ITER::value_type &)> f = [](typename ITER::value_type val){ return val; })
		{
			T min = std::numeric_limits<T>::max();
			T max = std::numeric_limits<T>::min();
			T sum = 0;
			size_t n = 0;

			for(ITER it = begin; it != end; ++it)
			{
				T val = f(*it);
				min = std::min(min, val);
				max = std::max(max, val);
				sum += val;
				n += 1;
			}

			return mpi_reduce_and_return(min, max, sum, n);
		}

		std::string to_string(const char * description)
		{
			char str[200];
			if constexpr(std::is_same_v<T,float>)                   snprintf(str, sizeof(str), " = %-55s    %8.2f <%8.2f - %8.2f> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,double>)             snprintf(str, sizeof(str), " = %-55s    %8.2f <%8.2f - %8.2f> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,int>)                snprintf(str, sizeof(str), " = %-55s    %8d <%8d - %8d> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned int>)       snprintf(str, sizeof(str), " = %-55s    %8u <%8u - %8u> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,long>)               snprintf(str, sizeof(str), " = %-55s    %8ld <%8ld - %8ld> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned long>)      snprintf(str, sizeof(str), " = %-55s    %8lu <%8lu - %8lu> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,long long>)          snprintf(str, sizeof(str), " = %-55s    %8lld <%8lld - %8lld> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned long long>) snprintf(str, sizeof(str), " = %-55s    %8llu <%8llu - %8llu> = \n", description, avg, min, max);
			else if constexpr(std::is_same_v<T,size_t>)             snprintf(str, sizeof(str), " = %-55s    %8zu <%8zu - %8zu> = \n", description, avg, min, max);
			else static_assert(sizeof(T*) == 0 /* aka false */, "minmaxavg is unable to print the type"); // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert, https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p2593r0.html
			return std::string(str);
		}

		std::string to_string_plain()
		{
			char str[200];
			if constexpr(std::is_same_v<T,float>)                   snprintf(str, sizeof(str), "%8.2f <%8.2f - %8.2f>", avg, min, max);
			else if constexpr(std::is_same_v<T,double>)             snprintf(str, sizeof(str), "%8.2f <%8.2f - %8.2f>", avg, min, max);
			else if constexpr(std::is_same_v<T,int>)                snprintf(str, sizeof(str), "%8d <%8d - %8d>", avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned int>)       snprintf(str, sizeof(str), "%8u <%8u - %8u>", avg, min, max);
			else if constexpr(std::is_same_v<T,long>)               snprintf(str, sizeof(str), "%8ld <%8ld - %8ld>", avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned long>)      snprintf(str, sizeof(str), "%8lu <%8lu - %8lu>", avg, min, max);
			else if constexpr(std::is_same_v<T,long long>)          snprintf(str, sizeof(str), "%8lld <%8lld - %8lld>", avg, min, max);
			else if constexpr(std::is_same_v<T,unsigned long long>) snprintf(str, sizeof(str), "%8llu <%8llu - %8llu>", avg, min, max);
			else if constexpr(std::is_same_v<T,size_t>)             snprintf(str, sizeof(str), "%8zu <%8zu - %8zu>", avg, min, max);
			else static_assert(sizeof(T*) == 0 /* aka false */, "minmaxavg is unable to print the type"); // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert, https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2022/p2593r0.html
			return std::string(str);
		}
	};
}


#endif /* BASIS_UTILITIES_MINMAXAVG_H_ */
