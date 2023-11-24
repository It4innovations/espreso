
#pragma once

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <type_traits>

#include "my_timer.hpp"






namespace espreso {




#ifndef MY_ABORT
#define MY_ABORT(message) do { my_abort(message, __FILE__, __LINE__); } while(false)
static void my_abort(const char * message, const char * file, int line)
{
    fprintf(stderr, "MY_ABORT: \"%s\". In file %s on line %d\n", message, file, line);
    fflush(stdout);
    fflush(stderr);
    exit(1);
}
#endif


template<typename T>
static std::vector<T> jv(const std::vector<T> & v1, const std::vector<T> & v2)
{
    std::vector<T> ret;
    ret.reserve(v1.size() + v2.size());
    ret.insert(ret.end(), v1.begin(), v1.end());
    ret.insert(ret.end(), v2.begin(), v2.end());
    return ret;
}

template<typename T>
static void av(std::vector<T> & v, const std::vector<T> & vin)
{
    v.insert(v.end(), vin.begin(), vin.end());
}





static void print_timer(const char * name, const my_timer & tm, int indent = 0)
{
    int number_width = 20;
    int number_precision = 6;
    int indent_size = 2;
    int label_width = 40;
    bool allow_irellevant_values = false;

    int spaces_before = indent * indent_size;
    int max_name_length = label_width - spaces_before - 1;
    int namelength = std::strlen(name);
    int spbef_name_length = spaces_before + std::min(namelength, max_name_length);
    int spaces_after = label_width - spbef_name_length - 1;

    printf("%*.*s:%*s", spbef_name_length, max_name_length, name, spaces_after, "");

    bool was_parallel = (tm.get_n_used_threads() > 1);
    bool print_total = allow_irellevant_values || !was_parallel;
    bool print_para = allow_irellevant_values || was_parallel;

    if(print_total) printf("%*.*f", number_width, number_precision, 1000.0 * tm.total_time());
    else printf("%*s", number_width, "");

    printf("%*.*f", number_width, number_precision, 1000.0 * tm.avg_time());

    if(print_para) printf("%*.*f", number_width, number_precision, 1000.0 * tm.avg_perthread_time());
    else printf("%*s", number_width, "");

    // printf("    ");
    // for(size_t i = 0; i < tm.get_n_used_threads(); i++)
    //     printf(" %*.*f", number_width, number_precision, 1000.0 * tm.perthread_time(i));

    printf("\n");
}



template<typename T>
static bool my_isnan(T val)
{
    static_assert(std::is_same_v<T,float> || std::is_same_v<T,double>);
    char str[100];
    sprintf(str, "%f", val);
    return (strstr(str, "nan") != nullptr);
}



template<typename F1, typename F2>
void loop_factorize_assemble(size_t n_iters, char waitcontpar, char joinsplit, my_timer & tm_factorize_outer, my_timer & tm_assemble_outer, F1 factorize, F2 assemble)
{
    if(joinsplit == 'J')
    {
        #pragma omp parallel for schedule(static,1) if(waitcontpar == 'P')
        for(size_t i = 0; i < n_iters; i++)
        {
            tm_factorize_outer.start();
            factorize(i);
            tm_factorize_outer.stop();
            tm_assemble_outer.start();
            assemble(i);
            tm_assemble_outer.stop();
        }
    }
    else
    {
        tm_factorize_outer.start();
        #pragma omp parallel for schedule(static,1) if(joinsplit == 'P')
        for(size_t i = 0; i < n_iters; i++) factorize(i);
        tm_factorize_outer.stop();
        tm_assemble_outer.start();
        #pragma omp parallel for schedule(static,1) if(waitcontpar == 'P')
        for(size_t i = 0; i < n_iters; i++) assemble(i);
        tm_assemble_outer.stop();
    }
}





template<typename T>
class my_stdallocator_wrapper
{
public:
    static constexpr bool is_data_host_accessible = true;
public:
    using value_type = T;
    std::allocator<T> a;
    my_stdallocator_wrapper() = default;
    T * allocate(size_t count)
    {
        return a.allocate(count);
    }
    void deallocate(T * ptr, size_t count)
    {
        a.deallocate(ptr, count);
    }
};





template<template<typename> typename A>
class my_limited_allocator_outer
{
public:
    template<typename T>
    class inner
    {
        // the capacity and size variables have to have larger or equal scope as all objects using this allocator
        // also this allocator wrapper does not account for overheads due to alignment or similar. Use just as an approximate limiter
    public:
        static constexpr bool is_data_host_accessible = A<T>::is_data_host_accessible;
    public:
        const size_t & capacity;
        size_t & current_size;
        A<T> inner_allocator;
        using value_type = T;
        inner(const size_t & capacity_, size_t & current_size_, A<T> inner_allocator_ = A<T>()) : capacity(capacity_), current_size(current_size_), inner_allocator(inner_allocator_) {}
        T * allocate(size_t count)
        {
            size_t new_size;

            #pragma omp atomic capture
            { current_size += count * sizeof(T); new_size = current_size; }
            if(new_size > capacity) MY_ABORT("Limited allocator exceeded its capacity");

            T * ptr = inner_allocator.allocate(count);
            return ptr;
        }
        void deallocate(T * ptr, size_t count)
        {
            if(ptr != nullptr)
            {
                #pragma omp atomic
                current_size -= count * sizeof(T);
            }

            inner_allocator.deallocate(ptr, count);
        }
    };
};

template<typename T, template<typename> typename A = my_stdallocator_wrapper>
using my_limited_allocator = typename my_limited_allocator_outer<A>::template inner<T>;





    // char strategy     = ms_get_strategy    (magicstring);
    // char method       = ms_get_method      (magicstring);
    // char device       = ms_get_device      (magicstring);
    // char joinsplit    = ms_get_joinsplit   (magicstring);
    // char waitcontpar  = ms_get_waitcontpar (magicstring);
    // char ordering     = ms_get_ordering    (magicstring);
    // char spdnfactor   = ms_get_spdnfactor  (magicstring);
    // char todensewhere = ms_get_todensewhere(magicstring);
    // char transcallman = ms_get_transcallman(magicstring);
    // char trsvtrsm     = ms_get_trsvtrsm    (magicstring);
    // char gemmsyrk     = ms_get_gemmsyrk    (magicstring);
    // char scaling      = ms_get_scaling     (magicstring);

static inline void check_magicstring_length(const char * magicstring) { if(strlen(magicstring) != 12) MY_ABORT("WRONG magicstring length"); }
static inline bool check_arg(char arg, const char * options) { return strchr(options, arg) != nullptr; }

static inline char ms_get_strategy    (const char * magicstring) { char arg = magicstring[ 0]; if(!check_arg(arg, "IE"))  MY_ABORT("magicstring strategy is WRONG");     return arg; }
static inline char ms_get_method      (const char * magicstring) { char arg = magicstring[ 1]; if(!check_arg(arg, "DPM")) MY_ABORT("magicstring method is WRONG");       return arg; }
static inline char ms_get_device      (const char * magicstring) { char arg = magicstring[ 2]; if(!check_arg(arg, "CG"))  MY_ABORT("magicstring device is WRONG");       return arg; }
static inline char ms_get_joinsplit   (const char * magicstring) { char arg = magicstring[ 3]; if(!check_arg(arg, "JSP")) MY_ABORT("magicstring joinsplit is WRONG");    return arg; }
static inline char ms_get_waitcontpar (const char * magicstring) { char arg = magicstring[ 4]; if(!check_arg(arg, "WCP")) MY_ABORT("magicstring waitcontpar is WRONG");  return arg; }
static inline char ms_get_ordering    (const char * magicstring) { char arg = magicstring[ 5]; if(!check_arg(arg, "AM"))  MY_ABORT("magicstring ordering is WRONG");     return arg; }
static inline char ms_get_spdnfactor  (const char * magicstring) { char arg = magicstring[ 6]; if(!check_arg(arg, "SDO")) MY_ABORT("magicstring spdnfactor is WRONG");   return arg; }
static inline char ms_get_todensewhere(const char * magicstring) { char arg = magicstring[ 7]; if(!check_arg(arg, "HD"))  MY_ABORT("magicstring todensewhere is WRONG"); return arg; }
static inline char ms_get_transcallman(const char * magicstring) { char arg = magicstring[ 8]; if(!check_arg(arg, "CM"))  MY_ABORT("magicstring transcallman is WRONG"); return arg; }
static inline char ms_get_trsvtrsm    (const char * magicstring) { char arg = magicstring[ 9]; if(!check_arg(arg, "VM"))  MY_ABORT("magicstring trsvtrsm is WRONG");     return arg; }
static inline char ms_get_gemmsyrk    (const char * magicstring) { char arg = magicstring[10]; if(!check_arg(arg, "GS"))  MY_ABORT("magicstring gemmsyrk is WRONG");     return arg; }
static inline char ms_get_scaling     (const char * magicstring) { char arg = magicstring[11]; if(!check_arg(arg, "CD"))  MY_ABORT("magicstring scaling is WRONG");      return arg; }

static void print_args_help()
{
    printf("mesh dimensionality                                              {2,3}\n");
    printf("element type         shape of element and number of vertices     {TRIANGLE3,SQUARE4;TETRA4,HEXA8}\n");
    printf("mesh size            number of elements along the edge\n");
    printf("n matrices           number of domains\n");
    printf("number of applys     number of F applications\n");
    printf("version string:\n");
    printf("  strategy           strategy                                    {Implicit, Explicit}\n");
    printf("  method             explicit method                             {Direct, sc Partial factorization, sc My}\n");
    printf("  device             sycl device                                 {Cpu, Gpu}\n");
    printf("  joinsplit          explicit direct factorize loop parallelism  {Join with main loop, split factorize Sequential, split factorize Parallel}\n");
    printf("  waitcontpar        main loop parallelism                       {sequential Wait, sequential Continue, Parallel}\n");
    printf("  ordering           cholmod permutation ordering algorithm      {Amd, Metis}\n");
    printf("  spdnfactor         factor sparsity in trs                      {Sparse, Dense, sparse with Optimize}\n");
    printf("  todensewhere       where to convert sparse to dense            {on Host, on Device}\n");
    printf("  transcallman       transposed matrices handling                {Call functions with trans flag, Manually transpose matrices}\n");
    printf("  trsvtrsm           right hand side of trs                      {trsV, trsM}\n");
    printf("  gemmsyrk           gemmsyrk                                    {trs+trs+Gemm, trs+Syrk}\n");
    printf("  scaling            mysc division                               {Constant, Decreasing}\n");
}





void run_dummy_parallel_region()
{
    #pragma omp parallel
    {

    }
}

}
