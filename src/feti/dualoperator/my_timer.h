
#ifndef SRC_FETI_DUALOPERATOR_MY_TIMER_H_
#define SRC_FETI_DUALOPERATOR_MY_TIMER_H_

#define ENABLE_DUALOP_EXPLICIT_GPU_TIMERS

#include <chrono>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "esinfo/mpiinfo.h"

namespace espreso {

#ifdef ENABLE_DUALOP_EXPLICIT_GPU_TIMERS

class my_timer
{
private:
    using clocktype = std::chrono::steady_clock;
    static constexpr double ticks_to_sec = ((double)clocktype::period::num) / clocktype::period::den; // duration of one tick
private:
    std::vector<clocktype::time_point> start_times;
    std::vector<long long> ticks_per_thread;
    std::vector<long long> ticks_per_thread_last_laps;
    long long total_ticks;
    long long n_stops;
    int nested_level_created;
    bool enabled;
public:
    my_timer(bool enabled_ = true) : enabled{enabled_}
    {
        if(!enabled) return;
        this->reset();
        start_times.resize(get_thread_count());
        ticks_per_thread_last_laps.resize(get_thread_count());
        ticks_per_thread.resize(get_thread_count(), 0LL);
        nested_level_created = omp_get_level();
    }
    void start()
    {
        if(!enabled) return;
        int idx;
        if(omp_get_level() == nested_level_created) idx = 0;
        else if(omp_get_level() == nested_level_created + 1) idx = get_thread_idx();
        else throw std::runtime_error("my timer can only be used in the same openmp nested level or one level deeper");
        this->start_times[idx] = clocktype::now();
    }
    void stop()
    {
        if(!enabled) return;
        int idx;
        if(omp_get_level() == nested_level_created) idx = 0;
        else if(omp_get_level() == nested_level_created + 1) idx = get_thread_idx();
        else throw std::runtime_error("my timer can only be used in the same openmp nested level or one level deeper");
        auto stop_time = clocktype::now();
        auto start_time = start_times[idx];
        long long n_ticks = (stop_time - start_time).count();
        this->ticks_per_thread[idx] += n_ticks;
        this->ticks_per_thread_last_laps[idx] = n_ticks;
#pragma omp atomic
        this->total_ticks += n_ticks;
#pragma omp atomic
        this->n_stops += 1;
    }
    void reset()
    {
        if(!enabled) return;
        this->total_ticks = 0LL;
        this->n_stops = 0LL;
    }
    double total_time() const
    {
        if(!enabled) return 0.0;
        return this->total_ticks * ticks_to_sec;
    }
    double avg_time() const
    {
        if(!enabled) return 0.0;
        return this->total_time() / n_stops;
    }
    double avg_perthread_time() const
    {
        if(!enabled) return 0.0;
        return this->total_time() / this->get_n_used_threads();
    }
    double perthread_time(int thread_idx) const
    {
        if(!enabled) return 0.0;
        return this->ticks_per_thread[thread_idx] * ticks_to_sec;
    }
    double lap_time(int thread_idx = -1)
    {
        if(!enabled) return 0.0;
        if(thread_idx < 0) thread_idx = get_thread_idx();
        return ticks_per_thread_last_laps[thread_idx] * ticks_to_sec;
    }
    int get_n_used_threads() const
    {
        if(!enabled) return 1;
        return std::count_if(start_times.begin(), start_times.end(), [](clocktype::time_point tp) { return tp.time_since_epoch().count() > 0LL; });
    }
    bool is_enabled() const
    {
        return enabled;
    }
private:
    int get_thread_idx() const
    {
        int idx = 0;
#ifdef _OPENMP
        idx = omp_get_thread_num();
#endif
        return idx;
    }
    int get_thread_count() const
    {
        int count = 1;
#ifdef _OPENMP
        count = omp_get_max_threads();
#endif
        return count;
    }
};

#else

class my_timer
{
public:
    my_timer() { }
    void start() { }
    void stop() { }
    void reset() { }
    double total_time() const { return 0.0; }
    double avg_time() const { return 0.0; }
    double avg_perthread_time() const { return 0.0; }
    double perthread_time(int) const { return 0.0; }
    double lap_time(int) { return 0.0; }
    int get_n_used_threads() const { return 1; }
};

#endif



inline void print_timer(const char * name, const my_timer & tm, int indent = 0)
{
#ifdef ENABLE_DUALOP_EXPLICIT_GPU_TIMERS
    if(!tm.is_enabled()) return;

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

    eslog::info("rank %4d %*.*s:%*s", info::mpi::rank, spbef_name_length, max_name_length, name, spaces_after, "");

    bool was_parallel = (tm.get_n_used_threads() > 1);
    bool print_total = allow_irellevant_values || !was_parallel;
    bool print_para = allow_irellevant_values || was_parallel;

    if(print_total) eslog::info("%*.*f", number_width, number_precision, 1000.0 * tm.total_time());
    else eslog::info("%*s", number_width, "");

    eslog::info("%*.*f", number_width, number_precision, 1000.0 * tm.avg_time());

    if(print_para) eslog::info("%*.*f", number_width, number_precision, 1000.0 * tm.avg_perthread_time());
    else eslog::info("%*s", number_width, "");

    // eslog::info("    ");
    // for(size_t i = 0; i < tm.get_n_used_threads(); i++)
    //     eslog::info(" %*.*f", number_width, number_precision, 1000.0 * tm.perthread_time(i));

    eslog::info("\n");
#endif
}

}

#endif /* SRC_FETI_DUALOPERATOR_MY_TIMER_H_ */
