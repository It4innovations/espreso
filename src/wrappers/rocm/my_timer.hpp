
#pragma once

#include <chrono>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace espreso {

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
public:
    my_timer()
    {
        this->reset();
        start_times.resize(get_thread_count());
        ticks_per_thread_last_laps.resize(get_thread_count());
        ticks_per_thread.resize(get_thread_count(), 0LL);
    }
    void start()
    {
        this->start_times[get_thread_idx()] = clocktype::now();
    }
    void stop()
    {
        auto stop_time = clocktype::now();
        auto start_time = start_times[get_thread_idx()];
        long long n_ticks = (stop_time - start_time).count();
        this->ticks_per_thread[get_thread_idx()] += n_ticks;
        this->ticks_per_thread_last_laps[get_thread_idx()] = n_ticks;
#pragma omp atomic
        this->total_ticks += n_ticks;
#pragma omp atomic
        this->n_stops += 1;
    }
    void reset()
    {
        this->total_ticks = 0LL;
        this->n_stops = 0LL;
    }
    double total_time() const
    {
        return this->total_ticks * ticks_to_sec;
    }
    double avg_time() const
    {
        return this->total_time() / n_stops;
    }
    double avg_perthread_time() const
    {
        return this->total_time() / this->get_n_used_threads();
    }
    double perthread_time(int thread_idx) const
    {
        return this->ticks_per_thread[thread_idx] * ticks_to_sec;
    }
    double lap_time(int thread_idx = -1)
    {
        if(thread_idx < 0) thread_idx = get_thread_idx();
        return ticks_per_thread_last_laps[thread_idx] * ticks_to_sec;
    }
    int get_n_used_threads() const
    {
        return std::count_if(start_times.begin(), start_times.end(), [](clocktype::time_point tp) { return tp.time_since_epoch().count() > 0LL; });
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

}
