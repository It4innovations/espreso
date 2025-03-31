
#ifndef SRC_BASIS_UTILITIES_STACKTIMER_H
#define SRC_BASIS_UTILITIES_STACKTIMER_H

#include <vector>
#include <stack>
#include <cstring>
#include <omp.h>

#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"



namespace espreso {



#ifdef ESPRESO_STACKTIMER_ENABLE

static constexpr size_t get_padding_size(size_t struct_size, size_t align) {
    size_t size = ((struct_size - 1) / align + 1) * align;
    return size - struct_size;
}

class stacktimer
{
private:
    struct stackitem
    {
        const char * name;
        double start_time;
    };
    struct stackwrapper
    {
        std::stack<stackitem> stk;
        char padding[get_padding_size(sizeof(std::stack<stackitem>),64)];
    };
private:
    void instance_init()
    {
        data.resize(omp_get_max_threads());
    }
    void instance_push(const char * name)
    {
        if(enabled <= 0) {
            return;
        }
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        int depth = stk.size();
        if(thread != 0) depth += main_thread_depth;
        int num_spaces = depth * indent;
        eslog::info("rank %3d thread %3d %*sstarted '%s'\n", info::mpi::rank, thread, num_spaces, "", name);
        double start_time = omp_get_wtime();
        stk.push(stackitem{name, start_time});
        if(thread == 0 && omp_get_num_threads() == 1) main_thread_depth++;
    }
    void instance_pop()
    {
        if(enabled <= 0) {
            return;
        }
        double stop_time = omp_get_wtime();
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        stackitem & item = stk.top();
        double time_ms = (stop_time - item.start_time) * 1000;
        int depth = stk.size() - 1;
        if(thread != 0) depth += main_thread_depth;
        int num_spaces = depth * indent;
        int width = 120;
        width -= num_spaces;
        width -= 8; // "rank 123"
        width -= 11; // " thread 123"
        width -= 9; // " finishd "
        width -= strlen(item.name) + 2; // "'%s' "
        width -= 16; // " %12.3f ms"
        width = std::max(width, 0);
        eslog::info("rank %3d thread %3d %*sfinishd '%s' %*s %12.3f ms\n", info::mpi::rank, thread, num_spaces, "", item.name, width, "", time_ms);
        stk.pop();
        if(thread == 0 && omp_get_num_threads() == 1) main_thread_depth--;
    }
    template<typename... Args>
    void instance_info(const char * fmt, Args... args)
    {
        if(enabled <= 0) {
            return;
        }
        char buffer[1024];
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        int depth = stk.size();
        if(thread != 0) depth += main_thread_depth;
        int num_spaces = depth * indent;
        snprintf(buffer, sizeof(buffer), "rank %3d thread %3d %*sinfo %s\n", info::mpi::rank, thread, num_spaces, "", fmt);
        eslog::info(buffer, args...);
    }
    void instance_finish()
    {
        data = std::vector<stackwrapper>();
    }
    void instance_enable()
    {
        #pragma omp atomic
        enabled++;
    }
    void instance_disable()
    {
        #pragma omp atomic
        enabled--;
    }
private:
    std::vector<stackwrapper> data;
    int main_thread_depth = 0;
    int indent = 2;
    int enabled = 0;
private:
    static stacktimer instance;
public:
    static void init()
    {
        instance.instance_init();
    }
    static void push(const char * name)
    {
        instance.instance_push(name);
    }
    static void pop()
    {
        instance.instance_pop();
    }
    template<typename... Args>
    static void info(const char * fmt, Args... args)
    {
        instance.instance_info(fmt, args...);
    }
    static void finish()
    {
        instance.instance_finish();
    }
    static void enable()
    {
        instance.instance_enable();
    }
    static void disable()
    {
        instance.instance_disable();
    }
};

#else

class stacktimer
{
private:
    static stacktimer instance;
public:
    static void init() {}
    static void push(const char * name) {}
    static void pop() {}
    template<typename... Args>
    static void info(const char * fmt, Args... args) {}
    static void finish() {}
    static void enable() {}
    static void disable() {}
};

#endif



}



#endif /* SRC_BASIS_UTILITIES_STACKTIMER_H */
