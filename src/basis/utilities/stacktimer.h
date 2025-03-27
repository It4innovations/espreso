
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
        if(!enabled) {
            return;
        }
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        eslog::info("rank %3d thread %3d %*sstarted '%s'\n", info::mpi::rank, thread, stk.size() * indent, "", name);
        double start_time = omp_get_wtime();
        stk.push(stackitem{name, start_time});
    }
    void instance_pop()
    {
        if(!enabled) {
            return;
        }
        double stop_time = omp_get_wtime();
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        stackitem & item = stk.top();
        double time_ms = (stop_time - item.start_time) * 1000;
        int width = 120;
        width -= (stk.size() - 1) * indent;
        width -= 8; // "rank 123"
        width -= 11; // " thread 123"
        width -= 9; // " finishd "
        width -= strlen(item.name) + 2; // "'%s' "
        width -= 16; // " %12.3f ms"
        width = std::max(width, 0);
        eslog::info("rank %3d thread %3d %*sfinishd '%s' %*s %12.3f ms\n", info::mpi::rank, thread, (stk.size() - 1) * indent, "", item.name, width, "", time_ms);
        stk.pop();
    }
    template<typename... Args>
    void instance_info(const char * fmt, Args... args)
    {
        if(!enabled) {
            return;
        }
        char buffer[1024];
        int thread = omp_get_thread_num();
        std::stack<stackitem> & stk = data[thread].stk;
        snprintf(buffer, sizeof(buffer), "rank %3d thread %3d %*sinfo %s\n", info::mpi::rank, thread, (int)(stk.size() * indent), "", fmt);
        eslog::info(buffer, args...);
    }
    void instance_finish()
    {
        data = std::vector<stackwrapper>();
    }
    void instance_enable()
    {
        enabled = true;
    }
    void instance_disable()
    {
        enabled = false;
    }
private:
    std::vector<stackwrapper> data;
    int indent = 2;
    bool enabled = false;
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
