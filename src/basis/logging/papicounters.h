
#ifndef SRC_BASIS_LOGGING_PAPICOUNTERS_H_
#define SRC_BASIS_LOGGING_PAPICOUNTERS_H_

#include "verbosity.h"
#include "basis/containers/allocators.h"
#include "esinfo/ecfinfo.h"
#include "wrappers/papi/w.papi.h"

#include <cstddef>
#include <vector>

namespace espreso {

class PAPICounters: public Verbosity<PAPICounters, 'p'> {
public:
    struct Event {
        const char* name;
        enum Type {
            START, CHECKPOINT, ACCUMULATED, END, DURATION, LOADSTEP
        } type;
    };

    void initOutput()
    {
        papi.init();
        if (papi.values) {
            _events.reserve(1000000);
            _values.resize(1000000 * papi.values);
            papi.read(_values.data());
        }
        for (size_t i = 1; i < _events.size(); ++i) {
            for (int v = 0; v < papi.values; ++v) {
                _values[papi.values * i + v] = _values[v];
            }
        }
    }

    void add(const char* region, Event::Type type)
    {
        if (_events.size() == _events.capacity()) {
            _events.reserve(2 * _events.capacity());
            _values.resize(2 * _values.capacity());
        }
        papi.read(_values.data() + papi.values * _events.size());
        _events.push_back(Event{ region, type });
    }

    void start(const char* region, const char* section)
    {
        add(section, Event::START);
    }

    void checkpoint(const char* region)
    {
        add(region, Event::CHECKPOINT);
    }

    void accumulated(const char* region)
    {
        add(region, Event::ACCUMULATED);
    }

    void end(const char* region)
    {
        add(region, Event::END);
        add(region, Event::DURATION);
    }

    void param(const char* name, const int &value)
    {
        // do nothing
    }

    void param(const char* name, const long &value)
    {
        // do nothing
    }

    void param(const char* name, const long unsigned int &value)
    {
        // do nothing
    }

    void param(const char* name, const double &value)
    {
        // do nothing
    }

    void param(const char* name, const char* value)
    {
        // do nothing
    }

    void ln()
    {
        // do nothing
    }

    void nextLoadStep(int step)
    {
        add(nullptr, Event::LOADSTEP);
    }

    void output(const char* msg, VerboseArg::COLOR color)
    {
        // do nothing
    }

    void error(const char* msg)
    {

    }

    void finish();

protected:
    std::vector<Event> _events;
    std::vector<long long, initless_allocator<long long> > _values;
    PAPI papi;
};

}




#endif /* SRC_BASIS_LOGGING_PAPICOUNTERS_H_ */
