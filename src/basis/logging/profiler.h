
#ifndef SRC_BASIS_LOGGING_PROFILER_H_
#define SRC_BASIS_LOGGING_PROFILER_H_

namespace espreso {
namespace profiler {

// in order to not be dependent on info::system
enum class BUILD_MODE {
    RELEASE,
    DEVEL,
    DEBUG,
    PROFILE
};

static constexpr BUILD_MODE buildmode()
{
#ifdef __ESMODE__
    return BUILD_MODE::__ESMODE__;
#else
    return BUILD_MODE::RELEASE;
#endif
}

static constexpr bool active()
{
    return buildmode() == BUILD_MODE::PROFILE;
}

template<bool active>
struct Event {

};

template<>
struct Event<false> {
    template <typename TType> void min(const TType &value) {}
    template <typename TType> void max(const TType &value) {}
};

template<>
struct Event<true> {
    const char* name;
    int level;

    enum TYPE: int {
        EMPTY      = 1 << 0,
        SYNC       = 1 << 1,
        START      = 1 << 2,
        CHECKPOINT = 1 << 3,
        END        = 1 << 4,
        INT        = 1 << 5,
        SIZE       = 1 << 6,
        DOUBLE     = 1 << 7,
        PARAM      = 1 << 8
    } type;

    union Data {
        double time;
        long   ivalue;
        long unsigned int svalue;
        double dvalue;
    } data;

    template <typename TType> void min(const TType &value);
    template <typename TType> void max(const TType &value);
};

inline Event<true>::TYPE operator|(Event<true>::TYPE t1, Event<true>::TYPE t2) { return static_cast<Event<true>::TYPE>(static_cast<int>(t1) | static_cast<int>(t2)); }

template <typename TType> inline Event<true>::TYPE type(const TType &type);
template <> inline Event<true>::TYPE type(const int &type) { return Event<true>::TYPE::INT; }
template <> inline Event<true>::TYPE type(const long &type) { return Event<true>::TYPE::INT; }
template <> inline Event<true>::TYPE type(const long unsigned int &type) { return Event<true>::TYPE::SIZE; }
template <> inline Event<true>::TYPE type(const float &type) { return Event<true>::TYPE::DOUBLE; }
template <> inline Event<true>::TYPE type(const double &type) { return Event<true>::TYPE::DOUBLE; }

template <typename TType> inline Event<true>::Data data(const TType &data);
template <> inline Event<true>::Data data(const int &data) { return Event<true>::Data { .ivalue = data }; }
template <> inline Event<true>::Data data(const long &data) { return Event<true>::Data { .ivalue = data }; }
template <> inline Event<true>::Data data(const long unsigned int &data) { return Event<true>::Data { .svalue = data }; }
template <> inline Event<true>::Data data(const float &data) { return Event<true>::Data { .dvalue = data }; }
template <> inline Event<true>::Data data(const double &data) { return Event<true>::Data { .dvalue = data }; }

template <typename TType> inline TType& datavalue(Event<true>::Data &data);
template <> inline long& datavalue(Event<true>::Data &data) { return data.ivalue; };
template <> inline long unsigned int& datavalue(Event<true>::Data &data) { return data.svalue; };
template <> inline double& datavalue(Event<true>::Data &data) { return data.dvalue; };

template <> inline void Event<true>::min(const int &value) { if (value < datavalue<long>(data)) datavalue<long>(data) = value; }
template <> inline void Event<true>::max(const int &value) { if (datavalue<long>(data) < value) datavalue<long>(data) = value; }
template <typename TType> inline void Event<true>::min(const TType &value) { if (value < datavalue<TType>(data)) datavalue<TType>(data) = value; }
template <typename TType> inline void Event<true>::max(const TType &value) { if (datavalue<TType>(data) < value) datavalue<TType>(data) = value; }

template<bool active>
struct Profiler {
};

template<>
class Profiler<false> {
public:
    static Profiler instance;
    static double time();

    Event<false>& syncstart(const char* name) { return dummy; }
    Event<false>& synccheckpoint(const char* name) { return dummy; }
    Event<false>& syncend(const char* name) { return dummy; }
    template<typename TValue>
    Event<false>& syncparam(const char* name, const TValue &value) { return dummy; }

    Event<false>& istart(const char* name) { return dummy; }
    Event<false>& icheckpoint(const char* name) { return dummy; }
    Event<false>& iend(const char* name) { return dummy; }
    template<typename TValue>
    Event<false>& iparam(const char* name, const TValue &value) { return dummy; }

    void print() {}
    static double initted;
    Event<false> dummy;
};

template<>
class Profiler<true> {
public:
    static Profiler instance;
    static double time();

    Event<true>& add(const char* name, int level, const Event<true>::TYPE &type, const Event<true>::Data &data)
    {
        events[size] = Event<true>{ name, level, type, data };
        return events[size++];
    }

    Event<true>& syncstart(const char* name) { return add(name, ++level, Event<true>::TYPE::START | Event<true>::TYPE::SYNC, Event<true>::Data{ .time = time() }); }
    Event<true>& synccheckpoint(const char* name) { return add(name, level, Event<true>::TYPE::CHECKPOINT | Event<true>::TYPE::SYNC, Event<true>::Data{ .time = time() }); }
    Event<true>& syncend(const char* name) { return add(name, level--, Event<true>::TYPE::END | Event<true>::TYPE::SYNC, Event<true>::Data{ .time = time() }); }
    template<typename TValue>
    Event<true>& syncparam(const char* name, const TValue &value) { return add(name, level, type(value) | Event<true>::TYPE::SYNC | Event<true>::TYPE::PARAM, data(value)); }

    Event<true>& istart(const char* name) { return add(name, ++level, Event<true>::TYPE::START, Event<true>::Data{ .time = time() }); }
    Event<true>& icheckpoint(const char* name) { return add(name, level, Event<true>::TYPE::CHECKPOINT, Event<true>::Data{ .time = time() }); }
    Event<true>& iend(const char* name) { return add(name, level--, Event<true>::TYPE::END, Event<true>::Data{ .time = time() }); }
    template<typename TValue>
    Event<true>& iparam(const char* name, const TValue &value) { return add(name, level, type(value) | Event<true>::TYPE::PARAM, data(value)); }


    Profiler();
    ~Profiler();
    void print();

    static double initted;
    int level, size;
    Event<true> *events;
};

inline Event<active()>& syncstart(const char* name)
{
    return Profiler<active()>::instance.syncstart(name);
}

inline Event<active()>& synccheckpoint(const char* name)
{
    return Profiler<active()>::instance.synccheckpoint(name);
}

inline Event<active()>& syncend(const char* name)
{
    return Profiler<active()>::instance.syncend(name);
}

template <typename TValue>
inline Event<active()>& syncparam(const char* name, const TValue &value)
{
    return Profiler<active()>::instance.syncparam(name, value);
}


inline Event<active()>& start(const char* name)
{
    return Profiler<active()>::instance.istart(name);
}

inline Event<active()>& checkpoint(const char* name)
{
    return Profiler<active()>::instance.icheckpoint(name);
}

inline Event<active()>& end(const char* name)
{
    return Profiler<active()>::instance.iend(name);
}

template <typename TValue>
inline Event<active()>& param(const char* name, const TValue &value)
{
    return Profiler<active()>::instance.iparam(name, value);
}

inline void print()
{
    Profiler<active()>::instance.print();
}

}
}

#endif /* SRC_BASIS_LOGGING_PROFILER_H_ */
