
#include "papicounters.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"
#include "wrappers/mpi/communication.h"

#include "omp.h"
#include <cstdio>
#include <cstring>
#include <ctime>

using namespace espreso;

void PAPICounters::finish()
{
    if (verbosity == 0 || papi.values == 0) {
        return;
    }
    if (papi.status) {
        eslog::info("\n ============================================================================================= \n");
        if (info::ecf->output.papi_events.size()) {
            eslog::info(" == INVALID PAPI EVENTS %67s == \n", info::ecf->output.papi_events.c_str());
        }
        if (info::ecf->output.papi_codes.size()) {
            eslog::info(" == INVALID PAPI CODES  %67s == \n", info::ecf->output.papi_codes.c_str());
        }
        eslog::info(" ============================================================================================= \n");
        return;
    }

    std::vector<long long> last(papi.values), prev, begin;
    papi.read(last.data());
    std::vector<long long, initless_allocator<long long> > min(_events.size() * papi.values), max(_events.size() * papi.values), sum(_events.size() * papi.values), measurement(_values.begin(), _values.begin() + _events.size() * papi.values);

    std::vector<size_t> block;
    for (size_t i = 0; i < _events.size(); i++) {
        switch (_events[i].type) {
        case Event::START:
            block.push_back(i);
            for (int v = 0; v < papi.values; ++v) {
                measurement[papi.values * i + v] = _values[papi.values * i + v];
                prev.push_back(_values[papi.values * i + v]);
                begin.push_back(_values[papi.values * i + v]);
            }
            break;
        case Event::ACCUMULATED:
        case Event::CHECKPOINT:
            for (int v = 0; v < papi.values; ++v) {
                measurement[papi.values * i + v] = _values[papi.values * i + v] - prev[papi.values * (block.size() - 1) + v];
                prev[papi.values * (block.size() - 1) + v] = _values[papi.values * i + v];
            }
            break;
        case Event::END:
            for (int v = 0; v < papi.values; ++v) {
                measurement[papi.values * i + v] = _values[papi.values * i + v] - prev[papi.values * (block.size() - 1) + v];
            }
            break;
        case Event::DURATION:
            for (int v = 0; v < papi.values; ++v) {
                measurement[papi.values * i + v] = _values[papi.values * i + v] - begin[papi.values * (block.size() - 1) + v];
            }
            for (int v = 0; v < papi.values; ++v) {
                prev.pop_back();
                begin.pop_back();
            }
            block.pop_back();
            break;
        default:
            break;
        }
    }

    { // synchronize across processes
        size_t eventsize = _events.size(), minsize, maxsize;
        Communication::allReduce(&eventsize, &minsize, 1, MPITools::getType<size_t>().mpitype, MPI_MIN);
        Communication::allReduce(&eventsize, &maxsize, 1, MPITools::getType<size_t>().mpitype, MPI_MAX);
        if (minsize == eventsize && maxsize == eventsize) {
            Communication::reduce(measurement.data(), min.data(), measurement.size(), MPITools::getType<long long>().mpitype, MPI_MIN, 0);
            Communication::reduce(measurement.data(), max.data(), measurement.size(), MPITools::getType<long long>().mpitype, MPI_MAX, 0);
            Communication::reduce(measurement.data(), sum.data(), measurement.size(), MPITools::getType<long long>().mpitype, MPI_SUM, 0);
        } else {
            min = max = sum = measurement;
            eslog::warning("Various number of time events (only root data are printed)\n");
        }
    }

    if (grank) {
        Communication::barrier();
        return;
    }

    auto print = [&] (size_t start, size_t end, int printeddepth) {
        auto replace = [] (const char *name) {
            std::string str = name;
            for (size_t i = 0; i < str.size(); ++i) {
                if (str[i] == ' ') { str[i] = '_'; }
            }
            str.resize(40, ' ');
            return str;
        };

        int depth = printeddepth - 1;
        std::vector<const char*> printed;
        for (size_t i = start; i <= end; i++) {
            switch (_events[i].type) {
            case Event::START:
                if (++depth == printeddepth) {
                    eslog::info("  PAPI::START::%s avg: ", replace(_events[i].name).c_str());
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", sum[papi.values * i + v] / info::mpi::size);
                    }
                    eslog::info("\n%55c min: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", min[papi.values * i + v]);
                    }
                    eslog::info("\n%55c max: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", max[papi.values * i + v]);
                    }
                    eslog::info("\n\n");
                }
                break;
            case Event::CHECKPOINT:
                if (depth == printeddepth) {
                    eslog::info("  PAPI::CHECK::%s avg: ", replace(_events[i].name).c_str());
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", sum[papi.values * i + v] / info::mpi::size);
                    }
                    eslog::info("\n%55c min: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", min[papi.values * i + v]);
                    }
                    eslog::info("\n%55c max: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", max[papi.values * i + v]);
                    }
                    eslog::info("\n\n");
                }
                break;
            case Event::ACCUMULATED:
                if (depth == printeddepth) {
                    if (std::find(printed.begin(), printed.end(), _events[i].name) == printed.end()) {
                        printed.push_back(_events[i].name);
                        int counter = 1;
                        std::vector<long long> amin(papi.values), amax(papi.values), asum(papi.values);
                        for (int v = 0; v < papi.values; ++v) {
                            amin[v] = min[papi.values * i + v];
                            amax[v] = max[papi.values * i + v];
                            asum[v] = sum[papi.values * i + v];
                        }
                        for (size_t j = i + 1; j <= end; ++j) {
                            if (_events[j].name == _events[i].name) {
                                for (int v = 0; v < papi.values; ++v) {
                                    amin[v] = std::min(min[papi.values * j + v], amin[v]);
                                    amax[v] = std::max(max[papi.values * j + v], amax[v]);
                                    asum[v] += sum[papi.values * j + v];
                                }
                                ++counter;
                            }
                        }
                        eslog::info("  PAPI::ACCUM::%s avg: ", replace(_events[i].name).c_str());
                        for (int v = 0; v < papi.values; ++v) {
                            eslog::info(" %15lu", asum[v] / info::mpi::size / counter);
                        }
                        eslog::info("\n%55c min: ", ' ' );
                        for (int v = 0; v < papi.values; ++v) {
                            eslog::info(" %15lu", amin[v]);
                        }
                        eslog::info("\n%55c max: ", ' ' );
                        for (int v = 0; v < papi.values; ++v) {
                            eslog::info(" %15lu", amax[v]);
                        }
                        eslog::info("\n%55c sum: ", ' ' );
                        for (int v = 0; v < papi.values; ++v) {
                            eslog::info(" %15lu", asum[v] / info::mpi::size);
                        }
                        eslog::info("\n\n");
                    }
                }
                break;
            case Event::END:
                if (depth == printeddepth) {
                    eslog::info("  PAPI::END::::%s avg: ", replace(_events[i].name).c_str());
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", sum[papi.values * i + v] / info::mpi::size);
                    }
                    eslog::info("\n%55c min: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", min[papi.values * i + v]);
                    }
                    eslog::info("\n%55c max: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", max[papi.values * i + v]);
                    }
                    eslog::info("\n\n");
                }
                break;
            case Event::DURATION:
                if (depth-- == printeddepth) {
                    eslog::info("  PAPI::TOTAL::%s avg: ", replace(_events[start].name).c_str());
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", sum[papi.values * i + v] / info::mpi::size);
                    }
                    eslog::info("\n%55c min: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", min[papi.values * i + v]);
                    }
                    eslog::info("\n%55c max: ", ' ' );
                    for (int v = 0; v < papi.values; ++v) {
                        eslog::info(" %15lu", max[papi.values * i + v]);
                    }
                    eslog::info("\n\n");
                }
                break;
            case Event::LOADSTEP:
                break;
            default:
                break;
            }
        }
        eslog::info(" ============================================================================================= \n");
    };

    eslog::info("\n ============================================================================================= \n");

    std::vector<size_t> begins;
    for (size_t i = 0; i < _events.size(); i++) {
        switch (_events[i].type) {
        case Event::START:
            begins.push_back(i);
            break;
        case Event::DURATION:
            if (begins.size() > 1) {
                print(begins.back(), i, begins.size() + 1);
            }
            begins.pop_back();
            break;
        default:
            break;
        }
    }

    eslog::info("\n ============================================================================================= \n");
    if (info::ecf->output.papi_events.size() && papi.status == 0) {
        eslog::info(" == PAPI EVENTS %75s == \n", info::ecf->output.papi_events.c_str());
    }
    if (info::ecf->output.papi_codes.size() && papi.status == 0) {
        eslog::info(" == PAPI CODES %76s == \n", info::ecf->output.papi_codes.c_str());
    }
    eslog::info(" ============================================================================================= \n");
    print(0, _events.size() - 1, 0);
    Communication::barrier();
}
