
#include "profiler.h"
#include "omp.h"

#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"

#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <fstream>

using namespace espreso;
using namespace espreso::profiler;

double Profiler<active()>::initted = Profiler<active()>::time();
Profiler<active()> Profiler<active()>::instance;

double Profiler<active()>::time()
{
    return omp_get_wtime();
}

Profiler<true>::Profiler()
: level(0), size(0), events(NULL)
{
    // be aware of the size of events
    events = new Event<true>[5000000];
}

Profiler<true>::~Profiler()
{
    delete[] events;
}

//static void gatherNameMap(Event<true>* events, int size, std::map<std::string, size_t> &n2index)
//{
//    std::vector<char> namemap;
//    if (info::mpi::rank == 0) {
//        for (int i = 0; i < size; ++i) {
//            size_t offset = (size_t)events[i].name;
//            int size = strlen(events[i].name);
//            namemap.insert(namemap.end(), reinterpret_cast<char*>(&offset), reinterpret_cast<char*>(&offset) + sizeof(offset));
//            namemap.insert(namemap.end(), reinterpret_cast<char*>(&size), reinterpret_cast<char*>(&size) + sizeof(size));
//            namemap.insert(namemap.end(), events[i].name, events[i].name + size);
//        }
//    }
//    Communication::broadcastUnknownSize(namemap);
//    for (size_t j = 0; j < namemap.size();) {
//        size_t offset = *reinterpret_cast<size_t*>(namemap.data() + j);
//        int size = *reinterpret_cast<int*>(namemap.data() + j + sizeof(size_t));
//        std::string name(namemap.data() + j + sizeof(size_t) + sizeof(int), namemap.data() + j + sizeof(size_t) + sizeof(int) + size);
//        n2index[name] = offset;
//        j += sizeof(size_t) + sizeof(int) + size;
//    }
//    namemap.clear();
//    for (int i = 0; i < size; ++i) {
//        if (events[i].type & Event<true>::TYPE::START) {
//            std::string name(events[i].name);
//            if (n2index.find(name) == n2index.end()) {
//                size_t offset = (size_t)events[i].name;
//                int size = strlen(events[i].name);
//                namemap.insert(namemap.end(), reinterpret_cast<char*>(&offset), reinterpret_cast<char*>(&offset) + sizeof(offset));
//                namemap.insert(namemap.end(), reinterpret_cast<char*>(&size), reinterpret_cast<char*>(&size) + sizeof(size));
//                namemap.insert(namemap.end(), events[i].name, events[i].name + size);
//            }
//        }
//    }
//    Communication::allGatherUnknownSize(namemap);
//    for (size_t j = 0; j < namemap.size(); ) {
//        size_t offset = *reinterpret_cast<size_t*>(namemap.data() + j);
//        int size = *reinterpret_cast<int*>(namemap.data() + j + sizeof(size_t));
//        std::string name(namemap.data() + j + sizeof(size_t) + sizeof(int), namemap.data() + j + sizeof(size_t) + sizeof(int) + size);
//        if (n2index.find(name) == n2index.end()) {
//            n2index[name] = offset;
//        }
//        j += sizeof(size_t) + sizeof(int) + size;
//    }
//}
//
//static void gatherIDs(Event<true>* events, int size, std::map<std::string, size_t> &n2index, std::vector<int> &asyncEvents, std::vector<std::vector<size_t> > &ids)
//{
//    std::vector<int> level, async, maxasync;
//    std::vector<size_t> offsets;
//    offsets.reserve(size);
//    for (int i = 0; i < size; ++i) {
//        if (events[i].type & Event<true>::TYPE::START) {
//            level.push_back(i);
//        }
//        offsets.push_back((size_t)(events[i].name - events[level.back()].name + n2index[events[level.back()].name]));
//        if (events[i].type & Event<true>::TYPE::SYNC) {
//            async.push_back(0);
//        } else {
//            if (async.back() == 0) {
//                ids.push_back({});
//            }
//            ++async.back();
//            ids.back().push_back(offsets.back());
//        }
//        if (events[i].type & Event<true>::TYPE::END) {
//            level.pop_back();
//        }
//    }
//
//    for (size_t i = 0; i < ids.size(); ++i) {
//        utils::sortAndRemoveDuplicates(ids[i]);
//    }
//
//    maxasync.resize(async.size());
//    Communication::allReduce(async.data(), maxasync.data(), async.size(), MPI_INT, MPI_MAX);
//    std::vector<int> asyncpoints;
//    for (size_t i = 0; i < maxasync.size(); ++i) {
//        if (maxasync[i]) {
//            if (async[i] == 0) {
//                ids.insert(ids.begin() + asyncpoints.size(), std::vector<size_t>());
//            }
//            asyncpoints.push_back(i);
//        }
//    }
//
//    std::vector<int> ssize(asyncpoints.size()), rsize(asyncpoints.size());
//    for (size_t i = 0; i < asyncpoints.size(); ++i) {
//        if (maxasync[asyncpoints[i]] && async[asyncpoints[i]]) {
//            ssize[i] = ids[i].size();
//        }
//    }
//    Communication::allReduce(ssize.data(), rsize.data(), asyncpoints.size(), MPI_INT, MPI_MAX);
//
//    std::vector<size_t> sids, rids;
//    for (size_t i = 0; i < asyncpoints.size(); ++i) {
//        if (ssize[i] == rsize[i]) {
//            sids.insert(sids.end(), ids[i].begin(), ids[i].end());
//        } else {
//            sids.insert(sids.end(), rsize[i], 0);
//        }
//    }
//    rids.resize(sids.size());
//    Communication::allReduce(sids.data(), rids.data(), sids.size(), MPITools::getType<size_t>().mpitype, MPI_MAX);
//
//    std::vector<int> scheck(asyncpoints.size()), rcheck(asyncpoints.size());
//    for (size_t i = 0, offset = 0; i < asyncpoints.size(); ++i) {
//        scheck[i] = std::includes(rids.begin() + offset, rids.begin() + offset + rsize[i], ids[i].begin(), ids[i].end());
//        offset += rsize[i];
//    }
//    Communication::allReduce(scheck.data(), rcheck.data(), asyncpoints.size(), MPI_INT, MPI_MIN);
//
//    for (size_t i = 0, offset = 0; i < rcheck.size(); ++i) {
//        if (rcheck[i] == 0) {
//            std::vector<size_t> tmp(ids[i].size());
//            auto it = std::set_difference(ids[i].begin(), ids[i].end(), rids.begin() + offset, rids.begin() + offset + rsize[i], tmp.begin());
//            tmp.resize(it - tmp.begin());
//            Communication::allGatherUnknownSize(tmp);
//            ids[i].assign(rids.begin() + offset, rids.begin() + offset + rsize[i]);
//            ids[i].insert(ids[i].end(), tmp.begin(), tmp.end());
//            utils::sortAndRemoveDuplicates(ids[i]);
//        } else {
//            ids[i].assign(rids.begin() + offset, rids.begin() + offset + rsize[i]);
//        }
//        offset += rsize[i];
//    }
//}

void Profiler<true>::print()
{
    return; // TODO
//    int finalsize = size;
//    std::map<std::string, size_t> n2index;
//    std::vector<int> asyncEvents, eventCounter;
//    std::vector<std::vector<size_t> > ids;
//
//    gatherNameMap(events, finalsize, n2index);
//    gatherIDs(events, finalsize, n2index, asyncEvents, ids);

}
