
#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"

#include <numeric>
#include "points.h"

using namespace espreso;

bool OpenFOAMPoints::readData(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
    size_t threads = info::env::OMP_NUM_THREADS;

    std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);

    std::vector<std::vector<Point> > points(threads);

    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        std::vector<Point> tpoints;

        const char *c = begin + tdistribution[t];
        while (c < end && *c != '(') { ++c; }
        while (c < begin + tdistribution[t + 1]) {
            tpoints.push_back({});

            c += 1; // skip '('
            tpoints.back().x = scaleFactor * readDouble(c);
            tpoints.back().y = scaleFactor * readDouble(c);
            tpoints.back().z = scaleFactor * readDouble(c);
            c += 2; // skip ')\n'
        }

        points[t].swap(tpoints);
    }

    for (size_t t = 0; t < threads; t++) {
        coordinates.insert(coordinates.end(), points[t].begin(), points[t].end());
    }

    nIDs.resize(coordinates.size());

    esint offset = coordinates.size();
    Communication::exscan(offset);
    std::iota(nIDs.begin(), nIDs.end(), offset);
    return true;
}




