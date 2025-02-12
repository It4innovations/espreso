
#include "faces.h"

#include "input/parsers/openfoam/openfoam.h"

#include "basis/containers/tarray.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"

#include <numeric>

using namespace espreso;

bool OpenFOAMFaces::readFaces(OpenFOAMData &data)
{
    size_t threads = info::env::OMP_NUM_THREADS;

    std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);

    std::vector<std::vector<esint> > fsize(threads), fnodes(threads);

    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        std::vector<esint> tsize, tnodes;

        const char *c = begin + tdistribution[t];
        while (c < end && *(c + 1) != '(') { ++c; }
        while (c < begin + tdistribution[t + 1]) {
            tsize.push_back(readInteger(c));
            c += 1; // skip '('
            for (esint i = 0; i < tsize.back(); ++i) {
                tnodes.push_back(readInteger(c));
            }
            c += 2; // skip ')\n'
        }

        fsize[t].swap(tsize);
        fnodes[t].swap(tnodes);
    }

    for (size_t t = 0; t < threads; t++) {
        data.fsize.insert(data.fsize.end(), fsize[t].begin(), fsize[t].end());
        data.fnodes.insert(data.fnodes.end(), fnodes[t].begin(), fnodes[t].end());
    }

    size_t offset = data.fsize.size();
    Communication::exscan(offset);
    data.fIDs.resize(data.fsize.size());
    std::iota(data.fIDs.begin(), data.fIDs.end(), offset);

    return true;
}

bool OpenFOAMFaces::readParents(std::vector<esint> &data)
{
    size_t threads = info::env::OMP_NUM_THREADS;

    std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, end - begin);

    std::vector<std::vector<esint> > elements(threads);

    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        std::vector<esint> telements;

        const char *c = begin + tdistribution[t];
        if (tdistribution[t]) {
            --c; toNext(c);
        }
        while (c < begin + tdistribution[t + 1]) {
            telements.push_back(readInteger(c));
            c += 1; // skip '\n'
        }

        elements[t].swap(telements);
    }

    for (size_t t = 0; t < threads; t++) {
        data.insert(data.end(), elements[t].begin(), elements[t].end());
    }

    return true;
}


