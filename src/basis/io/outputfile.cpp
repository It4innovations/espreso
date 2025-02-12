
#include "outputfile.h"
#include "writer.h"

#include "basis/containers/tarray.h"
#include "basis/containers/allocators.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

char OutputFile::buffer[bsize];

void OutputFile::_group()
{
    _distribution.push_back(_buffer.size());
}

OutputFilePack::~OutputFilePack()
{
    clear();
}

void OutputFilePack::clear()
{
    for (size_t i = 0; i < _files.size(); ++i) {
        delete _files[i];
    }
    _files.clear();
}

void OutputFilePack::groupData()
{
    _group();
}

void OutputFilePack::commitFile(const std::string &name)
{
    _files.push_back(new OutputFile());
    _files.back()->_name = name;
    _files.back()->_buffer.swap(_buffer);
    _files.back()->_distribution.swap(_distribution);
}

static int reorderedRank(int rank, int offset)
{
    int size = MPITools::subset->acrosssize * MPITools::subset->withinsize;
    return (rank - offset + size) % size;
}

void OutputFilePack::reorder()
{
    profiler::syncstart("reorder");
    std::vector<size_t> file(_files.size()), totalsize(_files.size()), groups(_files.size()), size, dataoffset, length, original;
    for (size_t i = 0; i < _files.size(); ++i) {
        file[i] = dataoffset.size();
        groups[i] = _files[i]->_distribution.size();
        for (size_t j = 0, prev = 0; j < _files[i]->_distribution.size(); prev = _files[i]->_distribution[j++]) {
            original.push_back(prev);
            dataoffset.push_back(_files[i]->_distribution[j] - prev);
        }
    }
    length = dataoffset;
    profiler::synccheckpoint("set_data_offset");
    Communication::exscan(size, dataoffset, MPITools::asynchronous);
    profiler::synccheckpoint("exscan");
    for (size_t i = 0; i < _files.size(); ++i) {
        for (size_t j = 0; j < _files[i]->_distribution.size(); ++j) {
            dataoffset[file[i] + j] += totalsize[i];
            totalsize[i] += size[file[i] + j];
        }
    }

    { // set a reasonable reorganization according to the stripe size setting (black magic)
        size_t stripe = info::ecf->output.stripe_size;
        size_t nwriters = MPITools::subset->acrosssize;
        size_t reduction = MPITools::subset->withinsize;
        size_t stripeoffset = 0;
        for (size_t i = 0; i < _files.size(); ++i) {
            size_t mult = totalsize[i] / (stripe * nwriters) + ((totalsize[i] % (stripe * nwriters)) ? 1 : 0);
            _files[i]->_distribution.resize(nwriters * reduction + 1, totalsize[i]);
            size_t fsize = 0, restsize = totalsize[i], restwriters = nwriters;
            _files[i]->_offset = stripeoffset;
            for (int r = 0; r < info::mpi::size; ++r) {
                if (r % reduction == 0) {
                    size_t size = std::min(mult * stripe, restsize);
                    _files[i]->_distribution[r] = fsize;
                    _files[i]->_distribution[r + 1] = fsize + size;
                    fsize += size;
                    restsize = totalsize[i] - fsize;
                    if (restsize) {
                        if (--restwriters && mult > 1) {
                            if (restsize < stripe * restwriters * (mult - 1)) {
                                --mult;
                            }
                        }
                    } else {
                        if (size) {
                            stripeoffset = (_files[i]->_offset + r + reduction >= (size_t)info::mpi::size) ? 0 : _files[i]->_offset + r + reduction;
                        }
                    }
                } else {
                    _files[i]->_distribution[r + 1] = _files[i]->_distribution[r];
                }
            }
        }
    }
    profiler::synccheckpoint("set_reorganization");

    std::vector<size_t> sBuffer, rBuffer;
    { // fill send buffer
        size_t buffers = 0;
        for (size_t i = 0; i < _files.size(); ++i) {
            buffers += _files[i]->_buffer.size() / sizeof(size_t) + sizeof(size_t);
        }
        sBuffer.reserve(_files.size() * 2 * info::mpi::size + 2 * info::mpi::size + buffers); // size, target, offset, size

        std::vector<size_t> ii(_files.size());
        for (int r = 0; r < info::mpi::size; ++r) {
            size_t prevsize = sBuffer.size();
            sBuffer.push_back(0); // total size
            sBuffer.push_back(r); // target

            if (r % MPITools::subset->withinsize == 0) { // target is writer
                for (size_t j = 0; j < _files.size(); ++j) {
                    size_t begin = _files[j]->_distribution[reorderedRank(r, _files[j]->_offset)];
                    size_t end = _files[j]->_distribution[reorderedRank(r, _files[j]->_offset) + 1];
                    ii[j] = 0;
                    while (ii[j] < groups[j] && dataoffset[file[j] + ii[j]] + length[file[j] + ii[j]] < begin) {
                        ++ii[j];
                    }

                    while (ii[j] < groups[j] && dataoffset[file[j] + ii[j]] < end) {
                        size_t databegin = std::max(dataoffset[file[j] + ii[j]], begin);
                        size_t dataend = std::min(dataoffset[file[j] + ii[j]] + length[file[j] + ii[j]], end);
                        size_t currentoffset = databegin - dataoffset[file[j] + ii[j]];
                        size_t datasize = dataend - databegin;
                        size_t bsize = datasize / sizeof(size_t) + ((datasize % sizeof(size_t)) ? 1 : 0);

                        sBuffer.push_back(j);
                        sBuffer.push_back(databegin);
                        sBuffer.push_back(datasize);
                        sBuffer.insert(sBuffer.end(),
                                reinterpret_cast<size_t*>(_files[j]->_buffer.data() + original[file[j] + ii[j]] + currentoffset),
                                reinterpret_cast<size_t*>(_files[j]->_buffer.data() + original[file[j] + ii[j]] + currentoffset) + bsize);

                        ++ii[j];
                    }
                }
            }
            sBuffer[prevsize] = sBuffer.size() - prevsize;
        }
    }

    profiler::synccheckpoint("sbuffer");

    if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer, 0, info::mpi::size, MPITools::asynchronous)) {
        eslog::internalFailure("cannot distribute permuted nodes.\n");
    }
    profiler::synccheckpoint("all_to_all");


    if (MPITools::subset->within.rank == 0) {
        for (size_t j = 0; j < _files.size(); ++j) {
            int rr = reorderedRank(info::mpi::rank, _files[j]->_offset);
            _files[j]->_buffer.resize(_files[j]->_distribution[rr + 1] - _files[j]->_distribution[rr]);
        }
        size_t offset = 0;
        for (int r = 0; r < info::mpi::size; r++) {
            size_t total = offset + rBuffer[offset];
            ++offset; // size
            ++offset; // me

            while (offset < total) {
                size_t j = rBuffer[offset++];
                size_t displacement = rBuffer[offset++];
                size_t size = rBuffer[offset++];
                memcpy(_files[j]->_buffer.data() + displacement - _files[j]->_distribution[reorderedRank(info::mpi::rank, _files[j]->_offset)], rBuffer.data() + offset, size);
                offset += size / sizeof(size_t) + ((size % sizeof(size_t)) ? 1 : 0);
            }
        }
    }
    profiler::synccheckpoint("rbuffer");
    profiler::syncend("reorder");
}

void OutputFilePack::write()
{
    profiler::syncstart("write");
    switch (info::ecf->output.writer) {
    case OutputConfiguration::WRITER::POSIX: eslog::internalFailure("POSIX writer does not work.\n"); break;
    case OutputConfiguration::WRITER::MPI:
    case OutputConfiguration::WRITER::MPI_COLLECTIVE:
        break;
    }

    for (size_t i = 0; i < _files.size(); ++i) {
        size_t chunk = _files[i]->_distribution[1];
        size_t chunkmax = INT32_MAX;
        size_t chunks = chunk / chunkmax + ((chunk % chunkmax) ? 1 : 0);
        size_t chunkoffset = _files[i]->_distribution[reorderedRank(info::mpi::rank, _files[i]->_offset)];
        size_t chunksize = _files[i]->_distribution[reorderedRank(info::mpi::rank, _files[i]->_offset) + 1] - chunkoffset;

        if (info::ecf->output.writer == OutputConfiguration::WRITER::MPI) {
            MPIWriter writer;
            if (MPITools::subset->within.rank == 0 && chunksize) {
                if (writer.open(MPITools::subset->across, _files[i]->_name)) {
                    eslog::error("WRITER: cannot create file '%s'\n", _files[i]->_name.c_str());
                }

                for (size_t c = 0; c < chunks; ++c) {
                    size_t size = std::min(chunkmax, _files[i]->_buffer.size() - c * chunkmax);
                    writer.store(_files[i]->_buffer.data() + chunkoffset - _files[i]->_distribution[reorderedRank(info::mpi::rank, _files[i]->_offset)], chunkoffset, size);
                    chunkoffset += size;
                }
                writer.close();
            }
        }
        if (info::ecf->output.writer == OutputConfiguration::WRITER::MPI_COLLECTIVE) {
            MPICollectiveWriter writer;
            if (MPITools::subset->within.rank == 0) {
                if (writer.open(MPITools::subset->across, _files[i]->_name)) {
                    eslog::error("WRITER: cannot create file '%s'\n", _files[i]->_name.c_str());
                }

                for (size_t c = 0; c < chunks; ++c) {
                    size_t size = std::min(chunkmax, _files[i]->_buffer.size() - c * chunkmax);
                    writer.store(_files[i]->_buffer.data() + chunkoffset - _files[i]->_distribution[reorderedRank(info::mpi::rank, _files[i]->_offset)], chunkoffset, size);
                    chunkoffset += size;
                }
                writer.close();
            }
        }
        profiler::synccheckpoint("write");
        profiler::syncparam("filesize", chunksize);
    }

    clear();
    profiler::syncend("write");
}
