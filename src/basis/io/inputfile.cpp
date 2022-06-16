
#include "inputfile.h"
#include "loader.h"

#include "basis/containers/tarray.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include <cstring>

using namespace espreso;

InputFile::InputFile()
: begin(NULL), end(NULL), hardend(NULL), maxchunk(0)
{

}

InputFilePack::InputFilePack(size_t minchunk, size_t overlap)
: fileindex((size_t)-1), minchunk(minchunk), overlap(overlap)
{

}

InputFilePack::~InputFilePack()
{
	clear();
}

void InputFilePack::clear()
{
        for (size_t i = 0; i < files.size(); ++i) {
                delete files[i];
        }
	files.clear();
}

void InputFilePack::commitFiles(const std::vector<std::string> &filepaths)
{
	paths = filepaths;
	for (size_t i = 0; i < paths.size(); ++i) {
		files.push_back(new InputFile());
	}
}

bool InputFilePack::next()
{
	if (fileindex != (size_t)-1) {
		files[fileindex]->begin = begin;
		files[fileindex]->end = end;
		files[fileindex]->hardend = hardend;
		files[fileindex]->distribution = distribution;
	}
	if (++fileindex < files.size()) {
		begin = files[fileindex]->begin;
		end = files[fileindex]->end;
		hardend = files[fileindex]->hardend;
		distribution = files[fileindex]->distribution;
		return true;
	}
	fileindex = (size_t)-1;
	return false;
}

void InputFilePack::prepare()
{
	profiler::syncstart("inputpack_prepare");
	std::vector<size_t> totalsize(paths.size());
	if (info::mpi::rank == 0) {
		for (size_t i = 0; i < paths.size(); ++i) {
			MPILoader loader;
			if (!loader.open(MPITools::subset->across, paths[i])) {
				totalsize[i] = loader.size();
				loader.close();
			} else {
				eslog::error("MESIO: file '%s' cannot be opened.\n", paths[i].c_str());
			}
		}
	}
	profiler::synccheckpoint("get_size");
	Communication::broadcast(totalsize.data(), totalsize.size(), MPITools::getType<size_t>().mpitype, 0);
	profiler::synccheckpoint("broadcast");

	{ // set a reasonable reorganization according to the stripe size setting (black magic)
		size_t stripe = info::ecf->input.stripe_size;
		size_t nloaders = MPITools::subset->acrosssize;
		size_t reduction = MPITools::subset->withinsize;
		size_t stripeoffset = 0;

		// increase the stripe size in order to read requested minimal amount of data for all processes
		if (stripe < minchunk * reduction) {
			stripe *= minchunk * reduction / stripe + ((minchunk * reduction % stripe) ? 1 : 0);
		}

		for (size_t i = 0; i < files.size(); ++i) {
			size_t mult = totalsize[i] / (stripe * nloaders) + ((totalsize[i] % (stripe * nloaders)) ? 1 : 0);
			files[i]->maxchunk = mult * stripe;
			files[i]->distribution.resize(info::mpi::size + 1);

			size_t fsize = 0, restsize = totalsize[i], restloaders = nloaders;
			int firstr = 0;
			if (totalsize[i] < stripe * nloaders / 2) { // file is read lesser than half of loaders
				size_t freesize = stripe * (nloaders - stripeoffset);
				size_t nstripes = totalsize[i] / stripe + ((totalsize[i] % stripe) ? 1 : 0);
				if (freesize < totalsize[i]) {
					if (stripeoffset < nloaders && nstripes <= nloaders) {
						stripeoffset = nloaders - nstripes;
					} else {
						stripeoffset = 0;
					}
				}
				if (nloaders <= stripeoffset) {
					stripeoffset = 0;
				}
				firstr = stripeoffset * reduction;
				stripeoffset += nstripes;
			}
			for (int r = firstr; r < info::mpi::size; ) {
				size_t maxreduction = std::min((int)reduction, info::mpi::size - r);
				size_t size = std::min(mult * stripe, restsize);
				size_t stepsize = std::max(size / maxreduction, std::min(size, minchunk));
				for (size_t rr = 0; rr < maxreduction; ++r, ++rr) {
					files[i]->distribution[r] = fsize + std::min(rr * stepsize, size);
					files[i]->distribution[r + 1] = fsize + std::min((rr + 1) * stepsize, size);
				}
				files[i]->distribution[r] = fsize + size; // fix the situation if size % reduction != 0
				fsize += size;
				restsize = totalsize[i] - fsize;
				if (restsize) {
					if (--restloaders && mult > 1) {
						if (restsize < stripe * restloaders * (mult - 1)) {
							--mult;
						}
					}
				}
			}

			if (MPITools::subset->within.rank == 0) {
				files[i]->data.reserve(files[i]->distribution[info::mpi::rank + MPITools::subset->within.size] - files[i]->distribution[info::mpi::rank] + overlap);
				files[i]->data.resize(files[i]->distribution[info::mpi::rank + MPITools::subset->within.size] - files[i]->distribution[info::mpi::rank]);
				files[i]->data.resize(files[i]->distribution[info::mpi::rank + MPITools::subset->within.size] - files[i]->distribution[info::mpi::rank] + overlap, 0);
			}
		}
	}
	profiler::syncend("inputpack_prepare");
}

void InputFilePack::read()
{
	eslog::startln("READER: STARTED", "READER");

	profiler::syncstart("inputpack_read");
	for (size_t i = 0; i < files.size(); ++i) {
		size_t chunk = files[i]->maxchunk;
		size_t chunkmax = INT32_MAX;
		size_t chunks = chunk / chunkmax + ((chunk % chunkmax) ? 1 : 0);
		size_t chunkoffset = files[i]->distribution[info::mpi::rank];
		size_t chunksize = files[i]->distribution[std::min(info::mpi::rank + MPITools::subset->within.size, info::mpi::size)] - chunkoffset;

		if (info::ecf->input.loader == InputConfiguration::LOADER::MPI) {
			MPILoader loader;
			if (MPITools::subset->within.rank == 0 && chunksize) {
				if (loader.open(MPITools::subset->across, paths[i])) {
					eslog::error("LOADER: cannot read file '%s'\n", paths[i].c_str());
				}
			}

			if (MPITools::subset->within.rank == 0 && chunksize) {
				for (size_t c = 0; c < chunks; ++c) {
					size_t size = std::min(chunkmax, files[i]->data.size() - overlap - c * chunkmax);
					loader.read(files[i]->data.data() + chunkoffset - files[i]->distribution[info::mpi::rank], chunkoffset, size);
					chunkoffset += size;
				}
				loader.close();
			}
		}
		if (info::ecf->input.loader == InputConfiguration::LOADER::MPI_COLLECTIVE) {
			MPICollectiveLoader loader;
			if (MPITools::subset->within.rank == 0) {
				if (loader.open(MPITools::subset->across, paths[i])) {
					eslog::error("LOADER: cannot read file '%s'\n", paths[i].c_str());
				}
			}

			if (MPITools::subset->within.rank == 0) {
				for (size_t c = 0; c < chunks; ++c) {
					size_t size = std::min(chunkmax, files[i]->data.size() - overlap - c * chunkmax);
					loader.read(files[i]->data.data() + chunkoffset - files[i]->distribution[info::mpi::rank], chunkoffset, size);
					chunkoffset += size;
				}
				loader.close();
			}
		}
		if (info::ecf->input.loader == InputConfiguration::LOADER::POSIX) {
			POSIXLoader loader;
			if (MPITools::subset->within.rank == 0) {
				if (loader.open(MPITools::subset->across, paths[i])) {
					eslog::error("LOADER: cannot read file '%s'\n", paths[i].c_str());
				}
			}

			if (MPITools::subset->within.rank == 0) {
				loader.read(files[i]->data.data(), chunkoffset, chunksize);
				loader.close();
			}
		}
		profiler::synccheckpoint("read");
		profiler::syncparam("filesize", chunksize);
	}
	profiler::syncend("inputpack_read");
	eslog::checkpointln("READER: FILE READ");

	profiler::syncstart("inputpack_postprocess");
	int reduction = MPITools::subset->within.size;

	if (reduction > 1) { // scatter data to non-readers
		std::vector<size_t> displacement(reduction + 1);
		std::vector<char> sBuffer, rBuffer;

		int writer = info::mpi::rank - MPITools::subset->within.rank;
		size_t totalsize = 0;
		for (size_t i = 0; i < files.size(); ++i) {
			totalsize += files[i]->distribution[writer + reduction] - files[i]->distribution[writer];
		}
		if (MPITools::subset->within.rank == 0) {
			sBuffer.reserve(totalsize);
		}
		for (int r = writer; r < writer + reduction; ++r) {
			for (size_t i = 0; i < files.size(); ++i) {
				size_t chunkoffset = files[i]->distribution[r] - files[i]->distribution[writer];
				size_t chunksize = files[i]->distribution[r + 1] - files[i]->distribution[r];
				displacement[r - writer + 1] += chunkoffset + chunksize;
				if (MPITools::subset->within.rank == 0 && chunksize) {
					sBuffer.insert(sBuffer.end(), files[i]->data.begin() + chunkoffset, files[i]->data.begin() + chunkoffset + chunksize);
				}
			}
		}
		profiler::synccheckpoint("sbuffer_data");

		Communication::scatterv(sBuffer, rBuffer, displacement, &MPITools::subset->within);
		profiler::synccheckpoint("scatterv");

		for (size_t i = 0, roffset = 0; i < files.size(); ++i) {
			size_t chunksize = files[i]->distribution[info::mpi::rank + 1] - files[i]->distribution[info::mpi::rank];
			files[i]->data.reserve(chunksize + overlap);
			files[i]->data.resize(chunksize);
			files[i]->data.resize(chunksize + overlap, 0); // fix valgrind warning since the initless allocator is used
			if (chunksize) {
				memcpy(files[i]->data.data(), rBuffer.data() + roffset, chunksize);
			}
			roffset += chunksize;
		}
		profiler::synccheckpoint("rbuffer_data");
		eslog::checkpointln("READER: DATA SCATTERED");
	}

	Communication::barrier();
	eslog::checkpointln("READER: SYNCHRONIZED");

	for (size_t i = 0; i < files.size(); ++i) {
		files[i]->begin = files[i]->data.data();
		files[i]->hardend = files[i]->begin + files[i]->data.size();
		files[i]->end = files[i]->hardend - overlap;
	}

	{ // overlap data
		std::vector<char, initless_allocator<char> > sBuffer(files.size() * overlap), rBuffer(files.size() * overlap);
		for (size_t i = 0; i < files.size(); ++i) {
			memcpy(sBuffer.data() + i * overlap, files[i]->data.data(), overlap);
		}
		profiler::synccheckpoint("sbuffer_overlap");
		Communication::receiveUpper(sBuffer, rBuffer);
		profiler::synccheckpoint("receive_upper");
		for (size_t i = 0; i < files.size(); ++i) {
			if (info::mpi::rank + 1 != info::mpi::size) {
				memcpy(files[i]->data.data() + files[i]->data.size() - overlap, rBuffer.data() + i * overlap, overlap);
			}
		}
		profiler::synccheckpoint("rbuffer_overlap");
	}
	profiler::syncend("inputpack_postprocess");
	eslog::endln("READER: ALIGNED");
}

void Metadata::read(const std::string &filename)
{
	profiler::syncstart("metadata_read");
	distribution.resize(2);
	if (info::mpi::rank == 0) {
		POSIXLoader loader;
		if (loader.open(MPITools::singleton->across, filename)) {
			eslog::error("MESIO error: cannot load metadata file '%s'\n", filename.c_str());
		}
		distribution = { 0, loader.size() };
		data.resize(distribution.back());
		profiler::synccheckpoint("prepare");
		loader.read(data.data(), 0, distribution.back());
		loader.close();
		profiler::synccheckpoint("read");
	}

	Communication::broadcast(&distribution.back(), sizeof(size_t), MPI_BYTE, 0);
	profiler::synccheckpoint("broadcast_size");
	data.resize(distribution.back());
	profiler::synccheckpoint("resize");
	Communication::broadcast(data.data(), data.size(), MPI_BYTE, 0);
	profiler::synccheckpoint("broadcast_data");
	begin = data.data();
	end = data.data() + data.size();
	hardend = end;
	profiler::syncend("metadata_read");
}
