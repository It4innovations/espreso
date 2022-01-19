
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
: begin(nullptr), end(nullptr), hardend(nullptr),
  overlap(0), totalSize(0), name{},
  maxchunk(0), loader(nullptr)
{

}

InputFile::InputFile(const std::string &name, size_t overlap)
: begin(nullptr), end(nullptr), hardend(nullptr),
  overlap(overlap), totalSize(0), name(name),
  maxchunk(0), loader(nullptr)
{

}

InputFile::~InputFile()
{
	if (loader) delete loader;
}

size_t InputFile::size(const std::string &file)
{
	MPILoader loader;
	if (!loader.open(MPITools::subset->across, file)) {
		return loader.size();
	}
	eslog::error("MESIO: file '%s' cannot be opened.\n", file.c_str());
	return 0;
}

void InputFile::setDistribution(const std::vector<size_t> &distribution)
{
	this->distribution = distribution;
	this->distribution.back() = totalSize;
	data.reserve(this->distribution[info::mpi::rank + MPITools::subset->within.size] - this->distribution[info::mpi::rank] + overlap);
	data.resize(this->distribution[info::mpi::rank + MPITools::subset->within.size] - this->distribution[info::mpi::rank]);
	data.resize(this->distribution[info::mpi::rank + MPITools::subset->within.size] - this->distribution[info::mpi::rank] + overlap, 0);
}

FilePack::FilePack(size_t minchunk, size_t overlap)
: minchunk(minchunk), overlap(overlap), fileindex(-1)
{

}

FilePack::FilePack(const std::vector<std::string> &filepaths, size_t minchunk, size_t overlap)
: minchunk(minchunk), overlap(overlap), fileindex(-1)
{
	for (size_t i = 0; i < filepaths.size(); ++i) {
		files.push_back(new InputFile(filepaths[i], overlap));
	}
}

FilePack::~FilePack()
{
	for (size_t i = 0; i < files.size(); ++i) {
		delete files[i];
	}
	files.clear();
}

InputFile* FilePack::add(const std::string &name)
{
	files.push_back(new InputFile(name, overlap));
	return files.back();
}

bool FilePack::next()
{
	if (fileindex != (size_t)-1) {
		this->swap(files[fileindex]);
	}
	if (++fileindex < files.size()) {
		this->swap(files[fileindex]);
		return true;
	}
	fileindex = -1;
	return false;
}

void FilePack::setTotalSizes()
{
	std::vector<size_t> totalsize(files.size());
	if (info::mpi::rank == 0) {
		while (next()) {
			totalsize[fileindex] = size(name);
		}
	}
	Communication::broadcast(totalsize.data(), totalsize.size(), MPITools::getType<size_t>().mpitype, 0);
	while (next()) {
		totalSize = totalsize[fileindex];
	}
}

InputFilePack::InputFilePack(size_t minchunk, size_t overlap)
: FilePack(minchunk, overlap)
{

}

InputFilePack::InputFilePack(const std::vector<std::string> &filepaths, size_t minchunk, size_t overlap)
: FilePack(filepaths, minchunk, overlap)
{

}

void InputFilePack::prepare()
{
	profiler::syncstart("inputpack_prepare");
	setTotalSizes();
	profiler::synccheckpoint("get_size");

	{ // set a reasonable reorganization according to the stripe size setting (black magic)
		size_t stripe = info::ecf->input.stripe_size;
		size_t nloaders = MPITools::subset->acrosssize;
		size_t reduction = MPITools::subset->withinsize;
		size_t stripeoffset = 0;

		// increase the stripe size in order to read requested minimal amount of data for all processes
		if (stripe < minchunk * reduction) {
			stripe *= minchunk * reduction / stripe + ((minchunk * reduction % stripe) ? 1 : 0);
		}

		while (next()) {
			size_t mult = totalSize / (stripe * nloaders) + ((totalSize % (stripe * nloaders)) ? 1 : 0);
			maxchunk = mult * stripe;
			distribution.resize(info::mpi::size + 1);

			size_t fsize = 0, restsize = totalSize, restloaders = nloaders;
			int firstr = 0;
			if (totalSize < stripe * nloaders / 2) { // file is read lesser than half of loaders
				size_t freesize = stripe * (nloaders - stripeoffset);
				size_t nstripes = totalSize / stripe + ((totalSize % stripe) ? 1 : 0);
				if (freesize < totalSize) {
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
					distribution[r] = fsize + std::min(rr * stepsize, size);
					distribution[r + 1] = fsize + std::min((rr + 1) * stepsize, size);
				}
				distribution[r] = fsize + size; // fix the situation if size % reduction != 0
				fsize += size;
				restsize = totalSize - fsize;
				if (restsize) {
					if (--restloaders && mult > 1) {
						if (restsize < stripe * restloaders * (mult - 1)) {
							--mult;
						}
					}
				}
			}

			if (MPITools::subset->within.rank == 0) {
				data.reserve(distribution[info::mpi::rank + MPITools::subset->within.size] - distribution[info::mpi::rank] + overlap);
				data.resize(distribution[info::mpi::rank + MPITools::subset->within.size] - distribution[info::mpi::rank]);
				data.resize(distribution[info::mpi::rank + MPITools::subset->within.size] - distribution[info::mpi::rank] + overlap, 0);
			}
		}
	}
	profiler::syncend("inputpack_prepare");
}

void InputFilePack::read()
{
	eslog::startln("READER: STARTED", "READER");

	profiler::syncstart("inputpack_read");
	while (next()) {
		size_t chunk = maxchunk;
		size_t chunkmax = INT32_MAX;
		size_t chunks = chunk / chunkmax + ((chunk % chunkmax) ? 1 : 0);
		size_t chunkoffset = distribution[info::mpi::rank];
		size_t chunksize = distribution[std::min(info::mpi::rank + MPITools::subset->within.size, info::mpi::size)] - chunkoffset;

		if (info::ecf->input.loader == InputConfiguration::LOADER::MPI) {
			MPILoader loader;
			if (MPITools::subset->within.rank == 0 && chunksize) {
				if (loader.open(MPITools::subset->across, name)) {
					eslog::error("LOADER: cannot read file '%s'\n", name.c_str());
				}
				for (size_t c = 0; c < chunks; ++c) {
					size_t size = std::min(chunkmax, data.size() - overlap - c * chunkmax);
					loader.read(data.data() + chunkoffset - distribution[info::mpi::rank], chunkoffset, size);
					chunkoffset += size;
				}
			}
		}
		if (info::ecf->input.loader == InputConfiguration::LOADER::MPI_COLLECTIVE) {
			MPICollectiveLoader loader;
			if (MPITools::subset->within.rank == 0) {
				if (loader.open(MPITools::subset->across, name)) {
					eslog::error("LOADER: cannot read file '%s'\n", name.c_str());
				}
				for (size_t c = 0; c < chunks; ++c) {
					size_t size = std::min(chunkmax, data.size() - overlap - c * chunkmax);
					loader.read(data.data() + chunkoffset - distribution[info::mpi::rank], chunkoffset, size);
					chunkoffset += size;
				}
			}
		}
		if (info::ecf->input.loader == InputConfiguration::LOADER::POSIX) {
			POSIXLoader loader;
			if (MPITools::subset->within.rank == 0) {
				if (loader.open(MPITools::subset->across, name)) {
					eslog::error("LOADER: cannot read file '%s'\n", name.c_str());
				}
				loader.read(data.data(), chunkoffset, chunksize);
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
		while (next()) {
			totalsize += distribution[writer + reduction] - distribution[writer];
		}
		if (MPITools::subset->within.rank == 0) {
			sBuffer.reserve(totalsize);
		}
		for (int r = writer; r < writer + reduction; ++r) {
			while(next()) {
				size_t chunkoffset = distribution[r] - distribution[writer];
				size_t chunksize = distribution[r + 1] - distribution[r];
				displacement[r - writer + 1] += chunkoffset + chunksize;
				if (MPITools::subset->within.rank == 0 && chunksize) {
					sBuffer.insert(sBuffer.end(), data.begin() + chunkoffset, data.begin() + chunkoffset + chunksize);
				}
			}
		}
		profiler::synccheckpoint("sbuffer_data");

		Communication::scatterv(sBuffer, rBuffer, displacement, &MPITools::subset->within);
		profiler::synccheckpoint("scatterv");

		size_t roffset = 0;
		while (next()) {
			size_t chunksize = distribution[info::mpi::rank + 1] - distribution[info::mpi::rank];
			data.reserve(chunksize + overlap);
			data.resize(chunksize);
			data.resize(chunksize + overlap, 0); // fix valgrind warning since the initless allocator is used
			if (chunksize) {
				memcpy(data.data(), rBuffer.data() + roffset, chunksize);
			}
			roffset += chunksize;
		}
		profiler::synccheckpoint("rbuffer_data");
		eslog::checkpointln("READER: DATA SCATTERED");
	}

	Communication::barrier();
	eslog::checkpointln("READER: SYNCHRONIZED");

	while (next()) {
		begin = data.data();
		hardend = begin + data.size();
		end = hardend - overlap;
	}

	{ // overlap data
		std::vector<char, initless_allocator<char> > sBuffer(files.size() * overlap), rBuffer(files.size() * overlap);
		while (next()) {
			memcpy(sBuffer.data() + fileindex * overlap, data.data(), overlap);
		}
		profiler::synccheckpoint("sbuffer_overlap");
		Communication::receiveUpper(sBuffer, rBuffer);
		profiler::synccheckpoint("receive_upper");
		while (next()) {
			if (info::mpi::rank + 1 != info::mpi::size) {
				memcpy(data.data() + data.size() - overlap, rBuffer.data() + fileindex * overlap, overlap);
			}
		}
		profiler::synccheckpoint("rbuffer_overlap");
	}
	profiler::syncend("inputpack_postprocess");
	eslog::endln("READER: ALIGNED");
}

AsyncFilePack::AsyncFilePack(size_t overlap)
: FilePack(0, overlap)
{

}

AsyncFilePack::AsyncFilePack(const std::vector<std::string> &filepaths, size_t overlap)
: FilePack(filepaths, 0, overlap)
{

}

void AsyncFilePack::iread(std::function<void(void)> callback)
{
	this->callback = callback;
	while (next()) {
		size_t chunkoffset = distribution[info::mpi::rank];
		size_t chunksize = distribution[std::min(info::mpi::rank + MPITools::subset->within.size, info::mpi::size)] - chunkoffset;
		size_t chunkmax = INT32_MAX;
		size_t chunks = chunksize / chunkmax + ((chunksize % chunkmax) ? 1 : 0);

		switch (info::ecf->input.loader) {
		case InputConfiguration::LOADER::MPI_COLLECTIVE: loader = new MPICollectiveLoader(); break;
		case InputConfiguration::LOADER::MPI: loader = new MPILoader(); break;
		case InputConfiguration::LOADER::POSIX: loader = new MPILoader(); break; // POSIX cannot be used to non-blocking read
		}

		if (chunksize) {
			if (loader->open(*MPITools::procs, name)) {
				eslog::error("LOADER: cannot read file '%s'\n", name.c_str());
			}
			if (chunks > 1) {
				eslog::error("LOADER: implement non-blocking multi-chunks reader\n");
			}
			for (size_t c = 0; c < chunks; ++c) {
				size_t size = std::min(chunkmax, data.size() - overlap - c * chunkmax);
				loader->iread(data.data() + chunkoffset - distribution[info::mpi::rank], chunkoffset, size);
				chunkoffset += size;
			}
		}
	}

	while (next()) {
		begin = data.data();
		hardend = begin + data.size();
		end = hardend - overlap;
	}
}

void AsyncFilePack::wait()
{
	while (next()) {
		loader->wait(); // waitall??
	}

	// overlap data
	std::vector<char, initless_allocator<char> > sBuffer(files.size() * overlap), rBuffer(files.size() * overlap);
	while (next()) {
		memcpy(sBuffer.data() + fileindex * overlap, data.data(), overlap);
	}
	Communication::receiveUpper(sBuffer, rBuffer);
	while (next()) {
		if (info::mpi::rank + 1 != info::mpi::size) {
			memcpy(data.data() + data.size() - overlap, rBuffer.data() + fileindex * overlap, overlap);
		}
	}

	callback();
	callback = [] () {};
}

void Metadata::read()
{
	profiler::syncstart("metadata_read");
	distribution.resize(2);
	if (info::mpi::rank == 0) {
		POSIXLoader loader;
		if (loader.open(MPITools::singleton->across, name)) {
			eslog::error("MESIO error: cannot load metadata file '%s'\n", name.c_str());
		}
		distribution = { 0, loader.size() };
		data.resize(distribution.back());
		profiler::synccheckpoint("prepare");
		loader.read(data.data(), 0, distribution.back());
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
