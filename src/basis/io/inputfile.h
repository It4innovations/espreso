
#ifndef SRC_BASIS_IO_INPUTFILE_H_
#define SRC_BASIS_IO_INPUTFILE_H_

#include "basis/containers/allocators.h"

#include <string>
#include <vector>

namespace espreso {

struct Loader;
struct MPISubset;
struct MPIGroup;

struct InputFile {
	friend class FilePack;
	static size_t size(const std::string &file);

	InputFile(const std::string &name);
	~InputFile();

	void swap(InputFile *other)
	{
		std::swap(begin, other->begin);
		std::swap(end, other->end);
		std::swap(hardend, other->hardend);
		data.swap(other->data);
		distribution.swap(other->distribution);
		std::swap(totalSize, other->totalSize);
		name.swap(other->name);
		std::swap(maxchunk, other->maxchunk);
		std::swap(loader, other->loader);
	}

	void setDistribution(const std::vector<size_t> &distribution);

	const char *begin, *end, *hardend;
	std::vector<char, initless_allocator<char> > data;
	std::vector<size_t> distribution;

	size_t totalSize;
	std::string name;
protected:
	InputFile();

	size_t maxchunk;
	Loader* loader;
};

struct Metadata: public InputFile {

	Metadata(const std::string &name): InputFile(name) {}

	void read(const std::string &filename);
};

struct FilePack: public InputFile {
	FilePack();
	FilePack(const std::vector<std::string> &filepaths);
	~FilePack();

	InputFile* add(const std::string &name);
	bool next();

	void setTotalSizes();

	size_t fileindex;
	std::vector<InputFile*> files;
};

struct InputFilePack: public FilePack {
	InputFilePack(size_t minchunk = 64 * 1024, size_t overlap = 1024);
	InputFilePack(const std::vector<std::string> &filepaths, size_t minchunk = 64 * 1024, size_t overlap = 1024);

	void prepare();
	void read();

	size_t minchunk, overlap;
};

struct AsyncFilePack: public FilePack {
	AsyncFilePack();
	AsyncFilePack(const std::vector<std::string> &filepaths);

	void iread();
	void wait();

	const size_t overlap = 0; // find better value
};

}

#endif /* SRC_BASIS_IO_INPUTFILE_H_ */
