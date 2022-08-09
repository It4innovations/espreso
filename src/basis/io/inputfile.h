
#ifndef SRC_BASIS_IO_INPUTFILE_H_
#define SRC_BASIS_IO_INPUTFILE_H_

#include "basis/containers/allocators.h"

#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct Loader;
struct MPISubset;
struct MPIGroup;

struct InputFile {
	friend class FilePack;
	static size_t size(const std::string &file);

	InputFile(const std::string &name, size_t overlap);
	~InputFile();

	void swap(InputFile *other)
	{
		std::swap(begin, other->begin);
		std::swap(end, other->end);
		std::swap(hardend, other->hardend);
		data.swap(other->data);
		distribution.swap(other->distribution);
		std::swap(overlap, other->overlap);
		std::swap(totalSize, other->totalSize);
		name.swap(other->name);
		std::swap(maxchunk, other->maxchunk);
		std::swap(loader, other->loader);

	}

	const char *begin, *end, *hardend;
	std::vector<char, initless_allocator<char> > data;
	std::vector<size_t> distribution;

	size_t overlap, totalSize;
	std::string name;
protected:
	InputFile(size_t overlap);

	size_t maxchunk;
	Loader* loader;
};

struct Metadata: public InputFile {

	Metadata(const std::string &name): InputFile(name, 0) { read(); }

	void read();
};

struct FilePack: public InputFile {
	FilePack(size_t minchunk, size_t overlap);
	FilePack(const std::vector<std::string> &filepaths, size_t minchunk, size_t overlap);
	~FilePack();

	InputFile* add(const std::string &name);
	bool next();

	void setTotalSizes(const MPIGroup &group);

	const size_t minchunk;
	size_t fileindex;
	std::vector<InputFile*> files;
};

struct InputFilePack: public FilePack {
	InputFilePack(size_t minchunk = 64 * 1024, size_t overlap = 1024);
	InputFilePack(const std::vector<std::string> &filepaths, size_t minchunk = 64 * 1024, size_t overlap = 1024);

	void prepare();
	void read();
};

struct AsyncFilePack: public FilePack {
	AsyncFilePack(size_t overlap = 1024);
	AsyncFilePack(const std::vector<std::string> &filepaths, size_t overlap = 1024);

	void iread(const MPIGroup &group);
	void wait(const MPIGroup &group);
};

}

#endif /* SRC_BASIS_IO_INPUTFILE_H_ */
