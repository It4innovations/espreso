
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
	friend class InputFilePack;

	InputFile();

	const char *begin, *end, *hardend;
	std::vector<char, initless_allocator<char> > data;
	std::vector<size_t> distribution;

protected:
	size_t maxchunk;
};

struct Metadata: public InputFile {

	void read(const std::string &filename);
};

struct InputFilePack: public InputFile {
	InputFilePack(size_t minchunk = 64 * 1024, size_t overlap = 1024);
	~InputFilePack();

	size_t size() { return files.size(); }

	void commitFiles(const std::vector<std::string> &filepaths);
	bool next();

	void prepare();
	void read();
	void clear();

	size_t fileindex;
	size_t minchunk, overlap;
	std::vector<std::string> paths;
	std::vector<InputFile*> files;
};

}

#endif /* SRC_BASIS_IO_INPUTFILE_H_ */
