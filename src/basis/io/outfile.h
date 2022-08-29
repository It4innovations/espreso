
#ifndef SRC_BASIS_IO_OUTFILE_H_
#define SRC_BASIS_IO_OUTFILE_H_

#include "basis/containers/allocators.h"
#include "writer.h"

#include <string>
#include <vector>

namespace espreso {

class MPIWriter;

class OutFile {
public:
	struct Block { size_t fileoffset = 0, size = 0, offset = 0;  };
	std::vector<Block> blocks;
	std::vector<char, initless_allocator<char> > data;

	void prepare();
	void open(const std::string &name);
	void store(size_t block);
	void close();

protected:
	MPIAsyncWriter _writer;
};

}



#endif /* SRC_BASIS_IO_OUTFILE_H_ */
