
#include "outfile.h"

#include "esinfo/eslog.hpp"

using namespace espreso;


void OutFile::prepare()
{
	size_t sum = 0;
	for (size_t i = 0; i < blocks.size(); ++i) {
		blocks[i].offset = sum;
		sum += blocks[i].size;
	}
	data.resize(sum);
}

void OutFile::open(const std::string &name)
{
	if (_writer.open(MPITools::singleton->within, name, blocks.size())) {
		eslog::error("WRITER: cannot create file '%s'\n", name.c_str());
	}
}

void OutFile::store(size_t block)
{
	_writer.store(data.data() + blocks[block].offset, blocks[block].fileoffset, blocks[block].size);
}

void OutFile::close()
{
	_writer.close();
}
