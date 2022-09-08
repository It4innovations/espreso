
#ifndef SRC_BASIS_IO_OUTPUTFILE_H_
#define SRC_BASIS_IO_OUTPUTFILE_H_

#include "basis/containers/allocators.h"

#include <string>
#include <vector>

namespace espreso {

class Writer;

class OutputFile {
friend class OutputFilePack;

public:
	void insert(int n)
	{
		_buffer.insert(_buffer.end(), buffer, buffer + n);
	}

	void insert(int n, char c)
	{
		_buffer.insert(_buffer.end(), n, c);
	}

	void insert(int size, const void *data)
	{
		_buffer.insert(_buffer.end(), reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data) + size);
	}

	// temporary buffer that is used during conversion to string by snprintf
	static const size_t bsize = 4 * 1024;
	static char buffer[bsize];

protected:
	OutputFile(): _offset(0) {}

	void _group();

	std::string _name;
	std::vector<char, initless_allocator<char> > _buffer;

	int _offset;
	std::vector<size_t> _distribution;
};

class OutputFilePack: public OutputFile {
public:
	~OutputFilePack();

	char* enlarge(size_t size)
	{
		size_t bsize = _buffer.size();
		_buffer.resize(bsize + size, 0);
		return _buffer.data() + bsize;
	}

	void groupData();
	void commitFile(const std::string &name);

	void reorder();
	void write();

	void directWrite();

protected:
	void clear();
	std::vector<OutputFile*> _files;
};

}

#endif /* SRC_BASIS_IO_OUTPUTFILE_H_ */
