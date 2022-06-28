
#ifndef SRC_INPUT_PARSERS_FILEBLOCK_H_
#define SRC_INPUT_PARSERS_FILEBLOCK_H_

#include "basis/io/inputfile.h"

#include <cstddef>

namespace espreso {

struct FileBlock {
	size_t bytesize;   // size of block in bytes
	size_t prevsize;   // size of block parsed by lower processes
	size_t size;       // end - begin / step
	size_t nextsize;   // size of block parsed by upper processes
	size_t begin, end; // points where parsing should be started / finished
	size_t step;

	void round(size_t &current, const size_t &start, const size_t &chunk)
	{
		if ((current - start) % chunk != 0) {
			current += chunk - (current - start) % chunk;
		}
	}

	FileBlock(const FilePack &file, size_t start, size_t size, size_t step, int rank)
	: prevsize(0),
	  size(0),
	  nextsize(0),
	  begin(std::max(start, file.distribution[rank])),
	  end(std::min(start + size, file.distribution[rank + 1])),
	  step(step)
	{
		round(begin, start, step);
		round(end, start, step);
		if (begin < end) {
			this->prevsize = begin - start;
			this->size = end - begin;
			this->nextsize = start + size - end;
			begin -= file.distribution[rank];
			end -= file.distribution[rank];
		} else {
			if (start + size < file.distribution[rank]) {
				this->prevsize = size;
				this->size = 0;
				this->nextsize = 0;
			} else {
				this->prevsize = 0;
				this->size = 0;
				this->nextsize = size;
			}
		}

		this->bytesize = this->size;
		this->prevsize /= step;
		this->size /= step;
		this->nextsize /= step;
	}
};

}

#endif /* SRC_INPUT_PARSERS_FILEBLOCK_H_ */
