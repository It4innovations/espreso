
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_

#include "foamfileheader.h"
#include "basis/io/inputfile.h"

#include "mpi.h"

namespace espreso {

struct RawParallelFoamFile {
	struct Distribution { size_t offset, size; const char *c; };

	RawParallelFoamFile(InputFile *input): input(input) {}

	std::vector<Distribution> distribution;
	InputFile *input;
};

struct ParallelFoamFile {
	FoamFileHeader header;
	size_t size = 0, begin = 0, end = 0;

	static void init();
	static void synchronize(const std::vector<ParallelFoamFile*> &files);
	static void finish();

	ParallelFoamFile& operator^=(const ParallelFoamFile &other) // used in MPI Reduce function
	{
		header ^= other.header;
		size = std::max(size, other.size);
		begin = std::max(begin, other.begin);
		end = std::max(end, other.end);
		return *this;
	}

protected:
	void scanFile(RawParallelFoamFile &file);

private:
	static MPI_Datatype oftype;
	static MPI_Op ofop;
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_ */
