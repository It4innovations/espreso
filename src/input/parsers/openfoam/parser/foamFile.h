
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_

#include "basis/io/inputfile.h"

#include "mpi.h"

namespace espreso {

#define MAX_CHAR_LENGTH 256L

struct FoamFileHeader {
	int version = -1, subversion = -1;

	enum class Format { UNKNOWN, ASCII, BINARY }
	format = Format::UNKNOWN;

	enum class Class { unknown, faceList, labelList, vectorField }
	foamClass = Class::unknown;

	char arch[MAX_CHAR_LENGTH] = { 0 };
	char note[MAX_CHAR_LENGTH] = { 0 };
	char location[MAX_CHAR_LENGTH] = { 0 };
	char object[MAX_CHAR_LENGTH] = { 0 };

	const char* read(InputFile *file);

	FoamFileHeader& operator^=(const FoamFileHeader &other) // used in MPI Reduce function
	{
		version = std::max(version, other.version);
		subversion = std::max(subversion, other.subversion);
		format = std::max(format, other.format);
		foamClass = std::max(foamClass, other.foamClass);
		for (long i = 0; i < MAX_CHAR_LENGTH; ++i) {
			arch[i] ^= other.arch[i];
			note[i] ^= other.note[i];
			location[i] ^= other.location[i];
			object[i] ^= other.object[i];
		}
		return *this;
	}
};

struct RawFoamFile {
	struct Distribution { size_t offset, size; const char *c; };

	RawFoamFile(InputFile *input): input(input) {}

	std::vector<Distribution> distribution;
	InputFile *input;
};

struct FoamFile {
	FoamFileHeader header;
	size_t size = 0, begin = 0, end = 0;

	static void init();
	static void synchronize(const std::vector<FoamFile*> &files);
	static void finish();

	FoamFile& operator^=(const FoamFile &other) // used in MPI Reduce function
	{
		header ^= other.header;
		size = std::max(size, other.size);
		begin = std::max(begin, other.begin);
		end = std::max(end, other.end);
		return *this;
	}

protected:
	void scanFile(RawFoamFile &file);

private:
	static MPI_Datatype oftype;
	static MPI_Op ofop;
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_CONTAINERS_FOAMFILE_H_ */
