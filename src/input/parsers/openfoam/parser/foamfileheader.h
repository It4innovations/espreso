
#ifndef SRC_INPUT_PARSERS_OPENFOAM_PARSER_FOAMFILEHEADER_H_
#define SRC_INPUT_PARSERS_OPENFOAM_PARSER_FOAMFILEHEADER_H_

#include "basis/io/inputfile.h"
#include <fstream>

namespace espreso {

#define MAX_CHAR_LENGTH 256L

struct FoamFileHeader {
	int version = -1, subversion = -1;

	enum class Format { unknown, ASCII, binary }
	format = Format::unknown;

	enum class Class {
		unknown,
		faceList, faceCompactList, labelList, vectorField, polyBoundaryMesh,
		pointScalarField, pointVectorField, volScalarField, volVectorField, surfaceScalarField }
	foamClass = Class::unknown;

	char arch[MAX_CHAR_LENGTH] = { 0 };
	char note[MAX_CHAR_LENGTH] = { 0 };
	char location[MAX_CHAR_LENGTH] = { 0 };
	char object[MAX_CHAR_LENGTH] = { 0 };
	int label, scalar; // data width

	const char* read(InputFile *file);
	void read(std::ifstream &is);

	int dimension()
	{
		switch (foamClass) {
		case FoamFileHeader::Class::labelList:          return 1;
		case FoamFileHeader::Class::vectorField:        return 3;
		case FoamFileHeader::Class::pointScalarField:   return 1;
		case FoamFileHeader::Class::pointVectorField:   return 3;
		case FoamFileHeader::Class::volScalarField:     return 1;
		case FoamFileHeader::Class::volVectorField:     return 3;
		case FoamFileHeader::Class::surfaceScalarField: return 1;
		default: return 0;
		}
	}

	FoamFileHeader(): label(4), scalar(8) {}

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
		label = std::max(label, other.label);
		scalar = std::max(scalar, other.scalar);
		return *this;
	}
};

}

#endif /* SRC_INPUT_PARSERS_OPENFOAM_PARSER_FOAMFILEHEADER_H_ */
