
#include "vectorField.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <fstream>

using namespace espreso;

template <typename T> static void _parse(const char* &c, T &value);

template <> void _parse(const char* &c, float &value)
{
	char *next;
	value = strtof(c, &next);
	c = next;
}

template <> void _parse(const char* &c, double &value)
{
	char *next;
	value = strtod(c, &next);
	c = next;
}

void OpenFOAMVectorField::parse(ivector<_Point<esfloat> > &coordinates)
{
	coordinates.resize(distribution.back().offset + distribution.back().size);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {

		const char *c = distribution[t].c;
		for (size_t i = distribution[t].offset; i < distribution[t].offset + distribution[t].size; ++i) {
			c += 1; // skip "("
			_parse(c, coordinates[i].x);
			_parse(c, coordinates[i].y);
			_parse(c, coordinates[i].z);
			c += 2; // skip ")\n"
		}
	}
}

esint OpenFOAMVectorField::load(const std::string &file, ivector<_Point<esfloat> > &coordinates)
{
	std::ifstream is(file);
	FoamFileHeader header(is);

	esint points;
	is >> points;
	coordinates.reserve(coordinates.size() + points);
	is.ignore(256, '(');
	if (header.format == FoamFileHeader::Format::ASCII) {
		for (esint p = 0; p < points; ++p) {
			is.ignore(16, '(');
			_Point<esfloat> coo;
			is >> coo.x >> coo.y >> coo.z;
			coordinates.push_back(coo);
		}
	} else {
		eslog::error("implement OpenFOAM binary reader.\n");
	}
	return points;
}
