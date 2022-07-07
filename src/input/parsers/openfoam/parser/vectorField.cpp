
#include "vectorField.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <fstream>

using namespace espreso;

template <int S> static esfloat _parse(const char* &c);

template <> esfloat _parse<4>(const char* &c)
{
	char *next;
	float value = strtof(c, &next);
	c = next;
	return value;
}

template <> esfloat _parse<8>(const char* &c)
{
	char *next;
	double value = strtof(c, &next);
	c = next;
	return value;
}

template <int T> static esfloat _parse(std::ifstream &is);

template <> esfloat _parse<4>(std::ifstream &is)
{
	float number;
	is.read(reinterpret_cast<char*>(&number), 4);
	return number;
}

template <> esfloat _parse<8>(std::ifstream &is)
{
	double number;
	is.read(reinterpret_cast<char*>(&number), 8);
	return number;
}

void OpenFOAMVectorField::parse(ivector<_Point<esfloat> > &coordinates)
{
	coordinates.resize(distribution.back().offset + distribution.back().size);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		const char *c = distribution[t].c;
		for (size_t i = distribution[t].offset; i < distribution[t].offset + distribution[t].size; ++i) {
			c += 1; // skip "("
			switch (header.label) {
			case 4:
				coordinates[i].x = _parse<4>(c);
				coordinates[i].y = _parse<4>(c);
				coordinates[i].z = _parse<4>(c);
				break;
			case 8:
				coordinates[i].x = _parse<8>(c);
				coordinates[i].y = _parse<8>(c);
				coordinates[i].z = _parse<8>(c);
				break;
			}
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
		for (esint p = 0; p < points; ++p) {
			coordinates.push_back({});
			switch (header.scalar) {
			case 4:
				coordinates.back().x = _parse<4>(is);
				coordinates.back().y = _parse<4>(is);
				coordinates.back().z = _parse<4>(is);
				break;
			case 8:
				coordinates.back().x = _parse<8>(is);
				coordinates.back().y = _parse<8>(is);
				coordinates.back().z = _parse<8>(is);
				break;
			}
		}
	}
	return points;
}
