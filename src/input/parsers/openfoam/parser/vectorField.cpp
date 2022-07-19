
#include "vectorField.h"

#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

#include <fstream>
#include <cstring>

using namespace espreso;

template <typename T, int S> static T _parseASCII(const char* &c)
{
	char *next;
	T value = strtof(c, &next);
	c = next;
	return value;
}

template <typename T, int S> static T _parseBinary(const char* &c)
{
	T number;
	memcpy(reinterpret_cast<char*>(&number), c, S);
	c += S;
	return number;
}

template <typename T, int S> static T _parse(std::ifstream &is)
{
	T number;
	is.read(reinterpret_cast<char*>(&number), S);
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
			switch (header.scalar) {
			case 4:
				coordinates[i].x = _parseASCII<float, 4>(c);
				coordinates[i].y = _parseASCII<float, 4>(c);
				coordinates[i].z = _parseASCII<float, 4>(c);
				break;
			case 8:
				coordinates[i].x = _parseASCII<double, 8>(c);
				coordinates[i].y = _parseASCII<double, 8>(c);
				coordinates[i].z = _parseASCII<double, 8>(c);
				break;
			}
			c += 2; // skip ")\n"
		}
	}
}

void OpenFOAMVectorField::load(InputFile *input, std::vector<esfloat> &data)
{
	FoamFileHeader header;
	input->begin = header.read(input);
	while (*input->begin != '(') {
		++input->begin;
	}
	// go back to the size of the list
	--input->begin;
	--input->begin;
	while (*input->begin != '\n') {
		--input->begin;
	}
	++input->begin;

	int dimension = header.dimension();
	// TODO: use threads
	const char *c = input->begin;
	esint size = strtol(c, nullptr, 10);
	while (*c++ != '(');
	if (header.format == FoamFileHeader::Format::ASCII) {
		while (*c++ != '\n');
		for (esint i = 0; i < size; ++i) {
			if (dimension > 1) {
				c += 1; // skip "("
			}
			switch (header.scalar) {
			case 4:
				for (int d = 0; d < dimension; ++d) {
					data.push_back(_parseASCII<float, 4>(c));
				}
				break;
			case 8:
				for (int d = 0; d < dimension; ++d) {
					data.push_back(_parseASCII<double, 8>(c));
				}
				break;
			}
			if (dimension > 1) {
				c += 2; // skip ")\n"
			}
		}
	} else {
		for (esint i = 0; i < size; ++i) {
			switch (header.scalar) {
			case 4:
				for (int d = 0; d < dimension; ++d) {
					data.push_back(_parseBinary<float, 4>(c));
				}
				break;
			case 8:
				for (int d = 0; d < dimension; ++d) {
					data.push_back(_parseBinary<double, 8>(c));
				}
				break;
			}
		}
	}
}

esint OpenFOAMVectorField::load(const std::string &file, ivector<_Point<esfloat> > &coordinates)
{
	std::ifstream is(file);
	FoamFileHeader header; header.read(is);

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
				coordinates.back().x = _parse<float, 4>(is);
				coordinates.back().y = _parse<float, 4>(is);
				coordinates.back().z = _parse<float, 4>(is);
				break;
			case 8:
				coordinates.back().x = _parse<double, 8>(is);
				coordinates.back().y = _parse<double, 8>(is);
				coordinates.back().z = _parse<double, 8>(is);
				break;
			}
		}
	}
	return points;
}
