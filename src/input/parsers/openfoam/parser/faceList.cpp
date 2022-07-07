
#include "faceList.h"

#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"

using namespace espreso;

template <int T> static esint _parse(std::ifstream &is);

template <> esint _parse<4>(std::ifstream &is)
{
	int number;
	is.read(reinterpret_cast<char*>(&number), 4);
	return number;
}

template <> esint _parse<8>(std::ifstream &is)
{
	long number;
	is.read(reinterpret_cast<char*>(&number), 8);
	return number;
}

void OpenFOAMFaceList::parse(ivector<Element::CODE> &type, ivector<esint> &enodes)
{
	type.resize(distribution.back().offset + distribution.back().size);
	std::vector<std::vector<esint> > nodes(info::env::OMP_NUM_THREADS);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::vector<esint> tnodes;
		tnodes.reserve(4 * distribution[t].size); // more space for general polygons?

		char* next;
		const char *c = distribution[t].c;
		for (size_t i = distribution[t].offset; i < distribution[t].offset + distribution[t].size; ++i) {
			int size = strtol(c, &next, 10); c = next + 1; // skip "("
			switch (size) {
			case 3: type[i] = Element::CODE::TRIANGLE3; break;
			case 4: type[i] = Element::CODE::SQUARE4; break;
			default:
				type[i] = Element::CODE::POLYGON; tnodes.push_back(size);
			}
			for (esint n = 0; n < size; ++n) {
				tnodes.push_back(strtol(c, &next, 10)); c = next;
			}
			c += 2; // skip ")\n"
		}

		nodes[t].swap(tnodes);
	}

	std::vector<esint> ncount(info::env::OMP_NUM_THREADS + 1);
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		ncount[t] = nodes[t].size();
	}
	utils::sizesToOffsets(ncount);
	enodes.resize(ncount.back());

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		std::copy(nodes[t].begin(), nodes[t].end(), enodes.begin() + ncount[t]);
	}
}

FoamFileHeader OpenFOAMFaceList::load(const std::string &file, ivector<Element::CODE> &type, ivector<esint> &enodes, esint offset)
{
	std::ifstream is(file);
	FoamFileHeader header(is);

	if (header.foamClass == FoamFileHeader::Class::faceList) {
		esint faces;
		is >> faces;
		type.reserve(type.size() + faces);
		enodes.reserve(enodes.size() + faces * 4.1);
		is.ignore(256, '(');
		if (header.format == FoamFileHeader::Format::ASCII) {
			for (esint f = 0; f < faces; ++f) {
				esint size, node;
				is >> size;
				is.ignore(16, '(');
				switch (size) {
				case 3: type.push_back(Element::CODE::TRIANGLE3); break;
				case 4: type.push_back(Element::CODE::SQUARE4); break;
				default: type.push_back(Element::CODE::POLYGON); enodes.push_back(size);
				}
				for (esint n = 0; n < size; ++n) {
					is >> node;
					enodes.push_back(node + offset);
				}
				is.ignore(16, ')');
			}
		} else {
			for (esint f = 0; f < faces; ++f) {
				esint size, node;
				switch (header.label) {
				case 4: size = _parse<4>(is); break;
				case 8: size = _parse<8>(is); break;
				}
				switch (size) {
				case 3: type.push_back(Element::CODE::TRIANGLE3); break;
				case 4: type.push_back(Element::CODE::SQUARE4); break;
				default: type.push_back(Element::CODE::POLYGON); enodes.push_back(size);
				}
				for (esint n = 0; n < size; ++n) {
					switch (header.label) {
					case 4: node = _parse<4>(is); break;
					case 8: node = _parse<8>(is); break;
					}
					enodes.push_back(node + offset);
				}
			}
		}
	}

	if (header.foamClass == FoamFileHeader::Class::faceCompactList) {
		esint faces;
		is >> faces;
		type.reserve(type.size() + faces - 1);
		enodes.reserve(type.size() + faces * 4.1);
		ivector<esint> dist(faces);
		is.ignore(256, '(');
		if (header.format == FoamFileHeader::Format::ASCII) {
			for (esint f = 0; f < faces; ++f) {
				is >> dist[f];
			}
			for (esint f = 1; f < faces; ++f) {
				switch (dist[f] - dist[f - 1]) {
				case 3: type.push_back(Element::CODE::TRIANGLE3); break;
				case 4: type.push_back(Element::CODE::SQUARE4); break;
				default: type.push_back(Element::CODE::POLYGON); enodes.push_back(dist[f] - dist[f - 1]);
				}
				for (esint n = dist[f - 1]; n < dist[f]; ++n) {
					esint node;
					is >> node;
					enodes.push_back(node + offset);
				}
			}
		} else {
			for (esint f = 0; f < faces; ++f) {
				switch (header.label) {
				case 4: dist[f] = _parse<4>(is); break;
				case 8: dist[f] = _parse<8>(is); break;
				}
			}
			is.ignore(256, '(');
			for (esint f = 1; f < faces; ++f) {
				switch (dist[f] - dist[f - 1]) {
				case 3: type.push_back(Element::CODE::TRIANGLE3); break;
				case 4: type.push_back(Element::CODE::SQUARE4); break;
				default: type.push_back(Element::CODE::POLYGON); enodes.push_back(dist[f] - dist[f - 1]);
				}
				for (esint n = dist[f - 1]; n < dist[f]; ++n) {
					switch (header.label) {
					case 4: enodes.push_back(_parse<4>(is) + offset); break;
					case 8: enodes.push_back(_parse<8>(is) + offset); break;
					}
				}
			}
		}
	}
	return header;
}
