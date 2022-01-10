
#include "foamFile.h"
#include "tokenizer.h"

#include "basis/containers/tarray.hpp"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/parsers/distributedscanner.h"

#include <cstdlib>

using namespace espreso;

const char* FoamFileHeader::read(InputFile *file)
{
	if (info::mpi::rank == 0) {
		const char *c = oftoken::toFoamFile(file->begin);
		while (*c != '}') {
			oftoken::dictionary dict = oftoken::dict(c);
			auto keyword = [dict] (const char *keyword) {
				return memcmp(dict.keyword.begin, keyword, std::min((size_t)(dict.keyword.end - dict.keyword.begin), strlen(keyword))) == 0;
			};
			auto value = [dict] (const char *value) {
				return memcmp(dict.value.begin, value, std::min((size_t)(dict.value.end - dict.value.begin), strlen(value))) == 0;
			};
			if (keyword("version")) {
				char *end;
				version = strtol(dict.value.begin, &end, 10);
				subversion = strtol(++end, NULL, 10);
			}
			if (keyword("format")) {
				if (value("ascii")) { format = Format::ASCII; }
				if (value("binary")) { format = Format::BINARY; }
			}
			if (keyword("arch")) {
				memcpy(arch, dict.value.begin + 1, std::min(dict.value.end - dict.value.begin - 3, MAX_CHAR_LENGTH));
			}
			if (keyword("class")) {
				if (value("faceList")) { foamClass = Class::faceList; }
				if (value("labelList")) { foamClass = Class::labelList; }
				if (value("vectorField")) { foamClass = Class::vectorField; }
			}
			if (keyword("location")) {
				memcpy(location, dict.value.begin + 1, std::min(dict.value.end - dict.value.begin - 3, MAX_CHAR_LENGTH));
			}
			if (keyword("object")) {
				memcpy(object, dict.value.begin, std::min(dict.value.end - dict.value.begin - 1, MAX_CHAR_LENGTH));
			}
			c = dict.value.end;
			while (*c++ != '\n');
		}

		if (format == Format::UNKNOWN) {
			eslog::internalFailure("unknown format of Open FOAM file (%s/%s).\n", location, object);
		}
		if (foamClass == Class::unknown) {
			eslog::internalFailure("unknown class of Open FOAM file (%s/%s).\n", location, object);
		}
		return c;
	}
	return file->begin;
}

void FoamFile::scan(InputFile *file, FoamFileDistribution &distribution)
{
	DistributedScanner::align(*file, "\n");

	const char *c = header.read(file);
	if (info::mpi::rank == 0) {
		while (*++c != '(');
		begin = c - file->begin;
		while (*c++ != '\n');
		const char *cc = file->begin + begin - 1;
		while (*--cc != '\n'); // go back to the size of the list
		size = strtol(cc + 1, NULL, 10);
	}

	distribution.distribution.resize(info::env::OMP_NUM_THREADS);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(info::env::OMP_NUM_THREADS, file->end - c);

	#pragma omp parallel for reduction(max:end)
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		const char* cend = c + tdistribution[t + 1];
		if (t && *(c - 1) != '\n') { while (c < cend && *c++ != '\n'); } // start at the new line

		FoamFileDistribution::Distribution dist{0, 0, c};
		while (c < cend) {
			if (*c == ')') { end = c - file->begin; break; }
			while (*c++ != '\n');
			++dist.size;
		}
		std::swap(distribution.distribution[t], dist);
	}

	end += file->distribution[info::mpi::rank];
	for (int t = 1; t < info::env::OMP_NUM_THREADS; ++t) {
		distribution.distribution[t].offset += distribution.distribution[t - 1].offset + distribution.distribution[t - 1].size;
	}
}




