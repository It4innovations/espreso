
#include "foamFile.h"
#include "tokenizer.h"

#include "basis/containers/tarray.hpp"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "input/parsers/distributedscanner.h"

using namespace espreso;

MPI_Datatype FoamFile::oftype;
MPI_Op FoamFile::ofop;

static void _ofreduce(void *in, void *out, int *len, MPI_Datatype *datatype)
{
	for (int i = 0; i < *len; i++) {
		*(reinterpret_cast<FoamFile*>(out) + i) ^= *(reinterpret_cast<FoamFile*>(in) + i);
	}
}

void FoamFile::init()
{
	MPI_Type_contiguous(sizeof(FoamFile), MPI_BYTE, &oftype);
	MPI_Type_commit(&oftype);

	MPI_Op_create(_ofreduce, 1, &ofop);
}

void FoamFile::synchronize(const std::vector<FoamFile*> &files)
{
	std::vector<FoamFile> ff(files.size());
	for (size_t i = 0; i < files.size(); ++i) {
		ff[i] = *files[i];
	}
	Communication::allReduce(ff.data(), NULL, ff.size(), FoamFile::oftype, FoamFile::ofop);

	for (size_t i = 0; i < files.size(); ++i) {
		*files[i] = ff[i];
	}
}

void FoamFile::finish()
{
	MPI_Op_free(&ofop);
	MPI_Type_free(&oftype);
}

void parse(FoamFileHeader &header, const oftoken::dictionary &dict)
{
	auto keyword = [dict] (const char *keyword) {
		return memcmp(dict.keyword.begin, keyword, std::min((size_t)(dict.keyword.end - dict.keyword.begin), strlen(keyword))) == 0;
	};
	auto value = [dict] (const char *value) {
		return memcmp(dict.value.begin, value, std::min((size_t)(dict.value.end - dict.value.begin), strlen(value))) == 0;
	};

	if (keyword("version")) {
		char *end;
		header.version = strtol(dict.value.begin, &end, 10);
		header.subversion = strtol(++end, NULL, 10);
	}
	if (keyword("format")) {
		if (value("ascii")) { header.format = FoamFileHeader::Format::ASCII; }
		if (value("binary")) { header.format = FoamFileHeader::Format::binary; }
	}
	if (keyword("arch")) {
		memcpy(header.arch, dict.value.begin + 1, std::min(dict.value.end - dict.value.begin - 2, MAX_CHAR_LENGTH));
		if (memcmp(header.arch, "LSB;", 4) != 0) {
			eslog::internalFailure("implement non LSB reader.\n");
		}
		if (memcmp(header.arch + 4, "label=32;", 9) == 0) {
			header.label = 4;
		}
		if (memcmp(header.arch + 4, "label=64;", 9) == 0) {
			header.label = 8;
		}
		if (memcmp(header.arch + 13, "scalar=32", 9) == 0) {
			header.scalar = 4;
		}
		if (memcmp(header.arch + 13, "scalar=64", 9) == 0) {
			header.scalar = 8;
		}
	}
	if (keyword("note")) {
		memcpy(header.note, dict.value.begin + 1, std::min(dict.value.end - dict.value.begin - 2, MAX_CHAR_LENGTH));
	}
	if (keyword("class")) {
		if (value("faceList")) { header.foamClass = FoamFileHeader::Class::faceList; }
		if (value("faceCompactList")) { header.foamClass = FoamFileHeader::Class::faceCompactList; }
		if (value("labelList")) { header.foamClass = FoamFileHeader::Class::labelList; }
		if (value("vectorField")) { header.foamClass = FoamFileHeader::Class::vectorField; }
	}
	if (keyword("location")) {
		memcpy(header.location, dict.value.begin + 1, std::min(dict.value.end - dict.value.begin - 3, MAX_CHAR_LENGTH));
	}
	if (keyword("object")) {
		memcpy(header.object, dict.value.begin, std::min(dict.value.end - dict.value.begin - 1, MAX_CHAR_LENGTH));
	}
}

const char* FoamFileHeader::read(InputFile *file)
{
	if (file->distribution[info::mpi::rank] == 0 && file->distribution[info::mpi::rank + 1] != 0) {
		const char *c = oftoken::toFoamFile(file->begin);
		while (*c != '}') {
			oftoken::dictionary dict = oftoken::dict(c);
			parse(*this, dict);
			c = dict.value.end;
			while (*c++ != '\n');
		}

		if (format == Format::unknown) {
			eslog::internalFailure("unknown format of Open FOAM file (%s/%s).\n", location, object);
		}
		if (foamClass == Class::unknown) {
			eslog::internalFailure("unknown class of Open FOAM file (%s/%s).\n", location, object);
		}
		return c;
	}
	return file->begin;
}

void FoamFileHeader::read(std::ifstream &is)
{
	oftoken::toFoamFile(is);
	char line[256];
	is.getline(line, 256);
	while(line[0] != '}') {
		oftoken::dictionary dict = oftoken::dict(line);
		parse(*this, dict);
		is.getline(line, 256);
	}
	is.getline(line, 256);

	if (format == Format::unknown) {
		eslog::internalFailure("unknown format of Open FOAM file (%s/%s).\n", location, object);
	}
	if (foamClass == Class::unknown) {
		eslog::internalFailure("unknown class of Open FOAM file (%s/%s).\n", location, object);
	}
}

void FoamFile::scanFile(RawFoamFile &file)
{
	DistributedScanner::align(*file.input, "\n");

	const char *c = header.read(file.input);
	if (file.input->distribution[info::mpi::rank] == 0 && file.input->distribution[info::mpi::rank + 1] != 0) {
		while (*++c != '(');
		begin = c - file.input->begin;
		while (*c++ != '\n');
		const char *cc = file.input->begin + begin - 1;
		while (*--cc != '\n'); // go back to the size of the list
		size = strtol(cc + 1, NULL, 10);
	}

	size_t coffset = c - file.input->begin;
	file.distribution.resize(info::env::OMP_NUM_THREADS);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(info::env::OMP_NUM_THREADS, file.input->end - c);

	#pragma omp parallel for reduction(max:end)
	for (int t = 0; t < info::env::OMP_NUM_THREADS; ++t) {
		const char* tc = file.input->begin + coffset + tdistribution[t];
		const char* cend = file.input->begin + coffset + tdistribution[t + 1];
		if (t && *(tc - 1) != '\n') { while (tc < cend && *tc++ != '\n'); } // start at the new line

		RawFoamFile::Distribution dist{0, 0, tc};
		while (tc < cend) {
			if (*tc == ')') { end = tc - file.input->begin; break; }
			while (tc < cend && *tc++ != '\n');
			++dist.size;
		}
		std::swap(file.distribution[t], dist);
	}

	end += file.input->distribution[info::mpi::rank];
	for (int t = 1; t < info::env::OMP_NUM_THREADS; ++t) {
		file.distribution[t].offset += file.distribution[t - 1].offset + file.distribution[t - 1].size;
	}
}




