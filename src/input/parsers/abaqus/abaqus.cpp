
#include "abaqus.h"

#include "parser/nblock.h"
#include "parser/eblock.h"
#include "parser/blockend.h"
#include "parser/nsel.h"
#include "parser/eset.h"
#include "parser/ssection.h"
#include "parser/material.h"
#include "parser/elemat.h"
#include "parser/nset.h"

#include "basis/containers/tarray.h"
#include "wrappers/mpi/communication.h"
#include "config/ecf/input/input.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include "parser/parser.h"
// msg remove this

using namespace espreso;

AbaqusLoader::AbaqusLoader(const InputConfiguration &configuration)
: _configuration(configuration), _pfile(32 * MAX_LINE_STEP * MAX_LINE_SIZE, MAX_LINE_STEP * MAX_LINE_SIZE)
{

}

void AbaqusLoader::load()
{
	eslog::start("ABAQUS: STARTED", "ABAQUS");
	eslog::param("database", _configuration.path.c_str());
	eslog::ln();

	_pfile.commitFiles({ _configuration.path });
	_pfile.prepare();
	eslog::checkpointln("ABAQUS: READING PREPARED");

	_pfile.read();
	eslog::checkpointln("ABAQUS: DATA READ");

	_pfile.next();
	AbaqusParser::offset = _pfile.distribution[info::mpi::rank];
	AbaqusParser::begin = _pfile.begin;
	AbaqusParser::end = _pfile.end;

	scan();
	eslog::checkpointln("ABAQUS: DATA SCANNED");

	parse();
	if (!_configuration.keep_material_sets) {
		std::fill(material.begin(), material.end(), 0);
	}
	eslog::endln("ABAQUS: DATA PARSED");
}

void AbaqusLoader::scan()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, _pfile.end - _pfile.begin);
	std::vector<std::vector<NList> > tNLists(threads);
	std::vector<std::vector<EList> > tELists(threads);
	std::vector<std::vector<BlockFinish> > tBlockFinishs(threads);
	std::vector<std::vector<Eset> > tEsets(threads);
	std::vector<std::vector<Nset> > tNsets(threads);
	std::vector<std::vector<SSection> > tSSections(threads);
	std::vector<std::vector<AbaqusMaterial> > tMaterials(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const char* tbegin = _pfile.begin + tdistribution[t];
		if (t && *(tbegin - 1) != '\n') {
			while (tbegin < _pfile.end && *tbegin++ != '\n'); // start at new line
		}

		while (tbegin < _pfile.begin + tdistribution[t + 1]) {
			

			if (memcmp(tbegin, BlockFinish::asterik, BlockFinish::nSize) == 0  &&
				memcmp(tbegin , BlockFinish::comment, BlockFinish::ncSize) != 0) {

				tBlockFinishs[t].push_back(BlockFinish().parse(tbegin ));
			}

			if (memcmp(tbegin, NList::upper, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
			}
			if (memcmp(tbegin, NList::lower, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
			}
			if (memcmp(tbegin, NList::sentence, NList::size) == 0) {
				tNLists[t].push_back(NList().parse(tbegin));
			}

			if (memcmp(tbegin, Eset::upper, Eset::size) == 0) {
				tEsets[t].push_back(Eset().parse(tbegin));
			}
			if (memcmp(tbegin, Eset::lower, Eset::size) == 0) {
				tEsets[t].push_back(Eset().parse(tbegin));
			}
			if (memcmp(tbegin, Eset::sentence, Eset::size) == 0) {
				tEsets[t].push_back(Eset().parse(tbegin));
			}

			if (memcmp(tbegin, EList::upper, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}
			if (memcmp(tbegin, EList::lower, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}
			if (memcmp(tbegin, EList::sentence, EList::size) == 0) {
				tELists[t].push_back(EList().parse(tbegin));
			}

			if (memcmp(tbegin, SSection::upper, SSection::size) == 0) {
				tSSections[t].push_back(SSection().parse(tbegin));
			}
			if (memcmp(tbegin, SSection::lower, SSection::size) == 0) {
				tSSections[t].push_back(SSection().parse(tbegin));
			}
			if (memcmp(tbegin, SSection::sentence, SSection::size) == 0) {
				tSSections[t].push_back(SSection().parse(tbegin));
			}

			if (memcmp(tbegin, AbaqusMaterial::upper, AbaqusMaterial::size) == 0) {
				tMaterials[t].push_back(AbaqusMaterial().parse(tbegin));
			}
			if (memcmp(tbegin, AbaqusMaterial::lower, AbaqusMaterial::size) == 0) {
				tMaterials[t].push_back(AbaqusMaterial().parse(tbegin));
			}
			if (memcmp(tbegin, AbaqusMaterial::sentence, AbaqusMaterial::size) == 0) {
				tMaterials[t].push_back(AbaqusMaterial().parse(tbegin));
			}

			if (memcmp(tbegin, Nset::upper, Nset::size) == 0) {
				tNsets[t].push_back(Nset().parse(tbegin));
			}
			if (memcmp(tbegin, Nset::lower, Nset::size) == 0) {
				tNsets[t].push_back(Nset().parse(tbegin));
			}
			if (memcmp(tbegin, Nset::sentence, Nset::size) == 0) {
				tNsets[t].push_back(Nset().parse(tbegin));
			}

			while (tbegin < _pfile.end && *tbegin++ != '\n');
		}
	}
	for (size_t t = 0; t < threads; t++) {
		_NLists.insert(_NLists.end(), tNLists[t].begin(), tNLists[t].end());
		_blockFinishs.insert(_blockFinishs.end(), tBlockFinishs[t].begin(), tBlockFinishs[t].end());
		_ELists.insert(_ELists.end(), tELists[t].begin(), tELists[t].end());
		_Esets.insert(_Esets.end(), tEsets[t].begin(), tEsets[t].end());
		_SSections.insert(_SSections.end(), tSSections[t].begin(), tSSections[t].end());
		_Materials.insert(_Materials.end(), tMaterials[t].begin() ,tMaterials[t].end());
		_Nsets.insert(_Nsets.end(), tNsets[t].begin() ,tNsets[t].end());
	}

	if (!Communication::allGatherUnknownSize(_blockFinishs)) {
		eslog::internalFailure("exchange ABAQUS block ends.\n");
	}

	for (size_t i = 0; i < _NLists.size(); i++) {
		_NLists[i].fillDistribution(_blockFinishs, _pfile.distribution);
	}

	for (size_t i = 0; i < _ELists.size(); i++) {
		_ELists[i].fillDistribution(_blockFinishs, _pfile.distribution);
	}

	for (size_t i = 0; i < _Esets.size(); i++) {
		_Esets[i].fillDistribution(_blockFinishs, _pfile.distribution);
	}

	for (size_t i = 0; i < _Nsets.size(); i++) {
		_Nsets[i].fillDistribution(_blockFinishs, _pfile.distribution);
	}

	if (!Communication::allGatherUnknownSize(_NLists)) {
		eslog::internalFailure("exchange ABAQUS NLists.\n");
	}

	if (!Communication::allGatherUnknownSize(_ELists)) {
		eslog::internalFailure("exchange ABAQUS ELists.\n");
	}

	if (!Communication::allGatherUnknownSize(_Esets)) {
		eslog::internalFailure("exchange ABAQUS ESets.\n");
	}

	if (!Communication::allGatherUnknownSize(_Nsets)) {
		eslog::internalFailure("exchange ABAQUS NSets.\n");
	}

	// fix distribution if EBlocks are across more processes and elements data have more lines
	for (size_t i = 0; i < _ELists.size(); i++) {
		_ELists[i].fixOffsets(_pfile.distribution);
		if (_pfile.begin != AbaqusParser::begin) {
			_pfile.distribution[info::mpi::rank] += _pfile.begin - AbaqusParser::begin;
		}
		if (_pfile.end != AbaqusParser::end) {
			_pfile.distribution[info::mpi::rank + 1] += _pfile.end - AbaqusParser::end;
		}
		_pfile.begin = AbaqusParser::begin;
		_pfile.end = AbaqusParser::end;
	}
}


void AbaqusLoader::parse()
{
	//_Elemats[0].create_dict(_SSections,_Materials);

	for (size_t i = 0; i < 1;i++){   //_NLists.size(); i++) {
		if (!_NLists[i].readData(nIDs, coordinates, 1)) {
			eslog::error("ABAQUS parser: something wrong happens while read NList.\n");
		}
	}

	for (size_t i = 0; i < _ELists.size(); i++) {
		if (!_ELists[i].readData(*this)) {
			eslog::error("ABAQUS parser: something wrong happens while read EList.\n");
		}
	}

	for (size_t i = 0; i < _Esets.size(); i++) {
		if (!_Esets[i].readData(eregions[_Esets[i].NAME])) {
			eslog::error("ABAQUS parser: something wrong happens while read ESets.\n");
		}
	}

	for (size_t i = 0; i < _Nsets.size(); i++) {
		if (!_Nsets[i].readData(nregions[_Nsets[i].NAME])) {
			eslog::error("ABAQUS parser: something wrong happens while read NSets.\n");
		}
	}
//	for (size_t i = 0; i < 1; i++) {
//		ESINFO(ERROR) << "Workbench parser: something wrong happens while read EBLOCK.";
//	}
}

