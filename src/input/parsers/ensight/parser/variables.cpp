
#include "variables.h"
#include "casefile.h"

#include "basis/io/inputfile.h"
#include "input/input.h"
#include "input/parsers/fileblock.h"
#include "input/parsers/distributedscanner.h"

using namespace espreso;

EnsightVariables::EnsightVariables(const EnsightCasefile &casefile, const EnsightGeometry &geofile, InputFilePack &variables, OrderedMeshDatabase &database)
: _casefile(casefile), _geofile(geofile), _variables(variables), _database(database)
{
	_keywords.header = _geofile._keywords.header;
}

void EnsightVariables::scan()
{
	DistributedScanner scanner;
	if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
		fillScanner<EnsightASCIIVariableKeywordParser>(_variables, scanner, _keywords);
	} else {
		fillScanner<EnsightBinaryVariableKeywordParser>(_variables, scanner, _keywords);
	}

	while (_variables.next()) {
		scanner.scan(_variables);
	}

	scanner.synchronize(_keywords.parts, _keywords.coordinates, _keywords.elements);
	scanner.sort(_keywords.parts, _keywords.coordinates, _keywords.elements);

	for (size_t c = 0, file = 0, index = 0; c < _keywords.coordinates.size(); ++c, ++index) {
		if (file != _keywords.coordinates[c].fileindex) {
			file = _keywords.coordinates[c].fileindex;
			index = 0;
		}
		_keywords.coordinates[c].nn = _geofile._keywords.coordinates[index].nn;
	}
	for (size_t e = 0, file = 0, index = 0; e < _keywords.elements.size(); ++e, ++index) {
		if (file != _keywords.elements[e].fileindex) {
			file = _keywords.elements[e].fileindex;
			index = 0;
		}
		_keywords.elements[e].ne = _geofile._keywords.elements[e].ne;
	}
}

void EnsightVariables::parse()
{
	size_t c = 0, e = 0;
	while (_variables.next()) {
		std::vector<OrderedMeshDatabase::Values> *_values;
		switch (_casefile.variables[_variables.fileindex].type) {
		case EnsightCasefile::Variable::Type::NODE:    _values = &_database.nvalues; break;
		case EnsightCasefile::Variable::Type::ELEMENT: _values = &_database.evalues; break;
		}
		_values->push_back({_casefile.variables[_variables.fileindex].name, _casefile.variables[_variables.fileindex].dimension});
		OrderedMeshDatabase::Values &values = _values->back();

		size_t vsize;
		if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
			vsize = 13;
		} else {
			vsize = sizeof(float);
		}

		std::vector<size_t> offset, size;
		for ( ; c < _keywords.coordinates.size() && _keywords.coordinates[c].fileindex == _variables.fileindex; ++c) {
			offset.push_back(_keywords.coordinates[c].offset);
			size.push_back(values.dimension * vsize * _keywords.coordinates[c].nn);
		}
		for ( ; e < _keywords.elements.size() && _keywords.elements[e].fileindex == _variables.fileindex; ++e) {
			offset.push_back(_keywords.elements[e].offset);
			size.push_back(values.dimension * vsize * _keywords.elements[e].ne);
		}
		for (size_t i = 0, voffset = 0; i < offset.size(); voffset += size[i++]) {
			std::vector<float> data;
			FileBlock block(_variables, offset[i], size[i], vsize, info::mpi::rank);
			if (block.size) {
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					data.reserve(block.size / vsize);
					for (const char *cc = _variables.begin + block.begin; cc < _variables.begin + block.end; cc += vsize) {
						data.push_back(atof(cc));
					}
				} else {
					data.resize(block.size / vsize);
					memcpy(data.data(), _variables.begin + block.begin, block.size);
				}

				esint coffset = voffset / vsize + block.prevsize / vsize;
				if (values.dimension == 1) {
					values.values.push_back({coffset, values.dimension});
					values.values.back().values.assign(data.begin(), data.end());
				} else {
					eslog::internalFailure("Implement parsing of multidimensional data.\n");
				}
			}
		}
	}
}
