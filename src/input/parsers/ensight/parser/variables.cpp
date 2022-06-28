
#include "variables.h"
#include "casefile.h"

#include "basis/containers/serializededata.h"
#include "basis/io/inputfile.h"
#include "input/input.h"
#include "input/parsers/fileblock.h"
#include "input/parsers/distributedscanner.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"

using namespace espreso;

EnsightVariables::EnsightVariables(const EnsightCasefile &casefile, const EnsightGeometry &geofile, AsyncFilePack &variables)
: _casefile(casefile), _geofile(geofile), _variables(variables)
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
		_keywords.elements[e].ne = _geofile._keywords.elements[index].ne;
	}
}

void EnsightVariables::parse(Mesh &mesh)
{
	std::vector<OrderedValues> values(_variables.files.size());

	size_t c = 0, e = 0;
	while (_variables.next()) {
		size_t vsize = _keywords.header.format == EnsightKeywords::Format::ASCII ? 13 : sizeof(float);

		std::vector<size_t> offset, size;
		for ( ; c < _keywords.coordinates.size() && _keywords.coordinates[c].fileindex == _variables.fileindex; ++c) {
			offset.push_back(_keywords.coordinates[c].offset);
			size.push_back(_casefile.variables[_variables.fileindex].dimension * vsize * _keywords.coordinates[c].nn);
		}
		for ( ; e < _keywords.elements.size() && _keywords.elements[e].fileindex == _variables.fileindex; ++e) {
			offset.push_back(_keywords.elements[e].offset);
			size.push_back(_casefile.variables[_variables.fileindex].dimension * vsize * _keywords.elements[e].ne);
		}

		size_t localSize = 0;
		for (size_t i = 0, voffset = 0; i < offset.size(); voffset += size[i++]) {
			FileBlock block(_variables, offset[i], size[i], vsize, info::mpi::rank);
			localSize += block.size;
		}
		values[_variables.fileindex].data.reserve(localSize);

		for (size_t i = 0, voffset = 0; i < offset.size(); voffset += size[i++]) {
			FileBlock block(_variables, offset[i], size[i], vsize, info::mpi::rank);
			if (block.size) {
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					for (const char *cc = _variables.begin + block.begin; cc < _variables.begin + block.end; cc += vsize) {
						values[_variables.fileindex].data.push_back(atof(cc));
					}
				} else {
					for (const float *cc = (float*)(_variables.begin + block.begin); cc < (float*)(_variables.begin + block.end); ++cc) {
						values[_variables.fileindex].data.push_back(*cc);
					}
				}
			}
		}
	}

	for (size_t i = 0; i < _variables.files.size(); ++i) {
		switch (_casefile.variables[i].type) {
		case EnsightCasefile::Variable::Type::NODE: {
			NodeData *data = mesh.nodes->appendData(1, NamedData::DataType::SCALAR, _casefile.variables[i].name);
			for (size_t n = 0; n < mesh.nodes->IDs->datatarray().size(); ++n) {
				data->data[n] = values[i].data[mesh.nodes->IDs->datatarray()[n]];
			}
		} break;
		case EnsightCasefile::Variable::Type::ELEMENT: {
			ElementData *data = mesh.elements->appendData(1, NamedData::DataType::SCALAR, _casefile.variables[i].name);
			data->data.assign(values[i].data.begin(), values[i].data.end());
		} break;
		}
	}
}
