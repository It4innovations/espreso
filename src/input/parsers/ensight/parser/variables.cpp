
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

void EnsightVariables::parse(Mesh &mesh, VariablesBlocks &variables)
{
	size_t c = 0, e = 0, cv = 0, ev = 0;
	while (_variables.next()) {
		size_t vsize = _keywords.header.format == EnsightKeywords::Format::ASCII ? 13 : sizeof(float);
		for ( ; c < _keywords.coordinates.size() && _keywords.coordinates[c].fileindex == _variables.fileindex; ++c, ++cv) {
			esint offset = _keywords.coordinates[c].offset;
			esint size = _casefile.variables[_variables.fileindex].dimension * vsize * _keywords.coordinates[c].nn;
			FileBlock block(_variables, offset, size, vsize, info::mpi::rank);
			if (block.size) {
				variables.nodes[cv].data.reserve(block.size);
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					for (const char *cc = _variables.begin + block.begin; cc < _variables.begin + block.end; cc += vsize) {
						variables.nodes[cv].data.push_back(atof(cc));
					}
				} else {
					for (const float *cc = (float*)(_variables.begin + block.begin); cc < (float*)(_variables.begin + block.end); ++cc) {
						variables.nodes[cv].data.push_back(*cc);
					}
				}
			}
		}
		for ( ; e < _keywords.elements.size() && _keywords.elements[e].fileindex == _variables.fileindex; ++e, ++ev) {
			esint offset = _keywords.elements[e].offset;
			esint size = _casefile.variables[_variables.fileindex].dimension * vsize * _keywords.elements[e].ne;
			FileBlock block(_variables, offset, size, vsize, info::mpi::rank);
			if (block.size) {
				variables.elements[ev].data.reserve(block.size);
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					for (const char *cc = _variables.begin + block.begin; cc < _variables.begin + block.end; cc += vsize) {
						variables.elements[ev].data.push_back(atof(cc));
					}
				} else {
					for (const float *cc = (float*)(_variables.begin + block.begin); cc < (float*)(_variables.begin + block.end); ++cc) {
						variables.elements[ev].data.push_back(*cc);
					}
				}
			}
		}
	}
}
