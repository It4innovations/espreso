
#include "ensightgold.h"
#include "basis/containers/serializededata.h"
#include "basis/io/outfile.h"
#include "basis/logging/profiler.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "mesh/mesh.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/clusterstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/contactinterfacestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/contactstore.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>

using namespace espreso;

EnSightGold::EnSightGold()
: _ftt(NULL), _withIDs(false), _step(-1)
{
	_geometry = _directory + _name + ".geo";
	_fixedDataPath = _directory;
}

EnSightGold::~EnSightGold()
{

}

EnSightGold::FTT::FTT(EnSightGold *parent)
: EnSightGold()
{
	_directory += std::to_string(step::outfrequency.current) + "/";
	_name = parent->_name + ".freq." + std::to_string(step::outfrequency.current);
	if (isRoot()) {
		createOutputDirectory();
	}
}

void EnSightGold::updateMesh()
{
	if (_measure) { eslog::startln("ENSIGHT: STORING STARTED", "ENSIGHT"); }
	profiler::syncstart("store_ensight");
	if (info::mpi::irank == 0) {
		geometry();
		if (info::ecf->output.store_decomposition) {
			decomposition();
		}
		casefile();
	}
	profiler::syncend("store_ensight");
	if (_measure) { eslog::endln("ENSIGHT: MESH STORED"); }
}

void EnSightGold::updateSolution()
{
	return;
	if (_measure) { eslog::startln("ENSIGHT RESULTS: STORING STARTED", "ENSIGHT RESULTS"); }
	EnSightGold *writer = this;

	if (step::outstep.type == step::TYPE::FTT && step::outftt.isFirst()) {
		_ftt = new FTT(this);
	}

	switch (step::outstep.type) {
	case step::TYPE::TIME:      _times.push_back(step::outtime.current); break;
	case step::TYPE::FREQUENCY: _times.push_back(step::outfrequency.current); break;
	case step::TYPE::FTT: _ftt->_times.push_back(step::outftt.time); writer = _ftt; break;
	}

	size_t nvars = 0;
	for (size_t di = 0; di < info::mesh->elements->data.size(); di++) {
		nvars += writer->edata(info::mesh->elements->data[di]);
	}
	for (size_t di = 0; di < info::mesh->nodes->data.size(); di++) {
		nvars += writer->ndata(info::mesh->nodes->data[di]);
	}
	if (nvars != _variables.size()) {
		_variables.clear();
	}

	if (_step != step::outstep.loadstep || info::mpi::grank == 0 || (step::outstep.type == step::TYPE::FTT && info::mpi::rank == 0)) {
		_step = step::outstep.loadstep;
		if (step::outduplicate.instances > 1 && step::outstep.substep - step::outduplicate.offset + 1 == step::outduplicate.size) {
			// TODO: generalize
			_times.clear();
			for (auto f = step::outfrequency.start; f < step::outfrequency.final; f += step::outfrequency.shift) {
				_times.push_back(f + step::outfrequency.shift);
			}
		}

		auto spaces = [] (const std::string &label, size_t size) -> std::string {
			if (size > label.size()) {
				return std::string(size - label.size(), ' ');
			}
			return "";
		};

		auto pushdata = [&] (std::vector<std::string> &variables, const NamedData *data, const std::string &var) {
			if (!storeData(data)) {
				return;
			}
			if (data->dataType == NamedData::DataType::VECTOR) {
				std::string name = dataname(data, 0);
				variables.push_back(
						"vector per " + var + ": " + spaces(var, 8) + "1 " + name + spaces(name, 30) + " " + writer->_directory + name + ".****"
				);
				return;
			}
			if (data->dimension == 1) {
				std::string name = dataname(data, 0);
				variables.push_back(
						"scalar per " + var + ": " + spaces(var, 8) + "1 " + name + spaces(name, 30) + " " + writer->_directory + name + ".****"
				);
				return;
			}
			for (int d = 0; d < data->dimension; d++) {
				std::string name = dataname(data, d);
				variables.push_back(
						"scalar per " + var + ": " + spaces(var, 8) + "1 " + name + spaces(name, 30) + " " + writer->_directory + name + ".****"
				);
			}
		};

		writer->_variables.clear();
		for (size_t i = 0; i < info::mesh->nodes->data.size(); i++) {
			pushdata(writer->_variables, info::mesh->nodes->data[i], "node");
		}
		for (size_t i = 0; i < info::mesh->elements->data.size(); i++) {
			pushdata(writer->_variables, info::mesh->elements->data[i], "element");
		}

		writer->casefile();
	}

	if (_measure) { eslog::checkpointln("ENSIGHT RESULTS: DATA INSERTED"); }

	writer->_writer.reorder();
	if (_measure) { eslog::checkpointln("ENSIGHT RESULTS: DATA REORDERED"); }

	writer->_writer.write();
	if (step::outstep.type == step::TYPE::FTT && step::outftt.isLast()) {
		delete _ftt;
	}

	Communication::barrier(MPITools::asynchronous);
	if (_measure) { eslog::endln("ENSIGHT RESULTS: DATA STORED"); }
}

std::string EnSightGold::dataname(const NamedData *data, int d)
{
	if (data->dataType == NamedData::DataType::VECTOR) {
		return data->name;
	}
	if (data->dimension == 1) {
		return data->name;
	}
	return data->name + data->suffix(d);
}

void EnSightGold::casefile()
{
	std::ofstream os(_path + _name + ".case");

	os << "\n# output from ESPRESO library (espreso.it4.cz)";
	os << "\n#";
	os << "\n";
	os << "\nFORMAT";
	os << "\ntype: ensight gold";
	os << "\n";
	os << "\nGEOMETRY";
	os << "\nmodel: " << _geometry;
	os << "\n";
	os << "\nVARIABLE";
	os << "\n";
	if (info::ecf->output.store_decomposition) {
		os << "scalar per element: BODY    " << _fixedDataPath << "BODY" << "\n";
		os << "scalar per element: DOMAIN  " << _fixedDataPath << "DOMAIN" << "\n";
		os << "scalar per element: CLUSTER " << _fixedDataPath << "CLUSTER" << "\n";
		os << "scalar per element: MPI     " << _fixedDataPath << "MPI" << "\n";
	}
	if (_times.size()) {
		for (size_t i = 0; i < _variables.size(); ++i) {
			os << "\n" << _variables[i];
		}
		os << "\n";
		os << "\nTIME";
		os << "\ntime set:               1";
		os << "\nnumber of steps:        " << _times.size();
		os << "\nfilename start numbers: 1";
		os << "\nfilename increment:     1";
		os << "\ntime values:";
		for (size_t i = 0; i < _times.size(); ++i) {
			if (i % 10 == 0) {
				os << "\n" << _times[i];
			} else {
				os << " " << _times[i];
			}
		}
	}
}

void EnSightGold::geometry()
{
	OutFile file;

	size_t fileHeader = 5 * 80;
	size_t partHeader = 2 * 80 + sizeof(int);
	size_t coorHeader = 1 * 80 + sizeof(int);
	size_t elemHeader = 1 * 80 + sizeof(int);

	std::vector<const RegionStore*> regions;
	regions.insert(regions.end(), info::mesh->elementsRegions.begin() + 1, info::mesh->elementsRegions.end());
	regions.insert(regions.end(), info::mesh->boundaryRegions.begin() + 1, info::mesh->boundaryRegions.end());
	regions.insert(regions.end(), info::mesh->contactInterfaces.begin(), info::mesh->contactInterfaces.end());

	size_t blocks = 0;
	for (size_t r = 0; r < regions.size(); ++r) {
		blocks += 3;
		for (size_t c = 0; c < regions[r]->distribution.code.size(); ++c) {
			if (regions[r]->distribution.code[c].totalSize) {
				switch (c) {
				case (size_t)Element::CODE::POLYGON: blocks += 2; break;
				case (size_t)Element::CODE::POLYHEDRON: blocks += 3; break;
				default: blocks += 1;
				}
			}
		}
	}
	file.blocks.resize(blocks);

	size_t block = 0, offset = fileHeader;
	std::vector<esint> rbegin(regions.size() + 1, blocks);
	for (size_t r = 0; r < regions.size(); ++r) {
		rbegin[r] = block;
		offset += partHeader + coorHeader;
		for (int d = 0; d < 3; ++d) {
			file.blocks[block + d].fileoffset = offset + sizeof(float) * info::mesh->elementsRegions[r]->nodeInfo.offset;
			file.blocks[block + d].size = sizeof(float) * info::mesh->elementsRegions[r]->nodeInfo.size;
			offset += sizeof(float) * info::mesh->elementsRegions[r]->nodeInfo.totalSize;
		}
		if (isRoot()) {
			file.blocks[block].fileoffset -= partHeader + coorHeader;
			file.blocks[block].size += partHeader + coorHeader;
		}
		block += 3;

		for (size_t c = 0; c < regions[r]->distribution.code.size(); ++c) {
			if (regions[r]->distribution.code[c].totalSize) {
				switch (c) {
				case (size_t)Element::CODE::POLYGON:
					if (isRoot()) {
						file.blocks[block].fileoffset = offset;
						file.blocks[block].size = elemHeader + sizeof(int) * regions[r]->distribution.code[c].size;
					} else {
						file.blocks[block].fileoffset = offset + elemHeader + sizeof(int) * regions[r]->distribution.code[c].offset;
						file.blocks[block].size = sizeof(int) * regions[r]->distribution.code[c].size;
					}
					offset += elemHeader + sizeof(int) * regions[r]->distribution.code[c].totalSize;
					block += 1;

					file.blocks[block].fileoffset = offset + sizeof(int) * regions[r]->distribution.polygonNodes.offset;
					file.blocks[block].size = sizeof(int) * regions[r]->distribution.polygonNodes.size;
					offset += sizeof(int) * regions[r]->distribution.polygonNodes.totalSize;
					block += 1;
					break;
				case (size_t)Element::CODE::POLYHEDRON:
					if (isRoot()) {
						file.blocks[block].fileoffset = offset;
						file.blocks[block].size = elemHeader + sizeof(int) * regions[r]->distribution.code[c].size;
					} else {
						file.blocks[block].fileoffset = offset + elemHeader + sizeof(int) * regions[r]->distribution.code[c].offset;
						file.blocks[block].size = sizeof(int) * regions[r]->distribution.code[c].size;
					}
					offset += elemHeader + sizeof(int) * regions[r]->distribution.code[c].totalSize;
					block += 1;

					file.blocks[block].fileoffset = offset + sizeof(int) * regions[r]->distribution.polyhedronFaces.offset;
					file.blocks[block].size = sizeof(int) * regions[r]->distribution.polyhedronFaces.size;
					offset += sizeof(int) * regions[r]->distribution.polyhedronFaces.totalSize;
					block += 1;

					file.blocks[block].fileoffset = offset + sizeof(int) * regions[r]->distribution.polyhedronNodes.offset;
					file.blocks[block].size = sizeof(int) * regions[r]->distribution.polyhedronNodes.size;
					offset += sizeof(int) * regions[r]->distribution.polyhedronNodes.totalSize;
					block += 1;
					break;
				default:
					if (isRoot()) {
						file.blocks[block].fileoffset = offset;
						file.blocks[block].size = elemHeader + sizeof(int) * regions[r]->distribution.code[c].size * Mesh::edata[c].nodes;
					} else {
						file.blocks[block].fileoffset = offset + elemHeader + sizeof(int) * regions[r]->distribution.code[c].offset * Mesh::edata[c].nodes;
						file.blocks[block].size = sizeof(int) * regions[r]->distribution.code[c].size * Mesh::edata[c].nodes;
					}
					offset += elemHeader + sizeof(int) * regions[r]->distribution.code[c].totalSize * Mesh::edata[c].nodes;
					block += 1;
				}
			}
		}
	}
	if (isRoot()) {
		file.blocks.front().fileoffset = 0;
		file.blocks.front().size += fileHeader;
	}

	file.prepare();
	if (_measure) { eslog::checkpointln("ENSIGHT: BLOCK ALLOCATED"); }

//	Communication::serialize([&] () {
//		for (size_t b = 0; b < file.blocks.size(); ++b) {
//			printf("<%9lu - %9lu> (%9lu %9lu)\n", file.blocks[b].fileoffset, file.blocks[b].fileoffset + file.blocks[b].size, file.blocks[b].size, (file.blocks[b].size - 84) / 4);
//		}
//	}, MPITools::asynchronous);

	file.open(_path + _geometry);
	if (_measure) { eslog::checkpointln("ENSIGHT: MESH FILE OPENED"); }

	auto check = [] (size_t region, const esint *mask) {
		esint maskOffset = region / (8 * sizeof(esint));
		esint bit = (esint)1 << (region % (8 * sizeof(esint)));
		return mask[maskOffset] & bit;
	};

	auto description = [] (char *data, const std::string &desc) {
		memcpy(data, desc.data(), desc.size());
		memset(data + desc.size(), ' ', 80 - desc.size());
		data[79] = '\n';
	};
	auto descriptionWithSize = [] (char *data, const std::string &desc, int size) {
		memcpy(data, desc.data(), desc.size());
		memset(data + desc.size(), ' ', 80 - desc.size());
		data[79] = '\n';
		memcpy(data + 80, &size, sizeof(int));
	};

	auto nranks = info::mesh->nodes->ranks->cbegin();
	auto bmask = info::mesh->nodes->bregions->cbegin();
	auto emask = info::mesh->nodes->eregions->cbegin();
	auto fblocks = file.blocks;

	if (isRoot()) {
		description(file.data.data()      , "C Binary");
		description(file.data.data() + 80 , "ESPRESO serialized output");
		description(file.data.data() + 160, "-------------------------");
		description(file.data.data() + 240, "node id off");
		description(file.data.data() + 320, "element id off");
		fblocks[0].offset += 400;
		for (size_t r = 0; r < regions.size(); ++r) {
			descriptionWithSize(file.data.data() + fblocks[rbegin[r]].offset, "part", r + 1);
			fblocks[rbegin[r]].offset += 80 + sizeof(int);
			description(file.data.data() + fblocks[rbegin[r]].offset, regions[r]->name);
			fblocks[rbegin[r]].offset += 80;
			descriptionWithSize(file.data.data() + fblocks[rbegin[r]].offset, "coordinates", regions[r]->nodeInfo.totalSize);
			fblocks[rbegin[r]].offset += 80 + sizeof(int);
		}
	}
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks, ++bmask, ++emask) {
		auto insert = [&file, &fblocks] (_Point<float> p, const size_t &bindex) {
			memcpy(file.data.data() + fblocks[bindex + 0].offset, &p.x, sizeof(float));
			memcpy(file.data.data() + fblocks[bindex + 1].offset, &p.y, sizeof(float));
			memcpy(file.data.data() + fblocks[bindex + 2].offset, &p.z, sizeof(float));
			fblocks[bindex + 0].offset += sizeof(float);
			fblocks[bindex + 1].offset += sizeof(float);
			fblocks[bindex + 2].offset += sizeof(float);
		};

		if (nranks->front() == info::mpi::rank) {
			size_t bindex = 0;
			for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++bindex) {
				if (check(r, emask->data())) {
					insert(info::mesh->nodes->coordinates->datatarray()[n], rbegin[bindex]);
				}
			}
			size_t ri = 1;
			for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++ri, ++bindex) {
				if (check(ri, bmask->data())) {
					insert(info::mesh->nodes->coordinates->datatarray()[n], rbegin[bindex]);
				}
			}
			for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++ri, ++bindex) {
				if (check(ri, bmask->data())) {
					insert(info::mesh->nodes->coordinates->datatarray()[n], rbegin[bindex]);
				}
			}
		}
	}

	for (size_t r = 0; r < regions.size(); ++r) {
		for (int b = 0; b < 3; ++b) {
			file.store(rbegin[r] + b);
		}
	}
	if (_measure) { eslog::checkpointln("ENSIGHT: COORDINATES SERIALIZED"); }

	std::vector<esint> output(info::mesh->nodes->outputOffset->boundarytarray().begin(), info::mesh->nodes->outputOffset->boundarytarray().end());

	auto eheader = [&file, &fblocks, &descriptionWithSize] (const RegionStore *region, const size_t &bindex, int code) {
		if (isRoot()) {
			descriptionWithSize(file.data.data() + fblocks[bindex].offset, EnsightOutputWriter::codetotype(code), region->distribution.code[code].totalSize);
			fblocks[bindex].offset += 80 + sizeof(int);
		}
	};

	auto node = [&file, &fblocks] (const size_t &bindex, int node) {
		memcpy(file.data.data() + fblocks[bindex].offset, &node, sizeof(int));
		fblocks[bindex].offset += sizeof(int);
	};

	auto eregion = [&] (const ElementsRegionStore *region, const size_t &bindex, Element::CODE code) {
		const auto &ep = info::mesh->elements->epointers->datatarray();
		auto enodes = info::mesh->elements->nodes->cbegin();
		esint prev = 0;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			if (ep[*e]->code == code) {
				enodes += *e - prev; prev = *e;
				for (auto n = enodes->begin(); n != enodes->end(); ++n) {
					node(bindex, info::mesh->nodes->outputOffset->datatarray()[output[*n]] + 1);
				}
			}
		}
	};

	auto epolygons = [&] (const ElementsRegionStore *region, const size_t &bindex) {
		const auto &ep = info::mesh->elements->epointers->datatarray();
		auto enodes = info::mesh->elements->nodes->cbegin();
		esint prev = 0;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			if (ep[*e]->code == Element::CODE::POLYGON) {
				enodes += *e - prev; prev = *e;
				node(bindex, enodes->front());
				for (auto n = enodes->begin() + 1; n != enodes->end(); ++n) {
					node(bindex + 1, info::mesh->nodes->outputOffset->datatarray()[output[*n]] + 1);
				}
			}
		}
	};

	auto epolyhedrons  = [&] (const ElementsRegionStore *region, const size_t &bindex) {
		const auto &ep = info::mesh->elements->epointers->datatarray();
		auto enodes = info::mesh->elements->nodes->cbegin();
		esint prev = 0;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			if (ep[*e]->code == Element::CODE::POLYHEDRON) {
				enodes += *e - prev; prev = *e;
				PolyElement poly(Element::decode(Element::CODE::POLYHEDRON, enodes->size()), enodes->begin());
				node(bindex, enodes->front());
				for (size_t n = 1; n < enodes->size(); ++n) {
					if (poly.isNode(n)) {
						node(bindex + 2, info::mesh->nodes->outputOffset->datatarray()[output[enodes->at(n)]] + 1);
					} else {
						node(bindex + 1, enodes->at(n));
					}
				}
			}
		}
	};

	auto bregion = [&file, &fblocks] (const BoundaryRegionStore *region) {

	};

	auto nextOutput = [] () {

	};

	size_t ri = 0;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++ri) {
		for (size_t c = 0, bindex = rbegin[ri] + 3; c < regions[ri]->distribution.code.size(); ++c) {
			if (regions[ri]->distribution.code[c].totalSize) {
				eheader(regions[ri], bindex, c);
				switch (c) {
				case (int)Element::CODE::POLYGON: epolygons(info::mesh->elementsRegions[r], bindex); bindex += 2; break;
				case (int)Element::CODE::POLYHEDRON: epolyhedrons(info::mesh->elementsRegions[r], bindex); bindex += 3; break;
				default: eregion(info::mesh->elementsRegions[r], bindex, (Element::CODE)c); bindex += 1;
				}
				nextOutput();
			}
		}
	}
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++ri) {
		for (size_t c = 0, bindex = rbegin[ri] + 3; c < regions[ri]->distribution.code.size(); ++c) {
			if (regions[ri]->distribution.code[c].totalSize) {
				eheader(regions[ri], bindex, c);
				bregion(info::mesh->boundaryRegions[r]);
				nextOutput();
				++bindex;
			}
		}
	}
	for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++ri) {
		for (size_t c = 0, bindex = rbegin[ri] + 3; c < regions[ri]->distribution.code.size(); ++c) {
			if (regions[ri]->distribution.code[c].totalSize) {
				eheader(regions[ri], bindex, c);
				bregion(info::mesh->contactInterfaces[r]);
				nextOutput();
				++bindex;
			}
		}
	}

	for (size_t r = 0; r < regions.size(); ++r) {
		for (esint b = rbegin[r] + 3; b < rbegin[r + 1]; ++b) {
			file.store(b);
		}
	}
	if (_measure) { eslog::checkpointln("ENSIGHT: ELEMENTS SERIALIZED"); }
	file.close();
}

int EnSightGold::ndata(const NamedData *data)
{
	if (!storeData(data)) {
		return 0;
	}

	auto niterator = [this] (esint size, esint *nodes, std::function<void(esint nindex)> callback) {
		for (esint n = 0; n < size; ++n) {
			callback(nodes[n]);
		}
		_writer.groupData();
	};

	if (data->dataType == NamedData::DataType::VECTOR) {
		std::stringstream file;
		file << _path + _directory + dataname(data, 0) + "." << std::setw(4) << std::setfill('0') << _times.size() + step::outduplicate.offset;

		if (isRoot()) {
			_writer.description(dataname(data, 0));
		}

		int part = 1;
		for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++part) {
			if (isRoot()) {
				_writer.description("part");
				_writer.int32ln(part);
				_writer.description("coordinates");
			}

			const ElementsRegionStore *region = info::mesh->elementsRegions[r];
			for (int d = 0; d < data->dimension; ++d) {
				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(data->store[nindex * data->dimension + d]);
				});
			}
			if (data->dimension == 2) {
				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(0);
				});
			}
		}

		auto boundary = [&] (const BoundaryRegionStore *region) {
			if (isRoot()) {
				_writer.description("part");
				_writer.int32ln(part);
				_writer.description("coordinates");
			}

			for (int d = 0; d < data->dimension; ++d) {
				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(data->store[nindex * data->dimension + d]);
				});
			}
			if (data->dimension == 2) {
				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(0);
				});
			}
		};

		for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++part) {
			boundary(info::mesh->boundaryRegions[r]);
		}
		for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++part) {
			boundary(info::mesh->contactInterfaces[r]);
		}
		_writer.commitFile(file.str());
	} else {
		for (int d = 0; d < data->dimension; ++d) {
			std::stringstream file;
			file << _path + _directory + dataname(data, d) + "." << std::setw(4) << std::setfill('0') << _times.size() + step::outduplicate.offset;

			if (isRoot()) {
				_writer.description(dataname(data, d));
			}

			int part = 1;
			for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r, ++part) {
				if (isRoot()) {
					_writer.description("part");
					_writer.int32ln(part);
					_writer.description("coordinates");
				}

				const ElementsRegionStore *region = info::mesh->elementsRegions[r];
				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(data->store[nindex * data->dimension + d]);
				});
			}

			auto boundary = [&] (const BoundaryRegionStore *region) {
				if (isRoot()) {
					_writer.description("part");
					_writer.int32ln(part);
					_writer.description("coordinates");
				}

				niterator(region->nodeInfo.size, region->nodes->datatarray().data() + region->nodeInfo.nhalo, [&] (esint nindex) {
					_writer.float32ln(data->store[nindex * data->dimension + d]);
				});
			};

			for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r, ++part) {
				boundary(info::mesh->boundaryRegions[r]);
			}
			for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r, ++part) {
				boundary(info::mesh->contactInterfaces[r]);
			}
			_writer.commitFile(file.str());
		}
	}
	return 1;
}

int EnSightGold::edata(const NamedData *data)
{
	if (!storeData(data)) {
		return 0;
	}

	auto eiterator = [&] (const ElementsRegionStore *region, int etype, std::function<void(const ElementsInterval &interval, esint eindex)> callback) {
		for (size_t i = 0; i < region->eintervals.size(); i++) {
			if (region->eintervals[i].code == etype) {
				for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e) {
					callback(region->eintervals[i], region->elements->datatarray()[e]);
				}
			}
		}
		_writer.groupData();
	};

	if (data->dataType == NamedData::DataType::VECTOR) {
		std::stringstream file;
		file << _path + _directory + dataname(data, 0) + "." << std::setw(4) << std::setfill('0') << _times.size() + step::outduplicate.offset;
		if (isRoot()) {
			_writer.description(dataname(data, 0));
		}

		for (size_t r = 1, part = 1; r < info::mesh->elementsRegions.size(); ++r, ++part) {
			if (isRoot()) {
				_writer.description("part");
				_writer.int32ln(part);
			}

			for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
				if (info::mesh->elementsRegions[r]->distribution.code[etype].totalSize) {
					if (isRoot()) {
						_writer.description(EnsightOutputWriter::codetotype(etype));
					}

					for (int d = 0; d < data->dimension; ++d) {
						eiterator(info::mesh->elementsRegions[r], etype, [&] (const ElementsInterval &interval, esint eindex) {
							_writer.float32ln(data->store[eindex * data->dimension + d]);
						});
					}
					if (data->dimension == 2) {
						eiterator(info::mesh->elementsRegions[r], etype, [&] (const ElementsInterval &interval, esint eindex) {
							_writer.float32ln(0);
						});
					}
				}
			}
		}
		_writer.commitFile(file.str());
	} else {
		for (int d = 0; d < data->dimension; ++d) {
			std::stringstream file;
			file << _path + _directory + dataname(data, d) + "." << std::setw(4) << std::setfill('0') << _times.size() + step::outduplicate.offset;

			if (isRoot()) {
				_writer.description(dataname(data, 0));
			}

			for (size_t r = 1, part = 1; r < info::mesh->elementsRegions.size(); ++r, ++part) {
				if (isRoot()) {
					_writer.description("part");
					_writer.int32ln(part);
				}
				for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
					if (info::mesh->elementsRegions[r]->distribution.code[etype].totalSize) {
						if (isRoot()) {
							_writer.description(EnsightOutputWriter::codetotype(etype));
						}
						eiterator(info::mesh->elementsRegions[r], etype, [&] (const ElementsInterval &interval, esint eindex) {
							_writer.float32ln(data->store[eindex * data->dimension + d]);
						});
					}
				}
			}
			_writer.commitFile(file.str());
		}
	}
	return 1;
}

void EnSightGold::decomposition()
{
	auto store = [&] (const std::string &name, std::function<double(const ElementsInterval &interval, esint eindex)> callback) {
		if (isRoot()) {
			_writer.description(name);
		}
		for (size_t r = 1, part = 1; r < info::mesh->elementsRegions.size(); ++r, ++part) {
			if (isRoot()) {
				_writer.description("part");
				_writer.int32ln(part);
			}

			for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
				if (info::mesh->elementsRegions[r]->distribution.code[etype].totalSize) {
					if (isRoot()) {
						_writer.description(EnsightOutputWriter::codetotype(etype));
					}
					for (size_t i = 0; i < info::mesh->elementsRegions[r]->eintervals.size(); i++) {
						if (info::mesh->elementsRegions[r]->eintervals[i].code == etype) {
							for (esint e = info::mesh->elementsRegions[r]->eintervals[i].begin; e < info::mesh->elementsRegions[r]->eintervals[i].end; ++e) {
								_writer.float32ln(callback(info::mesh->elementsRegions[r]->eintervals[i], info::mesh->elementsRegions[r]->elements->datatarray()[e]));
							}
						}
					}
					_writer.groupData();
				}
			}
		}
		_writer.commitFile(_path + _directory + name);
	};

	store("DOMAIN", [&] (const ElementsInterval &interval, esint eindex) {
		return interval.domain;
	});
	if (info::mesh->domains->cluster.size()) {
		esint cluster = info::mesh->clusters->offset;
		store("CLUSTER", [&] (const ElementsInterval &interval, esint eindex) {
			return info::mesh->domains->cluster[interval.domain - info::mesh->domains->offset] + cluster;
		});
	}
	store("MPI", [&] (const ElementsInterval &interval, esint eindex) {
		return info::mpi::rank;
	});
	store("BODY", [&] (const ElementsInterval &interval, esint eindex) {
		return info::mesh->elements->body->datatarray()[eindex];
	});
}

