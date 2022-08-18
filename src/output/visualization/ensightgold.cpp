
#include "ensightgold.h"
#include "basis/containers/serializededata.h"
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
	}
	profiler::synccheckpoint("serialize");
	if (_measure) { eslog::checkpointln("ENSIGHT: GEOMETRY SERIALIZED"); }
	if (info::mpi::grank == 0) {
		casefile();
	}
	Communication::barrier(MPITools::asynchronous);
	profiler::synccheckpoint("store_casefile");
	if (_measure) { eslog::checkpointln("ENSIGHT: CASEFILE STORED"); }

	_writer.reorder();
	profiler::synccheckpoint("reorder");
	if (_measure) { eslog::checkpointln("ENSIGHT: GEOMETRY DATA REORDERED"); }

	_writer.write();
	profiler::synccheckpoint("write");
	Communication::barrier(MPITools::asynchronous);
	profiler::syncend("store_ensight");
	if (_measure) { eslog::endln("ENSIGHT: GEOMETRY STORED"); }
}

void EnSightGold::updateSolution()
{
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
	if (isRoot()) {
		_writer.format();
		_writer.description("EnSight Gold geometry format");
		_writer.description("----------------------------");

		if (_withIDs) {
			_writer.description("node id given");
		} else {
			_writer.description("node id off");
		}
		_writer.description("element id off");
	}

	auto nodes = [&] (const RegionStore *store) {
		if (isRoot()) {
			_writer.description("coordinates");
			_writer.int32ln(store->nodeInfo.totalSize);
		}

		if (_withIDs) {
			for (esint n = 0, i = store->nodeInfo.nhalo; n < store->nodeInfo.size; ++n, ++i) {
				_writer.int32ln(info::mesh->nodes->IDs->datatarray()[store->nodes->datatarray()[i]]);
			}
			_writer.groupData();
		}

		for (int d = 0; d < 3; ++d) {
			for (esint n = 0, i = store->nodeInfo.nhalo; n < store->nodeInfo.size; ++n, ++i) {
				_writer.float32ln(info::mesh->nodes->coordinates->datatarray()[store->nodes->datatarray()[i]][d]);
			}
			_writer.groupData();
		}
	};

	auto forEachElement = [] (const ElementsRegionStore *region, int etype, std::function<void(serializededata<esint, esint>::const_iterator)> callback) {
		for (size_t i = 0; i < region->eintervals.size(); i++) {
			if (region->eintervals[i].code == etype) {
				for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e) {
					callback(info::mesh->elements->nodes->cbegin() + region->elements->datatarray()[e]);
				}
			}
		}
	};

	esint part = 0;
	for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
		const ElementsRegionStore *region = info::mesh->elementsRegions[r];
		if (isRoot()) {
			_writer.description("part");
			_writer.int32ln(++part);
			_writer.description(region->name);
		}

		nodes(region);
		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (region->distribution.code[etype].totalSize) {
				if (isRoot()) {
					_writer.description(EnsightOutputWriter::codetotype(etype));
					_writer.int32ln(region->distribution.code[etype].totalSize);
				}

				if (Mesh::element(etype).code == Element::CODE::POLYGON) {
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) { _writer.int32ln(it->front()); });
					_writer.groupData();
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) {
						for (auto n = it->begin() + 1; n != it->end(); ++n) {
							_writer.int32(region->getPosition(*n) + 1);
						}
						_writer.ln();
					});
					_writer.groupData();
				}
				if (Mesh::element(etype).code == Element::CODE::POLYHEDRON) {
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) {
						_writer.int32ln(it->front());
					});
					_writer.groupData();
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) {
						PolyElement poly(Element::decode(Element::CODE::POLYHEDRON, it->size()), it->begin());
						for (size_t n = 1; n < it->size(); ++n) {
							if (!poly.isNode(n)) {
								_writer.int32ln(it->at(n));
							}
						}
					});
					_writer.groupData();
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) {
						PolyElement poly(Element::decode(Element::CODE::POLYHEDRON, it->size()), it->begin());
						for (size_t n = 2; n < it->size(); ++n) {
							if (poly.isNode(n)) {
								_writer.int32(region->getPosition(it->at(n)) + 1);
							} else {
								_writer.ln();
							}
						}
						_writer.ln();
					});
					_writer.groupData();
				}
				if (Mesh::element(etype).code != Element::CODE::POLYGON && Mesh::element(etype).code != Element::CODE::POLYHEDRON) {
					forEachElement(region, etype, [&] (serializededata<esint, esint>::const_iterator it) {
						for (auto n = it->begin(); n != it->end(); ++n) {
							_writer.int32(region->getPosition(*n) + 1);
						}
						_writer.ln();
					});
					_writer.groupData();
				}
			}
		}
	}

	auto boundary = [&] (const BoundaryRegionStore *region) {
		if (isRoot()) {
			_writer.description("part");
			_writer.int32ln(++part);
			_writer.description(region->name);
		}

		nodes(region);
		if (region->dimension) {
			for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
				if (region->distribution.code[etype].totalSize) {
					if (isRoot()) {
						_writer.description(EnsightOutputWriter::codetotype(etype));
						_writer.int32ln(region->distribution.code[etype].totalSize);
					}

					if (Mesh::element(etype).code == Element::CODE::POLYGON) {
						for (size_t i = 0; i < region->eintervals.size(); i++) {
							if (region->eintervals[i].code == etype) {
								for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e) {
									auto element = region->elements->cbegin() + e;
									_writer.int32ln(element->front());
								}
							}
						}
						_writer.groupData();
						for (size_t i = 0; i < region->eintervals.size(); i++) {
							if (region->eintervals[i].code == etype) {
								for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e) {
									auto element = region->elements->cbegin() + e;
									for (auto n = element->begin() + 1; n != element->end(); ++n) {
										_writer.int32(region->getPosition(*n) + 1);
									}
									_writer.ln();
								}
							}
						}
						_writer.groupData();
					} else {
						for (size_t i = 0; i < region->eintervals.size(); i++) {
							if (region->eintervals[i].code == etype) {
								for (esint e = region->eintervals[i].begin; e < region->eintervals[i].end; ++e) {
									auto element = region->elements->cbegin() + e;
									for (auto n = element->begin(); n != element->end(); ++n) {
										_writer.int32(region->getPosition(*n) + 1);
									}
									_writer.ln();
								}
							}
						}
						_writer.groupData();
					}
				}
			}
		} else {
			if (isRoot()) {
				_writer.description(EnsightOutputWriter::codetotype(static_cast<int>(Element::CODE::POINT1)));
				_writer.int32ln(region->nodeInfo.totalSize);
			}

			for (esint i = 0; i < region->nodeInfo.size; ++i) {
				_writer.int32ln(region->nodeInfo.offset + i + 1);
			}
			_writer.groupData();
		}
	};

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		boundary(info::mesh->boundaryRegions[r]);
	}

	for (size_t r = 0; r < info::mesh->contactInterfaces.size(); ++r) {
		boundary(info::mesh->contactInterfaces[r]);
	}

	_writer.commitFile(_path + _geometry);
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

