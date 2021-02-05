
#include "eblock.h"
#include "et.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "wrappers/mpi/communication.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/element.h"
#include "input/parsers/ansyscdb/ansyscdb.h"

#include <algorithm>

using namespace espreso;

EBlock::EBlock()
: NUM_NODES(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineEndSize(1),
  valueSize(-1), valueLength(-1)
{

}

EBlock& EBlock::parse(const char* begin)
{
	auto error = [&] (std::string &line) {
		eslog::error("Workbench parse error: unknown format of EBLOCK: %s\n", line.c_str());
	};

	std::string commandLine = Parser::getLine(begin);
	std::string descriptionLine = Parser::getLine(begin + commandLine.size());

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> description = Parser::split(descriptionLine.substr(1, descriptionLine.size() - 2 - lineEndSize), "i");
	if (description.size() != 2) {
		error(descriptionLine);
	}

	switch (command.size()) {
	case 5:
		NDSEL = command[4].size() ? std::stol(command[4]) : -1;
	case 4:
		NDMAX = command[3].size() ? std::stol(command[3]) : -1;
	case 3:
		NUM_NODES = std::stoi(command[1]);
		Solkey = command[2].size();
		break;
	default:
		error(commandLine);
	}

	valueSize = std::stoi(description[0]);
	valueLength = std::stoi(description[1]);

	WorkbenchParser::fillIndices(begin, begin + commandLine.size() + descriptionLine.size());
	return *this;
}

bool EBlock::readData(const std::vector<ET> &et, AnsysCDBData &mesh)
{
	if (Solkey) {
		return solid(et, mesh);
	} else {
		return boundary(et, mesh);
	}
}

bool EBlock::solid(const std::vector<ET> &et, AnsysCDBData &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, last - first);

	std::vector<std::vector<esint> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tet(threads), tbody(threads), tmat(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<esint> nindices(20);

		std::vector<esint> esize, nodes, IDs;
		std::vector<int> type, ansystype, body, mat;

		auto data = first + tdistribution[t];
		if (tdistribution[t] && *(data - 1) != '\n') {
			while (data < first + tdistribution[t + 1] && *data++ != '\n'); // start at new line
		}

		{ // check if the first line contains an element description
			std::vector<esint> fields(11);
			const char *line = data;
			std::fill(fields.begin(), fields.end(), -1);
			for (int i = 0; i < 11; ++i) {
				if (*line != '\r' && *line != '\n') {
					memcpy(value.data(), line, valueLength);
					line += valueLength;
					fields[i] = atol(value.data());
				}
			};
			if (
					fields[8] > 20                 || // too much of nodes
					fields[1] > (esint)et.size()   || // invalid element type
					fields[5] < 0 || fields[5] > 1 || // birth flag (zero is treated later)
					std::any_of(fields.begin(), fields.end(), [] (int i) { return i == -1;})
			) {
				while (data < first + tdistribution[t + 1] && *data++ != '\n'); // start at the next line
			}
		}

		auto skip = [&] () {
			data += valueLength;
		};

		auto parse = [&] () {
			memcpy(value.data(), data, valueLength);
			skip();
			return atol(value.data());
		};

		while (data < first + tdistribution[t + 1]) {
			body.push_back(0);
			mat.push_back(parse() - 1); // material
			ansystype.push_back(parse() - 1); // etype
			type.push_back(0);
			skip(); // real constant
			skip(); // section ID
			skip(); // element coordinate system
			skip(); // birth / death
			skip(); // solid model reference number
			skip(); // element shape flag
			int nnodes = parse(); // number of nodes
			skip(); // not used
			IDs.push_back(parse() - 1); // element ID

			for (int i = 0; i < 8 && i < nnodes; i++) {
				nindices[i] = parse() - 1;
			}
			data += lineEndSize;

			auto readNextNodes = [&] () {
				for (int i = 0; i < nnodes - 8; i++) {
					nindices[i + 8] = parse() - 1;
				}
				data += lineEndSize;
			};

			switch (et[ansystype.back()].etype()) {
			case ET::ETYPE::D2SOLID_4NODES:
				if (nindices[2] == nindices[3]) { // triangle3
					esize.push_back(3);
					type.back() = (esint)Element::CODE::TRIANGLE3;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
				} else { // square4
					esize.push_back(4);
					type.back() = (esint)Element::CODE::SQUARE4;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 4);
				}
				break;
			case ET::ETYPE::D2SOLID_6NODES:
				esize.push_back(6);
				type.back() = (esint)Element::CODE::TRIANGLE6;
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 6);
				break;
			case ET::ETYPE::D2SOLID_8NODES:
				if (nindices[2] == nindices[3]) { // triangle6
					esize.push_back(6);
					type.back() = (esint)Element::CODE::TRIANGLE6;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
					nodes.insert(nodes.end(), nindices.begin() + 4, nindices.begin() + 6);
					nodes.push_back(nindices[7]);
				} else { // square8
					esize.push_back(8);
					type.back() = (esint)Element::CODE::SQUARE8;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 8);
				}
				break;
			case ET::ETYPE::D3SOLID_4NODES:
				esize.push_back(4);
				type.back() = (esint)Element::CODE::TETRA4;
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 4);
				break;
			case ET::ETYPE::D3SOLID_8NODES:
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra4
						esize.push_back(4);
						type.back() = (esint)Element::CODE::TETRA4;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
						nodes.push_back(nindices[4]);
					} else { // prisma6
						esize.push_back(6);
						type.back() = (esint)Element::CODE::PRISMA6;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
						nodes.insert(nodes.end(), nindices.begin() + 4, nindices.begin() + 7);
					}
				} else {
					if (nindices[4] == nindices[5]) { // pyramid5
						esize.push_back(5);
						type.back() = (esint)Element::CODE::PYRAMID5;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 5);
					} else { // hexa8
						esize.push_back(8);
						type.back() = (esint)Element::CODE::HEXA8;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 8);
					}
				}
				break;
			case ET::ETYPE::D3SOLID_10NODES:
				readNextNodes();
				esize.push_back(10);
				type.back() = (esint)Element::CODE::TETRA10;
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 10);
				break;
			case ET::ETYPE::D3SOLID_20NODES:
				readNextNodes();
				if (nindices[2] == nindices[3]) {
					if (nindices[4] == nindices[5]) { // tetra10
						esize.push_back(10);
						type.back() = (esint)Element::CODE::TETRA10;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
						nodes.push_back(nindices[4]);

						nodes.insert(nodes.end(), nindices.begin() + 8, nindices.begin() + 10);
						nodes.push_back(nindices[11]);
						nodes.insert(nodes.end(), nindices.begin() + 16, nindices.begin() + 19);
					} else { // prisma15
						esize.push_back(15);
						type.back() = (esint)Element::CODE::PRISMA15;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
						nodes.insert(nodes.end(), nindices.begin() + 4, nindices.begin() + 7);

						nodes.insert(nodes.end(), nindices.begin() + 8, nindices.begin() + 10);
						nodes.push_back(nindices[11]);
						nodes.insert(nodes.end(), nindices.begin() + 12, nindices.begin() + 14);
						nodes.push_back(nindices[15]);
						nodes.insert(nodes.end(), nindices.begin() + 16, nindices.begin() + 19);
					}
				} else {
					if (nindices[4] == nindices[5]) { // pyramid13
						esize.push_back(13);
						type.back() = (esint)Element::CODE::PYRAMID13;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 5);
						nodes.insert(nodes.end(), nindices.begin() + 8, nindices.begin() + 12);
						nodes.insert(nodes.end(), nindices.begin() + 16, nindices.begin() + 20);
					} else { // hexa20
						esize.push_back(20);
						type.back() = (esint)Element::CODE::HEXA20;
						nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 20);
					}
				}
				break;
			case ET::ETYPE::TARGET:
			case ET::ETYPE::CONTACT:
			case ET::ETYPE::PRETS:
			case ET::ETYPE::SKIP:
				esize.push_back(nnodes);
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + nnodes);
				type.back() = (esint)Element::CODE::NOT_SUPPORTED;
				break;
			default:
				eslog::error("ESPRESO Workbench parser: not implemented parsing of etype[%d] = %d\n", type.back() + 1, et[type.back()].type);
			}
		}
		tIDs[t].swap(IDs);
		tesize[t].swap(esize);
		tnodes[t].swap(nodes);
		ttype[t].swap(type);
		tet[t].swap(ansystype);
		tbody[t].swap(body);
		tmat[t].swap(mat);
	}

	for (size_t t = 0; t < threads; t++) {
		mesh.eIDs.insert(mesh.eIDs.end(), tIDs[t].begin(), tIDs[t].end());
		mesh.esize.insert(mesh.esize.end(), tesize[t].begin(), tesize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), tnodes[t].begin(), tnodes[t].end());
		mesh.etype.insert(mesh.etype.end(), ttype[t].begin(), ttype[t].end());
		mesh.et.insert(mesh.et.end(), tet[t].begin(), tet[t].end());
		mesh.body.insert(mesh.body.end(), tbody[t].begin(), tbody[t].end());
		mesh.material.insert(mesh.material.end(), tmat[t].begin(), tmat[t].end());
	}
	return true;
}

bool EBlock::boundary(const std::vector<ET> &et, AnsysCDBData &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, last - first);

	std::vector<std::vector<esint> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tet(threads), tbody(threads), tmat(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<esint> nindices(20);

		std::vector<esint> esize, nodes, IDs;
		std::vector<int> type, ansystype, body, mat;

		auto data = first + tdistribution[t];
		if (tdistribution[t] && *(data - 1) != '\n') {
			while (data < first + tdistribution[t + 1] && *data++ != '\n'); // start at new line
		}

		auto skip = [&] (const char* &data) {
			data += valueLength;
		};

		auto parse = [&] (const char* &data) {
			memcpy(value.data(), data, valueLength);
			skip(data);
			return atol(value.data());
		};

		while (data < first + tdistribution[t + 1]) {
			IDs.push_back(parse(data) - 1); // element ID
			body.push_back(0);
			type.push_back(0);
			ansystype.push_back(-1);
			skip(data); // section ID
			skip(data); // real constant
			mat.push_back(parse(data) - 1); // material
			skip(data); // element coordinate system

			int enodes = 0;
			while (*data != '\r' && *data != '\n') {
				nindices[enodes++] = parse(data) - 1;
			}
			data += lineEndSize;

			if (enodes == 2) { // line2
				esize.push_back(2);
				type.back() = (esint)Element::CODE::LINE2;
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 2);
			}

			if (enodes == 3) { // line3
				esize.push_back(3);
				type.back() = (esint)Element::CODE::LINE3;
				nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
			}

			if (enodes == 4) {
				if (nindices[2] == nindices[3]) { // triangle3
					esize.push_back(3);
					type.back() = (esint)Element::CODE::TRIANGLE3;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
				} else { // square4
					esize.push_back(4);
					type.back() = (esint)Element::CODE::SQUARE4;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 4);
				}
			}
			if (enodes == 8) {
				if (nindices[2] == nindices[3]) { // triangle6
					esize.push_back(6);
					type.back() = (esint)Element::CODE::TRIANGLE6;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 3);
					nodes.insert(nodes.end(), nindices.begin() + 4, nindices.begin() + 6);
					nodes.push_back(nindices[7]);
				} else { // square8
					esize.push_back(8);
					type.back() = (esint)Element::CODE::SQUARE8;
					nodes.insert(nodes.end(), nindices.begin(), nindices.begin() + 8);
				}
			}
		}

		tIDs[t].swap(IDs);
		tesize[t].swap(esize);
		tnodes[t].swap(nodes);
		ttype[t].swap(type);
		tet[t].swap(ansystype);
		tbody[t].swap(body);
		tmat[t].swap(mat);
	}

	for (size_t t = 0; t < threads; t++) {
		mesh.eIDs.insert(mesh.eIDs.end(), tIDs[t].begin(), tIDs[t].end());
		mesh.esize.insert(mesh.esize.end(), tesize[t].begin(), tesize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), tnodes[t].begin(), tnodes[t].end());
		mesh.etype.insert(mesh.etype.end(), ttype[t].begin(), ttype[t].end());
		mesh.et.insert(mesh.et.end(), tet[t].begin(), tet[t].end());
		mesh.body.insert(mesh.body.end(), tbody[t].begin(), tbody[t].end());
		mesh.material.insert(mesh.material.end(), tmat[t].begin(), tmat[t].end());
	}

	return true;
}

