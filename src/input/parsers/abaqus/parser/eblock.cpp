
#include "eblock.h"

#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "mesh/element.h"
#include "input/parsers/abaqus/abaqus.h"

using namespace espreso;

size_t EList::size = 8;
const char* EList::upper = "*ELEMENT";
const char* EList::lower = "*element";
const char* EList::sentence = "*Element";

EList::EList()
: NUM_NODES(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(0), elementSize(0), lineEndSize(-1),
  indexSize(-1), indexLength(-1), valueSize(-1), valueLength(-1),
  NUM_NODES_LINE1(0),NUM_NODE_LINE2(0)
{
	memset(TYPE, '\0', MAX_NAME_SIZE);
}

EList& EList::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);
	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	const char* linebegin = begin ;
	if (*(linebegin ) != '\n') {
		while (*linebegin++ != '\n'); // start at new line
	}

	std::string elementLine = Parser::getLine(linebegin);
	std::vector<std::string> elementData = Parser::split(elementLine, ",");

	if (*(linebegin ) != '\n') {
		while (*linebegin++ != '\n') {
			lineSize +=1;  // start at new line
		}
	}

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	std::vector<std::string> ele_type = Parser::split(command[1], "=", false);
	ele_type[1] = Parser::strip(ele_type[1]);
	memcpy(TYPE,ele_type[1].data(),ele_type[1].size());

	//std::istringstream iss(TYPE);
	if (strcmp(TYPE, "C3D20R") == 0 ){
		NUM_NODES = 20;
	}

	indexSize = 1;
	indexLength = elementData[0].size()+1;
	valueSize = elementData.size()-1;
	valueLength = elementData[1].size()+1;

	if ( atol(elementData.back().data()) == 0) {
		valueSize = valueSize -1;
	}

	NUM_NODES_LINE1 = (lineSize - indexLength) / valueLength;
	if (NUM_NODES_LINE1 < NUM_NODES){
		NUM_NODE_LINE2 = NUM_NODES - NUM_NODES_LINE1;
	}

	elementSize = lineSize + 1;
	lineSize = lineSize + 1;
	AbaqusParser::fillIndices(begin, begin + commandLine.size());
	return *this;
}


void EList::fixOffsets(std::vector<size_t> &dataOffsets)
{
	esint nodes;
	if (fRank == info::mpi::rank) {
		nodes = NUM_NODES;
	}
	Communication::broadcast(&nodes, sizeof(esint), MPI_BYTE, fRank);

	int size1 = lineSize; // first line
	int size2 = 0;
	if (NUM_NODE_LINE2 != 0) {
		size2 = NUM_NODE_LINE2 * valueLength -1 + lineEndSize + indexLength;
	}

	if (fRank != lRank && NUM_NODE_LINE2 != 0) {
		for (int rank = fRank + 1; rank <= lRank; rank++) {
			if ((dataOffsets[rank] - first) % (size1 + size2) != 0) {
				dataOffsets[rank] += size2;
				if (rank - 1 == info::mpi::rank) {
					end += size2;
				}
				if (rank == info::mpi::rank) {
					begin += size2;
					offset += size2;
				}
			}
		}
	}

	if (NUM_NODE_LINE2 == 0) {
		//valueSize = valueSize - 8 + nodes;
		//lineSize = valueSize * valueLength + lineEndSize;
		Communication::broadcast(&valueSize, sizeof(esint), MPI_BYTE, fRank);
		Communication::broadcast(&lineSize, sizeof(esint), MPI_BYTE, fRank);
		size1 = elementSize = lineSize;
	}
	if (NDSEL == -1) {
		NDSEL = (first - last) / (size1 + size2);
	}
	elementSize = size1 + size2;
}

bool EList::readData(MeshBuilder &mesh)
{
		return solid( mesh);
}

bool EList::solid(MeshBuilder &mesh)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();

	std::vector<std::vector<esint> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tet(threads), tbody(threads), tmat(threads);
	size_t size = (last - first) / elementSize;
//	if (  element_dict[TYPE].compare("D1LINE_2NODES")==0     ) {
	if (strcmp(TYPE, "B32") == 0 ) {
		std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> esize, nodes, IDs;
			std::vector<int> type, ansystype, body, mat;
			for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
				esize.push_back(3);
				type.push_back(0);
				type.back() = (esint)Element::CODE::LINE3;
				std::string element_line = Parser::getLine(element);
				std::vector<std::string> element_data = Parser::split(element_line, ",");
				IDs.push_back(atol(element_data[0].data()));
				for (size_t ii = 1; ii < element_data.size(); ii++){
					nodes.push_back(atol(element_data[ii].data()));
				}
				//element += elementSize;
				if (*(element ) != '\n') {  while (*element++ != '\n');  }  //element += lineEndSize;
			}
			tIDs[t].swap(IDs);
			tesize[t].swap(esize);
			tnodes[t].swap(nodes);
			ttype[t].swap(type);
		}
	}
	
	if (strcmp(TYPE, "C3D8R") == 0 ) {
//	else if (  element_dict[TYPE].compare("D3SOLID_8NODES")==0     ) {
		std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> esize, nodes, IDs;
			std::vector<int> type, ansystype, body, mat;
			for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
				esize.push_back(8);
				type.push_back(0);
				type.back() = (esint)Element::CODE::HEXA8;
				std::string element_line = Parser::getLine(element);
				std::vector<std::string> element_data = Parser::split(element_line, ",");
				IDs.push_back(atol(element_data[0].data()));
				for (size_t ii  =1; ii < element_data.size(); ii++){
					nodes.push_back(atol(element_data[ii].data()));
				}
				if (*(element ) != '\n') {  while (*element++ != '\n');  } // start at new line

			}
			tIDs[t].swap(IDs);
			tesize[t].swap(esize);
			tnodes[t].swap(nodes);
			ttype[t].swap(type);
		}
	}

	if (strcmp(TYPE, "C3D10") == 0 ){
//	else if(  element_dict [TYPE].compare("D3SOLID_10NODES")==0     ) {
		std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> esize, nodes, IDs;
			std::vector<int> type, ansystype, body, mat;
			for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
				esize.push_back(10);
				type.push_back(0);
				type.back() = (esint)Element::CODE::TETRA10;
				std::string element_line = Parser::getLine(element);
				std::vector<std::string> element_data = Parser::split(element_line, ",");
				IDs.push_back(atol(element_data[0].data()));
				for (size_t ii = 1; ii < element_data.size(); ii++){
					nodes.push_back(atol(element_data[ii].data()));
				}
				if (*(element ) != '\n') {  while (*element++ != '\n');  } // start at new line

			}  //for auto
		tIDs[t].swap(IDs);
		tesize[t].swap(esize);
		tnodes[t].swap(nodes);
		ttype[t].swap(type);
		}
	}

	if (strcmp(TYPE, "C3D20R") == 0 ){
//	else if(  element_dict [TYPE].compare("D3SOLID_10NODES")==0     ) {
		std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			std::vector<esint> esize, nodes, IDs;
			std::vector<int> type, ansystype, body, mat;
			for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
				esize.push_back(20);
				type.push_back(0);
				type.back() = (esint)Element::CODE::HEXA20;
				std::string element_line = Parser::getLine(element);
				std::vector<std::string> element_data = Parser::split(element_line, ",");
				IDs.push_back(atol(element_data[0].data()));
				for (size_t ii = 1; ii < element_data.size() - 1; ii++){
					nodes.push_back(atol(element_data[ii].data()));
				}
				if (*(element ) != '\n') {  while (*element++ != '\n');  } // start at new line

				element_line = Parser::getLine(element);
				element_data = Parser::split(element_line, ",");
				for (size_t ii = 0; ii < element_data.size(); ii++){
					nodes.push_back(atol(element_data[ii].data()));
				}
				if (*(element ) != '\n') {  while (*element++ != '\n');  } // start at new line


			}  //for auto
		tIDs[t].swap(IDs);
		tesize[t].swap(esize);
		tnodes[t].swap(nodes);
		ttype[t].swap(type);
		}
	}


/*
	size_t size = (last - first) / elementSize;

	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);


	std::vector<std::vector<esint> > tesize(threads), tnodes(threads), tIDs(threads);
	std::vector<std::vector<int> > ttype(threads), tet(threads), tbody(threads), tmat(threads);

	std::string element_line1 = Parser::getLine(first);

    #pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> esize, nodes, IDs;
			std::vector<int> type, body, mat;
		    if (strcmp(TYPE, "B32") == 0 ){
				for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
					esize.push_back(3);
					type.push_back(0);
					type.back() = (esint)Element::CODE::LINE3;
					std::string element_line = Parser::getLine(element);
					std::vector<std::string> element_data = Parser::split(element_line, ",");
					IDs.push_back(atol(element_data[0].data()));
                    for (int ii=1; ii<element_data.size();ii++){
                    	nodes.push_back(atol(element_data[ii].data()));
                    }
                    element += elementSize;
					element += lineEndSize;
				 }  //for auto
              } // if e type
		      if (strcmp(TYPE, "C3D8R") == 0 ){
				for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
					esize.push_back(8);
					type.push_back(0);
					type.back() = (esint)Element::CODE::HEXA8;
					std::string element_line = Parser::getLine(element);
					std::vector<std::string> element_data = Parser::split(element_line, ",");
					IDs.push_back(atol(element_data[0].data()));
					for (int ii=1; ii<element_data.size();ii++){
						nodes.push_back(atol(element_data[ii].data()));
					}
					//element += elementSize;
                                        if (*(element ) != '\n') {
                                                                while (*element++ != '\n'); // start at new line
                                                        }
					 }  //for auto
		      }// if e type C3D8R

		      if (strcmp(TYPE, "C3D10") == 0 ){
				for (auto element = first + elementSize * tdistribution[t]; element < first + elementSize * tdistribution[t + 1];) {
					esize.push_back(8);
					type.push_back(0);
					type.back() = (esint)Element::CODE::TETRA10;
					std::string element_line = Parser::getLine(element);
					std::vector<std::string> element_data = Parser::split(element_line, ",");
					IDs.push_back(atol(element_data[0].data()));
					for (int ii=1; ii<element_data.size();ii++){
						nodes.push_back(atol(element_data[ii].data()));
					}
					//element += elementSize;
                                        if (*(element ) != '\n') {
                                                                while (*element++ != '\n'); // start at new line
                                                        }
					}  //for auto
			  }// if e type C3D8R
        tIDs[t].swap(IDs);
		tesize[t].swap(esize);
		tnodes[t].swap(nodes);
		ttype[t].swap(type);
	} // for threads */
	for (size_t t = 0; t < threads; t++) {
		mesh.eIDs.insert(mesh.eIDs.end(), tIDs[t].begin(), tIDs[t].end());
		mesh.esize.insert(mesh.esize.end(), tesize[t].begin(), tesize[t].end());
		mesh.enodes.insert(mesh.enodes.end(), tnodes[t].begin(), tnodes[t].end());
		mesh.etype.insert(mesh.etype.end(), ttype[t].begin(), ttype[t].end());
	}
	mesh.body.resize(mesh.etype.size());
	mesh.material.resize(mesh.etype.size());
	return true;
}



