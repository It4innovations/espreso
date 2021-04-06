
#include "nblock.h"

#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"

#include "esinfo/envinfo.h"

using namespace espreso;

size_t NList::size = 5;
const char* NList::upper     = "*NODE";
const char* NList::lower     = "*node";
const char* NList::sentence  = "*Node";
NList::NList()
: NUMFIELD(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(0), lineEndSize(-1),
  indexSize(-1), indexLength(-1), valueSize(-1), valueLength(-1)
{

}

NList& NList::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);

	lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

	const char* linebegin = begin ;
	if (*(linebegin ) != '\n') {
		while (*linebegin++ != '\n'); // start at new line
	}

	if (*(linebegin ) != '\n') {
		while (*linebegin++ != '\n') {
			lineSize +=1;}; // start at new line
	}
	std::string nodeLine = Parser::getLine(linebegin);
	std::vector<std::string> nodeData = Parser::split(nodeLine, ",");

	indexSize = 1;
	indexLength = nodeData[0].size();
	valueSize = nodeData.size()-1;
	valueLength = nodeData[1].size();

	lineSize = lineSize+1;//lineEndSize;
	AbaqusParser::fillIndices(begin, begin + commandLine.size() );
	return *this;
}

bool NList::readData(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	return index_solid_line_x_y_z(nIDs, coordinates, scaleFactor);
}

bool NList::index_x_y_z(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	size_t size = (last - first) / lineSize;

	size_t offset = coordinates.size();
	nIDs.resize(offset + size);
	coordinates.resize(offset + size);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<char> value(valueLength + 1);
		std::vector<char> index(indexLength + 1);
		esint i = 0;
		double x, y, z;
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1]; ++i) {
			memcpy(index.data(), data, indexLength);
			data += indexLength;
			nIDs[offset + tdistribution[t] + i] = atol(index.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			x = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			y = scaleFactor * atof(value.data());

			memcpy(value.data(), data, valueLength);
			data += valueLength;
			z = scaleFactor * atof(value.data());

			data = data + lineEndSize +1;

			coordinates[offset + tdistribution[t] + i] = Point(x, y, z);
		}
	}
	return true;
}

bool NList::index_solid_line_x_y_z(std::vector<esint> &nIDs, std::vector<Point> &coordinates, double scaleFactor)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	esint size = (last - first) / lineSize;

	size_t offset = coordinates.size();
	nIDs.resize(offset + size);
	coordinates.resize(offset + size);
	std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint i = 0;
		double x, y, z;
		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1]; ++i) {

			std::string node_line = Parser::getLine(data);
			std::vector<std::string> node_data = Parser::split(node_line, ",");

			nIDs[offset + tdistribution[t] + i] = atol(node_data[0].data());

			x=atof(node_data[1].data());
			y=atof(node_data[2].data());
			z=atof(node_data[3].data());

			//element += elementSize;
			if (*(data ) != '\n') {
									while (*data++ != '\n'); // start at new line
							}
			coordinates[offset + tdistribution[t] + i] = Point(x, y, z);
		}
	}
	return true;
}
