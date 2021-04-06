
#include "ssection.h"

#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"

#include "esinfo/envinfo.h"

using namespace espreso;

size_t SSection::size = 14;
const char* SSection::upper     = "*SOLID SECTION";
const char* SSection::lower     = "*solid section";
const char* SSection::sentence  = "*Solid Section";
SSection::SSection()
{
	memset(Elset, '\0', MAX_NAME_SIZE);
	memset(Material, '\0', MAX_NAME_SIZE);
}

SSection& SSection::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	for (size_t i=0; i < command.size();i++)
	{
		std::vector<std::string> ssection_line = Parser::split(command[i], "=", false);
		ssection_line[0] = Parser::strip(ssection_line[0]);

		if ( StringCompare::caseInsensitiveEq(ssection_line[0], "ELSET")||
			 StringCompare::caseInsensitiveEq(ssection_line[0], "Elset")||
			 StringCompare::caseInsensitiveEq(ssection_line[0], "elset")) {
			ssection_line[1] = Parser::strip(ssection_line[1]);
			memcpy(Elset,ssection_line[1].data(),ssection_line[1].size());
		}

		if ( StringCompare::caseInsensitiveEq(ssection_line[0], "MATERIAL")||
			 StringCompare::caseInsensitiveEq(ssection_line[0], "material")||
			 StringCompare::caseInsensitiveEq(ssection_line[0], "Material")) {
			ssection_line[1] = Parser::strip(ssection_line[1]);
			memcpy(Material,ssection_line[1].data(),ssection_line[1].size());
		}
	}
	return *this;
}
/*
bool Eset::readData(std::vector<esint> &indices)
{

	size_t threads = environment->OMP_NUM_THREADS;

	const char *first = getFirst(), *last = getLast();
	esint size = (last - first) / lineSize;

	std::vector<size_t> tdistribution = tarray<esint>::distribute(threads, size);
	std::vector<std::vector<esint> > tindices(threads);


	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {

		for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1];) {

			std::string eleset_line = Parser::getLine(data);
			std::vector<std::string> element_data = Parser::split(eleset_line, ",");

			for (int ii=0; ii<element_data.size();ii++){
				if (atol(element_data[ii].data()) != 0) {
				tindices[t].push_back(atol(element_data[ii].data()));}
			}


			if (*(data ) != '\n') {
				while (*data++ != '\n'); // start at new line
			}
		}

	}

	for (size_t t = 0; t < threads; t++) {
		indices.insert(indices.end(), tindices[t].begin(), tindices[t].end());
	}

	std::sort(indices.begin(), indices.end());
	return true;

}
	size_t offset = coordinates.size();
	size_t adjustment = 0;
	nIDs.resize(offset + size);
	coordinates.resize(offset + size);
	std::vector<size_t> tdistribution = tarray<esint>::distribute(threads, size);

	const char *adjustment_data;
    if (first + lineSize * tdistribution[threads] < last) {
    	adjustment_data = first + lineSize * tdistribution[threads];
    	while (adjustment_data++ != last);
    	 adjustment++;
    }
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		const char *data_end =  first + lineSize * tdistribution[t + 1];
		if (t==threads) {data_end+=adjustment;}
		esint i = 0;
		double x, y, z;
		for (auto data = first + lineSize * tdistribution[t]; data < data_end; ++i) {

			std::string node_line = Parser::getLine(data);
			std::vector<std::string> node_data = Parser::split(node_line, ",");

			nIDs[offset + tdistribution[t] + i] = atol(node_data[0].data()) - 1;

			x=atof(node_data[1].data());
			y=atof(node_data[1].data());
			z=atof(node_data[1].data());

			//element += elementSize;
			if (*(data ) != '\n') {
									while (*data++ != '\n'); // start at new line
							}
			coordinates[offset + tdistribution[t] + i] = Point(x, y, z);
		}
	}
	return true;

*/



