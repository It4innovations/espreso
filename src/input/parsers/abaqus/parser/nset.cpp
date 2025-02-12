
#include "nset.h"

#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

#include <algorithm>

using namespace espreso;

size_t Nset::size = 5;
const char* Nset::upper     = "*NSET";
const char* Nset::lower     = "*nset";
const char* Nset::sentence  = "*Nset";
Nset::Nset()
: NUMFIELD(-1), Solkey(0), NDMAX(-1), NDSEL(-1),
  lineSize(0), lineEndSize(-1),
  indexSize(-1), indexLength(-1), valueSize(-1), valueLength(-1)
{

}

Nset& Nset::parse(const char* begin)
{
    std::string commandLine = Parser::getLine(begin);

    std::vector<std::string> command = Parser::split(commandLine, ",", false);
    for (size_t i = 0; i < command.size(); i++)
    {
        std::vector<std::string> nset_name = Parser::split(command[i], "=", false);
        nset_name[0] = Parser::strip(nset_name[0]);
        if (StringCompare::caseInsensitiveEq(nset_name[0], "NSET")||
            StringCompare::caseInsensitiveEq(nset_name[0], "Nset")||
            StringCompare::caseInsensitiveEq(nset_name[0], "nset")) {

            nset_name[1] = Parser::strip(nset_name[1]);
            memcpy(NAME,nset_name[1].data(),nset_name[1].size());
        }
    }

    lineEndSize = (*(commandLine.end() - 2) == '\r') ? 2 : 1;

    const char* linebegin = begin ;
    if (*(linebegin ) != '\n') {while (*linebegin++ != '\n');} // start at new line

    if (*(linebegin ) != '\n') {
        while (*linebegin++ != '\n') {
            lineSize +=1; // start at new line
        }
    }
    std::string nodeLine = Parser::getLine(linebegin);
    std::vector<std::string> nodeData = Parser::split(nodeLine, ",");

    indexSize = 1;
    indexLength = nodeData[0].size();
    valueSize = nodeData.size()-1;
    valueLength = nodeData[1].size();

    lineSize = lineSize+1;
    AbaqusParser::fillIndices(begin, begin + commandLine.size() );
    return *this;
}

bool Nset::readData(std::vector<esint> &indices)
{
    size_t threads = info::env::OMP_NUM_THREADS;

    const char *first = getFirst(), *last = getLast();
    esint size = (last - first) / lineSize;

    std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, size);
    std::vector<std::vector<esint> > tindices(threads);


    #pragma omp parallel for
    for (size_t t = 0; t < threads; t++) {
        for (auto data = first + lineSize * tdistribution[t]; data < first + lineSize * tdistribution[t + 1];) {

            std::string nset_line = Parser::getLine(data);
            std::vector<std::string> element_data = Parser::split(nset_line, ",");

            for (size_t ii = 0; ii < element_data.size(); ii++){
                if (atol(element_data[ii].data()) != 0) {
                    tindices[t].push_back(atol(element_data[ii].data()));
                }
            }


            if (*(data ) != '\n') {
                while (*data++ != '\n'); // start at new line
            }
        }//auto
        if (lRank == info::mpi::rank && t == threads - 1) {
            auto data = first + lineSize * tdistribution[t + 1];
            std::string eleset_line = Parser::getLine(data);
            std::vector<std::string> element_data = Parser::split(eleset_line, ",");

            for (size_t ii = 0; ii < element_data.size(); ii++){
                if (atol(element_data[ii].data()) != 0) {
                tindices[t].push_back(atol(element_data[ii].data()));}
            }
        }
    }

    for (size_t t = 0; t < threads; t++) {
        indices.insert(indices.end(), tindices[t].begin(), tindices[t].end());
    }

    std::sort(indices.begin(), indices.end());
    return true;

}

