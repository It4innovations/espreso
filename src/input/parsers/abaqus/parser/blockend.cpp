
/* #include "blockend.h"

using namespace espreso;

size_t BlockEnd::nSize = 2;
const char* BlockEnd::nUpper = "N,";
const char* BlockEnd::nLower = "n,";

size_t BlockEnd::unixSize = 3;
const char* BlockEnd::unixEnd = "-1\n";

size_t BlockEnd::winSize = 4;
const char* BlockEnd::winEnd = "-1\r\n";

BlockEnd::BlockEnd()
{

}

BlockEnd& BlockEnd::parse(const char* begin)
{
	while (*(--begin) != '\n');
	AbaqusParser::fillIndices(begin + 1, begin + 1);
	return *this;
}
*/

#include "blockend.h"

using namespace espreso;

size_t BlockFinish::nSize  = 1;
size_t BlockFinish::ncSize = 2;

const char* BlockFinish::asterik = "*";
const char* BlockFinish::comment = "**";


BlockFinish::BlockFinish()
{

}

BlockFinish& BlockFinish::parse(const char* begin)
{
	while (*(--begin) != '\n');
	AbaqusParser::fillIndices(begin + 1, begin + 1);
	return *this;
}

