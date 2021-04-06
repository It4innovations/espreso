
#include "blockend.h"

using namespace espreso;

BlockEnd::BlockEnd()
{

}

BlockEnd& BlockEnd::parse(const char* begin)
{
	while (*(--begin) != '\n');
	WorkbenchParser::fillIndices(begin + 1, begin + 1);
	return *this;
}

