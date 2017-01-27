
#include "parseerror.h"

#include "../../../basis/logging/logging.h"

using namespace espreso;

void input::ParseError::print() {
	ESINFO(ERROR) << "ParseError: " << message;
}


