#include "element.h"
#include <iomanip>
#include <iostream>


std::ostream& operator<<(std::ostream& os, const Element &e)
{
	for (size_t i = 0; i < e.size(); i++) {
		os << e.node(i) << " ";
	}
	return os;
}


