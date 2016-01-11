
#include "assembler.h"

namespace assembler {

template<>
size_t Assembler<FEM>::subdomains()
{
	return this->_input.mesh.parts();
}

template<>
size_t Assembler<BEM>::subdomains()
{
	return this->_input.surface.parts();
}

template<>
size_t Assembler<API>::subdomains()
{
	return 1;
}

template<>
size_t Assembler<API2>::subdomains()
{
	return 1;
}

}

