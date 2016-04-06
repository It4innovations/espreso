
#include "assembler.h"

namespace espreso {

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
	return this->_input.mesh->parts();
}

template<>
size_t Assembler<FEM>::DOFs()
{
	return this->_input.mesh.DOFs();
}

template<>
size_t Assembler<BEM>::DOFs()
{
	return this->_input.surface.DOFs();
}

template<>
size_t Assembler<API>::DOFs()
{
	return this->_input.mesh->DOFs();
}

}

