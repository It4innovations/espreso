
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

template<>
size_t Assembler<FEM>::rank()
{
	return this->_input.mesh.rank();
}

template<>
size_t Assembler<BEM>::rank()
{
	return this->_input.surface.rank();
}

template<>
size_t Assembler<API>::rank()
{
	return esconfig::MPIrank;
}

template<>
size_t Assembler<API2>::rank()
{
	return esconfig::MPIrank;
}

template<>
size_t Assembler<FEM>::size()
{
	return this->_input.mesh.size();
}

template<>
size_t Assembler<BEM>::size()
{
	return this->_input.surface.size();
}

template<>
size_t Assembler<API>::size()
{
	return esconfig::MPIsize;
}

template<>
size_t Assembler<API2>::size()
{
	return esconfig::MPIsize;
}

}

