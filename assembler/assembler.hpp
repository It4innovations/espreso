
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
	//return this->_surface.parts();
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
	return 1;
	//return this->_surface.parts();
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
	return 1;
	//return this->_surface.parts();
}
//
//template<>
//const mesh::Mesh& Assembler<FEM>::mesh() const
//{
//	return this->_mesh;
//}
//
//template<>
//const mesh::Mesh& Assembler<BEM>::mesh() const
//{
//	return this->_surface;
//}
//
//template<>
//const mesh::Mesh& Assembler<API>::mesh() const
//{
//	return this->_mesh;
//}

}

