
#include "assembler.h"

namespace assembler {

template<>
Assembler<FEM>::Assembler(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh.rank(), mesh.size()), _verbose(true) { }

template<>
Assembler<BEM>::Assembler(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh), _verbose(true)
{
	_surface.computeFixPoints(8);
}

template<>
Assembler<ELMER>::Assembler(const mesh::Mesh &mesh): _mesh(mesh), _surface(mesh.rank(), mesh.size()), _verbose(true) { }

template<>
size_t Assembler<FEM>::subdomains()
{
	return this->_mesh.parts();
}

template<>
size_t Assembler<BEM>::subdomains()
{
	return this->_surface.parts();
}

template<>
size_t Assembler<ELMER>::subdomains()
{
	return this->_surface.parts();
}

template<>
const mesh::Mesh& Assembler<FEM>::mesh() const
{
	return this->_mesh;
}

template<>
const mesh::Mesh& Assembler<BEM>::mesh() const
{
	return this->_surface;
}

template<>
const mesh::Mesh& Assembler<ELMER>::mesh() const
{
	return this->_mesh;
}

}

