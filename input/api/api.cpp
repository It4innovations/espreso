
#include "api.h"

using namespace esinput;


void API::elements(std::vector<mesh::Element*> &elements)
{
	// TODO: elements with various DOFS
	int DOFS = 3;

	elements.reserve(eIndices.size());
	eslocal indices[20], params[6];

	for (eslocal e = 0; e < eIndices.size(); e++) {
		for (eslocal i = 0; i < eIndices.size(); i += DOFS) {
			indices[i / DOFS] = eIndices[e][i] / DOFS;
		}
		switch(eIndices[e].size() / DOFS) {
		case Tetrahedron4NodesCount:
			elements.push_back(new mesh::Tetrahedron4(indices, eIndices[e].size() / DOFS, params));
			break;
		case Tetrahedron10NodesCount:
			elements.push_back(new mesh::Tetrahedron10(indices, eIndices[e].size() / DOFS, params));
			break;
		case Pyramid5NodesCount:
			elements.push_back(new mesh::Pyramid5(indices, eIndices[e].size() / DOFS, params));
			break;
		case Pyramid13NodesCount:
			elements.push_back(new mesh::Pyramid13(indices, eIndices[e].size() / DOFS, params));
			break;
		case Prisma6NodesCount:
			elements.push_back(new mesh::Prisma6(indices, eIndices[e].size() / DOFS, params));
			break;
		case Prisma15NodesCount:
			elements.push_back(new mesh::Prisma15(indices, eIndices[e].size() / DOFS, params));
			break;
		case Hexahedron8NodesCount:
			elements.push_back(new mesh::Hexahedron8(indices, eIndices[e].size() / DOFS, params));
			break;
		case Hexahedron20NodesCount:
			elements.push_back(new mesh::Hexahedron20(indices, eIndices[e].size() / DOFS, params));
			break;
		default:
			ESLOG(eslog::ERROR) << "Unknown element with " << eIndices[e].size() / DOFS << " indices.";
		}
	}
}



