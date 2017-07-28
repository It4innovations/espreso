
#include "../../basis/utilities/utils.h"

#include "../instance.h"

#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/elements/plane/planeelement.h"
#include "boundarybased3d.h"

using namespace espreso;

BoundaryBased3D::BoundaryBased3D()
: Physics("", NULL, NULL) // skipped because Physics is inherited virtually
{

}

void BoundaryBased3D::extractBoundaryNodes()
{
	_mesh->computeFacesOnDomainsSurface();
	_boundaryIndices.resize(_mesh->parts());

	#pragma omp parallel for
	for (size_t p = 0; p < _mesh->parts(); p++) {
		for (size_t i = 0; i < _mesh->faces().size(); i++) {
			if (std::find(_mesh->faces()[i]->domains().begin(), _mesh->faces()[i]->domains().end(), p) != _mesh->faces()[i]->domains().end()) {
				for (size_t n = 0; n < _mesh->faces()[i]->nodes(); n++) {
					_boundaryIndices[p].push_back(_mesh->faces()[i]->node(n));
				}
			}
		}
		std::sort(_boundaryIndices[p].begin(), _boundaryIndices[p].end());
		Esutils::removeDuplicity(_boundaryIndices[p]);
	}

	for (size_t p = 0; p < _mesh->parts(); p++) {
		for (size_t i = 0; i < _boundaryIndices[p].size(); i++) {
			Element *node = _mesh->nodes()[_boundaryIndices[p][i]];
			size_t d = std::lower_bound(node->domains().begin(), node->domains().end(), p) - node->domains().begin();
			for (size_t dof = 0; dof < pointDOFs().size(); dof++) {
				node->DOFsIndices()[pointDOFs().size() * d + dof] = pointDOFs().size() * i + dof;
			}
		}
		_instance->domainDOFCount[p] = pointDOFs().size() * _boundaryIndices[p].size();
	}
}


void BoundaryBased3D::boundaryTriangularization(std::vector<eslocal> &elements, std::vector<double> &coordinates, size_t domain)
{
	bool swap;

	std::vector<std::vector<eslocal> > triangles;

	for (size_t i = 0; i < _mesh->faces().size(); i++) {
		if (std::find(_mesh->faces()[i]->domains().begin(), _mesh->faces()[i]->domains().end(), domain) != _mesh->faces()[i]->domains().end()) {

			for (size_t p = 0; p < _mesh->faces()[i]->parentElements().size(); p++) {
				if (_mesh->faces()[i]->parentElements()[p]->domains().front() == domain) {
					swap = _mesh->faces()[i]->parentElements()[p]->isFaceSwapped(_mesh->faces()[i]);
					break;
				}
			}

			triangles = dynamic_cast<PlaneElement*>(_mesh->faces()[i])->triangularize();

			for (size_t t = 0; t < triangles.size(); t++) {
				for (size_t n = 0; n < 3; n++) {
					elements.push_back(std::lower_bound(_boundaryIndices[domain].begin(), _boundaryIndices[domain].end(), triangles[t][swap ? 2 - n : n]) - _boundaryIndices[domain].begin());
				}
			}
		}
	}

	for (size_t i = 0; i < _boundaryIndices[domain].size(); i++) {
		coordinates.push_back(_mesh->coordinates()[_boundaryIndices[domain][i]].x);
		coordinates.push_back(_mesh->coordinates()[_boundaryIndices[domain][i]].y);
		coordinates.push_back(_mesh->coordinates()[_boundaryIndices[domain][i]].z);
	}
}



