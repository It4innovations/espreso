
#include "boundarycondition.h"
#include "mesh.h"

using namespace espreso;

SurfaceCondition::~SurfaceCondition()
{
	for (size_t i = 0; i < faces.size(); i++) {
		delete faces[i];
	}
}


void ElementCondition::fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values)
{

}

void ElementCondition::fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain)
{

}

void SurfaceCondition::fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values)
{

}

void SurfaceCondition::fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain)
{

}

void NodeCondition::fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values)
{
	if (_type != ConditionType::DIRICHLET) {
		return;
	}

	std::vector<eslocal> offsets;
	for (size_t i = 0; i < DOFs.size(); i++) {
		if (_prescribed[static_cast<int>(DOFs[i])]) {
			offsets.push_back(i);
		}
	}

	indices.reserve(indices.size() + offsets.size() * nodes.size());
	values.reserve(indices.size() + offsets.size() * nodes.size());
	for (size_t i = 0; i < nodes.size(); i++) {
		const Point3D &p = mesh.coordinates()[nodes[i]];
		for (size_t j = 0; j < offsets.size(); j++) {
			indices.push_back(DOFs.size() * nodes[i] + offsets[j]);
			values.push_back(_expression[static_cast<int>(DOFs[offsets[j]])].evaluate(p.x, p.y, p.z));
		}
	}
}

void NodeCondition::fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain)
{
	if (_type != ConditionType::FORCES) {
		return;
	}

	std::vector<eslocal> offsets;
	for (size_t i = 0; i < DOFs.size(); i++) {
		if (_prescribed[static_cast<int>(DOFs[i])]) {
			offsets.push_back(i);
		}
	}

	const std::vector<eslocal> &l2c = mesh.coordinates().localToCluster(subdomain);

	for (size_t i = 0; i < nodes.size(); i++) {
		auto it = std::lower_bound(l2c.begin(), l2c.end(), nodes[i]);
		if (it != l2c.end() && *it == nodes[i]) {
			const Point3D &p = mesh.coordinates()[nodes[i]];
			for (size_t j = 0; j < offsets.size(); j++) {
				f[it - l2c.begin() + offsets[j]] = _expression[static_cast<int>(DOFs[offsets[j]])].evaluate(p.x, p.y, p.z);
			}
		}
	}
}

