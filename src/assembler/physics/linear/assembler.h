
#ifndef SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_
#define SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_

#include "../assembler.h"

namespace espreso {

struct LinearPhysics: public Physics {

	virtual void assemble()
	{
		ESINFO(PROGRESS2) << "Assemble matrices K and RHS";
		const std::map<eslocal, double> &forces_x = _mesh.coordinates().property(FORCES_X).values();
		const std::map<eslocal, double> &forces_y = _mesh.coordinates().property(FORCES_Y).values();
		const std::map<eslocal, double> &forces_z = _mesh.coordinates().property(FORCES_Z).values();

		K.resize(_mesh.parts());
		f.resize(_mesh.parts());
		for (size_t p = 0; p < _mesh.parts(); p++) {
			composeSubdomain(p);

			const std::vector<eslocal> &l2g = _mesh.coordinates().localToCluster(p);
			for (eslocal i = 0; i < l2g.size(); i++) {
				size_t n = _mesh.subdomainBoundaries()[l2g[i]].size();
				if (forces_x.find(l2g[i]) != forces_x.end()) {
					f[p][DOFs * i + 0] = forces_x.at(l2g[i]) / n;
				}
				if (forces_y.find(l2g[i]) != forces_y.end()) {
					f[p][DOFs * i + 1] = forces_y.at(l2g[i]) / n;
				}
				if (forces_z.find(l2g[i]) != forces_z.end()) {
					f[p][DOFs * i + 2] = forces_z.at(l2g[i]) / n;
				}
			}

			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
		ESINFO(PROGRESS2);
	}

	virtual void save()
	{
		ESINFO(PROGRESS2) << "Save matrices K and RHS";
		for (size_t p = 0; p < _mesh.parts(); p++) {
			std::ofstream osK(Logging::prepareFile(p, "K").c_str());
			osK << K[p];
			osK.close();

			std::ofstream osF(Logging::prepareFile(p, "f").c_str());
			osF << f[p];
			osF.close();
		}
	}

	LinearPhysics(const Mesh &mesh, size_t DOFs): Physics(mesh, DOFs) {};
	virtual ~LinearPhysics() {};

	std::vector<SparseMatrix> K, T; // T will be deleted
	std::vector<std::vector<double> > f;

protected:
	virtual void composeSubdomain(size_t subdomain) =0;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_LINEAR_ASSEMBLER_H_ */
