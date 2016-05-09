
#ifndef MESH_STRUCTURES_MATERIAL_H_
#define MESH_STRUCTURES_MATERIAL_H_

struct Material {

	Material():
		density(7850), youngModulus(2.1e11), poissonRatio(0.3),
		termalExpansion(1), termalCapacity(1), termalConduction(1) {};

	double density;

	double youngModulus;
	double poissonRatio;

	double termalExpansion;
	double termalCapacity;
	double termalConduction;
};



#endif /* MESH_STRUCTURES_MATERIAL_H_ */
