
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct ElasticMaterial: public ActionOperator {
	ElasticMaterial(int interval, const ParameterData &youngModulus, const ParameterData &poissonRatio, const ParameterData &shearModulus, ParameterData &elasticity)
	: youngModulus(youngModulus, interval),
	  poissonRatio(poissonRatio, interval),
	  shearModulus(shearModulus, interval),
	  elasticity(elasticity, interval)
	{

	}

	InputParameterIterator youngModulus, poissonRatio, shearModulus;
	OutputParameterIterator elasticity;

	void operator++()
	{
		++youngModulus; ++poissonRatio; ++shearModulus;
		++elasticity;
	}

	void move(int n)
	{
		youngModulus += n; poissonRatio += n; shearModulus += n;
		elasticity += n;
	}
};

template<size_t nodes, size_t gps>
struct ElasticityIsotropic2DPlaneStrain: public ElasticMaterial {
	using ElasticMaterial::ElasticMaterial;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double ex = youngModulus[gpindex];
			double mi = poissonRatio[gpindex];
			double k = ex * (1 - mi) / ((1 + mi) * (1 - 2 * mi));
			int ogp = 9 * gpindex;
			elasticity[ogp + 0] = k;                   elasticity[ogp + 1] = k * (mi / (1 - mi)); elasticity[ogp + 2] = 0;
			elasticity[ogp + 3] = k * (mi / (1 - mi)); elasticity[ogp + 4] = k;                   elasticity[ogp + 5] = 0;
			elasticity[ogp + 6] = 0;                   elasticity[ogp + 7] = 0;                   elasticity[ogp + 8] = k * ((1 - 2 * mi) / (2 * (1 - mi)));
		}
	}
};

template<size_t nodes, size_t gps>
struct ElasticityIsotropic2DPlaneStress: public ElasticMaterial {
	using ElasticMaterial::ElasticMaterial;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double ex = youngModulus[gpindex];
			double mi = poissonRatio[gpindex];
			double k = ex / (1 - mi * mi);
			int ogp = 9 * gpindex;
			elasticity[ogp + 0] = k;      elasticity[ogp + 1] = k * mi; elasticity[ogp + 2] = 0;
			elasticity[ogp + 3] = k * mi; elasticity[ogp + 4] = k;      elasticity[ogp + 5] = 0;
			elasticity[ogp + 6] = 0;      elasticity[ogp + 7] = 0;      elasticity[ogp + 8] = k * ((1 -  mi) / 2);
		}
	}
};

template<size_t nodes, size_t gps>
struct ElasticityIsotropic3D: public ElasticMaterial {
	using ElasticMaterial::ElasticMaterial;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double ex = youngModulus[gpindex];
			double mi = poissonRatio[gpindex];
			double ee = ex / ((1 + mi) * (1 - 2 * mi));
			int ogp = 36 * gpindex;
			elasticity[ogp +  0] = ee * (1.0 - mi); elasticity[ogp +  1] = ee * mi;         elasticity[ogp +  2] = ee * mi;         elasticity[ogp +  3] = 0;               elasticity[ogp +  4] = 0;               elasticity[ogp +  5] = 0;
			elasticity[ogp +  6] = ee * mi;         elasticity[ogp +  7] = ee * (1.0 - mi); elasticity[ogp +  8] = ee * mi;         elasticity[ogp +  9] = 0;               elasticity[ogp + 10] = 0;               elasticity[ogp + 11] = 0;
			elasticity[ogp + 12] = ee * mi;         elasticity[ogp + 13] = ee * mi;         elasticity[ogp + 14] = ee * (1.0 - mi); elasticity[ogp + 15] = 0;               elasticity[ogp + 16] = 0;               elasticity[ogp + 17] = 0;
			elasticity[ogp + 18] = 0;               elasticity[ogp + 19] = 0;               elasticity[ogp + 20] = 0;               elasticity[ogp + 21] = ee * (0.5 - mi); elasticity[ogp + 22] = 0;               elasticity[ogp + 23] = 0;
			elasticity[ogp + 24] = 0;               elasticity[ogp + 25] = 0;               elasticity[ogp + 26] = 0;               elasticity[ogp + 27] = 0;               elasticity[ogp + 28] = ee * (0.5 - mi); elasticity[ogp + 29] = 0;
			elasticity[ogp + 30] = 0;               elasticity[ogp + 31] = 0;               elasticity[ogp + 32] = 0;               elasticity[ogp + 33] = 0;               elasticity[ogp + 34] = 0;               elasticity[ogp + 35] = ee * (0.5 - mi);
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_H_ */
