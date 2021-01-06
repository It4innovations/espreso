
#ifndef SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_H_
#define SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_H_

#include "kernel.h"
#include "config/holders/expression.h"

#include <functional>
#include <map>

namespace espreso {

template <typename TEBoundaries, typename TEData> class serializededata;
class VectorsDense;
class MatrixDense;
struct MaterialBaseConfiguration;
struct HeatTransferLoadStepConfiguration;
struct Builder;

class HeatTransferKernel: public Kernel
{
public:
	static const int GP_TRIANGLE3 = 6;
	static const int GP_TRIANGLE6 = 6;
	static const int GP_SQUARE4   = 4;
	static const int GP_SQUARE8   = 9;

	static const int GP_TETRA4    = 4;
	static const int GP_TETRA10   = 15;
	static const int GP_PYRAMID5  = 8;
	static const int GP_PYRAMID13 = 14;
	static const int GP_PRISMA6   = 9;
	static const int GP_PRISMA15  = 9;
	static const int GP_HEXA8     = 8;
	static const int GP_HEXA20    = 8;

	static NodeData *initialTemperature, *temperature;
	static ElementData *translationMotion, *phase, *latentHeat, *gradient, *flux;
	static void createParameters();

	HeatTransferKernel(HeatTransferKernel *previous, HeatTransferLoadStepConfiguration &configuration);
	~HeatTransferKernel();

	void nextSubstep();
	void solutionChanged();

	bool boundaryWithSettings(size_t rindex);

	void processElements(const Builder &builder, InstanceFiller &filler);
	void processBoundary(const Builder &builder, size_t rindex, InstanceFiller &filler);
	void processSolution();


	const HeatTransferLoadStepConfiguration &configuration;

	// node, elements parameters
	ElementParameterInfo _N, _dN, _w;
	ElementParameterInfo _time, _timestep, _ncoordinates, _ntemperature, _ninitTemperature;
	// parameters in Gauss points
	ElementParameterInfo _coordinates;
	ElementParameterInfo _initTemperature, _temperature;
	ElementParameterInfo _thickness;
	ElementParameterInfo _invJ, _detJ, _dND;
	ElementParameterInfo _csCartesian, _csSpherical, _csCylindrical, _dens, _hc, _conductivity[4];
	ElementParameterInfo _K, _isotropicK;
	ElementParameterInfo _translationMotion, _u, _advection;

	ElementParameterInfo _stiffness, _mass;

	BoundaryParameterInfo _btemperature, _heatflow, _heatflux;

protected:
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);
	void setMaterials(const std::map<std::string, std::string> &settings);
	void printMaterials(const std::map<std::string, std::string> &settings);
	void examineMaterialParameter(const std::string &material, const std::string &name, const ECFExpression &settings, ParameterInfo &info, int dimension);
	template<class TSecond>
	void examineElementParameter(const std::string &name, const std::map<std::string, TSecond> &settings, ParameterInfo &info, int dimension, std::function<const Evaluator*(const TSecond &expr)> getevaluator);
	void examineElementParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ParameterInfo &info)
	{
		examineElementParameter<ECFExpression>(name, settings, info, 0, [] (const ECFExpression &expr) { return expr.evaluator; });
	}
	void examineElementParameter(const std::string &name, const std::map<std::string, ECFExpressionVector> &settings, ParameterInfo &info, int dimension)
	{
		examineElementParameter<ECFExpressionVector>(name, settings, info, dimension, [&] (const ECFExpressionVector &expr) { return expr.data[dimension].evaluator; });
	}

	void examineBoundaryParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ParameterInfo &info);
	void insertDependency(ParameterInfo &info, const std::vector<std::string> &parameters, esint interval, int dimension = 0);

	template<int code, int nodes, int gpsize>
	void isotropicK(const Builder &builder, InstanceFiller &filler);
	template<int code, int nodes, int gpsize>
	void generalK(const Builder &builder, InstanceFiller &filler);

	void initTemperature();
	void computeThickness();
	void computeCoodinates();
	void computeJacobian();
	void computeConductivity();
	void computeConductivity2D();
	void computeConductivity3D();
	void computeAdvection();

	void computeStiffness();

	void assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const;

};

}

#endif /* SRC_PHYSICS_KERNELS_HEATTRANSFER_KERNEL_H_ */
