
#ifndef SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_OPT_H_
#define SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_OPT_H_

#include "kernel.parameters.h"
#include "physics/assembler/operator.h"
#include "physics/kernels/basefunctions/basefunctions.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct Evaluator;
struct ECFExpression;
struct ECFExpressionVector;
struct ExpressionsToElements;
struct ExpressionsToBoundary;
class ConvectionConfiguration;
struct SolverDataProvider;

class KernelOpt
{
public:
	virtual ~KernelOpt() {}
	SolverDataProvider *solverDataProvider;

	virtual void nextSubstep() =0;
	virtual void solutionChanged() =0;
	virtual void processSolution() =0;

	virtual void updateStiffness(double *K, esint *perm, int interval) =0;

	virtual void updateStiffness(double *K, esint *perm, int region, int interval) =0;
	virtual void updateRHS(double *RHS, esint *perm, int region, int interval) =0;

protected:
	KernelOpt(SolverDataProvider *provider): solverDataProvider(provider)
	{
		BaseFunctions::setBaseFunctions();
	}

	void setMaterials(const std::map<std::string, std::string> &settings);
	void printMaterials(const std::map<std::string, std::string> &settings);
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

	void examineMaterialParameter(const std::string &material, const std::string &name, const ECFExpression &settings, ExpressionsToElements &builder, int dimension);

	template<class TSecond>
	void examineElementParameter(const std::string &name, const std::map<std::string, TSecond> &settings, ExpressionsToElements &builder, int dimension, std::function<const Evaluator*(const TSecond &expr)> getevaluator);
	void examineElementParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExpressionsToElements &builder);
	void examineElementParameter(const std::string &name, const std::map<std::string, ECFExpressionVector> &settings, ExpressionsToElements &builder, int dimension);

	void examineBoundaryParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExpressionsToBoundary &builder);
	void examineBoundaryParameter(const std::string &name, const std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection);

	std::vector<OperatorBuilder*> operators;
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_KERNELS_KERNEL_OPT_H_ */
