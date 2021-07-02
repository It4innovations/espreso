
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "parameters.h"
#include "analysis/assembler/operator.h"

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
class Vectors;

class Assembler
{
public:
	virtual ~Assembler() {}

	virtual bool hasKernel(int domain) =0;

	void addParameter(ParameterData &parameter)
	{
		parameters.push_back(&parameter);
	}

protected:

	void updateVersions();

	void printParamtereStats(const char* name, ParameterData &parameter);
	void printParamtereStats(const char* name, NamedData *data);

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

	std::vector<OperatorBuilder*> builders, results;
	std::vector<ParameterData*> parameters;
	int version;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */
