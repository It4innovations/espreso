
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "parameters.h"
#include "analysis/assembler/operator.h"

#include <memory>
#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct Evaluator;
struct ECFExpression;
struct ECFExpressionVector;
class ConvectionConfiguration;
struct SolverDataProvider;
class Vectors;

class Assembler
{
public:
	Assembler();
	virtual ~Assembler() {}

	virtual bool hasKernel(int domain) =0;

	void addParameter(ParameterData &parameter)
	{
		parameters.push_back(&parameter);
	}

	std::vector<std::vector<std::unique_ptr<ActionOperator> > > actionOps, actionRes;

protected:

	void updateVersions();

	void printParamtereStats(const char* name, ParameterData &parameter);
	void printParamtereStats(const char* name, NamedData *data);

	void setMaterials(const std::map<std::string, std::string> &settings);
	void printMaterials(const std::map<std::string, std::string> &settings);
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

	bool examineMaterialParameter(const std::string &material, const std::string &name, const ECFExpression &settings, ExternalValue &value, int dimension);

	template<class TSecond>
	bool examineElementParameter(const std::string &name, const std::map<std::string, TSecond> &settings, ExternalValue &value, int dimension, std::function<Evaluator*(const TSecond &expr)> getevaluator);
	bool examineElementParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExternalValue &value);
	bool examineElementParameter(const std::string &name, const std::map<std::string, ECFExpressionVector> &settings, ExternalValue &value, int dimension);

	void examineBoundaryParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ExternalValue &value);
	void examineBoundaryParameter(const std::string &name, const std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection);

	std::vector<OperatorBuilder*> builders, results;
	std::vector<ParameterData*> parameters;
	int version;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */
