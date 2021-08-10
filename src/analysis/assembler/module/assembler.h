
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "parameters.h"
#include "analysis/assembler/operator.h"

#include "analysis/composer/elementmapping.h"
#include "math2/primitives/vector_sparse.h"
#include "math2/primitives/matrix_info.h"
#include "math2/generalization/matrix_base.h"

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

	std::vector<std::vector<std::unique_ptr<ActionOperator> > > elementOps, elementRes;
	std::vector<std::vector<std::vector<std::unique_ptr<ActionOperator> > > > boundaryOps, boundaryRes;

protected:
	void iterate();

	void updateVersions();

	void printParamtereStats(const char* name, ParameterData &parameter);
	void printParamtereStats(const char* name, NamedData *data);

	void setMaterials(const std::map<std::string, std::string> &settings);
	void printMaterials(const std::map<std::string, std::string> &settings);
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

	bool examineMaterialParameter(const std::string &material, const std::string &name, ECFExpression &settings, ExternalElementValue &externalValue, int dimension);

	template<class TSecond>
	bool examineElementParameter(const std::string &name, std::map<std::string, TSecond> &settings, ExternalElementValue &externalValue, int dimension, std::function<ECFExpression*(TSecond &expr)> getExpr);
	bool examineElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalElementValue &externalValue);
	bool examineElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, ExternalElementValue &externalValue, int dimension);

	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalBoundaryValue &value);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection);

	std::vector<ParameterData*> parameters;
	int version;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */