
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "parameters.h"
#include "analysis/assembler/controller.h"
#include "analysis/assembler/operator.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct ECFExpression;
struct ECFExpressionVector;
class ConvectionConfiguration;
struct ImpedanceConfiguration;

class Assembler
{
public:
	Assembler();
	virtual ~Assembler() {}

	virtual bool hasKernel(int domain) =0;

	ParameterController controller;
	std::vector<std::vector<std::unique_ptr<ActionOperator> > > elementOps, elementRes;
	std::vector<std::vector<std::vector<std::unique_ptr<ActionOperator> > > > boundaryOps, boundaryRes;

protected:
	void iterate();

	void printParameterStats(const char* name, ParameterData &parameter);
	void printParameterStats(const char* name, NamedData *data);

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
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ImpedanceConfiguration> &settings, ExternalBoundaryValue &impedance);

	int version;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */
