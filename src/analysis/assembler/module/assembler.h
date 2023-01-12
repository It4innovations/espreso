
#ifndef SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_
#define SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_

#include "parameters.h"
#include "analysis/assembler/controller.h"
#include "analysis/assembler/operator.h"
#include "math/primitives/vector_sparse.h"

#include <vector>
#include <map>
#include <functional>

namespace espreso {

struct ECFExpression;
struct ECFExpressionVector;
struct ECFExpressionOptionalVector;
class ConvectionConfiguration;
struct ImpedanceConfiguration;
struct PointSourceConfiguration;

class Assembler
{
public:
	Assembler();
	virtual ~Assembler() {}

	ParameterController controller;
	std::vector<std::vector<std::unique_ptr<ActionOperator> > > elementOps, elementFiller, elementRes;
	std::vector<std::vector<std::vector<std::unique_ptr<ActionOperator> > > > boundaryOps, boundaryFiller, boundaryRes;

protected:
	void iterate();
	void fill();
	void results();

	void printElementVolume(const ElementParameter<egps> &weight, const ElementParameter<egps> &jacobian);
	void printBoundarySurface(const BoundaryParameter<egps> &weight, const BoundaryParameter<egps> &jacobian);
	void printParameterStats(ParameterData &parameter);
	void printParameterStats(NamedData *data);

	void printMaterials(const std::map<std::string, std::string> &settings);
	template<typename Ttype>
	void validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings);

	bool examineMaterialParameter(const std::string &material, const std::string &name, ECFExpression &settings, ExternalElementValue &externalValue, int dimension);

	template<class TSecond>
	bool examineElementParameter(const std::string &name, std::map<std::string, TSecond> &settings, ExternalElementValue &externalValue, int dimension, std::function<ECFExpression*(TSecond &expr)> getExpr);
	bool examineElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalElementValue &externalValue);
	bool examineElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, ExternalElementValue &externalValue, int dimension);

	template<class TSecond>
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, TSecond> &settings, ExternalBoundaryValue &externalValue, int dimension, std::function<ECFExpression*(TSecond &expr)> getExpr);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings, ExternalBoundaryValue &value);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ConvectionConfiguration> &settings, ParametersConvection &convection);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ImpedanceConfiguration> &settings, ExternalBoundaryValue &impedance);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, PointSourceConfiguration> &settings, ExternalBoundaryValue &point_source);

	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, ExternalBoundaryValue &value, int dimension);
	bool examineBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings, ExternalBoundaryValue &value, int dimension);

	int version;
};

template <typename T>
static void reset(T *t)
{
	if (t) {
		t->set(0);
	}
}

template <typename T, typename ...Other>
static void reset(T *t, Other... other)
{
	reset(t);
	reset(other...);
}

template <typename T>
static void update(T *t)
{
	if (t) {
		t->synchronize();
		t->touched = true;
	}
}

template <typename T, typename ...Other>
static void update(T *t, Other... other)
{
	update(t);
	update(other...);
}

}

#endif /* SRC_ANALYSIS_ASSEMBLER_MODULE_ASSEMBLER_H_ */
