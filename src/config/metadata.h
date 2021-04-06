
#ifndef SRC_CONFIG_METADATA_H_
#define SRC_CONFIG_METADATA_H_

#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct TensorConfiguration;
struct ECFAbstractCondition;

enum class ECFDataType {
	BOOL,
	STRING,
	INTEGER,
	POSITIVE_INTEGER,
	NONNEGATIVE_INTEGER,
	FLOAT,
	ENUM_FLAGS,
	OPTION,
	REGION,
	BOUNDARY_REGION,
	ELEMENTS_REGION,
	MATERIAL,
	LOAD_STEP,
	EXPRESSION,
	TENSOR,
	INTERVAL,
	SPACE,
	SEPARATOR,
	BEGINBLOCK,
	ENDBLOCK,
	BEGINCOLLAPSEBLOCK
};

enum class DIMENSION {
	D1,
	D2,
	D3,
	Z
};

struct ECFOption {
	std::string name;
	std::string description;
	std::function<bool(void)> isallowed;

	ECFOption& setname(const std::string &name) { this->name = name; return *this; }
	ECFOption& setdescription(const std::string &description) { this->description = description; return *this; }
	ECFOption& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }

	ECFOption() { isallowed = [] () { return true; }; }
};

struct SIUnit {
	int metre, kilogram, second, ampere, kelvin, mole, candela;

	SIUnit(): metre(0), kilogram(0), second(0), ampere(0), kelvin(0), mole(0), candela(0) {}

	SIUnit(int metre, int kilogram, int second, int ampere, int kelvin, int mole, int candela)
	: metre(metre), kilogram(kilogram), second(second), ampere(ampere), kelvin(kelvin), mole(mole), candela(candela)
	{}

	std::string unit() const;
};


struct ECFMetaData {
	std::string name;
	std::vector<std::string> description;
	std::vector<ECFDataType> datatype;
	std::vector<std::string> pattern;
	std::vector<ECFOption> options;
	std::vector<std::string> variables;
	bool visibleObjectName;
	TensorConfiguration *tensor;
	int tensor_row;
	int tensor_column;
	std::vector<TensorConfiguration*> tensors;
	SIUnit unit;

	ECFAbstractCondition *condition;

	std::function<bool(void)> isallowed;
	std::function<bool(void)> ismandatory;

	ECFMetaData& setname(const std::string &name) { this->name = name; return *this; }
	ECFMetaData& setdescription(const std::vector<std::string> &description) { this->description = description; return *this; }
	ECFMetaData& setdatatype(const std::vector<ECFDataType> &datatype) { this->datatype = datatype; return *this; }
	ECFMetaData& setpattern(const std::vector<std::string> &pattern) { this->pattern = pattern; return *this; }
	ECFMetaData& setvariables(const std::vector<std::string> &variables) { this->variables = variables; return *this; }
	ECFMetaData& settensor(TensorConfiguration &tensor) { this->tensor = &tensor; return *this; }
	ECFMetaData& settensor(TensorConfiguration &tensor, int row, int column) 
	{ this->tensor = &tensor; this->tensor_row = row; this->tensor_column = column; return *this; }
	ECFMetaData& registerTensor(TensorConfiguration& tensor)
	{ this->tensors.push_back(&tensor); return *this; }
	ECFMetaData& setunit(const SIUnit &unit) { this->unit = unit; return *this; }
	ECFMetaData& displayObjectName() { this->visibleObjectName = true; return *this; }
	ECFMetaData& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }
	ECFMetaData& mandatoryonly(std::function<bool(void)> ismandatory) { this->ismandatory = ismandatory; return *this; }
	ECFMetaData& addconstraint(const ECFAbstractCondition &condition);

	ECFMetaData& addoption(const ECFOption &option) { options.push_back(option); return *this; }

	static std::vector<std::string> getcoordinatevariables() { return { "X", "Y", "Z" }; }
	static std::vector<std::string> getboundaryconditionvariables() { return { "X", "Y", "Z", "INITIAL_TEMPERATURE", "TEMPERATURE", "TIME", "FREQUENCY" }; }
	static std::vector<std::string> getharmonicvariables() { return { "X", "Y", "Z", "TIME", "FREQUENCY" }; }
	static std::vector<std::string> getmaterialvariables() { return { "X", "Y", "Z", "TIME", "TEMPERATURE" }; }

	void checkdescription(const std::string &name, size_t size) const;
	void checkdatatype(const std::string &name, size_t size) const;
	void checkpattern(const std::string &name, size_t size) const;

	ECFMetaData suffix(size_t start) const;

	ECFMetaData();
	ECFMetaData(const ECFMetaData &other);
	ECFMetaData(ECFMetaData &&other);

	ECFMetaData& operator=(const ECFMetaData &other);
	~ECFMetaData();
};

}



#endif /* SRC_CONFIG_METADATA_H_ */
