
#include "metadata.h"
#include "conditions.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

std::string SIUnit::unit() const
{
	auto exponent = [] (const std::string &value, int exponent) {
		if (exponent == 0) {
			return std::string();
		}
		return value + std::to_string(exponent);
	};

	return
			exponent("m", metre) +
			exponent("kg", kilogram) +
			exponent("s", second) +
			exponent("A", ampere) +
			exponent("K", kelvin) +
			exponent("mol", mole) +
			exponent("cd", candela);
}

void ECFMetaData::checkdescription(const std::string &name, size_t size) const
{
	if (description.size() != size) {
		eslog::globalerror("ESPRESO internal error: '%s' has incorrect number of descriptions.\n", name.c_str());
	}
}

void ECFMetaData::checkdatatype(const std::string &name, size_t size) const
{
	if (datatype.size() != size) {
		eslog::globalerror("ESPRESO internal error: '%s' has incorrect number of datatypes.\n", name.c_str());
	}
}

void ECFMetaData::checkpattern(const std::string &name, size_t size) const
{
	if (pattern.size() != size) {
		eslog::globalerror("ESPRESO internal error: '%s' has incorrect number of pattern values.\n", name.c_str());
	}
}

ECFMetaData& ECFMetaData::addconstraint(const ECFAbstractCondition &condition) 
{ 
	delete this->condition; this->condition = condition.copy(); return *this; 
}

template <typename TType>
static std::vector<TType> getsuffix(size_t start, const std::vector<TType> &data)
{
	if (start < data.size()) {
		return std::vector<TType>(data.begin() + start, data.end());
	}
	return std::vector<TType>();
}

ECFMetaData ECFMetaData::suffix(size_t start) const
{
	ECFMetaData ret(*this);
	ret.setdescription(getsuffix(start, description));
	ret.setdatatype(getsuffix(start, datatype));
	ret.setpattern(getsuffix(start, pattern));
	ret.tensor = NULL;
	return ret;
}

ECFMetaData::ECFMetaData(): visibleObjectName(false), tensor(NULL), tensor_row(-1), tensor_column(-1)
{
	condition = new ECFCondition();
	isallowed = [] () { return true; };
	ismandatory = [] () { return true; };
}

ECFMetaData::ECFMetaData(const ECFMetaData &other)
: name(other.name), description(other.description),
  datatype(other.datatype), pattern(other.pattern),
  options(other.options), variables(other.variables),
  visibleObjectName(other.visibleObjectName),
  tensor(other.tensor), tensor_row(other.tensor_row),
  tensor_column(other.tensor_column), tensors(other.tensors),
  unit(other.unit), condition(other.condition->copy()),
  isallowed(other.isallowed), ismandatory(other.ismandatory)
{

}

ECFMetaData::ECFMetaData(ECFMetaData &&other)
: name(std::move(other.name)), description(std::move(other.description)),
  datatype(std::move(other.datatype)), pattern(std::move(other.pattern)),
  options(std::move(other.options)), variables(std::move(other.variables)),
  visibleObjectName(other.visibleObjectName),
  tensor(std::move(other.tensor)), tensor_row(other.tensor_row),
  tensor_column(other.tensor_column), tensors(other.tensors),
  unit(std::move(other.unit)), condition(other.condition),
  isallowed(std::move(other.isallowed)), ismandatory(std::move(other.ismandatory))
{
	other.condition = NULL;
}

ECFMetaData& ECFMetaData::operator=(const ECFMetaData &other)
{
	if (this != &other) {
		name = other.name;
		description = other.description;
		datatype = other.datatype;
		pattern = other.pattern;
		options = other.options;
		variables = other.variables;
		tensor = other.tensor;
		tensor_row = other.tensor_row;
		tensor_column = other.tensor_column;
		tensors = other.tensors;
		unit = other.unit;
		isallowed = other.isallowed;
		ismandatory = other.ismandatory;
		if (condition != NULL) {
			delete condition;
		}
		condition = other.condition->copy();
		visibleObjectName = other.visibleObjectName;
	}
	return *this;
}

ECFMetaData::~ECFMetaData()
{
	if (condition != NULL) {
		delete condition;
	}
}

