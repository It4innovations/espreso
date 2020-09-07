
#include "metadata.h"
#include "conditions.h"
#include "esinfo/eslog.hpp"

#include <sstream>

using namespace espreso;

Unit& Unit::add(UnitLibrary unit, int exponent)
{
	this->elements.push_back({unit, exponent});
	return *this;
}

void ECFMetaData::checkdescription(const std::string &name, size_t size) const
{
	if (description.size() != size) {
		eslog::internalFailure("'%s' has incorrect number of descriptions.\n", name.c_str());
	}
}

void ECFMetaData::checkdatatype(const std::string &name, size_t size) const
{
	if (datatype.size() != size) {
		eslog::internalFailure("'%s' has incorrect number of datatypes.\n", name.c_str());
	}
}

void ECFMetaData::checkpattern(const std::string &name, size_t size) const
{
	if (pattern.size() != size) {
		eslog::internalFailure("'%s' has incorrect number of pattern values.\n", name.c_str());
	}
}

ECFMetaData& ECFMetaData::addconstraint(const ECFAbstractCondition &condition) 
{ 
	delete this->condition; this->condition = condition.copy(); return *this; 
}

ECFMetaData& ECFMetaData::removeconstraint() 
{ 
	delete this->condition; this->condition = new ECFCondition(); return *this; 
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

ECFMetaData::ECFMetaData(): visibleObjectName(false), tensor(NULL), 
tensor_row(-1), tensor_column(-1), exporting(true), range_begin(0), 
range_end(0), gui_type(ECFGUIType::STATIC)
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
  unit(other.unit), exporting(other.exporting), range_begin(other.range_begin),
  range_end(other.range_end), gui_type(other.gui_type),
  pattern_name(other.pattern_name), pattern_item_name(other.pattern_item_name),
  condition(other.condition->copy()),
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
  unit(std::move(other.unit)), exporting(other.exporting), range_begin(other.range_begin),
  range_end(other.range_end), gui_type(other.gui_type), 
  pattern_name(other.pattern_name), pattern_item_name(other.pattern_item_name),
  condition(other.condition), isallowed(std::move(other.isallowed)), 
  ismandatory(std::move(other.ismandatory))
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
		exporting = other.exporting;
		range_begin = other.range_begin;
		range_end = other.range_end;
		gui_type = other.gui_type;
		pattern_name = other.pattern_name;
		pattern_item_name = other.pattern_item_name;
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

