
#include "json.h"

#include <sstream>
#include "config/conditions.h"
#include "config/ecf/material/tensor.h"

using namespace espreso;

ECFJSONExport::ECFJSONExport(ECFObject* obj, std::ostream& outputStream)
: m_stream(outputStream),
  m_rootValue(NULL),
  m_rootObject(obj),
  m_parent(NULL),
  m_tensorId(0)
{
	this->init();
}

ECFJSONExport::ECFJSONExport(ECFValue* value, ECFObject* parent, std::ostream& outputStream)
: m_stream(outputStream),
  m_rootValue(value),
  m_rootObject(NULL),
  m_parent(parent),
  m_tensorId(0)
{
	this->init();
}

void ECFJSONExport::exportToStream()
{
	ECFParameter* p = this->m_rootObject ? this->m_rootObject : this->m_rootValue;
	this->printParameter(p, NULL, 0, false);
}

std::string ECFJSONExport::spaces(const int indent)
{
	std::stringstream ss;
	for (int i = 0; i < indent; i++) ss << " ";
	return ss.str();
}

void ECFJSONExport::printKeyValuePair(const std::string& key,
	const std::string& val,
	const int indent,
	const char * valueEnclosed)
{
	this->m_stream
		<< this->spaces(indent)
		<< "\"" << key << "\""
		<< ": " << valueEnclosed
		<< val
		<< valueEnclosed;
}

void ECFJSONExport::printKeyArrayPair(const std::string& key,
	const std::vector<std::string>& arr, const int indent,
	const char* itemEnclosed)
{
	this->m_stream
		<< this->spaces(indent)
		<< "\"" << key << "\""
		<< ": "
		<< "[";

	std::vector<std::string> aux_arr;
	for (auto e = arr.begin(); e != arr.end(); ++e)
	{
		aux_arr.push_back(itemEnclosed);
		aux_arr.push_back(*e);
		aux_arr.push_back(itemEnclosed);
		aux_arr.push_back(", ");
	}
	if (aux_arr.size()) aux_arr.pop_back();
	for (auto e = aux_arr.begin(); e != aux_arr.end(); ++e)
	{
		this->m_stream << *e;
	}

	this->m_stream << "]";
}

void ECFJSONExport::printConstraint(ECFParameter* p, int indent)
{
	if (p->metadata.condition->isset())
	{
		this->m_stream << "," << std::endl;
		ECFObject* scope = m_rootObject ? m_rootObject : m_parent;
		std::string condition;
		scope->forEachParameters(
			[&] (const ECFParameter *conditionalParameter,
				const std::string& prefix
			) {
			if (p->metadata.condition->match(
					conditionalParameter->data()
			))
			{
				p->metadata.condition->bind(conditionalParameter, prefix);
				condition = p->metadata.condition->compose();
			}
		}, "", false, true);
		this->printKeyValuePair("condition", condition, indent);
	}
}

void ECFJSONExport::printRegisteredTensors(ECFParameter* p, int indent)
{
	if (p->metadata.tensors.size())
	{
		this->m_stream << "," << std::endl;

		std::vector<std::string> ids;
		for (auto it = p->metadata.tensors.cbegin();
			it != p->metadata.tensors.end();
			++it)
		{
			if ( m_tensors.find(*it) == m_tensors.end() )
			{ m_tensors[*it] = m_tensorId++; }
			ids.push_back(std::to_string(m_tensors[*it]));
		}
		this->printKeyArrayPair(
			"tensors",
			ids,
			indent
		);
	}
}

void ECFJSONExport::printTensor(ECFParameter* p, int indent)
{
	if (p->metadata.tensor != NULL)
	{
		if ( m_tensors.find(p->metadata.tensor) == m_tensors.end() )
		{ m_tensors[p->metadata.tensor] = m_tensorId++; }
		this->m_stream << "," << std::endl;
		auto tensor = [&] (int p_indent = 0) {
			this->printKeyValuePair(
				"id",
				std::to_string(m_tensors[p->metadata.tensor]),
				p_indent);
			this->m_stream << "," << std::endl;
			this->printKeyValuePair(
				"row",
				std::to_string(p->metadata.tensor_row),
				p_indent);
			this->m_stream << "," << std::endl;
			this->printKeyValuePair(
				"column",
				std::to_string(p->metadata.tensor_column),
				p_indent);
			this->m_stream << "," << std::endl;
			this->printKeyValuePair(
				"size",
				std::to_string(p->metadata.tensor->size),
				p_indent);
		};
		this->printKeyObjectPair("tensor", tensor, indent);
	}
}

void ECFJSONExport::printKeyObjectPair(const std::string& key,
	const std::function<void(int)>& printObjContent,
	const int indent)
{
	this->m_stream << this->spaces(indent);
	this->m_stream << "\"" << key << "\"" << ": " << "{" << std::endl;

	printObjContent(indent + 2);

	this->m_stream << std::endl << spaces(indent + 2) << "}" << std::endl;
}

void ECFJSONExport::printMetaData(ECFParameter* p, const int indent)
{
	ECFMetaData* md = &p->metadata;
	int ind = indent + 2;

	this->printKeyValuePair("name", md->name, ind);
	this->m_stream << "," << std::endl;

	this->printKeyArrayPair("description", md->description, ind);
	this->m_stream << "," << std::endl;

	auto datatypesToStringVector = [&](const std::vector<ECFDataType>& dts) -> std::vector<std::string>{
		std::vector<std::string> ret;
		for (auto d = dts.begin(); d != dts.end(); ++d)
		{
			ret.push_back(this->m_datatypes[*d]);
		}
		return ret;
	};
	this->printKeyArrayPair("datatype", datatypesToStringVector(md->datatype), ind);

	this->m_stream << "," << std::endl;

	this->printKeyArrayPair("pattern", md->pattern, ind);

	this->m_stream << "," << std::endl;

	auto printOptions = [&] () {
		std::vector<std::string> ret;
		for (auto option = md->options.begin(); option != md->options.end(); ++option)
		{
			std::stringstream ss;
			ss
				<< "{ \"name\": \""
				<< option->name
				<< "\", \"description\": \""
				<< option->description
				<< "\" }";
			ret.push_back(ss.str());
		}

		return ret;
	};
	this->printKeyArrayPair("options", printOptions(), ind, "");

	this->m_stream << "," << std::endl;

	this->printKeyArrayPair("variables", md->variables, ind);

	this->m_stream << "," << std::endl;

	this->printKeyValuePair("unit", md->unit.unit(), ind);

	this->printConstraint(p, ind);

	this->printRegisteredTensors(p, ind);

	this->printTensor(p, ind);

	this->m_stream << "," << std::endl;

	if (md->visibleObjectName)
		printKeyValuePair("visibleObjectName", "true", ind, "");
	else
		printKeyValuePair("visibleObjectName", "false", ind, "");
}

void ECFJSONExport::printValue(ECFParameter* val, const int indent)
{
	this->printKeyValuePair("value", val->getValue(), indent + 2);
	this->m_stream << "," << std::endl;

	auto printMetadata = [&](const int p_indent = 0) {
		printMetaData(val, p_indent);
	};
	this->printKeyObjectPair("metadata", printMetadata, indent + 2);
}

void ECFJSONExport::printParameter(ECFParameter* param,
	int * const id,
	const int indent,
	const bool printKey)
{
	this->m_stream
			<< std::endl
			<< this->spaces(indent);

	if (printKey)
	{
		this->m_stream << "\"";
		this->m_stream << param->name;
		if
		(
			param->metadata.datatype.size() &&
			(
				param->metadata.datatype[0] == ECFDataType::BEGINBLOCK ||
				param->metadata.datatype[0] == ECFDataType::ENDBLOCK ||
				param->metadata.datatype[0] == ECFDataType::BEGINCOLLAPSEBLOCK ||
				param->metadata.datatype[0] == ECFDataType::SPACE ||
				param->metadata.datatype[0] == ECFDataType::SEPARATOR
			)
		)
		{
			this->m_stream << "-" << *id;
			(*id)++;
		}
		this->m_stream << "\"" << ": ";
	}

	this->m_stream << "{" << std::endl;

	if (param->isValue())
	{
		this->printValue(param, indent);
	}
	else if (param->isObject())
	{
		int id = 0;

		ECFObject* obj = static_cast<ECFObject*>(param);

		auto printMetadata = [&](const int p_indent = 0) {
			printMetaData(obj, p_indent);
		};
		this->printKeyObjectPair("metadata", printMetadata, indent + 2);
		this->m_stream << "," << std::endl;

		ECFParameter* pattern = obj->getPattern();
		if (pattern)
		{
			pattern->name = "pattern";
			this->printParameter(pattern, &id, indent + 2);
		}

		auto p = obj->parameters.begin();
		if (p != obj->parameters.end()) {
			if (pattern) this->m_stream << "," << std::endl;
			this->printParameter(*p, &id, indent + 2);
			p++;
		}
		for (; p != obj->parameters.end(); ++p)
		{
			this->m_stream << "," << std::endl;
			this->printParameter(*p, &id, indent + 2);
		}
	}
	else
	{
		this->printValue(param, indent);
	}

	this->m_stream << std::endl << this->spaces(indent) << "}";
}

void ECFJSONExport::init()
{
	this->m_datatypes[ECFDataType::BOOL] = "bool";
	this->m_datatypes[ECFDataType::STRING] = "string";
	this->m_datatypes[ECFDataType::INTEGER] = "integer";
	this->m_datatypes[ECFDataType::POSITIVE_INTEGER] = "positive_integer";
	this->m_datatypes[ECFDataType::NONNEGATIVE_INTEGER] = "nonnegative_integer";
	this->m_datatypes[ECFDataType::FLOAT] = "float";
	this->m_datatypes[ECFDataType::ENUM_FLAGS] = "enum_flags";
	this->m_datatypes[ECFDataType::OPTION] = "option";
	this->m_datatypes[ECFDataType::REGION] = "region";
	this->m_datatypes[ECFDataType::BOUNDARY_REGION] = "boundary_region";
	this->m_datatypes[ECFDataType::ELEMENTS_REGION] = "elements_region";
	this->m_datatypes[ECFDataType::MATERIAL] = "material";
	this->m_datatypes[ECFDataType::LOAD_STEP] = "load_step";
	this->m_datatypes[ECFDataType::EXPRESSION] = "expression";
	this->m_datatypes[ECFDataType::TENSOR] = "tensor";
	this->m_datatypes[ECFDataType::INTERVAL] = "interval";
	this->m_datatypes[ECFDataType::SPACE] = "space";
	this->m_datatypes[ECFDataType::SEPARATOR] = "separator";
	this->m_datatypes[ECFDataType::BEGINBLOCK] = "beginblock";
	this->m_datatypes[ECFDataType::ENDBLOCK] = "endblock";
	this->m_datatypes[ECFDataType::BEGINCOLLAPSEBLOCK] = "begincollapseblock";
}
