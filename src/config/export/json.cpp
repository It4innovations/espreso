
#include "json.h"

#include <sstream>
#include "config/conditions.h"
#include "config/ecf/material/tensor.h"
#include "config/holders/valueholder.h"

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
<<<<<<< HEAD
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

void ECFJSONExport::printKeyValuePair(const std::string& key,
	const int val,
	const int indent,
	const char * valueEnclosed)
{
	this->m_stream
		<< this->spaces(indent)
		<< "\"" << key << "\""
		<< ": " << val;
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

void ECFJSONExport::printUnit(ECFParameter* p, const int indent)
{
	if (!p->metadata.unit.elements.size()) return;

	this->m_stream << "," << std::endl;

	std::vector<Unit::UnitElement>& units = p->metadata.unit.elements;

	std::stringstream ss;
	ss << "[";
	auto u = units.begin();
	ss << "{ \"unit\": \"" << m_units[u->unit] << "\", \"exponent\": \"";
	ss << std::to_string(u->exponent) << "\"}";
	for (++u; u != units.end(); ++u)
	{
		ss << ", { \"unit\": \"" << m_units[u->unit] << "\", \"exponent\": \"";
		ss << std::to_string(u->exponent) << "\"}";
	}
	ss << "]";

	this->printKeyValuePair("unit", ss.str(), indent, "");
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
				p->metadata.condition->bind(conditionalParameter, "currentForm");
				condition = p->metadata.condition->compose();
			}
		}, "", false, true);
		this->printKeyValuePair("conditions", condition, indent);
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

	if (md->description.size())
	{ 
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("description", md->description[0], ind); 
	}

	// this->printKeyArrayPair("pattern", md->pattern, ind);
	// this->m_stream << "," << std::endl;

	if (md->datatype.size() 
		&& (md->datatype[0] == ECFDataType::OPTION 
			|| md->datatype[0] == ECFDataType::ENUM_FLAGS))
	{
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("defaultValue", p->getValue(), ind);
		this->m_stream << "," << std::endl;
		auto printOptions = [&] (int indent) {
			auto option = md->options.begin();
			this->m_stream  << this->spaces(indent)
				<< "\"" << option->name << "\": "
				<< "\"" << option->description << "\"";
			for (++option; option != md->options.end(); ++option)
			{
				this->m_stream << "," << std::endl << this->spaces(indent)
				<< "\"" << option->name << "\": "
				<< "\"" << option->description << "\"";
			}
		};
		this->printKeyObjectPair("values", printOptions, ind);
	}

	if (md->datatype.size() && md->datatype[0] == ECFDataType::SEPARATOR)
	{
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("dashed", "true", ind);
	}
	if (md->datatype.size() && 
		(md->datatype[0] == ECFDataType::INTEGER
		|| md->datatype[0] == ECFDataType::POSITIVE_INTEGER
		|| md->datatype[0] == ECFDataType::NONNEGATIVE_INTEGER))
	{
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("min", md->range_begin, ind);
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("max", md->range_end, ind);
	}
	if (p->isObject() && md->gui_type == ECFGUIType::GROUP_COLLAPSED)
	{
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("type", "collapse", ind);
	}
	if (p->isObject() && md->gui_type == ECFGUIType::GROUP_BLOCK)
	{
		this->m_stream << "," << std::endl;
		this->printKeyValuePair("type", "border", ind);
	}
	// this->m_stream << "," << std::endl;

	// this->printKeyArrayPair("variables", md->variables, ind);

	// this->m_stream << "," << std::endl;

	this->printUnit(p, ind);

	this->printConstraint(p, ind);

	
	// this->printRegisteredTensors(p, ind);

	// this->printTensor(p, ind);

	// this->m_stream << "," << std::endl;

	// if (md->visibleObjectName)
	// 	printKeyValuePair("visibleObjectName", "true", ind, "");
	// else
	// 	printKeyValuePair("visibleObjectName", "false", ind, "");
=======
    std::stringstream ss;
    for (int i = 0; i < indent; i++) ss << " ";
    return ss.str();
}

void ECFJSONExport::printKeyValuePair(const std::string& key, 
    const std::string& val, 
    const int indent)
{
    this->m_stream 
        << this->spaces(indent) 
        << "\"" << key << "\""
        << ": \"" 
        << val 
        << "\"";
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

void ECFJSONExport::printKeyObjectPair(const std::string& key, 
    const std::function<void(int)>& printObjContent, 
    const int indent)
{
    this->m_stream << this->spaces(indent);
    this->m_stream << "\"" << key << "\"" << ": " << "{" << std::endl;
    
    printObjContent(indent + 2);

    this->m_stream << std::endl << spaces(indent + 2) << "}" << std::endl;
}

void ECFJSONExport::printMetaData(ECFMetaData *md, const int indent)
{
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
>>>>>>> ENH: [WEB-GUI] JSON export - SI unit metadata basic export
}

void ECFJSONExport::printValue(ECFParameter* val, const int indent)
{
	this->printKeyValuePair("type", 
		this->m_datatypes[val->metadata.datatype[0]], 
		indent + 2);

	this->m_stream << "," << std::endl;

	auto printMetadata = [&](const int p_indent = 0) {
		if (val->metadata.datatype.size() && 
			(val->metadata.datatype[0] == ECFDataType::INTEGER
			|| val->metadata.datatype[0] == ECFDataType::POSITIVE_INTEGER
			|| val->metadata.datatype[0] == ECFDataType::NONNEGATIVE_INTEGER))
		{
			this->printKeyValuePair("value", std::stoi(val->getValue()), p_indent + 2);
			this->m_stream << "," << std::endl;
		}
		else if(val->metadata.datatype.size() && val->metadata.datatype[0] == ECFDataType::BOOL)
		{
			std::string content = val->getValue();
			if (content.compare("FALSE") == 0)
				this->printKeyValuePair("value", 0, p_indent + 2);
			else
				this->printKeyValuePair("value", 1, p_indent + 2);
			this->m_stream << "," << std::endl;
		}
		else if(val->metadata.datatype.size() && 
			(val->metadata.datatype[0] != ECFDataType::OPTION) &&
			(val->metadata.datatype[0] != ECFDataType::ENUM_FLAGS))
		{
			this->printKeyValuePair("value", val->getValue(), p_indent + 2);
			this->m_stream << "," << std::endl;
		}
		
		printMetaData(val, p_indent);
	};
	this->printKeyObjectPair("metadata", printMetadata, indent + 2);
}

void ECFJSONExport::printParameter(ECFParameter* param,
	int * const id,
	const int indent,
	const bool printKey)
{
	if (!param->metadata.exporting) return;

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

		if (obj->metadata.gui_type == ECFGUIType::GROUP_BLOCK
			|| obj->metadata.gui_type == ECFGUIType::GROUP_COLLAPSED) {
			this->printObjectGroup(obj, &id, indent);
		}
		else {
			this->printObjectInTree(obj, &id, indent);
		}
	}
	else {
		this->printValue(param, indent);
	}

	this->m_stream << std::endl << this->spaces(indent) << "}" << std::endl;
}

void ECFJSONExport::printObjectInTree(ECFObject* object,
		int * const id,
		const int indent,
		const bool printKey)
{
	auto printGuiMetadata = [&](const int p_indent = 0) {
		if (object->metadata.gui_type == ECFGUIType::STATIC)
			printKeyValuePair("type", "static", indent + 2);
		else
			printKeyValuePair("type", "dynamic", indent + 2);
		this->m_stream << "," << std::endl;
		
		printKeyValuePair("name", object->metadata.name, indent + 2);
	};
	this->printKeyObjectPair("metadata", printGuiMetadata, indent + 2);
	this->m_stream << "," << std::endl;

	// Export parameters that should displayed in the object form
	std::vector<ECFParameter*> form_parameters; 
	for (auto p = object->parameters.begin(); 
		p != object->parameters.end(); ++p)
	{
		if (!(*p)->metadata.exporting) continue;
		if ((*p)->metadata.gui_type == ECFGUIType::FORM)
			form_parameters.push_back(*p);
	}
	if (form_parameters.size())
	{
		this->m_stream
		<< this->spaces(indent + 2)
		<< "\"form\": {" << std::endl;

		auto printMetadata = [&](const int p_indent = 0) {
			printMetaData(object, p_indent);
		};
		this->printKeyObjectPair("metadata", printMetadata, indent + 4);
		this->m_stream << "," << std::endl;

		this->m_stream
			<< this->spaces(indent + 4)
			<< "\"content\": {" << std::endl;

		auto p = form_parameters.begin();
		if (p != form_parameters.end()) 
		{
			this->printParameter(*p, id, indent + 6);
			p++;
		}
		for (; p != form_parameters.end(); ++p)
		{
			if (!(*p)->metadata.exporting) continue;
			this->m_stream << "," << std::endl;
			this->printParameter(*p, id, indent + 6);
		}

		this->m_stream << this->spaces(indent + 4) << "}" << std::endl;
		this->m_stream << std::endl << this->spaces(indent + 2) << "}," << std::endl;
	}

	// Export object's parameters depending on whether it is a static or 
	// a dynamic item in tree
	if (object->metadata.gui_type == ECFGUIType::STATIC)
		this->m_stream
		<< this->spaces(indent + 2)
		<< "\"content\": {" << std::endl;
	else
	{
		this->m_stream
		<< this->spaces(indent + 2)
		<< "\"createForm\": {" << std::endl;
		auto printMetadata = [&](const int p_indent = 0) {
			if (!object->metadata.pattern_name.empty())
				object->metadata.setname(object->metadata.pattern_name);
			printMetaData(object, p_indent);
		};
		this->printKeyObjectPair("metadata", printMetadata, indent + 4);
		this->m_stream << "," << std::endl;

		this->m_stream
			<< this->spaces(indent + 4)
			<< "\"content\": {" << std::endl;
	}

	this->printObjectContent(object, id, indent + 6, printKey);
	
	if (object->metadata.gui_type == ECFGUIType::DYNAMIC)
		this->m_stream << this->spaces(indent + 4) << "}" << std::endl;
	this->m_stream << std::endl << this->spaces(indent + 2) << "}" << std::endl;
}

void ECFJSONExport::printObjectContent(ECFObject* object,
	int * const id,
	const int indent,
	const bool printKey)
{
	ECFParameter* pattern = object->getPattern();
	if (pattern && !object->parameters.size() && object->metadata.datatype.size() == 1)
	{
		ECFObject* ptrn = static_cast<ECFObject*>(pattern);
		auto p = ptrn->parameters.begin();
		if (p != ptrn->parameters.end()) 
		{ this->printParameter(*p, id, indent); p++; }
		for (; p != ptrn->parameters.end(); ++p)
		{
			if (!(*p)->metadata.exporting) continue;
			this->m_stream << "," << std::endl;
			this->printParameter(*p, id, indent);
		}
		// pattern->name = "pattern";
		// this->printParameter(pattern, &id, indent + 2);
		// this->m_stream << 
		// 	this->spaces(indent + 2) << "," << std::endl;
	}
	else if (pattern && !object->parameters.size() && object->metadata.datatype.size() == 2)
	{
		std::string value = "";
		ECFValueHolder<std::string> val(value);
		val.name = "name";
		val.metadata.setname("Region")
					.setdescription({"Name of the region"})
					.setdatatype({ECFDataType::STRING});
		this->printParameter(&val, id, indent);
		this->m_stream << ",";
		pattern->setValue("");
		if (!object->metadata.pattern_item_name.empty())
			pattern->metadata.setname(object->metadata.pattern_item_name);
		this->printParameter(pattern, id, indent);
	}
	else {
		int parameters_printed = 0;
		auto p = object->parameters.begin();
		if (p != object->parameters.end())
		{
			if ((*p)->metadata.gui_type != ECFGUIType::FORM)
			{
				this->printParameter(*p, id, indent);
				parameters_printed++;
			}
			p++;
		}
		for (; p != object->parameters.end(); ++p)
		{
			if (!(*p)->metadata.exporting) continue;
			if ((*p)->metadata.gui_type == ECFGUIType::FORM) continue;
			if (parameters_printed) this->m_stream << "," << std::endl;
			this->printParameter(*p, id, indent);
			parameters_printed++;
		}
	}
}

void ECFJSONExport::printObjectGroup(ECFObject* object,
		int * const id,
		const int indent,
		const bool printKey)
{
	this->printKeyValuePair("type", 
		"group", 
		indent + 2);
	this->m_stream << "," << std::endl;
	auto printMetadata = [&](const int p_indent = 0) {
		printMetaData(object, p_indent);
	};
	this->printKeyObjectPair("metadata", printMetadata, indent + 2);
	this->m_stream << "," << std::endl;
	this->m_stream
		<< this->spaces(indent + 2)
		<< "\"content\": {" << std::endl;

	this->printObjectContent(object, id, indent + 4);

	this->m_stream << std::endl << this->spaces(indent + 2) << "}" << std::endl;
}

void ECFJSONExport::init()
{
	this->m_datatypes[ECFDataType::BOOL] = "switch";
	this->m_datatypes[ECFDataType::STRING] = "textbox";
	this->m_datatypes[ECFDataType::INTEGER] = "slider";
	this->m_datatypes[ECFDataType::POSITIVE_INTEGER] = "slider";
	this->m_datatypes[ECFDataType::NONNEGATIVE_INTEGER] = "slider";
	this->m_datatypes[ECFDataType::FLOAT] = "textbox";
	this->m_datatypes[ECFDataType::ENUM_FLAGS] = "combobox";
	this->m_datatypes[ECFDataType::OPTION] = "combobox";
	this->m_datatypes[ECFDataType::REGION] = "textbox";
	this->m_datatypes[ECFDataType::BOUNDARY_REGION] = "boundary_region";
	this->m_datatypes[ECFDataType::ELEMENTS_REGION] = "elements_region";
	this->m_datatypes[ECFDataType::MATERIAL] = "textbox";
	this->m_datatypes[ECFDataType::LOAD_STEP] = "load_step";
	this->m_datatypes[ECFDataType::EXPRESSION] = "mathbox";
	this->m_datatypes[ECFDataType::TENSOR] = "tensor";
	this->m_datatypes[ECFDataType::INTERVAL] = "interval";
	// this->m_datatypes[ECFDataType::SPACE] = "space";
	this->m_datatypes[ECFDataType::SEPARATOR] = "divider";

	this->m_units[Unit::UnitLibrary::AMPERE] = "A";
	this->m_units[Unit::UnitLibrary::METRE] = "m";
	this->m_units[Unit::UnitLibrary::KILOGRAM] = "kg";
	this->m_units[Unit::UnitLibrary::SECOND] = "s";
	this->m_units[Unit::UnitLibrary::KELVIN] = "K";
	this->m_units[Unit::UnitLibrary::MOLE] = "mol";
	this->m_units[Unit::UnitLibrary::CANDELA] = "cd";
	this->m_units[Unit::UnitLibrary::PASCAL] = "Pa";
}
