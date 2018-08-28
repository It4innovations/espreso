
#include "config/configuration.h"

#include <ostream>

namespace espreso
{
class ECFJSONExport
{
public:
	ECFJSONExport(ECFObject* obj,
		std::ostream& outputStream);
	ECFJSONExport(ECFValue* value,
		ECFObject* parent,
		std::ostream& outputStream);

	void exportToStream();

protected:
	std::ostream& m_stream;
	ECFParameter* m_rootValue;
	ECFObject* m_rootObject;
	ECFObject* m_parent;

	std::map<ECFDataType, std::string> m_datatypes;
	std::map<Unit::UnitLibrary, std::string> m_units;

	int m_tensorId;
	std::map<TensorConfiguration*, int> m_tensors;

	void init();

	void printParameter(ECFParameter* param,
		int * const id,
		const int indent = 0,
		const bool printKey = true);
	void printObjectContent(ECFObject* object,
		int * const id,
		const int indent = 0,
		const bool printKey = true);
	void printObjectGroup(ECFObject* object,
		int * const id,
		const int indent = 0,
		const bool printKey = true);
	void printObjectInTree(ECFObject* object,
		int * const id,
		const int indent = 0,
		const bool printKey = true);
	std::string spaces(const int indent = 0);
	void printKeyValuePair(const std::string& key,
		const std::string& val,
		const int indent = 0,
		const char * valueEnclosed = "\"");
	void printKeyValuePair(const std::string& key,
		const int val,
		const int indent = 0,
		const char * valueEnclosed = "\"");
	void printKeyArrayPair(const std::string& key,
		const std::vector<std::string>& arr,
		const int indent = 0,
		const char* itemEnclosed = "\"");
	void printKeyObjectPair(const std::string& key,
		const std::function<void(int)>& printObjContent,
		const int indent = 0);
	void printMetaData(ECFParameter* p, const int indent = 0);
	void printValue(ECFParameter* val, const int indent = 0);
	void printUnit(ECFParameter* p, const int indent = 0);
	void printConstraint(ECFParameter* p, int indent = 0);
	void printRegisteredTensors(ECFParameter* p, int indent = 0);
	void printTensor(ECFParameter* p, int indent = 0);
};
}
