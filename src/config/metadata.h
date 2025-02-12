
#ifndef SRC_CONFIG_METADATA_H_
#define SRC_CONFIG_METADATA_H_

#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct TensorConfiguration;
struct ECFAbstractCondition;

enum struct ECFDataType {
    BOOL,
    STRING,
    RANGE,
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
    SEPARATOR
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

struct Unit {

    enum class UnitLibrary {
        METRE,
        KILOGRAM,
        SECOND,
        AMPERE,
        KELVIN,
        MOLE,
        CANDELA,
        PASCAL
    };

    struct UnitElement {
        UnitLibrary unit;
        int exponent;
    };

    std::vector<UnitElement> elements;

    Unit() {}

    Unit& add(UnitLibrary unit, int exponent);
};

enum struct ECFGUIType {
    GROUP_BLOCK,
    GROUP_COLLAPSED,
    STATIC,
    DYNAMIC,
    FORM
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
    Unit unit;
    bool exporting;
    int range_begin;
    int range_end;
    ECFGUIType gui_type;
    std::string pattern_name;
    std::string pattern_item_name;

    ECFAbstractCondition *condition;

    std::function<bool(void)> isallowed;
    std::function<bool(void)> ismandatory;

    ECFMetaData& setname(const std::string &name) { this->name = name; return *this; }
    ECFMetaData& setdescription(const std::vector<std::string> &description) { this->description = description; return *this; }
    ECFMetaData& setdatatype(const std::vector<ECFDataType> &datatype) { this->datatype = datatype; return *this; }
    ECFMetaData& setpattern(const std::vector<std::string> &pattern) { this->pattern = pattern; return *this; }
    ECFMetaData& setvariables(const std::vector<std::string> &variables) { this->variables = variables; return *this; }
    ECFMetaData& setrange(int begin, int end) { this->range_begin = begin; this->range_end = end; return *this; }
    ECFMetaData& settensor(TensorConfiguration &tensor) { this->tensor = &tensor; return *this; }
    ECFMetaData& settensor(TensorConfiguration &tensor, int row, int column) 
    { this->tensor = &tensor; this->tensor_row = row; this->tensor_column = column; return *this; }
    ECFMetaData& registertensor(TensorConfiguration& tensor)
    { this->tensors.push_back(&tensor); return *this; }
    ECFMetaData& setunit(const Unit &unit) { this->unit = unit; return *this; }
    ECFMetaData& setcollapsed() { this->gui_type = ECFGUIType::GROUP_COLLAPSED; return *this; }
    ECFMetaData& displayobjectname() { this->visibleObjectName = true; return *this; }
    ECFMetaData& noexport() { this->exporting = false; return *this; }
    ECFMetaData& allowonly(std::function<bool(void)> isallowed) { this->isallowed = isallowed; return *this; }
    ECFMetaData& mandatoryonly(std::function<bool(void)> ismandatory) { this->ismandatory = ismandatory; return *this; }
    ECFMetaData& addconstraint(const ECFAbstractCondition &condition);
    ECFMetaData& removeconstraint();
    ECFMetaData& setdynamic() { this->gui_type = ECFGUIType::DYNAMIC; return *this; }
    ECFMetaData& setgroup() { this->gui_type = ECFGUIType::GROUP_BLOCK; return *this; }
    ECFMetaData& setform() { this->gui_type = ECFGUIType::FORM; return *this; }
    ECFMetaData& setpatternname(const std::string &name) { this->pattern_name = name; return *this; }
    ECFMetaData& setpatternitemname(const std::string &name) { this->pattern_item_name = name; return *this; }


    ECFMetaData& addoption(const ECFOption &option) { options.push_back(option); return *this; }

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
