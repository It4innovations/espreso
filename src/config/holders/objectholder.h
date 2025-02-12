
#ifndef SRC_CONFIG_HOLDERS_OBJECTHOLDER_H_
#define SRC_CONFIG_HOLDERS_OBJECTHOLDER_H_

#include "valueholder.h"
#include "config/configuration.h"

namespace espreso {

template <typename TParameter, typename TValue, typename... TArgs>
struct ECFValueMap: public ECFObject {
    std::map<TParameter, TValue> &value;
    std::tuple<TArgs...> args;

    ECFValueMap(std::map<TParameter, TValue> &value, TArgs... args): value(value), args(args...) {}
    ECFValueMap(std::map<TParameter, TValue> &value, std::tuple<TArgs...> &args): value(value), args(args) {}

    virtual ECFParameter* _getNamedParameter(const std::string &name)
    {
        if (ECFObject::_getNamedParameter(name) != NULL) {
            return ECFObject::_getNamedParameter(name);
        }
        TParameter key;
        if (ECFValueHolder<TParameter>(key).setValue(name)) {
            value.emplace(std::piecewise_construct, std::forward_as_tuple(key), args);
            return registerParameter(name, new ECFValueHolder<TValue>(value.at(key)), metadata.suffix(1));
        } else {
            return NULL;
        }
    }

    virtual void dropParameter(ECFParameter *parameter)
    {
        TParameter key;
        ECFValueHolder<TParameter>(key).setValue(parameter->name);
        value.erase(key);
        ECFObject::dropParameter(parameter);
    }

    virtual ECFParameter* getPattern() const
    {
        TParameter key{};
        if (_pattern.size() == 0) {
            _pattern.emplace(std::piecewise_construct, std::forward_as_tuple(key), args);
        }
        ECFParameter *parameter = new ECFValueHolder<TValue>(_pattern.at(key));
        parameter->setValue(metadata.pattern[1]);
        return registerPatternParameter(parameter);
    }

    virtual const void* data() const { return &value; }

private:
    mutable std::map<TParameter, TValue> _pattern;
};

template <typename TParameter, typename TValue>
struct ECFEnumMap: public ECFObject {
    std::map<TParameter, TValue> &value;

    ECFEnumMap(std::map<TParameter, TValue> &value): value(value) {}

    virtual ECFParameter* _getNamedParameter(const std::string &name)
    {
        if (ECFObject::_getNamedParameter(name) != NULL) {
            return ECFObject::_getNamedParameter(name);
        }
        TParameter key;
        if (ECFValueHolder<TParameter>(key).setValue(name)) {
            return registerParameter(name, new ECFEnumHolder<TValue>(value[key]), metadata.suffix(1));
        } else {
            return NULL;
        }
    }

    virtual void dropParameter(ECFParameter *parameter)
    {
        TParameter key;
        ECFValueHolder<TParameter>(key).setValue(parameter->name);
        value.erase(key);
        ECFObject::dropParameter(parameter);
    }

    virtual ECFParameter* getPattern() const
    {
        ECFParameter *parameter = new ECFEnumHolder<TValue>(_patternValue);
        parameter->metadata = metadata.suffix(1);
        parameter->setValue(metadata.pattern[1]);
        return registerPatternParameter(parameter);
    }

    virtual const void* data() const { return &value; }

private:
    mutable TValue _patternValue;
};

template <typename TParameter, typename TObject, typename... TArgs>
struct ECFObjectMap: public ECFObject {
    std::map<TParameter, TObject> &value;
    std::tuple<TArgs...> args;

    ECFObjectMap(std::map<TParameter, TObject> &value, TArgs... args): value(value), args(args...) {}
    ECFObjectMap(std::map<TParameter, TObject> &value, std::tuple<TArgs...> &args): value(value), args(args) {}

    virtual ECFParameter* _getNamedParameter(const std::string &name)
    {
        if (ECFObject::_getNamedParameter(name) != NULL) {
            return ECFObject::_getNamedParameter(name);
        }
        TParameter key;
        if (ECFValueHolder<TParameter>(key).setValue(name)) {
            auto it = value.emplace(std::piecewise_construct, std::forward_as_tuple(key), args);
            parameters.push_back(it.first->second.ecfdescription);
            parameters.back()->name = name;
            parameters.back()->metadata = metadata.suffix(1);
            return parameters.back();
        } else {
            return NULL;
        }
    }

    virtual void dropParameter(ECFParameter *parameter)
    {
        TParameter key;
        ECFValueHolder<TParameter>(key).setValue(parameter->name);
        value.erase(key);
        ECFObject::dropParameter(parameter);
    }

    virtual void dropAllParameters()
    {
        value.clear();
        ECFObject::dropAllParameters();
    }

    virtual ECFParameter* getPattern() const
    {
        TParameter key{};
        if (_pattern.size() == 0) {
            _pattern.emplace(std::piecewise_construct, std::forward_as_tuple(key), args);
        }
        ECFParameter *parameter = _pattern.at(key).ecfdescription;
        parameter->metadata = metadata.suffix(1);
        return parameter;
    }

    virtual const void* data() const { return &value; }

private:
    mutable std::map<TParameter, TObject> _pattern;
};

template <typename TParameter1, typename TParameter2, typename TValue>
struct ECFValueMapMap: public ECFObject {
    std::map<TParameter1, std::map<TParameter2, TValue> > &value;

    ECFValueMapMap(std::map<TParameter1, std::map<TParameter2, TValue> > &value): value(value) {}

    virtual ECFParameter* _getNamedParameter(const std::string &name)
    {
        if (ECFObject::_getNamedParameter(name) != NULL) {
            return ECFObject::_getNamedParameter(name);
        }
        TParameter1 key;
        if (ECFValueHolder<TParameter1>(key).setValue(name)) {
            return registerParameter(name, new ECFValueMap<TParameter2, TValue>(value[key]), metadata.suffix(1));
        } else {
            return NULL;
        }
    }

    virtual void dropParameter(ECFParameter *parameter)
    {
        TParameter1 key;
        ECFValueHolder<TParameter1>(key).setValue(parameter->name);
        value.erase(key);
        ECFObject::dropParameter(parameter);
    }

    virtual void dropAllParameters()
    {
        value.clear();
        ECFObject::dropAllParameters();
    }

    virtual ECFParameter* getPattern() const
    {
        return registerPatternParameter(new ECFValueMap<TParameter2, TValue>(_patternValue));
    }

    virtual const void* data() const { return &value; }

private:
    mutable std::map<TParameter2, TValue> _patternValue;
};


template <typename TParameter1, typename TParameter2, typename TObject, typename... TArgs>
struct ECFObjectMapMap: public ECFObject {
    std::map<TParameter1, std::map<TParameter2, TObject> > &value;
    std::tuple<TArgs...> args;

    ECFObjectMapMap(std::map<TParameter1, std::map<TParameter2, TObject> > &value, TArgs... args): value(value), args(args...) {}

    virtual ECFParameter* _getNamedParameter(const std::string &name)
    {
        if (ECFObject::_getNamedParameter(name) != NULL) {
            return ECFObject::_getNamedParameter(name);
        }
        TParameter1 key;
        if (ECFValueHolder<TParameter1>(key).setValue(name)) {
            return registerParameter(name, new ECFObjectMap<TParameter2, TObject, TArgs...>(value[key], args), metadata.suffix(1));
        } else {
            return NULL;
        }
    }

    virtual void dropParameter(ECFParameter *parameter)
    {
        TParameter1 key;
        ECFValueHolder<TParameter1>(key).setValue(parameter->name);
        value.erase(key);
        ECFObject::dropParameter(parameter);
    }

    virtual void dropAllParameters()
    {
        value.clear();
        ECFObject::dropAllParameters();
    }

    virtual ECFParameter* getPattern() const
    {
        return registerPatternParameter(new ECFObjectMap<TParameter2, TObject>(_patternValue));
    }

    virtual const void* data() const { return &value; }

private:
    mutable std::map<TParameter2, TObject> _patternValue;
};

}

#endif /* SRC_CONFIG_HOLDERS_OBJECTHOLDER_H_ */
