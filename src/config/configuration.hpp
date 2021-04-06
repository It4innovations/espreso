
#ifndef SRC_CONFIGURATION_CONFIGURATION_HPP_
#define SRC_CONFIGURATION_CONFIGURATION_HPP_

#include "holders/objectholder.h"
#include "holders/valueholder.h"
#include "holders/vectorholder.h"
#include <type_traits>

#include "configuration.h"

#define REGISTER(parameter, metadata, ...) ecfdescription->registerParameter(#parameter, parameter, metadata, ##__VA_ARGS__)
#define PNAME(parameter) #parameter

namespace espreso {

/////////// PARAMETER ///////////
/////////////////////////////////

// Child of ECFDescription
template <typename Ttype>
typename std::enable_if<std::is_class<Ttype>::value && std::is_base_of<ECFDescription, Ttype>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 0);
	metadata.checkpattern(name, 0);

	parameters.push_back(parameter.ecfdescription);
	parameters.back()->name = name;
	parameters.back()->metadata = metadata;
	parameters.back()->defaultName();
	return parameters.back();
}

// ENUM
template<typename Ttype>
typename std::enable_if<std::is_enum<Ttype>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 0);

	return registerParameter(name, new ECFEnumHolder<Ttype>(parameter), metadata);
}

// REST
template <typename Ttype>
typename std::enable_if<(!std::is_class<Ttype>::value && !std::is_enum<Ttype>::value) || (std::is_class<Ttype>::value && !std::is_base_of<ECFDescription, Ttype>::value), ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, Ttype &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 0);

	return registerParameter(name, new ECFValueHolder<Ttype>(parameter), metadata);
}

//////////// VECTOR /////////////
/////////////////////////////////

template<typename Ttype>
ECFParameter*
ECFObject::registerParameter(const std::string &name, std::vector<Ttype> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 1);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 0);

	return registerParameter(name, new ECFVectorHolder<Ttype>(parameter), metadata);
}

////////////// MAP //////////////
/////////////////////////////////

// TYPE2 = Child of ECFDescription
template<typename Ttype1, typename Ttype2, typename... TArgs>
typename std::enable_if<std::is_class<Ttype2>::value && std::is_base_of<ECFDescription, Ttype2>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 1);
	metadata.checkpattern(name, 1);

	return registerParameter(name, new ECFObjectMap<Ttype1, Ttype2, TArgs...>(parameter, args...), metadata);
}

// TYPE2 = ENUM
template<typename Ttype1, typename Ttype2>
typename std::enable_if<std::is_enum<Ttype2>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFEnumMap<Ttype1, Ttype2>(parameter), metadata);
}

// TYPE2 = REST
template<typename Ttype1, typename Ttype2>
typename std::enable_if<(!std::is_class<Ttype2>::value && !std::is_enum<Ttype2>::value) || (std::is_class<Ttype2>::value && !std::is_base_of<ECFDescription, Ttype2>::value), ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFValueMap<Ttype1, Ttype2>(parameter), metadata);
}

template<typename Ttype1, typename Ttype2, typename... TArgs>
typename std::enable_if<(!std::is_class<Ttype2>::value && !std::is_enum<Ttype2>::value) || (std::is_class<Ttype2>::value && !std::is_base_of<ECFDescription, Ttype2>::value), ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, Ttype2> &parameter, const ECFMetaData &metadata, TArgs... args)
{
	metadata.checkdescription(name, 2);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFValueMap<Ttype1, Ttype2, TArgs...>(parameter, args...), metadata);
}

//////////// MAP MAP ////////////
/////////////////////////////////

// TYPE3 = Child of ECFDescription
template<typename Ttype1, typename Ttype2, typename Ttype3, typename... TArgs>
typename std::enable_if<std::is_class<Ttype3>::value && std::is_base_of<ECFDescription, Ttype3>::value, ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata, TArgs... args)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 2);
	metadata.checkpattern(name, 2);

	return registerParameter(name, new ECFObjectMapMap<Ttype1, Ttype2, Ttype3, TArgs...>(parameter, args...), metadata);
}

// TYPE3 = REST
template<typename Ttype1, typename Ttype2, typename Ttype3>
typename std::enable_if<!std::is_class<Ttype3>::value || (std::is_class<Ttype3>::value && !std::is_base_of<ECFDescription, Ttype3>::value), ECFParameter*>::type
ECFObject::registerParameter(const std::string &name, std::map<Ttype1, std::map<Ttype2, Ttype3> > &parameter, const ECFMetaData &metadata)
{
	metadata.checkdescription(name, 3);
	metadata.checkdatatype(name, 3);
	metadata.checkpattern(name, 3);

	return registerParameter(name, new ECFValueMapMap<Ttype1, Ttype2, Ttype3>(parameter), metadata);
}

}


#endif /* SRC_CONFIGURATION_CONFIGURATION_HPP_ */
