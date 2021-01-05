
#ifndef SRC_INPUT_PARSERS_DISTRIBUTEDSCANNER_H_
#define SRC_INPUT_PARSERS_DISTRIBUTEDSCANNER_H_

#include "basis/logging/profiler.h"
#include "basis/containers/tarray.h"
#include "basis/io/inputfile.h"
#include "basis/utilities/packing.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/envinfo.h"

#include <string>
#include <vector>
#include <functional>

namespace espreso {

struct InputFile;

class DistributedScanner {
public:
	static bool check(const char *c, const std::string &key);
	static void align(InputFile &input, const std::string &chars);

	DistributedScanner()
	{
		memset(_filter, 0, 256);
		memset(_endFilter, 0, 256);
	}

	void add(const char* key, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip = [] (const char *c) { return 0; });
	void add(std::initializer_list<const char*> keys, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip = [] (const char *c) { return 0; });
	void addEnd(const char* key, std::function<void(const char *c)> constructor, std::function<size_t(const char *c)> skip = [] (const char *c) { return 0; });

	void scan(InputFile &input);
	void scanlines(InputFile &input);

	template<typename ...TArgs>
	void synchronize(TArgs& ...args)
	{
		profiler::syncstart("synchronize");
		_packedData.clear();
		_pack(0, args...);
		_clear(args...);
		_synchronize();
		const char *p = _packedData.data();
		while (p < _packedData.data() + _packedData.size()) {
			_unpack(0, p, args...);
		}
		profiler::syncend("synchronize");
	}

protected:
	struct Keyword {
		int offset, size;
		// accept a key and return minimal increment
		std::function<void(const char *c)> constructor;
		std::function<size_t(const char *c)> skip;
	};

	void _synchronize()
	{
		Communication::allGatherUnknownSize(_packedData);
	}

	template<typename TData>
	void _clear(TData &data)
	{
		data.clear();
	}

	template<typename TData, typename ...TOther>
	void _clear(TData &data, TOther& ...others)
	{
		_clear(data);
		_clear(others...);
	}

	template<typename TData>
	void _pack(int index, TData &data)
	{
		if (data.size()) {
			size_t size = utils::packedSize(index) + utils::packedSize(data);
			size_t offset = _packedData.size();
			_packedData.resize(offset + size);
			char *p = _packedData.data() + offset;
			utils::pack(index, p);
			utils::pack(data, p);
		}
	}

	template<typename TData, typename ...TOther>
	void _pack(int index, TData &data, TOther& ...others)
	{
		_pack(index++, data);
		_pack(index, others...);
	}

	template<typename TData>
	void _unpack(int index, const char* &p, TData &data)
	{
		int pindex;
		const char *pp = p;
		utils::unpack(pindex, pp);
		if (pindex == index) {
			TData tmp;
			utils::unpack(tmp, pp);
			data.insert(data.end(), tmp.begin(), tmp.end());
			p = pp;
		}
	}

	template<typename TData, typename ...TOther>
	void _unpack(int index, const char* &p, TData &data, TOther& ...others)
	{
		if (p < _packedData.data() + _packedData.size()) {
			_unpack(index++, p, data);
		}
		if (p < _packedData.data() + _packedData.size()) {
			_unpack(index, p, others...);
		}
	}

	std::string _keys;
	std::vector<Keyword> _keywords, _endKeywords;
	unsigned char _filter[256], _endFilter[256];
	std::vector<char> _packedData;
};

}

#endif /* SRC_INPUT_PARSERS_DISTRIBUTEDSCANNER_H_ */
