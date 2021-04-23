
#ifndef SRC_BASIS_CONTAINERS_PACKER_H_
#define SRC_BASIS_CONTAINERS_PACKER_H_

#include "basis/utilities/packing.h"

#include <cstddef>
#include <vector>

namespace espreso {

struct PackerBase {
	virtual size_t packedSize() { return 0; }
	virtual void pack(char* &p) { }
	virtual void unpack(const char* &p) { }

	virtual PackerBase* copy() const { return new PackerBase(*this); }
	virtual ~PackerBase() {}
};

template <typename... Other> struct PackerData: PackerBase  { };

template <typename T, typename... Other>
struct PackerData<T, Other...>: PackerBase {

	PackerData(T& data, Other& ... other): data(data), other(other...) {}
	PackerBase* copy() const { return new PackerData<T, Other...>(*this); }

	size_t packedSize() override
	{
		return utils::packedSize(data) + other.packedSize();
		return other.packedSize();
	}

	void pack(char* &p) override
	{
		utils::pack(data, p);
		other.pack(p);
	}

	void unpack(const char* &p) override
	{
		utils::unpack(data, p);
		other.unpack(p);
	}

	T &data;
	PackerData<Other...> other;
};

class Packer {
	char *_buffer;
	size_t _size;
	std::vector<PackerBase*> _data;

public:
	template <typename... Data>
	Packer(Data& ... data): _buffer(NULL), _size(0), _data{new PackerData<Data...>(data...)} {}
	~Packer()
	{
		if (_buffer) delete[] _buffer;
		for (size_t i = 0; i < _data.size(); ++i) {
			delete _data[i];
		}
	}

	Packer& operator+=(const Packer &other)
	{
		for (size_t i = 0; i < other._data.size(); ++i) {
			_data.push_back(other._data[i]->copy());
		}
		return *this;
	}

	size_t size() const
	{
		return _size;
	}

	const char* pack()
	{
		_size = 0;
		for (size_t i = 0; i < _data.size(); ++i) {
			_size += _data[i]->packedSize();
		}

		if (_buffer) delete[] _buffer;
		char *p = _buffer = new char[_size];

		for (size_t i = 0; i < _data.size(); ++i) {
			_data[i]->pack(p);
		}
		return _buffer;
	}

	void unpack() const
	{
		const char *p = _buffer;
		for (size_t i = 0; i < _data.size(); ++i) {
			_data[i]->unpack(p);
		}
	}
};

} // namespace espreso

#endif /* SRC_BASIS_CONTAINERS_PACKER_H_ */
