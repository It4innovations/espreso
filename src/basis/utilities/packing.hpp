
#include "packing.h"

namespace espreso {
namespace utils {

template<>
inline size_t packedSize(const std::string &data)
{
	return sizeof(size_t) + data.size();
}

template<typename Ttype>
inline size_t packedSize(const Ttype &data)
{
	return sizeof(Ttype);
}

template<>
inline size_t packedSize(const std::vector<std::string> &data)
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		size += packedSize(data[i]);
	}
	return size;
}

template<typename Ttype>
inline size_t packedSize(const std::vector<Ttype> &data)
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		size += packedSize(data[i]);
	}
	return size;
}

inline size_t packedSize(const VectorDense &data)
{
	return sizeof(data.size) + sizeof(double) * data.size;
}

template <typename TEBoundaries, typename TEData>
inline size_t packedSize(serializededata<TEBoundaries, TEData> *data)
{
	if (data != NULL) {
		return data->packedSize() + 1 + 3 * sizeof(size_t);
	}
	return 1;
}

template<>
inline void pack(const std::string &data, char* &p)
{
	size_t size = data.size();
	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	std::memcpy(p, data.data(), data.size());
	p += data.size();
}

inline void pack(VectorDense &data, char* &p)
{
	std::memcpy(p, &data.size, sizeof(data.size));
	p += sizeof(data.size);
	std::memcpy(p, data.vals, sizeof(double) * data.size);
	p += sizeof(double) * data.size;
}

template<typename Ttype>
inline void pack(const Ttype &data, char* &p)
{
	std::memcpy(p, &data, packedSize(data));
	p += packedSize(data);
}

template<>
inline void pack(const std::vector<std::string> &data, char* &p)
{
	size_t size = data.size();
	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	for (size_t i = 0; i < data.size(); i++) {
		pack(data[i], p);
	}
}

template<typename Ttype>
inline void pack(const std::vector<Ttype> &data, char* &p)
{
	size_t size = data.size();

	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	for (size_t i = 0; i < data.size(); i++) {
		pack(data[i], p);
	}
}

template <typename TEBoundaries, typename TEData>
inline void pack(serializededata<TEBoundaries, TEData> *data, char* &p)
{
	pack(data != NULL, p);
	if (data != NULL) {
		pack(data->threads(), p);
		pack(data->boundarytarray().size(), p);
		pack(data->datatarray().size(), p);
		data->pack(p);
	}
}

template<>
inline void unpack(std::string &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);
	data = std::string(p, size);
	p += size;
}

inline void unpack(VectorDense &data, const char* &p)
{
	std::memcpy(&data.size, p, sizeof(data.size));
	p += packedSize(data.size);
	data.resize(data.size);
	std::memcpy(data.vals, p, sizeof(double) * data.size);
	p += sizeof(double) * data.size;
}

template<typename Ttype>
inline void unpack(Ttype &data, const char* &p)
{
	std::memcpy(reinterpret_cast<void*>(&data), p, packedSize(data));
	p += packedSize(data);
}

template<>
inline void unpack(std::vector<std::string> &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	if (size) {
		data.resize(size);
		for (size_t i = 0; i < data.size(); i++) {
			unpack(data[i], p);
		}
	}
}

template<typename Ttype>
inline void unpack(std::vector<Ttype> &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	data.resize(size);
	for (size_t i = 0; i < data.size(); i++) {
		unpack(data[i], p);
	}
}

template <typename TEBoundaries, typename TEData>
inline void unpack(serializededata<TEBoundaries, TEData> *&data, const char* &p)
{
	if (data != NULL) {
		delete data;
		data = NULL;
	}

	bool notnull;
	unpack(notnull, p);
	if (notnull) {
		size_t threads, bsize, dsize;
		unpack(threads, p);
		unpack(bsize, p);
		unpack(dsize, p);
		data = new serializededata<TEBoundaries, TEData>(tarray<TEBoundaries>(threads, bsize), tarray<TEData>(threads, dsize));
		data->unpack(p);
	}
}

}
}
