
#ifndef SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_
#define SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_

#include "point.h"
#include "serializededata.h"

#include <functional>

namespace espreso {

/*
 * header:
 * [0 - 2] voxel difference, 000 for duplicated voxel with DIFF(1, 0, 0), 111 for FULL voxel
 * [3] - duplicated value
 */
struct VolumePacker {
	_Point<short> size;

	inline int bits(const short &diff) const
	{
		switch (diff) {
		case -31 ... -16: return 6;
		case -15 ...  -8: return 5;
		case  -7 ...  -4: return 4;
		case  -3 ...  -2: return 3;
		case          -1: return 2;
		case   0 ...   1: return 1;
		case           2: return 2;
		case   3 ...   4: return 3;
		case   5 ...   8: return 4;
		case   9 ...  16: return 5;
		case  17 ...  32: return 6;
		default: return 7; // full voxels will be stored
		}
	}

	inline int bits(const _Point<short> &diff) const
	{
		return std::max(bits(diff.x), std::max(bits(diff.y), bits(diff.z)));
	}

	VolumePacker(const _Point<short> &grid)
	{
		while (grid.x > (1 << ++size.x));
		while (grid.y > (1 << ++size.y));
		while (grid.z > (1 << ++size.z));
	}

	size_t analyze(serializededata<size_t, _Point<short> > *volumeIndices, size_t begin, size_t end) const
	{
		size_t bytes = 0;
		_Point<short> prev, diff;
		for (auto e = volumeIndices->cbegin() + begin; e != volumeIndices->cbegin() + end; ++e) {
			bool duplicate = false;
			for (auto v = e->begin(); v != e->end(); prev = *v++) {
				diff = *v - prev;
				if (duplicate && diff.x == 1 && diff.y == 0 && diff.z == 0) {
					bytes += 3;
				} else {
					int max = bits(diff);
					if (max == 7) {
						bytes += 4 + size.x + size.y + size.z;
					} else {
						bytes += 4 + 3 * max;
					}
				}
				duplicate = true;
			}
		}
		bytes = sizeof(size_t) + bytes / 8 + 1;
		return bytes;
	}

	void pack(serializededata<size_t, _Point<short> > *volumeIndices, size_t begin, size_t end, char *packed) const
	{
		struct {
			char *data; int empty;

			void insert(const char *value, int bits)
			{
				*data &= 255 << empty;
				if (bits <= empty) {
					_insert(*value, bits);
				} else {
					int half = bits - empty;
					_insert(*value >> half, empty);
					_insert(*value, half);
				}
			}
		private:
			void _insert(const char value, int bits)
			{
				empty -= bits;
				*data |= (value & (255 >> (8 - bits))) << empty;
				if (empty == 0) {
					data += 1;
					*data = 0;
					empty = 8;
				}
			}
		} stream{packed, 8};

		size_t count = (volumeIndices->cbegin() + end)->begin() - (volumeIndices->cbegin() + begin)->begin();
		for (size_t i = 0; i < sizeof(count); ++i) {
			stream.insert(reinterpret_cast<const char*>(&count) + i, 8);
		}

		_Point<short> prev, diff;
		for (auto e = volumeIndices->cbegin() + begin; e != volumeIndices->cbegin() + end; ++e) {
			int duplicate = 0;
			for (auto v = e->begin(); v != e->end(); prev = *v++) {
				diff = *v - prev;
				if (duplicate && diff.x == 1 && diff.y == 0 && diff.z == 0) {
					char header = 0;
					stream.insert(&header, 3);
				} else {
					char max = bits(diff);
					char header = (max << 1) | duplicate;
					stream.insert(&header, 4);
					if (max == 7) {
						if (size.x > 8) {
							stream.insert(reinterpret_cast<const char*>(&v->x), 8);
							stream.insert(reinterpret_cast<const char*>(&v->x) + 1, size.x - 8);
						} else {
							stream.insert(reinterpret_cast<const char*>(&v->x), size.x);
						}
						if (size.y > 8) {
							stream.insert(reinterpret_cast<const char*>(&v->y), 8);
							stream.insert(reinterpret_cast<const char*>(&v->y) + 1, size.y - 8);
						} else {
							stream.insert(reinterpret_cast<const char*>(&v->y), size.y);
						}
						if (size.z > 8) {
							stream.insert(reinterpret_cast<const char*>(&v->z), 8);
							stream.insert(reinterpret_cast<const char*>(&v->z) + 1, size.z - 8);
						} else {
							stream.insert(reinterpret_cast<const char*>(&v->z), size.z);
						}
					} else {
						if (max > 1) {
							diff += _Point<short>((1 << (max - 1)) - 1);
						}
						stream.insert(reinterpret_cast<const char*>(&diff.x), max);
						stream.insert(reinterpret_cast<const char*>(&diff.y), max);
						stream.insert(reinterpret_cast<const char*>(&diff.z), max);
					}
				}
				duplicate = 1;
			}
		}
	}

	void unpack(char *packed, std::function<void(const _Point<short>&, const size_t&)> callback) const
	{
		struct {
			char *data; int unread;

			char extract(int bits)
			{
				char value;
				extract(&value, bits);
				return value;
			}

			void extract(char *value, int bits)
			{
				*value = 0;
				if (bits <= unread) {
					*value = _extract(*data, bits);
				} else {
					int half = bits - unread;
					*value = _extract(*data, unread) << half;
					*value |= _extract(*data, half);
				}
			}
		private:
			char _extract(const char value, int bits)
			{
				unread -= bits;
				char result = (value >> unread) & (255 >> (8 - bits));
				if (unread == 0) {
					data += 1;
					unread = 8;
				}
				return result;
			}
		} stream{packed, 8};

		size_t count;
		for (size_t i = 0; i < sizeof(count); ++i) {
			stream.extract(reinterpret_cast<char*>(&count) + i, 8);
		}

		_Point<short> voxel;
		for (size_t v = 0, i = 0; v < count; ++v) {
			_Point<short> diff;
			char header = stream.extract(3);
			char duplicate = 0;
			switch (header) {
			case 0:
				duplicate = 1;
				diff = _Point<short>(1, 0, 0);
				voxel += diff;
				break;
			case 7:
				duplicate = stream.extract(1);
				if (size.x > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.x), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.x) + 1, size.x - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.x), size.x);
				}
				if (size.y > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.y), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.y) + 1, size.y - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.y), size.y);
				}
				if (size.z > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.z), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.z) + 1, size.z - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.z), size.z);
				}
				break;
			default:
				duplicate = stream.extract(1);
				stream.extract(reinterpret_cast<char*>(&diff.x), header);
				stream.extract(reinterpret_cast<char*>(&diff.y), header);
				stream.extract(reinterpret_cast<char*>(&diff.z), header);
				if (header > 1) {
					diff -= _Point<short>((1 << (header - 1)) - 1);
				}
				voxel += diff;
			}
			if (v && !duplicate) ++i;
			callback(voxel, i);
		}
	}
};

}



#endif /* SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_ */
