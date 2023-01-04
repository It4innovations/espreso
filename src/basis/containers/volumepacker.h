
#ifndef SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_
#define SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_

#include "point.h"
#include "serializededata.h"
#include "mesh/store/elementstore.h"

#include <functional>

namespace espreso {

/*
 * header:
 * [0 - 2] voxel difference, 000 for voxel with DIFF(1, 0, 0), 111 for FULL voxel
 * [3] - duplicated value
 */
struct VolumePacker {
	_Point<short> gridSize;
	size_t packSize, packedVoxels, variableSize;
	std::vector<std::string> names;
	char *packedData;

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
	: packSize(0), packedVoxels(0), variableSize(0), packedData(nullptr)
	{
		while (grid.x > (1 << ++gridSize.x));
		while (grid.y > (1 << ++gridSize.y));
		while (grid.z > (1 << ++gridSize.z));
	}

	~VolumePacker()
	{
		if (packedData) delete[] packedData;
	}

	void allocate()
	{
		packedData = new char[packSize];
	}

	void analyze(serializededata<size_t, _Point<short> > *volumeIndices, const std::vector<ElementData*> &variables)
	{
		size_t bytes = 0, i = 0;
		_Point<short> prev, diff;
		for (auto e = volumeIndices->cbegin(); e != volumeIndices->cend(); ++e, ++i) {
			if (e->size() == 0) continue;
			bool duplicate = false, filtered = false;
			for (size_t v = 0; v < variables.size(); ++v) {
				if (variables[v]->filter.size()) {
					filtered = true;
					for (size_t f = 0; f < variables[v]->filter.size(); ++f) {
						if (variables[v]->filter[f].first <= variables[v]->store[i] && variables[v]->store[i] <= variables[v]->filter[f].second) {
							filtered = false;
							break;
						}
					}
				}
			}
			if (!filtered) {
				packedVoxels += e->size();
				for (auto v = e->begin(); v != e->end(); prev = *v++) {
					diff = *v - prev;
					if (duplicate && diff.x == 1 && diff.y == 0 && diff.z == 0) {
						bytes += 3;
					} else {
						int max = bits(diff);
						if (max == 7) {
							bytes += 4 + gridSize.x + gridSize.y + gridSize.z;
						} else {
							bytes += 4 + 3 * max;
						}
					}
					if (!duplicate) {
						bytes += 8 * sizeof(float) * variables.size();
					}
					duplicate = true;
				}
			}
		}
		packSize = sizeof(size_t) + bytes / 8 + 1;
		variableSize = variables.size();
		for (size_t v = 0; v < variables.size(); ++v) {
			names.push_back(variables[v]->name);
		}
	}

	void pack(serializededata<size_t, _Point<short> > *volumeIndices, const std::vector<ElementData*> &variables)
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
		} stream{packedData, 8};

		for (size_t i = 0; i < sizeof(packedVoxels); ++i) {
			stream.insert(reinterpret_cast<const char*>(&packedVoxels) + i, 8);
		}

		size_t i = 0;
		_Point<short> prev, diff;
		for (auto e = volumeIndices->cbegin(); e != volumeIndices->cend(); ++e, ++i) {
			if (e->size() == 0) continue;
			bool duplicate = false, filtered = false;
			for (size_t v = 0; v < variables.size(); ++v) {
				if (variables[v]->filter.size()) {
					filtered = true;
					for (size_t f = 0; f < variables[v]->filter.size(); ++f) {
						if (variables[v]->filter[f].first <= variables[v]->store[i] && variables[v]->store[i] <= variables[v]->filter[f].second) {
							filtered = false;
							break;
						}
					}
				}
			}
			if (!filtered) {
				for (auto voxels = e->begin(); voxels != e->end(); prev = *voxels++) {
					diff = *voxels - prev;
					if (duplicate && diff.x == 1 && diff.y == 0 && diff.z == 0) {
						char header = 0;
						stream.insert(&header, 3);
					} else {
						char max = bits(diff);
						char header = (max << 1) | duplicate;
						stream.insert(&header, 4);
						if (max == 7) {
							if (gridSize.x > 8) {
								stream.insert(reinterpret_cast<const char*>(&voxels->x), 8);
								stream.insert(reinterpret_cast<const char*>(&voxels->x) + 1, gridSize.x - 8);
							} else {
								stream.insert(reinterpret_cast<const char*>(&voxels->x), gridSize.x);
							}
							if (gridSize.y > 8) {
								stream.insert(reinterpret_cast<const char*>(&voxels->y), 8);
								stream.insert(reinterpret_cast<const char*>(&voxels->y) + 1, gridSize.y - 8);
							} else {
								stream.insert(reinterpret_cast<const char*>(&voxels->y), gridSize.y);
							}
							if (gridSize.z > 8) {
								stream.insert(reinterpret_cast<const char*>(&voxels->z), 8);
								stream.insert(reinterpret_cast<const char*>(&voxels->z) + 1, gridSize.z - 8);
							} else {
								stream.insert(reinterpret_cast<const char*>(&voxels->z), gridSize.z);
							}
						} else {
							if (max > 1) {
								diff += _Point<short>((1 << (max - 1)) - 1); // diff is always non-negative
							}
							stream.insert(reinterpret_cast<const char*>(&diff.x), max);
							stream.insert(reinterpret_cast<const char*>(&diff.y), max);
							stream.insert(reinterpret_cast<const char*>(&diff.z), max);
						}
						if (!duplicate) {
							for (size_t v = 0; v < variables.size(); ++v) {
								float value = 0;
								if (variables[v]->dimension == 1) {
									value = variables[v]->store[i];
								} else {
									for (int d = 0; d < variables[v]->dimension; ++d) {
										value += variables[v]->store[variables[v]->dimension * i + d] * variables[v]->store[variables[v]->dimension * i + d];
									}
									value = std::sqrt(value);
								}
								for (size_t b = 0; b < sizeof(value); ++b) {
									stream.insert(reinterpret_cast<const char*>(&value) + b, 8);
								}
							}
						}
					}
					duplicate = 1;
				}
			}
		}
	}

	void unpack(std::function<void(const _Point<short>&, const size_t&, const float&)> callback) const
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
		} stream{packedData, 8};

		size_t count;
		for (size_t i = 0; i < sizeof(count); ++i) {
			stream.extract(reinterpret_cast<char*>(&count) + i, 8);
		}

		std::vector<float> values(variableSize);
		_Point<short> voxel;
		for (size_t i = 0; i < count; ++i) {
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
				if (gridSize.x > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.x), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.x) + 1, gridSize.x - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.x), gridSize.x);
				}
				if (gridSize.y > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.y), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.y) + 1, gridSize.y - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.y), gridSize.y);
				}
				if (gridSize.z > 8) {
					stream.extract(reinterpret_cast<char*>(&voxel.z), 8);
					stream.extract(reinterpret_cast<char*>(&voxel.z) + 1, gridSize.z - 8);
				} else {
					stream.extract(reinterpret_cast<char*>(&voxel.z), gridSize.z);
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
			if (!duplicate) {
				for (size_t v = 0; v < variableSize; ++v) {
					for (size_t b = 0; b < sizeof(float); ++b) {
						stream.extract(reinterpret_cast<char*>(&values[v]) + b, 8);
					}
				}
			}

			for (size_t v = 0; v < variableSize; ++v) {
				callback(voxel, v, values[v]);
			}
		}
	}
};

}



#endif /* SRC_BASIS_CONTAINERS_VOLUMEPACKER_H_ */
