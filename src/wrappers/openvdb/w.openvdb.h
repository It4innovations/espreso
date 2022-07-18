
#ifndef SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_
#define SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_

namespace espreso {

struct OpenVDBWrapperData;

template <typename TType> class _Point;
template <typename TEBoundaries, typename TEData> class serializededata;
struct ElementData;

struct OpenVDBWrapper {

	OpenVDBWrapper();
	~OpenVDBWrapper();

	void add_grid(const serializededata<esint, _Point<int> > *emap, const ElementData *data);
	void store_grids(const char *name);

protected:
	OpenVDBWrapperData *_data;
};

}

#endif /* SRC_WRAPPERS_OPENVDB_W_OPENVDB_H_ */
