
#ifndef SRC_WRAPPERS_HDF5_W_HDF5_H_
#define SRC_WRAPPERS_HDF5_W_HDF5_H_

namespace espreso {

struct MPIGroup;

class HDF5 {
	class H5File;
	struct H5Type;
	struct H5TypeWrapper { void *data; };

public:
	enum class MODE { READ, WRITE };

	static bool islinked();

	static H5TypeWrapper INT;
	static H5TypeWrapper LONG;
	static H5TypeWrapper ESINT;
	static H5TypeWrapper FLOAT;
	static H5TypeWrapper DOUBLE;

	HDF5(const char* file, MPIGroup &mpigroup, MODE mode);
	~HDF5();

	void append(
			const char* name, const H5TypeWrapper &source, const H5TypeWrapper &target,
			const void *data, esint esize, esint nelements, esint offset, esint totalsize);

	void read(
			const char* name, const H5TypeWrapper &target,
			void *data, esint esize, esint nelements, esint offset);

	H5File *_file;
};

}



#endif /* SRC_WRAPPERS_HDF5_W_HDF5_H_ */
