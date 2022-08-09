
#ifndef SRC_BASIS_IO_LOADER_H_
#define SRC_BASIS_IO_LOADER_H_

#include "wrappers/mpi/communication.h"
#include <cstdio>

#include "esinfo/mpiinfo.h"

namespace espreso {

struct Loader {
	Loader() {}
	virtual ~Loader() {}

	virtual int    open(const MPIGroup &group, const std::string &file) =0;
	virtual size_t size() =0;
	virtual void   read(char *data, size_t offset, size_t size) =0;
	virtual void   iread(char *data, size_t offset, size_t size) =0;

	virtual void   wait() =0;
};

struct POSIXLoader: public Loader {

	POSIXLoader(): f(nullptr) {}
	~POSIXLoader() { if (f) fclose(f); }

	int open(const MPIGroup &group, const std::string &file)
	{
		return (f = fopen(file.c_str(), "rb")) == NULL;
	}

	size_t size()
	{
		fseek(f, 0L, SEEK_END);
		return ftell(f);
	}

	void read(char *data, size_t offset, size_t size)
	{
		fseek(f, offset, SEEK_SET);
		size = fread(data, 1, size, f);
	}

	void iread(char *data, size_t offset, size_t size)
	{
		fseek(f, offset, SEEK_SET);
		size = fread(data, 1, size, f);
	}

	void wait()
	{

	}

protected:
	FILE *f;
};

struct MPILoader: public Loader {
	MPILoader(): MPIfile(nullptr), waiting{false} {}
	~MPILoader() { if (MPIfile) MPI_File_close(&MPIfile); }

	int open(const MPIGroup &group, const std::string &file)
	{
		return MPI_File_open(MPI_COMM_SELF, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile);
	}

	size_t size()
	{
		MPI_Offset size;
		MPI_File_get_size(MPIfile, &size);
		return size;
	}

	void read(char *data, size_t offset, size_t size)
	{
		MPI_File_read_at(MPIfile, offset, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
	}

	void iread(char *data, size_t offset, size_t size)
	{
		MPI_File_iread_at(MPIfile, offset, data, size, MPI_BYTE, &request);
		waiting = true;
	}

	void wait()
	{
		if (waiting) {
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			waiting = false;
		}
	}

protected:
	MPI_File MPIfile;
	MPI_Request request;
	bool waiting;
};

struct MPICollectiveLoader: public Loader {
	MPICollectiveLoader(): MPIfile(nullptr), waiting{false} {}
	~MPICollectiveLoader() { if (MPIfile) MPI_File_close(&MPIfile); }

	int open(const MPIGroup &group, const std::string &file)
	{
		return MPI_File_open(group.communicator, file.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &MPIfile);
	}

	size_t size()
	{
		MPI_Offset size;
		MPI_File_get_size(MPIfile, &size);
		return size;
	}

	void read(char *data, size_t offset, size_t size)
	{
		MPI_File_read_at_all(MPIfile, offset, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
	}

	void iread(char *data, size_t offset, size_t size)
	{
		MPI_File_iread_at_all(MPIfile, offset, data, size, MPI_BYTE, &request);
		waiting = true;
	}

	void wait()
	{
		if (waiting) {
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			waiting = false;
		}
	}

protected:
	MPI_File MPIfile;
	MPI_Request request;
	bool waiting;
};

}



#endif /* SRC_BASIS_IO_LOADER_H_ */
