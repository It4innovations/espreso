
#ifndef SRC_BASIS_IO_WRITER_H_
#define SRC_BASIS_IO_WRITER_H_

#include "wrappers/mpi/communication.h"

namespace espreso {

class Writer {
public:
	virtual ~Writer() {}

	virtual int  open(MPIGroup &group, const std::string &file) =0;
	virtual void store(const char* data, size_t offset, size_t size) =0;
	virtual void close() =0;
};

class MPIWriter: public Writer {
public:
	MPIWriter(): MPIfile(NULL) {}

	int open(MPIGroup &group, const std::string &file)
	{
		return MPI_File_open(MPI_COMM_SELF, file.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile);
	}

	void store(const char* data, size_t offset, size_t size)
	{
		MPI_File_write_at(MPIfile, offset, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
	}

	void close()
	{
		MPI_File_close(&MPIfile);
	}

protected:
	MPI_File MPIfile;
};

class MPICollectiveWriter: public Writer {
public:
	MPICollectiveWriter(): MPIfile(NULL) {}

	int open(MPIGroup &group, const std::string &file)
	{
		return MPI_File_open(group.communicator, file.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile);
	}

	void store(const char* data, size_t offset, size_t size)
	{
		MPI_File_write_at_all(MPIfile, offset, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
	}

	void close()
	{
		MPI_File_close(&MPIfile);
	}

protected:
	MPI_File MPIfile;
};

class MPIAsyncWriter: public Writer {
public:
	MPIAsyncWriter(): MPIfile(NULL), rindex(0) {}

	MPIAsyncWriter(const MPIAsyncWriter&) = delete;
	MPIAsyncWriter& operator=(const MPIAsyncWriter&) = delete;

	int open(MPIGroup &group, const std::string &file)
	{
		return MPI_File_open(MPI_COMM_SELF, file.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &MPIfile);
	}

	int open(MPIGroup &group, const std::string &file, size_t nrequests)
	{
		requests.resize(nrequests);
		return open(group, file);
	}

	void store(const char* data, size_t offset, size_t size)
	{
		MPI_File_iwrite_at(MPIfile, offset, data, size, MPI_BYTE, &requests[rindex++]);
	}

	void close()
	{
		MPI_Waitall(rindex, requests.data(), MPI_STATUSES_IGNORE);
		MPI_File_close(&MPIfile);
	}

protected:
	MPI_File MPIfile;
	int rindex;
	std::vector<MPI_Request> requests;
};

// POSIX does not work
//class POSIXWriter: public Writer {
//public:
//	POSIXWriter(): f(NULL) {}
//
//	bool only32b()
//	{
//		return false;
//	}
//
//	int open(MPIGroup &group, const std::string &file)
//	{
//		return (f = fopen(file.c_str(), "wb+")) == NULL;
//	}
//
//	void resize(size_t totalsize)
//	{
//		setbuf(f, NULL);
//		fseek(f, totalsize, SEEK_SET); // set the size of the file
//	}
//
//	void store(const char* data, size_t offset, size_t size)
//	{
//		fseek(f, offset, SEEK_SET);
//		fwrite(data + offset, sizeof(char), size, f);
//	}
//
//	void close()
//	{
//		fclose(f);
//	}
//
//protected:
//	FILE *f;
//};

}

#endif /* SRC_BASIS_IO_WRITER_H_ */
