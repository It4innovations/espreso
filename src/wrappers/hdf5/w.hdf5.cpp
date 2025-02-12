
#include "w.hdf5.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.hpp"

#include <string>
#include <cstring>
#include <cstddef>

#ifndef HAVE_HDF5

using namespace espreso;

HDF5::H5TypeWrapper HDF5::INT;
HDF5::H5TypeWrapper HDF5::LONG;
HDF5::H5TypeWrapper HDF5::ESINT;
HDF5::H5TypeWrapper HDF5::FLOAT;
HDF5::H5TypeWrapper HDF5::DOUBLE;

#else

#include "hdf5.h"
using namespace espreso;

struct HDF5::H5Type {
    static hid_t t_int;
    static hid_t t_long;
    static hid_t t_float;
    static hid_t t_double;
};

hid_t HDF5::H5Type::t_int = H5T_NATIVE_INT;
hid_t HDF5::H5Type::t_long = H5T_NATIVE_LONG;
hid_t HDF5::H5Type::t_float = H5T_NATIVE_FLOAT;
hid_t HDF5::H5Type::t_double = H5T_NATIVE_DOUBLE;

HDF5::H5TypeWrapper HDF5::INT = HDF5::H5TypeWrapper{ &HDF5::H5Type::t_int };
HDF5::H5TypeWrapper HDF5::LONG = HDF5::H5TypeWrapper{ &HDF5::H5Type::t_long };
HDF5::H5TypeWrapper HDF5::ESINT = HDF5::H5TypeWrapper{ sizeof(esint) == 4 ? &HDF5::H5Type::t_int :  &HDF5::H5Type::t_long };
HDF5::H5TypeWrapper HDF5::FLOAT = HDF5::H5TypeWrapper{ &HDF5::H5Type::t_float };
HDF5::H5TypeWrapper HDF5::DOUBLE = HDF5::H5TypeWrapper{ &HDF5::H5Type::t_double };

static void check(herr_t ret)
{
    if (ret < 0) { eslog::error("HDF5 error: %d\n", ret); }
}

static hid_t check(hid_t id)
{
    if (id < 0) { eslog::error("HDF5 error: %d\n", id); }
    return id;
}

class HDF5::H5File {
    hid_t handler;

public:
    H5File(const char* file, MPIGroup &mpigroup, HDF5::MODE mode)
    {
        hid_t fileacces = H5Pcreate(H5P_FILE_ACCESS);

        check(H5Pset_fapl_mpio(fileacces, mpigroup.communicator, MPI_INFO_NULL));
        switch (mode) {
        case MODE::WRITE:
            handler = check(H5Fcreate((std::string(file) + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fileacces));
            break;
        case MODE::READ:
            handler = check(H5Fopen((std::string(file)).c_str(), H5F_ACC_RDONLY, fileacces));
            break;
        }
        check(H5Pclose(fileacces));
    }

    ~H5File()
    {
        H5Fclose(handler);
    }

    void append(const char* name, const H5TypeWrapper &sourcetype, const H5TypeWrapper &targettype, const void *data, esint esize, esint nelements, esint offset, esint totalsize)
    {
        hsize_t h5start[2] = { (hsize_t)offset, 0 };
        hsize_t h5size[2] = { (hsize_t)nelements, (hsize_t)esize };
        hsize_t h5totalsize[2] = { (hsize_t)totalsize, (hsize_t)esize };

        // describe my data organisation
        hid_t h5memory = check(H5Screate_simple(esize == 1 ? 1 : 2, h5size, h5size));
        hid_t h5dataspace = check(H5Screate_simple(esize == 1 ? 1 : 2, h5totalsize, h5totalsize));

        // describe file data organisation
        hid_t h5dataset = check(H5Dcreate2(handler, name, *static_cast<hid_t*>(targettype.data), h5dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        hid_t h5filespace = check(H5Dget_space(h5dataset));
        check(H5Sselect_hyperslab(h5filespace, H5S_SELECT_SET, h5start, NULL, h5size, NULL));

        // describe transfer
        hid_t h5xfer = check(H5Pcreate(H5P_DATASET_XFER));
        check(H5Pset_dxpl_mpio(h5xfer, H5FD_MPIO_INDEPENDENT));

        // append data to file
        check(H5Dwrite(h5dataset, *static_cast<hid_t*>(sourcetype.data), h5memory, h5filespace, h5xfer, data));

        check(H5Sclose(h5memory));
        check(H5Pclose(h5xfer));
        check(H5Sclose(h5filespace));
        check(H5Dclose(h5dataset));
        check(H5Sclose(h5dataspace));
    }

    void read(const char* name, const H5TypeWrapper &targettype, void *data, esint esize, esint nelements, esint offset)
    {
        hsize_t h5start[2] = { (hsize_t)offset, 0 };
        hsize_t h5size[2] = { (hsize_t)nelements, (hsize_t)esize };

        // describe my data organization
        hid_t h5memory = check(H5Screate_simple(esize == 1 ? 1 : 2, h5size, h5size));

        // describe file data organization
        hid_t h5dataset = check(H5Dopen2(handler, name, H5P_DEFAULT));
        hid_t h5filespace = check(H5Dget_space(h5dataset));
        check(H5Sselect_hyperslab(h5filespace, H5S_SELECT_SET, h5start, NULL, h5size, NULL));

        // describe transfer
        hid_t h5xfer = check(H5Pcreate(H5P_DATASET_XFER));
        check(H5Pset_dxpl_mpio(h5xfer, H5FD_MPIO_INDEPENDENT));

        // read data from file
        check(H5Dread(h5dataset, *static_cast<hid_t*>(targettype.data), h5memory, h5filespace, h5xfer, data));

        check(H5Sclose(h5memory));
        check(H5Pclose(h5xfer));
        check(H5Sclose(h5filespace));
        check(H5Dclose(h5dataset));
    }
};

#endif

bool HDF5::islinked()
{
#ifndef HAVE_HDF5
    return false;
#else
    return true;
#endif
}

HDF5::HDF5(const char* file, MPIGroup &mpigroup, MODE mode)
: _file(NULL)
{
#ifndef HAVE_HDF5
    eslog::globalerror("ESPRESO run-time error: cannot call HDF5 library (the library is not linked).\n");
#else
    _file = new H5File(file, mpigroup, mode);
#endif
}

HDF5::~HDF5()
{
#ifdef HAVE_HDF5
    if (_file) { delete _file; }
#endif
}

void HDF5::append(
        const char* name, const H5TypeWrapper &source, const H5TypeWrapper &target,
        const void *data, esint esize, esint nelements, esint offset, esint totalsize)
{
#ifdef HAVE_HDF5
    _file->append(name, source, target, data, esize, nelements, offset, totalsize);
#endif
}

void HDF5::read(
        const char* name, const H5TypeWrapper &target,
        void *data, esint esize, esint nelements, esint offset)
{
#ifdef HAVE_HDF5
    _file->read(name, target, data, esize, nelements, offset);
#endif
}



