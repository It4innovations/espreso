
#include "progresslogger.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"

#include <cstdlib>
#include <cstring>
#include <string>

using namespace espreso;

ProgressFileLogger::ProgressFileLogger()
: out(NULL), err(NULL), bsize(0)
{

}

ProgressFileLogger::~ProgressFileLogger()
{
	if (out) {
		fclose(out);
	}
	if (err) {
		fclose(err);
	}
}

void ProgressFileLogger::setLogFile()
{
	if (rank == 0) {
		std::string outfile, errfile;
		if (info::mpi::isize == 1) {
			outfile = info::ecf->outpath + "/" + info::ecf->name + ".log";
			errfile = info::ecf->outpath + "/" + info::ecf->name + ".err";
		} else {
			outfile = info::ecf->outpath + "/" + info::ecf->name + "." + std::to_string(info::mpi::irank) + ".log";
			errfile = info::ecf->outpath + "/" + info::ecf->name + "." + std::to_string(info::mpi::irank) + ".err";
		}
		out = fopen(outfile.c_str(), "a");
		if (out == NULL) {
			eslog::error("Cannot create log file '%s'.\n", outfile.c_str());
		}
		setvbuf(out, NULL, _IOLBF, 256);
		fprintf(out, "%s", prebuffer);

		err = fopen(errfile.c_str(), "a");
		setbuf(err, NULL);
	}
}

void ProgressFileLogger::closeLogFile()
{
	if (out) {
		fclose(out);
		out = NULL;
	}
	if (err) {
		fclose(err);
		err = NULL;
	}
}
