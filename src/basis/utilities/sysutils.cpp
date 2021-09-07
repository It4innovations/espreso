
#include "sysutils.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/stepinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"

#include <execinfo.h>
#include <dirent.h>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/sysinfo.h>

namespace espreso {
namespace utils {

std::string getFileDirectory(const std::string &file)
{
	if (file.find("/") == std::string::npos) {
		return ".";
	}
	return file.substr(0, file.find_last_of('/'));
}

std::string createDirectory(const std::vector<std::string> &path)
{
	std::stringstream prefix;
	std::for_each(path.begin(), path.end(), [&] (const std::string &dir) { prefix << dir << "/"; });

	createDirectory(prefix.str());
	return prefix.str();
}

void createDirectory(const std::string &path)
{
	if (std::system(("mkdir -p " + path).c_str())) {
		eslog::error("Cannot create directory '%s'\n", path.c_str());
	}
}

bool exists(const std::string &path)
{
	return std::ifstream(path.c_str()).good();
}

void remove(const std::string &path)
{
	if (std::remove(path.c_str())) {
		eslog::error("Cannot remove '%s'\n", path.c_str());
	}
}

void createSymlink(const std::string &path, const std::string &link)
{
	if (symlink(path.c_str(), link.c_str())) {
		eslog::error("Cannot create link '%s' to '%s'\n", link.c_str(), path.c_str());
	}
}

void copyFile(const std::string &source, const std::string &destination)
{
	std::ifstream src(source.c_str(), std::ios::binary);
	if (!src.good()) {
		eslog::error("Cannot read file '%s'\n", source.c_str());
	}
	std::ofstream dst(destination, std::ios::binary);
	if (!dst.good()) {
		eslog::error("Cannot create file '%s'\n", destination.c_str());
	}
	dst << src.rdbuf();
}

std::string debugDirectory()
{
	std::stringstream path;
	path << info::ecf->outpath << "/DEBUG";
	path << "/loadstep" << step::step.loadstep;
	path << "/substep" << step::step.substep;
	path << "/iteration" << step::step.iteration;
	path << "/" << info::mpi::rank;
	return path.str();
}

std::string debugDirectory(step::Step &step)
{
	std::stringstream path;
	path << info::ecf->outpath << "/DEBUG";
	path << "/loadstep" << step.loadstep;
	path << "/substep" << step.substep;
	path << "/iteration" << step.iteration;
	path << "/" << info::mpi::rank;
	return path.str();
}

std::string prepareFile(const std::string &directory, const std::string &name, int domain)
{
	createDirectory(directory);
	if (domain != -1) {
		return std::string(directory + "/" + name + std::to_string(domain) + ".txt");
	} else {
		return std::string(directory + "/" + name + ".txt");
	}
}

std::string filename(const std::string &directory, const std::string &name)
{
	createDirectory(directory);
	return directory + "/" + name;
}

void listDirectory(const std::string &dir, std::vector<std::string> &files)
{
	DIR *d = opendir(dir.c_str());
	dirent *dp;
	while ((dp = readdir(d)) != NULL) {
		if (dp->d_type == DT_REG) {
			files.push_back(dp->d_name);
		}
	}
	closedir(d);
}

void callusleep(int usec)
{
	usleep(usec);
}

int nprocs()
{
	return get_nprocs();
}

void printStack()
{
	pid_t pid = getpid();
	std::string pstack = "pstack " + std::to_string(pid);
	if (system(pstack.c_str())) {
		// stack cannot be printed
	}
}

std::string getStack()
{
	std::vector<void*> stack(30);
	size_t size = backtrace(stack.data(), 30);
	char** functions = backtrace_symbols(stack.data(), size);

	std::stringstream command;
	command << "addr2line -sipfC -e " << info::ecf->exe;
	for (size_t i = 0; i < size; i++) {
		std::string function(functions[i]);
		size_t begin = function.find_last_of('[') + 1;
		size_t end = function.find_last_of(']');
		command << " " << function.substr(begin, end - begin);
	}
	command << "\n";
	free(functions);

	FILE *in;
	char buff[512];
	printf("popen\n");
	if(!(in = popen(command.str().c_str(), "r"))){
		return "Broken address to file lines command";
	}
	printf("popened\n");

	std::string message;
	while(fgets(buff, sizeof(buff), in) != NULL){
		message += buff;
	}
	pclose(in);

	return message;
}


}
}





