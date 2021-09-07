
#ifndef SRC_BASIS_UTILITIES_SYSUTILS_H_
#define SRC_BASIS_UTILITIES_SYSUTILS_H_

#include <string>
#include <vector>

namespace espreso {

namespace step { struct Step; }

namespace utils {

std::string getFileDirectory(const std::string &file);
std::string createDirectory(const std::vector<std::string> &path);
void createDirectory(const std::string &path);
void createSymlink(const std::string &path, const std::string &link);
void copyFile(const std::string &source, const std::string &destination);
bool exists(const std::string &path);
void remove(const std::string &path);
std::string debugDirectory();
std::string debugDirectory(step::Step &step);
std::string prepareFile(const std::string &directory, const std::string &name, int domain = -1);
std::string filename(const std::string &directory, const std::string &name);
void listDirectory(const std::string &dir, std::vector<std::string> &files);

void printStack();
std::string getStack();

void callusleep(int usec);
int nprocs();

}
}


#endif /* SRC_BASIS_UTILITIES_SYSUTILS_H_ */
