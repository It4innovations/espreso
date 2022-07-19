
#ifndef SRC_BASIS_UTILITIES_SYSUTILS_H_
#define SRC_BASIS_UTILITIES_SYSUTILS_H_

#include <string>
#include <vector>

namespace espreso {
namespace utils {

std::string getFileDirectory(const std::string &file);
std::string createDirectory(const std::vector<std::string> &path);
void createDirectory(const std::string &path);
void createSymlink(const std::string &path, const std::string &link);
void copyFile(const std::string &source, const std::string &destination);
bool exists(const std::string &path);
void remove(const std::string &path);
std::string debugDirectory();
std::string prepareFile(const std::string &directory, const std::string &name, int domain = -1);
void listDirectoryFiles(const std::string &dir, std::vector<std::string> &files);
void listDirectorySubdirectories(const std::string &dir, std::vector<std::string> &files);

void printStack();
std::string getStack();

void callusleep(int usec);
int nprocs();

}
}


#endif /* SRC_BASIS_UTILITIES_SYSUTILS_H_ */
