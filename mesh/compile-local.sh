rm -fr bin
mkdir bin

cp Makefile.Local bin/Makefile

cp -R src/* bin

make -C bin cleanall
make -C bin lib

cp bin/libesmesh.so ../libs

