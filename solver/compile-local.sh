rm -fr bin
mkdir bin

cp Makefile.Local bin/Makefile

cp src/* bin

make -C bin cleanall
make -C bin lib -j 16

