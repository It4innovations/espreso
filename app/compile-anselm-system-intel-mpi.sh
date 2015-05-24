#!/bin/bash

rm -fr bin
mkdir bin

cp Makefile.Anselm bin/Makefile

cp -R src/* bin

make -C bin cleanall
make -C bin opt -j 16

cp bin/espreso ../

