
./waf configure --skip-suitesparse --skip-mkl --skip-blas --skip-lapack
./waf $1
./waf configure --skip-suitesparse --skip-mkl --skip-blas --skip-lapack --intwidth=64
./waf $1
