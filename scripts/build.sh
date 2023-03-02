

declare -a COMPILERS=("INTEL" "GCC")


for COMPILER in "${COMPILERS[@]}"
do
    if [ "$COMPILER" = "INTEL" ]
    then
        ml purge
        if [ "$1" = "barbora" ]
        then
            source env/modules.barbora.icpc
        else
            source env/modules.karolina.icpc
        fi
        FLAGS=--cxxflags="-march=native -mtune=native"
    elif [ "$COMPILER" = "GCC" ]
    then
        ml purge
        if [ "$1" = "barbora" ]
        then
            source env/modules.barbora.gcc
        else
            source env/modules.karolina.gcc
        fi

        FLAGS=--cxxflags="-march=native -mtune=native"
    else
        echo "Incorrect compiler. Available compilers are INTEL, GCC"
        exit 1
    fi

    ./waf configure "$FLAGS" --out=./prebuild/${1}/${COMPILER}
    ./waf
    
done
