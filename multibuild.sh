#!/bin/bash

if [ "$#" -lt 1 ]
then
    echo "Usage:"
    echo "  ./multibuild.sh change <target_build_name>"
    echo "  ./multibuild.sh add <new_build_name>"
    echo "  ./multibuild.sh remove <build_name>"
    echo "  ./multibuild.sh rename <old_build_name> <new_build_name>"
    exit
fi



if ! find . -maxdepth 1 -type f -name ".multibuild_curr-*" | grep . >/dev/null
then
    touch .multibuild_curr-default
fi



current_build=""
for f in $(find . -maxdepth 1 -type f -name ".multibuild_curr-*")
do
    # f = "./.multibuild_curr-"
    current_build="${f:19}"
    break
done



action="${1}"
if [ "${action}" == "change" ]
then
    if [ "$#" -lt 2 ]
    then
        echo "not enough arguments"
        exit 2
    fi
    changeto="${2}"
    rm ".multibuild_curr-${current_build}"
    next_builddir="build-${changeto}"
    next_depdir="dependencies-${changeto}"
    next_lockfile=".lock-waf_linux_build-${changeto}"
    curr_builddir="build-${current_build}"
    curr_depdir="dependencies-${current_build}"
    curr_lockfile=".lock-waf_linux_build-${current_build}"
    mv build "${curr_builddir}"
    mv dependencies "${curr_depdir}"
    mv .lock-waf_linux_build "${curr_lockfile}"
    mv "${next_builddir}" build
    mv "${next_depdir}" dependencies
    mv "${next_lockfile}" .lock-waf_linux_build
    touch ".multibuild_curr-${changeto}"
elif [ "${action}" == "add" ]
then
    if [ "$#" -lt 2 ]
    then
        echo "not enough arguments"
        exit 2
    fi
    newname="${2}"
    mkdir "build-${newname}"
    mkdir "dependencies-${newname}"
elif [ "${action}" == "remove" ]
then
    if [ "$#" -lt 2 ]
    then
        echo "not enough arguments"
        exit 2
    fi
    oldname="${2}"
    if [ "${oldname}" == "${current_build}" ]
    then
        echo "cannot remove current build"
        exit 3
    else
        rm -rf "build-${oldname}"
        rm -rf "dependencies-${oldname}"
        rm -f ".lock-waf_linux_build-${oldname}"
    fi
elif [ "${action}" == "rename" ]
then
    if [ "$#" -lt 3 ]
    then
        echo "not enough arguments"
        exit 2
    fi
    oldname="${2}"
    newname="${3}"
    if [ "${oldname}" == "${current_build}" ]
    then
        mv ".multibuild_curr-${oldname}" ".multibuild_curr-${newname}"
    else
        mv "build-${oldname}" "build-${newname}"
        mv "dependencies-${oldname}" "dependencies-${newname}"
        mv ".lock-waf_linux_build-${oldname}" ".lock-waf_linux_build-${newname}"
    fi
else
    echo "wrong action"
    exit 1
fi



exit 0
