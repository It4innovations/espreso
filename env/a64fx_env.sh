#!/bin/bash

export MODULESHOME=/apps/aarch64/all/Lmod/8.7.8/lmod/lmod
export LMOD_PKG=/apps/aarch64/all/Lmod/8.7.8/lmod/lmod
export LMOD_CMD=/apps/aarch64/all/Lmod/8.7.8/lmod/lmod/libexec/lmod
export LMOD_DIR=/apps/aarch64/all/Lmod/8.7.8/lmod/lmod/libexec

if [ 1 ]; then

  # arch and flags
  ARCH=$(lscpu | grep Architecture | awk '{print $2}')
  AVX512=$(lscpu | grep avx512)
  MFLAGS=
  MPATH=

  if [ "$ARCH" == "aarch64" ]; then
    MFLAGS="aarch64"
  else
    if [ -z "$AVX512" ]; then
       MFLAGS="avx2"
    else
       MFLAGS="avx512"
    fi
  fi

  echo
  echo "  $MFLAGS modules + all modules"
  echo

  MPATH="/apps/$MFLAGS/modules/*"

  MODULEPATH=""

  for dir in $MPATH
  do
    # Exclude following directories
    if [[ (${dir##*/} == "all") ]]; then
        continue
    fi
    # In case that it's directory
    if [ -d "$dir" ]; then
        if [ -z "$MODULEPATH" ]; then
            MODULEPATH="$dir"
        else
            MODULEPATH="$MODULEPATH:$dir"
        fi
    fi
  done

  for dir in /apps/modules/*
  do
    # Exclude following directories
    if [[ (${dir##*/} == "all") ]]; then
        continue
    fi
    # In case that it's directory
    if [ -d "$dir" ]; then
        if [ -z "$MODULEPATH" ]; then
            MODULEPATH="$dir"
        else
            MODULEPATH="$MODULEPATH:$dir"
        fi
    fi
  done

  # export
  export MODULEPATH
  export CLUSTERNAME="CS"
  export ARCH
  export MFLAGS
fi
