#!/bin/bash

url="${1}"

if command -v wget >/dev/null 2>/dev/null
then
    wget "${url}"
elif command -v curl >/dev/null 2>/dev/null
then
    curl -sOL "${url}"
else
    echo "ERROR: no way to download a file"
    exit 1
fi
