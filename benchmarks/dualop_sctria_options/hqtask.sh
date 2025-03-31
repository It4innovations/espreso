#!/bin/bash

taskdir="${1}"

touch "${taskdir}/.started"

source env/mn5....sh

source "${taskdir}/config.sh"



# run timeout espreso
# > "${taskdir}/out.txt" 2> "${taskdir}/err.txt"



touch "${taskdir}/.finished"
