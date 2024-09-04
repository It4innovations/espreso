#!/bin/bash

find dependencies/ -maxdepth 1 -mindepth 1 -type d | xargs rm -rf
