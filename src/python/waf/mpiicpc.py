#! /usr/bin/env python
# encoding: utf-8
# WARNING! Do not edit! https://waf.io/book/index.html#_obtaining_the_waf_file

import os,sys
from waflib.Tools import ccroot,ar,gxx
from waflib.Configure import conf
@conf
def find_mpiicpc(conf):
	if sys.platform=='cygwin':
		conf.fatal('The Intel compiler does not work on Cygwin')
	cxx=conf.find_program('mpiicpc')
	conf.env.CXX = cxx
	conf.get_cc_version(cxx,icc=True)
	conf.env.CXX_NAME='icc'
def configure(conf):
	conf.find_mpiicpc()
	conf.find_ar()
	conf.gxx_common_flags()
	conf.gxx_modifier_platform()
	conf.cxx_load_tools()
	conf.cxx_add_flags()
	conf.link_add_flags()
