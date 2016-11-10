#! \file   Make.inc
#  \brief  source files in Driver/ directory and depending libraries
#  \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
#  \date   Jul. 12th 2015
#  \date   Feb. 29th 2016

# This file is part of Dissection
# 
# Dissection is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dissection is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

OBJS  = $(SOURCESCPP:%.cpp=%.o) $(SOURCESF77:%.for=%.o) $(SOURCESF90:%.f90=%.o)

all: $(OBJS) 

.SUFFIXES: .for .f90 .cpp

.cpp.o:
	$(CXX) $(CCFLAGS) -c $< -o $@

%.d: %.cpp
	$(CXX) $(MAKE_DEP_OPT) $(CCFLAGS) $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

include $(SOURCESCPP:%.cpp=%.d)


clean:
	echo "Cleaning the trash..."
	@rm -fr *~ *.o *.so core *.d *.mod *.a *.i90
