#! \file   Make.inc
#  \brief  fundamental setting of Makefile
#  \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
#  \date   Jul. 12th 2015

# This file is part of dissection solver
# 
# Dissection solver is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dissection solver is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with diisection solver.  If not, see <http://www.gnu.org/licenses/>.

OBJS       = $(SOURCESCPP:%.cpp=%.o)

all: $(OBJS) 

.SUFFIXES: .cpp

.cpp.o:
	$(CXX) $(CCFLAGS) -c $< -o $@

%.d: %.cpp
	$(CXX) $(MAKE_DEP_OPT) $(CCFLAGS) $< > $@.$$$$; \
        sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

include $(SOURCESCPP:%.cpp=%.d)

clean:
	echo "Cleaning the trash..."
	@rm -fr *~ *.o *.so core *.d *.a
