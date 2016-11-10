INT_WIDTH = 32 # Dissection does not support 64b int


SRC      = $(ROOT)/tools/dissection/src
INCLUDES = $(ROOT)/src/include
LIBS     = $(ROOT)/libs

CXX     = icpc

CCFLAGS = -O3 -restrict -DNO_TO_STRING -DCLOCK_GETTIME -DDD_REAL -DBLAS_MKL -DINT_WIDTH=$(INT_WIDTH) -I$(SRC) -I$(INCLUDES)

ifeq ($(LIBTYPE), SHARED)
  LD       = icpc
  CCFLAGS += -fPIC
  LDFLAGS  = -shared -L$(LIBS) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmetis -lqd -lscotch -lscotcherr -lpthread -lrt -o
  EXTLIB   = so
else
  LD      = ar
  LDFLAGS = rcs 
  EXTLIB  = a
endif

MAKE_DEP_OPT = -MM
