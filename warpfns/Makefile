# A Makefile for warpfns library

#include ${FSLCONFDIR}/default.mk

PROJNAME = warpfns
LIBS     = -lfsl-basisfield -lfsl-meshclass -lfsl-newimage -lfsl-miscmaths -lfsl-NewNifti -lfsl-utils

all: libfsl-warpfns.so

libfsl-warpfns.so: warpfns.o fnirt_file_reader.o fnirt_file_writer.o point_list.o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ ${LDFLAGS}

test_parallel_warpfns: test_parallel_warpfns.o warpfns.o fnirt_file_reader.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
