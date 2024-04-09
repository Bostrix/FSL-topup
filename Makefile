#Specify the default compiler
CXX = g++

#Specify the -fpic flag
CXXFLAGS += -fpic

#Define source files
SRCS = topup_matrices.cpp topup_file_io.cpp topup_costfunctions.cpp topupfns.cpp displacement_vector.cpp 
XSRCS = topup.cpp applytopup.cpp
#Define object files
OBJS = $(SRCS:.cc=.o) $(SRCS:.cpp=.o)
XOBJS = $(XSRCS:.cc=.o) $(XSRCS:.cpp=.o)

#Define library source files and directories
LIB_DIRS = warpfns basisfield miscmaths newimage NewNifti utils cprob znzlib meshclass
LIB_SRCS = $(foreach dir,$(LIB_DIRS),$(wildcard $(dir)/.cc))
LIB_OBJS = $(LIB_SRCS:.cc=.o) $(LIB_SRCS:.cpp=.o)

XFILES   = topup applytopup
FXFILES  = test_displacement_vector
SOFILES  = libfsl-topup.so


all:topup applytopup libfsl-topup.so

libfsl-topup.so: libraries topup_matrices.o topup_file_io.o topup_costfunctions.o topupfns.o displacement_vector.o $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -shared -o $@ $^ $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

topup: libraries topup.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

applytopup: libraries applytopup.o  libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

test_displacement_vector: libraries test_displacement_vector.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

#Phony target to build all libraries
.PHONY: libraries
libraries:
	@for dir in $(LIB_DIRS); do \
                $(MAKE) -C $$dir CXX=$(CXX) CXXFLAGS='$(CXXFLAGS)' LDFLAGS='$(LDFLAGS)' ; \
        done

#Clean rule
clean:
	rm -f $(XFILES) $(LIB_OBJS) $(shell find . -type f \( -name "*.o" -o -name "*.so" \))
