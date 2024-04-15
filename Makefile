#Specify the default compiler
CXX = g++

#Specify the -fpic flag
CXXFLAGS += -fpic

# Additional LDFLAGS for library
ADDED_LDFLAGS = -L${HOME}/FSL-topup/warpfns -L${HOME}/FSL-topup/meshclass -L${HOME}/FSL-topup/znzlib -L${HOME}/FSL-topup/miscmaths -L${HOME}/FSL-topup/basisfield -L${HOME}/FSL-topup/utils -lfsl-warpfns -lfsl-meshclass -lfsl-znz -lfsl-miscmaths -lfsl-basisfield -lfsl-utils

#Define source files
SRCS = topup_matrices topup_file_io topup_costfunctions topupfns displacement_vector topup applytopup 
#Define object files
OBJS = $(SRCS:.cc=.o) 

#Define library source files and directories
LIB_DIRS = warpfns basisfield miscmaths newimage NewNifti cprob utils znzlib meshclass
LIB_SRCS = $(foreach dir,$(LIB_DIRS),$(wildcard $(dir)/*.cc))
LIB_OBJS = $(LIB_SRCS:.cc=.o) 

XFILES   = topup applytopup
FXFILES  = test_displacement_vector
SOFILES  = libfsl-topup.so

USRLDFLAGS += -rdynamic

all:topup applytopup libfsl-topup.so

libfsl-topup.so: topup_matrices.o topup_file_io.o topup_costfunctions.o topupfns.o displacement_vector.o
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS} 

topup: libraries topup.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ topup.o libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} ${ADDED_LDFLAGS} -lblas -llapack -lz

applytopup: libraries applytopup.o  libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ applytopup.o libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} ${ADDED_LDFLAGS} -lblas -llapack -lz

test_displacement_vector: libraries test_displacement_vector.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS}  -lblas -llapack -lz

#Rule to build object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#Phony target to build all libraries
.PHONY: libraries
libraries:
	@for dir in $(LIB_DIRS); do \
	$(MAKE) -C $$dir CXX=$(CXX) CXXFLAGS='$(CXXFLAGS)' LDFLAGS='$(LDFLAGS)'; \
 	done
#Clean rule
clean:
	rm -f topup applytopup $(shell find . -type f \( -name "*.o" -o -name "*.so" \))

