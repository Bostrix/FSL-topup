#Specify the default compiler
CXX = g++

#Specify the -fpic flag
CXXFLAGS += -fpic

#Define source files
SRCS = topup_matrices topup_file_io topup_costfunctions topupfns displacement_vector topup applytopup
#Define object files
OBJS = $(SRCS:.cc=.o) 

#Define library source files and directories
LIB_DIRS = warpfns basisfield miscmaths newimage NewNifti utils cprob znzlib meshclass
LIB_SRCS = $(foreach dir,$(LIB_DIRS),$(wildcard $(dir)/*.cc))
LIB_OBJS = $(LIB_SRCS:.cc=.o) 

XFILES   = topup applytopup
FXFILES  = test_displacement_vector
SOFILES  = libfsl-topup.so


all:topup applytopup libfsl-topup.so

libfsl-topup.so: topup_matrices.o topup_file_io.o topup_costfunctions.o topupfns.o displacement_vector.o
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

topup:libraries topup.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

applytopup: libraries applytopup.o  libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

test_displacement_vector: libraries test_displacement_vector.o libfsl-topup.so $(LIB_OBJS)
	${CXX} ${CXXFLAGS} -o $@ $< libfsl-topup.so $(LIB_OBJS) ${LDFLAGS} -lblas -llapack -lz

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
	rm -f $(XFILES) $(LIB_OBJS) $(shell find . -type f \( -name "*.o" -o -name "*.so" \))
