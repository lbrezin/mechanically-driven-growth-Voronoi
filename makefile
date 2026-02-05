#standard places to find cuda files
CUDA_INC = $(SCC_CUDA_INCLUDE)
CUDA_LIB = $(SCC_CUDA_LIB)
#CUDA_LIB2 = /usr/local/cuda/lib

CXX := g++
CC := gcc
LINK := g++ #-fPIC
NVCC := nvcc

INCLUDES = -I. -I./src/ -I./ext_src/ -I./inc/ -I$(CUDA_INC) -I$(SCC_CGAL_INCLUDE) -I/usr/include
INCLUDES += -I$(SCC_BOOST_INCLUDE) -I$(SCC_NETCDF_INCLUDE) -I$(SCC_NETCDF_CXX_INCLUDE) -I$(SCC_EIGEN_INCLUDE)
LIB_CUDA = -L. -L$(CUDA_LIB) -lcuda -lcudart
LIB_CGAL = -L$(SCC_CGAL_LIB) -L/usr/lib64 -lCGAL -lCGAL_Core -lgmp -lmpfr
LIB_NETCDF = -lnetcdf -lnetcdf_c++ -L$(SCC_NETCDF_LIB) -L$(SCC_NETCDF_CXX_LIB)
#LIB_NETCDF = -lnetcdf_c++ -L$(SCC_NETCDF_LIB) -L$(SCC_NETCDF_CXX_LIB)

#common flags
COMMONFLAGS += -std=c++11 -DCGAL_DISABLE_ROUNDING_MATH_CHECK -O0
NVCCFLAGS += -g -G -arch=sm_35 -D_FORCE_INLINES $(COMMONFLAGS) -Wno-deprecated-gpu-targets #-Xptxas -fmad=false#-O0#-dlcm=ca#-G
CXXFLAGS += $(COMMONFLAGS) $(INCLUDES) 
CXXFLAGS += -g -w -frounding-math -Wabi-tag
CFLAGS += -g $(COMMONFLAGS) -frounding-math

CUOBJ_DIR=obj/cuobj
MODULES = databases models updaters utility
INCLUDES += -I./inc/databases -I./inc/models -I./inc/updaters -I./inc/utility -I./inc/analysis
OBJ_DIR=obj
SRC_DIR=src
BIN_DIR=.

.SECONDARY:

PROGS := $(wildcard *.cpp)
PROG_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.main.o,$(PROGS))
PROG_MAINS := $(patsubst %.cpp,$(BIN_DIR)/%.out,$(PROGS))


CPP_FILES := $(wildcard src/*/*.cpp)
CPP_FILES += $(wildcard src/*.cpp)
CU_FILES := $(wildcard src/*/*.cu)
CU_FILES += $(wildcard src/*.cu)

CLASS_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_FILES))
CU_OBJS := $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.cu.o,$(CU_FILES))


#cuda objects
$(OBJ_DIR)/%.cu.o : $(SRC_DIR)/%.cu 
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(LIB_CUDA)  -o $@ -c $<

#cpp class objects
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(LIB_CUDA) $(LIB_NETCDF) $(LIB_CGAL) -o $@ -c $<

#program objects
$(OBJ_DIR)/%.main.o: %.cpp
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(LIB_CUDA) $(LIB_NETCDF) $(LIB_CGAL) -o $@ -c $<

#Programs
%.out: $(OBJ_DIR)/%.main.o $(CLASS_OBJS) $(CU_OBJS)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(LIB_CUDA) $(LIB_CGAL) $(LIB_NETCDF) -o $@ $+

#target rules

all:build

float: CXXFLAGS += -DSCALARFLOAT
float: NVCCFLAGS += -DSCALARFLOAT
float: build

debug: CXXFLAGS += -g -DCUDATHREADSYNC -DDEBUGFLAGUP
debug: NVCCFLAGS += -g -G -lineinfo -Xptxas --generate-line-info -DDEBUGFLAGUP # -G note that in debug mode noise will always be reproducible
debug: build
build: $(CLASS_OBJS) $(CU_OBJS) $(PROG_MAINS)  $(PROGS)

clean: 
	rm -f $(PROG_OBJS) $(CLASS_OBJS) $(CU_OBJS) $(PROG_OBJS) $(PROG_MAINS)

print-%  : ; @echo $* = $($*)

