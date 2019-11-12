###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
CUDA_ROOT_DIR=/usr/tce/packages/cuda/cuda-9.2.148/

##########################################################

## CC COMPILER OPTIONS ##

# CC compiler options:
CC=g++
CC_FLAGS= -Wall -pedantic -std=c++11 -O0
# LA_FlAGS=-llapack -lblas
CC_LIBS=-framework Accelerate

##########################################################

## NVCC COMPILER OPTIONS ##

# NVCC compiler options:
NVCC=nvcc
NVCC_FLAGS=
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = ./src

# Object file directory:
OBJ_DIR = ./bin

# Include header file diretory:
INC_DIR = ./include

##########################################################

## Make variables ##

# Target executable name:
EXE = gp

# Object files:
OBJS = $(OBJ_DIR)/main.o \
		$(OBJ_DIR)/Vector.o \
		$(OBJ_DIR)/Matrix.o \
		$(OBJ_DIR)/kernel.o
		#  $(OBJ_DIR)/SetupJacobian.o

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -o $@

# Compile main .cpp file to object files:
$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(INC_DIR)/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

# Clean objects in object directory.
clean:
	$(RM) $(OBJ_DIR)/*.o $(EXE)
