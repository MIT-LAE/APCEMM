# Set here the desired C compiler and its optimization options
CC   = g++
COPT = -O3 -Wall  

ifndef ROOT
	# Set root directory
	ROOT := ..
endif

# Set src folder
SRC_DIR := ./src

CURR_DIR := SANDS

# Define any libraries to link into executable
LDFLAGS := -lm -lstdc++ -lfftw3

ifndef BUILD
	BUILD := $(ROOT)/build
endif
ifndef APP_DIR
	APP_DIR := $(BUILD)/apps
endif
OBJ_DIR := $(BUILD)/objects/$(CURR_DIR)

# Define output library:
LIB := libSands.a

# Define lib folder:
LIB_DIR := $(ROOT)/lib

# Define the executable file:
TARGET := SANDS

# Define any directories containing header files
INCLUDE := -Iinclude

HEADERS = $(wildcard ./include/*.h*)

# Define the C++ source files
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# Define the C++ object files
OBJ = ${SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o}


.PHONY: all build clean

all:    build $(APP_DIR)/$(TARGET) $(LIB_DIR)/$(LIB)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(@D)
	$(CC) $(COPT) $(INCLUDE) $(LDFLAGS) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(HEADERS) $(OBJ)
	@mkdir -p $(@D)

$(LIB_DIR)/$(LIB): $(OBJ)
	@mkdir -p $(LIB_DIR)
	@$(AR) cr $@ $^ 
	@echo "Library $(LIB) has been created"

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(LIB_DIR)

clean:
	@echo "---> Making clean in $(CURR_DIR) <---"
	-@rm -rv $(OBJ_DIR)/*
	-@rm -rv $(LIB_DIR)/$(LIB)

