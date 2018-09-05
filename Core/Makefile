#------------------------------------------------------------------------------
#      							C++ - APCEMM                                  #
#      A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)       #
#------------------------------------------------------------------------------
#
# !MODULE: Makefile (Main-level)
#
# !REMARKS:
#                                                                             .
# To display a complete list of options, type "make help".
# To clean the current folder of .o and excutable files, type "make clean".

# 'make'        build executable file 'APCEMM'
# 'make clean'  removes all .o and executable files

# Set default shell to bash
SHELL :=/bin/bash 

# Set root directory
ROOT := ..

# Set src folder
SRC_FOLDER := src

CURR_FOLDER := Main

#---------------------------------------------------
# Compiler settings
#---------------------------------------------------

# %%% OpenMP parallelization (off by default) %%%
ifndef OMP
	OMP :=no
endif

# Option to turn off OpenMP for testing
REGEXP :=(^[Nn]|^[Oo])
ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true) 
	USER_DEFS += -DNO_OMP
endif

# %%% Set default compiler %%%
#
# %%% If COMP is not defined, default to the $(CXX) variable, which %%%
# %%% is set in your .bashrc %%%
ifndef COMP
	COMP := $(CXX)
endif

ifeq ($(COMP),g++)

	# Base set of compiler flags
	CFLAGS := -Wall -Wextra -pedantic-errors #-Werror

	# Default optimization level for all routines (-O2)
	ifndef OPT
		OPT := -O2
		CFLAGS += $(OPT)
	endif

	CFLAGS += $(USER_DEFS)
	
endif

# -g    adds debugging information to the executable file
# -Wall turns on most compiler warnings
# -O3   turns on optimization mode #3

# Define any libraries to link into executable: use -llibname option
LDFLAGS := -lstdc++ -lm -lfftw3

ifndef BUILD
	BUILD := $(ROOT)/build
endif
ifndef APP_DIR
	APP_DIR := $(BUILD)/apps
endif
OBJ_DIR := $(BUILD)/objects/$(CURR_FOLDER)

# Define the executable file: TARGET
TARGET := APCEMM

# Define any directories containing header files
INCLUDE := -I./include

HEADERS = $(wildcard ./include/*.h) 

# Define the C++ source files
SRC = $(wildcard $(SRC_FOLDER)/*.cpp)

# Define the C++ object files
OBJ = ${SRC:$(SRC_FOLDER)/%.cpp=$(OBJ_DIR)/%.o}

all:	build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: $(SRC_FOLDER)/%.cpp
	@mkdir -p $(@D)
	$(COMP) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJ)
	@mkdir -p $(@D)
	$(COMP) $(CFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJ) $(LDFLAGS)
	@echo "---> Executable file $(TARGET) has been compiled! <---"

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CFLAGS += -DDEBUG -g
debug: all

release: CFLAGS += $(OPT)
release: all

clean:
	@echo "---> Making clean <---"
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*




