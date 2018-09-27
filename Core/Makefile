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
KPP  :=$(ROOT)/KPP
SANDS:=$(ROOT)/SANDS

# Set src folder
SRC_DIR := src

CURR_DIR := Core

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

USER_DEFS += -std=c++11

# %%% Set default compiler %%%
#
# %%% If COMP is not defined, default to the $(CXX) variable, which %%%
# %%% is set in your .bashrc %%%
ifndef COMP
	COMP := $(CXX)
endif

ifeq ($(COMP),g++)

	# Base set of compiler flags
	CFLAGS := -Wall -Wextra -pedantic-errors #-DUSE_NETCDF4 #-Werror

	# Default optimization level for all routines (-O3)
	ifndef OPT
		OPT := -O3
		CFLAGS += $(OPT)
	endif

	CFLAGS += $(USER_DEFS)
	
endif

# -g    adds debugging information to the executable file
# -Wall turns on most compiler warnings
# -O3   turns on optimization mode #3

# Define library folder
LIB_DIR := $(ROOT)/lib

# Add path to APCEMM's library folder
LDFLAGS := -L$(LIB_DIR)

# Add path/to/lib folder
LDFLAGS := $(LDFLAGS)

# Define any libraries to link into executable: use -llibname option
LDFLAGS := $(LDFLAGS) -lstdc++ -lm

# Add custom libraries
LDFLAGS := $(LDFLAGS) -lSands -lfftw3 -lKpp -lHeaders -lEpm -lnetcdf_c++

ifndef BUILD
	BUILD := $(ROOT)/build
endif
ifndef APP_DIR
	APP_DIR := $(BUILD)/apps
endif

OBJ_DIR := $(BUILD)/objects/$(CURR_DIR)

# Define the executable file: TARGET
TARGET := APCEMM

# Define any directories containing header files
INCLUDE := -Iinclude

HEADERS = $(wildcard include/*.h*) 

# Define the C++ source files
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# Define the C++ object files
OBJ = ${SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o}

all:	build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(COMP) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJ) $(KPPOBJ)
	@mkdir -p $(@D)
	$(COMP) $(CFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJ) $(LDFLAGS)
	@echo "---> Executable file $(TARGET) has been compiled! <---"

libkpp:
	@$(MAKE) -C $(KPP)

libsands:
	@$(MAKE) -C $(SANDS)

libhelp:
	@echo " --> In libKpp.a: "
	@$(AR) -t $(LIB_DIR)/libKpp.a
	@echo " --> In libSands.a: "
	@$(AR) -t $(LIB_DIR)/libSands.a

.PHONY: all libkpp libsands libhelp build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CFLAGS += -DDEBUG -g
debug: all

release: CFLAGS += $(OPT)
release: all

clean:
	@echo "---> Making clean in $(CURR_DIR) <---"
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*




