#------------------------------------------------------------------------------
#                 MIT Laboratory for Aviation and the Environment
#------------------------------------------------------------------------------
#                               C++ - APCEMM                                  #
#      A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)       #
#------------------------------------------------------------------------------
#
# !MODULE: Makefile (Core directory)
#
# !DESCRIPTION: This makefile compiles the various Core source files 
#
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the root        %%%
# %%% directory!                                                           %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the C-APCEMM code
# LIB        Specifies the directory where library files (*.a) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
#
# !REVISION HISTORY: 
# 04 Oct 2018 - T. Fritz - Initial version
# 09 Oct 2018 - T. Fritz - Added Make flags
#------------------------------------------------------------------------------

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT_DIR   := ../..
LIB_DIR    := $(ROOT_DIR)/lib
BUILD_DIR  := $(ROOT_DIR)/build
APP_DIR    := $(BUILD_DIR)/apps
INCL_DIR   := $(ROOT_DIR)/include/Core
SRC_DIR    := $(ROOT_DIR)/src
UTIL_DIR   := $(SRC_DIR)/Util
# (Kinetics Pre-Processor                  -> Chemistry)
KPP_DIR    := $(SRC_DIR)/KPP
# (Spectral Advection aNd Diffusion Solver -> Atmospheric Advection/Diffusion)
SANDS_DIR  := $(SRC_DIR)/SANDS
# (AIrcraft Microphysics                   -> Microphysics)
AIM_DIR    := $(SRC_DIR)/AIM
# (Early Plume Microphysics                -> Early Plume Microphysics)
EPM_DIR    := $(SRC_DIR)/EPM

# Include header file.  This returns GCC, SHELL, as well as the default 
# Makefile compilation rules for source code files.
include $(ROOT_DIR)/Makefile_header.mk

# List of source files: everything ending in .cpp
SOURCES    := $(wildcard *.cpp)

# List of object files (replace .cpp with .o)
OBJECTS    := ${SOURCES:%.cpp=%.o}

# List of header files
HEADERS    := $(wildcard $(INCL_DIR)/*.h*)

# List of library file
LIBRARY    := libCore.a

# Executable file
TARGET     := APCEMM.sh

###############################################################################
###                                                                         ###
###                         Makefile targets                                ###
###                                                                         ###
###############################################################################

.PHONY: all build lib exe clean realclean wipeout help debug
.PHONY: libUtil libKpp libSands libAim libEpm libCore

all:
	@$(MAKE) $(MAKEFLAGS) build
	@$(MAKE) $(MAKEFLAGS) lib
	@$(MAKE) $(MAKEFLAGS) exe

build:
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(APP_DIR)

lib:
	@$(MAKE) libUtil
	@$(MAKE) libKpp
	@$(MAKE) libSands
	@$(MAKE) libAim
	@$(MAKE) libEpm
	@$(MAKE) libCore

exe:
	$(LD) $(INCLUDE) -o $(TARGET) $(OBJECTS) $(LINK)
	@echo "** ===> Executable file $(TARGET) has been created! <=== **"
	mv -f $(TARGET) $(APP_DIR)

libCore: $(OBJECTS)

libUtil:
	@$(MAKE) -C $(UTIL_DIR)

libKpp:
	@$(MAKE) -C $(KPP_DIR)
	
libSands:
	@$(MAKE) -C $(SANDS_DIR)

libAim:
	@$(MAKE) -C $(AIM_DIR)

libEpm:
	@$(MAKE) -C $(EPM_DIR)

clean:
	@echo "===> Making clean in directory: Core <==="
	-@rm -vf *.o *.a $(TARGET)

realclean:
	@$(MAKE) clean
	@$(MAKE) -C $(UTIL_DIR) clean
	@$(MAKE) -C $(KPP_DIR) clean
	@$(MAKE) -C $(SANDS_DIR) clean
	@$(MAKE) -C $(AIM_DIR) clean
	@$(MAKE) -C $(EPM_DIR) clean
	@$(MAKE) wipeout

wipeout:
	rm -f $(LIB_DIR)/*.a

help:
	@$(MAKE) -C $(HELP) help

debug:
	@echo "ROOT      : $(ROOT_DIR)"
	@echo "KPP_DIR   : $(KPP_DIR)"
	@echo "SANDS_DIR : $(SANDS_DIR)"
	@echo "AIM_DIR   : $(AIM_DIR)"
	@echo "EPM_DIR   : $(EPM_DIR)"
	@echo "LIB       : $(LIB_DIR)"
	@echo "GCC       : $(GCC)"
	@echo "SOURCES   : $(SOURCES)"
	@echo "OBJECTS   : $(OBJECTS)"
	@echo "LIBRARY   : $(LIBRARY)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###                                                                         ###
###############################################################################

Advection.o         : Advection.cpp

Aircraft.o          : Aircraft.cpp

Ambient.o           : Ambient.cpp

BuildTime.o         : BuildTime.cpp

CallSolver.o        : CallSolver.cpp

Cluster.o           : Cluster.cpp

Combvec.o           : Combvec.cpp

Diffusion.o         : Diffusion.cpp

Emission.o          : Emission.cpp

Engine.o            : Engine.cpp

FileHandler.o       : FileHandler.cpp

Fuel.o              : Fuel.cpp

LiquidAer.o         : LiquidAer.cpp

Main.o              : Main.cpp

Mesh.o              : Mesh.cpp

PlumeModel.o        : PlumeModel.cpp

PrintMessage.o      : PrintMessage.cpp

ReadParameters.o    : ReadParameters.cpp

Reshape.o           : Reshape.cpp

Ring.o              : Ring.cpp

Save.o              : Save.cpp

Species.o           : Species.cpp

Structure.o         : Structure.cpp

SZA.o               : SZA.cpp

UpdateTime.o        : UpdateTime.cpp

Util.o              : Util.cpp

Vortex.o            : Vortex.cpp

