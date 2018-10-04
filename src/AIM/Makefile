#------------------------------------------------------------------------------
#                 MIT Laboratory for Aviation and the Environment
#------------------------------------------------------------------------------
#                               C++ - APCEMM                                  #
#      A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)       #
#------------------------------------------------------------------------------
#
# !MODULE: Makefile (AIM directory)
#
# !DESCRIPTION: This makefile compiles the various AIM source files 
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
#------------------------------------------------------------------------------

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

# Directories
ROOT_DIR  :=../..
LIB_DIR   :=$(ROOT_DIR)/lib
INCL_DIR  :=$(ROOT_DIR)/include/AIM

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
LIBRARY    := libAim.a

###############################################################################
###                                                                         ###
###                         Makefile targets                                ###
###                                                                         ###
###############################################################################

.PHONY: lib clean debug

lib: $(OBJECTS)
	@$(AR) crs $(LIBRARY) $(OBJECTS)
	mv $(LIBRARY) $(LIB_DIR)
	@echo "Library $(LIBRARY) has been created!"

clean:
	@echo "===> Making clean in directory: AIM <==="
	-@rm -vf *.o

debug:
	@echo "ROOT    : $(ROOT_DIR)"
	@echo "LIB     : $(LIB_DIR)"
	@echo "GCC     : $(GCC)"
	@echo "SOURCES : $(SOURCES)"
	@echo "OBJECTS : $(OBJECTS)"
	@echo "LIBRARY : $(LIBRARY)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###                                                                         ###
###############################################################################

Aerosol.o           : Aerosol.cpp

buildKernel.o       : buildKernel.cpp

Coagulation.o       : Coagulation.cpp

Nucleation.o        : Nucleation.cpp

