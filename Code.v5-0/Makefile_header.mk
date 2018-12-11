#------------------------------------------------------------------------------
#                 MIT Laboratory for Aviation and the Environment
#------------------------------------------------------------------------------
#      							C++ - APCEMM                                  #
#      A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)       #
#------------------------------------------------------------------------------
#
# !IROUTINE: Makefile_header.mk
#
# !DESCRIPTION: This sub-makefile defines the variables which specify 
# compilation options for the different support compiler/platform combinations.
# Also, the default makefile compilation rules are specified here.
#
# NOTE: We now use SHELL :=/bin/bash as the default Unix shell.
#
# !REVISION HISTORY:
# 03 Oct 2018 - T. Fritz - Initial version
# 10 Oct 2018 - T. Fritz - Added Makefile flags
# 28 Nov 2018 - T. Fritz - Fixed typo
# 						   Added optimization flags
# 						   Removed flags unrelevant to C/C++ code
#------------------------------------------------------------------------------

###############################################################################
###                                                                         ###
###  Set the default Unix shell and some error message variables            ###
###                                                                         ###
###############################################################################

# Set default shell to bash, for use with the Makefile conditionals
SHELL                :=/bin/bash

# Error message for bad COMPILER input
ERR_CMPLR            :="Select a compiler: COMPILER=g++, COMPILER=gcc"

# Error message for unknown compiler/OS combintation
ERR_OSCOMP           :="Makefile_header.mk not set up for this compiler/OS combination"

###############################################################################
###                                                                         ###
###  Set C-preprocessor switches representing user options.  These are not  ###
###  specific to any compiler, but are general options for the simulation.  ###
###                                                                         ###
###############################################################################

#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------

# %%%%% OpenMP parallelization (on by default) %%%%%
ifndef OMP
  OMP                :=0
endif

# %%%%% Set default compiler %%%%%

# %%%%% If COMPILER is not defined, default to the $(CXX) variable, which %%%%%
# %%%%% is set in your .bashrc, or when you load the compiler module      %%%%%
ifndef COMPILER
  COMPILER           :=$(CXX)
endif

# %%%%% Test if g++ is selected %%%%%
REGEXP               :=(^[Gg][+][+]|[Cc][Xx][Xx])
ifeq ($(shell [[ "$(COMPILER)" =~ $(REGEXP) ]] && echo true),true)
  COMPILE_CMD        :=$(CXX)
  USER_DEFS          += -DLINUX_CXX
endif

# %%%%% Test if g++ is selected %%%%%
REGEXP               :=(^[Gg][Cc][Cc]|^[Cc][Cc])
ifeq ($(shell [[ "$(COMPILER)" =~ $(REGEXP) ]] && echo true),true)
  COMPILE_CMD        :=$(CC)
  USER_DEFS          += -DLINUX_CC
endif

# %%%%% ERROR CHECK!  Make sure our COMPILER selection is valid! %%%%%
REGEXP               :=((-DLINUX_)?|CXX|CC)
ifneq ($(shell [[ "$(USER_DEFS)" =~ $(REGEXP) ]] && echo true),true)
  $(error $(ERR_CMPLR))
endif

#------------------------------------------------------------------------------
# Diagnostic settings
#------------------------------------------------------------------------------

# %%%%% NETCDF (for using new netCDF diagnostic output) %%%%%
REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(NC_DIAG)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DNC_DIAG
  IS_NC_DIAG         :=1
endif

#------------------------------------------------------------------------------
# Add test for mass conservation
#------------------------------------------------------------------------------

REGEXP               :=(^[Yy]|^[Yy][Ee][Ss])
ifeq ($(shell [[ "$(MASSCONS)" =~ $(REGEXP) ]] && echo true),true)
  USER_DEFS          += -DMASSCONS
endif

###############################################################################
###                                                                         ###
###  Set linker commands for local and external libraries (incl. netCDF)    ###
###                                                                         ###
###############################################################################

# netCDF Library include path.  

# Get the version number (e.g. "4130"=netCDF 4.1.3; "4200"=netCDF 4.2, etc.)
NC_VERSION := $(shell nc-config --version)
NC_VERSION := $(shell echo "$(NC_VERSION)" | sed 's|netCDF ||g')
NC_VERSION := $(shell echo "$(NC_VERSION)" | sed 's|\.||g')
NC_VERSION_LEN       :=$(shell perl -e "print length $(NC_VERSION)")

ifeq ($(NC_VERSION_LEN),3)
 NC_VERSION          :=$(NC_VERSION)0
endif
ifeq ($(NC_VERSION_LEN),2) 
 NC_VERSION          :=$(NC_VERSION)00
endif

# Add path/to/lib folder
LINK := $(LINK)

# Add path to APCEMM's library folder
LINK := -L$(LIB_DIR)

# Define any libraries to link into executable: use -llibname option
LINK := $(LINK) -lstdc++ -lm

LINK_FFTW := -lfftw3 -lfftw3f -lfftw3l
# -lfftw3_omp -lfftw3_threads

# Create linker command to create the APCEMM executable
LINK := $(LINK) -Wl,--start-group -lUtil -lSands $(LINK_FFTW) -lKpp -lEpm -lAim -lnetcdf_c++ -Wl,--end-group

###############################################################################
###                                                                         ###
###  Define settings for the GNU C++ COMPILER (aka g++)                     ###
###                                                                         ###
###############################################################################

ifeq ($(COMPILER),g++) 

  # Get the GNU Fortran version
  VERSIONTEXT        :=$(shell $(CXX) --version)
  VERSION            :=$(word 4, $(VERSIONTEXT))
  VERSION            :=$(subst .,,$(VERSION))

  # Base set of compiler flags
  CXXFLAGS            :=-std=c++11 -w

  # Default optimization level for all routines (-O3)
  ifndef OPT
    # Options of interest
    #  -limf                Intel math libraries - machine must have them
    #  -O3                  Highest safe optimization level
    #  -march=native        Make the binary machine-specific. If in doubt, 
    #                        use a specific architecture, eg...
    #  -march=corei7-avx    Binary uses optimizations for 
    #                        Intel Sandy-Bridge Xeon (e.g. E5-2680)
    #  -mfpmath=sse         Use SSE extensions
    #  -funroll-loops       Enable loop unrolling
    OPT              := -O3 -funroll-loops -limf -march=native -msse
  endif
  
  # Pick compiler options for debug run or regular run 
  REGEXP             := (^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(DEBUG)" =~ $(REGEXP) ]] && echo true),true)
    CXXFLAGS          += -Wall -Wextra -Wconversion
    CXXFLAGS          += -Wno-c++0x-compat
    CXXFLAGS          += -Wno-unused
    CXXFLAGS          += -g -O0
    USER_DEFS         += -DDEBUG
  else
    CXXFLAGS          += $(OPT)
  endif
  
  # Turn on OpenMP parallelization
  REGEXP             :=(^[Yy]|^[Yy][Ee][Ss])
  ifeq ($(shell [[ "$(OMP)" =~ $(REGEXP) ]] && echo true),true)
    CXXFLAGS          += -fopenmp
    USER_DEFS         += -DOMP
  endif
 
  # Append the user options in USER_DEFS to CXXFLAGS
  CXXFLAGS            += $(USER_DEFS)

  # Include options (i.e. for finding *.h* files)
  INCLUDE := -I$(ROOT_DIR)/include

endif

# Set the standard compiler variables
GCC  := $(COMPILE_CMD) $(CXXFLAGS) $(INCLUDE) $(LINK)
LD   := $(COMPILE_CMD) $(CXXFLAGS)

###############################################################################
###                                                                         ###
###  Specify pattern rules for compiliation                                 ###
###  (i.e. tell "make" how to compile files w/ different extensions)        ###
###                                                                         ###
###############################################################################

%.o : %.cpp
	$(GCC) -o $@ -c $< 
%.o : %.c
	$(GCC) -o $@ -c $<


###############################################################################
###                                                                         ###
###  Set Makefile flags representing user options                           ###
###                                                                         ###
###############################################################################

ifndef VERBOSE
  MAKEFLAGS += --no-print-directory
endif

###############################################################################
###                                                                         ###
###  Export global variables so that the main Makefile will see these       ###
###                                                                         ###
###############################################################################

export OMP
export GCC
export LD
export INCLUDE
export LINK
export SHELL
export MAKEFLAGS

###############################################################################
###                                                                         ###
###  Debug print output.  Normally you will leave the following lines       ###
###  commented out.  Uncomment these lines for testing.                     ###
###                                                                         ###
###############################################################################

#headerinfo:
#	@@echo '####### in Makefile_header.mk ########'
#	@@echo "COMPILER    : $(COMPILER)"
#	@@echo "DEBUG       : $(DEBUG)"
#	@@echo "CXX         : $(CXX)"
#	@@echo "CC          : $(CC)"
#	@@echo "INCLUDE     : $(INCLUDE)"
#	@@echo "LINK        : $(LINK)"
#	@@echo "USERDEFS    : $(USER_DEFS)"
#	@@echo "FLAGS       : $(CXXFLAGS)"

