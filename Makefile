#------------------------------------------------------------------------------
#                 MIT Laboratory for Aviation and the Environment
#------------------------------------------------------------------------------
#                               C++ - APCEMM                                  #
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

# Variable    Descripton
# --------    ----------
#
# !REVISION HISTORY:
# 04 Sep 2018 - T. Fritz - Initial version
# 27 Sep 2018 - T. Fritz - Added MAKE for all modules
# 03 Oct 2018 - T. Fritz - Architecture change
#------------------------------------------------------------------------------

# Directories
SRC_DIR    :=src
APCEMM_DIR :=$(SRC_DIR)/Core
LIB_DIR    :=lib
APP_DIR    :=object
###############################################################################
###                                                                         ###
###                         Makefile targets                                ###
###                                                                         ###
###############################################################################

.PHONY: all build exe clean realclean debug wipeout help headerinfo
.PHONY: lib libCore libUtil libKpp libSands libAim libEpm

all:
	@$(MAKE) -C $(APCEMM_DIR) all

build:
	@$(MAKE) -C $(APCEMM_DIR) build

lib:
	@$(MAKE) -C $(APCEMM_DIR) lib

libCore:
	@$(MAKE) -C $(APCEMM_DIR) libCore

libUtil:
	@$(MAKE) -C $(APCEMM_DIR) libUtil

libKpp:
	@$(MAKE) -C $(APCEMM_DIR) libKpp
	
libSands:
	@$(MAKE) -C $(APCEMM_DIR) libSands

libAim:
	@$(MAKE) -C $(APCEMM_DIR) libAim

libEpm:
	@$(MAKE) -C $(APCEMM_DIR) libEpm

exe:
	@$(MAKE) -C $(APCEMM_DIR) exe

clean:
	@$(MAKE) -C $(APCEMM_DIR) clean

realclean:
	@$(MAKE) -C $(APCEMM_DIR) realclean

debug:
	@$(MAKE) -C $(APCEMM_DIR) debug

wipeout:
	@$(MAKE) -C $(APCEMM_DIR) wipeout

help:
	@$(MAKE) -C $(APCEMM_DIR) help

headerinfo:
	@$(MAKE) -C $(APCEMM_DIR) headerinfo

