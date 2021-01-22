#------------------------------------------------------------------------------
#                 MIT Laboratory for Aviation and the Environment
#------------------------------------------------------------------------------
#                               C++ - APCEMM                                  #
#      A(ircraft) P(lume) C(hemistry) E(mission) M(icrophysics) M(odel)       #
#------------------------------------------------------------------------------
#
# !MODULE: Makefile
#
# !DESCRIPTION: Makefile for APCEMM run directories
# 
# !REMARKS:
#  (1) CODE_DIR (source code location), VERSION and LOG_DIR are automatically
#  	   passed to the make command using the default values defined in this file.
#
# !REVISION HISTORY:
# 11 Dec 2018 - T. Fritz - Initial version
# 27 Aug 2020 - T. Fritz - Rename C-APCEMM into APCEMM
#
#------------------------------------------------------------------------------

# Unix shell
SHELL          :=/bin/bash

###############################################################################
###                                                                         ###
###     CONFIGURABLE TOKENS: You can modify these for your environment      ###
###                                                                         ###
###############################################################################

# APCEMM run log filename prefix
ifndef VERSION
  VERSION      :=v05-00
endif

# Source code location
ifndef CODE_DIR
  CODE_DIR       :=/net/d13/data/aa681/statistical_global_contrails/Simulations/C-APCEMM/Code.$(VERSION)
  # CODE_DIR     :=$(HOME)/APCEMM/Code.$(VERSION)
endif

# Output log file destination (default is run directory)
ifndef LOG_DIR
  LOG_DIR      :=.
endif

###############################################################################
###                                                                         ###
###     DEFAULT COMPILATION OPTIONS: Set various switches if not defined    ###
###                                                                         ###
###############################################################################

ifndef TRACEBACK
  TRACEBACK    :=n
endif

ifndef BOUNDS
  BOUNDS       :=n
endif

ifndef DEBUG
  DEBUG        :=n
endif

ifndef TIMERS
  TIMERS       :=n
endif

ifndef RINGS
  RINGS        :=n
endif

###############################################################################
###                                                                         ###
###  RUN OPTIONS: Get various settings from the run directory name,         ###
###	 			  or from the type of simulation that is being used.        ###
###                                                                         ###
###############################################################################

# Run directory path
RUN_DIR        :=$(shell pwd)

# getRunInfo perl script location (default is run directory)
ifndef PERL_DIR
  PERL_DIR     :=$(RUN_DIR)
endif

# TODO
# Get run info... ~!@# To be implemented!
#

# General run information
TIMESTAMP    :=$(shell date +"%Y/%m/%d %H:%M")	

# Log files that will be written to the log directory
LOG_COMP     :="$(RUN_DIR)/lastbuild"
LOG_SP       :="$(LOG_DIR)/$(VERSION).log.sp"
LOG_MP       :="$(LOG_DIR)/$(VERSION).log.mp"

# Executables
EXE_SP       :=APCEMM
EXE_MP       :=APCEMM

# Get information about the Git version, because some features
# (like -C) are not available in older Git versions.  The git 
# version command returns "git version X.Y.Z", so we will just take
# the 3rd word and strip all the dots. (bmy, 12/21/16)
GIT_VERSION    :=$(subst .,,$(word 3, $(shell git --version)))
NEWER_THAN_185 :=$(shell perl -e "print ($(GIT_VERSION) gt 185)")

ifeq ($(NEWER_THAN_185),1)

 # Git version 1.8.5 and higher uses the -C option to look in
 # directories other than the current directory.  Use this to
 # get info about the last committed state of the code.
 CODE_BRANCH  :=$(shell git -C $(CODE_DIR) rev-parse --abbrev-ref HEAD)
 LAST_COMMIT  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%s") 
 COMMIT_DATE  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%cd") 
 COMMIT_USER  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%cn") 
 COMMIT_HASH  :=$(shell git -C $(CODE_DIR) log -n 1 --pretty=format:"%h") 
 CODE_STATUS  :=$(shell git -C $(CODE_DIR) status -s) 

else

 # Git versions older than 1.8.5 don't have the -C option,
 # so we have to manually change to the code dir to get information
 # about the last committed state of the code. (bmy, 12/21/16)
 CODE_BRANCH  :=$(shell cd $(CODE_DIR); git rev-parse --abbrev-ref HEAD; cd $(PWD))
 LAST_COMMIT  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%s"; cd $(PWD)) 
 COMMIT_DATE  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%cd"; cd $(PWD)) 
 COMMIT_USER  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%cn"; cd $(PWD)) 
 COMMIT_HASH  :=$(shell cd $(CODE_DIR); git log -n 1 --pretty=format:"%h"; cd $(PWD)) 
 CODE_STATUS  :=$(shell cd $(CODE_DIR); git status -s; cd $(PWD)) 

endif

# If COMPILER is not defined, default to the $(CXX) variable, which
# is set in your .bashrc, or when you load the compiler module
ifndef COMPILER
  COMPILER           :=$(CXX)
endif

###############################################################################
###                                                                         ###
###                            MAKEFILE TARGETS                             ###
###                                                                         ###
###############################################################################

# PHONY targets don't actually compile anything. They either are
# synonyms for other targets, they remove files, or they print output.
.PHONY: all     dataclean  logclean   fileclean  spclean
.PHONY: mpclean spexeclean mpexeclean superclean printallinfo

all: unittest

ut: unittest

unittest:
	@$(MAKE) realclean
	@$(MAKE) sp
	@$(MAKE) realclean
	@$(MAKE) mp
	

#########################
## Build and run       ##
#########################

sp:
	@$(MAKE) spclean
	@$(MAKE) -C $(CODE_DIR) OMP=no RINGS=$(RINGS) > $(LOG_SP)
	@$(MAKE) printbuildinfo > $(LOG_COMP).sp
	cp -f $(CODE_DIR)/build/apps/APCEMM.sh $(EXE_SP)
	./$(EXE_SP) >> $(LOG_SP)
	@$(MAKE) printruninfosp >> $(LOG_SP)
	@$(MAKE) printallinfosp

mp:
	@$(MAKE) mpclean
	@$(MAKE) -C $(CODE_DIR) OMP=yes RINGS=$(RINGS) > $(LOG_MP)
	@$(MAKE) printbuildinfo > $(LOG_COMP).mp
	cp -f $(CODE_DIR)/build/apps/APCEMM.sh $(EXE_MP)
	./$(EXE_MP) >> $(LOG_MP)
	@$(MAKE) printruninfomp >> $(LOG_MP)
	@$(MAKE) printallinfomp

#########################
## Build only          ##
#########################

spbuild:
	@$(MAKE) spexeclean
	@$(MAKE) splogclean
	@$(MAKE) -C $(CODE_DIR) OMP=no RINGS=$(RINGS)
	cp -f $(CODE_DIR)/build/apps/APCEMM.sh $(EXE_SP)
	@$(MAKE) printbuildinfo > $(LOG_COMP).sp
	@$(MAKE) printbuildinfo

mpbuild:
	@$(MAKE) mpexeclean
	@$(MAKE) mplogclean
	@$(MAKE) -C $(CODE_DIR) OMP=yes RINGS=$(RINGS)
	cp -f $(CODE_DIR)/build/apps/APCEMM.sh $(EXE_MP)
	@$(MAKE) printbuildinfo > $(LOG_COMP).mp
	@$(MAKE) printbuildinfo

#########################
## Run only            ##
#########################

sprun:
	@$(MAKE) spdataclean
	./$(EXE_SP) > $(LOG_SP)
	@$(MAKE) printruninfosp >> $(LOG_SP)
	@$(MAKE) printruninfosp

mprun:
	@$(MAKE) mpdataclean
	./$(EXE_MP) > $(LOG_MP)
	@$(MAKE) printruninfomp >> $(LOG_MP)
	@$(MAKE) printruninfomp

#########################
## Remove output data  ##
#########################

dataclean:

spdataclean:

mpdataclean:

########################
## Remove logs        ##
########################

logclean:
	@$(MAKE) splogclean
	@$(MAKE) mplogclean
	rm -f *.log

splogclean:
	rm -f $(LOG_SP)

mplogclean:
	rm -f $(LOG_MP)

########################
## Remove executables ##
########################

execlean:
	@$(MAKE) spexeclean
	@$(MAKE) mpexeclean

spexeclean:
	rm -f $(EXE_SP)
	rm -f $(LOG_COMP).sp
	
mpexeclean:
	rm -f $(EXE_MP)
	rm -f $(LOG_COMP).mp

########################
## Remove all         ##
########################

fileclean: dataclean logclean execlean

spclean: spdataclean splogclean spexeclean

mpclean: mpdataclean mplogclean mpexeclean

#########################
## Remove config files ##
#########################

extrafileclean:
	rm -f input.apcemm

########################
## Clean source       ##
########################

clean:
	@$(MAKE) -C $(CODE_DIR) clean

realclean:
	@$(MAKE) -C $(CODE_DIR) realclean
	@$(MAKE) execlean
	rm -f lastbuild.*

########################
## Clean and rm all   ##
########################

superclean: fileclean realclean

########################
## Print information  ##
########################

printruninfosp:
	@echo "RUN SETTINGS:"
	@echo "  Run directory       : $(RUN_DIR)"
	@echo "  Version             : $(VERSION)"
	@echo "  APCEMM Run Log File : $(LOG_SP)"

printruninfomp:
	@echo "RUN SETTINGS:"
	@echo "  Run directory       : $(RUN_DIR)"
	@echo "  Version             : $(VERSION)"
	@echo "  APCEMM Run Log File : $(LOG_MP)"

printbuildinfo:
	@echo "LAST BUILD INFORMATION:"
	@echo "  CODE_DIR     : $(CODE_DIR)"
	@echo "  CODE_BRANCH  : $(CODE_BRANCH)"
	@echo "  LAST_COMMIT  : $(LAST_COMMIT)"
	@echo "  COMMIT_DATE  : $(COMMIT_DATE)"
	@echo "  VERSION      : $(VERSION)"
	@echo "  TRACEBACK    : $(TRACEBACK)"
	@echo "  BOUNDS       : $(BOUNDS)"
	@echo "  DEBUG        : $(DEBUG)"
	@echo "  TIMERS       : $(TIMERS)"
	@echo "  COMPILER     : $(COMPILER) $(COMPILER_VERSION)"
	@echo "  Datetime     : $(TIMESTAMP)"

printallinfosp: 
	@$(MAKE) printbuildinfo
	@$(MAKE) printruninfosp

printallinfomp: 
	@$(MAKE) printbuildinfo
	@$(MAKE) printruninfomp

printcodeinfo:
	@echo -e "Code directory:  $(CODE_DIR)"
	@echo -e "Git branch:      $(CODE_BRANCH)"
	@echo -e "Last commit:"
	@echo -e "   Message:      $(LAST_COMMIT)"
	@echo -e "   Date:         $(COMMIT_DATE)"
	@echo -e "   Committer:    $(COMMIT_USER)"
	@echo -e "   Hash abbrev:  $(COMMIT_HASH)"
	@echo -e "Uncommitted files (if any):\n$(CODE_STATUS)"

###########################################################################
###                                                                     ###
###                             HELP SCREEN                             ###
###                                                                     ###
###########################################################################

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%     APCEMM Run Directory Makefile Help Screen      %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo
	@echo 'Usage: make -jN TARGET [OPTIONAL-FLAGS]'
	@echo ''
	@echo '-jN               Compiles N files at a time (reduces compile time)'
	@echo ''
	@echo '----------------------------------------------------------'
	@echo 'TARGET may be one of the following:'
	@echo '----------------------------------------------------------'
	@echo 'all               Default target (synonym for "unittest")'
	@echo ''
	@echo '%%% COMPILE AND RUN %%%'
	@echo 'unittest          Runs a APCEMM unit test:'
	@echo '                  (1) Builds & runs APCEMM with parallelization OFF;'
	@echo '                  (2) Builds & runs APCEMM with parallelzation ON;'
	@echo '                  (3) Checks to see if the output is identical.'
	@echo 'sp                Builds and runs APCEMM with parallelization OFF.'
	@echo 'mp                Builds and runs APCEMM with parallelization ON.'
	@echo ''
	@echo '%%% BUILD ONLY %%%'
	@echo 'spbuild           Just builds APCEMM with parallelization OFF.'
	@echo 'mpbuild           Just builds APCEMM with parallelization ON.'
	@echo ''
	@echo '%%% RUN ONLY %%%'
	@echo 'sprun             Just runs APCEMM with parallelization OFF.'
	@echo 'mprun             Just runs GEOS_Chem with parallelization ON.'  
	@echo ''
	@echo '%%% REMOVE DATA FILES %%%'
	@echo 'dataclean         Removes ALL APCEMM diagnostic and restart files'
	@echo 'spdataclean       Removes diagnostic and restart files from APCEMM'
	@echo '                    simulations with parallelization turned OFF.'
	@echo 'mpdataclean       Removes diagnostic and restart files from APCEMM'
	@echo '                    simulations with parallelization turned ON.'
	@echo ''
	@echo '%%% REMOVE LOG FILES %%%'
	@echo 'logclean          Removes all APCEMM and HEMCO output log files.'
	@echo 'splogclean        Removes APCEMM and HEMCO log files from APCEMM'
	@echo '                    simulations with parallelization turned OFF.'
	@echo 'mplogclean        Removes APCEMM and HEMCO log files from APCEMMs'
	@echo '                    simulations with parallelization turned ON.'
	@echo ''
	@echo '%%% REMOVE EXECUTABLE FILES %%%'
	@echo 'execlean          Removes all APCEMM executable files'
	@echo 'spexeclean        Removes APCEMM executable files for which the'
	@echo '                    parallelization has been turned OFF.'
	@echo 'mpexeclean        Removes APCEMM executable files for which the'
	@echo '                    parallelization has been turned ON.'
	@echo ''
	@echo '%%% REMOVE ALL FILES %%%'
	@echo 'fileclean         Synonym for: dataclean   logclean   execlean'
	@echo 'spclean           Synonym for: spdataclean splogclean spexeclean'
	@echo 'mpclean           Synonym for: mpdataclean mplogclean mpexeclean'
	@echo ''
	@echo '%%% REMOVE CONFIG FILES %%%'
	@echo 'extrafileclean  Removes input.geos and HEMCO_Config.rc files'
	@echo ''
	@echo '%%% CLEAN SOURCE CODE %%%'
	@echo 'clean             Makes "clean" in source code directory $CODE_DIR'
	@echo 'realclean         Makes "realclean" in the source code directory $CODE_DIR'
	@echo ''
	@echo '%%% CLEAN AND REMOVE ALL %%%'
	@echo 'superclean        Synonym for: fileclean realclean'
	@echo ''
	@echo '%%% PRINT INFORMATION %%%'
	@echo 'printruninfosp    Prints the run settings for APCEMM simulations'
	@echo '                    for which the parallelization is turned OFF.'
	@echo 'printruninfomp    Prints the run settings for APCEMM simulations'
	@echo '                    for which the parallelization is turned ON'
	@echo 'printbuildinfo    Prints the build settings for APCEMM simulations'
	@echo 'printallinfosp    Synonym for: printbuildinfosp printruninfosp'
	@echo 'printallinfomp    Synonym for: printbuildinfomp printruninfomp'
	@echo 'printcodeinfo     Print code directory git information'
	@echo ''
	@echo 'help              Prints this help screen'
	@echo ''
	@echo '----------------------------------------------------------'
	@echo 'OPTIONAL-FLAGS may be one of the following:'
	@echo '----------------------------------------------------------'
	@echo 'DEBUG=y           Compiles APCEMM with various debugging options'
	@echo 'BOUNDS=y          Compiles APCEMM with out-of-bounds error checks'
	@echo 'TRACEBACK=y       Compiles APCEMM with traceback error printout'
	@echo 'TIMERS=1          Compiles APCEMM with timers turned on'

