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

# Variable    Descripton
# --------    ----------
# APCEMMDIR	  Specifies the directory where APCEMM's core routines are found
# KPPDIR      Specifies the directory where the KPP routines are found
#
# !REVISION HISTORY:
# 4 Sep 2018 - T. Fritz - Initial version

#------------------------------------------------------------------------------

# Directories
APCEMMDIR := Core
KPPDIR    := KPP

###############################################################################
###																			###
###							Makefile targets								###
###																			###
###############################################################################


all: 
	@$(MAKE) -C $(KPPDIR)
	@$(MAKE) -C $(APCEMMDIR) all

clean:
	@$(MAKE) -C $(KPPDIR) clean
	@$(MAKE) -C $(APCEMMDIR) clean
