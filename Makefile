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
# 04 Sep 2018 - T. Fritz - Initial version
# 27 Sep 2018 - T. Fritz - Added MAKE for all modules
#------------------------------------------------------------------------------

# Directories
APCEMMDIR := Core
KPPDIR    := KPP    # (Kinetics Pre-Processor                  -> Chemistry)
SANDSDIR  := SANDS  # (Spectral Advection aNd Diffusion Solver -> Atmospheric Advection/Diffusion)
AIMDIR    := AIM    # (AIrcraft Microphysics                   -> Microphysics)
EPMDIR    := EPM    # (Early Plume Microphysics                -> Early Plume Microphysics)

###############################################################################
###																			###
###							Makefile targets								###
###																			###
###############################################################################


all: 
	@$(MAKE) -C $(KPPDIR)
	@$(MAKE) -C $(SANDSDIR)
	@$(MAKE) -C $(AIM)
	@$(MAKE) -C $(EPM)
	@$(MAKE) -C $(APCEMMDIR) all

clean:
	@$(MAKE) -C $(KPPDIR) clean
	@$(MAKE) -C $(SANDSDIR) clean
	@$(MAKE) -C $(AIM) clean
	@$(MAKE) -C $(EPM) clean
	@$(MAKE) -C $(APCEMMDIR) clean
