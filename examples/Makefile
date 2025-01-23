# BINS_PATH contains the executables to be tested; needs to be passed externally
#
PRJS=\
     Case01 \
     Case02 \
     Case03 \
     Case04 \
     Case05 \
     Case06 \
     Case07 \
     Case08 \
     Case10 \
     Case11 \
     Case12 \
     Case13 \
     Case14 \
     Case15 \
     Case16 \
     Case17 \
     Case18 \
     Case19 \
     Case20 \
     Case21 \
     Case24 \
     Case25 \
     Case26 \
     Case27 \
     Case28 \
     Case29 \
     Case30 \
     Case41 \
     Case43 \
     Case44 \
     Case46 \
     Case51 \
     Case52 \
     Case53 \
     Case70 \
     Case81 \
     Case82 \
     Case84 \
     Case85 \
     Case86 \
     Case87 \
     Case91 \
     Case92 \
     Case93

    #  Case42 
    #  Case34 

FILES=avg* grid* flow* scal* dns.out tlab.log dns.obs tlab.err tlab.war tlab.ini.bak tower* part* traj* residence* eps*

##########################################################################
#   Directives
##########################################################################
SHELL := /bin/bash

check:
	@echo -e "Running selected cases in examples directory."
	@echo -e "Reference incompressible values obtained with FFTW/ESTIMATE mode."
	@(for prj in $(PRJS); do (echo -n "[`date +%T`] Checking $$prj... "; cd $$prj; ../check_run.sh ../../$(BINS_PATH) ) done )

##########################################################################
clean: 
	@(for prj in $(PRJS); do (cd $$prj; rm -f $(FILES) ) done )

proper:
	-/bin/rm -rf *~
	@(for prj in $(PRJS); do (cd $$prj; rm -f $(FILES) ) done )

