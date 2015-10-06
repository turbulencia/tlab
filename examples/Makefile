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
     Case09 \
     Case10 \
     Case11 \
     Case12 \
     Case13 \
     Case14 \
     Case15 \
     Case18 \
     Case21 \
     Case22 \
     Case23 \
     Case24 \
     Case25 \
     Case26 \
     Case40

     # Case16 \
     # Case17 \

FILES=avg* grid* flow* scal* dns.out dns.log dns.err dns.war dns.def dns.ini.bak tower*

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

